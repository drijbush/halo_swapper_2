Module halo_parallel_module

  ! ToDo:
  ! 1) Rationalise error codes
  ! 2) Add flags for orthogonal cells which avoid unneccesary comms

  Use mpi_f08, Only : mpi_comm, mpi_request

  Use constants, Only : wp
  Use swap_module

  Use halo_setter_base_module, Only : halo_setter_base_class

  Implicit None

  Type, Private :: halo_comms
     Integer,    Dimension( 1:3 )                  :: remote_coord
     Integer                                       :: remote_rank
     Integer,    Dimension( 1:3 )                  :: comm_size
     Integer,    Dimension( 1:3 )                  :: comm_start
     Integer                                       :: tag
     Logical                                       :: is_local
     Real( wp ), Dimension( :, :, : ), Allocatable :: buffer
  End type halo_comms

  Type, Public, Extends( halo_setter_base_class ) :: halo_parallel_setter_2
     Private
     Type( mpi_comm ),                                   Private :: comm
     Integer,                                            Private :: halo_width
     Integer,             Dimension( 1:3 ),              Private :: local_size
     Integer,             Dimension( 1:3 ),              Private :: total_size
     Integer,             Dimension( 1:3 ),              Private :: first_point
     Integer,             Dimension( 1:3 ),              Private :: n_procs
     Integer,             Dimension( 1:3 ),              Private :: my_coords
     Logical,             Dimension( 1:3 ),              Private :: is_periodic
     Type( mpi_request ), Dimension(  :  ), Allocatable, Private :: msg_requests
     Type( halo_comms ),  Dimension(  :  ), Allocatable, Private :: recv_comms
     Type( halo_comms ),  Dimension(  :  ), Allocatable, Private :: send_comms
   Contains
     Generic,   Public  :: init => halo_parallel_init_old
     Generic,   Public  :: init => halo_parallel_init_f08
     Procedure, Public  :: fill => halo_fill
     Procedure, Private :: halo_parallel_init_old
     Procedure, Private :: halo_parallel_init_f08
  End Type halo_parallel_setter_2

  Private

Contains

  Subroutine halo_parallel_init_f08( H, local_size, halo_width, comm, error )

    ! NO OPTIMISATION DUE TO ORTHOG GRIDS
    ! ONLY WORKS for halo_width <= local_size - factor_comms needs more thought

    Use mpi_f08, Only : mpi_comm, mpi_topo_test, mpi_cart, mpi_cartdim_get, mpi_cart_get, &
         mpi_cart_sub, mpi_comm_free, mpi_comm_rank, &
         mpi_integer, mpi_allreduce, mpi_max, mpi_min, mpi_in_place

    ! GRIDS START AT ZERO

    Class( halo_parallel_setter_2 ),                   Intent( InOut ) :: H
    Integer,                         Dimension( 1:3 ), Intent( In    ) :: local_size
    Integer,                                           Intent( In    ) :: halo_width
    Type( mpi_comm ),                                  Intent( In    ) :: comm
    Integer,                                           Intent(   Out ) :: error

    Type one_d_comms
       Integer, Dimension( : ), Allocatable :: coords
       Integer, Dimension( : ), Allocatable :: sizes
       Integer, Dimension( : ), Allocatable :: starts
    End type one_d_comms

    Type( one_d_comms ), Dimension( 1:3 ) :: one_d_send, one_d_recv

    Type( mpi_comm ) :: axis_comm
    Type( mpi_comm ) :: plane_comm

    Integer :: comm_type
    Integer :: ndims
    Integer :: n_loc_max, n_loc_min
    Integer :: i

    Logical, Dimension( 1:3 ) :: is_this_axis
    Logical, Dimension( 1:3 ) :: is_this_orthog_plane

    error = 0

    ! Check it is a caretesian comunicator
    Call mpi_topo_test( comm, comm_type )
    If( comm_type /= mpi_cart ) Then
       error = 1
       Return
    End If

    ! Check is has 3 dimensions
    Call mpi_cartdim_get( comm, ndims )
    If( ndims /= 3 ) Then
       error = 2
       Return
    End If

    ! Check it is periodic all in directions
    Call mpi_cart_get( comm, ndims, H%n_procs, H%is_periodic, H%my_coords )
    If( .Not. All ( H%is_periodic ) ) Then
       error = 3
       Return
    End If

    ! Check consistency of sizes in each direction - across a given plane
    ! each process must have the same size as otherwise it is not a cartesian grid
    Do i = 1, 3
       ! Pick an axis
       is_this_axis = .False.
       is_this_axis( i ) = .True.
       ! Take the plane orthogonal to it
       is_this_orthog_plane = .Not. is_this_axis
       ! Create a communicator containg this process in that plane
       Call mpi_cart_sub( comm, is_this_orthog_plane, plane_comm )
       ! Find the maximum and minimum value of the local size in that plane
       Call mpi_allreduce( local_size( i ), n_loc_max, 1, mpi_integer, mpi_max, plane_comm )
       Call mpi_allreduce( local_size( i ), n_loc_min, 1, mpi_integer, mpi_min, plane_comm )
       ! If the max size is not the same as the min size that means not all procs
       ! in the plane have the same size
       If( n_loc_max /= n_loc_min ) Then
          error = 4
       End If
    End Do
    ! Check if an error occured
    Call mpi_allreduce( mpi_in_place, error, 1, mpi_integer, mpi_max, comm )
    If( error /= 0 ) Then
       Return
    End If

    H%comm       = comm
    H%local_size = local_size
    H%halo_width = halo_width

    ! In x, y, z directions work out the "factored" communications patters
    Do i = 1, 3
       ! Create a comunicator for the direction in question
       is_this_axis = .False.
       is_this_axis( i ) = .True.
       Call mpi_cart_sub( comm, is_this_axis, axis_comm )
       Call factor_comms( axis_comm, H%n_procs( i ), H%is_periodic( i ), &
            H%my_coords( i ), H%local_size( i ), H%halo_width, &
            H%total_size( i ), H%first_point( i ), one_d_recv( i ), one_d_send( i ) )
       Call mpi_comm_free( axis_comm, error )
    End Do

    Call combine_comms( H, one_d_recv, one_d_send )

  Contains

    Subroutine combine_comms( H, one_d_recv, one_d_send )

      Use mpi_f08, Only : mpi_cart_rank

      Class( halo_parallel_setter_2 ),                   Intent( InOut ) :: H
      Type( one_d_comms ),             Dimension( 1:3 ), Intent( In    ) :: one_d_recv
      Type( one_d_comms ),             Dimension( 1:3 ), Intent( In    ) :: one_d_send

      Integer, Dimension( 1:3 ) :: remote_coords
      Integer, Dimension( 1:3 ) :: comm_vec
      Integer, Dimension( 1:3 ) :: comm_size
      Integer, Dimension( 1:3 ) :: comm_start

      Integer :: remote_rank, my_rank
      Integer :: n_comms, i_comms
      Integer :: range_c
      Integer :: ix, iy, iz
      Integer :: i

      Call mpi_cart_rank( H%comm, H%my_coords, my_rank )

      ! Find the max abs coord difference for tag generation later
      range_c = -1
      Do  i = 1, 3
         range_c = Max( range_c, Maxval( Abs( one_d_recv( i )%coords( : ) - H%my_coords( i ) ) ) )
         range_c = Max( range_c, Maxval( Abs( one_d_send( i )%coords( : ) - H%my_coords( i ) ) ) )
      End Do

      ! MUST be same number of sends and recvs
      n_comms = Size( one_d_recv( 1 )%sizes ) * Size( one_d_recv( 2 )%sizes ) * Size( one_d_recv( 3 )%sizes )

      Allocate( H%recv_comms( 1:n_comms ) )
      Allocate( H%send_comms( 1:n_comms ) )

      ! First the recvs
      i_comms = 0
      Do iz = 1, Size( one_d_recv( 3 )%sizes )
         Do iy = 1, Size( one_d_recv( 2 )%sizes )
            Do ix = 1, Size( one_d_recv( 1 )%sizes )

               i_comms = i_comms + 1

               ! Note mpi_cart_rank respect periodic boundary conditions correctly
               remote_coords = [ one_d_recv( 1 )%coords( ix ), &
                    one_d_recv( 2 )%coords( iy ), &
                    one_d_recv( 3 )%coords( iz ) ]
               Call mpi_cart_rank( H%comm, remote_coords, remote_rank )

               ! comm vec is the data direction for a recv
               comm_vec = remote_coords - H%my_coords

               comm_size   = [ one_d_recv( 1 )%sizes( ix ),  one_d_recv( 2 )%sizes( iy ),  one_d_recv( 3 )%sizes( iz ) ]
               comm_start  = [ one_d_recv( 1 )%starts( ix ), one_d_recv( 2 )%starts( iy ), one_d_recv( 3 )%starts( iz ) ]

               H%recv_comms( i_comms )%remote_coord = remote_coords
               H%recv_comms( i_comms )%remote_rank  = remote_rank
               H%recv_comms( i_comms )%comm_size    = comm_size
               H%recv_comms( i_comms )%comm_start   = comm_start
               H%recv_comms( i_comms )%tag          = comm_vec_to_tag( range_c, comm_vec )
               H%recv_comms( i_comms )%is_local     = remote_rank == my_rank
               Allocate( H%recv_comms( i_comms )%buffer( 0:H%recv_comms( i_comms )%comm_size( 1 ) - 1, &
                    0:H%recv_comms( i_comms )%comm_size( 2 ) - 1, &
                    0:H%recv_comms( i_comms )%comm_size( 3 ) - 1 ) )

            End Do
         End Do
      End Do

      ! Now the sends
      i_comms = 0
      Do iz = 1, Size( one_d_send( 3 )%sizes )
         Do iy = 1, Size( one_d_send( 2 )%sizes )
            Do ix = 1, Size( one_d_send( 1 )%sizes )

               i_comms = i_comms + 1

               ! Note mpi_cart_rank respect periodic boundary conditions correctly
               remote_coords = [ one_d_send( 1 )%coords( ix ), &
                    one_d_send( 2 )%coords( iy ), &
                    one_d_send( 3 )%coords( iz ) ]
               Call mpi_cart_rank( H%comm, remote_coords, remote_rank )

               ! comm vec is the data direction for a recv, this is a send so negate it
               comm_vec = remote_coords - H%my_coords
               comm_vec = - comm_vec

               comm_size   = [ one_d_send( 1 )%sizes( ix ),  one_d_send( 2 )%sizes( iy ),  one_d_send( 3 )%sizes( iz ) ]
               comm_start  = [ one_d_send( 1 )%starts( ix ), one_d_send( 2 )%starts( iy ), one_d_send( 3 )%starts( iz ) ]

               H%send_comms( i_comms )%remote_coord = remote_coords
               H%send_comms( i_comms )%remote_rank  = remote_rank
               H%send_comms( i_comms )%comm_size    = comm_size
               H%send_comms( i_comms )%comm_start   = comm_start
               H%send_comms( i_comms )%tag          = comm_vec_to_tag( range_c, comm_vec )
               H%send_comms( i_comms )%is_local     = remote_rank == my_rank
               Allocate( H%send_comms( i_comms )%buffer( 0:H%send_comms( i_comms )%comm_size( 1 ) - 1, &
                    0:H%send_comms( i_comms )%comm_size( 2 ) - 1, &
                    0:H%send_comms( i_comms )%comm_size( 3 ) - 1 ) )

            End Do
         End Do
      End Do

      ! And finaly storgae for the mesage requests
      Allocate( H%msg_requests( 1:Size( H%recv_comms ) + Size( H%send_comms ) ) )

    End Subroutine combine_comms

    Pure Function comm_vec_to_tag( range_c, comm_vec ) Result( tag )

      Integer :: tag

      Integer,                   Intent( In ) :: range_c
      Integer, Dimension( 1:3 ), Intent( In ) :: comm_vec

      Integer, Dimension( 1:3 ) :: shifted_vec

      shifted_vec = comm_vec + range_c

      tag = 1 + shifted_vec( 1 ) + ( 2 * range_c + 1 ) * shifted_vec( 2 ) + &
           ( 2 * range_c + 1 ) * ( 2 * range_c + 1 ) * shifted_vec( 3 )

    End Function comm_vec_to_tag

    Subroutine factor_comms( comm, n_procs, is_periodic, my_coord, local_size, halo_width, total_size, first_point, &
         one_d_recv, one_d_send )

      Use mpi_f08, Only : mpi_allgather, mpi_integer

      Type( mpi_comm ),    Intent( In    ) :: comm
      Integer,             Intent( In    ) :: n_procs
      Logical,             Intent( In    ) :: is_periodic ! pass this for future support of open boundary cases
      Integer,             Intent( In    ) :: my_coord
      Integer,             Intent( In    ) :: local_size
      Integer,             Intent( In    ) :: halo_width
      Integer,             Intent(   Out ) :: total_size
      Integer,             Intent(   Out ) :: first_point
      Type( one_d_comms ), Intent(   Out ) :: one_d_recv
      Type( one_d_comms ), Intent(   Out ) :: one_d_send

      Integer, Dimension( : ), Allocatable :: axis_local_sizes
      Integer, Dimension( : ), Allocatable :: coord_to_recv_from, size_to_recv_from, recv_start
      Integer, Dimension( : ), Allocatable :: coord_to_send_to,   size_to_send_to,   send_start

      Integer :: left_first_proc, right_last_proc
      Integer :: points_remaining, size_this_proc
      Integer :: proc_periodic
      Integer :: start_comm, end_comm

      Allocate( axis_local_sizes( 0:n_procs - 1 ) )

      Call mpi_allgather( local_size, 1, mpi_integer, axis_local_sizes, 1, mpi_integer, comm )

      ! Now have the total size of the grid
      total_size = Sum( axis_local_sizes )

      ! Can work out the what my first point is (grids start at zero!)
      ! note the sum of a zero size array is zero by the standard
      first_point = Sum( axis_local_sizes( 0:my_coord - 1 ) )

      ! Now work out who we have to recv from
      Allocate( coord_to_recv_from( 0:-1 ) )
      Allocate( size_to_recv_from ( 0:-1 ) )
      Allocate( recv_start        ( 0:-1 ) )

      ! DOESN'T WORK for HALO_WIDTH > LOCAL_SIZE
      ! Need to come back and rethink

      ! First consider to the left
      points_remaining = halo_width
      left_first_proc = my_coord
      end_comm = -1
      Do While( points_remaining > 0 )
         left_first_proc = left_first_proc - 1
         ! Periodic Boundary conditions to find size, but note keep "real", unshifted coordinate
         ! of the proc for comms working - just is easier and MPI does most of the periodic stuff for us
         proc_periodic = Merge( Modulo( left_first_proc, n_procs ), left_first_proc, is_periodic )
         size_this_proc = Min( points_remaining, axis_local_sizes( proc_periodic ) )
         start_comm = end_comm - size_this_proc + 1
         coord_to_recv_from = [ coord_to_recv_from, left_first_proc ]
         size_to_recv_from  = [ size_to_recv_from,  size_this_proc  ]
         recv_start         = [ recv_start,         start_comm      ]
         points_remaining = points_remaining - size_this_proc
         end_comm = start_comm - 1
      End Do

      ! Now the middle - i.e. comms orthogonal to this axis, so just a local size length bit
      coord_to_recv_from = [ coord_to_recv_from, my_coord   ]
      size_to_recv_from  = [ size_to_recv_from,  local_size ]
      recv_start         = [ recv_start,         0          ]

      ! And now to the right
      points_remaining = halo_width
      right_last_proc = my_coord
      start_comm = local_size
      Do While( points_remaining > 0 )
         right_last_proc = right_last_proc + 1
         proc_periodic = Merge( Modulo( right_last_proc, n_procs ), right_last_proc, is_periodic )
         size_this_proc = Min( points_remaining, axis_local_sizes( proc_periodic ) )
         coord_to_recv_from = [ coord_to_recv_from, right_last_proc ]
         size_to_recv_from  = [ size_to_recv_from,  size_this_proc  ]
         recv_start         = [ recv_start,         start_comm      ]
         points_remaining = points_remaining - size_this_proc
         start_comm = start_comm + size_this_proc
      End Do

      one_d_recv%coords   = coord_to_recv_from
      one_d_recv%sizes    = size_to_recv_from
      one_d_recv%starts   = recv_start

      ! Now work out who we have to send to
      ! On second thoughts I think the send and recv will always be the same!
      ! But keep them different for the moment just in case I have missed something,
      ! and as start places different it just makes the code easier to follow
      ! for a minimal investment in memory
      Allocate( coord_to_send_to( 0:-1 ) )
      Allocate( size_to_send_to ( 0:-1 ) )
      Allocate( send_start      ( 0:-1 ) )

      ! Note we need to order the sends in the opposite order to the recvs
      ! So first to the right
      points_remaining = halo_width
      right_last_proc = my_coord
      start_comm = local_size - halo_width
      Do While( points_remaining > 0 )
         right_last_proc = right_last_proc + 1
         ! Periodic Boundary conditions
         proc_periodic = Merge( Modulo( right_last_proc, n_procs ), right_last_proc, is_periodic )
         size_this_proc = Min( points_remaining, axis_local_sizes( proc_periodic ) )
         coord_to_send_to = [ coord_to_send_to, right_last_proc ]
         size_to_send_to  = [ size_to_send_to,  size_this_proc  ]
         send_start       = [ send_start,       start_comm      ]
         points_remaining = points_remaining - size_this_proc
         start_comm = start_comm + size_this_proc
      End Do

      ! Now the middle - i.e. comms orthogonal to this axis, so just a local size length bit
      coord_to_send_to = [ coord_to_send_to, my_coord   ]
      size_to_send_to  = [ size_to_send_to,  local_size ]
      send_start       = [ send_start,       0          ]

      ! Last consider to the left
      points_remaining = halo_width
      left_first_proc = my_coord
      end_comm = halo_width - 1
      Do While( points_remaining > 0 )
         left_first_proc = left_first_proc - 1
         proc_periodic = Merge( Modulo( left_first_proc, n_procs ), left_first_proc, is_periodic )
         size_this_proc = Min( points_remaining, axis_local_sizes( proc_periodic ) )
         start_comm = end_comm - size_this_proc + 1
         coord_to_send_to = [ coord_to_send_to, left_first_proc ]
         size_to_send_to  = [ size_to_send_to,  size_this_proc  ]
         send_start       = [ send_start,         start_comm    ]
         points_remaining = points_remaining - size_this_proc
         end_comm = start_comm - 1
      End Do

      one_d_send%coords = coord_to_send_to
      one_d_send%sizes  = size_to_send_to
      one_d_send%starts = send_start

    End Subroutine factor_comms

  End Subroutine halo_parallel_init_f08

  Subroutine halo_parallel_init_old( H, local_size, halo_width, comm, error )

    Use mpi_f08, Only : mpi_comm

    ! GRIDS START AT ZERO

    Class( halo_parallel_setter_2 ),                   Intent( InOut ) :: H
    Integer,                         Dimension( 1:3 ), Intent( In    ) :: local_size
    Integer,                                           Intent( In    ) :: halo_width
    Integer,                                           Intent( In    ) :: comm
    Integer,                                           Intent(   Out ) :: error

    Type( mpi_comm ) :: comm_f08

    comm_f08%mpi_val = comm
    Call halo_parallel_init_f08( H, local_size, halo_width, comm_f08, error )

  End Subroutine halo_parallel_init_old

  Subroutine halo_fill( H, halo_width, hdlb, gin, hout, error )

    Use mpi_f08, Only : mpi_waitall, mpi_statuses_ignore, mpi_isend, mpi_irecv, &
         mpi_double_precision, mpi_request_null

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Class( halo_parallel_setter_2 ),                             Intent( InOut ) :: H
    Integer,                                                     Intent( In    ) :: halo_width
    Integer,    Dimension( 1:3 ),                                Intent( In    ) :: hdlb
    Real( wp ), Dimension( 0:, 0:, 0: ),                         Intent( In    ) :: gin
    Real( wp ), Dimension( - halo_width:, - halo_width:, - halo_width: ), Intent(   Out ) :: hout
    Integer,                                                     Intent(   Out ) :: error

    Integer, Dimension( 1:3 ) :: s, e, ss, es

    Integer :: i_comms

    error = 0

    If( halo_width /= H%halo_width ) Then
       error = 10
       Return
    End If

    Call H%inc_n_calls()

    ! Loops over the recvs
    Do i_comms = 1, Size( H%recv_comms )
       ! Post the recv
       ! FIX MPI_DOUBLE_PRECISON to make portable
       If( .Not. H%recv_comms( i_comms )%is_local ) Then
          Call mpi_irecv( H%recv_comms( i_comms )%buffer, Size( H%recv_comms( i_comms )%buffer ), mpi_double_precision, &
               H%recv_comms( i_comms )%remote_rank, H%recv_comms( i_comms )%tag, H%comm, &
               H%msg_requests( i_comms ) )
       Else
          H%msg_requests( i_comms ) = mpi_request_null
       End If
    End Do

    ! Loops over the sends
    Do i_comms = 1, Size( H%recv_comms )
       ! Copy the relevant part into the buffer
       If( .Not. H%send_comms( i_comms )%is_local ) Then
          s = H%send_comms( i_comms )%comm_start
          e = H%send_comms( i_comms )%comm_start + H%send_comms( i_comms )%comm_size - 1
          H%send_comms( i_comms )%buffer = gin( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) )
          ! Post the send
          Call mpi_isend( H%send_comms( i_comms )%buffer, Size( H%send_comms( i_comms )%buffer ), mpi_double_precision, &
               H%send_comms( i_comms )%remote_rank, H%send_comms( i_comms )%tag, H%comm, &
               H%msg_requests( i_comms + Size( H%recv_comms ) ) )
       Else
          H%msg_requests( i_comms + Size( H%recv_comms )  ) = mpi_request_null
       End If
    End Do

    ! Wait on the async comms
    Call mpi_waitall( Size( H%msg_requests ), H%msg_requests, mpi_statuses_ignore )

    Do i_comms = 1, Size( H%recv_comms )
       ! Copy the buffer into the relevant parts
       s = H%recv_comms( i_comms )%comm_start
       e = H%recv_comms( i_comms )%comm_start + H%recv_comms( i_comms )%comm_size - 1
       If( .Not. H%recv_comms( i_comms )%is_local ) Then
          hout( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) ) = H%recv_comms( i_comms )%buffer
       Else
          ss = H%send_comms( i_comms )%comm_start
          es = H%send_comms( i_comms )%comm_start + H%send_comms( i_comms )%comm_size - 1
          hout( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) ) = &
               gin( ss( 1 ):es( 1 ), ss( 2 ):es( 2 ), ss( 3 ):es( 3 ) )
       End If
    End Do

  End Subroutine halo_fill

End Module halo_parallel_module
