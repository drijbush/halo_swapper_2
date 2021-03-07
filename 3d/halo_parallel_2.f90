Module halo_parallel_module

  ! ToDo:
  ! 1) Rationalise error codes
  ! 2) Add flags for orthogonal cells which avoid unneccesary comms

  Use mpi_f08, Only : mpi_comm, mpi_request

  Use constants  , Only : wp
  Use swap_module, Only : halo_dim_plan_type

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
       Call mpi_comm_free( axis_comm, error )
    End Do

  Contains

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

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Class( halo_parallel_setter_2 ),                             Intent( InOut ) :: H
    Integer,                                                     Intent( In    ) :: halo_width
    Integer,    Dimension( 1:3 ),                                Intent( In    ) :: hdlb
    Real( wp ), Dimension( 0:, 0:, 0: ),                         Intent( In    ) :: gin
    Real( wp ), Dimension( - halo_width:, - halo_width:, - halo_width: ), Intent(   Out ) :: hout
    Integer,                                                     Intent(   Out ) :: error

    ! Dummy code to kill warnings while developing elsewhere
    error = 0
    hout = gin
    Write( *, * ) H%halo_width, halo_width, hdlb
    
!!$    Integer, Dimension( 1:3 ) :: s, e, ss, es
!!$
!!$    Integer :: i_comms
!!$
!!$    error = 0
!!$
!!$    If( halo_width /= H%halo_width ) Then
!!$       error = 10
!!$       Return
!!$    End If
!!$
!!$    Call H%inc_n_calls()
!!$
!!$    ! Loops over the recvs
!!$    Do i_comms = 1, Size( H%recv_comms )
!!$       ! Post the recv
!!$       ! FIX MPI_DOUBLE_PRECISON to make portable
!!$       If( .Not. H%recv_comms( i_comms )%is_local ) Then
!!$          Call mpi_irecv( H%recv_comms( i_comms )%buffer, Size( H%recv_comms( i_comms )%buffer ), mpi_double_precision, &
!!$               H%recv_comms( i_comms )%remote_rank, H%recv_comms( i_comms )%tag, H%comm, &
!!$               H%msg_requests( i_comms ) )
!!$       Else
!!$          H%msg_requests( i_comms ) = mpi_request_null
!!$       End If
!!$    End Do
!!$
!!$    ! Loops over the sends
!!$    Do i_comms = 1, Size( H%recv_comms )
!!$       ! Copy the relevant part into the buffer
!!$       If( .Not. H%send_comms( i_comms )%is_local ) Then
!!$          s = H%send_comms( i_comms )%comm_start
!!$          e = H%send_comms( i_comms )%comm_start + H%send_comms( i_comms )%comm_size - 1
!!$          H%send_comms( i_comms )%buffer = gin( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) )
!!$          ! Post the send
!!$          Call mpi_isend( H%send_comms( i_comms )%buffer, Size( H%send_comms( i_comms )%buffer ), mpi_double_precision, &
!!$               H%send_comms( i_comms )%remote_rank, H%send_comms( i_comms )%tag, H%comm, &
!!$               H%msg_requests( i_comms + Size( H%recv_comms ) ) )
!!$       Else
!!$          H%msg_requests( i_comms + Size( H%recv_comms )  ) = mpi_request_null
!!$       End If
!!$    End Do
!!$
!!$    ! Wait on the async comms
!!$    Call mpi_waitall( Size( H%msg_requests ), H%msg_requests, mpi_statuses_ignore )
!!$
!!$    Do i_comms = 1, Size( H%recv_comms )
!!$       ! Copy the buffer into the relevant parts
!!$       s = H%recv_comms( i_comms )%comm_start
!!$       e = H%recv_comms( i_comms )%comm_start + H%recv_comms( i_comms )%comm_size - 1
!!$       If( .Not. H%recv_comms( i_comms )%is_local ) Then
!!$          hout( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) ) = H%recv_comms( i_comms )%buffer
!!$       Else
!!$          ss = H%send_comms( i_comms )%comm_start
!!$          es = H%send_comms( i_comms )%comm_start + H%send_comms( i_comms )%comm_size - 1
!!$          hout( s( 1 ):e( 1 ), s( 2 ):e( 2 ), s( 3 ):e( 3 ) ) = &
!!$               gin( ss( 1 ):es( 1 ), ss( 2 ):es( 2 ), ss( 3 ):es( 3 ) )
!!$       End If
!!$    End Do

  End Subroutine halo_fill

End Module halo_parallel_module
