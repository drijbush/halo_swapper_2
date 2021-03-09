Module halo_parallel_module

  ! ToDo:
  ! 1) Rationalise error codes
  ! 2) Add flags for orthogonal cells which avoid unneccesary comms

  Use swap_module, Only : halo_dim_plan_type, FILL_X, FILL_Y, FILL_Z

  Use mpi_f08, Only : mpi_comm, mpi_request

  Use constants  , Only : wp

  Use halo_setter_base_module, Only : halo_setter_base_class

  Implicit None

  Type, Public, Extends( halo_setter_base_class ) :: halo_parallel_setter_2
     Private
     Type( mpi_comm ),                             Private :: comm
     Integer,                                      Private :: halo_width
     Integer,                    Dimension( 1:3 ), Private :: local_size
     Integer,                    Dimension( 1:3 ), Private :: total_size
     Integer,                    Dimension( 1:3 ), Private :: first_point
     Integer,                    Dimension( 1:3 ), Private :: last_point
     Integer,                    Dimension( 1:3 ), Private :: n_procs
     Integer,                    Dimension( 1:3 ), Private :: my_coords
     Logical,                    Dimension( 1:3 ), Private :: is_periodic
     Logical                                     , Private :: corners = .True.
     Type( halo_dim_plan_type ), Dimension( 1:3 ), Private :: dim_plans
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
         mpi_integer, mpi_allreduce, mpi_max, mpi_min, mpi_in_place, mpi_comm_free

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

    ! Check it is periodic all in directions, and get other useful data,
    ! size of proc grid and where I am in the proc grid
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
       Call mpi_comm_free( plane_comm, error )
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

    ! In x, y, z directions work out the communication plans
    Do i = 1, 3
       ! Create a comunicator for the direction in question
       is_this_axis = .False.
       is_this_axis( i ) = .True.
       Call mpi_cart_sub( comm, is_this_axis, axis_comm )
       ! Now get the plans
       Call H%dim_plans( i )%init( local_size( i ), halo_width, axis_comm, error )
       ! And inquire of the plans some useful data
       Call H%dim_plans( i )%inquire( i_start = H%first_point( i ), i_end = H%last_point( i ), n_tot = H%total_size( i ) )
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

    Use constants, Only : wp

    Use swap_module, Only : halo_dim_plan_type

    Class( halo_parallel_setter_2 ),                             Intent( InOut ) :: H
    Integer,                                                     Intent( In    ) :: halo_width
    Integer,    Dimension( 1:3 ),                                Intent( In    ) :: hdlb
    Real( wp ), Dimension( 0:, 0:, 0: ),                         Intent( In    ) :: gin
    Real( wp ), Dimension( - halo_width:, - halo_width:, - halo_width: ), Intent(   Out ) :: hout
    Integer,                                                     Intent(   Out ) :: error

    ! Subsequent calling tree allocates to right size unlike earlier efforts
    ! Need to rationalise this but let's get message passing right first
    Real( wp ), Dimension( :, :, : ), Allocatable :: temp1
    Real( wp ), Dimension( :, :, : ), Allocatable :: temp2

    error = 0

    ! Just x sorted so far - let's get that going then look at subsequent
    Call H%dim_plans( 1 )%fill( FILL_X, H%corners, gin  , temp1 )
    Call H%dim_plans( 2 )%fill( FILL_Y, H%corners, temp1, temp2 )
    Call H%dim_plans( 3 )%fill( FILL_Z, H%corners, temp2, temp1  )
    Hout = temp1
    
  End Subroutine halo_fill

End Module halo_parallel_module
