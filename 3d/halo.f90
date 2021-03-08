Program halo3

  Use constants, Only : wp
  
  Use, Intrinsic :: iso_fortran_env, Only : output_unit

  Use mpi_f08, Only : mpi_comm, mpi_comm_world, mpi_init, mpi_comm_size, mpi_comm_rank, mpi_finalize, &
       mpi_allreduce, mpi_in_place, mpi_integer, mpi_sum, mpi_bcast, mpi_barrier, mpi_cart_create, &
       mpi_cart_coords, mpi_dims_create, mpi_cart_sub, mpi_comm_free

  Use swap_module, Only : halo_dim_plan_type

  Use halo_parallel_module, Only : halo_parallel_setter_2
  
  Implicit None

  Type( halo_parallel_setter_2 ) :: H

  Type( mpi_comm ) :: cart_comm
  Type( mpi_comm ) :: plane_comm
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: data_3d
  Real( wp ), Dimension( :, :, : ), Allocatable :: data_3d_with_halo
  
  Real, Dimension( 1:3 ) :: rtmp_3d

  Real :: rtmp

  Integer, Dimension( :, : ), Allocatable :: n_data_all_3d
  Integer, Dimension( :, : ), Allocatable :: i_start_3d, i_end_3d

  Integer, Dimension( : ), Allocatable :: data
  Integer, Dimension( : ), Allocatable :: n_data_all
  Integer, Dimension( : ), Allocatable :: i_start, i_end

  Integer, Dimension( 1:3 ) :: n_data_3d
  Integer, Dimension( 1:3 ) :: n_3d
  Integer, Dimension( 1:3 ) :: np_grid, p_coords

  Logical, Dimension( 1:3 ) :: is_this_axis, is_this_orthog_plane
  
  Integer :: n
  Integer :: n_data
  Integer :: rank, nprc, me_cart, me_plane
  Integer :: n_halo
  Integer :: out
  Integer :: check_val
  Integer :: i
  Integer :: ix, iy, iz

  Integer :: error

  Logical :: worked_size, worked_data

  Type( halo_dim_plan_type ) :: dim_plan
  
  Call mpi_init( error )
  Call mpi_comm_size( mpi_comm_world, nprc, error )
  Call mpi_comm_rank( mpi_comm_world, rank, error )

  out = 10 + rank  

  If( rank == 0 ) Then
     Call Random_number( rtmp )
     n_halo = 1 + Int( 20.0 * rtmp )
  End If
  Call mpi_bcast( n_halo, 1, mpi_integer, 0, mpi_comm_world, error )

  Call Random_number( rtmp )
  n_data = 1 + Int( 15.0 * rtmp )

  Allocate( n_data_all( 0:nprc - 1 ) )
  n_data_all = 0
  n_data_all( rank ) = n_data

  Call mpi_allreduce( mpi_in_place, n_data_all, Size( n_data_all ), mpi_integer, mpi_sum, mpi_comm_world, error )
  n = Sum( n_data_all )

  If( rank == 0 ) Then
     Write( output_unit, * ) 'n =      ', n
     Write( output_unit, * ) 'n_halo = ', n_halo
     Write( output_unit, * ) 'n proc = ', n_data_all
     Flush( output_unit )
  End If
  Call mpi_barrier( mpi_comm_world, error )

  Allocate( i_start( 0:nprc - 1 ) )
  Allocate( i_end  ( 0:nprc - 1 ) )
  Do i = 0, nprc - 1
     i_start( i ) = Sum( n_data_all( 0:i - 1 ) )
  End Do
  i_end = i_start + n_data_all - 1

  Write( out, * ) 'n, n_halo ', n, n_halo
  Write( out, * ) 'n proc ', n_data_all

  Call mpi_barrier( mpi_comm_world, error )
  Write( out, * ) '!!!!!!!!!!!!!!!!!'
  If( rank == 0 ) Write( *, * ) '!!!!!!!!!!!!!!!!!'
  Call mpi_barrier( mpi_comm_world, error )
  Call dim_plan%init( n_data, n_halo, mpi_comm_world, error )
  Call dim_plan%report( out )
  Allocate( data( i_start( rank ):i_end( rank ) ) )
  Do i = Lbound( data, Dim = 1 ), Ubound( data, Dim = 1 )
     data( i ) = i
  End Do
  Call dim_plan%fill( data, .True. )
  ! Finally check everything has worked
  ! First check bounds
  worked_size = Ubound( data, Dim = 1 ) == i_end( rank ) + n_halo
  worked_size = worked_size .And. Lbound( data, Dim = 1 ) == i_start( rank ) - n_halo
  ! Now the data
  worked_data = .True.
  Do i = Ubound( data, Dim = 1 ), Lbound( data, Dim = 1 ), -1
     worked_data = worked_data .And. Modulo( i, n ) == data( i )
  End Do
  Write( output_unit, * ) 'Checking rank ', rank, worked_size, worked_data
  If( ( .Not. worked_size ) .Or. ( .Not. worked_data ) ) Write( *, * ) rank, ' BUSTED', n_halo, n_data_all
  Flush( output_unit )

  ! END OF 1D

  ! Now 3d
  np_grid = 0
  Call mpi_dims_create( nprc, 3, np_grid )
  Call mpi_cart_create( mpi_comm_world, 3, np_grid, [ .True., .True., .True. ], .True., &
       cart_comm )
  ! Get the coordinates for this processor in the process grid
  Call mpi_comm_rank( cart_comm, me_cart,    error )
  Call mpi_cart_coords( cart_comm, me_cart, 3, p_coords )

  ! Generate grid dims
  Do i = 1, 3
     is_this_axis = .False.
     is_this_axis( i ) = .True.
     ! Take the plane orthogonal to it
     is_this_orthog_plane = .Not. is_this_axis
     ! Create a communicator containg this process in that plane
     Call mpi_cart_sub( cart_comm, is_this_orthog_plane, plane_comm )
     ! Zero proc in plane decides dimension
     Call mpi_comm_rank( plane_comm, me_plane, error )
     If( me_plane == 0 ) Then
        Call Random_number( rtmp_3d )
        n_data_3d = 1 + Int( 9.0 * rtmp_3d )
     End If
     Call mpi_bcast( n_data_3d, Size( n_data_3d ), mpi_integer, 0, plane_comm, error )
     Call mpi_comm_free( plane_comm, error )
  End Do

  Allocate( n_data_all_3d( 1:3, 0:nprc - 1 ) )
  n_data_all_3d = 0
  n_data_all_3d( :, rank ) = n_data_3d

  Call mpi_allreduce( mpi_in_place, n_data_all_3d, Size( n_data_all_3d ), mpi_integer, mpi_sum, mpi_comm_world, error )
  n_3d = Sum( n_data_all_3d, Dim = 2 )

  If( rank == 0 ) Then
     Write( *, * ) '!!!!!!!!!!'
     Write( *, * ) 
     Write( *, * )  '3D'
     Write( output_unit, * ) 'n =      ', n_3d
     Write( output_unit, * ) 'n_halo = ', n_halo
     Write( output_unit, * ) 'n proc = ', n_data_all_3d
     Flush( output_unit )
  End If
  Call mpi_barrier( mpi_comm_world, error )
  
  Allocate( i_start_3d( 1:3, 0:nprc - 1 ) )
  Allocate( i_end_3d  ( 1:3, 0:nprc - 1 ) )
  Do i = 0, nprc - 1
     i_start_3d( :, i ) = Sum( n_data_all_3d( :, 0:i - 1 ), Dim = 2 )
  End Do
  i_end_3d = i_start_3d + n_data_all_3d - 1

  Call H%init( n_data_3d, n_halo, cart_comm, error )
  If( error /= 0 ) Then
     Write( *, * ) 'Error = ', error, rank
     Call mpi_finalize( error )
     Stop
  End If

  Allocate( data_3d( i_start_3d( 1, rank ):i_end_3d( 1, rank ), &
       i_start_3d( 2, rank ):i_end_3d( 2, rank ), &
       i_start_3d( 3, rank ):i_end_3d( 3, rank ) ) )
  
  Do iz = i_start_3d( 3, rank ), i_end_3d( 3, rank )
     Do iy = i_start_3d( 2, rank ), i_end_3d( 2, rank )
        Do ix = i_start_3d( 1, rank ), i_end_3d( 1, rank )
           data_3d( ix, iy, iz ) = ix + 10 * iy + 100 * iz
        End Do
     End Do
  End Do

!!$  Allocate( data_3d_with_halo, Source = data_3d )
  Allocate( data_3d_with_halo( i_start_3d( 1, rank ) - n_halo:i_end_3d( 1, rank ) + n_halo, &
       i_start_3d( 2, rank ):i_end_3d( 2, rank ), &
       i_start_3d( 3, rank ):i_end_3d( 3, rank ) ) )
  Call H%fill( n_halo, n_data_3d, data_3d, data_3d_with_halo, error )
  
  Call mpi_barrier( mpi_comm_world, error )
  worked_size = Ubound( data_3d_with_halo, Dim = 1 ) == i_end_3d( 1, rank ) + n_halo
  worked_size = worked_size .And. Lbound( data_3d_with_halo, Dim = 1 ) == i_start_3d( 1, rank ) - n_halo
  worked_data = .True.
  Do iz = i_start_3d( 3, rank ), i_end_3d( 3, rank )
     Do iy = i_start_3d( 2, rank ), i_end_3d( 2, rank )
        Do ix = i_start_3d( 1, rank ), i_end_3d( 1, rank )
           check_val = Modulo( ix, n_3d( 1 ) )
           check_val = check_val + 10  * Modulo( iy, n_3d( 2 ) )
           check_val = check_val + 100 * Modulo( iz, n_3d( 3 ) )
           worked_data = worked_data .And. Nint( data_3d_with_halo( ix, iy, iz ) ) == check_val
        End Do
     End Do
  End Do
  Write( output_unit, * ) 'Checking rank ', rank, worked_size, worked_data
  If( ( .Not. worked_size ) .Or. ( .Not. worked_data ) ) Write( output_unit, * ) rank, ' BUSTED', n_halo, n_data_all_3d
  
  Call mpi_finalize( error )

Contains

End Program halo3
