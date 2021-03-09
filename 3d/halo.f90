Program halo3

  Use constants, Only : wp
  
  Use, Intrinsic :: iso_fortran_env, Only : output_unit

  Use mpi_f08, Only : mpi_comm, mpi_comm_world, mpi_init, mpi_comm_size, mpi_comm_rank, mpi_finalize, &
       mpi_allreduce, mpi_in_place, mpi_integer, mpi_sum, mpi_bcast, mpi_barrier, mpi_cart_create, &
       mpi_cart_coords, mpi_dims_create, mpi_cart_sub, mpi_comm_free, mpi_cart_get

  Use swap_module, Only : halo_dim_plan_type

  Use halo_parallel_module, Only : halo_parallel_setter_2
  
  Implicit None

  Integer, Parameter :: base = 10

  Type dims_data
     Integer, Dimension( : ), Allocatable :: n
  End Type dims_data

  Type( dims_data ), Dimension( 1:3 ) :: n_data_all_3d
  Type( dims_data ), Dimension( 1:3 ) :: i_start_3d, i_end_3d
  
  Type( halo_parallel_setter_2 ) :: H

  Type( mpi_comm ) :: cart_comm
  Type( mpi_comm ) :: plane_comm
  Type( mpi_comm ) :: axis_comm
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: data_3d
  Real( wp ), Dimension( :, :, : ), Allocatable :: data_3d_with_halo
  
  Real :: rtmp

  Integer, Dimension( : ), Allocatable :: data
  Integer, Dimension( : ), Allocatable :: n_data_all
  Integer, Dimension( : ), Allocatable :: i_start, i_end

  Integer, Dimension( 1:3 ) :: n_data_3d
  Integer, Dimension( 1:3 ) :: n_3d
  Integer, Dimension( 1:3 ) :: np_grid, p_coords, n_coords, g_ranks

  Logical, Dimension( 1:3 ) :: is_this_axis, is_this_orthog_plane, is_periodic
  
  Integer :: n
  Integer :: n_data
  Integer :: rank, nprc, me_cart, me_plane
  Integer :: n_halo
  Integer :: out
  Integer :: check_val
  Integer :: sx, sy, sz
  Integer :: ex, ey, ez
  Integer :: i, j
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
  Call mpi_cart_get( cart_comm, 3, n_coords, is_periodic, p_coords, error )
  Write( output_unit, '( "grid dat ", 7( i2, 1x ), 3( l1, 1x ) )' ) rank, n_coords, p_coords, is_periodic
  Flush( output_unit )

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
        Call Random_number( rtmp )
        n_data_3d( i ) = 1 + Int( 9.0 * rtmp )
     End If
     Call mpi_bcast( n_data_3d( i ), 1, mpi_integer, 0, plane_comm, error )
     Call mpi_comm_free( plane_comm, error )
  End Do
  Do i = 1, 3
     is_this_axis = .False.
     is_this_axis( i ) = .True.
     ! Take the plane orthogonal to it
     ! Create a communicator containg this process in that plane
     Call mpi_cart_sub( cart_comm, is_this_axis, axis_comm )
     Call mpi_comm_size( axis_comm, n_coords( i ), error )
     Call mpi_comm_rank( axis_comm, g_ranks( i ), error )
     Allocate( n_data_all_3d( i )%n( 0:n_coords( i ) - 1 ) )
     n_data_all_3d( i )%n = 0
     n_data_all_3d( i )%n( g_ranks( i ) ) = n_data_3d( i )
     Call mpi_allreduce( mpi_in_place, n_data_all_3d( i )%n, Size( n_data_all_3d( i )%n ), &
          mpi_integer, mpi_sum, axis_comm, error )
     Call mpi_comm_free( axis_comm, error )
     n_3d( i ) = Sum( n_data_all_3d( i )%n )
  End Do


  If( rank == 0 ) Then
     Write( *, * ) '!!!!!!!!!!'
     Write( *, * ) 
     Write( *, * )  '3D'
     Write( output_unit, * ) 'n =      ', n_3d
     Write( output_unit, * ) 'n_halo = ', n_halo
     Do i = 1, 3
        Write( output_unit, * ) 'dir ', i, ' n proc = ', n_data_all_3d( i )%n
     End Do
     Flush( output_unit )
  End If
  Call mpi_barrier( mpi_comm_world, error )
  
  Do i = 1, 3
     Allocate( i_start_3d( i )%n( 0:n_coords( i ) - 1 ) )
     Allocate( i_end_3d( i )%n( 0:n_coords( i ) - 1 ) )
     Do j = 0,  n_coords( i ) - 1
        i_start_3d( i )%n( j ) = Sum( n_data_all_3d( i )%n( 0:j - 1 ) )
     End Do
     i_end_3d( i )%n = i_start_3d( i )%n + n_data_all_3d( i )%n - 1
  End Do

  sx = i_start_3d( 1 )%n( g_ranks( 1 ) )
  sy = i_start_3d( 2 )%n( g_ranks( 2 ) )
  sz = i_start_3d( 3 )%n( g_ranks( 3 ) )
  
  ex = i_end_3d( 1 )%n( g_ranks( 1 ) )
  ey = i_end_3d( 2 )%n( g_ranks( 2 ) )
  ez = i_end_3d( 3 )%n( g_ranks( 3 ) )
  
  Call H%init( n_data_3d, n_halo, cart_comm, error )
  If( error /= 0 ) Then
     Write( *, * ) 'Error = ', error, rank
     Call mpi_finalize( error )
     Stop
  End If

  Allocate( data_3d( sx:ex, sy:ey, sz:ez ) )
  
  Do iz = sz, ez
     Do iy = sy, ey
        Do ix = sx, ex
           data_3d( ix, iy, iz ) = ix + base * iy + base * base * iz
        End Do
     End Do
  End Do

  Allocate( data_3d_with_halo( sx - n_halo:ex + n_halo, sy:ey, sz:ez ) )
  
  Call H%fill( n_halo, n_data_3d, data_3d, data_3d_with_halo, error )
  
  Call mpi_barrier( mpi_comm_world, error )
  worked_size = Ubound( data_3d_with_halo, Dim = 1 ) == ex + n_halo
  worked_size = worked_size .And. Lbound( data_3d_with_halo, Dim = 1 ) == sx - n_halo
!!$  Write( out, * )
!!$  Write( out, * ) np_grid, g_ranks, p_coords
!!$  Write( out, * ) n_halo, n_3d
!!$  Write( out, * ) sx, sy, sz, ex, ey, ez
!!$  Write( out, * ) Lbound( data_3d, Dim = 1 ), Ubound( data_3d, Dim = 1 )
!!$  Write( out, * ) Lbound( data_3d_with_halo, Dim = 1 ), Ubound( data_3d_with_halo, Dim = 1 )
!!$  Write( out, * )
  worked_data = .True.
  Do iz = Lbound( data_3d_with_halo, Dim = 3 ), Ubound( data_3d_with_halo, Dim = 3 )
     Do iy = Lbound( data_3d_with_halo, Dim = 2 ), Ubound( data_3d_with_halo, Dim = 2 )
        Do ix = Lbound( data_3d_with_halo, Dim = 1 ), Ubound( data_3d_with_halo, Dim = 1 )
           check_val =                           Modulo( ix, n_3d( 1 ) )
           check_val = check_val + base *        Modulo( iy, n_3d( 2 ) )
           check_val = check_val + base * base * Modulo( iz, n_3d( 3 ) )
           worked_data = worked_data .And. Nint( data_3d_with_halo( ix, iy, iz ) ) == check_val
           If( Nint( data_3d_with_halo( ix, iy, iz ) ) /= check_val ) Then
              Write( out, * ) ix, iy, iz, check_val, data_3d_with_halo( ix, iy, iz )
           End If
        End Do
     End Do
  End Do
!!$  Write( out, * )
!!$  Do iz = Lbound( data_3d, Dim = 3 ), Ubound( data_3d, Dim = 3 )
!!$     Write( out, * ) iz
!!$     Write( out, '( 2( 2( i3, 1x ) / ) )' ) Nint( data_3d( :, :, iz ) )
!!$  End Do
!!$  Write( out, * ) 'with halo'
!!$  Do iz = Lbound( data_3d_with_halo, Dim = 3 ), Ubound( data_3d_with_halo, Dim = 3 )
!!$     Write( out, * ) iz
!!$     Write( out, '( 2( 4( i3, 1x ) / ) )' ) Nint( data_3d_with_halo( :, :, iz ) )
!!$  End Do
  Write( output_unit, * ) 'Checking rank ', rank, worked_size, worked_data
  If( ( .Not. worked_size ) .Or. ( .Not. worked_data ) ) Write( output_unit, * ) rank, ' BUSTED', n_halo
  
  Call mpi_finalize( error )

Contains

End Program halo3
