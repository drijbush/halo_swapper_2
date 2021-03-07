Program halo3

  Use, Intrinsic :: iso_fortran_env, Only : output_unit

  Use mpi_f08, Only : mpi_comm_world, mpi_init, mpi_comm_size, mpi_comm_rank, mpi_finalize, &
       mpi_allreduce, mpi_in_place, mpi_integer, mpi_sum, mpi_bcast, mpi_barrier

  Use swap_module, Only : halo_dim_plan_type
  
  Implicit None

  Real :: rtmp

  Integer, Dimension( : ), Allocatable :: data
  Integer, Dimension( : ), Allocatable :: n_data_all
  Integer, Dimension( : ), Allocatable :: i_start, i_end

  Integer :: n
  Integer :: n_data
  Integer :: rank, nprc
  Integer :: prev, next
  Integer :: n_halo
  Integer :: out
  Integer :: i

  Integer :: error

  Logical :: worked_size, worked_data

  Type( halo_dim_plan_type ) :: dim_plan
  
  Call mpi_init( error )
  Call mpi_comm_size( mpi_comm_world, nprc, error )
  Call mpi_comm_rank( mpi_comm_world, rank, error )

  prev = Modulo( rank - 1, nprc )
  next = Modulo( rank + 1, nprc )

  out = 10 + rank  

  If( rank == 0 ) Then
     Call Random_number( rtmp )
     n_halo = 1 + Int( 20.0 * rtmp )
  End If
  Call mpi_bcast( n_halo, 1, mpi_integer, 0, mpi_comm_world, error )

  ! HACK
!!$  n_halo = 3

  Call Random_number( rtmp )
  n_data = 1 + Int( 15.0 * rtmp )

  Allocate( n_data_all( 0:nprc - 1 ) )
  n_data_all = 0
  n_data_all( rank ) = n_data

  !HACK
!!$  n_data_all = 2
  
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
  Allocate( data( i_start( rank ):i_end( rank ) ) )
  Do i = Lbound( data, Dim = 1 ), Ubound( data, Dim = 1 )
     data( i ) = i
  End Do
  Call dim_plan%fill( data, .True. )
  ! Finally check everything has worked
  ! First check bounds
  worked_size = Ubound( data, Dim = 1 ) == i_end( rank ) + n_halo
  ! Now the data
  worked_data = .True.
  Do i = Ubound( data, Dim = 1 ), Lbound( data, Dim = 1 ), -1
     worked_data = worked_data .And. Modulo( i, n ) == data( i )
  End Do
  Write( output_unit, * ) 'Checking rank ', rank, worked_size, worked_data
  If( ( .Not. worked_size ) .Or. ( .Not. worked_data ) ) Write( *, * ) rank, ' BUSTED', n_halo, n_data_all

  Call mpi_finalize( error )

Contains

End Program halo3
