Module swap_module

  Use mpi_f08, Only : mpi_comm, mpi_datatype
  
  Implicit None

  Integer, Parameter, Public :: FILL_X = 1
  Integer, Parameter, Public :: FILL_Y = 2
  Integer, Parameter, Public :: FILL_Z = 3
  
  Type, Private :: halo_step_type
     Private
     Integer                   :: n_want, n_wanted
     Integer, Dimension( 1:2 ) :: can_give
     Integer, Dimension( 1:2 ) :: got
  End Type halo_step_type

  Type, Private :: halo_plan_type
     Private
     Type( mpi_comm )                                    :: comm
     Integer                                             :: rank
     Integer                                             :: prev, next
     Integer                                             :: direction
     Integer                                             :: i_start, i_end
     Integer                                             :: n_tot
     Integer                                             :: n_halo
     Type( halo_step_type ), Dimension( : ), Allocatable :: steps
     Type( mpi_Datatype   )                              :: real_handle
   Contains
     Procedure :: plan_left       => plan_halo_swap_left
     Procedure :: plan_right      => plan_halo_swap_right
     Procedure :: report          => report_plan
     Procedure :: swap_left_1d    => swap_int_1d_left
     Procedure :: swap_right_1d   => swap_int_1d_right
     Procedure :: swap_left_3d_x  => swap_real_3d_left_x
     Procedure :: swap_right_3d_x => swap_real_3d_right_x
     Procedure :: swap_left_3d_y  => swap_real_3d_left_y
     Procedure :: swap_right_3d_y => swap_real_3d_right_y
     Procedure :: swap_left_3d_z  => swap_real_3d_left_z
     Procedure :: swap_right_3d_z => swap_real_3d_right_z
     Procedure :: swap_1d         => swap_int_1d
     Procedure :: swap_3d_x       => swap_real_3d_x
     Procedure :: swap_3d_y       => swap_real_3d_y
     Procedure :: swap_3d_z       => swap_real_3d_z
  End Type halo_plan_type

  Type, Public :: halo_dim_plan_type
     Private
     Type( halo_plan_type ) :: left
     Type( halo_plan_type ) :: right
   Contains
     Procedure, Public :: init    => halo_dim_plan_init
     Procedure, Public :: inquire => halo_dim_plan_inquire
     Procedure, Public :: fill_1d => halo_dim_plan_fill_1d
     Procedure, Public :: fill_3d => halo_dim_plan_fill_3d
     Procedure, Public :: report  => halo_dim_plan_report
     Procedure, Public :: free    => halo_dim_plan_free
     Generic :: fill => fill_1d
     Generic :: fill => fill_3d
     Final :: halo_dim_plan_final
  End Type halo_dim_plan_type

  Private

  Integer, Parameter :: LEFT  = -1
  Integer, Parameter :: RIGHT = +1

  Integer, Parameter :: PLAN_WANT_TAG = 10
  Integer, Parameter :: PLAN_GIVE_TAG = 20
  Integer, Parameter :: SWAP_1D_LEFT  = 30
  Integer, Parameter :: SWAP_1D_RIGHT = 40
  Integer, Parameter :: SWAP_3D_LEFT  = 50
  Integer, Parameter :: SWAP_3D_RIGHT = 60
  
  
Contains

  Subroutine halo_dim_plan_inquire( plan, i_start, i_end, n_tot )
    ! Inquire properties of a plan - NEEDS COMPLETING

    Implicit None

    Class( halo_dim_plan_type ), Intent( InOut )         :: plan
    Integer                    , Intent( Out ), Optional :: i_start
    Integer                    , Intent( Out ), Optional :: i_end
    Integer                    , Intent( Out ), Optional :: n_tot

    If( Present( i_start ) ) Then
       i_start = plan%left%i_start
    End If
    
    If( Present( i_end ) ) Then
       i_end = plan%left%i_end
    End If
    
    If( Present( n_tot ) ) Then
       n_tot = plan%left%n_tot
    End If
    
  End Subroutine halo_dim_plan_inquire

  Subroutine halo_dim_plan_init( plan, local_size, halo_width, comm, error )

    Use mpi_f08, Only : mpi_comm, mpi_allreduce, mpi_scan, mpi_integer, mpi_sum, mpi_type_create_f90_real

    Use constants, Only : wp
    
    Implicit None

    Class( halo_dim_plan_type ), Intent( InOut ) :: plan
    Integer,                     Intent( In    ) :: local_size
    Integer,                     Intent( In    ) :: halo_width
    Type( mpi_comm ),            Intent( In    ) :: comm
    Integer,                     Intent(   Out ) :: error

    Integer :: i_start, i_end
    Integer :: n_tot

    error = 0

    ! Generate the global size and local start and end points
    Call mpi_allreduce( local_size, n_tot, 1, mpi_integer, mpi_sum, comm, error )
    Call mpi_scan( local_size, i_end, 1, mpi_integer, mpi_sum, comm, error )
    i_end = i_end - 1
    i_start = i_end - local_size + 1

    ! Plans for comms in each direction
    Call plan%left%plan_left  ( comm, n_tot, i_start, i_end, halo_width )
    Call plan%right%plan_right( comm, n_tot, i_start, i_end, halo_width )

    ! And set up a handle for a real type to use in the communication
    Call mpi_type_create_f90_real( Precision( 1.0_wp ), Range( 1.0_wp ), plan%left%real_handle, error )
    plan%right%real_handle = plan%left%real_handle
    
  End Subroutine halo_dim_plan_init

  Subroutine halo_dim_plan_fill_1d( plan, data, report )

    Class( halo_dim_plan_type ),                 Intent( In    )              :: plan
    Integer                    , Dimension( : ), Intent( InOut ), Allocatable :: data
    Logical                    , Optional      , Intent( In    )              :: report

    Call plan%left%swap_left_1d ( data, report )
    Call plan%right%swap_right_1d( data, report )

  End Subroutine halo_dim_plan_fill_1d

  Subroutine halo_dim_plan_fill_3d( plan, which, lbd, data, data_with_halo )

    Use constants, Only : wp
    
    Class( halo_dim_plan_type ),                                               Intent( In    )              :: plan
    Integer                                                                  , Intent( In    )              :: which
    Integer                    , Dimension( 1:3                             ), Intent( In    )              :: lbd
    Real( wp )                 , Dimension( lbd( 1 ):, lbd( 2 ):, lbd( 3 ): ), Intent( In    )              :: data
    Real( wp )                 , Dimension(         :,         :,         : ), Intent(   Out ), Allocatable :: data_with_halo

    Real( wp ), Dimension( :, :, : ), Allocatable :: temp

    Integer, Dimension( 1:3 ) :: lbd_h, ubd_h
    
    lbd_h = lbd
    ubd_h = Ubound( data )

    Select Case( which )

    Case( FILL_X )
       lbd_h( 1 ) = lbd_h( 1 ) - plan%left%n_halo
       Allocate( temp( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%left%swap_3d_x( lbd, data, lbd_h, temp )
       
       ubd_h( 1 ) = ubd_h( 1 ) + plan%right%n_halo
       Allocate( data_with_halo( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%right%swap_3d_x( lbd_h, temp, lbd_h, data_with_halo )

    Case( FILL_Y )
       lbd_h( 2 ) = lbd_h( 2 ) - plan%left%n_halo
       Allocate( temp( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%left%swap_3d_y( lbd, data, lbd_h, temp )
       
       ubd_h( 2 ) = ubd_h( 2 ) + plan%right%n_halo
       Allocate( data_with_halo( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%right%swap_3d_y( lbd_h, temp, lbd_h, data_with_halo )

    Case( FILL_Z )
       lbd_h( 3 ) = lbd_h( 3 ) - plan%left%n_halo
       Allocate( temp( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%left%swap_3d_z( lbd, data, lbd_h, temp )
       
       ubd_h( 3 ) = ubd_h( 3 ) + plan%right%n_halo
       Allocate( data_with_halo( lbd_h( 1 ):ubd_h( 1 ), lbd_h( 2 ):ubd_h( 2 ), lbd_h( 3 ):ubd_h( 3 ) ) )
       Call plan%right%swap_3d_z( lbd_h, temp, lbd_h, data_with_halo )

    End Select

  End Subroutine halo_dim_plan_fill_3d

  Subroutine halo_dim_plan_report( plan, unit )

    Class( halo_dim_plan_type ), Intent( In ) :: plan
    Integer                    , Intent( In ) :: unit
    
    Call plan%left%report ( unit )
    Call plan%right%report( unit )

  End Subroutine halo_dim_plan_report

  Subroutine plan_halo_swap_left( halo_plan, comm, n_tot, i_start, i_end, n_halo )

    ! Plan the comunications to complete a left halo - i.e. the lower halo

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_integer, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Implicit None

    Class( halo_plan_type ), Intent( InOut ) :: halo_plan
    Type ( mpi_comm )      , Intent( In    ) :: comm
    Integer                , Intent( In    ) :: n_tot
    Integer                , Intent( In    ) :: i_start
    Integer                , Intent( In    ) :: i_end
    Integer                , Intent( In    ) :: n_halo

    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Integer, Dimension( 1:2 ) :: want, wanted
    Integer, Dimension( 1:2 ) :: can_give, got
    Integer, Dimension( 1:2 ) :: dom_wanted

    Integer :: rank, nprc
    Integer :: prev, next
    Integer :: i_hold_lo, i_hold_hi
    Integer :: n_want, n_wanted
    Integer :: shift
    Integer :: error
    Integer :: i

    Type( halo_step_type ) :: this_step

    ! Work out my neighbouring procs
    Call mpi_comm_size( comm, nprc, error )
    Call mpi_comm_rank( comm, rank, error )
    prev = Modulo( rank - 1, nprc )
    next = Modulo( rank + 1, nprc )

    halo_plan%direction = LEFT
    
    halo_plan%comm = comm
    halo_plan%rank = rank
    halo_plan%prev = prev
    halo_plan%next = next

    halo_plan%i_start = i_start
    halo_plan%i_end   = i_end
    halo_plan%n_tot   = n_tot
    halo_plan%n_halo  = n_halo

    ! List of what indices need communicating at each step
    If( Allocated( halo_plan%steps ) ) Deallocate( halo_plan%steps )
    Allocate( halo_plan%steps( 1:0 ) )
    
    ! The current range of data I hold as the process progresses
    i_hold_lo = i_start
    i_hold_hi = i_end

    ! We will always have some halo to set up on the first step
    n_want   = Huge( n_want )
    n_wanted = Huge( n_wanted )
    want   = Huge( want )
    wanted = Huge( wanted )

    ! Communication step counter
    i = 0
    
    Do 
       i = i + 1

       ! Send to neighbouring proc the indices that I still need to fill in to complete the halo
       If( n_want > 0 ) Then
          want = [ i_start - n_halo, i_hold_lo - 1 ]
          !Send what is wanted    to   prev proc
          Call mpi_isend( want  , Size( want   ), mpi_integer, prev, 10, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       ! Recv what indices are wanted from next proc
       If( n_wanted > 0 ) Then
          Call mpi_irecv( wanted, Size( wanted ), mpi_integer, next, 10, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       ! Complete comms
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

       ! Work out how many data items I still want, and how many my neighbour wants
       n_want   = Max( 0, want  ( 2 ) - want  ( 1 ) + 1 )
       n_wanted = Max( 0, wanted( 2 ) - wanted( 1 ) + 1 )

       ! If neither I or my neighbour wants anything wqe are all done 
       If( n_want == 0 .And. n_wanted == 0 ) Exit

       ! If my neighbour needs some data ...
       If( n_wanted > 0 ) Then
          ! Given PBCs find where the indices wanted by my neighbour maps into our indexing
          ! Do this by finding where the upper bound maps into our indexing range, and then
          ! as it is a linear list of indices we can find the lower bound from the number of elements

          ! Do it painfully obviously as I am finding this tricky

          ! Find upper bound - note as I am sending left this should be next to my domain and so definitely in range
          dom_wanted( 2 ) = wanted( 2 )
          If( wanted( 2 ) < i_hold_lo ) Then
             Do
                dom_wanted( 2 ) = dom_wanted( 2 ) + n_tot
                If( dom_wanted( 2 ) >= i_hold_lo .And. dom_wanted( 2 ) <= i_hold_hi ) Exit
             End Do
          Else If( wanted( 2 ) > i_hold_hi ) Then
             Do
                dom_wanted( 2 ) = dom_wanted( 2 ) - n_tot
                If( dom_wanted( 2 ) >= i_hold_lo .And. dom_wanted( 2 ) <= i_hold_hi ) Exit
             End Do
          End If
          ! And hence the lower bound
          dom_wanted( 1 ) = dom_wanted( 2 ) - ( wanted( 2 ) - wanted( 1 ) )
          ! And given what is wanted expressed in our local indexing work out what we can actually give
          ! which is constrained by the current size of what we hold
          can_give = [ Max( dom_wanted( 1 ), i_hold_lo ), Min( i_end, dom_wanted( 2 ) ) ]
          !Send what can be given to   next proc
          Call mpi_isend( can_give, Size( can_give ), mpi_integer, next, 20, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If

       ! If I need some more data
       If( n_want > 0 ) Then
          ! Recv what can be given from prev_proc
          Call mpi_irecv( got     , Size( got      ), mpi_integer, prev, 20, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       
       ! Complete comms
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

       If( n_want > 0 ) Then
          ! The data sent to me will be in the indexing of the remote proc
          ! Shift back into the index range I want
          shift = want( 2 ) - got( 2 )
          got = got + shift
          i_hold_lo = got( 1 )
       End If

       ! Add this step of the plan to the list of what is needed to be done
       this_step%n_want   = n_want
       this_step%n_wanted = n_wanted
       this_step%can_give = can_give
       this_step%got      = got
       halo_plan%steps = [ halo_plan%steps, this_step ]

    End Do

  End Subroutine plan_halo_swap_left

  Subroutine plan_halo_swap_right( halo_plan, comm, n_tot, i_start, i_end, n_halo )

    ! Plan the comunications to complete a right halo - i.e. the upper halo

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_integer, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Implicit None

    Class( halo_plan_type ), Intent( InOut ) :: halo_plan
    Type ( mpi_comm )      , Intent( In    ) :: comm
    Integer                , Intent( In    ) :: n_tot
    Integer                , Intent( In    ) :: i_start
    Integer                , Intent( In    ) :: i_end
    Integer                , Intent( In    ) :: n_halo

    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Integer, Dimension( 1:2 ) :: want, wanted
    Integer, Dimension( 1:2 ) :: can_give, got
    Integer, Dimension( 1:2 ) :: dom_wanted

    Integer :: rank, nprc
    Integer :: prev, next
    Integer :: i_hold_lo, i_hold_hi
    Integer :: n_want, n_wanted
    Integer :: shift
    Integer :: error
    Integer :: i

    Type( halo_step_type ) :: this_step

    ! Work out my neighbouring procs
    Call mpi_comm_size( comm, nprc, error )
    Call mpi_comm_rank( comm, rank, error )
    prev = Modulo( rank - 1, nprc )
    next = Modulo( rank + 1, nprc )

    halo_plan%direction = RIGHT
    
    halo_plan%comm = comm
    halo_plan%rank = rank
    halo_plan%prev = prev
    halo_plan%next = next

    halo_plan%i_start = i_start
    halo_plan%i_end   = i_end
    halo_plan%n_tot   = n_tot
    halo_plan%n_halo  = n_halo

    ! List of what indices need communicating at each step
    If( Allocated( halo_plan%steps ) ) Deallocate( halo_plan%steps )
    Allocate( halo_plan%steps( 1:0 ) )
    
    ! The current range of data I hold as the process progresses
    i_hold_lo = i_start
    i_hold_hi = i_end

    ! We will always have some halo to set up on the first step
    n_want   = Huge( n_want )
    n_wanted = Huge( n_wanted )
    want   = Huge( want )
    wanted = Huge( wanted )

    ! Communication step counter
    i = 0
    
    Do 
       i = i + 1

       ! Send to neighbouring proc the indices that I still need to fill in to complete the halo
       If( n_want > 0 ) Then
          want = [ i_hold_hi + 1, i_end + n_halo ]
          !Send what is wanted    to   next proc
          Call mpi_isend( want  , Size( want   ), mpi_integer, next, PLAN_WANT_TAG, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       ! Recv what indices are wanted from prev proc
       If( n_wanted > 0 ) Then
          Call mpi_irecv( wanted, Size( wanted ), mpi_integer, prev, PLAN_WANT_TAG, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       ! Complete comms
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

       ! Work out how many data items I still want, and how many my neighbour wants
       n_want   = Max( 0, want  ( 2 ) - want  ( 1 ) + 1 )
       n_wanted = Max( 0, wanted( 2 ) - wanted( 1 ) + 1 )

       ! If neither I or my neighbour wants anything we are all done 
       If( n_want == 0 .And. n_wanted == 0 ) Exit

       ! If my neighbour needs some data ...
       If( n_wanted > 0 ) Then
          ! Given PBCs find where the indices wanted by my neighbour maps into our indexing
          ! Do this by finding where the upper bound maps into our indexing range, and then
          ! as it is a linear list of indices we can find the lower bound from the number of elements

          ! Do it painfully obviously as I am finding this tricky

          ! Find lower bound - note as I am sending right this should be next to my domain and so definitely in range
          dom_wanted( 1 ) = wanted( 1 )
          If( wanted( 1 ) < i_hold_lo ) Then
             Do
                dom_wanted( 1 ) = dom_wanted( 1 ) + n_tot
                If( dom_wanted( 1 ) >= i_hold_lo .And. dom_wanted( 1 ) <= i_hold_hi ) Exit
             End Do
          Else If( wanted( 1 ) > i_hold_hi ) Then
             Do
                dom_wanted( 1 ) = dom_wanted( 1 ) - n_tot
                If( dom_wanted( 1 ) >= i_hold_lo .And. dom_wanted( 1 ) <= i_hold_hi ) Exit
             End Do
          End If
          ! And hence the lower bound
          dom_wanted( 2 ) = dom_wanted( 1 ) + ( wanted( 2 ) - wanted( 1 ) )
          ! And given what is wanted expressed in our local indexing work out what we can actually give
          ! which is constrained by the current size of what we hold
          can_give = [ Max( dom_wanted( 1 ), i_start ), Min( i_hold_hi, dom_wanted( 2 ) ) ]
          !Send what can be given to prev proc
          Call mpi_isend( can_give, Size( can_give ), mpi_integer, prev, PLAN_GIVE_TAG, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If

       ! If I need some more data
       If( n_want > 0 ) Then
          ! Recv what can be given from next proc
          Call mpi_irecv( got     , Size( got      ), mpi_integer, next, PLAN_GIVE_TAG, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       
       ! Complete comms
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

       If( n_want > 0 ) Then
          ! The data sent to me will be in the indexing of the remote proc
          ! Shift back into the index range I want
          shift = want( 1 ) - got( 1 )
          got = got + shift
          i_hold_hi = got( 2 )
       End If

       ! Add this step of the plan to the list of what is needed to be done
       this_step%n_want   = n_want
       this_step%n_wanted = n_wanted
       this_step%can_give = can_give
       this_step%got      = got
       halo_plan%steps = [ halo_plan%steps, this_step ]

    End Do

  End Subroutine plan_halo_swap_right

  Subroutine swap_int_1d( plan, data, report )

    ! Complete a left halo update with integer data. Mainly for debugging.

    Class( halo_plan_type ),                 Intent( In    )              :: plan
    Integer                , Dimension( : ), Intent( InOut ), Allocatable :: data
    Logical                , Optional      , Intent( In    )              :: report

    Select Case( plan%direction )
    Case( LEFT )
       Call plan%swap_left_1d( data, report )
    Case( RIGHT )
       Call plan%swap_right_1d( data, report )
    End Select

  End Subroutine swap_int_1d

  Subroutine swap_real_3d_x( plan, lbd, data, lbd_h, data_with_halo )

    Use constants, Only : wp
    
    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Select Case( plan%direction )
    Case( LEFT )
       Call plan%swap_left_3d_x( lbd, data, lbd_h, data_with_halo )
    Case( RIGHT )
       Call plan%swap_right_3d_x( lbd, data, lbd_h, data_with_halo )
    End Select

  End Subroutine swap_real_3d_x

  Subroutine swap_real_3d_y( plan, lbd, data, lbd_h, data_with_halo )

    Use constants, Only : wp
    
    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Select Case( plan%direction )
    Case( LEFT )
       Call plan%swap_left_3d_y( lbd, data, lbd_h, data_with_halo )
    Case( RIGHT )
       Call plan%swap_right_3d_y( lbd, data, lbd_h, data_with_halo )
    End Select

  End Subroutine swap_real_3d_y

  Subroutine swap_real_3d_z( plan, lbd, data, lbd_h, data_with_halo )

    Use constants, Only : wp
    
    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Select Case( plan%direction )
    Case( LEFT )
       Call plan%swap_left_3d_z( lbd, data, lbd_h, data_with_halo )
    Case( RIGHT )
       Call plan%swap_right_3d_z( lbd, data, lbd_h, data_with_halo )
    End Select

  End Subroutine swap_real_3d_z

  Subroutine swap_int_1d_left( plan, data, report )

    ! Complete a left halo update with integer data. Mainly for debugging.

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_integer, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Class( halo_plan_type ),                 Intent( In    )              :: plan
    Integer                , Dimension( : ), Intent( InOut ), Allocatable :: data
    Logical                , Optional      , Intent( In    )              :: report

    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( : ), Allocatable :: data_with_halo

    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give
  
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: out
    Integer :: error
    Integer :: i, j

    Logical :: do_report

    If( Present( report ) ) Then
       do_report = report
    Else
       do_report = .False.
    End If

    out = 10 + plan%rank

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    Allocate( data_with_halo( Lbound( data, Dim = 1 ) - n_halo: Ubound( data, Dim = 1 ) ) )
    data_with_halo( Lbound( data, Dim = 1 ):Ubound( data, Dim = 1 ) ) = data

    If( do_report ) Then
       Write( out, * )
       Write( out, * ) "Initial"
       Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
       Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
       Write( out, * )
       Flush( out )
    End If
 
    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got

       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call mpi_isend( data_with_halo( can_give( 1 ) ), can_give( 2 ) - can_give( 1 ) + 1, mpi_integer, next, &
               SWAP_1D_LEFT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          Call mpi_irecv( data_with_halo( got( 1 ) ), got( 2 ) - got( 1 ) + 1, mpi_integer, prev, &
               SWAP_1D_LEFT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

    End Do

    Call move_alloc( data_with_halo, data )

    If( do_report ) Then
       Write( out, * )
       Write( out, * ) "Final"
       Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
       Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
       Flush( out )
    End If

  End Subroutine swap_int_1d_left

  Subroutine swap_int_1d_right( plan, data, report )

    ! Complete a right halo update with integer data. Mainly for debugging.

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_integer, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Class( halo_plan_type ),                 Intent( In    )              :: plan
    Integer                , Dimension( : ), Intent( InOut ), Allocatable :: data
    Logical                , Optional      , Intent( In    )              :: report

    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( : ), Allocatable :: data_with_halo

    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give
  
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: out
    Integer :: error
    Integer :: i, j

    Logical :: do_report

    If( Present( report ) ) Then
       do_report = report
    Else
       do_report = .False.
    End If

    out = 10 + plan%rank

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    Allocate( data_with_halo( Lbound( data, Dim = 1 ): Ubound( data, Dim = 1 ) + n_halo ) )
    data_with_halo( Lbound( data, Dim = 1 ):Ubound( data, Dim = 1 ) ) = data

    If( do_report ) Then
       Write( out, * )
       Write( out, * ) "Initial"
       Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
       Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
       Write( out, * )
       Flush( out )
    End If

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got

       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call mpi_isend( data_with_halo( can_give( 1 ) ), can_give( 2 ) - can_give( 1 ) + 1, mpi_integer, prev, &
               SWAP_1D_RIGHT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          Call mpi_irecv( data_with_halo( got( 1 ) ), got( 2 ) - got( 1 ) + 1, mpi_integer, next, &
               SWAP_1D_RIGHT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )

    End Do

    Call Move_alloc( data_with_halo, data )

    If( do_report ) Then
       Write( out, * )
       Write( out, * ) "Final"
       Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
       Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
       Flush( out )
    End If

  End Subroutine swap_int_1d_right

  Subroutine swap_real_3d_left_x( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: nx, n_loc_y, n_loc_z
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    ! Corners makes no difference in the x case
    n_loc_y = Size( data, Dim = 2 )
    n_loc_z = Size( data, Dim = 3 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data

    Allocate( buffer_send( 1:n_halo * n_loc_y * n_loc_z ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_y * n_loc_z ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( can_give, [ lbd( 2 ), ubd( 2 ) ], [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          nx = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, next, SWAP_3D_LEFT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          nx = got( 2 ) - got( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, prev, SWAP_3D_LEFT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          nx = got( 2 ) - got( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call copy_out( got, [ lbd( 2 ), ubd( 2 ) ], [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If
          
    End Do
    
  End Subroutine swap_real_3d_left_x
  
  Subroutine swap_real_3d_right_x( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: nx, n_loc_y, n_loc_z
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    ! Corners makes no difference in the x case
    n_loc_y = Size( data, Dim = 2 )
    n_loc_z = Size( data, Dim = 3 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data

    Allocate( buffer_send( 1:n_halo * n_loc_y * n_loc_z ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_y * n_loc_z ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( can_give, [ lbd( 2 ), ubd( 2 ) ], [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          nx = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, prev, SWAP_3D_RIGHT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          nx = got( 2 ) - got( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, next, SWAP_3D_RIGHT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          nx = got( 2 ) - got( 1 ) + 1
          n_msg = nx * n_loc_y * n_loc_z
          Call copy_out( got, [ lbd( 2 ), ubd( 2 ) ], [ lbd( 3 ), ubd( 3 )  ], &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If

    End Do
    
  End Subroutine swap_real_3d_right_x
  
  Subroutine swap_real_3d_left_y( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: ny, n_loc_x, n_loc_z
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )
    
    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    n_loc_x = Size( data, Dim = 1 )
    n_loc_z = Size( data, Dim = 3 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data

    Allocate( buffer_send( 1:n_halo * n_loc_x * n_loc_z ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_x * n_loc_z ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( [ lbd( 1 ),ubd( 1 ) ], can_give, [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          ny = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, next, SWAP_3D_LEFT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          ny = got( 2 ) - got( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, prev, SWAP_3D_LEFT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          ny = got( 2 ) - got( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call copy_out( [ lbd( 1 ), ubd( 1 ) ], got, [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If
          
    End Do
    
  End Subroutine swap_real_3d_left_y

  Subroutine swap_real_3d_right_y( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: ny, n_loc_x, n_loc_z
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )
    
    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    n_loc_x = Size( data, Dim = 1 )
    n_loc_z = Size( data, Dim = 3 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data
    
    Allocate( buffer_send( 1:n_halo * n_loc_x * n_loc_z ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_x * n_loc_z ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( [ lbd( 1 ), ubd( 1 ) ], can_give, [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          ny = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, prev, SWAP_3D_RIGHT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          ny = got( 2 ) - got( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, next, SWAP_3D_RIGHT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          ny = got( 2 ) - got( 1 ) + 1
          n_msg = ny * n_loc_x * n_loc_z
          Call copy_out( [ lbd( 1 ), ubd( 1 ) ], got, [ lbd( 3 ), ubd( 3 ) ], &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If
          
    End Do
    
  End Subroutine swap_real_3d_right_y

  Subroutine swap_real_3d_left_z( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo
    
    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: nz, n_loc_x, n_loc_y
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    n_loc_x = Size( data, Dim = 1 )
    n_loc_y = Size( data, Dim = 2 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data

    Allocate( buffer_send( 1:n_halo * n_loc_x * n_loc_y ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_x * n_loc_y ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( [ lbd( 1 ), ubd( 1 ) ], [ lbd( 2 ), ubd( 2 ) ], can_give, &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          nz = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, next, SWAP_3D_LEFT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          nz = got( 2 ) - got( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, prev, SWAP_3D_LEFT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          nz = got( 2 ) - got( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call copy_out( [ lbd( 1 ), ubd( 1 ) ], [ lbd( 2 ), ubd( 2 ) ], got, &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If
          
    End Do
    
  End Subroutine swap_real_3d_left_z

  Subroutine swap_real_3d_right_z( plan, lbd, data, lbd_h, data_with_halo )

    Use mpi_f08, Only : mpi_comm, mpi_request_null, mpi_request, mpi_statuses_ignore, &
         mpi_comm_size, mpi_comm_rank, mpi_isend, mpi_irecv, mpi_waitall

    Use constants, Only : wp

    Class( halo_plan_type ),                                                     Intent( In    ) :: plan
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd
    Real( wp )             , Dimension( lbd( 1 )  :, lbd( 2 )  :, lbd( 3 )  : ), Intent( In    ) :: data
    Integer                , Dimension( 1:3                                   ), Intent( In    ) :: lbd_h
    Real( wp )             , Dimension( lbd_h( 1 ):, lbd_h( 2 ):, lbd_h( 3 ): ), Intent(   Out ) :: data_with_halo

    Real( wp ), Dimension( : ), Allocatable :: buffer_send
    Real( wp ), Dimension( : ), Allocatable :: buffer_recv
    
    Type( mpi_request ), Dimension( 1:2 ) :: requests

    Type( mpi_comm ) :: comm

    Integer, Dimension( 1:3 ) :: ubd
    Integer, Dimension( 1:2 ) :: got
    Integer, Dimension( 1:2 ) :: can_give

    Integer :: nz, n_loc_x, n_loc_y
    Integer :: n_msg
    Integer :: prev, next
    Integer :: n_want, n_wanted
    Integer :: n_halo
    Integer :: error
    Integer :: i

    ubd = Ubound( data )

    comm   = plan%comm
    prev   = plan%prev
    next   = plan%next
    n_halo = plan%n_halo

    n_loc_x = Size( data, Dim = 1 )
    n_loc_y = Size( data, Dim = 2 )

    data_with_halo( lbd( 1 ):ubd( 1 ), lbd( 2 ):ubd( 2 ), lbd( 3 ):ubd( 3 ) ) = data

    Allocate( buffer_send( 1:n_halo * n_loc_x * n_loc_y ) )
    Allocate( buffer_recv( 1:n_halo * n_loc_x * n_loc_y ) )

    Do i = 1, Size( plan%steps )

       n_want = plan%steps( i )%n_want
       n_wanted = plan%steps( i )%n_wanted
       can_give = plan%steps( i )%can_give
       got = plan%steps( i )%got
       
       ! Send out data and Recieve new data
       If( n_wanted > 0 ) Then
          Call copy_in( [ lbd( 1 ), ubd( 1 ) ], [ lbd( 2 ), ubd( 2 ) ], can_give, &
               Lbound( data_with_halo ), data_with_halo, buffer_send )
          nz = can_give( 2 ) - can_give( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call mpi_isend( buffer_send, n_msg, plan%real_handle, prev, SWAP_3D_RIGHT, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          nz = got( 2 ) - got( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call mpi_irecv( buffer_recv, n_msg, plan%real_handle, next, SWAP_3D_RIGHT, comm, requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       If( n_want > 0 ) Then
          nz = got( 2 ) - got( 1 ) + 1
          n_msg = nz * n_loc_x * n_loc_y
          Call copy_out( [ lbd( 1 ), ubd( 1 ) ], [ lbd( 2 ), ubd( 2 ) ], got, &
               Lbound( data_with_halo ), buffer_recv, data_with_halo )
       End If
          
    End Do
    
  End Subroutine swap_real_3d_right_z

  Impure Elemental Subroutine halo_dim_plan_free( plan )

    Use mpi_f08, Only : mpi_type_free

    Class( halo_dim_plan_type ), Intent( InOut ) :: plan

    Integer :: error

    Call mpi_type_free( plan%left%real_handle, error )

    Deallocate( plan%right%steps )
    Deallocate( plan%left%steps )
    
  End Subroutine halo_dim_plan_free
  
  Impure Elemental Subroutine halo_dim_plan_final( plan )

    Type( halo_dim_plan_type ), Intent( InOut ) :: plan

    Call plan%free()

  End Subroutine halo_dim_plan_final

  Subroutine report_plan( plan, unit )

    ! Report a halo plan

    Class( halo_plan_type ), Intent( In ) :: plan
    Integer                , Intent( In ) :: unit

    Integer :: i

    Write( unit, * )
    Write( unit, * ) 'Steps in plan for rank ', plan%rank, ' neighbours ', plan%prev, plan%next
    Write( unit, * ) 'Plan is for a ', Merge( 'Left ', 'Right', plan%direction == LEFT ), ' swap'
    Write( unit, '( 3( a7 ), 2( a14 ) )' ) 'Step   ', 'Want   ', 'Wanted ', '  Can give    ', '  Got         ' 
    Do i = 1, Size( plan%steps )
       Write( unit, '( 7( i4, 3x ) )' ) &
            i, plan%steps( i )%n_want, plan%steps( i )%n_wanted, plan%steps( i )%can_give, plan%steps( i )%got
    End Do
    Write( unit, * )

  End Subroutine report_plan

  Pure Subroutine copy_in( xb, yb, zb, lb, in, buff )

    Use constants, Only : wp

    Integer   , Dimension( 1:2     ), Intent( In    ) :: xb
    Integer   , Dimension( 1:2     ), Intent( In    ) :: yb
    Integer   , Dimension( 1:2     ), Intent( In    ) :: zb
    Integer   , Dimension( 1:3     ), Intent( In    ) :: lb
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: in
    Real( wp ), Dimension( :       ), Intent(   Out ) :: buff

    Integer :: ib
    Integer :: ix, iy, iz

    ib = 0
    Do iz = zb( 1 ), zb( 2 )
       Do iy = yb( 1 ), yb( 2 )
          Do ix = xb( 1 ), xb( 2 )
             ib = ib + 1
             buff( ib ) = in( ix, iy, iz )
          End Do
       End Do
    End Do
       
  End Subroutine copy_in
    
  Pure Subroutine copy_out( xb, yb, zb, lb, buff, out )

    Use constants, Only : wp

    Integer   , Dimension( 1:2     ), Intent( In    ) :: xb
    Integer   , Dimension( 1:2     ), Intent( In    ) :: yb
    Integer   , Dimension( 1:2     ), Intent( In    ) :: zb
    Integer   , Dimension( 1:3     ), Intent( In    ) :: lb
    Real( wp ), Dimension( :       ), Intent( In    ) :: buff
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: out

    Integer :: ib
    Integer :: ix, iy, iz

    ib = 0
    Do iz = zb( 1 ), zb( 2 )
       Do iy = yb( 1 ), yb( 2 )
          Do ix = xb( 1 ), xb( 2 )
             ib = ib + 1
             out( ix, iy, iz ) = buff( ib )
          End Do
       End Do
    End Do
       
  End Subroutine copy_out
  
End Module swap_module
