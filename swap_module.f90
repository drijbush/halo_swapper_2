Module swap_module

  Use mpi_f08, Only : mpi_comm
  
  Implicit None

  Type, Private :: halo_step_type
     Private
     Integer                   :: n_want, n_wanted
     Integer, Dimension( 1:2 ) :: can_give
     Integer, Dimension( 1:2 ) :: got
  End Type halo_step_type

  Type, Public :: halo_plan_type
     Private
     Type( mpi_comm )                                    :: comm
     Integer                                             :: rank
     Integer                                             :: prev, next
     Integer                                             :: direction
     Integer                                             :: i_start, i_end
     Integer                                             :: n_tot
     Integer                                             :: n_halo
     Type( halo_step_type ), Dimension( : ), Allocatable :: steps
   Contains
     Procedure, Public :: plan_left  => plan_halo_swap_left
     Procedure, Public :: plan_right => plan_halo_swap_right
     Procedure, Public :: report     => report_plan
     Procedure, Public :: swap_left  => swap_int_1d_left
     Procedure, Public :: swap_right => swap_int_1d_right
     Procedure, Public :: swap       => swap_int_1d
  End Type halo_plan_type

  Type, Public :: halo_dim_plan_type
     Private
     Type( halo_plan_type ) :: left
     Type( halo_plan_type ) :: right
   Contains
     Procedure, Public :: init => halo_dim_plan_init
     Procedure, Public :: fill => halo_dim_plan_fill
  End Type halo_dim_plan_type

  Private

  Integer, Parameter :: LEFT  = -1
  Integer, Parameter :: RIGHT = +1

Contains

  Subroutine halo_dim_plan_init( plan, local_size, halo_width, comm, error )

    Use mpi_f08, Only : mpi_comm, mpi_allreduce, mpi_scan, mpi_integer, mpi_sum
    
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

    Call plan%left%plan_left  ( comm, n_tot, i_start, i_end, halo_width )
    Call plan%right%plan_right( comm, n_tot, i_start, i_end, halo_width )
    
  End Subroutine halo_dim_plan_init

  Subroutine halo_dim_plan_fill( plan, data, report )

    Class( halo_dim_plan_type ),                 Intent( In    )              :: plan
    Integer                    , Dimension( : ), Intent( InOut ), Allocatable :: data
    Logical                    , Optional      , Intent( In    )              :: report

    Call plan%left%swap ( data, report )
    Call plan%right%swap( data, report )

  End Subroutine halo_dim_plan_fill

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
!!$          want = [ i_start - n_halo, i_hold_lo - 1 ]
          want = [ i_hold_hi + 1, i_end + n_halo ]
          !Send what is wanted    to   next proc
          Call mpi_isend( want  , Size( want   ), mpi_integer, next, 10, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       ! Recv what indices are wanted from prev proc
       If( n_wanted > 0 ) Then
          Call mpi_irecv( wanted, Size( wanted ), mpi_integer, prev, 10, comm, requests( 2 ), error )
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
          Call mpi_isend( can_give, Size( can_give ), mpi_integer, prev, 20, comm, requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If

       ! If I need some more data
       If( n_want > 0 ) Then
          ! Recv what can be given from next proc
          Call mpi_irecv( got     , Size( got      ), mpi_integer, next, 20, comm, requests( 2 ), error )
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
       Call plan%swap_left( data, report )
    Case( RIGHT )
       Call plan%swap_right( data, report )
    End Select

  End Subroutine swap_int_1d

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

    comm = plan%comm
    prev = plan%prev
    next = plan%next

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
          Call mpi_isend( data( can_give( 1 ) ), can_give( 2 ) - can_give( 1 ) + 1, mpi_integer, next, 30, comm, &
               requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
          Allocate( data_with_halo( got( 1 ):Ubound( data, Dim = 1 ) ) )
          data_with_halo = -100
          data_with_halo( Lbound( data, Dim = 1 ): ) = data
          Call mpi_irecv( data_with_halo( got( 1 ) ), got( 2 ) - got( 1 ) + 1, mpi_integer, prev, 30, comm, &
               requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       ! Put new data in resized old array
       If( n_want > 0 ) Then
          Call move_alloc( data_with_halo, data )
       End If

       ! Report current status of data
       If( do_report ) Then
          Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
          Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
          Write( out, * )
          Flush( out )
       End If

    End Do

    If( do_report ) Then
       Write( out, * )
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

    comm = plan%comm
    prev = plan%prev
    next = plan%next

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
          Call mpi_isend( data( can_give( 1 ) ), can_give( 2 ) - can_give( 1 ) + 1, mpi_integer, prev, 30, comm, &
               requests( 1 ), error )
       Else
          requests( 1 ) = mpi_request_null
       End If
       If( n_want > 0 ) Then
!!$          Allocate( data_with_halo( got( 1 ):Ubound( data, Dim = 1 ) ) )
          Allocate( data_with_halo( Lbound( data, Dim = 1 ):got( 2 ) ) )
          data_with_halo = -100
          data_with_halo( Lbound( data, Dim = 1 ):got( 1 ) - 1 ) = data
          Call mpi_irecv( data_with_halo( got( 1 ) ), got( 2 ) - got( 1 ) + 1, mpi_integer, next, 30, comm, &
               requests( 2 ), error )
       Else
          requests( 2 ) = mpi_request_null
       End If
       Call mpi_waitall( Size( requests ), requests, mpi_statuses_ignore, error )
       ! Put new data in resized old array
       If( n_want > 0 ) Then
          Deallocate( data )
          Allocate( data, Mold = data_with_halo )
          data = data_with_halo
          Deallocate( data_with_halo )
       End If

       ! Report current status of data
       If( do_report ) Then
          Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
          Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
          Write( out, * )
          Flush( out )
       End If

    End Do

    If( do_report ) Then
       Write( out, * )
       Write( out, * )
       Write( out, * ) "Final"
       Write( out, '( "Indices: ", 3x, 200( i3, 1x ) )' ) ( j, j = Lbound( data, 1 ), Ubound( data, 1 ) )
       Write( out, '( "Current: ", i3, 200( i3, 1x ) )' ) Size( data ), data
       Flush( out )
    End If

  End Subroutine swap_int_1d_right

  Subroutine report_plan( plan )

    ! Report a halo plan

    Class( halo_plan_type ), Intent( In ) :: plan

    Integer :: out
    Integer :: i

    out = 10 + plan%rank

    Write( out, * )
    Write( out, * ) 'Steps in plan for rank ', plan%rank, ' neighbours ', plan%prev, plan%next
    Write( out, * ) 'Plan is for a ', Merge( 'Left ', 'Right', plan%direction == LEFT ), ' swap'
    Write( out, '( 3( a7 ), 2( a14 ) )' ) 'Step   ', 'Want   ', 'Wanted ', '  Can give    ', '  Got         ' 
    Do i = 1, Size( plan%steps )
       Write( out, '( 7( i4, 3x ) )' ) &
            i, plan%steps( i )%n_want, plan%steps( i )%n_wanted, plan%steps( i )%can_give, plan%steps( i )%got
    End Do

  End Subroutine report_plan
  
End Module swap_module
