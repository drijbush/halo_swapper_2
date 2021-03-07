Module halo_setter_base_module

  Use constants, Only : wp
  Implicit None

  Type, Public, Abstract :: halo_setter_base_class
     Private
     Integer :: n_calls
   Contains
     Generic,                     Public  :: init        => base_init
     Procedure,                   Public  :: allocate    => halo_allocate
     Procedure,                   Public  :: free        => halo_free
     Procedure,                   Public  :: inc_n_calls => halo_inc_n_calls
     Procedure( fill ), Deferred, Public  :: fill
     Procedure,                   Private :: base_init
  End type halo_setter_base_class

  Abstract Interface

     Subroutine fill( H, halo_width, hdlb, gin, hout, error )
       Use constants, Only : wp
       Import :: halo_setter_base_class
       Class( halo_setter_base_class ),                             Intent( InOut ) :: H
       Integer,                                                     Intent( In    ) :: halo_width
       Integer,    Dimension( 1:3 ),                                Intent( In    ) :: hdlb
       Real( wp ), Dimension( 0:, 0:, 0: ),                         Intent( In    ) :: gin
       Real( wp ), Dimension( hdlb( 1 ):, hdlb( 2 ):, hdlb( 3 ): ), Intent(   Out ) :: hout
       Integer,                                                     Intent(   Out ) :: error
     End Subroutine fill

  End Interface

Contains

  Subroutine base_init( H, error )

    Class( halo_setter_base_class      ), Intent( InOut ) :: H
    Integer,                              Intent(   Out ) :: error

    H%n_calls = 0

    error = 0

  End Subroutine base_init

  Pure Subroutine halo_allocate( H, glb, gub, halo_width, hout )

    Class( halo_setter_base_class ),               Intent( InOut ) :: H
    Integer,    Dimension( 1:3     ),              Intent( In    ) :: glb
    Integer,    Dimension( 1:3     ),              Intent( In    ) :: gub
    Integer,                                       Intent( In    ) :: halo_width
    Real( wp ), Dimension( :, :, : ), Allocatable, Intent(   Out ) :: hout

    Integer, Dimension( 1:3 ) :: hlb, hub

    H%n_calls = 0

    hlb       = glb - halo_width
    hub       = gub + halo_width
    Allocate( hout( hlb( 1 ):hub( 1 ), hlb( 2 ):hub( 2 ), hlb( 3 ):hub( 3 ) ) )

  End Subroutine halo_allocate

  Pure Subroutine halo_free( H, hout )

    Class( halo_setter_base_class ),               Intent( InOut ) :: H
    Real( wp ), Dimension( :, :, : ), Allocatable, Intent(   Out ) :: hout

    H%n_calls = 0

    Deallocate( hout )

  End Subroutine halo_free

  Subroutine halo_inc_n_calls( H )

    Class( halo_setter_base_class ), Intent( InOut ) :: H

    H%n_calls = H%n_calls + 1

  End Subroutine halo_inc_n_calls

End Module halo_setter_base_module
