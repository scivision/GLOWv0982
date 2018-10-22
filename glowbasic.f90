program glowbasic

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file Glowlicense.txt.
! For more information see the file Glow.txt.

! Version 0.982, 2/2018

! Adapted from glowdriver by Stan Solomon, 2/2016

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
! Runs GLOW for designated inputs once, or multiple times.
! MPI and netCDF libraries not required.

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.

! Other definitions:
! f107p    Solar 10.7 cm flux for previous day
! ap       Ap index of geomagnetic activity
! stl      Local Solar Time, hours
! z        altitude array, km
! ef       auroral electron maxwellian flux, erg/cm2/s
! ec       auroral electron maxwellian characteristic energy, eV
! fmono    auroral electron monoenergetic flux, erg/cm2/s
! emono    auroral electron monoenergetic energy, eV
! pfmax    auroral proton maxwellian energy flux, erg/cm2/s
! pemax    auroral proton maxwellian characteristic energy, eV
! pfmono   auroral proton monoenergetic flux, erg/cm2/s
! pemono   auroral proton monoenergetic energy, eV
! itail    Switch for including low-energy tail of Maxwellian distributions
! pedcond  Pederson conductivities in S/m (mho)
! hallcond Hall conductivities in S/m (mho)
! zun,zvn  neutral wind components (not in use)

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission

  use cglow,only: jmax,nbins,npbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: pflux,pener,pdel,pia
  use cglow,only: cglow_init
  use cglow,only: data_dir

  implicit none

  character(len=1024) :: iri90_dir

  real,allocatable :: z(:)
  real,allocatable :: zun(:), zvn(:)
  real,allocatable :: pedcond(:), hallcond(:)
  real :: stl,fmono,emono,pfmax,pemax,pfmono,pemono
  integer :: i,j,itail,instance,iostatus

!
! Initialize standard switches:
!
  iscale=1
  xuvfac=3.
  kchem=4
  jlocal=0
  itail=0
  fmono=0.
  emono=0.
  pfmono=0.
  pemono=0.
!
! Set data directories:
!
  data_dir    = 'data/'
  iri90_dir   = 'data/iri90/'
!
! Set number of altitude levels:
!
  jmax = 102
!
! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(pedcond(jmax))
  allocate(hallcond(jmax))
!
! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
  call cglow_init
!
! Call EGRID to set up electron energy grid:
!
  call egrid (ener, del, nbins)
!
! Call PEGRID to set up proton energy grid:
!
  call pegrid (pener, pdel, npbins)
!
! Loop to call GLOW for designated inputs until end-of-file or any character on standard input:
!
  do instance=1,10000
!
! Get input values:
!
    write(6,"('Enter date, UT, lat, lon, F107a, F107, F107p, Ap, Ef, Ec, PFmax, PEmax')")
    read(5,*,iostat=iostatus) idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,pfmax,pemax
    if (iostatus /= 0) stop
!
! Calculate local solar time:
!
    stl = ut/3600. + glong/15.
    if (stl < 0.) stl = stl + 24.
    if (stl >= 24.) stl = stl - 24.
!
! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
!
    call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                 z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Call MAXT to put auroral electron flux specified by input into phitop array:
!
    phitop(:) = 0.
    if ((ef>.001 .and. ec>1.) .or. (fmono>.001 .and. emono>1.)) &
      call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
!
! Call MAXT to put auroral proton flux specified by input into pflux array:
!
    pflux(:) = 0.
    if ((pfmax>.001 .and. pemax>1.) .or. (pfmono>.001 .and. pemono>1.)) &
      call maxt (pfmax,pemax,pener,pdel,npbins,itail,pfmono,pemono,pflux)
!
! Fill altitude array, converting to cm:
!
    zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!
    call glow
!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
    do j=1,jmax
      call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                    zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
                    pedcond(j), hallcond(j))
    enddo
!
! Output section:
!
    write(6,"(1x,i7,10f8.1,f10.1)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,pfmax,pemax
    write(6,"('   Z     Tn       O        N2        NO      Ne(in)    Ne(out)  Ionrate      O+       O2+      NO+       N(2D)    Pederson   Hall')")
    write(6,"(1x,0p,f5.1,f6.0,1p,12e10.2)") (z(j),ztn(j),zo(j),zn2(j),zno(j),ze(j), &
      ecalc(j),tir(j),zxden(3,j),zxden(6,j),zxden(7,j),zxden(10,j),pedcond(j),hallcond(j),j=1,jmax)
    write(6,"('   Z      3371    4278    5200    5577    6300    7320   10400    3644    7774    8446    3726    LBH     1356    1493    1304')")
    write(6,"(1x,f5.1,15f8.2)")(z(j),(zeta(i,j),i=1,15),j=1,jmax)

  enddo

stop

end program glowbasic
