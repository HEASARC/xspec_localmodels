c idre - Integrating Disc Ring Emission
c Fortran77 integration subroutine for XSPEC
c
c ref. Dovciak M., Karas V., Yaqoob T. (2004)
c -------------------------------------------
c REFERENCES:
c
c Dovciak M., Karas V. & Yaqoob, T. (2004). An extended scheme for fitting X-ray
c data with accretion disk spectra in the strong gravity regime. ApJS, 153, 205.
c
c Dovciak M., Karas V., Martocchia A., Matt G. & Yaqoob T. (2004). XSPEC model
c to explore spectral features from black hole sources. In Proc. of the workshop
c on processes in the vicinity of black holes and neutron stars. S.Hledik &
c Z.Stuchlik, Opava. In press. [astro-ph/0407330]
c
c Dovciak M. (2004). Radiation of accretion discs in strong gravity. Faculty of
c Mathematics and Physics, Charles University, Prague. PhD thesis.
c [astro-ph/0411605]
c -------------------------------------------
c
c This subroutine integrates local emission for axisymmetric and stationary
c emission of the accretion disc near rotating (Kerr) black hole (characterized
c by angular momentum a/M) for an observer with inclination angle theta_o.
c
c Axially integrated function comprising of g-factor, lensing, cosine of
c emission angle and limb darkening/brightening law is needed for the
c integration. This function differs for different a/M and theta_o. It is read
c and interpolated for particular a/M and theta_o from the fits file called
c 'KBHlineNN.fits', where NN is '00', '01', '02',... The format of this fits
c file is described in detail below.
c
c By 'local' it is meant 'with respect to the local inertial frame
c connected with the fluid in the accretion disc' everywhere in this code
c (local emission, local angle of emitted ray, local disc normal, ...).
c The local frame should be defined with each KBHlineNN.fits because some of
c the transfer functions used for computation of the tables are calculated
c with respect to it (g-factor, cosine).
c
c
c The subroutine 'idre' has 8 parameters:
c ear - array of energy bins (same as ear for local models in XSPEC)
c ne  - number of energy bins (same as ne for local models in XSPEC)
c photar(ne) - array of photon number density flux per bin
c              (same as photar for local models in XSPEC)
c idre_param - 10 more parameters needed for integration (explained below)
c cmodel     - 32-byte string with a base name of a FITS file with tables for
c              axisymmetric emission (e.g. "KBHline" for KBHlineNN.fits)
c ne_loc  - number of points (in energies) where local photon number density
c           flux (per keV) is defined in emissivity subroutine
c ear_loc - array of the local energies where the local photon flux far_loc is
c           defined
c far_loc - array of the local photon flux (per keV)
c
c For an example on how to use 'idre' and 'emissivity' subroutines in local
c models see e.g. gaussian line model "xskyrline.f".
c
c -------------------------------------------
c idre_param:
c
c idre_param(1): am      - black hole angular momentum (0 <= am <= 1)
c idre_param(2): thetaO  - observer inclination in degrees (0-pole, 90-equator)
c idre_param(3): rin     - inner edge of non-zero disc emissivity in GM/c^2
c int(idre_param(4)): ms - if r_in < r_ms, whether to integrate from r_in (ms=0)
c                          or from marginally stable orbit (MSO)
c                          (ms=1, i.e. zero emissivity below MSO)
c idre_param(5): rout    - outer edge of non-zero disc emissivity in GM/c^2
c                        - if rin > rout => zero total flux
c int(idre_param(6)): smooth  - whether to smooth the resulting spectrum
c                               (0-no, 1-yes)
c idre_param(7): normal - how to normalize the final spectrum
c                         = 0. - normalization to unity (total flux=1.)
c                                (e.g. for line)
c                         > 0. - normalization to the flux at 'normal' keV
c                                (usually used for continuum)
c                         < 0. - final spectrum is not normalized in idre
c idre_param(8):      zshift - overall Doppler shift
c int(idre_param(9)): ntable - table of integrated transfer function
c                              (defines fits file with tables), 0<= ntable <= 99
c idre_param(10):     alpha  - the inner radial power-law index used below boundary
c                              radius rb
c idre_param(11):     beta   - the outer radial power-law index used above boundary
c                              radius rb
c idre_param(12):     rb     - the boundary radius between inner and outer radial
c                              power-law emissivity (in units of GM/c^2)
c
c NOTE: Unlike the subroutine ide, the subroutine idre does not need any
c       external emissivity subroutine
c
c -------------------------------------------
c Integrated transfer function in the fits files KBHlineNN.fits
c
c Values of integrated transfer function are stored here as tables for different
c values of horizon (r_horizon, not a/M!) and observer inclination angle
c theta_o. Each table consists of values of the function for different radii r
c and g-factors (i.e. observed energy). Particular r_horizon, theta_o, r and
c g-factor, where the function is evaluated, are defined at the beginning of the
c fits file as vectors.
c
c Definition of the file KBHlineNN.fits:
c 0. All of the extensions defined below are binary.
c 1. The first extension contains one row with three columns that define bins in
c    the g-factor:
c    - integer in the first column defines the width of the bins (0 - constant,
c      1 - exponentially growing),
c    - real number in the second column defines the lower boundary of the first
c      bin (minimum of the g-factor),
c    - real number in the third column defines the upper boundary of the last
c      bin (maximum of the g-factor).
c 2. The second extension contains a vector of the values of the radius relative
c    to the horizon, r-r_horizon, in GM/c2 .
c 3. The third extension contains a vector of the horizon values in GM/c2
c    (1.00 <= r_horizon <= 2.00).
c 4. The fourth extension contains a vector of the values of the observer's
c    inclination angle theta_o in degrees (0 <= theta_o <= 90, 0-pole, 90-disc).
c 5. All the previous vectors have to have values sorted in an increasing order
c 6. In the following extensions the integrated transfer function is stored,
c    each extension is for a particular value of r_horizon and theta_o. The
c    values of r_horizon and theta_o are changing with each extension in the
c    following order:
c    r_horizon[1] x theta_o[1],
c    r_horizon[1] x theta_o[2],
c    r_horizon[1] x theta_o[3],
c    ...
c    ...
c    r_horizon[2] x theta_o[1],
c    r_horizon[2] x theta_o[2],
c    r_horizon[2] x theta_o[3],
c    ...
c    ...
c    Each extension has one column.
c 7. Each row corresponds to a particular value of r-r_horizon (see 2. above).
c 8. Each element corresponding to a particular column and row is a vector. Each
c    element of this vector corresponds to a value of the function for a
c    particular bin in the g-factor. This bin can be calculated from the number
c    of the corresponding element of the vector and data from the first
c    extension (see 1. above).
c
c For an example of a fits file with integrated transfer function see
c KBHline00.fits.
c
c -----
c KBHline00.fits, KBHline01.fits, KBHline02.fits
c
c These tables are computed for an optically thick and geometrically thin
c accretion disc near Kerr black hole. The medium between the disc and the
c observer is supposed to be optically thin for the wavelengths one is
c interested in. Therefore ray tracing in vacuum Kerr space-time could be used.
c
c When calculating the integrated transfer function, it was supposed that the
c fluid in the disc rotates on stable circular (free) orbits above marginally
c stable orbit (MSO). Below MSO the fluid is freely falling and has the same
c energy and angular momentum as the matter which is on the MSO.
c
c The observer is placed in the direction phi = pi/2. The black hole rotates
c counter-clockwise.
c
c Three sets of tables for different limb darkening/brightening
c laws are precalculated. All of them were computed from tables in the
c KBHtables00.fits file and therefore these tables are calculated for the same
c values of the black-hole horizon and observer's inclination:
c - values of r_horizon are: 1.00, 1.05, 1.10, 1.15, ..., 1.90, 1.95, 2.00
c   (21 elements)
c - values of theta_o are: 0.1, 1, 5, 10, 15, 20, 25, ..., 75, 80, 85, 89
c   (20 elements)
c All of these tables have equidistant bins in the g-factor which fall in the
c interval 0.001-1.7 .
c Three sets of tables are available:
c - KBHline00.fits for isotropic emission
c - KBHline01.fits for Laor's limb darkening (flux ~ (1+2.06*mu); mu = cosine of
c   emission angle),
c - KBHline02.fits for Haardt's limb brightening (flux ~ ln (1+1/mu)).
c All of these tables have 300 bins in the g -factor and 500 values of the
c radius r-r_horizon which are exponentially increasing from 0 to 999.
c
c  1.6.2007 => we have changed format of the data file with transfer function
c 14.5.2008 => added broken radial power-law emissivity
c           => inner and outer edge is no longer relative to the horizon
c
c===============================================================================
      subroutine idre(ear,ne,photar,idre_param,cmodel,
     $                ne_loc,ear_loc,far_loc)
c-----------------------variables declarations--------------------------------
      implicit real(a-h,o-z)

      INTERFACE
         SUBROUTINE idre_ini(dc,r_vec,nr,am,thetaO,ntable,
     $                       cmodel,edivision,min,max,nbins_final)
           real, dimension(:), allocatable, intent(inout) :: dc(:)
           real, dimension(:), allocatable, intent(inout) :: r_vec(:)
           integer, intent(inout) :: nr
           real, intent(in) :: am
           real, intent(in) :: thetaO
           integer, intent(inout) :: ntable
           character(len=*), intent(in) :: cmodel
           integer, intent(out) :: edivision
           real, intent(out) :: min
           real, intent(out) :: max
           integer, intent(out) :: nbins_final
         END SUBROUTINE idre_ini
      END INTERFACE

      real ear(0:ne),idre_param(*),photar(ne)
      real ear_loc(0:ne_loc),far_loc(ne_loc)
      double precision fc(ne),dfc,photar_int,rbb,ralpha,photar1(ne)
      integer edivision,smooth,ntable,ms,nbins
      real am_old,thetaO_old,normal
      real min,max,Enorm,gmin,gmax,gmin_loc,gmax_loc
      real alpha,beta,rb
      integer ntable_old,nbins_old
      character(len=32) cmodel
      character*(32) kyrh,kyrin,kyrms
      character*(74) errortxt
      character*(80) errortext
      character*(128) fgmstr,pname,pkyrin,pkyrms,pkyrh
      character*(255) kydir

c     pointers to arrays used in integration...
      external fgmstr,fpmstr

      real, dimension(:), allocatable :: dc(:)
      real, dimension(:), allocatable :: r_vec(:)
      real(8), dimension(:), allocatable :: emission(:)

      save

      data pname,pkyrh,pkyrin,pkyrms/'KYDIR','KYRH','KYRIN','KYRMS'/
      data kydir/''/
      data am_old,thetaO_old/2*-1./
      data ntable_old,nbins_old/2*-1/
      data pi/3.1415926535897932/
      data pi2/6.2831853071795865/

c------------------------execution starts here--------------------------------
c let's check if input parameters have reasonable values
      am=idre_param(1)
      if(am.lt.0..or.am.gt.1.)then
       ierr=1
       errortxt='a/M must be >= 0 and <= 1'
       goto 9999
      endif
      am2=am*am
c outer and inner horizons
      r_plus=1.+sqrt(1.-am2)
c      r_minus=1.-sqrt(1.-am2)
      thetaO=idre_param(2)
      if(thetaO.lt.0..or.thetaO.gt.90.)then
       ierr=1
       errortxt='theta_o must be >= 0 and <= 90'
       goto 9999
      endif
      rin=idre_param(3)-r_plus
      if(rin.lt.0.)rin=0.
      ms=int(idre_param(4))
      if(ms.ne.0.and.ms.ne.1)then
       ierr=1
       errortxt='ms must be 0 or 1'
       goto 9999
      endif
      rout=idre_param(5)-r_plus
      if(rout.lt.0.)rout=0.
      smooth=int(idre_param(6))
      if(smooth.ne.0.and.smooth.ne.1)then
       ierr=1
       errortxt='smooth must be 0 or 1'
       goto 9999
      endif
      normal=idre_param(7)
      zshift=idre_param(8)
      if(zshift.le.-1.)then
       ierr=1
       errortxt='zshift must be larger than -1'
       goto 9999
      endif
      ntable=int(idre_param(9))
      if(ntable.lt.0.or.ntable.gt.99)then
       ierr=1
       errortxt='ntable must be >= 0 and <= 99'
       goto 9999
      endif
      alpha=idre_param(10)
      beta=idre_param(11)
      rb=idre_param(12)
c      rbb=rb**(beta-alpha)
      rbb=DBLE(rb)**(-alpha+1.d0)
      if(ne.lt.1)then
       ierr=1
       errortxt='ne must be larger than 0'
       goto 9999
      endif
c initialize some variables:
      do i=1,ne
       fc(i)=0.d0
      enddo
c zzshift - multiplication factor for gfac from zshift
      zzshift=1.0/(1.0+zshift)
      zzshift2=zzshift**2
c marginally stable orbit in the Kerr geometry:
      pom1=(1.+am)**(0.33333333)
      pom2=(1.-am)**(0.33333333)
      pom3=(1.-am2)**(0.33333333)
      pom=1.+pom3*(pom1+pom2)
      pom1=sqrt(3.*am2+pom**2)
      rms=3.+pom1-sqrt((3.-pom)*(3.+pom+2.*pom1))
c let's check whether rin and rout are reasonable
      if(ms.ne.0.and.rin.lt.(rms-r_plus))rin=rms-r_plus
      if(rin.ge.rout)then
       do i=1,ne
        photar(i)=0.
       enddo
       return
      endif
c      write(text,*)'r_horizon = ',r_plus
c      call xwrite(text,5)
c      write(text,*)'r_ms = ',rms
c      call xwrite(text,5)
      write(kyrh,*)r_plus
      call fpmstr(pkyrh,kyrh)
      write(kyrin,*)rin+r_plus
      call fpmstr(pkyrin,kyrin)
      write(kyrms,*)rms
      call fpmstr(pkyrms,kyrms)
c initialize data tables dc() (interpolated for given am, thetaO):
c we read the tables if we haven't yet or if we have changed them or am,
c thetaO or if we have changed the KYDIR director...
c so that we can get new dc()
      if(ntable_old.eq.-1.or.kydir.ne.fgmstr(pname).or.
     $   ntable_old.ne.ntable.or.am_old.ne.am.or.thetaO_old.ne.thetaO)
     $ then
       call idre_ini(dc,r_vec,nr,am,thetaO,ntable,cmodel,
     $               edivision,min,max,nbins)
       if(ntable.eq.-1)then
        do i=1,ne
         photar(i)=0.
        enddo
        return
       endif
c      We allocate memory for emission only if we read the tables for the first
c      time or if we have changed number of bins (in case the tables have
c      changed)
       if(ntable_old.eq.-1.or.nbins.ne.nbins_old)then
c       free memory from emission if we have already read the tables
c       and number of bins has changed (nbins)
        if(ntable_old.ne.-1.and.nbins.ne.nbins_old)then
         ierr=0
         errortxt='failed to free memory from emission...'
         deallocate(emission, stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Allocate memory for emission...
        ierr=0
        errortxt='failed to allocate memory for emission...'
        allocate(emission(nbins), stat=ierr)
        if(ierr.ne.0)goto 9999
       endif
       am_old=am
       thetaO_old=thetaO
       ntable_old=ntable
       nbins_old=nbins
       kydir=fgmstr(pname)
      endif
c now let's check whether rin and rout fall in the range covered by data tables
      if(rin.lt.r_vec(1))then
       ierr=1
       errortxt='r out of range covered by data tables'
       goto 9999
      endif
      if(rout.gt.r_vec(nr))then
       ierr=1
       errortxt='r out of range covered by data tables'
       goto 9999
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Let's write the functions to the file for TESTing purposes...
C      open(51,file='test2.dat',status='unknown')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c let's clear emission
      emission=0.d0
c given r, find corresponding indices in r_vec:
      ir0=1
      irn=nr
      do i=2,nr
       if(rin.ge.(r_vec(i)+r_vec(i-1))/2.)ir0=i
       if(rout.le.(r_vec(nr-i+2)+r_vec(nr-i+1))/2.)irn=nr-i+1
      enddo
c let's integrate...loop over the disc;
c rnow is polar coord. - radius in units of GM/c^2:
      do ir=ir0,irn
c Let's multiply the emission by factor r^(-alpha+1)
c +1 because dS=dr*dphi (i.e. without r)
       if((r_vec(ir)+r_plus).le.rb)then
        ralpha=DBLE(r_vec(ir)+r_plus)**(-alpha+1.d0)
       else
c        ralpha=rbb*(memr(r_vec+ir-1)+r_plus)**(-beta+1)
        ralpha=rbb*DBLE((r_vec(ir)+r_plus)/rb)**(-beta+1.d0)
       endif
       if(ir0.eq.irn)then
        drad=rout-rin
       else
        if(ir.eq.ir0)drad=(r_vec(ir0+1)+r_vec(ir0))/2.-rin
        if(ir.eq.irn)drad=rout-(r_vec(irn)+r_vec(irn-1))/2.
        if(ir.gt.ir0.and.ir.lt.irn)drad=
     $                           (r_vec(ir+1)-r_vec(ir-1))/2.
       endif
       do ibin=1,nbins
        emission(ibin) = emission(ibin) +
     $   dc(ir+nr*(ibin-1))*ralpha*drad
       enddo
      enddo
c Let's multiply the emission by zzshift^2 (due to overall redshift)
      do ibin=1,nbins
       emission(ibin) = emission(ibin)*zzshift2
      enddo
c Let's rebin the flux
      do ie=1,ne
       do ie_loc=1,ne_loc
        Enorm=(ear_loc(ie_loc)+ear_loc(ie_loc-1))/2.
        gmin=ear(ie-1)/zzshift/Enorm
        gmax=ear(ie)/zzshift/Enorm
c      let's find lower and upper index of local energy -> energy intervals
c      with lower edge between this lower and upper value are shifted and fall
c      into the energy bin ear(ie)-ear(ie-1) (at least partially)
c      we use real variables xlow, xhigh for these indices so that they don't
c      overflow for very small energy bins
        if(edivision.eq.0)then
         xlow=(gmin-min)/(max-min)*nbins
         xhigh=(gmax-min)/(max-min)*nbins
        else
         if(edivision.eq.1)then
          xlow=(alog10(gmin)-alog10(min))/(alog10(max)-alog10(min))*
     $          nbins
          xhigh=(alog10(gmax)-alog10(min))/
     $          (alog10(max)-alog10(min))*nbins
         endif
        endif
c      if the lower index is smaller than 0 => the whole first
c      interval of local energies falls into this bin; we set lower
c      index to 1 (real value to 0, we round it to 1 later on)
c      the upper energy must be larger than ear(ie-1) in this case...(otherwise
c      no part of energy bin falls into this bin)
        if((xlow.lt.0.).and.(xhigh.ge.0.))xlow=0.
c      if the upper index is larger than nbins => the last interval of
c      energies also falls into this bin; we set upper index to
c      nbins (real value to (nbins-1), we round it to nbins later on)
c      the lower energy must be smaller than ear(ie) in this case...(otherwise
c      no part of energy bin falls into this bin)
        if((xhigh.ge.nbins).and.(xlow.lt.nbins))xhigh=nbins-1.
c      otherwise no energy bin falls into this bin:
        if((xlow.lt.0.).or.(xhigh.ge.nbins))goto 300
c      let's change real values of indices to integer values
        ilow=nint(xlow+0.5)
        ihigh=nint(xhigh+0.5)
        do ibin=ilow,ihigh
         if(edivision.eq.0)then
          gmin_loc=min+(max-min)*(ibin-1)/nbins
          gmax_loc=min+(max-min)*ibin/nbins
         endif
         if(edivision.eq.1)then
          gmin_loc=min*(max/min)**(real(ibin-1)/nbins)
          gmax_loc=min*(max/min)**(real(ibin)/nbins)
         endif
         dfc=0.d0
c       the whole interval of local energies gmax_loc - gmin_loc
c       falls into the bin gmax - gmin
         if((gmin_loc.ge.gmin).and.(gmax_loc.le.gmax))then
          dfc=1.d0
         else
c	       only the higher part of the local energy interval falls into the
c        bin gmax - gmin => we have to consider only part of bin_loc
          if((gmin_loc.lt.gmin).and.(gmax_loc.le.gmax).and.
     $       (gmax_loc.gt.gmin))
     $     dfc=(gmax_loc-gmin)/(gmax_loc-gmin_loc)
c	       only the lower part of the local energy interval falls into the
c        bin gmax - gmin => we have to consider only part of bin_loc
          if((gmin_loc.ge.gmin).and.(gmin_loc.lt.gmax).and.
     $       (gmax_loc.gt.gmax))
     $     dfc=(gmax-gmin_loc)/(gmax_loc-gmin_loc)
c	       only the middle part of the local energy interval falls into the
c        bin gmax - gmin => we have to consider only part of bin_loc
          if((gmin_loc.lt.gmin).and.(gmax_loc.gt.gmax))
     $      dfc=(gmax-gmin)/(gmax_loc-gmin_loc)
         endif
         fc(ie)=fc(ie)+dfc*emission(ibin)*far_loc(ie_loc)
     $          *(ear_loc(ie_loc)-ear_loc(ie_loc-1))
        enddo
300     continue
       enddo
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Let's write the functions to the file for TESTing purposes...
C      close(51)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c let's smooth the spectra a little bit if it is required
      if(smooth.eq.0)then
       do i=1,ne
         photar1(i)=fc(i)
       enddo
      else
       photar1(1)=fc(1)
       photar1(ne)=fc(ne)
       do i=2,ne-1
        photar1(i)=(fc(i-1)/(ear(i-1)-ear(i-2))
     $          +2.0*fc(i)/(ear(i)-ear(i-1))
     $          +fc(i+1)/(ear(i+1)-ear(i)))/4.0*(ear(i)-ear(i-1))
       enddo
      endif
c let's normalize the spectrum if it is required (i.e. normal >= 0.)
      if(normal.ge.0.)then
c      if normal is 0 or if normal is out of energy range ear(0)-ear(ne)
c      we normalize the spectrum to total flux equal to unity
       if((ear(0).gt.normal.or.ear(ne).lt.normal).or.normal.eq.0.)then
        photar_int=0.d0
        do i=1,ne
         photar_int=photar_int+photar1(i)
        enddo
       else
c       otherwise we normalize to the maximum at normal keV (i.e. the photon
c       flux density per keV is 1 at normal keV)
c       let's find the energy bin where normal falls
        do i=1,ne
         if(ear(i-1).le.normal.and.ear(i).ge.normal)goto 100
        enddo
100     continue
c       let's calculate flux at normal keV - we interpolate two neighbouring
c       values (that are at the middle of energy bins)
        if(((ear(i)+ear(i-1))/2..ge.normal.and.i.gt.1).or.
     $    ((ear(i)+ear(i-1))/2..lt.normal.and.i.eq.ne))i=i-1
        photar_int=photar1(i)/(ear(i)-ear(i-1))+
     $   (2.*normal-(ear(i)+ear(i-1)))/(ear(i+1)-ear(i-1))*
     $   (photar1(i+1)/(ear(i+1)-ear(i))-photar1(i)/(ear(i)-ear(i-1)))
c       if the flux at normal keV is too small (zero) we normalize the spectrum
c       to total flux equal to unity
        if(photar_int.lt.1e-30)then
         photar_int=0.d0
         do i=1,ne
          photar_int=photar_int+photar1(i)
         enddo
        endif
       endif
       if(photar_int.gt.1e-30)then
        do i=1,ne
         photar(i)=photar1(i)/photar_int
        enddo
       endif
      else
       do i=1,ne
        photar(i)=photar1(i)
       enddo
      endif
9999  continue
      if(ierr.ne.0)then
       errortext='idre: '//errortxt
       call xwrite(errortext,5)
       STOP
      endif
      return
      end
c==============================================================================
      subroutine idre_ini(dc,r_vec,nr,am,thetaO,ntable,cmodel,
     $                    edivision,min,max,nbins_final)
c read the data tables from the fits file and
c initialize data table dc(:,:,:) for given ntable and parameters am and thetaO
c-----------------------variables declarations---------------------------------
      implicit real(a-h,o-z)
      
      integer, intent(inout) :: nr
      integer, intent(out) :: edivision
      real, intent(in) :: am
      real, intent(in) :: thetaO
      integer, intent(inout) :: ntable
      real, intent(out) :: min
      real, intent(out) :: max
      integer, intent(out) :: nbins_final
c maximum number of different values for am, thetaO in each data set:
      character(len=32), intent(in) :: cmodel
c     pointers to arrays from idre...
      real, dimension(:), allocatable, intent(inout) :: dc(:)
      real, dimension(:), allocatable, intent(inout) :: r_vec(:)
      
      character*(255) KBHfl,kydir
      character*(74) errortxt
      character*(80) errortext
      integer nhorizon,nincl,ierr,nbins
      real gmin,gmax,thetaOO
      integer ntable_old
      integer nhorizon_old,nincl_old,nr_old
c     pointers to tmp arrays...
      real, dimension(:), allocatable :: amvec(:)
      real, dimension(:), allocatable :: rmvec(:)
      real, dimension(:), allocatable :: thvec(:)
      real, dimension(:), allocatable :: gminimum(:)
      real, dimension(:), allocatable :: gmaximum(:)
      real, dimension(:), allocatable :: emission0(:)
      real, dimension(:), allocatable :: dc0(:)
      real(8), dimension(:), allocatable :: dc1(:)
c these are needed to work with the fits file...
      integer status,unit,blocksize,hdutype
      integer frow,felem,nelems,colnum
      real nulle
      double precision tscal,tzero
      character*(1) datatype
      character*(8) ttype,tunit,snull,tdisp
      logical anynull
      integer lenn
      character*(128) fgdatd,pname,fgmstr
      external fgdatd,fgmstr

      save

      data pname/'KYDIR'/
      data kydir/''/
      data ntable_old/-1/
      data nhorizon_old,nincl_old,nr_old/3*0/
      data pi/3.1415926535897932/
      data pi2/6.2831853071795865/

c--------------subroutine execution starts here-------------------------------
c let us initialize some of the parameters
       edivision=0
       nbins_final=500
       min=0.001
       max=1.66
       gdelta_final=(max-min)/nbins_final
c reading table files including vector r_vec -- execute on the
c first pass only or when we choose different tables
      if(kydir.ne.fgmstr(pname).and.ntable_old.ne.-1)ntable_old=-2
      if(ntable_old.eq.-1.or.ntable.ne.ntable_old)then
       write(*,*)'idre: initializing data tables, please wait...'
       write(*,*)'      Ref.: Dovciak M., Karas V. & Yaqoob T.'
       write(*,*)'            ApJS July 2004, Volume 153, Issue 1, '//
     $           'pp. 205-221'
       write(*,*)'      -------------------------------------------'//
     $           '-----------'

c  The status parameter must always be initialized.
       status=0
c  Get an unused Logical Unit Number to use to open the FITS file.
       call ftgiou(unit,status)
       if(status.ne.0)goto 999
c  Open the FITS file for readonly access
c  - if set try KYDIR directory, otherwise look in the working directory
c    or in the xspec directory where tables are usually stored...
       ierr=0
       KBHfl = fgmstr(pname)
       lenn = lenact(KBHfl)
       if(lenn.ne.0)then
        KBHfl = KBHfl(1:lenn)//'/'//cmodel
        call ftopen(unit,KBHfl,0,blocksize,status)
        if(status.eq.0)kydir=fgmstr(pname)
       else
        call ftopen(unit,cmodel,0,blocksize,status)
        if(status.ne.0)then
         status=0
         KBHfl = fgdatd()
         lenn = lenact(KBHfl)
         KBHfl = KBHfl(1:lenn)//cmodel
         call ftopen(unit,KBHfl,0,blocksize,status)
        endif
       endif
       if(status.ne.0)then
        call printerr(status)
c       let's free the unit number
        call ftfiou(unit, status)
        ntable=-1
        errortxt=
     $   'idre: set the KYDIR to the directory with the KY tables'
        call xwrite(errortxt,5)
        return
       endif
c  Let's read tables (binary tables => hdutype=2)
       hdutype=2
       colnum=1
       frow=1
       felem=1
       nulle=-1.
c Move to the extension 'r_vector' and read r_vec values
       call ftmrhd(unit,1,hdutype,status)
       if(ntable_old.eq.-1.or.ntable.ne.ntable_old)then
c       read number of r_vector values
        call ftgnrw(unit,nr,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        write(*,*)
c        write(*,'("Number of r_vec values: ",i3)')nr
ccccccccccccccccccccccccccccccc
c       We allocate memory for arrays only if we read the tables for the first
c       time or if we have changed the dimensions of the arrays (in case the
c       tables have changed)
        if(ntable_old.eq.-1.or.nr.ne.nr_old)then
c        Firstly we have to free allocated memory for these arrays if we have
c        already allocated it and the dimensions of the arrays have changed
         if((ntable_old.ne.-1).and.(nr.ne.nr_old))then
c        Free memory from tmp arrays...
          ierr=0
          errortxt='failed to free memory from tmp arrays...'
          deallocate(r_vec, stat=ierr)
          if(ierr.ne.0)goto 9999
         endif
c        Allocate memory for r_vec...
         ierr=0
         errortxt='failed to allocate memory for tmp arrays...'
         allocate(r_vec(nr), stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Read the data in the 'r_vector' table
        nelems=nr
c       FTGCVE reads the VALUES from the first column.
        call ftgcve(unit,colnum,frow,felem,nelems,nulle,r_vec,
     $              anynull,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        do i=1,nr
c         write(*,'(i4,3x,f11.6)')i,memr(r_vec+i-1)
c        enddo
cccccccccccccccccccccccccccccccc
       endif
c Move to the extension 'r_horizon' and read r_horizon values
       call ftmrhd(unit,1,hdutype,status)
       if(ntable_old.eq.-1.or.ntable.ne.ntable_old)then
c       read number of r_horizon values
        call ftgnrw(unit,nhorizon,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        write(*,'("Number of r_horizon values: ",i3)')nhorizon
ccccccccccccccccccccccccccccccc
c       We allocate memory for arrays only if we read the tables for the first
c       time or if we have changed the dimensions of the arrays (in case the
c       tables have changed)
        if(ntable_old.eq.-1.or.nhorizon.ne.nhorizon_old)then
c        Firstly we have to free allocated memory for these arrays if we have
c        already allocated it and the dimensions of the arrays have changed
         if((ntable_old.ne.-1).and.(nhorizon.ne.nhorizon_old))then
c        Free memory from tmp arrays...
          ierr=0
          errortxt='failed to free memory from tmp arrays...'
          deallocate(amvec, stat=ierr)
          if(ierr.ne.0)goto 9999
          deallocate(rmvec, stat=ierr)
          if(ierr.ne.0)goto 9999
         endif
c        Allocate memory for amvec and rmvec...
         ierr=0
         errortxt='failed to allocate memory for tmp arrays...'
         allocate(amvec(nhorizon), stat=ierr)
         if(ierr.ne.0)goto 9999
         allocate(rmvec(nhorizon), stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Read the data in the 'r_horizon' table
        nelems=nhorizon
c       FTGCVE reads the VALUES from the first column.
        call ftgcve(unit,colnum,frow,felem,nelems,nulle,rmvec(1),
     $              anynull,status)
        if(status.ne.0)goto 999
c       let's calculate am values for r_horizon values
        do i=1,nhorizon
         amvec(i)=sqrt(1.0-(rmvec(i)-1.0)*(rmvec(i)-1.0))
ccccccccccccccccccccccccccccccc
c         write(*,'(i4,3x,f4.2,3x,f8.6)')i,memr(rmvec+i-1),
c          memr(amvec+i-1)
ccccccccccccccccccccccccccccccc
        enddo
c       check if current value of am is in the range covered by data tables
        if((am.lt.amvec(nhorizon)).or.(am.gt.amvec(1)))then
         ierr=1
         errortxt='a/M out of table range'
         goto 9999
        endif
       endif
c Move to the extension 'inclination' and read inclination values
       call ftmrhd(unit,1,hdutype,status)
       if(ntable_old.eq.-1.or.ntable.ne.ntable_old)then
c       read number of inclination values
        call ftgnrw(unit,nincl,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        write(*,*)
c        write(*,'("Number of inclination values: ",i3)')nincl
ccccccccccccccccccccccccccccccc
c       We allocate memory for arrays only if we read the tables for the first
c       time or if we have changed the dimensions of the arrays (in case the
c       tables have changed)
        if(ntable_old.eq.-1.or.nincl.ne.nincl_old)then
c        Firstly we have to free allocated memory for these arrays if we have
c        already allocated it and the dimensions of the arrays have changed
         if((ntable_old.ne.-1).and.(nincl.ne.nincl_old))then
c        Free memory from tmp arrays...
          ierr=0
          errortxt='failed to free memory from tmp arrays...'
          deallocate(thvec, stat=ierr)
          if(ierr.ne.0)goto 9999
         endif
c        Allocate memory for thvec...
         ierr=0
         errortxt='failed to allocate memory for tmp arrays...'
         allocate(thvec(nincl), stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Read the data in the 'inclination' table
        nelems=nincl
c       FTGCVE reads the VALUES from the first column.
        call ftgcve(unit,colnum,frow,felem,nelems,nulle,thvec(1),
     $              anynull,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        do i=1,nincl
c         write(*,'(i4,3x,f4.1)')i,memr(thvec+i-1)
c        enddo
cccccccccccccccccccccccccccccccc
c       check if current value of thetaO is in the range covered by data tables
c       if thetaO is "almost" 0 and fits tables are defined for higher values
c       then thetaO but from close to 0 we change thetaO to the lowest value
c       in the fits tables...(for convenience we put 0.1 equal to 0.0 ...)
        thetaOO=thetaO
        if(((thvec(1)-thetaO).le.0.1).and.(thetaO.le.0.1))
     $   thetaOO=thvec(1)
        if((thetaOO.lt.thvec(1)).or.
     $   (thetaOO.gt.thvec(nincl)))then
         ierr=1
         errortxt='theta_o out of table range'
         goto 9999
        endif
       endif
c Let's read the transfer functions...

c      Let's find out how many energy bins there are in the tables (nbins)
        call ftmrhd(unit,1,hdutype,status)
        call ftgbcl(unit,3,ttype,tunit,datatype,nbins,tscal,tzero,
     $   snull,tdisp,status)
        if(status.ne.0)goto 999
ccccccccccccccccccccccccccccccc
c        write(*,'("Number of r_horizon values: ",i3)')nhorizon
ccccccccccccccccccccccccccccccc
c      We allocate memory for arrays only if we read the tables for the first
c      time or if we have changed the dimensions of the arrays (in case the
c      tables have changed)
       if(ntable_old.eq.-1.or.nhorizon.ne.nhorizon_old.or.
     $   nincl.ne.nincl_old.or.nr.ne.nr_old)then
c       Firstly we have to free allocated memory for these arrays if we have
c       already allocated it and the dimensions of the arrays have changed
        if((ntable_old.ne.-1).and.(nhorizon.ne.nhorizon_old.or.
     $    nincl.ne.nincl_old.or.nr.ne.nr_old))then
c        Free memory from tmp arrays...
         ierr=0
         errortxt='failed to free memory from tmp arrays...'
         deallocate(gminimum, stat=ierr)
         if(ierr.ne.0)goto 9999
         deallocate(gmaximum, stat=ierr)
         if(ierr.ne.0)goto 9999
         deallocate(emission0, stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Allocate memory for gminimum, gmaximum and emission...
        ierr=0
        errortxt='failed to allocate memory for tmp arrays...'
        allocate(gminimum(nr*nhorizon*nincl), stat=ierr)
        if(ierr.ne.0)goto 9999
        allocate(gmaximum(nr*nhorizon*nincl), stat=ierr)
        if(ierr.ne.0)goto 9999
        allocate(emission0(nr*nbins*nhorizon*nincl), stat=ierr)
        if(ierr.ne.0)goto 9999
       endif
c  Let's read gminimum, gmaximum and emission...
       nelems=nr*nbins
       do ihorizon=1,nhorizon
        do iincl=1,nincl
         if(ihorizon.ne.1.or.iincl.ne.1)then
          call ftmrhd(unit,1,hdutype,status)
         endif
         if(status.ne.0)goto 999
c        to read the file only once we have to read in blocks (all columns
c        from the extension are put to buffer together)
c        let's find out how many rows are going to be read into the buffer
         call ftgrsz(unit,nrow,status)
         if(status.ne.0)goto 999
         nelements=nrow
         do irow=1,nr,nrow
          
c         the last block to read may be smaller:
          if((nr-irow+1).lt.nrow)nelements=(nr-irow+1)
     
          call ftgcve(unit,1,irow,1,nelements,nulle,
     $     gminimum((irow-1)+nr*(ihorizon-1)
     $     +nr*nhorizon*(iincl-1)+1),anynull,status)
     
          call ftgcve(unit,2,irow,1,nelements,nulle,
     $     gmaximum((irow-1)+nr*(ihorizon-1)
     $     +nr*nhorizon*(iincl-1)+1),anynull,status)
          if(status.ne.0)goto 999
          
          call ftgcve(unit,3,irow,1,nelements*nbins,nulle,
     $     emission0((irow-1)*nbins+nr*nbins*(ihorizon-1)
     $     +nr*nbins*nhorizon*(iincl-1)+1),anynull,status)
          if(status.ne.0)goto 999
          
         enddo
        enddo
       enddo
ccccccccccccccccccccccccccccccc
c      open(66,file='emission.txt',status='unknown')
c      do i=1,nr
c       write(66,'(300g15.7)')(memr(emission0+j-1+nbins*(i-1)+
c    $	 nr*nbins*2+nr*nbins*nhorizon*8),j=1,nbins)
c      enddo
c      close(66)
cccccccccccccccccccccccccccccccc
c The FITS file must always be closed before exiting the program.
c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
       call ftclos(unit, status)
       call ftfiou(unit, status)
c Check for any error, and if so print out error messages.
c The PRINTERR subroutine is listed near the end of this file.
999    if (status .ne. 0)then
        call printerr(status)
        STOP
       endif
c      Let's create array dc() for interpolated transfer functions:
c      We allocate memory for arrays only if we read the tables for the first
c      time or if we have changed the dimensions of the arrays (in case the
c      tables have changed)
       if(ntable_old.eq.-1.or.nr.ne.nr_old)then
c       Firstly we have to free allocated memory for these arrays if we have
c       already allocated it and the dimensions of the arrays have changed
        if((ntable_old.ne.-1).and.(nr.ne.nr_old))then
c        Free memory from tmp arrays...
         ierr=0
         errortxt='failed to free memory from tmp arrays...'
         deallocate(dc, stat=ierr)
         if(ierr.ne.0)goto 9999
        endif
c       Allocate memory for dc...
        ierr=0
        errortxt='failed to allocate memory for tmp arrays...'
        allocate(dc(nr*nbins_final), stat=ierr)
        if(ierr.ne.0)goto 9999
       endif
       nhorizon_old=nhorizon
       nincl_old=nincl
       nr_old=nr
       ntable_old=ntable
       write(*,*)'       ...initializing finished'
      endif
c end of reading table fits file................................................
c     check if current value of thetaO is in the range covered by data tables
c     if thetaO is "almost" 0 and fits tables are defined for higher values
c     then thetaO but from close to 0 we change thetaO to the lowest value
c     in the fits tables...(for convenience we put 0.1 equal to 0.0 ...)
      thetaOO=thetaO
      if(((thvec(1)-thetaO).le.0.1).and.(thetaO.le.0.1))
     $  thetaOO=thvec(1)
      iam0=0
      ith0=0
c given am and thetaO, look up the stored tables that are needed now:
      do i=1,nhorizon
       if(am.le.amvec(i))then
        iam0=i
       else
        goto 1003
       endif
      enddo
1003  continue
      do i=1,nincl
       if(thetaOO.ge.thvec(i))then
        ith0=i
       else
        goto 1004
       endif
      enddo
1004  continue
      if((iam0.eq.nhorizon).and.(am.eq.amvec(nhorizon)))
     $ iam0=iam0-1
      if((ith0.eq.nincl).and.(thetaOO.eq.thvec(nincl)))
     $ ith0=ith0-1
c of all tables, four sets are needed now;
c 1) a_low/th0_low, 2) a_low/th0_high, 3) a_high/th0_low, 4) a_high/th0_high:
c interpolate data by bilinear interpolation to desired values of am and thetaO;
c result is in dc():
      rm=1.+sqrt(1.-am*am)
      ttmp=(rm-rmvec(iam0))/(rmvec(iam0+1)-rmvec(iam0))
      utmp=(thetaOO-thvec(ith0))/(thvec(ith0+1)-thvec(ith0))
      ttmp1=1.-ttmp
      utmp1=1.-utmp
      do j=1,nr
c gminimum
       y1=gminimum(j+nr*(iam0-1)+nr*nhorizon*(ith0-1))
       y2=gminimum(j+nr*iam0+nr*nhorizon*(ith0-1))
       y3=gminimum(j+nr*iam0+nr*nhorizon*ith0)
       y4=gminimum(j+nr*(iam0-1)+nr*nhorizon*ith0)
       gmin=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4)
c gmaximum
       y1=gmaximum(j+nr*(iam0-1)+nr*nhorizon*(ith0-1))
       y2=gmaximum(j+nr*iam0+nr*nhorizon*(ith0-1))
       y3=gmaximum(j+nr*iam0+nr*nhorizon*ith0)
       y4=gmaximum(j+nr*(iam0-1)+nr*nhorizon*ith0)
       gmax=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4)
c      Allocate memory for dc0...
       ierr=0
       errortxt='failed to allocate memory for tmp arrays...'
       allocate(dc0(nbins), stat=ierr)
       if(ierr.ne.0)goto 9999
       allocate(dc1(nbins_final), stat=ierr)
       if(ierr.ne.0)goto 9999
       do i=1,nbins
c emission
        y1=emission0(i+nbins*(j-1)
     $     +nr*nbins*(iam0-1)+nr*nbins*nhorizon*(ith0-1))
        if(y1.lt.0.)y1=0.
        y2=emission0(i+nbins*(j-1)
     $     +nr*nbins*iam0+nr*nbins*nhorizon*(ith0-1))
        if(y2.lt.0.)y2=0.
        y3=emission0(i+nbins*(j-1)
     $     +nr*nbins*iam0+nr*nbins*nhorizon*ith0)
        if(y3.lt.0.)y3=0.
        y4=emission0(i+nbins*(j-1)
     $     +nr*nbins*(iam0-1)+nr*nbins*nhorizon*ith0)
        if(y4.lt.0.)y4=0.
        dc0(i)=utmp1*(ttmp1*y1+ttmp*y2)+utmp*(ttmp*y3+ttmp1*y4)
       enddo
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c The following is made for gdivision=0. if you want to change to gdivision=1
c you have to make appropriate changes
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       dc1=0.d0
       if((gmin.ne.0.).and.(gmax.ne.0.))then
        ibin=0
        gdelta=(gmax-gmin)/nbins
        glow=gmin-0.5*gdelta
        ghigh=gmin+0.5*gdelta
        if(glow.lt.min)then
         ibin_final=1
        else
         ibin_final=int((glow-min)/gdelta_final)+1
        endif
        glow_final=min+gdelta_final*(ibin_final-1)
        ghigh_final=min+gdelta_final*ibin_final
        do while ((ibin_final.le.nbins_final).and.(ibin.le.nbins))
         if(ibin.eq.0)then
          flower=0.
         else
          flower=dc0(ibin)
         endif
         if(ibin.eq.nbins)then
          fupper=0.
         else
          fupper=dc0(ibin+1)
         endif
         if(ghigh.le.ghigh_final)then
          if(ghigh.gt.glow_final)then
           if(glow.lt.glow_final)then
            dc1(ibin_final)=dc1(ibin_final)+
     $       0.5*(fupper+(flower+
     $       (glow_final-glow)*(fupper-flower)/gdelta))*
     $       (ghigh-glow_final)
           else
            dc1(ibin_final)=dc1(ibin_final)+
     $       gdelta*(flower+fupper)/2.
           endif
          endif
          ibin=ibin+1
          glow=gmin+gdelta*0.5*(2.*ibin-1)
          ghigh=gmin+gdelta*0.5*(2.*ibin+1)
         else
          if(glow.lt.glow_final)then
           dc1(ibin_final)=dc1(ibin_final)+
     $      (flower+0.5*(ghigh_final+glow_final-2.*glow)*
     $      (fupper-flower)/gdelta)*
     $      (ghigh_final-glow_final)
          else
           if(glow.lt.ghigh_final)
     $      dc1(ibin_final)=dc1(ibin_final)+
     $       0.5*(flower+(flower+
     $       (ghigh_final-glow)*(fupper-flower)/gdelta))*
     $       (ghigh_final-glow)
          endif
          ibin_final=ibin_final+1
          glow_final=min+gdelta_final*(ibin_final-1)
          ghigh_final=min+gdelta_final*ibin_final
         endif
        enddo
       endif
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do i=1,nbins_final
         dc(j+nr*(i-1))=real(dc1(i))
       enddo
c      Free memory from tmp arrays...
       ierr=0
       errortxt='failed to free memory from tmp arrays...'
       deallocate(dc0, stat=ierr)
       if(ierr.ne.0)goto 9999
       deallocate(dc1, stat=ierr)
       if(ierr.ne.0)goto 9999
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Let's write the functions to the file for TESTing purposes...
C      open(52,file='test1.dat',status='unknown')
C      do i=1,nbins_final
C       write(52,'(501g14.6)')min+(i-1)*gdelta_final,
C     $  (memr(dc+j-1+nr*(i-1)),j=1,nr)
C      enddo
C      close(52)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
9999  continue
      if(ierr.ne.0)then
       errortext='idre: '//errortxt
       call xwrite(errortext,5)
       STOP
      endif
      return
      end
c==============================================================================
      subroutine printerr(status)

c  This subroutine prints out the descriptive text corresponding to the
c  error status value and prints out the contents of the internal
c  error message stack generated by FITSIO whenever an error occurs.

      integer status
      character*(30) errtext
      character*(80) errmessage,text

c  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return

c  The FTGERR subroutine returns a descriptive 30-character text string that
c  corresponds to the integer error status number.  A complete list of all
c  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
c      print *,'FITSIO Error Status =',status,': ',errtext
      write(text,*)'idre: FITSIO Error Status =',status,': ',errtext
      call xwrite(text,5)

c  FITSIO usually generates an internal stack of error messages whenever
c  an error occurs.  These messages provide much more information on the
c  cause of the problem than can be provided by the single integer error
c  status value.  The FTGMSG subroutine retrieves the oldest message from
c  the stack and shifts any remaining messages on the stack down one
c  position.  FTGMSG is called repeatedly until a blank message is
c  returned, which indicates that the stack is empty.  Each error message
c  may be up to 80 characters in length.  Another subroutine, called
c  FTCMSG, is available to simply clear the whole error message stack in
c  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
c          print *,errmessage
          call xwrite(errmessage,5)
          call ftgmsg(errmessage)
      end do
      return
      end
