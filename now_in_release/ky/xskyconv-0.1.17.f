c kyc1onv - convolution general relativistic model - non-axisymmetric version
c Fortran77 model subroutine for XSPEC
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
c This convolution subroutine takes input photar() array as a definition of the
c local flux across the accretion disc around a black hole and adds a power-law
c radial dependence and limb darkening/brightening law to it. The output is 
c total spectrum of an accretion disc. All relativistic effects are taken into 
c account. This model calls subroutine idre() for integrating local emission 
c over the disc and uses the fits file 'KBHlineNN.fits' defining the transfer 
c function needed for integration. For details on idre() and the fits file see 
c the subroutine idre() in xsidre.f.
c
c There are several restrictions that arise from the fact that we use existing 
c XSPEC models for definition of the local flux:
c - by local XSPEC models only the energy dependence of the photon flux can be 
c   defined, 
c - only a certain type of radial dependence of the local photon flux can be 
c   imposed - we have chosen to use a power-law radial dependence, 
c - there is no intrinsic azimuthal dependence of the local photon flux, the
c   only azimuthal dependence comes through limb darkening/brightening law 
c   (emission angle depends on azimuth)
c - local flux can highly depend on the energy resolution, i.e. on the energy 
c   binning used, if the energy resolution is not high enough. This is because 
c   the flux is defined in the centre of each bin. A large number of bins is 
c   needed for highly varying local flux. 
c
c For emissivities that cannot be defined by existing XSPEC models, or where the
c limitations mentioned above are too restrictive, one has to add a new 
c user-defined model to XSPEC (by adding a new subroutine to XSPEC). This method
c is more flexible and faster than when using this convolution model, and hence
c it is recommended even for cases when this model could be used. In any new 
c model for XSPEC one can use the common ray-tracing driver for relativistic 
c smearing of the local emission: ide() for non-axisymmetric models and idre()
c or idre2() for axisymmetric ones. 
c
c par1  ... a/M     - black hole angular momentum (0 <= a/M <= 1)
c par2  ... theta_o - observer inclination in degrees (0-pole, 90-disc)
c par3  ... rin     - inner edge of non-zero disc emissivity in GM/c^2
c par4  ... ms  - 0 - we integrate from inner edge rin-rh
c                 1 - if the inner edge of the disc is below marginally stable 
c                     orbit then we integrate emission above MSO only
c par5  ... rout    - outer edge of non-zero disc emissivity in GM/c^2
c par6  ... alpha - inner power-law index for radial dependence of emissivity,
c                   scales as r^(-alpha) below boundary radius rb
c par7  ... beta  - outer power-law index for radial dependence of emissivity,
c                   scales as rb^(beta-alpha) * r^(-beta) above boundary 
c                   radius rb
c par8  ... rb    - boundary radius between inner and outer power-law radial 
c                   dependence of emissivity (in units of GM/c^2)
c par9  ... zshift - overall Doppler shift
c par10 ... ntable - defines fits file with tables (0 <= ntable <= 99)
c                    0 -> for isotropic emission (flux ~ 1)
c                    1 -> for Laor's limb darkening (flux ~ (1+2.06*mu))
c                    2 -> for Haardt's limb brightening (flux ~ ln (1+1/mu))
c par11 ... ne_loc - number of grid points in local energy (energy resolution of
c                    local flux, the grid is equidistant in logarithmic scale)
c par12 ... normal - how to normalize the spectra
c                    = 0 - normalization to unity (total flux=1, e.g. for line)
c                    > 0 - normalization so that the flux=1 at 'normal' keV 
c                          (e.g. for continuum)
c                    < 0 - spectrum is not normalized
c
c===============================================================================
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Uncomment next lines if you want a stand-alone program outside XSPEC
Cc you will need library with idre() subroutine
C      program ky_main
C      implicit real(a-h,o-z)
C      parameter(ne=200,e_min=0.01,e_max=15.,nparam=12)
C      real ear(0:ne),param(nparam),photar(ne),photer(ne)
Cc parameters:     a/M    thetaO  rin  ms  rout alpha beta rb
C      data param/0.9982,   70.,   1., 0., 999.,  3.,  3., 400.,
C     $    0.,   0.,  100.,  1./
Cc       zshift limb ne_loc normal
Cc-----------------------------------------------------------------------------
C      ifl=1
C      do i=0,ne
Cc       ear(i)=e_min+i*(e_max-e_min)/ne
C       ear(i)=e_min*(e_max/e_min)**(float(i)/ne)
C      enddo
C      do i=1,ne
C       photar(i)=2./(ear(i)+ear(i-1))
C      enddo
C      call kyconv(ear,ne,param,ifl,photar,photer)
C      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c===============================================================================
      subroutine kyconv(ear,ne,param,ifl,photar,photer)
      implicit real(a-h,o-z)
      real ear(0:ne),param(*),photar(ne),photer(ne)
      real idre_param(12)
      integer ne_loc,ierr,ne_loc_old
      real, dimension(:), allocatable :: ear_local(:)
      real, dimension(:), allocatable :: far_local(:)
      character*(80) errortxt
      character*(32) cmodel

      save ne_loc_old,ear_local,far_local

      data ne_loc_old/-1/

c Let's initialize parameters for subroutine idre()
      cmodel="KBHline"
c a/M - black hole angular momentum
      idre_param(1)=param(1)
c theta_o - observer inclination
      idre_param(2)=param(2)
c rin - inner edge of non-zero disc emissivity
      idre_param(3)=param(3)
c ms - whether to integrate from rin or rms
      idre_param(4)=param(4)
c rout - outer edge of non-zero disc emissivity
      idre_param(5)=param(5)
c alpha - index of inner radial powerlaw dependence
      idre_param(10)=param(6)
c beta - index of outer radial powerlaw dependence
      idre_param(11)=param(7)
c rb - boundary radius between inner and outer radial powerlaw dependence
      idre_param(12)=param(8)
c zshift - overall Doppler shift
      idre_param(8)=param(9)
c ntable - table model (defines fits file with tables)
      idre_param(9)=param(10)
      ntable=int(param(10))
c smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
      idre_param(6)=0.
c normal - how to normalize the final spectrum
      idre_param(7)=param(12)
c number of points in local energy (resolution of the local flux)
      ne_loc=int(param(11))
      if(ne_loc_old.eq.-1.or.ne_loc.ne.ne_loc_old)then
       if(ne_loc_old.ne.-1.and.ne_loc.ne.ne_loc_old)then
        ierr=0
        errortxt='Failed to free memory for tmp arrays.'
        deallocate(ear_local, stat=ierr)
        if(ierr.ne.0)goto 999
        deallocate(far_local, stat=ierr)
        if(ierr.ne.0)goto 999
       endif
       ierr=0
       errortxt='Failed to allocate memory for tmp arrays.'
       allocate(ear_local(0:ne_loc), stat=ierr)
       if(ierr.ne.0)goto 999
       allocate(far_local(ne_loc), stat=ierr)
       if(ierr.ne.0)goto 999
999    continue
       if(ierr.ne.0)then
        call xwrite(errortxt,5)
        STOP
       endif
       ne_loc_old=ne_loc
      endif
c     assigning ear_local
      do i=0,ne_loc
       ear_local(i)=(ear(0)+ear(1))/2.0
     $  *((ear(ne-1)+ear(ne))/(ear(0)+ear(1)))**(float(i)/(ne_loc-1))
      enddo
c     assigning far_local
      do i=1,ne_loc
       do j=1,ne-1
        if(((ear_local(i-1)+ear_local(i))/2..gt.
     $      (ear(j)+ear(j-1))/2.).and.
     $     ((ear_local(i-1)+ear_local(i))/2..le.
     $       (ear(j+1)+ear(j))/2.))then
         far_local(i)=photar(j)/(ear(j)-ear(j-1))+
     $    ((ear_local(i-1)+ear_local(i))-(ear(j)+ear(j-1)))/
     $     (ear(j+1)-ear(j-1))*
     $     (photar(j+1)/(ear(j+1)-ear(j))-photar(j)/(ear(j)-ear(j-1)))
         goto 1000
        endif
       enddo
c      if j-loop finished before assigning anything to far_local, do it now
       far_local(i)=photar(j)/(ear(j)-ear(j-1))
1000   continue
      enddo
c     leave the last value as-is
      do i=1,ne_loc-1
       far_local(i)=(far_local(i)+far_local(i+1))/2.
      enddo
c model number (defines the set of data tables to be used; 0<=ntabel<=99):
       if(ntable.le.9)then
        write(cmodel,'(a7,2i1,a5)')'KBHline',0,ntable,'.fits'
       else
        write(cmodel,'(a7,i2,a5)')'KBHline',ntable,'.fits'
       endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Let's write input parameters to a text file
C      open(52,file='kyconv.txt',status='unknown')
C      write(52,'(a10,1f12.6)')'a/M       ',param(1)
C      write(52,'(a10,1f12.6)')'theta_o   ',param(2)
C      write(52,'(a10,1f12.6)')'rin       ',param(3)
C      write(52,'(a10,1i12)')'ms        ',int(param(4))
C      write(52,'(a10,1f12.6)')'rout      ',param(5)
C      write(52,'(a10,1f12.6)')'alpha     ',param(6)
C      write(52,'(a10,1f12.6)')'beta      ',param(7)
C      write(52,'(a10,1f12.6)')'rb        ',param(8)
C      write(52,'(a10,1f12.6)')'zshift    ',param(9)
C      write(52,'(a10,1i12)')'ntable    ',int(param(10))
C      write(52,'(a10,1i12)')'ne_loc    ',int(param(11))
C      write(52,'(a10,1f12.6)')'normal    ',param(12)
C      write(52,'(a10,1i12)')'smooth    ',int(idre_param(6))
C      close(52)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call idre(ear,ne,photar,idre_param,cmodel,ne_loc,
     $          ear_local,far_local)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc final spectrum output -- write ear() and photar() into file:
C      open(61,file='kyc_photar.dat',status='unknown')
C      do i=1,ne
C       write(61,'(2g14.6)')0.5*(ear(i-1)+ear(i)),
C     $                     photar(i)/(ear(i)-ear(i-1))
C      enddo
C      close(61)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      return
      end
