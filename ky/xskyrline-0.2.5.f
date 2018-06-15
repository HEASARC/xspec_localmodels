c kyrline - relativistic line - axisymmetric version
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
c This subroutine takes local Gaussian line emission and gives total spectrum
c of an accretion disc taking into account all relativistic effects. It calls
c subroutine idre() for integrating local emission over the disc and uses the
c fits file 'KBHlineNN.fits' defining the transfer function needed for
c integration. For details on idre() and the fits file see the subroutine idre().
c in xsidre.f
c
c par1 ... a/M     - black hole angular momentum (0 <= a/M <= 1)
c par2 ... theta_o - observer inclination in degrees (0-pole, 90-disc)
c par3 ... rin     - inner edge of non-zero disc emissivity in GM/c^2
c par4 ... ms - 0 - we integrate from inner edge rin
c               1 - if the inner edge of the disc is below marginally stable
c                   orbit then we integrate emission above MSO only
c par5 ... rout   - outer edge of non-zero disc emissivity in GM/c^2
c par6 ... Erest  - rest energy of the line (keV)
c par7 ... alpha  - inner power-law index for radial dependency of emissivity,
c                   scales as r^(-alpha) below boundary radius rb
c par8 ... beta   - outer power-law index for radial dependency of emissivity,
c                   scales as rb^(beta-alpha) * r^(-beta) above boundary radius
c par9 ... rb  - boundary radius between inner and outer power-law radial
c                dependency of emissivity (in units of GM/c^2)
c par10 ... zshift - overall Doppler shift
c par11 ... limb   - defines fits file with tables (0 <= ntable <= 99)
c                   0 -> for isotropic emission (flux ~ 1)
c                   1 -> for Laor's limb darkening (flux ~ (1+2.06*mu))
c                   2 -> for Haardt's limb brightening (flux ~ ln (1+1/mu))
c
c===============================================================================
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc Uncomment next lines if you want a stand-alone program outside XSPEC
Cc you will need library with idre() subroutine
C      program ky_main
C      implicit real(a-h,o-z)
C      parameter(ne=200,e_min=0.,e_max=10.,nparam=11)
C      real ear(0:ne),param(nparam),photar(ne),photer(ne)
Cc parameters:    a/M  thetaO rin  ms rout erest alpha beta rb
C      data param/0.9982, 70., 1., 1., 400,  6.4,   3.,  3., 400.,
C     $  0.,    0./
Cc     zshift limb
C      ifl=1
C      do i=0,ne
C       ear(i)=e_min+i*(e_max-e_min)/ne
Cc       ear(i)=e_min*(e_max/e_min)**(float(i)/ne)
C      enddo
C      call kyrline(ear,ne,param,ifl,photar,photer)
C      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c===============================================================================
      subroutine kyrline(ear,ne,param,ifl,photar,photer)
      implicit real(a-h,o-z)
      real ear(0:ne),param(*),photar(ne),photer(ne)
      parameter(ne_loc=1)
      real ear_loc(0:ne_loc),far_loc(ne_loc)
      real idre_param(12)
      character*(32) cmodel
      character*(80) errortext

c Let's initialize parameters for subroutine idre()
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
c Erest - normalization of local energy
      ear_loc(0)=param(6)-1e-3
      ear_loc(1)=param(6)+1e-3
      far_loc(1)=1./(ear_loc(1)-ear_loc(0))
      if(param(6).le.1e-1)then
       errortext='kyrline: Erest has to be larger than 0.1'
       call xwrite(errortext,5)
       do i=1,ne
        photar(i)=0.
       enddo
       return
      endif
c alpha - index of inner radial powerlaw dependence
      idre_param(10)=param(7)
c beta  - index of outer radial powerlaw dependence
      idre_param(11)=param(8)
c rb - boundary radius between inner and outer radial powerlaw
c      dependence
      idre_param(12)=param(9)
c zshift - overall Doppler shift
      idre_param(8)=param(10)
c ntable - table model (defines fits file with tables)
      idre_param(9)=param(11)
      ntable=int(param(11))
      if(ntable.lt.0.or.ntable.gt.99)then
       errortext='kyrline: ntable must be >= 0 and <= 99'
       call xwrite(errortext,5)
       do i=1,ne
        photar(i)=0.
       enddo
       return
      endif
c model number (defines the set of data tables to be used; 0<=ntabel<=99):
       if(ntable.le.9)then
        write(cmodel,'(a7,2i1,a5)')'KBHline',0,ntable,'.fits'
       else
        write(cmodel,'(a7,i2,a5)')'KBHline',ntable,'.fits'
       endif
c smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
      idre_param(6)=0.
c normal - how to normalize the final spectrum
      idre_param(7)=0.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc let's write the input parameters to a file
C      open(52,file='kyrline.txt',status='unknown')
C      write(52,'(a10,1f12.6)')'a/M       ',param(1)
C      write(52,'(a10,1f12.6)')'theta_o   ',param(2)
C      write(52,'(a10,1f12.6)')'rin       ',param(3)
C      write(52,'(a10,1i12)')'ms        ',int(param(4))
C      write(52,'(a10,1f12.6)')'rout      ',param(5)
C      write(52,'(a10,1f12.6)')'Erest     ',param(6)
C      write(52,'(a10,1f12.6)')'alpha     ',param(7)
C      write(52,'(a10,1f12.6)')'beta      ',param(8)
C      write(52,'(a10,1f12.6)')'rb        ',param(9)
C      write(52,'(a10,1f12.6)')'zshift    ',param(10)
C      write(52,'(a10,1i12)')'ntable    ',int(param(11))
C      write(52,'(a10,1i12)')'smooth    ',int(idre_param(6))
C      write(52,'(a10,1f12.6)')'normal    ',idre_param(7)
C      close(52)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c Let's integrate local emission over the accretion disc
      call idre(ear,ne,photar,idre_param,cmodel,ne_loc,ear_loc,far_loc)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cc final spectrum output -- write ear() and photar() into file:
C      open(61,file='kyr_photar.dat',status='unknown')
C      do i=1,ne
C       write(61,'(2g14.6)')0.5*(ear(i-1)+ear(i)),
C     $                     photar(i)/(ear(i)-ear(i-1))
C      enddo
C      close(61)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      return
      end
c=============================================================================
