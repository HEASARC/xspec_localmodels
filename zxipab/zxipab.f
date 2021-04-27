      SUBROUTINE zxipab(ear, ne, param, ifl, photar, photer)

      implicit none

      INTEGER ne, ifl
      REAL ear(0:ne),param(5),photar(ne),photer(ne)

C---
C XSPEC model subroutine for redshifted "power law covering fraction
C ionized absorption".  Koji is making this up based on pwab and zxipcf models
C---
C number of model parameters: 5
C       1-3     nHmin, nHmax, beta from pwab
C       4       xi
C       5       REDSHIFT
C---
C Arguments :
C    ear      r         i: energy ranges
c    ne       i         i: number of energy ranges
c    param    r         i: model parameter values
c    ifl      i         i: the dataset
c    photar   r         r: the transmission fraction
c    photer   r         r: (output error array)
C---

      INTEGER LNE
      PARAMETER(LNE=100)
      
      REAL eparam(4), e0, e(0:LNE), ph(LNE), phe(LNE), zfac
      REAL nhmin, nhmax, beta, a, abslog, dnh, lognh, nh
      INTEGER ie,i12, inh, nnh

      nhmin=param(1)
      nhmax=param(2)
      beta=param(3)

C shift energies to the emitter frame
      zfac = 1.0 + param(5)
      DO ie = 0, ne
         ear(ie) = ear(ie) * zfac
      ENDDO

c     nh - use 1e22 in preparation for later calculations
      eparam(1) = 1.0
c     logxi
      eparam(2) = param(4)
c     cover fraction = 1
      eparam(3)= 1.0
c     redshift = 0
      eparam(4) = param(5)
      

c     get transmission at 12 keV - little energy array from 10-14 keV
      do ie=0,lne
         e(ie)=10.0+float(ie)*4.0/lne
      end do
      call zxipcf(e, lne, eparam, ifl, ph, phe)
      i12=lne/2
      e0=12.0/((-log(ph(i12)))**(-0.33333))

c     calculate ionized absorber
      call zxipcf(ear, ne, eparam, ifl, photar, photer)

c     do sum so total covering fraction is unity
c     loop round values of Nh. lognh is initially set outside the loop
c     to avoid an undefined variable in the case when nnh=0.

      dnh=0.1
      nnh = 1 + INT((log10(nhmax)-log10(nhmin))/dnh)

      if (beta.ne.-1.0) then
         a=(beta+1.0)/((nhmax**(beta+1))-(nhmin**(beta+1)))
      else
         a=1.0/log(nhmax/nhmin)
      end if
      do ie = 1, ne
         if (ear(ie).gt.12.0) then
            abslog = -(ear(ie)/e0)**(-3.)
         else
            abslog = log( photar( ie ) )
         end if
         photar( ie ) = 0.0
         lognh = log10(nhmin)
         do inh = 1, nnh
            nh=10**(lognh+dnh/2.0)
            photar( ie ) = photar( ie ) + a*(nh**beta)
     &          *exp(abslog*nh)*(10**(lognh+dnh)-10**(lognh))
            lognh = lognh + dnh
         end do
         photar( ie ) = photar( ie ) + a*(nhmax**beta)
     &          *exp(abslog*nhmax)*(nhmax-10**(lognh-dnh))
      end do

C shift energies back to the observed frame

      DO ie = 0, ne
         ear(ie) = ear(ie) / zfac
      ENDDO

      RETURN
      END
