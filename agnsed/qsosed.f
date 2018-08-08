
      Subroutine qsosed(ear,ne,param,ifl,photar,photer)

      implicit none

      integer NPAR
      parameter(NPAR=6)

      integer ne, ifl
      real ear(0:ne), photar(ne), param(NPAR), photer(ne)

!------------------------------------------------------------------------------------ 
! This subroutine checks whether any of the parameters have changed, if they have, 
! then it creates a new energy binning scheme with 5000 bins increasing exponentially
! in size.  The input parameters and new binning scheme are used in a call to the 
! subroutine 'asdiskf'.  The returned fluxes are then corrected for redshift and
! shifted back into the original energy binning scheme.
!------------------------------------------------------------------------------------ 

!     ear = Energy array
!     ne = size of energy array (1 more than flux array as gives boundaries)
!     param = array containing model parameters
!     ifl = spectrum number
!     photar = output flux array (ne-1 parameters)

      integer NN
      parameter(NN=5000)

      integer n,i
      real e(0:NN),ph(NN)
      real dloge,zfac,oldpar(NPAR),newemin,newemax
      logical parchange, echange

      REAL fstart(ne), fend(ne)
      integer istart(ne), iend(ne)

      save oldpar,ph,e

c     param(1) mass in solar
c     param(2) luminosity distance
c     param(3) log mass accretion rate in L/LEdd
c     param(4) astar
c     param(5) cosi 0.5 is isotripic  (14)
c     param(6) redshift (11)

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      zfac = 1.0 + param(6)

!     test whether any parameters have changed
      parchange=.false.
      do i=1,NPAR,1
         if (param(i).ne.oldpar(i)) parchange=.true.
      end do

!     test whether own energy grid needs to be changed
      if (ear(0).eq.0.0) then 
         newemin=ear(1)-(ear(2)-ear(1))/10.0
      else
         newemin=min(1.0e-5,ear(0)*zfac)
      end if
      newemax=max(1.0e3,ear(ne)*zfac)

      if ((e(0).ne.newemin).or.(e(NN).ne.newemax)) then
         echange=.true.
      end if

!     if necessary recalculate the energy grid

      if (echange) then

         dloge=log10(newemax/newemin)/float(NN)
         e(0) = newemin
         do n=1,NN,1
!           populate the new energy array and zero flux
            e(n)=10**(log10(e(0))+dloge*float(n))
         end do

      endif
         
!     call the irradiated disc if parameters or energy grid have changed
      if (parchange .or. echange) then 

         call qsoamydiskf(e,NN,param,ifl,ph)

!        now redshift the energy bins back
         do n=0,NN,1
            e(n)=e(n)/zfac
         end do

!        and now correct the flux for redshift
         do n=1,NN,1
            ph(n)=ph(n)/zfac
         end do

!        save the new parameters into the oldpar array
         do i=1,NPAR,1
            oldpar(i)=param(i)
         end do

      end if

!     rebin the calculated fluxes back onto original energy grid
      CALL inibin(NN, e, ne, ear, istart, iend, fstart, fend, 0)
      CALL erebin(NN, ph, ne, istart, iend, fstart, fend, photar)
     
      return
      end
      subroutine qsoamydiskf(ear,ne,param,ifl,photar)


c     program to integrate the disk equations from shakura-sunyaev disk
c     as given by Novikov and Thorne

      implicit none
      double precision Htmax,tseedpl
      double precision m,mdot,rsg
      double precision rep,displ0 ! add add add :g is covering fraction
      double precision astar,z1,z2,rms
      double precision pow,rcor,amytemp,amytempreph
      double precision seeddis,lumipl!add
      double precision rpow !add add add
      double precision trepd,trepdis ! add add add
      double precision tnt ! add add add
      double precision rgcm,pi,t,t0,r,dr,dlogr,dflux,tprev
      double precision dfluxint,dfluxseed !add add add
      double precision en,kkev,h,kevhz,d0,d,logrout
      double precision dllth,dldiskint !add add add
      integer i,icor,iout,imax,n,ne,ifl
      logical first,firstdisk,firstdisk2,firstdisk3,firstpl
      double precision flux(ne),displ
      double precision mdotedd,alpha,eff
      double precision fluxint(ne),fluxseed(ne) ! add add add
      double precision Ht, a !scale height and albedo add add add
      double precision cosi
      double precision gammah,gammas,tauh,taus,teh,tes,ys,yh !add
      double precision kteh,ktes
      double precision ysb,tausb,yhb,tauhb
      double precision theta0 !geometry
      double precision ledd
      double precision rmax,tsg !critical temperature[K] for soft compton
      real ear(0:ne),photar(ne),param(*)
      real lpar(5),lphot(ne),lphote(ne),hpar(5),hphot(ne),hphote(ne)
      real photarint(ne),photarseed(ne)  !add add add
      real lphotall(ne),hphotall(ne)  !add add add


c     constants
      h=6.62617d-27 ! [erg s]
      kkev=1.16048d7  ![K/keV]
      kevhz=2.417965d17 ![Hz/keV]
      pi=4.0*atan(1.0)      

c     system parameters
      tseedpl=-1.0d0
      m=dble(param(1))         !in solar units
      Htmax=100.0d0  !upper limit of scaleheight
      d0=dble(param(2))     !in Mpc
      cosi=dble(param(5))
      mdotedd=dble(10**(param(3)))      !in L/Ledd if -ve plot disc
      astar=dble(param(4))     
      alpha=0.1
      rep=1.0d0 !@L>Le reprocess is not considered
c     corona parameters
      gammas=2.5d0 !passive disc
      ktes=0.2d0
      tes=ktes/511.0
      kteh=100.0d0
      teh=kteh/511.0
      taus=(3.0/(tes*((gammas+0.5)**2.0-2.25))+2.25)**0.5-1.5
      ys=4.0*(tes+4.0*tes**2)*(taus+taus**2)
      ysb=(gammas*4.0/9.0)**(-4.5)
      tausb=(ysb/(4.0*tes))**0.5
      ledd=1.39d38*m
      displ0=0.020d0*ledd
      a=0.30d0!albedo
      d=d0*1d6*3.085677d18 !mod
      rgcm=1.477d5*m  ![cm]
c     get rms
      z1=((1-astar**2)**(1./3.))
      z1=z1*(((1+astar)**(1./3.))+((1-astar)**(1./3.)))
      z1=1+z1
      z2=sqrt(3*astar*astar+z1*z1)
      if(astar.ge.0.0d0)then
        rms=3.+z2-sqrt((3.-z1)*(3.+z1+2.0*z2))
      else
        rms=3.+z2+sqrt((3.-z1)*(3.+z1+2.0*z2))  !mod
      endif
      eff=1.0-sqrt(1.0-2.0/(3.0*rms))  !0.057 for a=0

c     mdot in g/s Ledd=1.39e38 m and c2=c*c=8.99e20[cm2/s2] 
      mdot=mdotedd*1.39d18*m/(8.98755*eff) ![g/s]

c     iv -ve then calculate rout=rsg frmo laor & netzer 1989
         logrout=((m/1e9)**(-2.0/9.0)) * ((mdotedd)**(4.0/9.0))
         logrout=2150*logrout* (alpha**(2.0/9.0))
         logrout=log10(logrout)
         rsg=10**logrout !rout


c---------

c     initialise 
      do n=1,ne,1
         photar(n)=0.0
         flux(n)=0.0
         photarint(n)=0.0
         photarseed(n)=0.0
         fluxint(n)=0.0
         fluxseed(n)=0.0
         lphotall(n)=0.0 !add add add 
         hphotall(n)=0.0 !add add add 
      end do
      
      imax=2000
      iout=1000
      icor=10
      first=.true.
      firstdisk=.true.
      firstdisk2=.true.
      firstdisk3=.true.
      firstpl=.true.
      displ=0.0d0
      seeddis=0.0d0
      t0=0.0d0
      rmax=0.0d0
      tprev=0.0d0
      dflux=0.0d0
      dfluxint=0.0d0
      dfluxseed=0.0d0

c calculate  Rpow@L<Ledd first run

      do i=1,imax,1  !outer sc + disk is calculated
         dlogr=log10(rsg/rms)/float(imax)
         r=10**(log10(rms)+float(i-1)*dlogr+dlogr/2.0)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         t=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
        if(firstpl)then
            displ=displ+2*2*pi*r*dr*rgcm*rgcm*5.670367e-5*t**4
            rpow=r
            if (displ.ge.displ0)then
                firstpl=.false.
            endif
        endif
      enddo
c---------

c calculation for intersept seed photons without considerint reprocess
      rpow=dmin1(rpow,rsg)
      rpow=dmax1(rpow,rms)
      rcor=2*rpow
      rcor=dmin1(rcor,rsg)
      print *,"-------setup for qsosed------"
      write(*,*) 'Gamma_warm=',gammas,'kTe_warm=',ktes
      print *,'Gamma_hard=internal calculation  ',"kTe_hot=",kteh
      print *,"albedo=",a
      print *,"Rwarm/Rhot=",rcor/rpow
      print *,"Htmax=",Htmax

      Ht=dmin1(rpow,Htmax)

!2nd run for T(R)=slim disk and Tmax for all L
      tsg=amytemp(m,astar,mdot,rms,rsg)
      do i=1,imax,1
c            dlogr=log10(rsg/rms)/float(imax)
            dlogr=log10(rsg/rms)/float(imax)
            r=10**(log10(rms)+float(i-1)*dlogr+dlogr/2.0)
            dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
            tnt=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
            t=tnt

c     setting Rcor for soft compton
            tprev=t
c            rcor=dmax1(rcor,rpow)
            if (r.ge.rpow)then
               if (rep.eq.1.0d0)then
                  trepdis=amytempreph(m,astar,mdot,rms,r,
     &                 rep,displ,Ht,a) ! outer disk temperature including rep of Ldis only
               elseif(rep.ge.0.0d0)then
                  trepdis=t
               endif
               theta0=asin(Ht/r)                           
               seeddis=seeddis+2.0*2.0*pi*r*dr*rgcm*rgcm*5.670367e-5
     &              *(trepdis**4)*(theta0-0.5*sin(2*theta0))/pi !intersept seed photon at radius r
c               seedi=seedi+2.0*2.0*pi*r*dr*rgcm*rgcm
c     &              *5.670367e-5*(t**4)*(theta0-0.5*sin(2*theta0))/pi !intersept seed photon at radius r
            endif
      enddo                     !end of r

      lumipl=displ+seeddis ![erg/s]
c------- set rcor------
      if (rcor.le.rpow)then
        icor=0
        rcor=rpow
        print *, "==========================================="
        print *, "Rwarm is smaller than Rhot <reset as Rwarm=Rhot>"
        print *, "==========================================="
      endif
c--------------------




c     print *, fsc,mdot,mdotscd,trepd,trepsc

      do i=1,icor+iout,1
         if (i.le.icor) then
            dlogr=log10(rcor/rpow)/float(icor)
            r=10**(log10(rpow)+float(i-1)*dlogr+dlogr/2.0)
         else
            dlogr=log10(rsg/rcor)/float(iout)
            r=10**(log10(rcor)+float(i-icor-1)*dlogr+dlogr/2.0)
         end if
c     geometry factor
         theta0=asin(Ht/r)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         tnt=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
         t=tnt

         if (rep.gt.0) then     !if reprocess is included
            trepd=amytempreph(m,astar,mdot,rms,r,
     &           rep,lumipl,Ht,a) ! outer disk temperature including rep
         else                   ! no reprocess    
            trepd=amytemp(m,astar,mdot,rms,r)
         endif
         
      
c setting seed phton temperature for PL
        if (first) then
            if (rcor.gt.rpow) then
                t0=exp(ysb)*trepd/kkev !with soft comptonmmd
            else
                t0=trepd/kkev !if there is no low compton region
            endif
            first=.false.
        end if

c------------------------
         t=t/kkev ! [keV] 
         trepd=trepd/kkev  ![keV]

c        go over each photon energy - midpoint of bin
         do n=1,ne,1
            en=dble(log10(ear(n))+log10(ear(n-1)))
            en=en/2.0
            en=10**en

c           do blackbody spectrum  @r>rcor
c            if ((en.lt.30.0*trepd).and.(r.gt.rcor)) then
            if ((en.lt.30.0*trepd).and.(r.gt.rcor)) then
               dflux=pi*2.0*h*((en*kevhz)**3)/8.98755d20
               dflux=dflux*4.0*pi*r*dr*rgcm*rgcm
               dflux=dflux/(exp(en/(trepd))-1.0)
            else
               dflux=0.0d0
            end if

            flux(n)=flux(n)+(dflux)

            

c            do black body @r>rpow
            if ((en.lt.30.0*t).and.(r.gt.rpow)) then
               dfluxint=pi*2.0*h*((en*kevhz)**3)/8.98755d20
               dfluxint=dfluxint*4.0*pi*r*dr*rgcm*rgcm
               dfluxint=dfluxint/(exp(en/(t))-1.0)
               else
               dfluxint=0.0d0
            end if
            fluxint(n)=fluxint(n)+(dfluxint)

c     do seed photon for hard compton
            if ((en.lt.30.0*trepd).and.(r.gt.rpow)) then
               dfluxseed=pi*2.0*h*((en*kevhz)**3)/8.98755d20
               dfluxseed=dfluxseed*4.0*pi*r*dr*rgcm*rgcm
               dfluxseed=dfluxseed/(exp(en/(trepd))-1.0)
               dfluxseed=dfluxseed*(theta0-0.5*sin(2.0*theta0))/pi
               else
               dfluxseed=0.0d0
            end if
            fluxseed(n)=fluxseed(n)+(dfluxseed)

         end do 


c  start ---- mod mod mod mod slice nthcomp
         dllth=0.0
         dldiskint=2.0*2.0*pi*r*dr*rgcm*rgcm*5.670367e-5
     &        *(( trepd*kkev)**4) !mod mod mod
         dldiskint=dldiskint/(4.0*pi*d*d) ! erg/s/cm^2

         if (i.le.icor) then
            lpar(1)=sngl(gammas)
            lpar(2)=sngl(ktes) !kTe if - plot comp
            lpar(3)=sngl(trepd) !set to reprocessed effective temperature mod mod
            lpar(4)=0.0         ! int type
            lpar(5)=0.0

               ifl=1
c thcomp is called but scale may be wrong yet               
               call donthcomp(ear,ne,lpar,ifl,lphot,lphote)
               do n=1,ne,1
                  dllth=dllth+lphot(n)*ear(n)*kevhz*h !photons/s/cm2/Hz *Hz *keV
               end do
         end if
c         print *,r,dldiskint,dllth,t*fcol,trep*fcol,tss/kkev
         do n=1,ne,1
            if (dllth.eq.0) then
               lphot(n)=0.0
            else 
               lphot(n)=lphot(n)*sngl(dldiskint/dllth)
            endif 
            lphotall(n)=lphotall(n)+lphot(n)               
         end do 

c  end ----- mod mod mod mod slice nthcomp

      end do                    !end of r

c==== powerlaw

        gammah=7.0/3.0*((displ+seeddis)/seeddis-1.0)**(-0.10)

c      if (gammah.le.1.4d0)then
c         gammah=1.4
c      endif

      tauh=(3.0/(teh*((gammah+0.5)**2.0-2.25))+2.25)**0.5-1.5
      yh=4.0*(teh+4.0*teh**2)*(tauh+tauh**2)
      yhb=(gammah*4.0/9.0)**(-4.5)
      tauhb=(1+4.0*yhb/(4.0*(teh+4.0*teh*teh)))**0.5-1
      tauhb=tauhb/2



c      tot=0.0
      do n=1,ne,1
c        this is ergs cm^2 s-1 Hz^-1 
         flux(n)=flux(n)/(4.0*pi*d*d)
         fluxint(n)=fluxint(n)/(4.0*pi*d*d) 
         fluxseed(n)=fluxseed(n)/(4.0*pi*d*d) 

c        photons is flux/hv - photons cm^2 s-1 Hz^-1
         flux(n)=flux(n)/(h*kevhz*ear(n))
         fluxint(n)=fluxint(n)/(h*kevhz*ear(n))!intrinsic disk flux at r>rpow
         fluxseed(n)=fluxseed(n)/(h*kevhz*ear(n))!intrinsic disk flux at r>rpow

c        now multiply by energy band in Hz
         photar(n) = sngl(flux(n)*kevhz) * (ear(n)-ear(n-1)) !kept 
         photarint(n) = sngl(fluxint(n)*kevhz) * (ear(n)-ear(n-1)) !kept 
         photarseed(n) = sngl(fluxseed(n)*kevhz) * (ear(n)-ear(n-1)) !kept 

      end do

      displ=displ/(4.0*pi*d*d)!add add add
      lumipl=lumipl/(4.0*pi*d*d)!add add add
      displ0=displ0/(4.0*pi*d*d)!add add add
c     now add high temp compton emission

      hpar(1)=sngl(gammah)
      hpar(2)=sngl(kteh)     !Te
      hpar(3)=SNGL(t0)
      hpar(4)=0.0
      hpar(5)=0.0



         call donthcomp(ear,ne,hpar,ifl,hphot,hphote)
         pow=0.0
         do n=1,ne,1
            pow=pow+hphot(n)*ear(n)*kevhz*h !hard compton luminosity
         end do 

      print *,"                         "
      print *,"------hot compton-----"
      print *,"Gamma_hot=",gammah,"Rhot=",rpow
      print *,"T(Rhot)nt=",amytempreph(m,astar,mdot,rms,rpow,rep,
     &        lumipl,Ht,a),"Tseed=",t0*kkev
      print *,"Ldis,hot/Ledd=",displ*4.0*pi*d*d/ledd,"Lhot/Ledd=",
     &        lumipl*4.0*pi*d*d/ledd
      print *,"                         "
      print *,"------warm compton-----"
      print *,"T(Rw)=",amytempreph(m,astar,mdot,rms,rcor,rep,
     &        lumipl,Ht,a),"Rwarm=",rcor
      print *,"                         "
      print *,"------Rout-----"
      print *,"rout(rsg)=", rsg
      print *,"tout=", tsg
      print *,"-----------"

      do n=1,ne,1
            photar(n)=SNGL(photar(n)*cosi/0.5
     &               +lphotall(n)*sngl(cosi/0.5)
     &               +hphot(n)*sngl((lumipl/pow)))

      end do
           

      return
      end

