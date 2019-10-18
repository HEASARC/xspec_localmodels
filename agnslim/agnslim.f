      Subroutine agnslim(ear,ne,param,ifl,photar,photer)

      implicit none

      integer NPAR
      parameter(NPAR=14)

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
c     param(2) distance
c     param(3) log mass accretion rate in L/LEdd
c     param(4) spin parameter
c     param(5) cosi 0.5 is isotripic
c     param(6) Te of hot compton
c     param(7) Te of warm compton
c     param(8) gamma for hot compton
c     param(9) gamma for warm compton   if -ve plot only warm comp
c     param(10) R_hot
c     param(11) R_warm
c     param(12) log10 rout/rg
c     param(13) rin
c     param(14) redshift

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      zfac = 1.0 + param(14)

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

         call amydiskslim(e,NN,param,ifl,ph)

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

      subroutine amydiskslim(ear,ne,param,ifl,photar)


c     program to integrate the disk equations from shakura-sunyaev disk
c     as given by Novikov and Thorne

      implicit none
      double precision rin,rin0
      double precision m,mdot,rsg,mdotstart
      double precision astar,z1,z2,rms,rh
      double precision rcor,amytemp
      double precision seeddis,lumipl!add
      double precision rpow !add add add
      double precision trepd
      double precision tslim,tnt ! add add add
      double precision rgcm,pi,t,r,dr,dlogr,dflux,tprev
      double precision en,kkev,h,kevhz,d0,d,logrout
      double precision dllth,dlhth,dldiskint !add add add
      integer i,ipow,icor,iout,imax,n,ne,ifl
      logical first,firstdisk,firstdisk2,firstdisk3,firstpl
      double precision flux(ne),displ
      double precision mdotedd,alpha,eff,grad
      double precision cosi
      double precision gammah,gammas,teh,tes
      double precision ledd
      double precision tedd,tsg !critical temperature[K] for soft compton
      double precision trepdedd !critical temperature[K] for soft compton
      double precision redd, fcrit,t4r2,t4r2prev!radius at which T>Tedd and radius at Tmax
      double precision tedd0
      double precision factor
      character(255) comment
      real ear(0:ne),photar(ne),param(*)
      real lpar(5),lphot(ne),lphote(ne),hpar(5),hphot(ne),hphote(ne)
      real lphotall(ne),hphotall(ne)  !add add add

      imax=2000 !integrate Rin(or Rms) to Rsg
      iout=1000 !integrage Rout to Rsg
      icor=10
      ipow=100

c     constants
      h=6.62617d-27 ! [erg s]
      kkev=1.16048d7  ![K/keV]
      kevhz=2.417965d17 ![Hz/keV]
      pi=4.0*atan(1.0)      

c     system parameters
      m=dble(param(1))         !in solar units
      d0=dble(param(2))     !in Mpc
      cosi=dble(param(5))
      mdotedd=dble(10**(param(3)))      !in L/Ledd if -ve plot disc
      astar=dble(param(4))
      mdotstart=6.0d0
      alpha=0.1
      factor=2.39d0
c     corona parameters
      gammah=dble(abs(param(8)))
      gammas=dble(abs(param(9)))
      tes=dble(abs(param(7)))/511.0
      teh=dble(abs(param(6)))/511.0
      rcor=dble(abs(param(11)))

      d=d0*1d6*3.085677d18 !mod
      rgcm=1.477d5*m  ![cm]
      tedd=2.27651d5*(1.0d7/m)**0.25
      trepdedd=0.0d0
      tedd0=0.0d0
c     get rms
      z1=((1-astar**2)**(1./3.))
      z1=z1*(((1+astar)**(1./3.))+((1-astar)**(1./3.)))
      z1=1+z1
      z2=sqrt(3*astar*astar+z1*z1)
      rms=3.+z2-sqrt((3.-z1)*(3.+z1+2.0*z2))
      eff=1.0-sqrt(1.0-2.0/(3.0*rms))  !0.057 for a=0
      rh=1+sqrt(1-astar*astar)
c     mdot in g/s Ledd=1.39e38 m and c2=c*c=8.99e20[cm2/s2]
      mdot=mdotedd*1.39d18*m/(8.98755*eff) ![g/s]
      ledd=1.39d38*m
      fcrit=3.749d23*(1.0d7/m) * factor

c   definition of Rin by Watarai et al.
      if(mdotedd.le.mdotstart)then  !L<2Le*factor
         rin0=rms
      elseif(mdotedd.le.100.0d0)then !L<6Le*factor
c         rin0=rms*sqrt(10.0d0/mdotedd) !=2Rg @9Ledd
         rin0=rms*(mdotedd/mdotstart)
     &        **(log10(rh/rms)/log10(100.0d0/mdotstart))
         rin0=dmax1(rin0,rh) ! lowerlimit
      else
         rin0=rh
      endif

      if (param(13).eq.-1.0d0)then  !rin=-1
         rin=rin0               !fix from Watarai el
      elseif (mdotedd.ge.mdotstart)then !mdot>6
         if (param(13).le.rh)then
            rin=rh
         else  
            rin=dble(param(13))
         endif
      else !mdot <6
         if (param(13).le.rms)then
            rin=rms
         else
            rin=dble(param(13))
         endif
      endif
         


      logrout=dble(param(12))   !log10 Rout/Rg

c     iv -ve then calculate rout=rsg frmo laor & netzer 1989
      if (param(12).lt.0.0) then
         logrout=((m/1e9)**(-2.0/9.0)) * ((mdotedd)**(4.0/9.0))
         logrout=2150*logrout* (alpha**(2.0/9.0))
         logrout=log10(logrout)
      end if
      rsg=10**logrout !rout
      tsg=amytemp(m,astar,mdot,rms,rsg)

c------------aaaaa
      rpow=dmin1(dble(param(10)),rsg)
      rpow=dmax1(rpow,rin)
      if (rpow.le.rin)then
         ipow=0.0d0
      endif
      if(param(10).gt.rsg)then
         print *, "==========================================="
         print *, "Rhot is larger than Rout <reset>"
         print *, "Rhot=",rpow
         print *, "==========================================="
      elseif(param(10).lt.rin)then
        print *, "==========================================="
        print *, "Rhot is smaller than Rin <reset>"
        print *, "Rhot=",rpow
        print *, "==========================================="
      endif
c---------

c     initialise 
      do n=1,ne,1
         photar(n)=0.0
         flux(n)=0.0
c         photarint(n)=0.0
c         photarseed(n)=0.0
c         fluxint(n)=0.0
c        fluxseed(n)=0.0
         lphotall(n)=0.0 !add add add 
         hphotall(n)=0.0 !add add add 
      end do
      
      first=.true.
      firstdisk=.true.
      firstdisk2=.true.
      firstdisk3=.true.
      firstpl=.true.
      displ=0.0d0
      seeddis=0.0d0
      tprev=0.0d0
      t4r2prev=0.0d0
      dflux=0.0d0
c      dfluxint=0.0d0
c      dfluxseed=0.0d0
      redd=0.0d0

c calculate Redd@L>Ledd
      do i=1,imax,1  !
         dlogr=log10(rsg/rms)/float(imax)
         r=10**(log10(rms)+float(i-1)*dlogr+dlogr/2.0)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         t=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
c---------

c         if (param(3).gt.log10(factor))then !change
         t4r2 = t**4*r*r
         if ((t4r2.le.fcrit).and.(t4r2prev.gt.fcrit)) then
            if(firstdisk3)then
                      redd=r
               tedd=t
               firstdisk3=.false.
            endif
c------?????????
         elseif ((t4r2.gt.fcrit).and.(i.eq.imax)) then
            redd=rsg*t4r2/fcrit
c------?????????
         endif
         t4r2prev=t4r2
c        endif
      enddo

c calculation for intersept seed photons without considerint reprocess
c      rcor=rin

      if(param(3).gt.0.0d0)then
         tedd0=amytemp(m,astar,mdot,rms,redd)
      endif

!2nd run for T(R)=slim disk and Tmax for all L
      do i=1,imax,1
c            dlogr=log10(rsg/rms)/float(imax)
            dlogr=log10(rsg/rin)/float(imax)
            r=10**(log10(rin)+float(i-1)*dlogr+dlogr/2.0)
            dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
            if (r.gt.rms)then
               tnt=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
            else
               tnt=0.0d0
            endif
            t=tnt
            if ((r.le.rsg).and.(r.lt.redd)) then
               tslim=tedd0*(r/redd)**(-0.5)
               if(mdotedd.ge.mdotstart)then
                  t=tslim
               else !mdotedd<10Ledd
c                  t=dmin1(tnt,tslim)
                  tnt=dmin1(tnt,tslim)
                 grad=4.0*log10(tslim/tnt)/log10(mdotstart/factor)
                 t=4.0*log10(tnt)+log10(mdotedd/factor)*grad
                 t=10**t
                 t=t**0.25
               endif
            endif
c calculate Ldiss,hot
        if(rpow.gt.rin)then
            if(firstpl)then
                displ=displ+2*2*pi*r*dr*rgcm*rgcm*5.670367e-5*t**4
                if (r.ge.rpow)then
                    firstpl=.false.
                endif
            endif
        else
            displ=0.0d0
        endif

c     setting Rcor for soft compton
            tprev=t
c            rcor=dmax1(rcor,rpow)
      enddo                     !end of r

        lumipl=displ
c------- set rcor------
      if (rcor.le.rpow)then
         icor=0
         rcor=rpow
      print *, "==========================================="
      print *, "Rwarm is smaller than Rhot <reset as Rwarm=Rhot>"
      print *, "Rwarm=", rcor
      print *, "==========================================="
      endif
c--------------------
c------- set rcor------0109_rev
      if (rcor.gt.rsg)then
         rcor=rsg
      print *, "==========================================="
      print *, "Rwarm is larger than Rout <reset as Rwarm=Rout>"
      print *, "Rwarm=", rcor
      print *, "==========================================="
      endif
c--------------------
c     print *, fsc,mdot,mdotscd,trepd,trepsc
      do i=1,ipow+icor+iout,1
         if (i.le.ipow) then
            dlogr=log10(rpow/rin)/float(ipow)
            r=10**(log10(rin)+float(i-1)*dlogr+dlogr/2.0)
         elseif (i.le.ipow+icor) then
            dlogr=log10(rcor/rpow)/float(icor)
            r=10**(log10(rpow)+float(i-ipow-1)*dlogr+dlogr/2.0)
         else
            dlogr=log10(rsg/rcor)/float(iout)
            r=10**(log10(rcor)+float(i-icor-ipow-1)*dlogr+dlogr/2.0)
         end if


c     geometry factor

         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         if(r.gt.rms)then
            tnt=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
         else
            tnt=0.0d0
         endif
         t=tnt
            if ((r.le.rsg).and.(r.lt.redd))then
               tslim=tedd0*(r/redd)**(-0.5)
c               if((tslim.lt.t).or.(mdotedd.gt.10.0d0))then  ! 1.6d0 is checked by eye on Mdot-L
               if(mdotedd.ge.mdotstart)then
                  t=tslim
               else !mdotedd<10Ledd
                  tnt=dmin1(tnt,tslim)
                 t=4.0*log10(tslim/tnt)/log10(mdotstart/factor)
                 t=4.0*log10(tnt)+log10(mdotedd/factor)*t
                 t=10**t
                 t=t**0.25
               endif
            endif

        trepd=t



c------------------------
         t=t/kkev ! [keV] 
         trepd=trepd/kkev  ![keV]

c        go over each photon energy - midpoint of bin
         do n=1,ne,1
            en=dble(log10(ear(n))+log10(ear(n-1)))
            en=en/2.0
            en=10**en

c           do blackbody spectrum  @r>rcor
            if ((en.lt.30.0*trepd).and.(r.gt.rcor).and.(r.gt.rin)) then
               dflux=pi*2.0*h*((en*kevhz)**3)/8.98755d20
               dflux=dflux*4.0*pi*r*dr*rgcm*rgcm
               dflux=dflux/(exp(en/(trepd))-1.0)
            else
               dflux=0.0d0
            end if

            flux(n)=flux(n)+(dflux)
         end do 

c powerlaw


c  start ---- slice nthcomp
         dllth=0.0
         dlhth=0.0
         dldiskint=2.0*2.0*pi*r*dr*rgcm*rgcm*5.670367e-5
     &        *(( trepd*kkev)**4) !mod mod mod
         dldiskint=dldiskint/(4.0*pi*d*d) ! erg/s/cm^2

c   ----- slithe ncomp for hard comp geometry =1

        if (i.le.ipow) then
            ifl=1
            hpar(1)=sngl(gammah)
            hpar(2)=abs(param(6)) !kTe if - plot comp
            hpar(3)=sngl(trepd) !set to reprocessed effective temperature mod mod
            hpar(4)=0.0         ! int type
            hpar(5)=0.0

c thcomp is called but scale may be wrong yet
            call donthcomp(ear,ne,hpar,ifl,hphot,hphote)
            do n=1,ne,1
            dlhth=dlhth+hphot(n)*ear(n)*kevhz*h !photons/s/cm2/Hz *Hz *keV
            end do
        end if

         if ((i.le.ipow+icor).and.(i.gt.ipow)) then
            lpar(1)=sngl(gammas)
            lpar(2)=abs(param(7)) !kTe if - plot comp
            lpar(3)=sngl(trepd) !set to reprocessed effective temperature mod mod
            lpar(4)=0.0         ! int type
            lpar(5)=0.0
c thcomp is called but scale may be wrong yet               
               ifl=1
               call donthcomp(ear,ne,lpar,ifl,lphot,lphote)
               do n=1,ne,1
                  dllth=dllth+lphot(n)*ear(n)*kevhz*h !photons/s/cm2/Hz *Hz *keV
               end do
         end if

         do n=1,ne,1
            if (dllth.eq.0) then
               lphot(n)=0.0
            else 
               lphot(n)=lphot(n)*sngl(dldiskint/dllth)
            endif 
            lphotall(n)=lphotall(n)+lphot(n)
            if (dlhth.eq.0) then
                hphot(n)=0.0
            else
                hphot(n)=hphot(n)*sngl(dldiskint/dlhth)
            endif
            hphotall(n)=hphotall(n)+hphot(n)
         end do 

c  end ----- mod mod mod mod slice nthcomp
c----------- end of slice pow

      end do                    !end of r

c==== powerlaw


c      tot=0.0
      do n=1,ne,1
c        this is ergs cm^2 s-1 Hz^-1 
         flux(n)=flux(n)/(4.0*pi*d*d)
c        photons is flux/hv - photons cm^2 s-1 Hz^-1
         flux(n)=flux(n)/(h*kevhz*ear(n))
c        now multiply by energy band in Hz
         photar(n) = sngl(flux(n)*kevhz) * (ear(n)-ear(n-1)) !kept 
      end do

      displ=displ/(4.0*pi*d*d)!add add add
      lumipl=lumipl/(4.0*pi*d*d)!add add add
c     now add high temp compton emission

c      write(comment,*) '------hot compton-----'
c      CALL xwrite(comment, 20)
c      WRITE(comment,*) 'Lhot/Ledd=',lumipl*4.0*pi*d*d/ledd
c      CALL xwrite(comment, 20)
      write(comment,'(3(a,1pg13.4))') 'rin=',rin, '( rin_calc=',rin0,')'
      CALL xwrite(comment, 20)
      WRITE(comment,'(3(a,1pg13.4))') 'r_hot=', rpow,
     &  'r_warm=',rcor, 'r_edd=', redd
      CALL xwrite(comment, 20)
      WRITE(comment,'(3(a,1pg13.4))') 'Tcri(slim)=',
     & tedd0*(rcor/redd)**(-0.5)
      CALL xwrite(comment, 20)
      write(comment, '(3(a,1pg13.4))') 'log_rout= ',logrout
      CALL xwrite(comment, 20)
      write(comment,'(3(a,1pg13.4))') 'tout=', tsg
      CALL xwrite(comment, 20)
      write(comment,*) '------system-----'
      CALL xwrite(comment, 20)
      write(comment,'(3(a,1pg9.3))') 'fcrit=',
     & fcrit*4*3.14*5.67e-5*rgcm**2/ledd,
     & '(critical flux is the Eddington flux)'
      CALL xwrite(comment, 20)
      write(comment,'(3(a,1pg13.4))') 'efficiency=', eff
      CALL xwrite(comment, 20)
      write(comment,'(3(a,1pg13.4))') 'r_H=', rh
      CALL xwrite(comment, 20)


      do n=1,ne,1
         if ((param(6).lt.0.0).or.(param(7).lt.0.0).or.
     &          (param(9).lt.0.0)) then
             if (param(6).lt.0.0)photar(n)=
     &                  SNGL(hphotall(n)*cosi/0.5)
             if (param(7).lt.0.0)
     &             photar(n)=SNGL(lphotall(n)*cosi/0.5)
             if (param(9).lt.0.0) photar(n)=photar(n)
     &             *sngl(cosi/0.5)
          else
             photar(n)=SNGL(photar(n)*cosi/0.5
     &             +lphotall(n)*sngl(cosi/0.5)
     &             +hphotall(n)*sngl(cosi/0.5))
          endif
      end do
           
      return
      end

