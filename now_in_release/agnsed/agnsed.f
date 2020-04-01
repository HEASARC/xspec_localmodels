
      Subroutine agnsed(ear,ne,param,ifl,photar,photer)

      implicit none

      integer NPAR
      parameter(NPAR=15)

      integer ne, ifl
      real ear(0:ne), photar(ne), param(NPAR), photer(ne)

!------------------------------------------------------------------------------------ 
! This subroutine checks whether any of the parameters have changed, if they have, 
! then it creates a new energy binning scheme with 5000 bins increasing exponentially
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
c     param(6) Te of hard compton (5)
c     param(7) Te of warm compton
c     param(8) gamma hot       (9)
c     param(9) gamma warm      (8)
c     param(10) Rhot
c     param(11) Rwarm (15)
c     param(12) log10 rout/rg (6)
c     param(13) Htmax is upper limit of the corona scale height (16)
c     param(14) reprocess 0=ff 1=on
c     param(15) redshift (11)

c     param(12) reprocess L/L_ave 1=average reprocess, 0:no reprocess <0: intrinsic disk emission from r>rpow 
c this model has no errors
c      param(13) albedo
c      param(17) Tseed if <0 Tseed is set to the disk/sc innermost temperature




      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      zfac = 1.0 + param(15)

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

         call amydiskf(e,NN,param,ifl,ph)

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



      subroutine amydiskf(ear,ne,param,ifl,photar)


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
      integer i,icor,imax,iout,n,ne,ifl
      logical first,firstdisk,firstdisk2,firstdisk3,firstpl
      double precision flux(ne),displ
      double precision mdotedd,alpha,eff
      double precision fluxint(ne),fluxseed(ne) ! add add add
      double precision Ht, a !scale height and albedo add add add
      double precision cosi
      double precision gammah,gammas,tauh,taus,teh,tes,ys,yh !add
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
      Htmax=dble(param(13))  !upper limit of scaleheight
      d0=dble(param(2))     !in Mpc
      cosi=dble(param(5))
      mdotedd=dble(10**(param(3)))      !in L/Ledd if -ve plot disc
      astar=dble(param(4))     
      alpha=0.1
      rep=dble(param(14)) !@L>Le reprocess is not considered
c     corona parameters
      gammas=dble(abs(param(9)))
      tes=dble(abs(param(7)))/511.0
      teh=dble(abs(param(6)))/511.0
      rcor=dble(abs(param(11)))
      taus=(3.0/(tes*((gammas+0.5)**2.0-2.25))+2.25)**0.5-1.5
      ys=4.0*(tes+4.0*tes**2)*(taus+taus**2)
      ysb=(gammas*4.0/9.0)**(-4.5)
      tausb=(ysb/(4.0*tes))**0.5


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
      ledd=1.39d38*m
      print *, 'rms=',rms,'eff=',eff


      logrout=dble(param(12))   !log10 Rout/Rg
c     get reprocess parameter
      a=0.30d0!albedo
            
c     iv -ve then calculate rout=rsg frmo laor & netzer 1989
      if (param(12).lt.0.0) then
         logrout=((m/1e9)**(-2.0/9.0)) * ((mdotedd)**(4.0/9.0))
         logrout=2150*logrout* (alpha**(2.0/9.0))
c         write(*,*) 'rout= ',logrout
         logrout=log10(logrout)
      end if
      rsg=10**logrout !rout
      tsg=amytemp(m,astar,mdot,rms,rsg)

c------------aaaaa
      rpow=dmin1(dble(param(10)),rsg)
      rpow=dmax1(rpow,rms)
c      print *,"Rhot=",rpow
      if(param(10).gt.rsg)then
         print *, "==========================================="
         print *, "Rhot is larger than Rout <reset>"
         print *, "==========================================="
      endif
      if(param(10).lt.rms)then
        print *, "==========================================="
        print *, "Rhot is smaller than Rms <reset>"
        print *, "==========================================="
        param(10)=sngl(rms)
      endif
      Ht=dmin1(rpow,Htmax)

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

c calculate Redd@L>Ledd and Rpow@L<Ledd first run
      do i=1,imax,1  !outer sc + disk is calculated
         dlogr=log10(rsg/rms)/float(imax)
         r=10**(log10(rms)+float(i-1)*dlogr+dlogr/2.0)
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)
         t=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
c---------
        if(rpow.gt.rms)then
            if(firstpl)then
                displ=displ+2*2*pi*r*dr*rgcm*rgcm*5.670367e-5*t**4
                if (r.ge.rpow)then
                    firstpl=.false.
                endif
            endif
        endif
        tnt=amytemp(m,astar,mdot,rms,r) !intrinsic effective temperature [K]
        t=tnt

c     setting Rcor for soft compton
            tprev=t
c            rcor=dmax1(rcor,rpow)
            if (r.ge.rpow)then
               if (rep.eq.1.0d0)then
                  trepdis=amytempreph(m,astar,mdot,rms,r,
     &                 rep,displ,Ht,a) ! outer disk temperature including rep of Ldis only
               elseif (rep.eq.0.0d0)then
                  trepdis=amytemp(m,astar,mdot,rms,r)
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
c------- set rcor------0109_rev
      if (rcor.gt.rsg)then
         rcor=rsg
        print *, "==========================================="
        print *, "Rwarm is larger than Rout <reset as Rwarm=Rout>"
        print *, "==========================================="
      endif
c--------------------



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
     &        *((trepd*kkev)**4) !mod mod mod
         dldiskint=dldiskint/(4.0*pi*d*d) ! erg/s/cm^2

         if (i.le.icor) then
            lpar(1)=abs(param(9))
            lpar(2)=abs(param(7)) !kTe if - plot comp
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

      if (param(8).gt.0) then
        gammah=dble(param(8))
      else
        gammah=7.0/3.0*((displ+seeddis)/seeddis-1.0)**(-0.10)
      endif

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
      hpar(2)=abs(param(6))     !Te
      hpar(3)=SNGL(t0)
      hpar(4)=0.0
      hpar(5)=0.0



         call donthcomp(ear,ne,hpar,ifl,hphot,hphote)
         pow=0.0
         do n=1,ne,1
            pow=pow+hphot(n)*ear(n)*kevhz*h !hard compton luminosity
         end do 

      print *,"------hot compton-----"
      print *,"gamma_hot=",gammah
      print *,"T(Rhot)nt=",amytempreph(m,astar,mdot,rms,rpow,rep,
     &        lumipl,Ht,a)
      print *,"Ldis,hot/Ledd=",displ*4.0*pi*d*d/ledd,"Lhot/Ledd=",
     &        lumipl*4.0*pi*d*d/ledd
      print *,"                         "
      print *,"------warm compton-----"
      print *,"T(Rwarm)=",amytempreph(m,astar,mdot,rms,rcor,rep,
     &        lumipl,Ht,a)
      print *,"                         "
      print *,"------Rout-----"
      print *,"rout=", rsg, "tout=", tsg
      print *,"-----------"

      do n=1,ne,1
         if (param(14).eq.1.0)then
            if ((param(6).lt.0.0).or.(param(7).lt.0.0).or.
     &              (param(9).lt.0.0)) then
                if (param(6).lt.0.0)  photar(n)=
     &                 SNGL(hphot(n)*(lumipl/pow)*rep)
                if (param(7).lt.0.0)
     &                 photar(n)=SNGL(lphotall(n)*cosi/0.5)
                if (param(9).lt.0.0) photar(n)=photar(n)
     &                 *sngl(cosi/0.5)
            else
                photar(n)=SNGL(photar(n)*cosi/0.5
     &                +lphotall(n)*sngl(cosi/0.5)
     &                +hphot(n)*sngl((lumipl/pow)*rep))
            end if
         elseif(param(14).eq.0)then
            if ((param(6).lt.0.0).or.(param(7).lt.0.0).or.
     &           (param(9).lt.0.0)) then
               if (param(6).lt.0.0)  photar(n)=
     &              SNGL(hphot(n)*(lumipl/pow))
               if (param(7).lt.0.0)
     &              photar(n)=SNGL(lphotall(n)*cosi/0.5)
               if (param(9).lt.0.0) photar(n)=photar(n)
     &              *sngl(cosi/0.5)
            else
               photar(n)=SNGL(photar(n)*cosi/0.5
     &              +lphotall(n)*sngl(cosi/0.5)
     &              +hphot(n)*sngl((lumipl/pow)))
            end if
         endif
      end do
           

      return
      end

      function amytemp(m0,astar,mdot0,rms,r0)
c     find temperature as a function of mass, spin and r in units of rg

        implicit none
        double precision m0,astar,mdot0,r,r0,rgcm,y,yms
        double precision pi,y1,y2,y3,part1,part2,part3,rms,c,b
        double precision amytemp



        pi=4.0*atan(1.0)
        rgcm=1.477d5*m0

        r=r0
        y=sqrt(r0)
        yms=sqrt(rms)
        y1=2.0*cos((acos(astar)-pi)/3.0)
        y2=2.0*cos((acos(astar)+pi)/3.0)
        y3=-2.0*cos((acos(astar)/3.0))

        part3=3.0*((y3-astar)**2)*log((y-y3)/(yms-y3))
        part3=part3/(y*y3*(y3-y1)*(y3-y2))
        part2=3.0*((y2-astar)**2)*log((y-y2)/(yms-y2))
        part2=part2/(y*y2*(y2-y1)*(y2-y3))
        part1=3.0*((y1-astar)**2)*log((y-y1)/(yms-y1))
        part1=part1/(y*y1*(y1-y2)*(y1-y3))
        c=1.0-yms/y-(3.0*astar/(2.0*y))*log(y/yms)-part1-part2-part3
        b=1.0-3.0/r+2.0*astar/(r**1.5)
        amytemp=3.0*6.67e-8*1.989d33*m0*mdot0
        amytemp=amytemp/(8.0*pi*5.670367e-5*(r*rgcm)**3)
        amytemp=(amytemp*c/b)**0.25
        return
        end

        function amytempreph(m0,astar,mdot0,rms,r0,rep,discu,H,a) !reprocess is included by pl=0.02Le

        implicit none
        double precision m0,astar,mdot0,r,r0,rgcm,y,yms
        double precision pi,y1,y2,y3,part1,part2,part3,rms,c,b
        double precision rep,rep0 !0=no reprocess, 1=reprocess, 1~1.** variable
        double precision amytempreph,factor,discu
        double precision a,H !a=albedo, H is scaleheight

        pi=4.0*atan(1.0)
        rgcm=1.477d5*m0
        r=r0
        y=sqrt(r0)
        yms=sqrt(rms)
        y1=2.0*cos((acos(astar)-pi)/3.0)
        y2=2.0*cos((acos(astar)+pi)/3.0)
        y3=-2.0*cos((acos(astar)/3.0))
c      factor = 4.0*discu/(mdot0*9.0e20)*rep
        rep0=dble(H/rms*(1-a))
        rep0=dble(rep0*(1+H*H/(r*r))**(-1.50))

        factor = 4.0*discu/(mdot0*8.9875e20)*rep0*rep*0.5 !factor 0.5 for one side radiation

        part3=3.0*((y3-astar)**2)*log((y-y3)/(yms-y3))
        part3=part3/(y*y3*(y3-y1)*(y3-y2))
        part2=3.0*((y2-astar)**2)*log((y-y2)/(yms-y2))
        part2=part2/(y*y2*(y2-y1)*(y2-y3))
        part1=3.0*((y1-astar)**2)*log((y-y1)/(yms-y1))
        part1=part1/(y*y1*(y1-y2)*(y1-y3))
        c=1.0-yms/y-(3.0*astar/(2.0*y))*log(y/yms)-part1-part2-part3
        b=1.0-3.0/r+2.0*astar/(r**1.5)
        amytempreph=3.0*6.67e-8*1.989d33*m0*mdot0
        amytempreph=amytempreph/(8.0*pi*5.670367e-5*(r*rgcm)**3)
        amytempreph=amytempreph*(c/b+factor)
        amytempreph=amytempreph**0.25
        return
        end





