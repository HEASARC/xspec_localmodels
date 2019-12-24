!This file contains the convolution version of the thcomp model.
!Convolution method was implemented on the basis of the simplcut model
!by James F. Steiner (see http://jfsteiner.com/simplcut/)

subroutine thcompf(ear,ne,param,IFL,photar,photer) 
  IMPLICIT NONE
  integer ifl,ne,i,j

  real ear(0:ne),param(3),photar(ne),photer(ne),kte,tau,gam_tau,z_red
  double precision theta0,theta,alp,tauthin,tauthick,a,b,c,q,delta
  real thpar(3),tmparrt(ne)
  real psum,tsum

  ifl=0
  photer(1)=0.0

  gam_tau=param(1)
  kte=param(2)
  z_red=param(3)

  if(gam_tau.gt.0.0)then
     IF (gam_tau.eq.1.0) THEN
        gam_tau=1.001
     ENDIF

     alp=gam_tau-1.
     theta0=kte/511.
     theta=theta0*(1.+1.*theta0+3.*theta0**2.)
     tauthin=1.2/(1.+1.*theta0+5.*theta0**2.)
     tauthick=0.25/(1.+1.*theta0+3.*theta0**2.)
     
     a=tauthick*theta!>0
     b=tauthin*theta !>0
     c=-1./((1.5+alp)**2.-2.25)!<0
     delta=b**2.-4.*a*c
     q=-0.5*(b+sqrt(delta)) !<0 from numerical recipes: q=-0.5*(b+sgn(b)*sqrt(delta))
     tau=real(c/q) !>0 (square root, the second one, tau=q/a is <0)

!     write(*,*)"Gamma parametrization used. Tau=", tau
  else
     tau=-gam_tau
!     write(*,*)"Tau parametrization used"
  endif

  thpar(1) = tau
  thpar(2) = kte
  thpar(3) = z_red

  !     Initialize arrays
  !     NTHCOMP CALL

  do i=1,ne
     tmparrt(i) = photar(i)/(ear(i)-ear(i-1)) ! initialize  TO PHOTON DENSITY!
  end do

  call msrunthcomp(ear,ne,thpar,IFL,tmparrt)

  psum = 0.
  tsum = 0.
  do j=1,ne                 ! SUM THE PHOTONS  (*** INDEX RUN MAY BE TRUNCATED TO AVOID BAD LAST CHAN ****)  
     psum = psum+photar(j)
     tsum = tsum+tmparrt(j)
  enddo

  !     PHOTAR UNITS ARE:
  !     (Photons/cm^2/s per BIN) (i.e, flx integrated over the bin width)

  photar = tmparrt/tsum*psum ! normalize the output

  RETURN
END subroutine thcompf


SUBROUTINE msrunthcomp(Ear,Ne,Param,Ifl,Spec)
  !     driver for the Comptonization code solving Kompaneets equation
  !     seed photons - arbitrary XSPEC model

  !     number of model parameters: 3
  !     1: Thomson optical depth
  !     2: plasma temperature in keV
  !     3: redshift
  
  implicit none
  INTEGER Ne,Ifl
  REAL Param(*),Ear(0:Ne),Spec(Ne),seedspec(Ne),xear(Ne)

  INTEGER np,i,j,jl
  REAL xn , normfac , z_red
  REAL xth(0:Ne) , spt(0:Ne), xprim(0:Ne)

  LOGICAL IFLAG

  Ifl=0
  seedspec = spec

  z_red = param(3)
!  z_red = 0. ! previously param(4) in non-convolving thcomp

  np = 3
  xear = ear/511.

  call thcompton_fun(seedspec,param(2)/511.,param(1),xear,xth,Ne,spt,IFLAG)

  xn = (1+z_red)/511.            

  normfac=1.0

  !     zero final array
  do i=1,ne
     spec(i) = 0
  end do

  do i=0,ne
     xprim(i) = 0
  end do

  !     put primary into final array only if scale >= 0.

  if (IFLAG.or.(z_red.gt.0.0)) then
     j = 1
     do i=0,ne
        do while ((j .lt. ne).and.511.*xth(j) .lt. ear(i)*(1+z_red))
           j = j + 1
        end do
        if (j .lt. ne) then
           if (j .gt. 1) then
              jl = j - 1
              xprim(i) = spt(jl)+(ear(i)/511.*(1+z_red)-xth(jl))*(spt(jl+1)-spt(jl))/(xth(jl+1)-xth(jl))
           else
              xprim(i)=spt(1)
           endif
        else
           xprim(i)=spt(j)
        endif
     end do
  else
     xprim=spt
  endif

  xprim(Ne) = spt(Ne)

  do i=1,ne
     spec(i) = 0.5*(xprim(i)/ear(i)**2+xprim(i-1)/ear(i-1)**2)*(ear(i)-ear(i-1))*normfac
  end do

  RETURN
END SUBROUTINE msrunthcomp


subroutine thcompton_fun(specin,theta0,tau,xin,x,nen,sptot,IFLAG) 
  !     version: January 96
  !     
  !     Thermal Comptonization; solves Kompaneets eq. with some
  !     relativistic corrections. See Lightman & Zdziarski (1987), ApJ
  !     The seed spectrum is blackbody.
  IMPLICIT NONE
  integer nen
  real delta,xmin,xmax,deltal,xnr,xr,taukn,arg,flz,pi
  real sptot(0:nen),x(0:nen),xin(0:nen),dphesc(nen),dphdot(nen)
  real rel(nen),bet(nen),c0(nen),c1(nen),c2(nen),specin(nen)
  integer i,j,jl,jnr,jrel
  real(kind=8)::w,w1,z1,z2,z3,z4,z5,z6
  real tau,alpha,theta0,tauthin,tauthick,escape
  !     input parameters:
  real theta
  real fnorm,dif1
  logical IFLAG

  pi=3.1415927 
  !     clear arrays (important for repeated calls)
  sptot(0)=0.
  do j=1,nen
     dphesc(j)=0.
     dphdot(j)=0.
     rel(j)=0.
     bet(j)=0.
     c0(j)=0.
     c1(j)=0.
     c2(j)=0.
     sptot(j)=0.
  enddo

  !     a phenomenological relativistic correction to the temperature
  !     an increase of internal T at higher values
  theta=theta0*(1+1*theta0+3*theta0**2)

!  write(*,*)'internal kT=',theta*511

  !     delta is the 10-log interval of the photon array.

  xmin=xin(0)
  xmax=xin(nen)

  delta = log10(xmax/xmin)/(nen)
  deltal=delta*log(10.)

  !     
  !     X - ARRAY FOR PHOTON ENERGIES
  !     
  do j=0,nen
     x(j)=xmin*10.**(j*delta)
  enddo

  !     compute c0(x), and rel(x) arrays
  !     c2(x) is the relativistic correction to Kompaneets equation
  !     rel(x) is the Klein-Nishina cross section divided by the Thomson one
  do j=1,nen-1
     w=x(j)
     !     c2 is the correction to the diffusion coefficient calculated at w1
     !     w1 is x(j+1/2) (x(i) defined up to jmax+1)
     !     
     w1=sqrt(x(j)*x(j+1))
     c0(j)= real(w1**4/(1+4.6*w1+1.1*w1**2))
     c1(j)= theta
     c2(j)= 1
     if (w.le.0.05) then
        !     use asymptotic limit for rel(x) for x less than 0.05
        rel(j)=real(1-2*w+26*w*w/5)
     else
        z1=(1+w)/w**3
        z2=1+2*w
        z3=log(z2)
        z4=2*w*(1+w)/z2
        z5=z3/2/w
        z6=(1+3*w)/z2/z2
        rel(j)=real(0.75*(z1*(z4-z3)+z5-z6))
     end if
  enddo

  dif1 = 0.
  do i=0,nen
     dif1 = dif1+abs(xin(i)-x(i))
  enddo

  IFLAG = .TRUE.            ! INTERP NEEDED IF TRUE                                 
  if (dif1 .lt. 1e-3) then
     dphdot = specin
     IFLAG = .FALSE.        ! SAME GRID so NO INTERP                                     
  else                      ! interpolate to get a new array                        
     j=1
     do i=1,nen
        do while (j .le. (nen-1) .and. xin(j) .lt. x(i))
           j=j+1
        end do
        if (j .lt. nen) then
           if (j .gt. 1) then
              jl = j-1
              dphdot(i) = specin(jl)+(x(i)-xin(jl))*(specin(jl+1)-specin(jl))/(xin(jl+1)-xin(jl))
           else
              dphdot(i)=specin(1)
           endif
        else
           dphdot(i)=specin(j)
        endif
     end do
  endif

  dphdot(nen)=specin(nen)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !CCCCCCCCCCCCCCCCCCCCCCCCCccc

  !     
  !     compute beta array, the probalility of escape per Thomson time.
  !     bet evaluated for spherical geometry and nearly uniform sources.
  !     Between x=0.1 and 1.0, a function flz modifies beta to allow
  !     the increasingly large energy change per scattering to gradually
  !     eliminate spatial diffusion

  jnr=INT(log10(0.1/xmin)/delta+1)
  jnr=min(jnr,nen-1)
  !  in jrel epsilon_R changed from 1 to 2, AAZ 30.07.2019
  !  jrel=INT(log10(1/xmin)/delta+1)
  jrel=INT(log10(2/xmin)/delta+1)
  jrel=min(jrel,nen)
  xnr=x(jnr)
  xr=x(jrel)
  !     tauthin is a phenomenological factor for the optically-thin part of beta
  tauthin=1.2/(1+1*theta0+5*theta0**2)
  !     tauthick is a correction factor for the optically-thick part of beta
  !     with a phenomenological dependence on theta0. In the limit of low kT, 
  !     this correction factor should be 0.25
  tauthick=0.25/(1+1*theta0+3*theta0**2)
  do j=1, jnr-1
     taukn=tau*rel(j)
     bet(j)= 1/(tau*(tauthin+taukn*tauthick))
  end do
  do j=jnr,jrel
     taukn=tau*rel(j)
     arg=(x(j)-xnr)/(xr-xnr)
     flz=1-arg
     bet(j)= 1/(tau*(tauthin+taukn*tauthick*flz))
  end do
  do j=jrel+1,nen
     taukn=tau*rel(j)
     bet(j)=1/(tau*tauthin)
  end do

  !     Note that tau in thermlc affects only the overall norm, and it is 
  !     actually not needed at all.

  call mythermlc(1.,theta,deltal,x,nen,dphesc,dphdot,bet,c0,c1,c2)
  !     the unscattered fraction in a uniform sphere with uniform 
  !     production of photons 
  escape=(2*tau**2 - 1 + exp(-2*tau)*(2*tau+1))*3/(8*tau**3) 
  !     the geometrical average with escape for central production of photons 
  escape=Sqrt(escape*Exp(-tau))
 ! write(*,*) 'escaped fraction=',escape

  !  alpha=sqrt(2.25+bet(1)/theta)-1.5
  alpha=sqrt(2.25+1/(theta*(tau*tauthin+(tau**2.)*tauthick)))-1.5 !to get agreement between alpha computed from rel(1)=1
!  write(*,*) 'alpha = ',alpha

  !     the spectrum in E F_E
  do j=1,nen-1
     !     add the unscattered fraction
     dphesc(j)=dphesc(j)+escape*dphdot(j)
     sptot(j)=dphesc(j)*x(j)**2
  enddo
  !     normalize the spectrum to unit number of photons.
  !     sptot is nu F_nu
  fnorm = 0
  do j = 1, nen-1
     fnorm = fnorm + sptot(j)/x(j)
  enddo
  fnorm=fnorm*log(x(2)/x(1))
  do j = 1, nen-1
     sptot(j) = sptot(j)/fnorm
  enddo

  return  
end subroutine thcompton_fun


subroutine mythermlc(tau,theta,deltal,x,nen,dphesc,dphdot,bet,c0,c1,c2)
  !     This program computes the effects of Comptonization by
  !     nonrelativistic thermal electrons in a sphere including escape, and
  !     relativistic corrections up to photon energies of 1 MeV.
  !     the dimensionless photon energy is x=hv/(m*c*c)
  !     
  !     The input parameters and functions are:
  !     dphdot(x), the photon production rate 
  !     tau, the Thomson scattering depth
  !     theta, the temperature in units of m*c*c
  !     c2(x), and bet(x), the coefficients in the K-equation and the
  !     probability of photon escape per Thomson time, respectively,
  !     including Klein-Nishina corrections
  !     The output parameters and functions are:
  !     dphesc(x), the escaping photon density
  implicit none 
  integer nen
  real tau,theta,deltal
  real x(0:nen),dphesc(nen),dphdot(nen),bet(nen),c0(nen),c1(nen)
  integer j,jj
  real a(nen),b(nen),c(nen),c20,w1,w2,t1,t2,t3,x32,aa,c2(nen) 
  real d(nen),alp(nen),g(nen),gam(nen),u(nen)

  !     u(x) is the dimensionless photon occupation number
  c20=tau/deltal
  !     
  !     determine u
  !     define coefficients going into equation
  !     a(j)*u(j+1)+b(j)*u(j)+c(j)*u(j-1)=d(j)
  do j=2,nen-1
     w1=sqrt(x(j)*x(j+1))
     w2=sqrt(x(j-1)*x(j))
     !     w1 is x(j+1/2)
     !     w2 is x(j-1/2)
     a(j)=-c20*c0(j)*(c1(j)/deltal/w1+0.5*c2(j))
     t1=-c20*c0(j)*(0.5*c2(j)-c1(j)/deltal/w1)
     t2=c20*c0(j-1)*(c1(j-1)/deltal/w2+0.5*c2(j-1))
     t3=x(j)**3*(tau*bet(j))
     b(j)=t1+t2+t3
     c(j)=c20*c0(j-1)*(0.5*c2(j-1)-c1(j-1)/deltal/w2)
     d(j)=x(j)*dphdot(j)
  enddo

  !     define constants going into boundary terms
  !     u(1)=aa*u(2) (zero flux at lowest energy)
  !     u(jx2) given from region 2 above
  x32=sqrt(x(1)*x(2))
  aa=(theta/deltal/x32+0.5)/(theta/deltal/x32-0.5)
  !     
  !     zero flux at the highest energy
  u(nen)=0
  !     
  !     invert tridiagonal matrix
  alp(2)=b(2)+c(2)*aa
  gam(2)=a(2)/alp(2)
  do j=3,nen-1
     alp(j)=b(j)-c(j)*gam(j-1)
     gam(j)=a(j)/alp(j)
  enddo
  g(2)=d(2)/alp(2)
  do j=3,nen-2
     g(j)=(d(j)-c(j)*g(j-1))/alp(j)
  enddo
  g(nen-1)=(d(nen-1)-a(nen-1)*u(nen)-c(nen-1)*g(nen-2))/alp(nen-1)
  u(nen-1)=g(nen-1)
  do j=3,nen-1
     jj=nen+1-j
     u(jj)=g(jj)-gam(jj)*u(jj+1)
  enddo
  u(1)=aa*u(2)
  
  !     compute new value of dph(x) and new value of dphesc(x)

  do j=1,nen-1
     dphesc(j)=x(j)*x(j)*u(j)*bet(j)*tau
  enddo

  return
end subroutine mythermlc

