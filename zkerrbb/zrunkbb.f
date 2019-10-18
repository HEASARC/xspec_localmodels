
      SUBROUTINE ZRUNKBB(EAR, NE, ETA, ASTAR, THETA, MBH, MDD, DBH, 
     &                   FCOL, RFLAG, LFLAG, ZBH, PHOTAR, FLUXE, NEX1, 
     &                   NEX2)

      IMPLICIT NONE

      Integer NE
      Real ETA, ASTAR, THETA, MBH, MDD, DBH, FCOL, ZBH
      Real EAR(0:NE), PHOTAR(NE), FLUXE(NE,2)
      Integer RFLAG, LFLAG
      Integer nex1(NE,2), nex2(NE,2)

C
C     Written by: Li-Xin Li
C              Harvard-Smithsonian Center for Astrophysics
C              60 Garden St., Cambridge, MA 02138
C                      
C     Version:         October, 2004
C
C     Reference:  Li, Zimmerman, Narayan, & McClintock 2004, astro-ph/0411583
C
C     The program computes the black-body emission spectrum from a thin 
C     Keplerian accretion disk around a Kerr black hole. All relativistic 
C     effects are taken into account, including frame-dragging, Doppler boost, 
C     gravitational redshift, and bending of light by the gravity of the black 
C     hole. In particular, self-irradiation of the disk as a result of light 
C     deflection is included. The inner boundary of the disk is fixed at the 
C     marginally stable orbit. Torque at the inner boundary of the disk 
C     is allowed to be nonzero. However, when this program is applied to data
C     reduction, a zero torque (eta = 0) is recommended since we have found
C     that the effect of a nonzero torque on the spectrum can, to a good 
C     approximation, be absorbed into a zero torque model by adjusting the
C     mass accretion rate and the normalization.
C
C     The program reads in the precalculated data file kerrbb.fits and 
C     uses it to fit spectral data by linear interpolation. The first
C     extension (GRIDVALS) contains 4 rows specifying:
C     - a grid of black hole spin (a/M);
C     - a grid of disk inclination angles (in degrees);
C     - a grid of the torque at the inner boundary of the disk (eta, dimensionless);
C     - a grid of spectrum energy (in keV)
C     The second extension (FLUX) contains the flux density at each grid 
C     point specified by the above four parameters, in photons/keV/cm^2/sec.
C     [There are four columns in this file: FLUX1, when both self-irradiation 
C     and limb-darkening are turned off; FLUX2, when self-irradiation is off 
C     but limb-darkening is on; FLUX3, when self-irradiation is on but 
C     limb-darkening is off; FLUX4, when both self-irradiation and 
C     limb-darkening are on.]
C
C     The program requires the following input arguments
C     - "EAR(0:NE)": array of observed energy (frequency) bins (in units of 
C        keV);
C     - "NE": number of observed energy bins (size of the energy array);
C     - "eta" : ratio of the disk power produced by a
C         torque at the disk inner boundary to the disk power arising 
C         from accretion. It must be >= 0 and <=1. When eta = 0, the 
C         solution corresponds to that of a standard Keplerian disk with 
C         zero torque at the inner boundary;
C     - "astar" : specific angular momentum of the black hole in 
C        units of the black hole mass M (geometrized units G=c=1). a 
C         should be >= -1 and < 1;
C     - "theta" : disk's inclination angle (the angle between the 
C        axis of the disk and the line of sight). It is expressed in 
C         degrees. i=0 is for a "face-on" accretion disk. i should be <=  
C         85 degree;
C     - "Mbh" : the mass of the black hole in units of the solar mass;
C     - "Mdd" : the "effective" mass accretion rate of the disk 
C        in units of 10^18 g/sec. When eta = 0 (zero torque at the inner 
C        boundary), this is just the mass accretion rate of the disk. When 
C         eta is nonzero, the effective mass accretion rate = (1+eta) times 
C         the true mass accretion rate of the disk. The total disk 
C         luminosity is then "epsilon" times "the effective mass accretion 
C         rate" times "c^2", where epsilon is the radiation efficiency of a 
C         standard accretion disk around the Kerr black hole;
C     - "Dbh" : the distance from the observer to the black 
C         hole in units of kpc;
C     - "fcol" : spectral hardening factor, T_col/T_eff. It should
C         be greater than 1.0, and considered to be 1.5-1.9 for accretion 
C         disks around a stellar-mass black hole. See, e.g., Shimura and 
C         Takahara 1995, ApJ, 445, 780;
C     - "rflag" : a flag to switch on/off the effect of 
C         self-irradiation (never allowed to be free). Self-irradiation is 
C         included when rflag is > 0. Self-irradiation is not included when 
C         rflag is <= 0;
C     - "lflag" : a flag to switch on/off the effect of limb-
C         darkening (never allowed to be free). The disk emission is assumed 
C         to be limb-darkened when lflag is > 0. The disk emission is assumed 
C         to be isotropic when lflag is <= 0.
C     - "ZBH": the redshift of the source. This is only used to modify the
C         variable ca (change by Jack Steiner). The redshift will presumably
C         be appropriate for the distance DBH input
C
C
C     The program outputs an array of the observed photon flux in each bin 
C     of energy: "PHOTAR(NE)", in units of photons/cm^2/second. For
C     example, PHOTAR(i) gives the observed photon number flux with 
C     observed photon energy between EAR(i) and EAR(i+1).

      INTEGER NG, NS, NTH, NENER
      PARAMETER(NG=6, NS=46, NTH=18, NENER=601)

      Real FLUX0(NS,NTH,NG,NENER), FLUX(NENER), E(NENER), GI(NG) 
      Real AI(NS), THETAI(NTH), ear0(NENER)

      Integer IREAD
      Integer lflagsav, rflagsav
      Real difx, dify, difz, diftx, difty, diftz, dx, dy, dz, 
     $     dife, difte, lslope
      Real ca, cb, cc, t, s, w, tc, sc, wc
      Real de, earm
      Real cflux1, cflux2, cflux3, cflux4, cflux5, cflux6, cflux7, 
     $     cflux8, iflux1, iflux2, iflux12, iflux0
      Integer nx, ny, nz, nxb, nyb, nzb, nx1, nx2, ny1, ny2, nz1, nz2
      Integer nex, nexb
      Integer i, j, k, n1, n2, n3, n4
      Integer ilun, block, status, hdutyp, icol

      Logical qanyf

      character(256) fgmodf, datdir, contxt
      character(256) filenm
      integer lenact
      external lenact, fgmodf

      save iread, lflagsav, rflagsav
      data iread/0/

      ca = Mdd**0.25*Mbh**(-0.5)*fcol
      cb = Mbh*Mbh/(Dbh*Dbh)*fcol**(-4.)
      cc = cb*ca**2.

C Modification from Jack Steiner for extragalactic sources
      ca = ca / (1.0+zbh)
      
      if ((iread .ne. 99) .or. (lflagsav .ne. lflag) .or. 
     $       (rflagsav .ne. rflag)) then
         datdir = fgmodf()
         filenm = datdir(:lenact(datdir))//'kerrbb.fits'

         status = 0

C Open FITS input file

         CALL getlun(ilun)
         CALL ftopen(ilun, filenm, 0, block, status)
         contxt = 'Failed to open '//filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Move to first (GRIDVALS) extension

         CALL ftmrhd(ilun, 1, hdutyp, status)
         contxt = 'Failed to move to first extension of '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Read four sets of grid values

         CALL ftgcve(ilun, 1, 1, 1, ng, 0.0, gi, qanyf, status)
         CALL ftgcve(ilun, 1, 2, 1, ns, 0.0, ai, qanyf, status)
         CALL ftgcve(ilun, 1, 3, 1, nth, 0.0, thetai, qanyf, status)
         CALL ftgcve(ilun, 1, 4, 1, nener, 0.0, e, qanyf, status)
         contxt = 'Failed to read GRIDVALS data from '
     &            //filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Move to second (FLUX) extension
          
         CALL ftmrhd(ilun, 1, hdutyp, status)
         contxt = 'Failed to move to second extension of '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

C Read the appropriate flux data
         
         if ((rflag .le. 0).and.(lflag .le. 0)) then
            icol = 1
         else  if ((rflag .le. 0).and.(lflag .gt. 0)) then
            icol = 2
         else  if ((rflag .gt. 0).and.(lflag .le. 0)) then
            icol = 3
         else
            icol = 4
         endif

         CALL ftgcve(ilun, icol, 1, 1, ng*ns*nth*nener, 0.0, 
     &               flux0, qanyf, status)
         contxt = 'Failed to read FLUX data from '//
     &            filenm(:lenact(filenm))
         IF ( status .NE. 0 ) GOTO 999

       CALL ftclos(ilun, status)
         CALL frelun(ilun)
                
       iread = 99
         lflagsav = lflag
         rflagsav = rflag

 999     CONTINUE
         IF ( status .NE. 0 ) THEN
            CALL xwrite(contxt, 10)
            WRITE(contxt, '(a,i6)') 'RUNKBB: Status = ', status
            CALL xwrite(contxt, 10)
         ENDIF
      
      endif

      nx = 1
      difx = abs(astar-ai(1))
      do i = 1,ns
       diftx = abs(astar-ai(i))
       if (diftx .lt. difx) then
          nx = i
          difx = diftx
       end if
      enddo

      ny = 1
      dify = abs(theta-thetai(1))
      do i = 1,nth
       difty = abs(theta-thetai(i))
       if (difty .lt. dify) then
          ny = i
          dify = difty
       end if
      enddo

      nz = 1
      difz = abs(eta-gi(1))
      do i = 1,ng
       diftz = abs(eta-gi(i))
       if (diftz .lt. difz) then
          nz = i
          difz = diftz
       end if
      enddo

      dx = astar - ai(nx)
      dy = theta-thetai(ny)
      dz = eta - gi(nz)
      nxb = nx + INT(sign(1.,dx))
      nyb = ny + INT(sign(1.,dy))
      nzb = nz + INT(sign(1.,dz))

      nx1 = min(nx,nxb)
      nx2 = max(nx,nxb)
      ny1 = min(ny,nyb)
      ny2 = max(ny,nyb)
      nz1 = min(nz,nzb)
      nz2 = max(nz,nzb)

c Trap special case of values being at top of ranges

      IF ( nx2 .GT. NS ) THEN
         nx2 = NS
         nx1 = nx2 - 1
      ENDIF
      IF ( ny2 .GT. NTH ) THEN
         ny2 = NTH
         ny1 = ny2 - 1
      ENDIF
      IF ( nz2 .GT. NG ) THEN
         nz2 = NG
         nz1 = nz2 - 1
      ENDIF


      t = (astar-ai(nx1))/(ai(nx2)-ai(nx1))
      tc = 1. - t
      s = (theta-thetai(ny1))/(thetai(ny2)-thetai(ny1))
      sc = 1. - s
      w = (eta-gi(nz1))/(gi(nz2)-gi(nz1))
      wc = 1. - w

      do i = 1,nener
       ear0(i) = e(i)*ca

       cflux1 = tc*sc*wc*flux0(nx1,ny1,nz1,i)
       cflux2 = t*sc*wc*flux0(nx2,ny1,nz1,i)
       cflux3 = t*s*wc*flux0(nx2,ny2,nz1,i)
       cflux4 = tc*s*wc*flux0(nx1,ny2,nz1,i)
       cflux5 = tc*sc*w*flux0(nx1,ny1,nz2,i)
       cflux6 = t*sc*w*flux0(nx2,ny1,nz2,i)
       cflux7 = t*s*w*flux0(nx2,ny2,nz2,i)
       cflux8 = tc*s*w*flux0(nx1,ny2,nz2,i)

       flux(i) = cc*(cflux1 + cflux2 + cflux3 + cflux4 
     $           + cflux5 + cflux6 + cflux7 + cflux8)
            
      enddo

      do i = 1,NE
         do k=1,2
            earm =ear(i+k-2)

C     When bin energy earm < ear0(1), calculate the specific flux at
C     earm according to law: specific flux proportional to energy^(-2/3)

            if (earm .lt. ear0(1)) then
               fluxe(i,k) = flux(1)*(ear0(1)/earm)**(2./3.)
               nex1(i,k) = 1
               nex2(i,k) = 1

C     When bin energy earm > ear0(nener), cutoff the flux density

            else if (earm .gt. ear0(nener)) then
               fluxe(i,k) = 1.e-20
               nex1(i,k) = nener
               nex2(i,k) = nener

C     Use linear interpolation when ear0(1) < earm < ear0(nener)

            else
               nex = 1
               dife = abs(earm-ear0(1))
               do j = 1,nener
                  difte = abs(earm-ear0(j))
                  if (difte .lt. dife) then
                     nex = j
                     dife = difte
                  end if
               enddo
         
               de = earm - ear0(nex)
               nexb = nex + INT(sign(1.,de))
               nex1(i,k) = min(nex,nexb)
               nex2(i,k) = max(nex,nexb)

               lslope = (log10(flux(nex2(i,k)))-log10(flux(nex1(i,k))))
     $               /(log10(ear0(nex2(i,k)))-log10(ear0(nex1(i,k))))
               fluxe(i,k) = flux(nex1(i,k))*10.**(lslope*(log10(earm)
     $               -log10(ear0(nex1(i,k)))))
            endif
         enddo

         n1=nex1(i,1)
         n2=nex2(i,1)
         n3=nex1(i,2)
         n4=nex2(i,2)

         iflux1=0.5*(fluxe(i,1)+flux(n2))*(ear0(n2)-ear(i-1))
         iflux2=0.5*(fluxe(i,2)+flux(n3))*(ear(i)-ear0(n3))
         iflux12=0.5*(fluxe(i,1)+fluxe(i,2))*(ear(i)-ear(i-1))

         if (n3 .gt. n2) then
            iflux0=0.0
            do k=n2,n3-1
               iflux0=iflux0+0.5*(flux(k)+flux(k+1))*(ear0(k+1)-ear0(k))
            enddo

            PHOTAR(i) = iflux0 + iflux1 + iflux2

         else if (n3 .eq. n2) then
            PHOTAR(i) = iflux1 + iflux2

         else
            PHOTAR(i) = iflux12

         endif
      
      enddo
      
      RETURN
      END

