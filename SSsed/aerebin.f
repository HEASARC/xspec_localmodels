      SUBROUTINE aerebin(nin, inar, nout, istart, istop, fstart, fstop,
     &                  outar)
c		rashafer 1 april 1987
c	subroutine to perform a rebinning of an array, using simple
c	block interpolation.
c	note that many of the routines must be initialized by INIBIN
c
      INTEGER nin
c                                I: The dimension of the input array
      REAL inar(nin)
c                                I: input array
      INTEGER nout
c                                I: The dimension of the output array
      INTEGER istart(nout)
c                                I: The first bin
      INTEGER istop(nout)
c                                I: The last bin
      REAL fstart(nout)
c                                I: The fraction of the first bin
      REAL fstop(nout)
c                                I: The fraction of the last bin
      REAL outar(nout)
c                                R: The rebinned values.
c
      INTEGER i, j
c
      DO i = 1, nout
         IF (istart(i).GT.0) THEN
            outar(i) = fstart(i)*inar(istart(i))
            IF (istop(i).GT.istart(i)) THEN
               DO j = istart(i) + 1, istop(i) - 1
                  outar(i) = outar(i) + inar(j)
               ENDDO
               outar(i) = outar(i) + fstop(i)*inar(istop(i))
            ENDIF
         ELSE
            outar(i) = 0.
         ENDIF
      ENDDO
      RETURN
      END
