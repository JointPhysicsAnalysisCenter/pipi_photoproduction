c------------------------------------------------------------------------
      SUBROUTINE UPDOWN (N,XI,FI,X,F,DF)
      implicit double precision(a-h,o-z)
C
C Subroutine performing the Lagrange interpolation with the
C upward and downward correction method.  F: interpolated
C value.  DF: error estimated.
C Copyright (c) Tao Pang 1997.
C
      PARAMETER (NMAX=27) 
      DIMENSION XI(N),FI(N),DP(NMAX,NMAX),DM(NMAX,NMAX)
C
      IF (N.GT.NMAX) STOP 'Dimension too large'
      DX = ABS(XI(N)-XI(1))
      DO     100 I = 1, N
        DP(I,I) = FI(I)
        DM(I,I) = FI(I)
        DXT = ABS(X-XI(I))
        IF (DXT.LT.DX) THEN
           I0 = I
           DX = DXT
        ELSE
        END IF
  100 CONTINUE
      J0 = I0
C
C Evaluate correction matrices
C
      DO     200 I = 1, N-1
        DO   150 J = 1, N-I
          K = J+I
          DT =(DP(J,K-1)-DM(J+1,K))/(XI(K)-XI(J))
          DP(J,K) = DT*(XI(K)-X)
          DM(J,K) = DT*(XI(J)-X)
  150   CONTINUE
  200 CONTINUE
C
C Update the approximation
C
      F = FI(I0)
      IT = 0
      IF(X.LT.XI(I0)) IT = 1
      DO     300 I = 1, N-1
        IF ((IT.EQ.1).OR.(J0.EQ.N)) THEN
          I0 = I0 - 1
          DF = DP(I0,J0)
          F  = F + DF
          IT = 0
          IF (J0.EQ.N) IT = 1
        ELSE IF ((IT.EQ.0).OR.(I0.EQ.1)) THEN
          J0 = J0 + 1
          DF = DM(I0,J0)
          F = F + DF
          IT = 1
          IF (I0.EQ.1) IT = 0
        END IF
  300 CONTINUE
      DF = ABS(DF)
      RETURN
      END
! Interpolacja motodą funkcji sklejanych
c------------------------------------------------------------------------
      SUBROUTINE ZSPL3(N,T,Y,H,B,U,V,Z)
      implicit double precision(a-h,o-z)
      DIMENSION  T(N),Y(N),H(N),B(N),U(N),V(N),Z(N)       
      DO 2 I = 1,N-1
        H(I) = T(I+1) - T(I)
        B(I) = (Y(I+1) -Y(I))/H(I)    
   2  CONTINUE    
      U(2) = 2.0*(H(1) + H(2))
      V(2) = 6.0*(B(2) - B(1))
      DO 3 I = 3,N-1
        U(I) = 2.0*(H(I) + H(I-1)) - H(I-1)**2/U(I-1)     
        V(I) = 6.0*(B(I) - B(I-1)) - H(I-1)*V(I-1)/U(I-1) 
   3  CONTINUE    
      Z(N) = 0.0  
      DO 4 I = N-1,2,-1     
        Z(I) = (V(I) - H(I)*Z(I+1))/U(I)
   4  CONTINUE    
      Z(1) = 0.0
      RETURN
      END 
c------------------------------------------------------------------------
      FUNCTION SPL3(N,T,Y,Z,X)
      implicit double precision(a-h,o-z)
      DIMENSION  T(N),Y(N),Z(N)       
      DO 2 I = N-1,1,-1     
        DIFF = X - T(I)     
        IF(DIFF .GE. 0.0)  GO TO 3    
   2  CONTINUE    
      I = 1 
   3  H = T(I+1) - T(I)     
      B = (Y(I+1) - Y(I))/H - H*(Z(I+1) + 2.0*Z(I))/6.0 
      P = 0.5*Z(I) + DIFF*(Z(I+1) - Z(I))/(6.0*H) 
      P = B + DIFF*P
      SPL3 = Y(I) + DIFF*P  
      RETURN      
      END 
c------------------------------------------------------------------------
      complex*16 function funkcjakulista(l,m,costh,phi)
      implicit double precision (a-h,o-z)
      complex*16 CI
      PI=3.141592653589793116d0
      CI=DCMPLX(0.D0,1.D0)
      sinth=DSQRT(1.d0-costh**2)
      if (l.eq.0) then
        funkcjakulista=1.d0/(2.d0*DSQRT(Pi))
      else if (l.eq.1) then
        if(m.eq.-1) then
        funkcjakulista=DSQRT(3.d0/(2.d0*Pi))*(CDEXP(-CI*phi)
     &  *sinth)/2.d0
        else if(m.eq.0) then
        funkcjakulista=(DSQRT(3.d0/Pi)*costh)/2.d0
        else if (m.eq.1) then
        funkcjakulista=-DSQRT(3.d0/(2.d0*Pi))*(CDEXP(CI*phi)
     &  *sinth)/2.d0
        endif
      else if (l.eq.2) then
        if (m.eq.-2) then
        funkcjakulista=(CDEXP(-2.d0*CI*phi)*DSQRT(15.d0/(2.d0*Pi))
     &  *sinth**2)/4.d0
        else if (m.eq.-1) then
        funkcjakulista=(CDEXP(-CI*phi)*DSQRT(15.d0/(2.d0*Pi))*
     &  costh*sinth)/2.d0
        else if (m.eq.0) then
        funkcjakulista=(DSQRT(5.d0/Pi)*(-1.d0+3.d0*costh**2))/
     &  4.d0
        else if (m.eq.1) then
        funkcjakulista=-(CDEXP(CI*phi)*DSQRT(15.d0/(2.d0*Pi))*
     &  costh*sinth)/2.d0
        else if (m.eq.2) then
        funkcjakulista=(CDEXP(2.d0*CI*phi)*DSQRT(15.d0/(2.d0*Pi))
     &  *sinth**2)/4.d0
        endif
      else if (l.eq.3) then
        if (m.eq.-3) then
        funkcjakulista=(CDEXP(-3.d0*CI*phi)*DSQRT(35.d0/Pi)*sinth**3)
     &  /8.d0
        else if (m.eq.-2) then
        funkcjakulista=(CDEXP(-2.d0*CI*phi)*DSQRT(105.d0/(2.d0*PI))*
     &  costh*sinth**2)/4.d0
        else if (m.eq.-1) then
        funkcjakulista=CDEXP(-CI*phi)*DSQRT(21.d0/Pi)*
     &  (-1.d0+5.d0*costh**2)*sinth/8.d0
        else if (m.eq.0) then
        funkcjakulista=(DSQRT(7.d0/Pi)*(-3.d0*costh+5.d0*
     &  costh**3))/4.d0
        else if (m.eq.1) then
        funkcjakulista=-(CDEXP(CI*phi)*DSQRT(21.d0/Pi)*
     &  (-1.d0+5.d0*costh**2)*sinth)/8.d0
        else if (m.eq.2) then
        funkcjakulista=(CDEXP(2.d0*CI*phi)*DSQRT(105.d0/(2.d0*Pi))
     &  *costh*sinth**2)/4.d0
        else if (m.eq.3) then
        funkcjakulista=(-CDEXP(3.d0*CI*phi)*DSQRT(35.d0/Pi)*sinth**3)
     &  /8.d0
        endif
      else if (l.eq.4) then
        if (m.eq.-4) then
        funkcjakulista=3.d0/16.d0*CDEXP(-4.d0*CI*phi)*
     &  DSQRT(35.d0/(2.d0*PI))*sinth**4
        else if (m.eq.-3) then
        funkcjakulista=3.d0/8.d0*CDEXP(-3.d0*CI*phi)*DSQRT(35.d0/Pi)*
     &  costh*sinth**3
        else if (m.eq.-2) then
        funkcjakulista=3.d0/8.d0*CDEXP(-2.d0*CI*phi)*
     &  DSQRT(5.d0/(2.d0*Pi))*(-1.d0+7.d0*costh**2)*sinth**2
        else if (m.eq.-1) then
        funkcjakulista=3.d0/8.d0*CDEXP(-CI*phi)*DSQRT(5.d0/Pi)*costh*
     &  (-3.d0+7.d0*costh**2)*sinth
        else if (m.eq.0) then
        funkcjakulista=(3.d0*(3.d0-30.d0*costh**2+35.d0*
     &  costh**4))/(16.d0*DSQRT(Pi))
        else if (m.eq.1) then
        funkcjakulista=(-3.d0*CDEXP(CI*phi)*DSQRT(5.d0/Pi)*costh*
     &  (-3.d0+7.d0*costh**2)*sinth)/8.d0
        else if (m.eq.2) then
        funkcjakulista=(3.d0*CDEXP(2.d0*CI*phi)*DSQRT(5.d0/(2.d0*Pi))*
     &  (-1.d0+7.d0*costh**2)*sinth**2)/8.d0
        else if (m.eq.3) then
        funkcjakulista=-3.d0/8.d0*CDEXP(3.d0*CI*phi)*DSQRT(35.d0/Pi)*
     &  costh*sinth**3
        else if (m.eq.4) then
        funkcjakulista=3.d0/16.d0*CDEXP(4.d0*CI*phi)*
     &  DSQRT(35.d0/(2.d0*Pi))*sinth**4
        endif
      else if (l.eq.5) then
        if(m.eq.-5) then
        funkcjakulista=3.d0/32.d0*CDEXP(-5.d0*CI*phi)*DSQRT(77.d0/Pi)*
     &  sinth**5
        else if(m.eq.-4) then
        funkcjakulista=3.d0/16.d0*CDEXP(-4.d0*CI*phi)*
     &  DSQRT(385.d0/(2.d0*Pi))*costh*sinth**4
        else if(m.eq.-3) then
        funkcjakulista=1.d0/32.d0*CDEXP(-3.d0*CI*phi)*DSQRT(385.d0/Pi)*
     &  (-1.d0+9.d0*costh**2)*sinth**3
        else if(m.eq.-2) then
        funkcjakulista=1.d0/8.d0*CDEXP(-2.d0*CI*phi)*
     &  DSQRT(1155.d0/(2.d0*Pi))*costh*(-1.d0+3.d0*costh**2)*sinth**2
        else if(m.eq.-1) then
        funkcjakulista=1.d0/16.d0*CDEXP(-CI*phi)*DSQRT(165.d0/(2.d0*Pi))
     &  *(1.d0-14.d0*costh**2+21.d0*costh**4)*sinth
        else if(m.eq.0) then
        funkcjakulista=1.d0/16.d0*DSQRT(11.d0/Pi)*
     &  (15.d0*costh-70.d0*costh**3+63.d0*costh**5)
        else if(m.eq.1) then
        funkcjakulista=-1.d0/16.d0*CDEXP(CI*phi)*
     &  DSQRT(165.d0/(2.d0*Pi))*(1.d0-14.d0*costh**2+21.d0*costh**4)*
     &  sinth
        else if(m.eq.2) then
        funkcjakulista=1.d0/8.d0*CDEXP(2.d0*CI*phi)*
     &  DSQRT(1155.d0/(2.d0*Pi))*costh*(-1.d0+3.d0*costh**2)*sinth**2
        else if(m.eq.3) then
        funkcjakulista=-1.d0/32.d0*CDEXP(3.d0*CI*phi)*DSQRT(385.d0/Pi)*
     &  (-1.d0+9.d0*costh**2)*sinth**3
        else if(m.eq.4) then
        funkcjakulista=3.d0/16.d0*CDEXP(4.d0*CI*phi)*
     &  DSQRT(385.d0/(2.d0*Pi))*costh*sinth**4
        else if(m.eq.5) then
        funkcjakulista=-3.d0/32.d0*CDEXP(5.d0*CI*phi)*
     &  DSQRT(77.d0/Pi)*sinth**5
        endif
      else if(l==6) then
        if(m==-6) then
        funkcjakulista=1.d0/64.d0*dsqrt(3003.d0/pi)*cdexp(-6.d0*ci*phi)*
     &  sinth**6
        else if(m==-5) then
        funkcjakulista=3.d0/32.d0*dsqrt(1001.d0/pi)*cdexp(-5.d0*ci*phi)*
     &  sinth**5*costh
        else if(m==-4) then
        funkcjakulista=3.d0/32.d0*dsqrt(91.d0/2.d0/pi)*
     &  cdexp(-4.d0*ci*phi)*sinth**4*(11.d0*costh**2-1.d0)
        else if(m==-3) then
        funkcjakulista=dsqrt(1365.d0/pi)/32.d0*cdexp(-3.d0*ci*phi)*
     &  sinth**3*costh*(11.d0*costh**2-3.d0)
        else if(m==-2) then
        funkcjakulista=dsqrt(1365.d0/pi)/64.d0*cdexp(-2.d0*ci*phi)*
     &  sinth**2*(33.d0*costh**4-18.d0*costh**2+1.d0)
        else if(m==-1) then
        funkcjakulista=dsqrt(273.d0/2.d0/pi)/16.d0*cdexp(-ci*phi)*
     &  sinth*costh*(33.d0*costh**4-30.d0*costh**2+5.d0)
        else if(m==0) then
        funkcjakulista=dsqrt(13.d0/pi)/32.d0*(231.d0*costh**6-
     &  315.d0*costh**4+105.d0*costh**2-5.d0)
        else if(m==1) then
        funkcjakulista=-dsqrt(273.d0/2.d0/pi)/16.d0*cdexp(ci*phi)*
     &  sinth*costh*(33.d0*costh**4-30.d0*costh**2+5.d0)
        else if(m==2) then
        funkcjakulista=dsqrt(1365.d0/pi)/64.d0*cdexp(2.d0*ci*phi)*
     &  sinth**2*(33.d0*costh**4-18.d0*costh**2+1.d0)
        else if(m==3) then
        funkcjakulista=-dsqrt(1365.d0/pi)/32.d0*cdexp(3.d0*ci*phi)*
     &  sinth**3*costh*(11.d0*costh**2-3.d0)
        else if(m==4) then
        funkcjakulista=3.d0/32.d0*dsqrt(91.d0/2.d0/pi)*
     &  cdexp(4.d0*ci*phi)*sinth**4*(11.d0*costh**2-1.d0)
        else if(m==5) then
        funkcjakulista=-3.d0/32.d0*dsqrt(1001.d0/pi)*cdexp(5.d0*ci*phi)*
     &  sinth**5*costh
        else if(m==6) then
        funkcjakulista=1.d0/64.d0*dsqrt(3003.d0/pi)*cdexp(6.d0*ci*phi)*
     &  sinth**6
        endif
      else if(l==7) then
        if(m==-7) then
        funkcjakulista=3.d0/64.d0*sqrt(715.d0/(2.d0*pi))*
     &  exp(-7.d0*ci*phi)*sinth**7
        else if(m==-6) then
        funkcjakulista=3.d0/64.d0*sqrt(5005.d0/pi)*exp(-6.d0*ci*phi)
     &  *sinth**6*costh
        else if(m==-5) then
        funkcjakulista=3.d0/64.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(-5.d0*ci*phi)*sinth**5*(13.d0*costh**2 - 1.d0)
        else if(m==-4) then
        funkcjakulista=3.d0/32.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(-4.d0*ci*phi)*sinth**4*costh*(13.d0*costh**2 - 3.d0)
        else if(m==-3) then
        funkcjakulista=3.d0/64.d0*sqrt(35.d0/(2.d0*pi))*
     &  exp(-3.d0*ci*phi)*sinth**3*
     &  (143.d0*costh**4-66.d0*costh**2+3.d0)
        else if(m==-2) then
        funkcjakulista=3.d0/64.d0*sqrt(35.d0/pi)*exp(-2.d0*ci*phi)
     &  *sinth**2*costh*(143.d0*costh**4 - 110.d0*costh**2 + 15.d0)
        else if(m==-1) then
        funkcjakulista=1.d0/64.d0*sqrt(105.d0/(2.d0*pi))*exp(-ci*phi)*
     &  sinth*(429.d0*costh**6-495.d0*costh**4+135.d0*costh**2 - 5.d0)
        else if(m==0) then
        funkcjakulista=1.d0/32.d0*sqrt(15.d0/pi)*
     &  (429.d0*costh**7-693.d0*costh**5 +315.d0*costh**3-35.d0*costh)
        else if(m==1) then
        funkcjakulista=-1.d0/64.d0*sqrt(105.d0/(2.d0*pi))*exp(ci*phi)*
     &  sinth*(429.d0*costh**6 - 495.d0*costh**4+135.d0*costh**2-5.d0)
        else if(m==2) then
        funkcjakulista=3.d0/64.d0*sqrt(35.d0/pi)*exp(2.d0*ci*phi)*
     &  sinth**2*costh*(143.d0*costh**4 - 110.d0*costh**2 + 15.d0)
        else if(m==3) then
        funkcjakulista=-3.d0/64.d0*sqrt(35.d0/(2.d0*pi))*
     &  exp(3.d0*ci*phi)*sinth**3*(143.d0*costh**4-66.d0*costh**2+3.d0)
        else if(m==4) then
        funkcjakulista=3.d0/32.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(4.d0*ci*phi)*sinth**4*costh*(13.d0*costh**2 - 3.d0)
        else if(m==5) then
        funkcjakulista=-3.d0/64.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(5.d0*ci*phi)*sinth**5*(13.d0*costh**2 - 1.d0)
        else if(m==6) then
        funkcjakulista=3.d0/64.d0*sqrt(5005.d0/pi)*exp(6.d0*ci*phi)*
     &  sinth**6*costh
        else if(m==7) then
        funkcjakulista=-3.d0/64.d0*sqrt(715.d0/(2.d0*pi))*
     &  exp(7.d0*ci*phi)*sinth**7
        endif
      else if(l==8) then
        if(m==-8) then
        funkcjakulista=3.d0/256.d0*sqrt(12155.d0/(2.d0*pi))*
     &  exp(-8.d0*ci*phi)*sinth**8
        else if(m==-7) then
        funkcjakulista=3.d0/64.d0*sqrt(12155.d0/(2.d0*pi))*
     &  exp(-7.d0*ci*phi)*sinth**7*costh
        else if(m==-6) then
        funkcjakulista=1.d0/128.d0*sqrt(7293.d0/pi)*exp(-6.d0*ci*phi)
     &  *sinth**6*(15.d0*costh**2 - 1.d0)
        else if(m==-5) then
        funkcjakulista=3.d0/64.d0*sqrt(17017.d0/(2.d0*pi))*
     &  exp(-5.d0*ci*phi)*sinth**5*costh* (5.d0*costh**2 - 1.d0)
        else if(m==-4) then
        funkcjakulista=3.d0/128.d0*sqrt(1309.d0/(2.d0*pi))
     &  *exp(-4.d0*ci*phi)*sinth**4*(65.d0*costh**4-26.d0*costh**2
     &  +1.d0)
        else if(m==-3) then
        funkcjakulista=1.d0/64.d0*sqrt(19635.d0/(2.d0*pi))*
     &  exp(-3.d0*ci*phi)*sinth**3*costh*(39.d0*costh**4-26.d0*costh**2
     &  +3.d0)
        else if(m==-2) then
        funkcjakulista=3.d0/128.d0*sqrt(595.d0/pi)*exp(-2.d0*ci*phi)*
     &  sinth**2*(143.d0*costh**6-143.d0*costh**4+33.d0*costh**2-1.d0)
        else if(m==-1) then
        funkcjakulista=3.d0/64.d0*sqrt(17.d0/(2.d0*pi))*exp(-ci*phi)*
     &  sinth*costh*(715.d0*costh**6-1001.d0*costh**4+385.d0*costh**2
     &  -35.d0)
        else if(m==0) then
        funkcjakulista=1.d0/256.d0*sqrt(17.d0/pi)*(6435.d0*costh**8-
     &  12012.d0*costh**6+6930.d0*costh**4-1260.d0*costh**2+35.d0)
        else if(m==1) then
        funkcjakulista=-3.d0/64.d0*sqrt(17.d0/(2.d0*pi))*exp(ci*phi)
     &  *sinth*costh*(715.d0*costh**6-1001.d0*costh**4+385.d0*costh**2
     &  -35.d0)
        else if(m==2) then
        funkcjakulista=3.d0/128.d0*sqrt(595.d0/pi)*exp(2.d0*ci*phi)*
     &  sinth**2*(143.d0*costh**6-143.d0*costh**4+33.d0*costh**2-1.d0)
        else if(m==3) then
        funkcjakulista=-1.d0/64.d0*sqrt(19635.d0/(2.d0*pi))*
     &  exp(3.d0*ci*phi)*sinth**3*costh*(39.d0*costh**4-26.d0*costh**2
     &  +3.d0)
        else if(m==4) then
        funkcjakulista=3.d0/128.d0*sqrt(1309.d0/(2.d0*pi))*
     &  exp(4.d0*ci*phi)*sinth**4*(65.d0*costh**4-26.d0*costh**2+1.d0)
        else if(m==5) then
        funkcjakulista=-3.d0/64.d0*sqrt(17017.d0/(2.d0*pi))*
     &  exp(5.d0*ci*phi)*sinth**5*costh* (5.d0*costh**2 - 1.d0)
        else if(m==6) then
        funkcjakulista=1.d0/128.d0*sqrt(7293.d0/pi)*exp(6.d0*ci*phi)*
     &  sinth**6*(15.d0*costh**2 - 1.d0)
        else if(m==7) then
        funkcjakulista=-3.d0/64.d0*sqrt(12155.d0/(2.d0*pi))*
     &  exp(7.d0*ci*phi)*sinth**7*costh
        else if(m==8) then
        funkcjakulista=3.d0/256.d0*sqrt(12155.d0/(2.d0*pi))*
     &  exp(8.d0*ci*phi)*sinth**8
        endif
      else if(l==9) then
        if(m==-9) then
        funkcjakulista=1.d0/512.d0*sqrt(230945.d0/pi)*exp(-9.d0*ci*phi)
     &  *sinth**9
        else if(m==-8) then
        funkcjakulista=3.d0/256.d0*sqrt(230945.d0/(2.d0*pi))*
     &  exp(-8.d0*ci*phi)*sinth**8*costh
        else if(m==-7) then
        funkcjakulista=3.d0/512.d0*sqrt(13585.d0/pi)*exp(-7.d0*ci*phi)*
     &  sinth**7*(17.d0*costh**2 - 1.d0)
        else if(m==-6) then
        funkcjakulista=1.d0/128.d0*sqrt(40755.d0/pi)*exp(-6.d0*ci*phi)
     &  *sinth**6*costh*(17.d0*costh**2 - 3.d0)
        else if(m==-5) then
        funkcjakulista=3.d0/256.d0*sqrt(2717.d0/pi)*exp(-5.d0*ci*phi)*
     &  sinth**5*(85.d0*costh**4-30.d0*costh**2+1.d0)
        else if(m==-4) then
        funkcjakulista=3.d0/128.d0*sqrt(95095.d0/(2.d0*pi))*
     &  exp(-4.d0*ci*phi)*sinth**4*costh*(17.d0*costh**4-10.d0*costh**2
     &  + 1.d0)
        else if(m==-3) then
        funkcjakulista=1.d0/256.d0*sqrt(21945.d0/pi)*exp(-3.d0*ci*phi)
     &  *sinth**3*(221.d0*costh**6-195.d0*costh**4+39.d0*costh**2 -1.d0)
        else if(m==-2) then
        funkcjakulista=3.d0/128.d0*sqrt(1045.d0/pi)*exp(-2.d0*ci*phi)*
     &  sinth**2*costh*(221.d0*costh**6-273.d0*costh**4+91.d0*costh**2-
     &  7.d0)
        else if(m==-1) then
        funkcjakulista=3.d0/256.d0*sqrt(95.d0/(2.d0*pi))*exp(-ci*phi)*
     &  sinth*(2431.d0*costh**8-4004.d0*costh**6+2002.d0*costh**4-
     &  308.d0*costh**2+7.d0)
        else if(m==0) then
        funkcjakulista=1.d0/256.d0*sqrt(19.d0/pi)*(12155.d0*costh**9-
     &  25740.d0*costh**7+18018.d0*costh**5-4620.d0*costh**3+315.d0
     &  *costh)
        else if(m==1) then
        funkcjakulista=-3.d0/256.d0*sqrt(95.d0/(2.d0*pi))*exp(ci*phi)*
     &  sinth*(2431.d0*costh**8-4004.d0*costh**6+2002.d0*costh**4-
     &  308.d0*costh**2+7.d0)
        else if(m==2) then
        funkcjakulista=3.d0/128.d0*sqrt(1045.d0/pi)*exp(2.d0*ci*phi)*
     &  sinth**2*costh*(221.d0*costh**6-273.d0*costh**4+91.d0*costh**2-
     &  7.d0)
        else if(m==3) then
        funkcjakulista=-1.d0/256.d0*sqrt(21945.d0/pi)*exp(3.d0*ci*phi)
     &  *sinth**3*(221.d0*costh**6-195.d0*costh**4+39.d0*costh**2-1.d0)
        else if(m==4) then
        funkcjakulista=3.d0/128.d0*sqrt(95095.d0/(2.d0*pi))*
     &  exp(4.d0*ci*phi)*sinth**4*costh*(17.d0*costh**4-10.d0*costh**2
     &  +1.d0)
        else if(m==5) then
        funkcjakulista=-3.d0/256.d0*sqrt(2717.d0/pi)*exp(5.d0*ci*phi)*
     &  sinth**5*(85.d0*costh**4-30.d0*costh**2+1.d0)
        else if(m==6) then
        funkcjakulista=1.d0/128.d0*sqrt(40755.d0/pi)*exp(6.d0*ci*phi)*
     &  sinth**6*costh*(17.d0*costh**2 - 3.d0)
        else if(m==7) then
        funkcjakulista=-3.d0/512.d0*sqrt(13585.d0/pi)*exp(7.d0*ci*phi)*
     &  sinth**7*(17.d0*costh**2 - 1.d0)
        else if(m==8) then
        funkcjakulista=3.d0/256.d0*sqrt(230945.d0/(2.d0*pi))*
     &  exp(8.d0*ci*phi)*sinth**8*costh
        else if(m==9) then
        funkcjakulista=-1.d0/512.d0*sqrt(230945.d0/pi)*
     &  exp(9.d0*ci*phi)*sinth**9
        endif
      else if(l==10) then
        if (m==-10) then
        funkcjakulista=(sqrt(969969.d0/pi)*
     &  exp(-10.d0*ci*phi)*sinth**10)/1024.d0
        else if(m==-9) then
        funkcjakulista=1.d0/512.d0*sqrt(4849845.d0/pi)*
     &  exp(-9.d0*ci*phi)*sinth**9*costh
        else if(m==-8) then
        funkcjakulista=1.d0/512.d0*sqrt(255255.d0/(2.d0*pi))*
     &  exp(-8.d0*ci*phi)*sinth**8*(19.d0*costh**2-1.d0)
        else if(m==-7) then
        funkcjakulista=3.d0/512.d0*sqrt(85085.d0/pi)*exp(-7.d0*ci*phi)
     &  *sinth**7*costh*(19.d0*costh**2-3.d0)
        else if(m==-6) then
        funkcjakulista=(3.d0*sqrt(5005.d0/pi)*exp(-6.d0*ci*phi)*
     &  sinth**6*(323.d0*costh**4-102.d0*costh**2+3.d0))/1024.d0
        else if(m==-5) then
        funkcjakulista=3.d0/256.d0*sqrt(1001.d0/pi)*exp(-5.d0*ci*phi)*
     &  sinth**5*costh*(323.d0*costh**4-170.d0*costh**2+15.d0)
        else if(m==-4) then
        funkcjakulista=3.d0/256.d0*sqrt(5005.d0/(2.d0*pi))*
     &  exp(-4.d0*ci*phi)*sinth**4*(323.d0*costh**6-255.d0*costh**4+
     &  45.d0*costh**2-1.d0)
        else if(m==-3) then
        funkcjakulista=3.d0/256.d0*sqrt(5005.d0/pi)*exp(-3.d0*ci*phi)*
     &  sinth**3*costh*(323.d0*costh**6-357.d0*costh**4+105.d0*costh**2
     &  -7.d0)
        else if(m==-2) then
        funkcjakulista=3.d0/512.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(-2.d0*ci*phi)*sinth**2*(4199.d0*costh**8-6188.d0*costh**6+
     &  2730.d0*costh**4-364.d0*costh**2+7.d0)
        else if(m==-1) then
        funkcjakulista=1.d0/256.d0*sqrt(1155.d0/(2.d0*pi))*exp(-ci*phi)
     &  *sinth*costh*(4199.d0*costh**8-7956.d0*costh**6+
     &  4914.d0*costh**4-1092.d0*costh**2+63.d0)
        else if(m==0) then
        funkcjakulista=1.d0/512.d0*sqrt(21.d0/pi)*(46189.d0*costh**10-
     &  109395.d0*costh**8+90090.d0*costh**6-30030.d0*costh**4+
     &  3465.d0*costh**2-63.d0)
        else if(m==1) then
        funkcjakulista=-1.d0/256.d0*sqrt(1155.d0/(2.d0*pi))*exp(ci*phi)
     &  *sinth*costh*(4199.d0*costh**8-7956.d0*costh**6+
     &  4914.d0*costh**4-1092.d0*costh**2+63.d0)
        else if(m==2) then
        funkcjakulista=3.d0/512.d0*sqrt(385.d0/(2.d0*pi))*
     &  exp(2.d0*ci*phi)*sinth**2*(4199.d0*costh**8-6188.d0*costh**6+
     &  2730.d0*costh**4-364.d0*costh**2+7.d0)
        else if(m==3) then
        funkcjakulista=-3.d0/256.d0*sqrt(5005.d0/pi)*exp(3.d0*ci*phi)*
     &  sinth**3*costh*(323.d0*costh**6-357.d0*costh**4+105.d0*costh**2
     &  -7.d0)
        else if(m==4) then
        funkcjakulista=3.d0/256.d0*sqrt(5005.d0/(2.d0*pi))*
     &  exp(4.d0*ci*phi)*sinth**4*(323.d0*costh**6-255.d0*costh**4+
     &  45.d0*costh**2-1.d0)
        else if(m==5) then
        funkcjakulista=-3.d0/256.d0*sqrt(1001.d0/pi)*exp(5.d0*ci*phi)*
     &  sinth**5*costh*(323.d0*costh**4-170.d0*costh**2+15.d0)
        else if(m==6) then
        funkcjakulista=(3.d0*sqrt(5005.d0/pi)*exp(6.d0*ci*phi)*
     &  sinth**6*(323.d0*costh**4-102.d0*costh**2+3.d0))/1024.d0
        else if(m==7) then
        funkcjakulista=-3.d0/512.d0*sqrt(85085.d0/pi)*exp(7.d0*ci*phi)*
     &  sinth**7*costh*(19.d0*costh**2-3.d0)
        else if(m==8) then
        funkcjakulista=1.d0/512.d0*sqrt(255255.d0/(2.d0*pi))*
     &  exp(8.d0*ci*phi)*sinth**8*(19.d0*costh**2-1.d0)
        else if(m==9) then
        funkcjakulista=-1.d0/512.d0*sqrt(4849845.d0/pi)*
     &  exp(9.d0*ci*phi)*sinth**9*costh
        else if(m==10) then
        funkcjakulista=(sqrt(969969.d0/pi)*exp(10.d0*ci*phi)*sinth**10)
     &  /1024.d0
        endif
      else
      print*,"l-value beyond the scope..."
      endif
      end 
c------------------------------------------------------------------------
! pochodne wielomianów Legendre'a
      double precision function LegendrePrim(n,x)
      implicit double precision (a-h,o-z)
      if (n.eq.0) then
      LegendrePrim=0.d0
      else if (n.eq.1) then
      LegendrePrim=1.d0
      else if (n.eq.2) then
      LegendrePrim=3.d0*x
      else if (n.eq.3) then
      LegendrePrim=0.5d0*(15.d0*x**2-3.d0)
      else if (n.eq.4) then
      LegendrePrim=1.d0/8.d0*(140.d0*x**3-60.d0*x)
      else if (n.eq.5) then
      LegendrePrim=(15.d0-210.d0*x**2+315.d0*x**4)/8.d0
      else if (n.eq.6) then
      LegendrePrim=(210.d0*x - 1260.d0*x**3 + 1386.d0*x**5)/16.d0
      else if (n.eq.7) then
      LegendrePrim=(-35.d0 + 945.d0*x**2 - 3465.d0*x**4 + 3003.d0*x**6)
     &/16.d0
      endif
      end
c------------------------------------------------------------------------
! funkcje Legendre'a drugiego rodzaju
      double precision function LegendreQ(n,x)
      integer n
      double precision x
      if (n.eq.0) then
      LegendreQ=-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x)
      else if (n.eq.1) then
      LegendreQ=-1.d0+x*(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      else if (n.eq.2) then
      LegendreQ=-3.d0*x/2.d0+0.5d0*(-1.d0+3.d0*x**2)*
     &(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      else if (n.eq.3) then
      LegendreQ=2.d0/3.d0-5.d0/2.d0*x**2-0.5d0*x*(3.d0-5.d0*x**2)*
     &(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      else if (n.eq.4) then
      LegendreQ=55.d0/24.d0*x-35.d0/8.d0*x**3+1.d0/8.d0*
     &(3.d0-30.d0*x**2+35.d0*x**4)*
     &(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      else if (n.eq.5) then
      LegendreQ=-8.d0/15.d0+49.d0/8.d0*x**2-63.d0/8.d0*x**4+
     &1.d0/8.d0*x*(15.d0-70.d0*x**2+63.d0*x**4)*
     &(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      else if (n.eq.6) then
      LegendreQ=-231.d0/80.d0*x+119.d0/8.d0*x**3-231.d0/16.d0*x**5+
     &1.d0/16.d0*(-5.d0+105.d0*x**2-315.d0*x**4+231.d0*x**6)*
     &(-0.5d0*dlog(x-1.d0)+0.5d0*dlog(1.d0+x))
      endif
      end

      SUBROUTINE CHGM(A,B,X,HG)
C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,... )
C                x  --- Argument
C       Output:  HG --- M(a,b,x)
C       Routine called: GAMMA for computing â(x)
C       ===================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
        IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
           HG=1.0D+300
        ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
           HG=1.0D0
        ELSE IF (A.EQ.-1.0D0) THEN
           HG=1.0D0-X/B
        ELSE IF (A.EQ.B) THEN
           HG=DEXP(X)
        ELSE IF (A-B.EQ.1.0D0) THEN
           HG=(1.0D0+X/B)*DEXP(X)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           HG=(DEXP(X)-1.0D0)/X
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           R=1.0D0
           HG=1.0D0
           DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
        ENDIF
        IF (HG.NE.0.0D0) RETURN
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        IF (A.LT.2.0D0) NL=0
        IF (A.GE.2.0D0) THEN
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
           ELSE
              CALL GAMMA(A,TA)
              CALL GAMMA(B,TB)
              XG=B-A
              CALL GAMMA(XG,TBA)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 SUM1=SUM1+R1
20               SUM2=SUM2+R2
              HG1=TB/TBA*X**(-A)*DCOS(PI*A)*SUM1
              HG2=TB/TA*DEXP(X)*X**(A-B)*SUM2
              HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
35            A=A+1.0D0
        ENDIF
        IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
        A=A1
        X=X0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
c---------------------------------------------------
