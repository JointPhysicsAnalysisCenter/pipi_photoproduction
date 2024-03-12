      module parameters
      implicit none
      real (8), parameter :: mp     = 0.93827203d0
      real (8), parameter :: mpi    = 0.13957018d0
      real (8), parameter :: pi     = 4.d0*datan(1.d0)
!      real (8), parameter :: eps    = 0.d0
      real (8), parameter :: hbc    = 0.19732d0
!     parameters of the Regge parametrization
      real (8), parameter, dimension (0:19) :: param  = (/
     c     0.4895696553113608d0, ! alpha rho
     c     0.9425449473756354d0,
     c     0.4895696553113608d0, ! alpha f
     c     0.9425449473756354d0,
     c     1.07473280868023d0,   ! alpha Pomeron
     c     0.4343313339478808d0,
     c     0.16164225485369113d0,
     c     5.00657330231329d0,   ! C(0,1,2) rho
     c     10.0982220619908d0,
     c     0.5730288256207643d0,
     c     128.86955863248116d0, ! D(0,1) rho
     c     1.3840421677412782d0,
     c     71.35385175012772d0,  ! C(0,1) f
     c     3.1790121902525574d0,
     c     71.35385175012772d0,  ! D(0,1) f
     c     3.1790121902525574d0,
     c     23.894196195009776d0, ! C(0,1) Pomeron
     c     2.205005669899481d0,
     c     23.894196195009776d0, ! D(0,1) Pomeron
     c     2.205005669899481d0 /)
      complex (8), parameter :: xr = dcmplx(1.d0,0.d0)
      complex (8), parameter :: xi = dcmplx(0.d0,1.d0)
      complex (8), parameter :: xzero = dcmplx(0.d0,0.d0)
      end module parameters

!      include 'cgamma.f90'

!       program main
!       use parameters
!       implicit none
!       integer :: i
!       integer, parameter        :: unitinput = 10
!       integer, parameter        :: nSAID = 30
!       character (25), parameter :: inputfilename = 'SAID_PW/list.txt'
!       character (25), dimension (nSAID) :: SAIDfile
!       real (8) :: s, plab, Elab, t
!       real (8), dimension (0:1) :: sigtot
!       complex (8), dimension (0:1,0:2) :: ABC
!       complex (8), dimension (0:4,0:9,0:249) :: PW
! !     Read SAID file names
! !      open(unit=unitinput, file=inputfilename, status='unknown')
! !      do i = 1, nSAID
! !         read(unitinput,*) SAIDfile(i)
! !         write (0,*) i, SAIDfile(i)
! !      enddo
! 
!       call SAID_PiN_PW(PW)
!       
!       open(unit=1,file='sigtot.dat',status='unknown')
!       
!       s = 1.2d0
!       t = 0.0d0
!       do i=1, 6000
!          call sigtotal (s,PW,sigtot)
!          !call Full_ABC (s,t,PW,ABC)
!          Elab = (s - mp*mp - mpi*mpi)/(2.d0*mp)
!          plab = sqrt(Elab*Elab - mpi*mpi)
! !          write (*,*) s, plab, sigtot(1), sigtot(0)
!          write (1,*) s, plab, sigtot(1), sigtot(0)
!          s = s + 5.d-3
!          print*,plab
!       enddo
!       end program main
      

      module piNamps
      contains

      subroutine observables (s,t,tpi,PW,obs)
      use parameters
      implicit none
      real (8) :: s, t, tpi
      real (8), dimension (0:2,0:1) :: obs
      complex (8), dimension (0:4,0:9,0:249) :: PW
      real (8) ::  plab, Elab, qcm, z, sint, fac, tau, sth, spth
      complex (8), dimension (0:1,0:2) :: ABC
      complex (8) :: B, C
      sth  = (mp + mpi)**2
      spth = (mp - mpi)**2
      Elab = (s - mp*mp - mpi*mpi)/(2.d0*mp)
      plab = dsqrt(Elab*Elab - mpi*mpi)
      qcm  = dsqrt((s-sth)*(s-spth))/(2.d0*dsqrt(s))
      z    = 1.d0 + t/(2.d0*qcm*qcm)
      sint = dsqrt(1.d0 - z*z)
      fac  = 1.d0/pi/s * (mp/4.d0/qcm)**2
      tau  = t/4.d0/mp/mp

      call Full_ABC(s,t,tpi,PW,ABC)

!     pi^- p --> pi^- p
      C = ABC(0,2) + ABC(1,2)
      B = ABC(0,1) + ABC(1,1)
      
      obs(0,0) = fac*( (1.d0-tau)*abs(C)*abs(C) - 
     &     tau*( (s*tau + plab*plab)/(1.d0 - tau) )*abs(B)*abs(B) )
      obs(0,1) = -sint/(16.d0*pi*dsqrt(s))*aimag(C*dconjg(B))/obs(0,0)
      obs(0,0) = obs(0,0)*389.352d0

!     pi^+ p --> pi^+ p      
      C = ABC(0,2) - ABC(1,2)
      B = ABC(0,1) - ABC(1,1)

      obs(1,0) = fac*( (1.d0-tau)*abs(C)*abs(C) -
     &     tau*( (s*tau + plab*plab)/(1.d0 - tau) )*abs(B)*abs(B) )
      obs(1,1) = -sint/(16.d0*pi*dsqrt(s))*aimag(C*dconjg(B))/obs(1,0)
      obs(1,0) = obs(1,0)*389.352d0

!     pi^- p --> pi^- p !!!!!!!!which reaction is this?
      C = -dsqrt(2.d0)*ABC(1,2)
      B = -dsqrt(2.d0)*ABC(1,1)
      
      obs(2,0) = fac*( (1.d0-tau)*abs(C)*abs(C) - 
     &     tau*( (s*tau + plab*plab)/(1.d0-tau) )*abs(B)*abs(B) )
      obs(2,1) = -sint/(16.d0*pi*dsqrt(s))*aimag(C*dconjg(B))/obs(2,0)
      obs(2,0) = obs(2,0)*389.352d0
      return
      end subroutine observables


      subroutine sigtotal (s,PW,sigtot)
      use parameters
      implicit none
      real (8) :: s,tpi
      real (8), dimension (0:1) :: sigtot
      complex (8), dimension (0:4,0:9,0:249) :: PW
      real (8), parameter :: t = 0.d0
      real (8) :: plab, Elab
      complex (8), dimension (0:1,0:2) :: ABC
      ABC = xzero
      Elab = (s - mp*mp - mpi*mpi)/(2.d0*mp)
      plab = dsqrt(Elab*Elab - mpi*mpi)
      call Full_ABC(s,t,mpi**2,PW,ABC)
!       print*,s,ABC,plab
      sigtot(0) = 0.389352d0*aimag(ABC(0,2) + ABC(1,2))/plab
      sigtot(1) = 0.389352d0*aimag(ABC(0,2) - ABC(1,2))/plab
      return
      end subroutine sigtotal


      subroutine Full_ABC(s,t,tpi,PW,ABC)
      use parameters
      implicit none
      real (8), intent(in):: s, t, tpi
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), intent(out), dimension (0:1,0:2) :: ABC
      complex (8), dimension (0:1,0:2) :: ABCSaid, ABCRegge
      real (8) :: Elab, eps
      real (8), parameter :: thres1 = 1.5d0, thres2 = 2.1d0
      Elab = ( s - mp*mp - mpi*mpi)/(2.d0*mp)
! modification to account for pion virtuality
c        Elab = ( s - mp*mp - tpi)/(2.d0*mp)
! modification to account for pion virtuality
      ABC = xzero
      ABCSaid = xzero
      ABCRegge = xzero
      if (Elab .lt. thres1) then
         call SAID_ABC(s,t,tpi,PW,ABC) !version w-o threshold factor
c           call SAID_ABC_thcor(s,t,tpi,PW,ABC) !version with threshold factor
      elseif (Elab .gt. thres2) then
         call Regge_ABC(s,t,tpi,ABC)
      else! ( (Elab .ge. thres1) .and. (Elab .le. thres2) then
c           call SAID_ABC_thcor(s,t,tpi,PW,ABCSaid)
         call SAID_ABC(s,t,tpi,PW,ABCSaid)
         call Regge_ABC(s,t,tpi,ABCRegge)
         eps = (Elab - thres1) / (thres2 - thres1)
         ABC(0,0) = (1.d0-eps)*ABCSaid(0,0) + eps*ABCRegge(0,0)
         ABC(0,1) = (1.d0-eps)*ABCSaid(0,1) + eps*ABCRegge(0,1)
         ABC(0,2) = (1.d0-eps)*ABCSaid(0,2) + eps*ABCRegge(0,2)

         ABC(1,0) = (1.d0-eps)*ABCSaid(1,0) + eps*ABCRegge(1,0)
         ABC(1,1) = (1.d0-eps)*ABCSaid(1,1) + eps*ABCRegge(1,1)
         ABC(1,2) = (1.d0-eps)*ABCSaid(1,2) + eps*ABCRegge(1,2)
      endif
      return
      end subroutine Full_ABC


      subroutine SAID_ABC (s,t,tpi,PW,ABC,lfix)
      use parameters
      use, intrinsic :: omp_lib
      implicit none
      real (8)                          :: s, t, tpi,xLbd
      complex (8), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2)  :: ABC
      integer, optional:: lfix
      integer :: L, k, i
      real (8) :: w, E, qcm, z, Tlab, Elab, plab,E1,E2,p1,p2,Zpl1,Zpl2,
     &Zmi1,Zmi2
      real (8) :: eps, sth, spth
      real (8), dimension (0:9) :: pol, pol1p
      complex (8) :: fp,fm,gp,gm,f1p,f1m,f2p,f2m,f1,g1,f3,g3
      complex (8), dimension (0:9) :: PW1m, PW1p, PW3m, PW3p

      ABC = xzero

      w = dsqrt(s)
!       E = (s + mp*mp - mpi*mpi)/(2.d0*w)
! modification to account for pion virtuality
 !     E = (s + mp*mp - tpi)/(2.d0*w)
! modification to account for pion virtuality
      sth = (mp + mpi)**2
      spth = (mp - mpi)**2
      
      E1=(s+mp**2-tpi)/(2.d0*w)
      E2=(s+mp**2-mpi**2)/(2.d0*w)
      p1=dsqrt(E1**2-mp**2)
      p2=dsqrt(E2**2-mp**2)
      
      Zpl1=dsqrt(E1+mp)
      Zmi1=dsqrt(E1-mp)
      Zpl2=dsqrt(E2+mp)
      Zmi2=dsqrt(E2-mp)
      
      
!       qcm = dsqrt((s-sth)*(s-spth))/(2.d0*w)
      qcm=dsqrt(p1*p2)
!       z = 1.d0 + t/(2.d0*qcm*qcm) ! z = cos(theta)
! modification to account for pion virtuality
      z=(2.d0*s*(t-2.d0*mp**2)+(s+mp**2-tpi)*(s+mp**2-mpi**2))/
     &dsqrt(xlbd(s,mp**2,tpi)*xlbd(s,mp**2,mpi**2))
! modification to account for pion virtuality
       Tlab = (s - mp*mp - mpi*mpi)/(2.d0*mp)-mpi !KE of pion in lab
! modification to account for pion virtuality
!      Tlab = (s - mp*mp - tpi)/(2.d0*mp)-mpi !KE of pion in lab
! modification to account for pion virtuality
      Elab = Tlab + mpi
      
      plab = dsqrt(Elab*Elab - mpi*mpi)
! modification to account for pion virtuality
 !     plab = dsqrt(Elab*Elab - tpi)
! modification to account for pion virtuality
      if (plab.gt.2.5d0) return
            
      PW1m = xzero
      PW1p = xzero
      PW3m = xzero
      PW3p = xzero
      
!       do i=0,100
!       print*,i
!       print*,PW(0,1,i)
!       print*,PW(1,1,i)
!       print*,PW(2,1,i)
!       print*,PW(3,1,i)
!       print*,PW(4,1,i)
!       enddo
!       stop

      k = floor(plab/0.025d0)
c       print*,s,plab,omp_get_thread_num()
      eps = (plab - PW(0,0,k)) / 0.025d0

      do L = 0, 9
         PW1m (L) = (1.d0 - eps)*PW(1,L,k) + eps*PW(1,L,k+1)
         PW1p (L) = (1.d0 - eps)*PW(2,L,k) + eps*PW(2,L,k+1)
         PW3m (L) = (1.d0 - eps)*PW(3,L,k) + eps*PW(3,L,k+1)
         PW3p (L) = (1.d0 - eps)*PW(4,L,k) + eps*PW(4,L,k+1)
      enddo
      
      fp  = xzero
      fm  = xzero
      gp  = xzero
      gm  = xzero
      f1p = xzero
      f1m = xzero
      f2p = xzero
      f2m = xzero
      f1  = xzero
      g1  = xzero
      f3  = xzero
      g3  = xzero

      pol = 0.d0
      pol1p = 0.d0

      call legendreP(z,pol,pol1p)

      if(.not. present(lfix)) then
        do L = 0, 7
           f1 = f1 + ( (L+1.d0)*PW1p(L) + L*PW1m(L) )*pol(L)
           f3 = f3 + ( (L+1.d0)*PW3p(L) + L*PW3m(L) )*pol(L)
           g1 = g1 + ( PW1p(L) - PW1m(L) )*pol1p(L)
           g3 = g3 + ( PW3p(L) - PW3m(L) )*pol1p(L)
        enddo
      else
!         print*,"lfix set to ",lfix
        do L = 0, 7
           if(L/=lfix) cycle
           f1 = f1 + ( (L+1.d0)*PW1p(L) + L*PW1m(L) )*pol(L)
           f3 = f3 + ( (L+1.d0)*PW3p(L) + L*PW3m(L) )*pol(L)
           g1 = g1 + ( PW1p(L) - PW1m(L) )*pol1p(L)
           g3 = g3 + ( PW3p(L) - PW3m(L) )*pol1p(L)
        enddo    
      endif
      f1 = f1/qcm
      f3 = f3/qcm
      g1 = g1/qcm
      g3 = g3/qcm

      fp = (f1 + 2.d0*f3)/3.d0
      gp = (g1 + 2.d0*g3)/3.d0
      fm = (f1 - f3)/3.d0
      gm = (g1 - g3)/3.d0

      f2p = -gp
      f1p = fp - z*f2p
      f2m = -gm
      f1m = fm - z*f2m

      ABC(0,0) =4.d0*pi*((w+mp)/(Zpl1*Zpl2)*f1p-(w-mp)/(Zmi1*Zmi2)*f2p)
      ABC(1,0) =4.d0*pi*((w+mp)/(Zpl1*Zpl2)*f1m-(w-mp)/(Zmi1*Zmi2)*f2m)
      ABC(0,1) =4.d0*pi*( 1.d0/(Zpl1*Zpl2)*f1p + 1.d0/(Zmi1*Zmi2)*f2p )
      ABC(1,1) =4.d0*pi*( 1.d0/(Zpl1*Zpl2)*f1m + 1.d0/(Zmi1*Zmi2)*f2m )
      ABC(0,2) = ABC(0,0) + 
     &     (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*ABC(0,1)
      ABC(1,2) = ABC(1,0) +
     &     (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*ABC(1,1)

      return
      end subroutine SAID_ABC

      subroutine SAID_ABC_thcor (s,t,tpi,PW,ABC,lfix)
      use parameters
      implicit none
      real (8)                          :: s, t, tpi,xLbd
      complex (8), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2)  :: ABC
      integer, optional:: lfix
      integer :: L, k, i
      real (8) :: w, E, qcm, z, Tlab, Elab, plab,E1,E2,p1,p2,Zpl1,Zpl2,
     &Zmi1,Zmi2
      real (8) :: eps, sth, spth
      real (8), dimension (0:9) :: pol, pol1p
      complex (8) :: fp,fm,gp,gm,f1p,f1m,f2p,f2m,f1,g1,f3,g3
      complex (8), dimension (0:9) :: PW1m, PW1p, PW3m, PW3p
      
      real:: poff,pon,ratio, lambda=1.d0

      ABC = xzero

      w = dsqrt(s)
!       E = (s + mp*mp - mpi*mpi)/(2.d0*w)
! modification to account for pion virtuality
 !     E = (s + mp*mp - tpi)/(2.d0*w)
! modification to account for pion virtuality
      sth = (mp + mpi)**2
      spth = (mp - mpi)**2
      
      E1=(s+mp**2-tpi)/(2.d0*w)
      E2=(s+mp**2-mpi**2)/(2.d0*w)
      p1=dsqrt(E1**2-mp**2)
      p2=dsqrt(E2**2-mp**2)
      
      Zpl1=dsqrt(E1+mp)
      Zmi1=dsqrt(E1-mp)
      Zpl2=dsqrt(E2+mp)
      Zmi2=dsqrt(E2-mp)
      
      
!       qcm = dsqrt((s-sth)*(s-spth))/(2.d0*w)
      qcm=dsqrt(p1*p2)
!       z = 1.d0 + t/(2.d0*qcm*qcm) ! z = cos(theta)
! modification to account for pion virtuality
      z=(2.d0*s*(t-2.d0*mp**2)+(s+mp**2-tpi)*(s+mp**2-mpi**2))/
     &dsqrt(xlbd(s,mp**2,tpi)*xlbd(s,mp**2,mpi**2))
! modification to account for pion virtuality
       Tlab = (s - mp*mp - mpi*mpi)/(2.d0*mp)-mpi !KE of pion in lab
! modification to account for pion virtuality
!      Tlab = (s - mp*mp - tpi)/(2.d0*mp)-mpi !KE of pion in lab
! modification to account for pion virtuality
      Elab = Tlab + mpi
      
       plab = dsqrt(Elab*Elab - mpi*mpi)
c        print*,"pl=",plab,Elab,Tlab,s
! modification to account for pion virtuality
 !     plab = dsqrt(Elab*Elab - tpi)
! modification to account for pion virtuality
      if (plab.gt.2.5d0) return
            
      PW1m = xzero
      PW1p = xzero
      PW3m = xzero
      PW3p = xzero
      

      k = floor(plab/0.025d0)
      eps = (plab - PW(0,0,k)) / 0.025d0

      do L = 0, 9
         PW1m (L) = (1.d0 - eps)*PW(1,L,k) + eps*PW(1,L,k+1)
         PW1p (L) = (1.d0 - eps)*PW(2,L,k) + eps*PW(2,L,k+1)
         PW3m (L) = (1.d0 - eps)*PW(3,L,k) + eps*PW(3,L,k+1)
         PW3p (L) = (1.d0 - eps)*PW(4,L,k) + eps*PW(4,L,k+1)
      enddo
      
      fp  = xzero
      fm  = xzero
      gp  = xzero
      gm  = xzero
      f1p = xzero
      f1m = xzero
      f2p = xzero
      f2m = xzero
      f1  = xzero
      g1  = xzero
      f3  = xzero
      g3  = xzero

      pol = 0.d0
      pol1p = 0.d0

      call legendreP(z,pol,pol1p)
      
      poff=p1
      pon=p2
      
!       ratio=poff/pon*(lambda+pon)/(lambda+poff)
      ratio=poff/pon
      
      if(.not. present(lfix)) then
        do L = 0, 7
          f1 = f1 + ( (L+1.d0)*PW1p(L) + L*PW1m(L) )*pol(L)*
     &    ratio**(L+0.5)
          f3 = f3 + ( (L+1.d0)*PW3p(L) + L*PW3m(L) )*pol(L)*
     &    ratio**(L+0.5)
          g1 = g1 + ( PW1p(L) - PW1m(L) )*pol1p(L)*
     &    ratio**(L+0.5)
          g3 = g3 + ( PW3p(L) - PW3m(L) )*pol1p(L)*
     &    ratio**(L+0.5)
        enddo
      else
!         print*,"lfix set to ",lfix
        do L = 0, 7
           if(L/=lfix) cycle
           f1 = f1 + ( (L+1.d0)*PW1p(L) + L*PW1m(L) )*pol(L)*
     &     ratio**(L+0.5)
           f3 = f3 + ( (L+1.d0)*PW3p(L) + L*PW3m(L) )*pol(L)*
     &     ratio**(L+0.5)
           g1 = g1 + ( PW1p(L) - PW1m(L) )*pol1p(L)*
     &     ratio**(L+0.5)
           g3 = g3 + ( PW3p(L) - PW3m(L) )*pol1p(L)*
     &     ratio**(L+0.5)
        enddo    
      endif

!       do L = 0, 7
!          f1 = f1 + ( (L+1.d0)*PW1p(L) + L*PW1m(L) )*pol(L)*ratio**L
!          f3 = f3 + ( (L+1.d0)*PW3p(L) + L*PW3m(L) )*pol(L)*ratio**L
!          g1 = g1 + ( PW1p(L) - PW1m(L) )*pol1p(L)*ratio**(L-1)
!          g3 = g3 + ( PW3p(L) - PW3m(L) )*pol1p(L)*ratio**(L-1)
!       enddo
      f1 = f1/qcm
      f3 = f3/qcm
      g1 = g1/qcm
      g3 = g3/qcm

      fp = (f1 + 2.d0*f3)/3.d0
      gp = (g1 + 2.d0*g3)/3.d0
      fm = (f1 - f3)/3.d0
      gm = (g1 - g3)/3.d0

      f2p = -gp
      f1p = fp - z*f2p
      f2m = -gm
      f1m = fm - z*f2m

      ABC(0,0) =4.d0*pi*((w+mp)/(Zpl1*Zpl2)*f1p-(w-mp)/(Zmi1*Zmi2)*f2p)
      ABC(1,0) =4.d0*pi*((w+mp)/(Zpl1*Zpl2)*f1m-(w-mp)/(Zmi1*Zmi2)*f2m)
      ABC(0,1) =4.d0*pi*( 1.d0/(Zpl1*Zpl2)*f1p + 1.d0/(Zmi1*Zmi2)*f2p )
      ABC(1,1) =4.d0*pi*( 1.d0/(Zpl1*Zpl2)*f1m + 1.d0/(Zmi1*Zmi2)*f2m )
      ABC(0,2) = ABC(0,0) + 
     &     (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*ABC(0,1)
      ABC(1,2) = ABC(1,0) +
     &     (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*ABC(1,1)

      return
      end subroutine SAID_ABC_thcor
      
      subroutine Regge_ABC(s,t,tpi,ABC)
      use parameters
      implicit none
      real (8) :: s, t, tpi
      complex (8), dimension (0:1,0:2) :: ABC
      integer  :: mo = 0
      real (8) :: a0rho, aprho, a0f, apf, a0P, apP, appP
      real (8) :: C0rho, C1rho, C2rho, D0rho, D1rho, C0f, C1f
      real (8) :: D0f, D1f, C0P, C1P, D0P, D1P
      real (8) :: Elab, som, u, nu
      complex (8) :: alpRho, alpF, alpPom
      complex (8) :: R1rho, R2rho, R1f, R2f, R1Pom, R2Pom
      complex (8) :: Crho, Brho, Cf, Bf, CPom, BPom
      complex (8) :: Arho, Af, APom
      complex (8), external :: cgamma

      ABC = xzero

      a0rho = param(0)
      aprho = param(1)
      a0f   = param(2)
      apf   = param(3)
      a0P   = param(4)
      apP   = param(5)
      appP  = param(6)
      C0rho = param(7)
      C1rho = param(8)
      C2rho = param(9)
      D0rho = param(10)
      D1rho = param(11)
      C0f   = param(12)
      C1f   = param(13)
      D0f   = param(14) 
      D1f   = param(15)
      C0P   = param(16)
      C1P   = param(17)
      D0P   = param(18)
      D1P   = param(19)
! modification to account for pion virtuality
!       Elab = (s - mp*mp - mpi*mpi)/(2.d0*mp)
      Elab = (s - mp*mp - tpi)/(2.d0*mp)
!       som  = 2.d0*mp*mp + 2.d0*mpi*mpi
      som  = 2.d0*mp*mp + mpi*mpi+tpi
! modification to account for pion virtuality      
      som  = 2.d0*mp*mp + mpi*mpi+tpi
      u    = -s - t + som
      nu   = (s-u)/(4.d0*mp)

      alpRho = a0rho + t*aprho
      alpF   = a0f + t*apf
      alpPom = a0P + t*apP + t*t*appP



      R1rho = cgamma(mo,0.d0-alpRho)/2.d0*(1.d0 - 
     &     cdexp(-xi*pi*alpRho))*nu**alpRho
      R2rho = cgamma(mo,1.d0-alpRho)/2.d0*(1.d0 - 
     &     cdexp(-xi*pi*alpRho))*nu**(alpRho-1)
      R1f = cgamma(mo,1.d0-alpF)/2.d0*(1.d0 + 
     &     cdexp(-xi*pi*alpF))*nu**(alpF)
      R2f = cgamma(mo,1.d0-alpF)/2.d0*(1.d0 +
     &     cdexp(-xi*pi*alpF))*nu**(alpF-1)
      R1Pom = cgamma(mo,1.d0-alpPom)/2.d0*(1.d0 + 
     &     cdexp(-xi*pi*alpPom))*nu**alpPom
      R2Pom = cgamma(mo,1.d0-alpPom)/2.d0*(1.d0+
     &     cdexp(-xi*pi*alpPom))*nu**(alpPom-1)


      Crho = -C0rho*((1.d0 + C2rho)*dexp(C1rho*t) - C2rho)*R1rho
      Brho = D0rho*dexp(D1rho*t)*R2rho
      Cf   = -C0f*dexp(C1f*t)*R1f
      Bf   = -D0f*dexp(D1f*t)*R2f
      CPom = -C0P*dexp(C1P*t)*R1Pom
      BPom = -D0P*dexp(D1P*t)*R2Pom


      Arho = Crho - (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*Brho
      Af   = Cf   - (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*Bf 
      APom = CPom - (Elab + t/(4.d0*mp))/(1.d0 - t/(4.d0*mp*mp))*BPom

      ABC(0,0) = APom + Af
      ABC(0,1) = BPom + Bf
      ABC(0,2) = CPom + Cf
      
      ABC(1,0) = Arho
      ABC(1,1) = Brho
      ABC(1,2) = Crho

      return
      end subroutine Regge_ABC

      
      subroutine SAID_PiN_PW(PW)
      use parameters
      implicit none
      complex (8), intent(out), dimension (0:4,0:9,0:249) :: PW
      integer, parameter :: maxbin = 100
      real (8), dimension (0:100) :: plab
      real (8), dimension (0:4,0:9,0:100) :: RePW, ImPW
      real (8), dimension (0:5) :: tmp

      integer :: i, L
      integer, parameter        :: unitinput = 10
      integer, parameter        :: nSAID = 30
      character (25), parameter :: inputfilename = 'SAID_PW/list.txt'
      character (35), dimension (1:nSAID) :: SAIDfile
      integer, dimension (1:nSAID) :: SAIDnum

      PW = xzero

!     Read SAID file names
      open(unit=unitinput, file=inputfilename, status='old')
      do i = 1, nSAID
         read(unitinput,*) SAIDfile(i)
!         write (0,*) i,SAIDfile(i)
         SAIDnum(i) = 300+i+1
         SAIDfile(i) = 'SAID_PW/' // SAIDfile(i)
      enddo

!      do i=1, nSAID
!      write (0,*) i, SAIDnum(i), SAIDfile(i)
!      enddo
      
      plab = 0.d0
      RePW = 0.d0
      ImPW = 0.d0
      tmp  = 0.d0

      do i=1,nSAID
         open (SAIDnum(i), file=SAIDfile(i), status='old')
      enddo
      do i = 0, maxbin
         read (SAIDnum(1),*) plab(i), tmp(0), tmp(1), tmp(2), tmp(3), !loading S-wave amps -2 IJ combinations
     &        RePW(2,0,i), ImPW(2,0,i), tmp(4), tmp(5)
         read (SAIDnum(2),*) plab(i), tmp(0), tmp(1), tmp(2), tmp(3),
     &        RePW(4,0,i), ImPW(4,0,i),tmp(4),tmp(5)
         PW(2,0,i) = RePW(2,0,i) + xi*ImPW(2,0,i)
         PW(4,0,i) = RePW(4,0,i) + xi*ImPW(4,0,i)
         do L = 1, 7
            read (SAIDnum(4*L+0-1),*) plab(i), tmp(0), tmp(1), tmp(2), !loadin higher waves - 4 IJ combinations
     &           tmp(3), RePW(1,L,i), ImPW(1,L,i), tmp(4), tmp(5)
            read (SAIDnum(4*L+1-1),*) plab(i), tmp(0), tmp(1), tmp(2),
     &           tmp(3), RePW(2,L,i), ImPW(2,L,i), tmp(4), tmp(5)
            read (SAIDnum(4*L+2-1),*) plab(i), tmp(0), tmp(1), tmp(2),
     &           tmp(3), RePW(3,L,i), ImPW(3,L,i), tmp(4), tmp(5)
            read (SAIDnum(4*L+3-1),*) plab(i), tmp(0), tmp(1), tmp(2),
     &           tmp(3), RePW(4,L,i), ImPW(4,L,i), tmp(4), tmp(5)
            PW(1,L,i) = RePW(1,L,i) + xi*ImPW(1,L,i)
            PW(2,L,i) = RePW(2,L,i) + xi*ImPW(2,L,i)
            PW(3,L,i) = RePW(3,L,i) + xi*ImPW(3,L,i)
            PW(4,L,i) = RePW(4,L,i) + xi*ImPW(4,L,i)
            PW(0,0,i) = plab(i)/1000.d0
         enddo
      enddo

      do i=1,nSAID
         close(SAIDnum(i))
      enddo
      return
      end subroutine SAID_PiN_PW


      subroutine legendreP(z,pol,pol1p)
      implicit none
      real (8) :: z
      real (8), dimension (0:9) :: pol, pol1p
      pol(0) = 1.d0
      pol(1) = z
      pol(2) = (3*z*z - 1.d0)/2.d0
      pol(3) = (-3*z + 5*z*z*z)/2.d0
      pol(4) = (3.d0 - 30*z*z + 35*z**4)/8.d0
      pol(5) = (15*z - 70*z*z*z + 63*z**5)/8.d0
      pol(6) = (-5.d0 +105*z*z - 315*z**4 + 231*z**6)/16.d0
      pol(7) = (-35*z + 315*z*z*z - 693*z**5 + 429*z**7 )/16.d0
      pol(8) = (35.d0 - 1260*z*z + 6934*z**4 - 12012*z**6 + 
     &     6435*z**8)/128.d0
      pol(9) = (315*z - 4620*z*z*z + 18018*z**5 - 25740*z**7 + 
     &     12155*z**9 )/129.d0

      pol1p(0) = 0.d0
      pol1p(1) = 1.d0
      pol1p(2) = 3.d0*z
      pol1p(3) = (-3.d0 + 15*z*z)/2.d0
      pol1p(4) = (-60*z + 140*z*z*z)/8.d0
      pol1p(5) = (15.d0 - 210*z*z + 315*z**4)/8.d0
      pol1p(6) = (210*z - 1260*z*z*z + 1386*z**5)/16.d0
      pol1p(7) = (-35.d0 + 945*z*z - 3465*z**4 + 3003*z**6)/16.d0
      pol1p(8) = (-2520*z + 27720*z*z*z - 72072*z**5 + 
     &     51480*z**7)/128.d0
      pol1p(9) = (315.d0 - 13860*z*z + 90090*z**4 - 180180*z**6 + 
     &     109395*z**8)/128.d0
      return
      end subroutine legendreP
      
      !       LB envelope routines
      subroutine AMP12_SAID(s,t,tpi,PW,A12,B12,lfix)
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A12,B12
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      integer, intent(in):: lfix
      
      call SAID_ABC(s,t,tpi,PW,ABC,lfix)
      A12=ABC(0,0)+2.d0*ABC(1,0)
      B12=ABC(0,1)+2.d0*ABC(1,1)
      end
      
      subroutine AMP32_SAID(s,t,tpi,PW,A32,B32,lfix)
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A32,B32
      complex (8),intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      integer, intent(in):: lfix
      
      call SAID_ABC(s,t,tpi,PW,ABC,lfix)
      A32=ABC(0,0)-ABC(1,0)
      B32=ABC(0,1)-ABC(1,1)
      end
      
      subroutine AMP12_SAID_thcor(s,t,tpi,PW,A12,B12,lfix)
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A12,B12
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      integer, intent(in):: lfix
      
      call SAID_ABC_thcor(s,t,tpi,PW,ABC)
      A12=ABC(0,0)+2.d0*ABC(1,0)
      B12=ABC(0,1)+2.d0*ABC(1,1)
      end
      
      subroutine AMP32_SAID_thcor(s,t,tpi,PW,A32,B32,lfix)
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A32,B32
      complex (8),intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      integer, intent(in):: lfix
      
      call SAID_ABC_thcor(s,t,tpi,PW,ABC,lfix)
      A32=ABC(0,0)-ABC(1,0)
      B32=ABC(0,1)-ABC(1,1)
      end
      
      end module
