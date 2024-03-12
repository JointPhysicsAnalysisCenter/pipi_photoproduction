module dwave

real(8), parameter:: mt=1.2755d0,gamma_t=.1867d0

contains

complex(8) function directD(epsl,q,p1,p2,k1,lambda1,lambda2)
    use parameters
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    complex(8),dimension(1:4,1:4) :: w1sl,w2sl
    complex(8),dimension(1:4) :: spin1,spin2,u2
    integer, intent(in) :: lambda1,lambda2
    real(8) :: fst,s,s12,t,A,alpha,f1,fv,t0,mu0,ilo_sk_4d,qdotp1_p2
    complex(8) :: bw,regge,cilo_sk_4d,w1(0:3),w2(0:3),epsdotk1_k2,d(0:3)
    integer :: ii,jj
    real(8) :: b1,theta1,phi1,theta2,phi2,p1mod,p2mod,v3abs,k1_k2sq

    real(8) :: a1,a2,delta,k2(0:3)

    k2=q+p1-k1-p2

    s=ilo_sk_4d(p1+q,p1+q)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s12=ilo_sk_4d(k1+k2,k1+k2)

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)

    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(spin1,p1mod,theta1,phi1,lambda1)
    call bispinor(u2,p2mod,theta2,phi2,lambda2)
    call dirconj(u2,spin2)

     k1_k2sq=ilo_sk_4d(k1-k2,k1-k2)
     qdotp1_p2=ilo_sk_4d(q,p1-p2)
     epsdotk1_k2=cilo_sk_4d(epsl,dcmplx(k1)-dcmplx(k2))

     w1=-q*(epsdotk1_k2*ilo_sk_4d(p1-p2,k1-k2)-cilo_sk_4d(epsl,dcmplx(p1-p2))*k1_k2sq/3.d0&
     +ilo_sk_4d(p1-p2,k1+k2)*cilo_sk_4d(epsl,dcmplx(k1)+dcmplx(k2))*k1_k2sq/s12/3.d0)&
     +epsl*(ilo_sk_4d(q,k1-k2)*ilo_sk_4d(p1-p2,k1-k2)-&
     qdotp1_p2*k1_k2sq/3.d0+ilo_sk_4d(p1-p2,k1+k2)*ilo_sk_4d(q,k1+k2)*k1_k2sq/(3.d0*s12))

     w2=qdotp1_p2*(epsdotk1_k2*(k1-k2)-epsl*k1_k2sq/3.d0+(k1+k2)*cilo_sk_4d(epsl,dcmplx(k1)+dcmplx(k2))*k1_k2sq/(3.d0*s12))-&
     cilo_sk_4d(epsl,dcmplx(p1)-dcmplx(p2))*(ilo_sk_4d(q,k1-k2)*(k1-k2)-q*k1_k2sq/3.d0+&
     (k1+k2)*ilo_sk_4d(q,k1+k2)*k1_k2sq/(3.d0*s12))

    bw=1.d0/(mt**2-s12-xi*mt*gamma_t)

    call v4xgam(w1,w1sl)
    call v4xgam(w2,w2sl)

    ! Sum over dirac indices
    directD=dcmplx(0.0d0)
    do jj=1,4
       do ii=1,4
          directD=directD+spin2(ii)*(w1sl(ii,jj)+w2sl(ii,jj))*spin1(jj)
       end do
    end do

    alpha=0.55d0+0.8d0*t
    regge=(-1.0d0+exp(-xi*pi*alpha))*(s/1.0d0)**alpha /(sin(pi*alpha))

!     directD=1.0/
    directD=regge/s*directD*bw

  end function directD

  function tau(v,M) result(val)
  use parameters
  implicit none
  complex(8), intent(in):: v(0:3)
  integer, intent(in):: M
  complex(8):: val(1:3)
  select case (M)
  case (2)
    val=[v(1)-xi*v(2),-(xi*v(1)+v(2)),dcmplx(0.d0)]*dsqrt(2.d0*pi/15.d0)
  case (-2)
    val=[v(1)+xi*v(2),(xi*v(1)-v(2)),dcmplx(0.d0)]*dsqrt(2.d0*pi/15.d0)
  case (1)
    val=[-v(3),xi*v(3),-(v(1)-xi*v(2))]*dsqrt(2.d0*pi/15.d0)
  case (-1)
    val=[v(3),xi*v(3),(v(1)+xi*v(2))]*dsqrt(2.d0*pi/15.d0)
  case (0)
    val=[-v(1)*dsqrt(4.d0*pi/5.d0)/3.d0,-v(2)*dsqrt(4.d0*pi/5.d0)/3.d0,v(3)*dsqrt(16.d0*pi/5.d0)/3.d0]
  case default
    print*,"Wrong selector index"
    stop
  end select
  end function

  complex(8) function W1M(k,v,dp,M) result(val)
  use parameters
  implicit none
  real(8), intent(in):: k
  complex(8), intent(in):: v(0:3),dp(0:3)
  integer, intent(in):: M
  select case (M)
  case (2)
    val=(v(1)*dp(1)-xi*v(1)*dp(2)-xi*v(2)*dp(1)-v(2)*dp(2))*dsqrt(2.d0*pi/15.d0)
  case (-2)
    val=(v(1)*dp(1)+xi*v(1)*dp(2)+xi*v(2)*dp(1)-v(2)*dp(2))*dsqrt(2.d0*pi/15.d0)
  case (1)
    val=(-v(1)*dp(3)-v(3)*dp(1)+xi*v(2)*dp(3)+xi*v(3)*dp(2))*dsqrt(2.d0*pi/15.d0)
  case (-1)
    val=(v(1)*dp(3)+xi*v(2)*dp(3)+v(3)*dp(1)+xi*v(3)*dp(2))*dsqrt(2.d0*pi/15.d0)
  case (0)
    val=-(v(1)*dp(1)+v(2)*dp(2)-2.d0*v(3)*dp(3))*dsqrt(4.d0*pi/5.d0)/3.d0
  case default
    print*,"Wrong selector index"
    stop
  end select
    val=val*4.d0*k**2
  end function

  complex(8) function directDM(epsl,q,p1,p2,k1,lambda1,lambda2,M)
    use parameters
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    complex(8),dimension(1:4,1:4) :: w1sl,w2sl
    complex(8),dimension(1:4) :: spin1,spin2,u2
    integer, intent(in) :: lambda1,lambda2,M
    real(8) :: fst,s,s12,t,A,alpha,f1,fv,t0,mu0,ilo_sk_4d,k
    complex(8) :: bw,regge,cilo_sk_4d,w1(0:3),w2(0:3),dp(0:3),W2vec(1:3)
    integer :: ii,jj
    real(8) :: b1,theta1,phi1,theta2,phi2,p1mod,p2mod,v3abs

    real(8) :: a1,a2,delta,k2(0:3)

    k2=q+p1-k1-p2

    s=ilo_sk_4d(p1+q,p1+q)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s12=ilo_sk_4d(k1+k2,k1+k2)

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k=v3abs(k1)
    dp=dcmplx(p1-p2)

    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(spin1,p1mod,theta1,phi1,lambda1)
    call bispinor(u2,p2mod,theta2,phi2,lambda2)
    call dirconj(u2,spin2)

     w1=-(q*W1M(k,epsl,dp,M)-epsl*W1M(k,dcmplx(q),dp,M))

     W2vec=-4.d0*k**2*(ilo_sk_4d(q,p1-p2)*tau(epsl,M)-cilo_sk_4d(epsl,dcmplx(p1)-dcmplx(p2))*tau(dcmplx(q),M))

     w2=[dcmplx(0.d0),W2vec]

!      print*,W2vec
!      print*,w2

    bw=1.d0/(mt**2-s12-xi*mt*gamma_t)

    call v4xgam(w1,w1sl)
    call v4xgam(w2,w2sl)

    ! Sum over dirac indices
    directDM=dcmplx(0.0d0)
    do jj=1,4
       do ii=1,4
          directDM=directDM+spin2(ii)*(w1sl(ii,jj)+w2sl(ii,jj))*spin1(jj)
       end do
    end do

    alpha=0.55d0+0.8d0*t
    regge=(-1.0d0+exp(-xi*pi*alpha))*(s/1.0d0)**alpha /(sin(pi*alpha))

!     directD=1.0/
    directDM=regge/s*directDM*bw

  end function directDM

  complex(8) function directD2pl(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      real(8):: theta,phi,costh
      complex(8):: funkcjakulista

      call v3angles(k1,theta,phi)
      costh=cos(theta)
      directD2pl=directDM(eps,q,p1,p2,k1,l1,l2,2)*funkcjakulista(2,2,costh,phi)
  end function

  complex(8) function directD2mi(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      real(8):: theta,phi,costh
      complex(8):: funkcjakulista

      call v3angles(k1,theta,phi)
      costh=cos(theta)
      directD2mi=directDM(eps,q,p1,p2,k1,l1,l2,-2)*funkcjakulista(2,-2,costh,phi)
  end function

  complex(8) function directD1pl(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      real(8):: theta,phi,costh
      complex(8):: funkcjakulista

      call v3angles(k1,theta,phi)
      costh=cos(theta)
      directD1pl=directDM(eps,q,p1,p2,k1,l1,l2,1)*funkcjakulista(2,1,costh,phi)
  end function

  complex(8) function directD1mi(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      real(8):: theta,phi,costh
      complex(8):: funkcjakulista

      call v3angles(k1,theta,phi)
      costh=cos(theta)
      directD1mi=directDM(eps,q,p1,p2,k1,l1,l2,-1)*funkcjakulista(2,-1,costh,phi)
  end function

  complex(8) function directD0(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      real(8):: theta,phi,costh
      complex(8):: funkcjakulista

      call v3angles(k1,theta,phi)
      costh=cos(theta)
      directD0=directDM(eps,q,p1,p2,k1,l1,l2,0)*funkcjakulista(2,0,costh,phi)
  end function

  complex(8) function dwave_anal_LM(s,t,mpipi,l,l1,l2,M)
  use frames
  implicit none
  real(8), intent(in):: s,t,mpipi
  integer, intent(in):: l,l1,l2,M
  real(8):: phi=.2d0,costh=.2d0
  real(8):: q(0:3),p1(0:3),k1(0:3),p2(0:3)
  complex(8):: eps(0:3),funkcjakulista

   call Helicity4Vectors(s,t,mpipi**2,costh,phi,l,eps,q,p1,k1,p2)

   dwave_anal_LM=directDM(eps,q,p1,p2,k1,l1,l2,M)
   end function
end module dwave

module dwave_test
real(8), private:: ps,pt,pspipi,pcosth
integer, private:: pl,pl1,pl2,pJ,pM
real(8), private, parameter:: pi=acos(-1.)
contains

complex(8) function dwave_complete(s,t,mpipi,costh,phi,l,l1,l2)
use dwave
use frames
implicit none
real(8), intent(in):: s,t,mpipi,costh,phi
integer, intent(in):: l,l1,l2

real(8):: q(0:3),p1(0:3),k1(0:3),p2(0:3)
complex(8):: eps(0:3)

   call Helicity4Vectors(s,t,mpipi**2,costh,phi,l,eps,q,p1,k1,p2)

   dwave_complete=directD(eps,q,p1,p2,k1,l1,l2)
end function


complex(8) function dwave_LM(s,t,mpipi,l,l1,l2,J,M)
implicit none
real(8), intent(in):: s,t,mpipi
integer, intent(in):: l,l1,l2,J,M
complex(8):: gaussc

   ps=s
   pt=t
   pspipi=mpipi**2
   pl=l
   pl1=l1
   pl2=l2
   pJ=J
   pM=M
   dwave_LM=gaussc(dwave_costh,-1.d0,1.d0,1.d-3)
end function

complex(8) function dwave_costh(costh)
implicit none
real(8), intent(in):: costh
complex(8):: gaussc
   pcosth=costh
!    print*,pi
   dwave_costh=gaussc(dwave_phi,0.d0,2.d0*pi,1.d-3)
end function

complex(8) function dwave_phi(phi)
use frames
use dwave
implicit none
real(8), intent(in):: phi
real(8):: q(0:3),p1(0:3),k1(0:3),p2(0:3)
complex(8):: eps(0:3),funkcjakulista

   call Helicity4Vectors(ps,pt,pspipi,pcosth,phi,pl,eps,q,p1,k1,p2)

   dwave_phi=directD(eps,q,p1,p2,k1,pl1,pl2)*dconjg(funkcjakulista(pJ,pM,pcosth,phi))

end function

subroutine doDwaveTest(unit,s,t,mpipi,Jmax)
use dwave
implicit none
real(8), intent(in):: s,t,mpipi
integer, intent(in):: Jmax,unit
integer:: l,m,lgam,l1,l2
complex(8):: a,b
   do lgam=-1,1,2
   do l=1,Jmax
      print*,l
      do m=-l,l
         do l1=-1,1,2
            do l2=-1,1,2
              a=dwave_LM(s,t,mpipi,lgam,l1,l2,l,M)
              write(unit,"(""n: lgam="",I2,2X,""l="",I2,2X,""m="",I2,2X,""l1="",I2,2X,""l2="",I2,2X,""Re="",F12.3,""Im="",F12.3)")&
              lgam,l,m,l1,l2,dble(a),dimag(a)
              if(l.eq.2) then
              b=dwave_anal_LM(s,t,mpipi,lgam,l1,l2,M)
              write(unit,"(""a: lgam="",I2,2X,""l="",I2,2X,""m="",I2,2X,""l1="",I2,2X,""l2="",I2,2X,""Re="",F12.3,""Im="",F12.3)")&
              lgam,l,m,l1,l2,dble(b),dimag(b)
              endif
              flush(unit)
            enddo
         enddo
      enddo
      write(unit,*)"-----"
   enddo
   enddo
end subroutine
end module dwave_test
