module swave

real(8), parameter:: msigma=0.6d0, gamma_sigma=.45d0, mf0980=0.990d0,gamma_f0980=0.055d0, mf01370=1.35d0, &
gamma_f01370=0.35d0

contains


complex(8) function directS(epsl,q,p1,p2,k1,lambda1,lambda2)
    use parameters
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    complex(8),dimension(4,4) :: matrix
    complex(8),dimension(4) :: spin1,spin2,u2
    integer, intent(in) :: lambda1,lambda2
    real(8) :: fst,s,s12,t,A,alpha,f1,fv,t0,mu0,ilo_sk_4d
    complex(8) :: w(0:3),bw,regge,cilo_sk_4d
    integer :: ii,jj
    real(8) :: theta1,phi1,theta2,phi2,p1mod,p2mod,v3abs

    real(8) :: a1,a2,delta,k2(0:3),d(0:3)

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

    bw=1.d0/(msigma**2-s12-xi*msigma*gamma_sigma)+1d0/(mf01370**2-s12-xi*mf01370*gamma_f01370)

    d=p1-p2
    w=(q*cilo_sk_4d(epsl,dcmplx(d))-epsl*ilo_sk_4d(q,d))
    call v4xgam(w,matrix)
!     matrix=sl(w)!A*identity+B*sl(k1+k2)/2.0_dp
!     spin1=spinor(p1,mN,lambda1,1)
!     spin2=spinor(p2,mN,lambda2,2)

    ! Sum over dirac indices
    directS=0.0d0
    do jj=1,4
       do ii=1,4
          directS=directS+spin2(ii)*matrix(ii,jj)*spin1(jj)
       end do
    end do

    alpha=0.55+0.8*t
    regge=(-1+exp(-xi*pi*alpha))*(s/1.0d0)**alpha !/(sin(pi_dp*alpha))

    directS=regge/s*directS*bw


  end function directS

  complex(8) function directSRes(epsl,q,p1,p2,k1,lambda1,lambda2,RM,RG)
    use parameters
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    integer, intent(in) :: lambda1,lambda2
    real(8), intent(in):: RM,RG
    complex(8),dimension(4,4) :: matrix
    complex(8),dimension(4) :: spin1,spin2,u2
    real(8) :: fst,s,s12,t,A,alpha,f1,fv,t0,mu0,ilo_sk_4d
    complex(8) :: w(0:3),bw,regge,cilo_sk_4d
    integer :: ii,jj
    real(8) :: theta1,phi1,theta2,phi2,p1mod,p2mod,v3abs

    real(8) :: a1,a2,delta,k2(0:3),d(0:3)

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

    bw=1.d0/(RM**2-s12-xi*RM*RG)

    d=p1-p2
    w=(q*cilo_sk_4d(epsl,dcmplx(d))-epsl*ilo_sk_4d(q,d))
    call v4xgam(w,matrix)
!     matrix=sl(w)!A*identity+B*sl(k1+k2)/2.0_dp
!     spin1=spinor(p1,mN,lambda1,1)
!     spin2=spinor(p2,mN,lambda2,2)

    ! Sum over dirac indices
    directSRes=0.0d0
    do jj=1,4
       do ii=1,4
          directSRes=directSRes+spin2(ii)*matrix(ii,jj)*spin1(jj)
       end do
    end do

    alpha=0.55+0.8*t
    regge=(-1.+exp(-xi*pi*alpha))*(s/1.0d0)**alpha !/(sin(pi_dp*alpha))

    directSRes=regge/s*directSRes*bw


  end function directSRes

    complex(8) function directS500(epsl,q,p1,p2,k1,lambda1,lambda2)
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    integer, intent(in) :: lambda1,lambda2

    directS500=directSRes(epsl,q,p1,p2,k1,lambda1,lambda2,msigma,gamma_sigma)
    end function

    complex(8) function directS980(epsl,q,p1,p2,k1,lambda1,lambda2)
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    integer, intent(in) :: lambda1,lambda2

    directS980=directSRes(epsl,q,p1,p2,k1,lambda1,lambda2,mf0980,gamma_f0980)
    end function

    complex(8) function directS1370(epsl,q,p1,p2,k1,lambda1,lambda2)
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    integer, intent(in) :: lambda1,lambda2

    directS1370=directSRes(epsl,q,p1,p2,k1,lambda1,lambda2,mf01370,gamma_f01370)
    end function

    complex(8) function directS_nr(epsl,q,p1,p2,k1,lambda1,lambda2)
    use parameters
    implicit none
    real(8),dimension(0:3), intent(in) :: p1,q,p2,k1
    complex(8),dimension(0:3), intent(in) :: epsl
    complex(8),dimension(4,4) :: matrix
    complex(8),dimension(4) :: spin1,spin2,u2
    integer, intent(in) :: lambda1,lambda2
    real(8) :: fst,s,s12,t,A,alpha,f1,fv,t0,mu0,ilo_sk_4d,s1min,s1max,factor
    complex(8) :: w(0:3),bw,regge,cilo_sk_4d
    integer :: ii,jj
    real(8) :: theta1,phi1,theta2,phi2,p1mod,p2mod,v3abs

    real(8) :: a1,a2,delta,k2(0:3),d(0:3)

    k2=q+p1-k1-p2

    s=ilo_sk_4d(p1+q,p1+q)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s12=ilo_sk_4d(k1+k2,k1+k2)

    s1min=4.*mpi**2
    s1max=(sqrt(s)-mp)**2
    factor=(s1min-s12)*(s1max-s12)

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)

    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(spin1,p1mod,theta1,phi1,lambda1)
    call bispinor(u2,p2mod,theta2,phi2,lambda2)
    call dirconj(u2,spin2)

!     bw=1d0/(msigma**2-s12-xi*msigma*gamma_sigma)+1d0/(mf01370**2-s12-xi*mf01370*gamma_f01370)

    d=p1-p2
    w=(q*cilo_sk_4d(epsl,dcmplx(d))-epsl*ilo_sk_4d(q,d))
    call v4xgam(w,matrix)
!     matrix=sl(w)!A*identity+B*sl(k1+k2)/2.0_dp
!     spin1=spinor(p1,mN,lambda1,1)
!     spin2=spinor(p2,mN,lambda2,2)

    ! Sum over dirac indices
    directS_nr=0.0d0
    do jj=1,4
       do ii=1,4
          directS_nr=directS_nr+spin2(ii)*matrix(ii,jj)*spin1(jj)
       end do
    end do

    alpha=0.55+0.8*t
    regge=(-1.+exp(-xi*pi*alpha))*(s/1.0d0)**alpha !/(sin(pi_dp*alpha))

    directS_nr=regge/s*directS_nr*factor


  end function directS_nr

  complex(8) function directSRes_MHel(s,t,mpipi,l,l1,l2,RM,RG)
    use parameters
    implicit none
    real(8),intent(in) :: s,t,mpipi
    integer, intent(in) :: l,l1,l2
    complex(8),dimension(0:3):: epsl

    real(8), intent(in):: RM,RG
    complex(8),dimension(4,4) :: matrix
    complex(8),dimension(4) :: spin1,spin2,u2
    real(8) :: fst,s1,A,alpha,f1,fv,t0,mu0,ilo_sk_4d
    complex(8) :: w(0:3),bw,regge,cilo_sk_4d
    integer :: ii,jj
    real(8) :: p1mod,p2mod,qmod,costhq,sinthq,costh1,E1,E2

    real(8) :: a1,a2,delta,q(0:3),d(0:3),p1(0:3),p2(0:3),k(0:3)


    s1=mpipi**2
    qmod=(s1-t)/(2.d0*dsqrt(s1))
    E1=(s-mp**2+t)/(2.d0*dsqrt(s1))
    E2=(s-mp**2-s1)/(2.d0*dsqrt(s1))
    p1mod=dsqrt(E1**2-mp**2)
    p2mod=dsqrt(E2**2-mp**2)
    costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
    sinthq=dsqrt(1.d0-costhq**2)
    costh1=(mp**2-E1*E2-t/2.d0)/(p1mod*p2mod)

    call bispinor(spin1,p1mod,dacos(costh1),0.d0,l1)
    call bispinor(u2,p2mod,pi,0.d0,l2)
    call dirconj(u2,spin2)

    bw=1.d0/(RM**2-s1-xi*RM*RG)

    q=qmod*[1.d0,-sinthq,0.d0,costhq]
    epsl=-dble(l)/sqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),xi*dble(l),dcmplx(sinthq)]
    p2=[E2,0.d0,0.d0,-p2mod]
    k=[mpipi,0.d0,0.d0,0.d0]
    p1=k+p2-q

    d=p1-p2
    w=(q*cilo_sk_4d(epsl,dcmplx(d))-epsl*ilo_sk_4d(q,d))
    call v4xgam(w,matrix)


    ! Sum over dirac indices
    directSRes_MHel=0.0d0
    do jj=1,4
       do ii=1,4
          directSRes_MHel=directSRes_MHel+spin2(ii)*matrix(ii,jj)*spin1(jj)
       end do
    end do

    alpha=0.55+0.8*t
    regge=(-1.+exp(-xi*pi*alpha))*(s/1.0d0)**alpha !/(sin(pi_dp*alpha))

    directSRes_MHel=regge/s*directSRes_MHel*bw


  end function directSRes_MHel

  complex(8) function directS500_MHel(s,t,mpipi,l,l1,l2)
  implicit none
   real(8),intent(in) :: s,t,mpipi
   integer, intent(in) :: l,l1,l2
   directS500_MHel=directSRes_MHel(s,t,mpipi,l,l1,l2,msigma,gamma_sigma)
  end function

  complex(8) function directS980_MHel(s,t,mpipi,l,l1,l2)
  implicit none
   real(8),intent(in) :: s,t,mpipi
   integer, intent(in) :: l,l1,l2
   directS980_MHel=directSRes_MHel(s,t,mpipi,l,l1,l2,mf0980,gamma_f0980)
  end function

  complex(8) function directS1370_MHel(s,t,mpipi,l,l1,l2)
  implicit none
   real(8),intent(in) :: s,t,mpipi
   integer, intent(in) :: l,l1,l2
   directS1370_MHel=directSRes_MHel(s,t,mpipi,l,l1,l2,mf01370,gamma_f01370)
  end function

  complex(8) function directS_nr_MHel(s,t,mpipi,l,l1,l2)
   use parameters
    implicit none
    real(8),intent(in) :: s,t,mpipi
    integer, intent(in) :: l,l1,l2
    complex(8),dimension(0:3):: epsl
    complex(8),dimension(4,4) :: matrix
    complex(8),dimension(4) :: spin1,spin2,u2
    real(8) :: fst,s1,A,alpha,f1,fv,t0,mu0,ilo_sk_4d,s1max,s1min
    complex(8) :: w(0:3),bw,regge,cilo_sk_4d
    integer :: ii,jj
    real(8) :: p1mod,p2mod,qmod,costhq,sinthq,costh1,E1,E2,factor

    real(8) :: a1,a2,delta,q(0:3),d(0:3),p1(0:3),p2(0:3),k(0:3)


    s1=mpipi**2
    qmod=(s1-t)/(2.d0*dsqrt(s1))
    E1=(s-mp**2+t)/(2.d0*dsqrt(s1))
    E2=(s-mp**2-s1)/(2.d0*dsqrt(s1))
    p1mod=dsqrt(E1**2-mp**2)
    p2mod=dsqrt(E2**2-mp**2)
    costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
    sinthq=dsqrt(1.d0-costhq**2)
    costh1=(mp**2-E1*E2-t/2.d0)/(p1mod*p2mod)

    call bispinor(spin1,p1mod,dacos(costh1),0.d0,l1)
    call bispinor(u2,p2mod,pi,0.d0,l2)
    call dirconj(u2,spin2)

    q=qmod*[1.d0,-sinthq,0.d0,costhq]
    epsl=-dble(l)/sqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),xi*dble(l),dcmplx(sinthq)]
    p2=[E2,0.d0,0.d0,-p2mod]
    k=[mpipi,0.d0,0.d0,0.d0]
    p1=k+p2-q

    d=p1-p2
    w=(q*cilo_sk_4d(epsl,dcmplx(d))-epsl*ilo_sk_4d(q,d))
    call v4xgam(w,matrix)

    s1min=4.*mpi**2
    s1max=(sqrt(s)-mp)**2
    factor=(s1min-s1)*(s1max-s1)

    ! Sum over dirac indices
    directS_nr_MHel=0.0d0
    do jj=1,4
       do ii=1,4
          directS_nr_MHel=directS_nr_MHel+spin2(ii)*matrix(ii,jj)*spin1(jj)
       end do
    end do

    alpha=0.55+0.8*t
    regge=(-1.+exp(-xi*pi*alpha))*(s/1.0d0)**alpha !/(sin(pi_dp*alpha))

    directS_nr_MHel=regge/s*directS_nr_MHel*factor


  end function directS_nr_MHel


end module swave
