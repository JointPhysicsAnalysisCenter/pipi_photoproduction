module pwave
    implicit none
    real(8):: pt,ps,pmpipi,pcosth,pphi
    real(8),parameter:: eps=10.d0**(-4)
    private:: pt,ps,pmpipi,pcosth,pphi,eps
    logical:: free
contains
complex(8) function directPPom(eps,q,p1,p2,k1,l1,l2)
    implicit none
    complex(8),intent(in)::eps(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3)
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0P=1.23d0
    
    integer::i
    
    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(eps,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    k=k1+k2
    v=cmplx(k)*epsxk1_k2-cmplx(qxk1_k2)*eps

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)
    
    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    result=-result*f(s,t,BP)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,aP,aPprime)/s
!     result=-result*F1(t,t0P)*Fv(t,mrho,mu0)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,aP,aPprime)/s
!     print*,result
    directPpom=result
end function

complex(8) function directPPom_M(eps,q,p1,p2,k1,l1,l2,M)
use parameters
    implicit none
    complex(8),intent(in)::eps(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2,M
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1,costhq,qmod,k1mod
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3)
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0P=1.23d0

    integer::i,l

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k1mod=v3abs(k1)
    qmod=v3abs(q)
    costhq=q(3)/qmod
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(eps,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    l=-nint(sqrt(2.)*real(eps(1))/costhq)
    k=k1+k2
    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*eps)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)

    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    result=-result*f(s,t,BP)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,aP,aPprime)/s
    directPPom_M=result
end function

complex(8) function directPPom_MHel(s,t,mpipi,l,l1,l2,M)
use parameters
use frames
    implicit none
    real(8), intent(in)::s,t,mpipi
    integer,intent(in):: l,l1,l2,M
    complex(8)::eps(0:3)
    real(8)::costhq,sinthq,k(0:3),s1,qmod,k1mod,omega,E1,E2,costh1
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0P=1.23d0

    integer::i

    s1=mpipi**2
    qmod=(s1-t)/(2.d0*dsqrt(s1))
    E1=(s-mp**2+t)/(2.d0*dsqrt(s1))
    E2=(s-mp**2-s1)/(2.d0*dsqrt(s1))
    p1mod=dsqrt(E1**2-mp**2)
    p2mod=dsqrt(E2**2-mp**2)
    omega=dsqrt(s1)/2.d0
    k1mod=dsqrt(omega**2-mpi**2)
    costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
    sinthq=dsqrt(1.d0-costhq**2)
    costh1=(mp**2-E1*E2-t/2.d0)/(p1mod*p2mod)

    call bispinor(u1,p1mod,dacos(costh1),0.d0,l1)
    call bispinor(u2,p2mod,pi,0.d0,l2)
    call dirconj(u2,u2bar)
    k=[sqrt(s1),0.d0,0.d0,0.d0]
    eps=-dble(l)/sqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),xi*dble(l),dcmplx(sinthq)]

    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*eps)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)

    result=-result*f(s,t,BP)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,aP,aPprime)/s
    directPPom_MHel=result
end function

complex(8) function directPf2(epsl,q,p1,p2,k1,l1,l2)
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.d0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1,alfa
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3)
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4
    complex(8),parameter::ci=cmplx(0.d0,1.d0)
    real(8),parameter:: pi=dacos(-1.d0)
    
    integer::i
    
    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(epsl,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    k=k1+k2
    v=cmplx(k)*epsxk1_k2-cmplx(qxk1_k2)*epsl

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)
    
    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    if (.not. free) then
        result=-result*f(s,t,Bf2)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,af2,af2prime)/s
    else
        alfa=af2+af2prime*t
        result=-result*BWM(s1,.775d0,.149d0)*(1.+exp(-ci*pi*alfa))*(s/s0)**alfa/s
    endif
    directPf2=result
end function

complex(8) function BW(s,M,GAM)
    implicit none
    real(8),intent(in)::s,M,GAM
    BW=1.d0/(M**2-s-cmplx(0.d0,1.d0)*M*GAM)
end function

complex(8) function BWM(s,M,GAM)
    use parameters
    implicit none
    real(8),intent(in)::s,M,GAM
    real(8)::GAMM,qm,qs
    
    GAMM=GAM*(q(sqrt(s))/q(M))**3*2.d0/(1.d0+(q(sqrt(s))/q(M))**2)
    
    BWM=1.d0/(M**2-s-cmplx(0.d0,1.d0)*M*GAMM)
    
!     BWM=1.d0/(M**2-s-cmplx(0.d0,1.d0)*M*GAM)
    
    
    contains
        real(8) function q(m)
        implicit none
            real(8), intent(in)::m
            q=sqrt(m**2/4.d0-mpi**2)
        end function
end function

real(8) function f(s,t,b)
    implicit none
    real(8),intent(in)::s,t,b
    
!     f=exp(A*t/2.d0)
    f=exp(b*t)
end function

real(8) function F1(t,t0)
use parameters
implicit none
real(8), intent(in):: t,t0
    F1=(4.d0*mp**2-2.8d0*t)/(4.d0*mp**2-t)/(1.d0-t/t0)**2
end function

real(8) function FV(t,mv,mu0)
implicit none
real(8), intent(in):: t,mv,mu0
    FV=(2.d0*mu0**2+mv**2)/(1.d0-t/mv**2)/(2.d0*mu0**2+mv**2-t)
end function

complex(8) function R(s,t,a,aprime)
    implicit none
    real(8),intent(in)::s,t,a,aprime
    complex(8),parameter::ci=cmplx(0.d0,1.d0)
    real(8),parameter:: pi=dacos(-1.d0),s0=1.d0
    real(8)::alfa
    alfa=a+aprime*t
    R=(1.d0+exp(-ci*pi*alfa))/(sin(pi*alfa))*(s/s0)**alfa
end function

complex(8) function RN(s,t,a,aprime)
    implicit none
    real(8),intent(in)::s,t,a,aprime
    complex(8),parameter::ci=cmplx(0.d0,1.d0)
    real(8),parameter:: pi=dacos(-1.d0),s0=1.d0
    real(8)::alfa
    alfa=a+aprime*t
    RN=(alfa/a)*(1.d0+exp(-ci*pi*alfa))/(sin(pi*alfa))*(s/s0)**alfa
end function

real(8) function dsdt(s,t)
    use parameters
    implicit none
    real(8),intent(in):: s,t
    real(8):: mmin,mmax,gauss
    pt=t
    ps=s
!     mmin=sqrt(4.d0*mpi**2)
!     mmax=sqrt((mp**2*t + s*t)/(2.d0*mp**2) + Sqrt((-4.d0*mp**6*t + 8.d0*mp**4*s*t - 4.d0*mp**2*s**2*t + &
!             mp**4*t**2 - 2.d0*mp**2*s*t**2 + s**2*t**2)/mp**4)/2.d0)
    mmin=0.4d0
    mmax=1.2d0
    dsdt=gauss(dsdtdm,mmin,mmax,eps)
end function

real(8) function dsdtdm(mpipi)
    use parameters
    implicit none
    real(8),intent(in):: mpipi
    real(8):: gauss
    pmpipi=mpipi
    dsdtdm=gauss(dsdtdmdcosth,-1.d0,1.d0,eps)
end function

real(8) function dsdtdmdcosth(costh)
    use parameters
    implicit none
    real(8),intent(in):: costh
    real(8)::gauss
    pcosth=costh
    dsdtdmdcosth=gauss(dsdtdmdcosthdphi,0.d0,2.d0*pi,eps)
end function
    
real(8) function dsdtdmdcosthdphi(phi) result(res)
    use parameters
    use frames
    implicit none
    real(8),intent(in):: phi
    real(8)::s1,q(0:3),p1(0:3),k1(0:3),p2(0:3),k
    complex(8)::epsl(0:3)
    integer:: l1,l2
    res=0.d0
    s1=pmpipi**2
    k=sqrt(s1/4.d0-mpi**2)
    
    call Helicity4Vectors(ps,pt,s1,pcosth,phi,+1,epsl,q,p1,k1,p2)
    do l1=-1,1,2
        do l2=-1,1,2
            res=res+389.d0*k/((2.d0*pi)**3*4.d0*pi*(ps-mp**2)**2*8.d0*2.d0)* &
            abs(directPpom(epsl,q,p1,p2,k1,l1,l2)+directPf2(epsl,q,p1,p2,k1,l1,l2))**2
!               abs(directPpom(epsl,q,p1,p2,k1,l1,l2))**2
!             print*,directP(epsl,q,p1,p2,k1,l1,l2)
!             if (isnan(res)) then
!                 print*,"s1=",s1,"t=",pt,"phi=",phi,"costh=",pcosth,directP(epsl,q,p1,p2,k1,l1,l2)
!             endif
        enddo
    enddo
end function

real(8) function dsdtdmdcosthdphi4D(mpipi,t,costh,phi) result(res)
    use parameters
    use frames
    implicit none
    real(8),intent(in):: mpipi,t,costh,phi
    real(8)::s1,q(0:3),p1(0:3),k1(0:3),p2(0:3),k
    complex(8)::epsl(0:3)
    integer:: l1,l2
    real(8), parameter:: Egam=5.d0
    res=0.d0
    ps=mp**2+2*Egam*mp
    s1=mpipi**2
    k=sqrt(s1/4.d0-mpi**2)
    
    call Helicity4Vectors(ps,t,s1,costh,phi,+1,epsl,q,p1,k1,p2)
    do l1=-1,1,2
        do l2=-1,1,2
            res=res+389.d0*k/((2.d0*pi)**3*4.d0*pi*(ps-mp**2)**2*8.d0*2.d0)* &
            abs(directPpom(epsl,q,p1,p2,k1,l1,l2)+directPf2(epsl,q,p1,p2,k1,l1,l2))**2
!               abs(directPpom(epsl,q,p1,p2,k1,l1,l2))**2
!             print*,directP(epsl,q,p1,p2,k1,l1,l2)
!             if (isnan(res)) then
!                 print*,"s1=",s1,"t=",pt,"phi=",phi,"costh=",pcosth,directP(epsl,q,p1,p2,k1,l1,l2)
!             endif
        enddo
    enddo
end function

complex(8) function directPTMD(epsl,q,p1,p2,k1,l1,l2)
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,spipi
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    OmegaWu1(1:4),u1baru2,k1mk2(0:3),OmegaW(1:4,1:4)
!     couplings
    real(8),parameter:: g=5.96d0,mrho=0.77526d0,t0f2=0.1d-4,mf2=1.2755d0,mp=.93827,frho=5.33d0,&
    GfNN=2.12d0,FfNN=0.d0,fvfgam=0.d0,gvfgam=e/frho*5.76d0
    integer::i,alfa,beta,rho,sigma
    real, dimension(0:3,0:3):: gmunu
    
    gmunu=reshape((/1.,0.,0.,0.,&
                    0.,-1.,0.,0.,&
                    0.,0.,-1.,0.,&
                    0.,0.,0.,-1./),shape(gmunu))
    
    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    k=k1+k2
    
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=dcmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(epsl,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    spipi=ilo_sk_4d(k1+k2,k1+k2)
    directPTMD=dcmplx(0.d0)
    OmegaW=dcmplx(0.d0)
    do alfa=0,3
        do beta=0,3
            do rho=0,3
                do sigma=0,3
                OmegaW=OmegaW+OmegaTensor(alfa,beta,GfNN,FfNN)*W(rho,sigma)*gmunu(alfa,rho)*gmunu(beta,sigma)
                enddo
            enddo
        enddo
    enddo
    
    OmegaWu1=matmul(OmegaW,u1)
    
    directPTMD=-g*f(s,t,4.19d0)*BWM(spipi,.775d0,.149d0)*PRT(s,t)*u1baru2(u2bar,OmegaWu1)*0.02d0
    
    contains
    
    real(8) function gbar(mu,nu)
    implicit none
    integer, intent(in):: mu,nu
    gbar=-gmunu(mu,nu)+(q(mu)-k(mu))*(q(nu)-k(nu))/mf2**2
    end function
    
!     real(8) function P(alfa,beta,rho,sigma)
!     implicit none
!     integer, intent(in)::alfa,beta,rho,sigma
!     P=0.5*(gbar(alfa,rho)*gbar(beta,sigma)+gbar(alfa,sigma)*gbar(beta,rho))-gbar(alfa,beta)*gbar(rho,sigma)/3.d0
!     end function
    
    complex(8) function W(rho,sigma)
    implicit none
    integer, intent(in)::rho,sigma
    
    W=fvfgam/mf2**4*&
    (-cilo_sk_4d(k1mk2,epsl)*ilo_sk_4d(q,k)+cilo_sk_4d(epsl,dcmplx(k))*ilo_sk_4d(q,k1-k2))*(q(rho)+k(rho))*(q(sigma)+k(sigma))+&
    gvfgam/mf2*&
    (cilo_sk_4d(k1mk2,epsl)*(q(rho)+k(rho))*(q(sigma)+k(sigma))-cilo_sk_4d(dcmplx(k),epsl)*(k1(rho)-k2(rho))*(q(sigma)+k(sigma))-&
    cilo_sk_4d(dcmplx(k),epsl)*(q(rho)+k(rho))*(k1(sigma)-k2(sigma))-ilo_sk_4d(k1-k2,q)*epsl(rho)*(q(sigma)+k(sigma))&
    -ilo_sk_4d(k1-k2,q)*(k(rho)+q(rho))*epsl(sigma)+&
    2.*ilo_sk_4d(q,k)*(epsl(rho)*(k1(sigma)-k2(sigma))+epsl(sigma)*(k1(rho)-k2(rho)))&
    )
    end function
    
    function gamma_mu(mu) result(res)
    implicit none
    integer, intent(in):: mu
    complex(8), dimension(1:4,1:4):: res,d0,d1,d2,d3,d5,unit
    common/gammas/d0,d1,d2,d3,d5,unit
    select case (mu)
    case (0)
        res=d0
    case (1)
        res=d1
    case (2)
        res=d2
    case (3)
        res=d3
    case (4)
        res=unit
    case default
        print*,"Lorentz index for gamma array out of range"
    end select
    end function
    
    function OmegaTensor(rho,sigma,G,F) result(res)
    implicit none
    real(8), intent(in):: G,F
    complex(8), dimension(1:4,1:4):: res,p1_p2slash,p1plp2slash
    integer, intent(in):: rho,sigma
   
    call v4xgam(dcmplx(p1-p2),p1_p2slash)
    call v4xgam(dcmplx(p1+p2),p1plp2slash)
    res=G/mp*(&
    -(p1(rho)+p2(rho))*(-gamma_mu(sigma)+p1_p2slash*(p1(sigma)-p2(sigma))/mf2**2)&
    -(p1(sigma)+p2(sigma))*(-gamma_mu(rho)+p1_p2slash*(p1(rho)-p2(rho))/mf2**2)&
    )&
    +(-gmunu(rho,sigma)+(p1(rho)-p2(rho))*(p1(sigma)-p2(sigma))/mf2**2)&
    *(2.*G/(3.*mp)*p1plp2slash+F/(3.*mp**2)*ilo_sk_4d(p1+p2,p1+p2)*gamma_mu(4))&
    +F/mp**2*(p1(rho)+p2(rho))*(p1(sigma)+p2(sigma))*gamma_mu(4)
    end function

    complex(8) function PRT(s,t)
    use parameters
    implicit none
    real(8), intent(in):: s,t
    real(8),parameter:: s0=1.d0
    real(8):: alfaT
    alfaT=af2+af2prime*t
    PRT=pi*af2prime/gamma(alfaT-1.)*(1.+exp(-xi*pi*alfaT))/(2.*sin(pi*alfaT))*(s/s0)*(alfaT-2.0)
    end function
    
end function

! 2023-02-14 helicity projected f2 exchange amplitudes

real(8) function rho(costhq,M)
implicit none
real(8), intent(in):: costhq
integer, intent(in):: M
real(8):: sinthq

    sinthq=sqrt(1.d0-costhq**2)
    rho=0.
    select case(M)
        case (1)
            rho=sinthq/sqrt(2.d0)
        case (0)
            rho=costhq
        case (-1)
            rho=-sinthq/sqrt(2.d0)
        case default
        print*,"Wring M value"
        stop
    end select

end function

real(8) function b(costhq,lbd,M)
implicit none
real(8), intent(in):: costhq
integer, intent(in):: M,lbd
real(8):: sinthq,delta_lbdM
    sinthq=sqrt(1.d0-costhq**2)
    delta_lbdM=0.d0
    if(lbd==M) delta_lbdM=1.d0
    select case (M)
        case (1,-1)
            b=delta_lbdM+(costhq-1.d0)/2.*(-1.d0)**(dble(lbd-M)/2.d0)
        case (0)
            b=-lbd*sinthq/sqrt(2.d0)
    end select
end function

complex(8) function directPf2_M(epsl,q,p1,p2,k1,l1,l2,M)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2,M
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1,k1mod,qmod,costhq,theta,costh,phi,alfa
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3),funkcjakulista
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4
    
    integer::i,l
    
    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(epsl,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    k=k1+k2
    k1mod=v3abs(k1)
    qmod=v3abs(q)
    costhq=q(3)/qmod
    call v3angles(k1,theta,phi)
    costh=cos(theta)
    l=-nint(sqrt(2.)*real(epsl(1))/costhq)
    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*epsl)*funkcjakulista(1,M,costh,phi)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)
    
    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    if (.not. free) then
        result=-result*f(s,t,Bf2)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,af2,af2prime)/s
    else
        alfa=af2+af2prime*t
        result=-result*BWM(s1,.775d0,.149d0)*((1.+exp(-xi*pi*alfa))*(s/s0)**alfa)/s
    endif

    directPf2_M=result
end function

complex(8) function directPf2pw_M(epsl,q,p1,p2,k1,l1,l2,M)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2,M
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1,k1mod,qmod,costhq,theta,costh,phi,alfa
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3),funkcjakulista
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4

    integer::i,l

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(epsl,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    k=k1+k2
    k1mod=v3abs(k1)
    qmod=v3abs(q)
    costhq=q(3)/qmod
    l=-nint(sqrt(2.)*real(epsl(1))/costhq)
    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*epsl)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)

    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    if (.not. free) then
        result=-result*f(s,t,Bf2)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,af2,af2prime)/s
    else
        alfa=af2+af2prime*t
        result=-result*BWM(s1,.775d0,.149d0)*((1.+exp(-xi*pi*alfa))*(s/s0)**alfa)/s
    endif

    directPf2pw_M=result
end function

complex(8) function directPf2pw_MHel(s,t,mpipi,l,l1,l2,M)
    use parameters
    implicit none
    real(8), intent(in)::s,t,mpipi
    integer,intent(in):: l,l1,l2,M
    complex(8)::epsl(0:3)
    real(8):: theta1,theta2,phi1,phi2,p1mod,p2mod,E1,E2,omega,costh1
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::k(0:3),s1,k1mod,qmod,costhq,sinthq,alfa
    complex(8)::u1(1:4),u2(1:4),u2bar(1:4),v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4

    integer::i

    s1=mpipi**2
    qmod=(s1-t)/(2.d0*dsqrt(s1))
    E1=(s-mp**2+t)/(2.d0*dsqrt(s1))
    E2=(s-mp**2-s1)/(2.d0*dsqrt(s1))
    p1mod=dsqrt(E1**2-mp**2)
    p2mod=dsqrt(E2**2-mp**2)
    omega=dsqrt(s1)/2.d0
    k1mod=dsqrt(omega**2-mpi**2)
    costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
    sinthq=dsqrt(1.d0-costhq**2)
    costh1=(mp**2-E1*E2-t/2.d0)/(p1mod*p2mod)

    call bispinor(u1,p1mod,dacos(costh1),0.d0,l1)
    call bispinor(u2,p2mod,pi,0.d0,l2)
    call dirconj(u2,u2bar)

    k=[2.d0*omega,0.d0,0.d0,0.d0]
    epsl=-dble(l)/sqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),xi*dble(l),dcmplx(sinthq)]

    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*epsl)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)


    alfa=af2+af2prime*t
    result=-result*BWM(s1,.775d0,.149d0)*((1.+exp(-xi*pi*alfa))*(s/s0)**alfa)/s

    directPf2pw_MHel=result
end function


complex(8) function directPf2_M_nr(epsl,q,p1,p2,k1,l1,l2,M)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2,M
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,ilo_sk_4d,k2(0:3)
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::qxk1_k2,k(0:3),s,t,s1,k1mod,qmod,costhq,theta,costh,phi,alfa,factor,s1min,s1max
    complex(8)::cilo_sk_4d,u1(1:4),u2(1:4),u2bar(1:4),epsxk1_k2,v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2,k1mk2(0:3),funkcjakulista
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4

    integer::i,l

    p1mod=v3abs(p1)
    p2mod=v3abs(p2)
    k2=q+p1-k1-p2
    call v3angles(p1,theta1,phi1)
    call v3angles(p2,theta2,phi2)
    call bispinor(u1,p1mod,theta1,phi1,l1)
    call bispinor(u2,p2mod,theta2,phi2,l2)
    call dirconj(u2,u2bar)
    k1mk2=cmplx(k1-k2)
    epsxk1_k2=cilo_sk_4d(epsl,k1mk2)
    qxk1_k2=ilo_sk_4d(q,k1-k2)
    k=k1+k2
    k1mod=v3abs(k1)
    qmod=v3abs(q)
    costhq=q(3)/qmod
    call v3angles(k1,theta,phi)
    costh=cos(theta)
    l=-nint(sqrt(2.)*real(epsl(1))/costhq)
    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*epsl)*funkcjakulista(1,M,costh,phi)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)

    s1=ilo_sk_4d(k1+k2,k1+k2)
    t=ilo_sk_4d(p1-p2,p1-p2)
    s=ilo_sk_4d(p1+q,p1+q)
    s1min=4.*mpi**2
    s1max=(sqrt(s)-mp)**2
    factor=(s1min-s1)*(s1max-s1)
    if (.not. free) then
        result=-result*f(s,t,Bf2)*BWM(s1,.775d0,.149d0)*g*bgamropom*RN(s,t,af2,af2prime)/s
    else
        alfa=af2+af2prime*t
        result=-result*factor*((1.+exp(-xi*pi*alfa))*(s/s0)**alfa)/s
    endif

    directPf2_M_nr=result
end function

complex(8) function directPf2_M_nrHel(s,t,mpipi,l,l1,l2,M)
    use parameters
    implicit none
    real(8), intent(in)::s,t,mpipi
    integer,intent(in):: l,l1,l2,M
    complex(8)::epsl(0:3)
    real(8):: theta1,theta2,phi1,phi2,v3abs,p1mod,p2mod,E1,E2,costh1
    real(8),parameter:: e=dsqrt(4.d0*dacos(-1.d0)/137.d0),s0=1.0
    real(8),parameter::A=6.0d0 !Value from Ballam et al. Phys.Rev. D5, 1972, page 574
    real(8),parameter::BP=3.6d0,aP=1.08d0,aPprime=0.2d0 ! pomeron parameters
    real(8),parameter::Bf2=0.55d0, af2=0.5d0,af2prime=0.9d0 !f2 parameters from Vincent paper on rho photoproduction
    real(8)::s1,k1mod,qmod,costhq,sinthq,alfa,factor,s1min,s1max,k(0:3),omega
    complex(8)::u1(1:4),u2(1:4),u2bar(1:4),v(0:3), &
    vgamma(1:4,1:4),vgammau1(1:4),result,u1baru2
!     couplings
    real(8),parameter:: g=5.96d0,bgamropom=2.506d0,mrho=0.77526d0,mu0=10000.d0,t0f2=0.1d-4

    integer::i

    s1=mpipi**2
    qmod=(s1-t)/(2.d0*dsqrt(s1))
    E1=(s-mp**2+t)/(2.d0*dsqrt(s1))
    E2=(s-mp**2-s1)/(2.d0*dsqrt(s1))
    p1mod=dsqrt(E1**2-mp**2)
    p2mod=dsqrt(E2**2-mp**2)
    omega=dsqrt(s1)/2.d0
    k1mod=dsqrt(omega**2-mpi**2)
    costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
    sinthq=dsqrt(1.d0-costhq**2)
    costh1=(mp**2-E1*E2-t/2.d0)/(p1mod*p2mod)

    call bispinor(u1,p1mod,dacos(costh1),0.d0,l1)
    call bispinor(u2,p2mod,pi,0.d0,l2)
    call dirconj(u2,u2bar)

    k=[2.d0*omega,0.d0,0.d0,0.d0]
    epsl=-dble(l)/sqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),xi*dble(l),dcmplx(sinthq)]

    v=-2.*k1mod*sqrt(4.*pi/3.)*(cmplx(k)*b(costhq,l,M)-qmod*rho(costhq,M)*epsl)

    call v4xgam(v,vgamma)
    vgammau1=matmul(vgamma,u1)
    result=u1baru2(u2bar,vgammau1)

    s1min=4.*mpi**2
    s1max=(sqrt(s)-mp)**2
    factor=(s1min-s1)*(s1max-s1)

    alfa=af2+af2prime*t
    result=-result*factor*((1.+exp(-xi*pi*alfa))*(s/s0)**alfa)/s

    directPf2_M_nrHel=result
end function

complex(8) function directPf2pl_nr(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPf2pl_nr=directPf2_M_nr(epsl,q,p1,p2,k1,l1,l2,+1)
end

complex(8) function directPf20_nr(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPf20_nr=directPf2_M_nr(epsl,q,p1,p2,k1,l1,l2,0)
end

complex(8) function directPf2mi_nr(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPf2mi_nr=directPf2_M_nr(epsl,q,p1,p2,k1,l1,l2,-1)
end


complex(8) function directPPompl(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPPompl=directPPom_M(epsl,q,p1,p2,k1,l1,l2,+1)
end

complex(8) function directPPom0(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPPom0=directPPom_M(epsl,q,p1,p2,k1,l1,l2,0)
end

complex(8) function directPpommi(epsl,q,p1,p2,k1,l1,l2)
    use parameters
    implicit none
    complex(8),intent(in)::epsl(0:3)
    real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
    integer,intent(in):: l1,l2
    directPpommi=directPPom_M(epsl,q,p1,p2,k1,l1,l2,-1)
end

      complex(8) function directPf2pl(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      directPf2pl=directPf2_M(eps,q,p1,p2,k1,l1,l2,+1)
      end function

      complex(8) function directPf20(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      directPf20=directPf2_M(eps,q,p1,p2,k1,l1,l2,0)
      end function

      complex(8) function directPf2mi(eps,q,p1,p2,k1,l1,l2)
      implicit none
      complex(8),intent(in)::eps(0:3)
      real(8),intent(in) ::q(0:3),p1(0:3),p2(0:3),k1(0:3)
      integer,intent(in):: l1,l2
      directPf2mi=directPf2_M(eps,q,p1,p2,k1,l1,l2,-1)
      end function

end module pwave

module pwavetest

contains
subroutine pomtest(s,t,mpipi,costh,phi,l,l1,l2)
use pwave
use frames
implicit none
real(8), intent(in):: s,t,costh,phi,mpipi
integer, intent(in):: l,l1,l2
real(8), dimension(0:3):: q,p1,k1,p2
complex(8):: eps(0:3),funkcjakulista,orygin,suma

call Helicity4Vectors(s,t,mpipi**2,costh,phi,l,eps,q,p1,k1,p2)
orygin=directPPom(eps,q,p1,p2,k1,l1,l2)
suma=directPPom_M(eps,q,p1,p2,k1,l1,l2,+1)*funkcjakulista(1,1,costh,phi)+&
directPPom_M(eps,q,p1,p2,k1,l1,l2,0)*funkcjakulista(1,0,costh,phi)+&
directPPom_M(eps,q,p1,p2,k1,l1,l2,-1)*funkcjakulista(1,-1,costh,phi)
print*,orygin
print*,suma
print*,orygin/suma
end subroutine

subroutine f2test(s,t,mpipi,costh,phi,l,l1,l2)
use pwave
use frames
implicit none
real(8), intent(in):: s,t,costh,phi,mpipi
integer, intent(in):: l,l1,l2
real(8), dimension(0:3):: q,p1,k1,p2
complex(8):: eps(0:3),orygin,suma,funkcjakulista

call Helicity4Vectors(s,t,mpipi**2,costh,phi,l,eps,q,p1,k1,p2)

print*,directPf2pl(eps,q,p1,p2,k1,l1,l2)
print*,directPf2pw_M(eps,q,p1,p2,k1,l1,l2,+1)*funkcjakulista(1,1,costh,phi)

print*,directPf20(eps,q,p1,p2,k1,l1,l2)
print*,directPf2pw_M(eps,q,p1,p2,k1,l1,l2,0)*funkcjakulista(1,0,costh,phi)

print*,directPf2mi(eps,q,p1,p2,k1,l1,l2)
print*,directPf2pw_M(eps,q,p1,p2,k1,l1,l2,-1)*funkcjakulista(1,-1,costh,phi)

end subroutine

end module
