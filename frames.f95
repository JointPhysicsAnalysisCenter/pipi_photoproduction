module frames
implicit none
integer, parameter:: H=1, GJ=2, GJpipp=3, GJpimp=4
contains
      subroutine Helicity4Vectors(s,t,s1,costh,phi,l,eps,q,p1,k1,p2)
      implicit none
      integer,intent(in)::l
      real(8),intent(in)::s,t,s1,costh,phi
      real(8),intent(out)::q(0:3),p1(0:3),k1(0:3),p2(0:3)
      complex(8),intent(out)::eps(0:3)
      real(8) :: qmod,Eq,costhq,sinthq,costh1,sinth1,p1mod,p2mod,E1,E2, &
      omega,kmod,sinth
      real(8):: pm,pim,rom,omegm
      complex(8),parameter::ci=dcmplx(0.d0,1.d0)
      common/masses/pm,pim,rom,omegm
      qmod=(s1-t)/(2.d0*dsqrt(s1))
      Eq=qmod
      E1=(s-pm**2+t)/(2.d0*dsqrt(s1))
      E2=(s-pm**2-s1)/(2.d0*dsqrt(s1))
      p1mod=dsqrt(E1**2-pm**2)
      p2mod=dsqrt(E2**2-pm**2)
      omega=dsqrt(s1)/2.d0
      kmod=dsqrt(omega**2-pim**2)
      costhq=(E1**2-E2**2-qmod**2)/(2.d0*qmod*p2mod)
      sinthq=dsqrt(1.d0-costhq**2)
      costh1=(pm**2-E1*E2-t/2.d0)/(p1mod*p2mod)
      sinth1=dsqrt(1.d0-costh1**2)
      eps=-dble(l)/dsqrt(2.d0)*[dcmplx(0.d0),dcmplx(costhq),ci*dble(l),dcmplx(sinthq)]
      q=[Eq,-qmod*sinthq,0.d0,qmod*costhq]
      p1=[E1,p1mod*sinth1,0.d0,p1mod*costh1]
      sinth=dsqrt(1.d0-costh**2)
      k1=[omega,kmod*sinth*dcos(phi),kmod*sinth*dsin(phi),kmod*costh]
      p2=[E2,0.d0,0.d0,-p2mod]
      end subroutine
      
      subroutine GJ4Vectors(s,t,s1,costh,phi,l,eps,q,p1,k1,p2)
      implicit none
      integer, intent(in)::l
      real(8),intent(in)::s,t,s1,costh,phi
      real(8),intent(out)::q(0:3),p1(0:3),k1(0:3),p2(0:3)
      complex(8),intent(out)::eps(0:3)
      real(8)::qmod,Eq,E1,p1mod,E2,p2mod,omega,kmod,costh2,sinth2,xlbd, &
      sinth,k2(0:3)
      real(8):: pm,pim,rom,omegm
      complex(8),parameter::ci=dcmplx(0.d0,1.d0)
      common/masses/pm,pim,rom,omegm
      qmod=(s1-t)/(2.d0*dsqrt(s1))
      Eq=qmod
      E1=(s-pm**2+t)/(2.d0*dsqrt(s1))
      p1mod=dsqrt(E1**2-pm**2)
      E2=(s-s1-pm**2)/(2.d0*dsqrt(s1))
      p2mod=dsqrt(E2**2-pm**2)
      omega=dsqrt(s1)/2.d0
      kmod=dsqrt(omega**2-pim**2)
      costh2=(2.d0*s1*(s1-s-t+pm**2)+(s1-t)*(s-s1-pm**2))/((s1-t)*dsqrt(xlbd(s,s1,pm**2)))
      sinth2=dsqrt(1.d0-costh2**2)
      sinth=dsqrt(1.d0-costh**2)
      q=[Eq,0.d0,0.d0,qmod]
      p2=[E2,p2mod*sinth2,0.d0,p2mod*costh2]
      k1=[omega,kmod*sinth*dcos(phi),kmod*sinth*dsin(phi),kmod*costh]
      k2=[omega,-kmod*sinth*dcos(phi),-kmod*sinth*dsin(phi),-kmod*costh]
      p1=k1+k2+p2-q
      eps=-dble(l)/dsqrt(2.d0)*[dcmplx(0.d0),dcmplx(1.d0),ci*dble(l),dcmplx(0.d0)]
      end subroutine
      
      subroutine GJpipp4Vectors(s,t2,s1,costh,phi,l,eps,q,p1,k1,p2)
      use parameters
      implicit none
      integer, intent(in)::l
      real(8),intent(in)::s,t2,s1,costh,phi
      real(8),intent(out)::q(0:3),p1(0:3),k1(0:3),p2(0:3)
      complex(8),intent(out)::eps(0:3)
      real(8)::qmod,Eq,E1,p1mod,E2,pmod,omega1,omega2,k2mod,costh2,sinth2,xlbd, &
      sinth,k2(0:3),costhq,sinthq
      real(8):: pm,pim,rom,omegm
      complex(8),parameter::ci=dcmplx(0.d0,1.d0)
      common/masses/pm,pim,rom,omegm
      qmod=(s+t2-pm**2-pim**2)/(2.d0*dsqrt(s1))
      Eq=qmod
      E1=(s1+pm**2-t2)/(2.d0*dsqrt(s1))
      p1mod=dsqrt(E1**2-pm**2)
      E2=(s1-pim**2+pm**2)/(2.d0*dsqrt(s1))
      omega1=(s1-pm**2+pim**2)/(2.d0*sqrt(s1))
      pmod=dsqrt(E2**2-pm**2)
      omega2=(s-pim**2-s1)/(2.d0*sqrt(s1))
      k2mod=dsqrt(xLbd(s,pim**2,s1)/(4.d0*s1))
      costh2=(-2.d0*s1*(s-s1+t2)+(s1-t2+pm**2)*(s-s1-pim**2))/dsqrt(xLbd(s1,pm**2,t2)*xlbd(s,s1,pim**2))
      costhq=(pm**2+2.d0*qmod*E1-s)/(2.d0*qmod*p1mod)
      sinth2=dsqrt(1.d0-costh2**2)
      sinth=dsqrt(1.d0-costh**2)
      sinthq=dsqrt(1.d0-costhq**2)
      q=[Eq,qmod*sinthq,0.d0,qmod*costhq]
      p2=[E2,pmod*sinth*cos(pi+phi),pmod*sinth*sin(pi+phi),-pmod*costh]
      k1=[omega1,pmod*sinth*dcos(phi),pmod*sinth*dsin(phi),pmod*costh]
      k2=[omega2,k2mod*sinth2,0.d0,k2mod*costh2]
      p1=k1+k2+p2-q
      eps=ci/dsqrt(2.d0)*[dcmplx(0.d0),-dcmplx(l)*costhq,dcmplx(0.d0),dcmplx(l)*sinthq]
      end subroutine
      
end module frames
