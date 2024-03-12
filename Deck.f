c       real(8) function tmax(s,m1,m2,m3,m4)
c       implicit none
c       real(8), intent(in)::s,m1,m2,m3,m4
c       real(8):: xLbd
c         tmax=m1**2+m3**2-((s-m2**2+m1**2)*(s-m4**2+m3**2)-
c      &  sqrt(xLbd(s,m1**2,m2**2)*xLbd(s,m3**2,m4**2)))/(2.d0*s)
c       end function
      
      complex(8) function DeckPL(eps,q,p1,p2,k1,l1,l2,PW)
      use piNamps
      implicit none
     
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex(8)::eps(0:3),u1(1:4),u2(1:4),u2bar(1:4),A32,B32,Apl,Bpl,v
      real(8)::q(0:3),p1(0:3),p2(0:3),k1(0:3),theta1,theta2,phi1,phi2
      real(8):: v3abs,p1mod,p2mod,ilo_sk_4d,s2,t,
     &tpi,k2(0:3),e,tmx,tmax,s,ff
      integer::i,j,l1,l2
      complex(8):: Qu(0:3),Qslash(1:4,1:4),ck2(0:3),cilo_sk_4d,
     &cp1plp2(0:3),ck1(0:3)
      complex(8) d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),
     &d3(1:4,1:4),d5(1:4,1:4),unit(1:4,1:4),ABC(0:1,0:2)
      real(8):: pm,pim,rom,omegm
      common/masses/pm,pim,rom,omegm
      common/gammas/d0,d1,d2,d3,d5,unit
      
!       e=0.30282d0
      e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
      p1mod=v3abs(p1)
      p2mod=v3abs(p2)
      do i=0,3
        k2(i)=q(i)+p1(i)-k1(i)-p2(i)
      enddo
      call v3angles(p1,theta1,phi1)
      call v3angles(p2,theta2,phi2)
      s2=ilo_sk_4d(p2+k1,p2+k1)
      t=ilo_sk_4d(p1-p2,p1-p2)
      tpi=ilo_sk_4d(q-k2,q-k2)
      s=ilo_sk_4d(q+p1,q+p1)
      tmx=tmax(s,0.d0,pm,pim,sqrt(s2))
c       print*,s2
      call AMP32(s2,t,tpi,PW,A32,B32)
!       call AMP32(s2,t,pim**2,PW,A32,B32)

c       call Full_ABC(s2,t,tpi,PW,ABC)
      Apl=A32
      Bpl=B32
c       Apl=ABC(0,0)-ABC(1,0)
c       Bpl=ABC(0,1)-ABC(1,1)
      call bispinor(u1,p1mod,theta1,phi1,l1)
      call bispinor(u2,p2mod,theta2,phi2,l2)
      call dirconj(u2,u2bar)
      call complexify(k1,Qu)
      call v4xgam(Qu,Qslash)
      call complexify(k2,ck2)
      call complexify(k1,ck1)
      call complexify(p1+p2,cp1plp2)
      
      v=e*(-cilo_sk_4d(eps,ck2)/ilo_sk_4d(q,k2)
     &+cilo_sk_4d(eps,cp1plp2)/ilo_sk_4d(q,p1+p2)
     &)
     
      ff=exp((tpi-tmx)/.9d0**2)
      DeckPL=dcmplx(0.d0)
      do i=1,4
      do j=1,4
      DeckPL=DeckPL+v*u2bar(i)*
     &(Apl*unit(i,j)+Qslash(i,j)*Bpl)*u1(j)*ff
      enddo
      enddo

      end

      complex(8) function DeckMI(eps,q,p1,p2,k1,l1,l2,PW)
      use piNamps
      implicit none
      
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex(8)::eps(0:3),u1(1:4),u2(1:4),u2bar(1:4),A12,B12,A32,B32,
     &Ami,Bmi
      real(8)::q(0:3),p1(0:3),p2(0:3),k1(0:3),k2(0:3),p2plk2(0:3),
     &p1minp2(0:3),tpi,qmink1(0:3),e,tmx,tmax,s,ff
      real(8):: v3abs,p1mod,p2mod,ilo_sk_4d,s2,t,theta1,theta2,phi1,
     &phi2
      integer::i,j,l1,l2
      complex(8):: Qu(0:3),Qslash(1:4,1:4),ck1(0:3),v,cilo_sk_4d,
     &cp1plp2(0:3)
      complex(8) d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),
     &d3(1:4,1:4),d5(1:4,1:4),unit(1:4,1:4),ABC(0:1,0:2)
      real(8):: pm,pim,rom,omegm
      common/masses/pm,pim,rom,omegm
      common/gammas/d0,d1,d2,d3,d5,unit
!       e=0.30282d0
      e=dsqrt(4.d0*dacos(-1.d0)/137.d0)
      p1mod=v3abs(p1)
      p2mod=v3abs(p2)
      call v3angles(p1,theta1,phi1)
      call v3angles(p2,theta2,phi2)
      do i=0,3
        k2(i)=q(i)+p1(i)-k1(i)-p2(i)
      enddo
      call v4plusv4(k2,p2,p2plk2)
      call v4plusv4(p1,-p2,p1minp2)
      call v4plusv4(q,-k1,qmink1)
      s2=ilo_sk_4d(p2+k2,p2+k2)
      t=ilo_sk_4d(p1minp2,p1minp2)
      tpi=ilo_sk_4d(qmink1,qmink1)
      s=ilo_sk_4d(q+p1,q+p1)
      tmx=tmax(s,0.d0,pm,pim,sqrt(s2))
      call AMP12(s2,t,tpi,PW,A12,B12)
      call AMP32(s2,t,tpi,PW,A32,B32)
!       call AMP12(s2,t,pim**2,PW,A12,B12)
!       call AMP32(s2,t,pim**2,PW,A32,B32)
      Ami=(A32+2.d0*A12)/3.d0
      Bmi=(B32+2.d0*B12)/3.d0
c       call Full_ABC(s2,t,tpi,PW,ABC)
c       Ami=ABC(0,0)+ABC(1,0)
c       Bmi=ABC(0,1)+ABC(1,1)
      call bispinor(u1,p1mod,theta1,phi1,l1)
      call bispinor(u2,p2mod,theta2,phi2,l2)
      call dirconj(u2,u2bar)
      
      call complexify(k2,Qu)
      call v4xgam(Qu,Qslash)
            
      call complexify(k1,ck1)
      call complexify(p1+p2,cp1plp2)
      
      v=e*(cilo_sk_4d(eps,ck1)/ilo_sk_4d(q,k1)-cilo_sk_4d(eps,cp1plp2)/
     &ilo_sk_4d(q,p1+p2))
      
      ff=exp((tpi-tmx)/.9d0**2)
      DeckMI=dcmplx(0.d0)
      do i=1,4
      do j=1,4
      DeckMI=DeckMI+v*u2bar(i)*
     &(Ami*unit(i,j)+Qslash(i,j)*Bmi)*u1(j)*ff
      enddo
      enddo
      end
      
      subroutine AMP12(s,t,tpi,PW,A12,B12)
      use piNamps
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A12,B12
      complex (8), intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      
      call Full_ABC(s,t,tpi,PW,ABC)
      A12=ABC(0,0)+2.d0*ABC(1,0)
      B12=ABC(0,1)+2.d0*ABC(1,1)
      end
      
      subroutine AMP32(s,t,tpi,PW,A32,B32)
      use piNamps
      implicit none
      real (8), intent(in):: s,t,tpi
      complex (8), intent(out):: A32,B32
      complex (8),intent(in), dimension (0:4,0:9,0:249) :: PW
      complex (8), dimension (0:1,0:2) :: ABC
      call Full_ABC(s,t,tpi,PW,ABC)
      A32=ABC(0,0)-ABC(1,0)
      B32=ABC(0,1)-ABC(1,1)
      end
      
      complex(8) function Amin_hel(s,t,l1,l2,PW) result(val)
      implicit none
      real(8):: s,t,pmod,theta,costh,sinth,xLbd,k1(0:3),k2(0:3),omega,
     &Q(0:3)
      integer:: i,l1,l2
      complex(8) A12,A32,B12,B32,Amin,Bmin,u1(1:4),u2(1,4),u2bar(1:4),
     &Qslash(1:4,1:4),Amat(1:4,1:4),Bmat(1:4,1:4),cQ(0:3),A(1:4,1:4),
     &utmp(1:4)
      complex (8),intent(in), dimension (0:4,0:9,0:249) :: PW
      real(8):: pm,pim,rom,omegm
      complex(8) d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),
     &d3(1:4,1:4),d5(1:4,1:4),unit(1:4,1:4)
      common/masses/pm,pim,rom,omegm
      common/gammas/d0,d1,d2,d3,d5,unit
      call AMP12(s,t,pim**2,PW,A12,B12)
      call AMP32(s,t,pim**2,PW,A32,B32)
      Amin=(A12-A32)/3.d0
      Bmin=(B12-B32)/3.d0
      pmod=dsqrt(xLbd(s,pm**2,pim**2)/(4.d0*s))
      omega=dsqrt(pmod**2+pim**2)
      costh=1.d0+t/(2.d0*pmod**2)
      sinth=dsqrt(1.d0-costh**2)
      theta=dacos(costh)
      k1(0)=omega
      k1(1)=0.d0
      k1(2)=0.d0
      k1(3)=-pmod
      k2(0)=omega
      k1(1)=-pmod*sinth
      k1(2)=0.d0
      k1(3)=-pmod*costh
      Q=(k1+k2)/2.d0
      call bispinor(u1,pmod,0.d0,0.d0,l1)
      call bispinor(u2,pmod,theta,0.d0,l2)
      call dirconj(u2,u2bar) 
      Amat=unit
      Amat=Amat*Amin
      call complexify(Q,cQ)
      call v4xgam(cQ,Qslash)
      Bmat=Qslash*Bmin
      A=Amat+Bmat
      utmp=matmul(A,u1)
      val=dcmplx(0.d0,0.d0)
      do i=1,4
        val=val+u2bar(i)*utmp(i)
      enddo
      end function Amin_hel
      
      real(8) function strength_amin(s,t,PW) result(val)
      implicit none
      complex (8),intent(in), dimension (0:4,0:9,0:249) :: PW
      real(8):: p,costh,s,t,xLbd
      real(8):: pm,pim,rom,omegm
      common/masses/pm,pim,rom,omegm
      complex(8):: Amin_hel
      integer:: i,j
      val=0.d0
      p=sqrt(xLbd(s,pim**2,pm**2)/(4.d0*s))
      costh=1.d0+t/(2.d0*p**2)
      do i=-1,1,2
      do j=-1,1,2
        val=val+abs(Amin_hel(s,t,i,j,PW))**2
      enddo
      enddo
      end function strength_amin
