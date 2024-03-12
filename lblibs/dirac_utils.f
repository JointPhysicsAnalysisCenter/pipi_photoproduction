      subroutine initGammas()
      implicit none
      complex*16 d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),d3(1:4,1:4),
     &d5(1:4,1:4),unit(1:4,1:4),cu,ci,zero
    
      parameter (cu=dcmplx(1.d0,0.d0))
      parameter (ci=dcmplx(0.d0,1.d0))
      parameter (zero=dcmplx(0.d0,0.d0))
      common/gammas/d0,d1,d2,d3,d5,unit
      data d0/16*zero/
      data d1/16*zero/
      data d2/16*zero/
      data d3/16*zero/
      data d5/16*zero/
      data unit/16*zero/
      d0(1,1)=cu
      d0(2,2)=cu
      d0(3,3)=-cu
      d0(4,4)=-cu
      
      d1(1,4)=cu
      d1(2,3)=cu
      d1(3,2)=-cu
      d1(4,1)=-cu
      
      d2(1,4)=-ci
      d2(2,3)=ci
      d2(3,2)=ci
      d2(4,1)=-ci
      
      d3(1,3)=cu
      d3(2,4)=-cu
      d3(3,1)=-cu
      d3(4,2)=cu
      
      d5(1,3)=cu
      d5(2,4)=cu
      d5(3,1)=cu
      d5(4,2)=cu
      
      unit(1,1)=cu
      unit(2,2)=cu
      unit(3,3)=cu
      unit(4,4)=cu
      end
      
      subroutine printd(d)
      complex*16 d(4,4)
      
      do i=1,4
        write(*,12)(d(i,j),j=1,4)
      enddo
12    format("(",F4.1,",",F4.1,")(",F4.1,",",F4.1,")(",F4.1,",",F4.1,
     &")(",F4.1,",",F4.1,")")  
      end
      
      subroutine bispinor(b,pmod,theta,phi,ihel)
c       ihel =+/1 correspond to +/- hslash/2
      implicit none
      complex*16, intent(out):: b(1:4)
      real*8, intent(in):: pmod,theta,phi
      integer, intent(in):: ihel
      complex*16 ci
      real*8 costh,sinth_2,costh_2,Ep
      real*8 pm,pim,rom,omegm
      common/masses/pm,pim,rom,omegm
      
      costh=dcos(theta)
      
      costh_2=dsqrt(1.d0+costh)/dsqrt(2.d0)

      ci=dcmplx(0.d0,1.d0)
      Ep=dsqrt(pmod**2+pm**2)
      sinth_2=dsqrt(1.d0-costh_2**2)
c       print*,sinth_2,costh_2
      if(ihel.eq.1) then
      b(1)=dsqrt(Ep+pm)*costh_2
      b(2)=dsqrt(Ep+pm)*sinth_2*cdexp(ci*phi)
      b(3)=dsqrt(Ep+pm)*pmod/(Ep+pm)*costh_2
      b(4)=dsqrt(Ep+pm)*pmod/(Ep+pm)*sinth_2*cdexp(ci*phi)
      return
      endif
      if(ihel.eq.-1) then
      b(1)=-dsqrt(Ep+pm)*sinth_2*cdexp(-ci*phi)
      b(2)=dsqrt(Ep+pm)*costh_2
      b(3)=dsqrt(Ep+pm)*pmod/(Ep+pm)*sinth_2*cdexp(-ci*phi)
      b(4)=-dsqrt(Ep+pm)*pmod/(Ep+pm)*costh_2
      return
      endif
      end
      
      subroutine dirconj(bin,bout)
      implicit none
      complex*16, intent(in):: bin(1:4)
      complex*16, intent(out):: bout(1:4)
      integer i,j
      complex*16 d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),
     &d3(1:4,1:4),d5(1:4,1:4),unit(1:4,1:4),zero
      parameter (zero=dcmplx(0.d0,0.d0))
      common/gammas/d0,d1,d2,d3,d5,unit
c       print*,"d0=",d0
      do i=1,4
        bout(i)=zero
      enddo
      do j=1,4
        do i=1,4
        bout(j)=bout(j)+dconjg(bin(i))*d0(i,j)
        enddo
      enddo
      end
      
      subroutine v4xgam(v4,v4xg)
      implicit none
      complex*16, intent(in):: v4(0:3)
      complex*16, intent(out):: v4xg(1:4,1:4)
      complex*16 d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),d3(1:4,1:4),
     &d5(1:4,1:4),unit(1:4,1:4)
      integer :: i,j
      common/gammas/d0,d1,d2,d3,d5,unit
      call zeromatrix(v4xg)
          do i=1,4
            do j=1,4
              v4xg(i,j)=v4xg(i,j)+v4(0)*d0(i,j)
            enddo
          enddo
          
          do i=1,4
            do j=1,4
              v4xg(i,j)=v4xg(i,j)-v4(1)*d1(i,j)
            enddo
          enddo

          do i=1,4
            do j=1,4
              v4xg(i,j)=v4xg(i,j)-v4(2)*d2(i,j)
            enddo
          enddo

          do i=1,4
            do j=1,4
              v4xg(i,j)=v4xg(i,j)-v4(3)*d3(i,j)
            enddo
          enddo
      end
      
      subroutine scalarxgam(s,gin,gout)
      implicit none
      complex*16, intent(in) :: s,gin(1:4,1:4)
      complex*16, intent(out) :: gout(1:4,1:4)
      integer i,j
      do i=1,4
        do j=1,4
          gout(i,j)=s*gin(i,j)
        enddo
      enddo
      end
      
      subroutine zeromatrix(cm)
      implicit none
      complex*16, intent(out):: cm(1:4,1:4)
      integer i,j
      do i=1,4
        do j=1,4
          cm(i,j)=dcmplx(0.d0)
        enddo
      enddo
      end
      
      complex*16 function u1baru2(u1bar,u2)
      implicit none
      complex*16, intent(in):: u1bar(4),u2(4)
      integer i
      u1baru2=dcmplx(0.d0)
      
      do i=1,4
      u1baru2=u1baru2+u1bar(i)*u2(i)
      enddo
    
      end
      
      complex*16 function trace(gam)
      complex*16, intent(in):: gam(4,4)
      trace=dcmplx(0.d0)
      do i=1,4
      trace=trace+gam(i,i)
      enddo
      end
      
      subroutine hermconj(in,out)
      implicit none
      complex*16, intent(in):: in(1:4,1:4)
      complex*16, intent(out):: out(1:4,1:4)
      integer i,j
      do i=1,4
        do j=1,4
          out(i,j)=dconjg(in(j,i))
        enddo
      enddo
      end
      
      subroutine matmul(in1,in2,out)
      implicit none
      complex*16, intent(in):: in1(1:4,1:4),in2(1:4,1:4)
      complex*16, intent(out):: out(1:4,1:4)
      integer i,j,k
      call zeromatrix(out)
      do i=1,4
      do j=1,4
        do k=1,4
        out(i,j)=out(i,j)+in2(i,k)*in1(k,j)
        enddo
      enddo
      enddo
      end
      
      subroutine matdirconj(in,out)
      implicit none
      complex*16,intent(in):: in(1:4,1:4)
      complex*16,intent(out):: out(1:4,1:4)
      complex*16 hin(1:4,1:4),hind0(1:4,1:4)
      complex*16 d0(1:4,1:4),d1(1:4,1:4),d2(1:4,1:4),d3(1:4,1:4),
     &d5(1:4,1:4),unit(1:4,1:4)
      common/gammas/d0,d1,d2,d3,d5,unit
      call hermconj(in,hin)
      call matmul(d0,hin,hind0)
      call matmul(hind0,d0,out)
      end
      
      
! bispinor with cosinetheta rather than theta as argument
      subroutine bispinor_costh(b,pmod,costh,phi,ihel)
c       ihel =+/1 correspond to +/- hslash/2
      implicit none
      complex*16, intent(out):: b(1:4)
      real*8, intent(in):: pmod,phi
      integer, intent(in):: ihel
      complex*16 ci,costh,sinth_2,costh_2
      real*8 Ep
      real*8 pm,pim,rom,omegm
      common/masses/pm,pim,rom,omegm
      
      costh_2=cdsqrt(1.d0+costh)/dsqrt(2.d0)

      ci=dcmplx(0.d0,1.d0)
      Ep=dsqrt(pmod**2+pm**2)
      sinth_2=cdsqrt(1.d0-costh_2**2)
c       print*,sinth_2,costh_2
      if(ihel.eq.1) then
      b(1)=dsqrt(Ep+pm)*costh_2
      b(2)=dsqrt(Ep+pm)*sinth_2*cdexp(ci*phi)
      b(3)=dsqrt(Ep+pm)*pmod/(Ep+pm)*costh_2
      b(4)=dsqrt(Ep+pm)*pmod/(Ep+pm)*sinth_2*cdexp(ci*phi)
      return
      endif
      if(ihel.eq.-1) then
      b(1)=-dsqrt(Ep+pm)*sinth_2*cdexp(-ci*phi)
      b(2)=dsqrt(Ep+pm)*costh_2
      b(3)=dsqrt(Ep+pm)*pmod/(Ep+pm)*sinth_2*cdexp(-ci*phi)
      b(4)=-dsqrt(Ep+pm)*pmod/(Ep+pm)*costh_2
      return
      endif
      end
