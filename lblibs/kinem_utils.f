      real*8 function xLbd(x1,x2,x3)
      implicit none
      real*8 x1,x2,x3
      xLbd=x1**2+x2**2+x3**2-2.d0*(x1*x2+x1*x3+x2*x3)
      end
      
      subroutine complexify(re4v,cmplx4v)
      implicit none
      real*8, intent(in):: re4v(0:3)
      complex*16,intent(out):: cmplx4v(0:3)
      integer i
      do i=0,3
      cmplx4v(i)=dcmplx(re4v(i))
      enddo
      end
      
      subroutine v4plusv4(v1,v2,v1plv2)
      implicit none
      real*8, intent(in):: v1(0:3),v2(0:3)
      real*8, intent(out):: v1plv2(0:3)
      integer i
      do i=0,3
        v1plv2(i)=v1(i)+v2(i)
      enddo
      end
      
      subroutine cv4plusv4(v1,v2,v1plv2)
      implicit none
      complex*16, intent(in):: v1(0:3),v2(0:3)
      complex*16, intent(out):: v1plv2(0:3)
      integer i
      do i=0,3
        v1plv2(i)=v1(i)+v2(i)
      enddo
      end

      real*8 function ilo_sk_3d(x1,x2)
      implicit none
      real*8, intent(in):: x1(1:3),x2(1:3)
      integer i
      ilo_sk_3d=0.d0
      do i=1,3
      ilo_sk_3d=ilo_sk_3d+x1(i)*x2(i)
      enddo
      end

      real(8) function ilo_sk_4d(x1,x2)
      implicit none
      real(8), intent(in):: x1(0:3),x2(0:3)
      integer i
      ilo_sk_4d=x1(0)*x2(0)
      do i=1,3
      ilo_sk_4d=ilo_sk_4d-x1(i)*x2(i)
      enddo
      end
      
      complex(8) function cilo_sk_4d(x1,x2)
      implicit none
      complex(8), intent(in):: x1(0:3),x2(0:3)
      integer i
      cilo_sk_4d=x1(0)*x2(0)
      do i=1,3
      cilo_sk_4d=cilo_sk_4d-x1(i)*x2(i)
      enddo
      end
      
      double precision function g(mu,nu)
      implicit none
      integer :: mu,nu
      if(mu.eq.nu) then
        if(mu.eq.0) then
        g=1.d0
        else
        g=-1.d0
        endif
      else
      g=0.d0
      endif  
      end
      
      function v3abs(v4) result(res)
      implicit none
      real(8):: v4(0:3),res
      integer:: i
      res=0.d0
      do i=1,3
      res=res+v4(i)**2
      enddo
      res=dsqrt(res)
      end function v3abs
      
      subroutine v3angles(v4,theta,phi)
      implicit none
      real(8), intent(in) :: v4(0:3)
      real(8), intent(out):: theta,phi
      real(8)::v3abs,vabs
      vabs=v3abs(v4)
      theta=dacos(v4(3)/vabs)
      phi=datan2(v4(2),v4(1))
      end subroutine v3angles
      
      real(8) function kibblefun(s,t,m1sq,m2sq,m3sq,m4sq) result(res)
      implicit none
      real(8), intent(in):: s,t,m1sq,m2sq,m3sq,m4sq
      real(8):: u
      u=m1sq+m2sq+m3sq+m4sq-s-t
      
      res=s*t*u-s*(m1sq-m3sq)*(m2sq-m4sq)-t*(m1sq-m2sq)*(m3sq-m4sq)-
     &(m1sq*m4sq-m3sq*m2sq)*(m1sq+m4sq-m3sq-m2sq)
      
      end function kibblefun
      
      real(8) function tmin(s,m1,m2,m3,m4)
      implicit none
      real(8), intent(in):: s,m1,m2,m3,m4
      real(8):: xlbd
        tmin = m1**2+m3**2-((s-m2**2+m1**2)*(s-m4**2+m3**2)
     &  +sqrt(xlbd(s,m1**2,m2**2)*xlbd(s,m3**2,m4**2)))/(2*s)
      end function

      real(8) function tmax(s,m1,m2,m3,m4)
      implicit none
      real(8), intent(in):: s,m1,m2,m3,m4
      real(8):: xlbd
        tmax = m1**2+m3**2-((s-m2**2+m1**2)*(s-m4**2+m3**2)
     &  -sqrt(xlbd(s,m1**2,m2**2)*xlbd(s,m3**2,m4**2)))/(2*s)
      end function


