module cuba_wrappers


contains

subroutine cuhre_int(ndim, ncomp, lower,upper,integral,neval,userdata) 
implicit none
! real(8), intent(in):: a
integer:: ndim, ncomp, nvec, last, mineval, maxeval
real(8):: epsrel, epsabs, userdata
! parameter (ndim = 2)
! parameter (ncomp = 1)
parameter (nvec = 1)
parameter (epsrel = 1D-3)
parameter (epsabs = 1D-9)
parameter (last = 0)
parameter (mineval = 0)
parameter (maxeval = 50000)

integer:: spin
character*(*) statefile
parameter (statefile = "")
parameter (spin = -1)

integer:: key
parameter (key = 0)

! c integer ScaledIntegrand
! c external ScaledIntegrand


real(8), intent(out):: integral(ncomp)
real(8):: error(ncomp), prob(ncomp)
real(8), intent(in):: lower(ndim),upper(ndim)
integer:: verbose=0, nregions, neval, fail

integer:: maxdim,i
parameter (maxdim = 16)


real(8):: cupper(maxdim)
common /ubound/ cupper

real(8):: clower(maxdim)
common /lbound/ clower


do i=1,ndim
    clower(i)=lower(i)
    cupper(i)=upper(i)
enddo

call cuhre(ndim, ncomp, ScaledIntegrand, userdata, nvec, epsrel, epsabs, verbose + last, &
mineval, maxeval, key,statefile, spin, nregions, neval, fail, integral, error, prob)

end subroutine

integer function ScaledIntegrand(ndim, x, ncomp, result,userdata)
implicit none
integer:: ndim, ncomp
real(8):: x(ndim), result(ncomp), userdata

integer:: maxdim
parameter (maxdim = 16)

! integer:: gaussIntegrand
! external gaussIntegrand

real(8):: upper(maxdim)
common /ubound/ upper

real(8):: lower(maxdim)
common /lbound/ lower

integer:: dim, comp
real(8):: range, jacobian, scaledx(maxdim)
        
        
jacobian = 1
do dim = 1, ndim
    range = upper(dim) - lower(dim)
    jacobian = jacobian*range
    scaledx(dim) = lower(dim) + x(dim)*range
!     print*,scaledx(dim)
enddo

    ScaledIntegrand = gaussIntegrand(ndim, scaledx, ncomp, result,userdata)

do comp = 1, ncomp
    result(comp) = result(comp)*jacobian
enddo
        
end function

integer function gaussIntegrand(ndim, x, ncomp, f,userdata)
integer:: ndim, ncomp
real(8):: x(ndim), f(ncomp),userdata

real(8):: x1,x2,xkwexp,xexp

x1=x(1)
x2=x(2)
a=userdata

f(1)=xkwexp(x1,x2,a)
f(2)=xexp(x1,x2,a)
print*,"1:",x1,x2,a,xkwexp(x1,x2,a)
print*,"2:",x1,x2,a,xexp(x1,x2,a)
! if(f(1)>0 .or. f(2)>0) print*,"f",f(1),f(2)

gaussIntegrand=0
end function

end module cuba_wrappers





