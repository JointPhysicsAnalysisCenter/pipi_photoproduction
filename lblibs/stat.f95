module stat
implicit none
contains

real(8) function uniform(a,b) result(val)
implicit none
real(8), intent(in):: a,b
real(8):: x
call random_number(x)
val=a+(b-a)*x
end

real(8) function expdist(lbd) result(val)
implicit none
real(8), intent(in):: lbd
real(8):: x,R
    x=uniform(0.d0,20.d0)
    R=1.d0-exp(-lbd*x)
    val= -1.d0/lbd*log(R)
end

real(8) function circdist(mu,tau) result(val)
implicit none
real(8), intent(in):: mu,tau
real(8), parameter:: pi=acos(-1.d0)
real(8):: u,e,alfa
    u=uniform(-pi/2.d0,pi/2.d0)
    e=expdist(.8d0)
    alfa=2.d0
    val=mod(tau*sin(alfa*u)/cos(u)**(1.d0/alfa)*(cos((1.d0-alfa)*u)/e)**((1.d0-alfa)/alfa)+mu,2.d0*pi)

end

subroutine AvgAndStdev(data,Mean,StdDev)
        implicit none
        real(8), intent(in):: data(:)
        real(8), intent(out):: Mean,StdDev
        real(8):: Variance
        integer:: i,SimCount
        SimCount=size(data)
        Mean = 0.0                           ! compute mean
        do i = 1, SimCount
            Mean = Mean + Data(i)
        end do
        Mean = Mean / SimCount

        Variance = 0.0                       ! compute variance
        DO i = 1, SimCOunt
            Variance = Variance + (Data(i) - Mean)**2
        END DO
        Variance = Variance / (SimCount - 1)
        StdDev   = SQRT(Variance)
end subroutine

subroutine SampleDist4D(fun,low,high,i,sample,g)
    implicit none
    real(8):: fun
    real(8), intent(in):: low(0:3),high(0:3)
    integer, intent(in):: i
    real(8), dimension(0:i-1,0:3),intent(out):: sample
    logical, optional:: g
    integer:: lsize, hsize,allocstatus
    real(8), dimension(0:3):: stepsizes
    
    integer, parameter:: meshdensity=100
    integer:: j,k,l,m
    real(8):: maxfun=0.,x0,x1,x2,x3,f
    logical:: runend=.false.,condition=.true.
    
    do k=0,3
        stepsizes(k)=(high(k)-low(k))/meshdensity
    enddo
    
    do j=0,meshdensity-1
        x0=low(0)+j*stepsizes(0)
        do k=0,meshdensity-1
            x1=low(1)+j*stepsizes(1)
            do l=0,meshdensity-1
                x2=low(2)+j*stepsizes(2)
                do m=0,meshdensity-1
                    x3=low(3)+j*stepsizes(3)
                    if (present(g)) then
                    if (.not. g(x0,x1,x2,x3)) then
                      cycle
                    endif
                    endif
                    f=fun(x0,x1,x2,x3)
                    if(f>maxfun) then
                        maxfun=f
                    end if
                enddo
            enddo
        enddo
    enddo

    print*,maxfun
    
    j=0
    call random_seed()
    do while(.not. runend)
        x0=uniform(low(0),high(0))
        x1=uniform(low(1),high(1))
        x2=uniform(low(2),high(2))
        x3=uniform(low(3),high(3))
        f=uniform(0.d0,maxfun)
        if(present(g))then
            condition=g(x0,x1,x2,x3)
        endif
        if(f<fun(x0,x1,x2,x3) .and. condition) then
            print*,j
            sample(j,0)=x0
            sample(j,1)=x1
            sample(j,2)=x2
            sample(j,3)=x3
            j=j+1
        endif
        if(j==i) then
            runend=.true.
        endif
    enddo
end subroutine

function urandom_seed(n, stat) result(seed)
        !! Returns a seed array filled with random values from `/dev/urandom`.
        integer, intent(in)            :: n
        integer, intent(out), optional :: stat
        integer                        :: seed(n)
        integer                        :: fu, rc

        open (access='stream', action='read', file='/dev/urandom', &
            form='unformatted', iostat=rc, newunit=fu)
        if (present(stat)) stat = rc
        if (rc == 0) read (fu) seed
        close (fu)
end function urandom_seed
end module
