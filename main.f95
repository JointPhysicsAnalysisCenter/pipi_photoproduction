program main
use piNamps
use swave
use dwave
use pwave
use parameters
use frames
implicit none
complex (8), dimension (0:4,0:9,0:249) :: PW
common/pwcommon/PW
real(8), parameter:: egam=3.7d0,t=-0.55,s1=0.77d0**2,costh=.5d0,phi=.2d0
integer, parameter:: l=+1,l1=+1,l2=-1 !Note l1,l2=+/-1 correspond to nucleon's +/- 1/2 helicities
complex(8):: eps(0:3),DeckPL,DeckMI
real(8),dimension(0:3):: q,p1,k1,p2
real(8):: s=mp**2+2.d0*egam*mp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! s - p-gamma energy
! t - target-recoil 4-momentum transfer squared
! s1 - pipi energy
! costh,phi - cosine of theta and phi in the helicity frame
!
! Call GJ4Vectors(s,t,s1,costh,phi,l,eps,q,p1,k1,p2) to obtain 4-vectors in the Gottfried-Jackson frame
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call initCommons()
call SAID_PiN_PW(PW)
call Helicity4Vectors(s,t,s1,costh,phi,l,eps,q,p1,k1,p2)

print*,"S-wave f0(980)"
print*, directS980(eps,q,p1,p2,k1,l1,l2)
print*, "--------------------------------"
print*,"P-wave rho(770) /f2-exchange/"
print*, directPf2(eps,q,p1,p2,k1,l1,l2)
print*, "--------------------------------"
print*,"D-wave f2(1270)"
print*, directD(eps,q,p1,p2,k1,l1,l2)
print*, "--------------------------------"
print*, "Deck (for Deck we additionally need to pass pi-p partial waves)"
print*, DeckPL(eps,q,p1,p2,k1,l1,l2,PW)
print*, DeckMI(eps,q,p1,p2,k1,l1,l2,PW)
print*, "--------------------------------"
end program main

