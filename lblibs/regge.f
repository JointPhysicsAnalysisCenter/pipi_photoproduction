c 2019.09.09
C------------------------------------------------------------
	complex*16 function Regge(xmv,s,t)
	
	implicit real*8 (a-h,o-z)
	double complex xr,xi	
	
	xr = (1.d0,0.d0)
	xi = (0.d0,1.d0)
	
	pi = dacos(-1.d0)
	
	a0 = 1.d0 
	ap = 0.9d0
	
	alpha = a0 + ap*(t - xmv**2)
	arg = 1.d0 - alpha
	
	Regge = -0.5d0*(1.d0 - xr*dcos(pi*alpha)  
     >                      + xi*dsin(pi*alpha) )*
     >          Gamaf(arg)*( (ap*s)**alpha )/DEXP(A0*DLOG(S))   
      
      	return
      	end
C--------------------------------------------------------
      	
	double precision function Gamaf(z)
	
	implicit real*8 (a-h,o-z)
	real*8 c(26)
	integer i
	
c	Note: only good for |z| < 3.5 

      c(1) =  1.d0
      c(2) =  0.5772156649015329
      c(3) = -0.6558780715202538
      c(4) = -0.0420026350340952
      c(5) =  0.1665386113822915
      c(6) = -0.0421977345555443
  	  c(7) = -0.0096219715278770
  	  c(8) =  0.0072189432466630
  	  c(9) = -0.0011651675918591
  	  c(10)= -0.0002152416741149
  	  c(11)=  0.0001280502823882
  	  c(12)= -0.0000201348547807
  	  c(13)= -0.0000012504934821
  	  c(14)=  0.0000011330272320
  	  c(15)= -0.0000002056338417
  	  c(16)=  0.0000000061160950
  	  c(17)=  0.0000000050020075
  	  c(18)= -0.0000000011812746
  	  c(19)=  0.0000000001043427
  	  c(20)=  0.0000000000077823
  	  c(21)= -0.0000000000036968
  	  c(22)=  0.0000000000005100
  	  c(23)= -0.0000000000000206
  	  c(24)= -0.0000000000000054
  	  c(25)=  0.0000000000000014
  	  c(26)=  0.0000000000000001      

	 Gaminv = 0.d0  	
  	 do i=1,26

	 Gaminv = Gaminv + c(i)*(z**i)  	

	 enddo
	
	 Gamaf = 1.d0/Gaminv
	
	 return
	 end	
 
