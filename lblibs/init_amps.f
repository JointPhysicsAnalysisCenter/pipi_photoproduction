      subroutine initCommons()
      call initMasses()
!       call initKinematics()
      call initPropag() !must be called before initConstants()
      call initConstants()
      call initGammas()
      end
      
      subroutine initConstants()
      implicit real*8 (a-h,o-z)
      common/stale/pi,e,gpipiro,GrhoV,GrhoT,GomV,GomT,groompi,gropigam,
     &gompigam,CONV
      common/kinematyka/omega,t,empipi,s
      common/masses/pm,pim,rom,omegm
      common/propagator/iprop
      complex*16 regge
      
      pi=dacos(-1.d0)
      e=0.30282d0
      gpipiro=6.05d0
      GrhoV=2.27d0
      GrhoT=13.85d0
      GomV=11.54d0
      GomT=0.d0
c   groompi has dimension GeV^-1
      groompi=14.d0 !dolna wartość
      gropigam=0.75d0*E/rom
      gompigam=1.82d0*E/omegm
      CONV=389.38d0
!       pm2=pm**2

      if(iprop.eq.1) then
!       S0=PM2+2.d0*PM*OMEGA
!       S0=2.9d0**2*2.d0
      S0=2.5d0
      GrhoV=GrhoV/rom**2/CDABS(REGGE(rom,S0/2.d0,0.D0))
      GrhoT=GrhoT/rom**2/CDABS(REGGE(rom,S0/2.d0,0.D0))
      GomV=GomV/omegm**2/CDABS(REGGE(omegm,S0/2.d0,0.D0))
      endif
     
      end
      
      subroutine initKinematics()
      implicit real*8(a-h,o-z)
      common/kinematyka/omega,t,spipi,s
      common/masses/pm,pim,rom,omegm
      
!       omega=3.3d0
      omega=4.8596d0
      t=-.7d0
      spipi=1.0
      s=pm**2+2.d0*omega*pm
      end
      
      subroutine initMasses()
      implicit real*8(a-h,o-z)
      common/masses/pm,pim,rom,omegm
      
      pm=0.93827203d0
      pim=.13957018d0
      rom=0.77549d0
      omegm=0.78265d0
      end
      
      subroutine initPropag()
      integer iprop
      common/propagator/iprop
      iprop=0
      end
