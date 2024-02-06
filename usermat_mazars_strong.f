*deck,usermat      USERDISTRIB  parallel                       gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,statev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7, var8)
c******************************************************************

#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), statev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c


      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7, var8

c
c************************ User-defined section*****************************
c
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,EMOD,ENU,EA,EB,EC
      DOUBLE PRECISION Et0,Ec0
      DOUBLE PRECISION aat,bbt,aac,bbc,b
      double precision  Damage
      double precision  eqstrain,pStrain(3),meqstrain,meqstraintemp
      INTEGER      i,j,K1,K2,matsta
      double PRECISION invK(ncomp,2*ncomp),pstress(3)
      double PRECISION yz,pec
      external symeqn
      integer symeqn,op
      double PRECISION tt(2*ncomp,ncomp),SesSrn(ncomp,ncomp)
      double PRECISION tt2(ncomp,ncomp)
      double PRECISION r,Y0,sup1,sup2,Yt,Yc,s,Y,jStrain,iStrain,fm
      
c-----------------------Calculate element stiffness coefficient---------------------------

      PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0)
c-------------EMOD(Elastic modulus)£¬ENU(Poisson's ratio)--------------------------
c-------------Et0,Ec0(Tensile and compressive evolution thresholds)
c-------------aat, bbt, aac, bbc,b represent model parameters
c-------------s represents the substitution rate of Aeolian sand(0-100)-------
      EMOD=prop(1)
      ENU=prop(2)
      Et0=prop(3)
      Ec0=prop(4)
      s=prop(5)
      aac=prop(6)
      aat=prop(7)
      bbc=prop(8)
      bbt=prop(9)
      

      
c     

C------------------------------------------------------------------
C--------------------------Updated strain--------------------------------
C------------------------------------------------------------------
      do i=1,ncomp
          Strain(i)=Strain(i)+dStrain(i)
      end do
      

c     Call 'getDamage' to calculate the damage factor
c
	Damage=statev(1)

      call getr(ncomp,Strain,r,pStrain,stress,pstress,Damage)
      Y0=r*Et0+(1-r)*Ec0
      iStrain=pStrain(1)+pStrain(2)+pStrain(3)
      jStrain=(pStrain(1)-pStrain(2))**2
      jStrain=jStrain+(pStrain(2)-pStrain(3))**2
      jStrain=jStrain+(pStrain(1)-pStrain(3))**2
      jStrain=0.5*jStrain
      sup1=iStrain/(2*(1-2*ENU))+sqrt(jStrain)/(2*(1+ENU))
      sup2=iStrain/(5*(1-2*ENU))+6*sqrt(jStrain)/(5*(1+ENU))
      Yt=max(Et0,sup1)
      Yc=max(Ec0,sup2)
      Y=r*Yt+(1-r)*Yc
      
      
      if(Y.gt.Y0) then
      call getDamage(
     &           Damage,aat,bbt,aac,bbc,r,s,Y,Y0,fm)
      end if
      
      
      
      if (Damage.lt.statev(1)) then
          Damage=statev(1)
      end if
      

      
      if(Damage.GE.1) then
          Damage=0.999999
      end if
      
      do i=1,ncomp
          do j=1,2*ncomp
              invK(i,j)=0.0
              tt(j,i)=0.0
              tt2(j,i)=0.0
          end do
      end do
      
      
      EC=(ONE+ENU)/(ONE-Damage)/EMOD
      EA=-ENU/(ONE-Damage)/EMOD
      EB=(ONE+ENU)/(ONE-Damage)/EMOD+EA
      do i=nShear+1,ncomp
          invK(i,i)=EC
      end do
      
      do i=1,nDirect
          do j=1,nDirect
              invK(i,j)=EA
          end do
      end do
      
      do i=1,3
          invK(i,i)=invK(i,i)+EB
      end do
      
      do i=1,ncomp
          do j=1,2*ncomp
              tt2(j,i)=invK(i,j)
          end do
      end do
      
      op=symeqn(invK,ncomp,ncomp,-ncomp,0)
      
      do i=1,ncomp
          do j=1,2*ncomp
              tt(j,i)=invK(i,j)
          end do
      end do


c     
C     Consistent tangent operator matrix
C
      DO K1=1, ncomp
	     DO K2=1, ncomp
	      dsdePl(K1,K2)=invK(K1,ncomp+K2)
	     END DO
      END DO
      

C------------------------------------------------------------------
C--------------------------Updata Stress--------------------------------
C------------------------------------------------------------------

      do i=1,ncomp
          stress(i)=ZERO
      end do
      
      DO K1=1, ncomp
	     DO K2=1, ncomp
	          stress(K2)=stress(K2)+dsdePl(K2,K1)*Strain(K1)
	     END DO
      END DO
          
      
      

      
C------------------------------------------------------------------
C------------------------Update state variables - damage factor and other parameters----------------------
C------------------------------------------------------------------
      statev(1)=Damage
      statev(2)=meqstrain
      statev(3)=r
      statev(4)=Y
      statev(5)=r
      RETURN
      END

c------------------------------------------------------------------
c------------------Get damage factor D subroutine 'getDamage'--------------------
c------------------------------------------------------------------ 
      SUBROUTINE getDamage(
     &           Damage,aat,bbt,aac,bbc,r,s,Y,Y0,fm)
#include "impcom.inc"
      INTEGER ncomp
      double precision  Damage
      integer I,J
      DOUBLE PRECISION aat,bbt,aac,bbc,r,Y,Y0
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,EMOD,ENU,EA,EB,EC
      PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0)
      character*2 ch1,ch2,ch3,ch4
      DOUBLE PRECISION A,B,pi,xz,s,fs2,fs,fm,fps,rstrain,A1
      real p

      
      p=1
      pi=4*atan(p)
c------------Calculate A------------   
      A=aat*(2*r**2-r)+aac*(2*r**2-3*r+1)
      fs2=-4*10**(-5)*s**2+0.005*s+1
      A=A*fs2

c------------Calculate B-------------   
      B=r**(r**2-2*r+2)*bbt+(1-r**(r**2-2*r+2))*bbc
      
c------------Calculate correction coefficient-------------  
			rstrain=28.98/(1+exp(0.1569*s-4.805))     
			rstrain=28.98-rstrain 
      fps=(100-rstrain)/100*(1/B)
      A1=0.000125*s**2+0.0045*s+0.19
      fm=2*A1/pi*atan(5*B*(Y-1.5*fps))
      fm=fm-2*A1/pi*atan(5*B*(Y0-1.5*fps))+1
      

c-------------Calculate damage factor---------
      Damage=(1-A)*Y0/Y+A*exp(-B*(Y-Y0)*fm)
      Damage=1-Damage
      RETURN
      END
      
c-------------------------------------------------------
c----------Calculate triaxial factor r and principal strain-----------
c-------------------------------------------------------
      SUBROUTINE getr(ncomp,Strain,r,pStrain,stress,pstress,Damage)
      integer ncomp
      double PRECISION Strain(ncomp),eqstrain,pStrain(3),r
      external prinst
      double PRECISION svar(11),stress(ncomp),Damage,pstress(3)
      double PRECISION absstress,kstress
      
      
      do i=1, 11
         svar(i) =0.0
      end do
      do i=1,ncomp
          svar(i)=Strain(i)
      end do
      call prinst(svar)
      pStrain(1)=svar(7)
      pStrain(2)=svar(8)
      pStrain(3)=svar(9)
      
      do i=1, 11
         svar(i) =0.0
      end do
      do i=1,ncomp
          svar(i)=stress(i)/(1-Damage)
      end do
      call prinst(svar)
      pstress(1)=svar(7)
      pstress(2)=svar(8)
      pstress(3)=svar(9)
      
      absstress=0.0
      kstress=0.0
      
      do i=1,3
          absstress=absstress+abs(pstress(i))
          kstress=kstress+(abs(pstress(i))+pstress(i))/2
	end do
	
      r=kstress/absstress
      if (r.lt.0) then
          r=0.0
      end if
	
	      
      RETURN
      END



