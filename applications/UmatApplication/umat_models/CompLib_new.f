!     UMAT - Eleni Gerolymatou
      SUBROUTINE CompLib_New(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!

      USE tensors
!      USE nlmodule   !! non local?

      IMPLICIT NONE
!
!--------------------------------------------------------------- 
!     Declarating UMAT vaiables and constants 
      CHARACTER*80 CMNAME
      INTEGER NTENS,NSTATEV,NPROPS,LAYER,NOEL,NPT,KSPT,KSTEP,KINC,
     1 NDI,NSHR
      DOUBLE PRECISION  STRESS(NTENS),STATEV(NSTATEV),SSE,SPD,SCD,
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 MAT(2,2)
! ---------------------------------------------------------------      
      DOUBLE PRECISION RPL, DRPLDT, DTIME,
     1 TEMP, PNEWDT, CELENT, DTEMP 
! ---------------------------------------------------------------      
      DOUBLE PRECISION T(3,3), DEPS(3,3), JAC(3,3,3,3),
     1 evoid, epsP(3,3), pc, pt,
     2 ELMOD(3,3,3,3), Idev(3,3,3,3), TrialSTRESS(3,3)
     3 , TrialTh, uu,param(2),derve(6),dervp(6),dervs(6),dervh(6) 
     4 , TrialGth,Trialp,Trialq,TrialF,A(13,13),B(13),dum,xx(13)
     5 ,theta,TolG, TolF,epsG,epsF
     6 ,normaldev(3,3) ,Depsp(3,3), NormDepsp, Isym(3,3,3,3)      
      integer maxiter, kcounter,indx(13)
!     parameters for the material
      DOUBLE PRECISION E,nu,alpha,beta,Mf,cc,m,rhos,ksis,rhom,ksim,pc0
      DOUBLE PRECISION dep,deq,ps,pm,pss,pmm,denom,rrm
      DOUBLE PRECISION dfds(3,3),dfdep(3,3),aux(3,3),Hplast(3,3,3,3)
      DOUBLE PRECISION DSTRANV(6)
!     preparing the information

	  

	  
      CALL Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)      
! ---------------------------------------------------------------  
!     Defining material parameters 
      E     =  PROPS(1) !young modulus
      nu    =  PROPS(2) !poisson modulus 
!     other material properties: alpha,beta...	
      alpha =  PROPS(3) !alpha variable of the exponential
      beta  =  PROPS(4) !beta variable of the exponential	
      Mf    =  PROPS(5) !Mf
      cc    =  PROPS(6) !c for the Lode angle
      m     =  PROPS(7) !m value for pt
      rhos  =  PROPS(8) !rhos
      ksis  =  PROPS(9) !ksis
      rhom  =  PROPS(10)!rhom 
      ksim  =  PROPS(11)!ksim 
      pc0   =  PROPS(12)!dummy variable to make units fit
! ---------------------------------------------------------------     
!     State variables 
!     Void ratio as state variable (scalar)
      evoid=STATEV(1)
!     Plastic strains as state variable (tensor R3xR3)
      CALL Vectortomatrix(STATEV(2:7),epsP) 
!     Compression Strength
      pc= STATEV(8)+STATEV(9)
      pt= m*STATEV(9)
!     ps
      param(1) = STATEV(8)
!     pm
      param(2) = STATEV(9)
! ---------------------------------------------------------------      
!     Evaluating the elastic stiffness tensor
      CALL ELMOD1(E, nu, ELMOD) 
! ---------------------------------------------------------------    
!     Trial Elastic trial step
!     Trial stress 
      TrialSTRESS=(ELMOD.double.DEPS)+T
!     the invariants need to be evaluated next
!     Deviator fourth unit tensor/ (ds/dsigma)
      CALL Idev1(Idev) 
!     Trial deviator stress  
!      TrialDevSTRESS=Idev.double.TrialSTRESS   
!     Trial yield function F
      call Fyield(TrialSTRESS,STATEV,NSTATEV,NTENS,
     &  PROPS,NPROPS,param,TrialF)
! ---------------------------------------------------------------       
      if (TrialF<0.0d0) then
!        Elastic Steps, Trial step is ok
         T=TrialSTRESS
         JAC=ELMOD
         evoid=evoid+(1.0d0+evoid)*(Tr(Deps))
         STATEV(1)=evoid
!        No need to save state variables 
!        They are already correct
!         write(6,"(E10.5)") DSTRAN
! ---------------------------------------------------------------       
      else ! Plastic corrector step

!         write(6,*) 'plastic', TrialF
! ---------------------------------------------------------------
!     Defining tolerances for G and F and maximum iterations
         TolG=1.0E-4               ! Tolerance for residuals
         TolF=1.0E-9               ! Tolerance for increments
         maxiter=50                ! maximum iterations
!     Initializing variables
!     Elastic strain increment
         xx(1:6) = DSTRAN
         xx(4:6) = xx(4:6)/2.0
		 DSTRANV = xx(1:6)
!     Plastic strain increment
         xx(7:12)= 0.0
!     Plastic multiplier
         xx(13)  = 0.0

         if (abs( STATEV(8))  <1.0d-8) write(*,*) 'ps not set '
         if (abs( STATEV(9))  <1.0d-8) write(*,*) 'pm not set '
!     ps
      param(1) = STATEV(8)
!     pm
      param(2) = STATEV(9)
! ---------------------------------------------------------------
!        Cycle for plastic increment 
         do kcounter=1,maxiter
! ---------------------------------------------------------------
!            calling the yield function
             call Fyield(TrialSTRESS,xx(7:12),NSTATEV,NTENS,
     &            PROPS,NPROPS,param,TrialF)
!            calling the derivative evaluation
             call DerivF(TrialSTRESS,xx,NTENS,
     &            PROPS,NPROPS,param,derve,dervp,dervs)

!            defining the matrices A and B
!            initializing
             A=0.0
             B=0.0
             do i=1,6
                 A(i,i)     = 1.0
                 A(i,i+6)   = 1.0
                 A(i+6,i+6) = 1.0
                 A(i+6,13)  =-dervs(i)
                 A(13,i)    = derve(i)
                 A(13,i+6)  = dervp(i)
                 B(i)  = DSTRANV(i)-xx(i)-xx(i+6)
                 B(i+6)= abs(xx(13))*dervs(i)-xx(i+6) 
             enddo
             B(13)=-TrialF

!            evaluate residuals
             epsG = sqrt(dot3(B,B,13))  

             call ludcmp(A,13,13,indx,dum)
             call lubksb(A,13,13,indx,B)
!            update the solution
             xx = xx + B
			 
!             write(6,*) B(1), B(2), B(3)
!             write(6,*) B(4), B(5), B(6)
!             write(6,*) B(7), B(8), B(9)
!             write(6,*) B(10), B(11), B(12)
!             write(6,*) B(13)

!			 stop
!            evaluate residuals
             epsF = sqrt(dot3(B,B,13)) 
!             write(6,*) TolG, TolF
             DEPS(1,1)=1.0d0*xx(1)
             DEPS(2,2)=1.0d0*xx(2)
             DEPS(3,3)=1.0d0*xx(3)
             DEPS(1,2)=1.0d0*xx(4)
             DEPS(2,1)=1.0d0*xx(4)
             DEPS(1,3)=1.0d0*xx(5)
             DEPS(3,1)=1.0d0*xx(5)
             DEPS(2,3)=1.0d0*xx(6)
             DEPS(3,2)=1.0d0*xx(6)
             TrialSTRESS = T + (ELMOD.double.DEPS)
			 
!            Evaluating ep and eq
             dep = xx(7)+xx(8)+xx(9)
             deq = (xx(7)-dep/3.)**2+(xx(8)-dep/3.)**2+
     &             (xx(9)-dep/3.)**2+2.0*xx(10)**2+
     &             2.0*xx(11)**2+2.0*xx(12)**2
             deq = sqrt(2./3.*deq)
!            Getting ps
             ps = STATEV(8)*(1.+rhos*(-dep+ksis*deq))
!             ps = STATEV(8)*(1.+rhos*(-dep))
!            Getting pm
             pm = STATEV(9)*(1.+rhom*(-abs(dep)+ksim*deq))
             pm=min(pm,0.0)			 
             evoid=STATEV(1)+(1.0+STATEV(1))*
     &      (xx(1)+xx(2)+xx(3)+xx(7)+xx(8)+xx(9))
!             pm = STATEV(9)
			 
!            Correcting pc
             param(1) = ps
!            Correcting pt
             param(2) = pm
			 
!             write(6,*) NOEL,NPT,KSTEP,KINC
!             write(6,*) kcounter, epsF, epsG
       
! --------------------------------------------------------------- 
! --------------------------------------------------------------- 
             if (epsF<TolF.AND.epsG<TolG) then
                 exit         
             endif 
             if ( kcounter == maxiter-3) then
                 write(*,*) 'MaxIteration: ', epsF, epsG
            endif 
         enddo  
!         write(6,*) kcounter
!     End of iterations 
! ---------------------------------------------------------------
!     Increment of plastic strains DepsP
! ---------------------------------------------------------------   

! ---------------------------------------------------------------

! --------------------------------------------------------------- 
!     Saving State variables 
!     Alpha isotropic as state variable (scalar)
         STATEV(1)=evoid
!     Plastic strains as state variable (tensor R3xR3)
!         CALL Matrixtovector(STATEV(2:7),epsp)
         STATEV(2:4) = STATEV(2:4) + xx(7:9)
         STATEV(5:7) = STATEV(5:7) + 2.0*xx(10:12)
         STATEV(8)   = ps
         STATEV(9)   = pm
         STATEV(10)  = pc
		 STATEV(11:16) = dervs
! ---------------------------------------------------------------

! ---------------------------------------------------------------  
!     Consistent elastoplastic modulus
         T=TrialSTRESS
!     Calling the derivatives
         call DerivF(TrialSTRESS,xx,NTENS,
     &            PROPS,NPROPS,param,derve,dervp,dervs)
!     Starting with the denominator
!     Will be needing dfds a lot...
         call Vectortomatrix(dervs,dfds) 
!     and dfdep
         call Vectortomatrix(dervp,dfdep) 
         aux = ELMOD.double.dfds
         denom = dfds.double.aux 
!         write(6,*) denom  
         dum=dfdep.double.dfds
         denom = denom - dum
!         write(6,*) dervs(1),dervp(1),dfdep.double.dfds,denom
         dfdep = dfds.double.ELMOD

         Hplast = aux.dyad.dfdep
!         write(6,*) Hplast(1,1,1,1)

         JAC=ELMOD - Hplast/denom
!         write(6,*) xx(13)		 
          
!         JAC=ELMOD  
      endif
! ---------------------------------------------------------------
      Call Solution(NTENS, NDI, NSHR, T, STRESS, JAC
     1 , DDSDDE) 
	 
      !call give_nl(STATEV,NSTATEV,NOEL,NPT)
      ! give nonlocal behavior. matteo told me to comment it
! --------------------------------------------------------------
      END SUBROUTINE CompLib_New




      SUBROUTINE Fyield(STRESS,STATEV,NSTATEV,NTENS,
     &  PROPS,NPROPS,param,ff)
C
      USE tensors
      implicit none
C --------------------------------------------------------------- 
      integer ntens, ndi, nshr, nstatev, nprops, noel, npt,
     & layer, kspt, kstep, kinc
      double precision stress(3,3), statev(nstatev),props(nprops),
     & param(2), ff  
C --------------------------------------------------------------- 
C     Declarating other vaiables and constants    
      double precision p,q,uu,TrialTh,TrialGth,Mf,m,cc,alpha,beta,
     & pc,pt
C ---------------------------------------------------------------        

!     other material properties: alpha,beta...	
      alpha =  PROPS(3) !alpha variable of the exponential
      beta  =  PROPS(4) !beta variable of the exponential	
      Mf    =  PROPS(5) !Mf
      cc    =  PROPS(6) !c for the Lode angle
      m     =  PROPS(7)  

      pc    = param(1)+param(2)
      pt    = m*param(2)

      call pq(STRESS,p,q)
!     Lode angle 
      TrialTh  = cos3theta1(STRESS)  
      TrialGth = 2.0*cc/((1.0+cc)-(1.0-cc)*TrialTh)

!     auxilliary variable for the yield surface
      uu = exp(-((p-pt)/(pc-pt)-alpha)**2/beta)
      ff= q**2.0d0+Mf**2*TrialGth**2*uu*(p-pt)*(p-pc)
      ff=ff/PROPS(12)

      return
! ---------------------------------------------------------------           
      END SUBROUTINE Fyield    
	  
	  
	  
	  
	  
	  
	  
	  
	  
      SUBROUTINE DerivF(STRESS,xx,NTENS,
     &  PROPS,NPROPS,param,derve,dervp,dervs)
C
      USE tensors
      implicit none
C --------------------------------------------------------------- 
      integer ntens, ndi, nshr, nstatev, nprops, noel, npt,
     & layer, kspt, kstep, kinc
      double precision stress(3,3), xx(13),props(nprops),
     & param(2), ff,derve(6),dervp(6),dervs(6),ders(3,3) 
C --------------------------------------------------------------- 
C     Declarating other variables and constants    
      double precision p,q,J2s,J3s,uu,TrialTh,TrialGth,STRAINDEV(3,3),
     & Mf,cc,alpha,beta,pc,pt,derinv(6),ELMOD(3,3,3,3),Idev(3,3,3,3),
     & STRESSDEV(3,3),temp2(3,3),Iunit(3,3,3,3),tempder(3,3)
      double precision dfdps,dfdpm,ps,pm,ep,eq,Iun(3,3),strain(3,3),TINY
      DOUBLE PRECISION E,nu,m,rhos,ksis,rhom,ksim,pc0	  
C ---------------------------------------------------------------        
!     initializing parameters
      TINY  =  1.0D-14 
!     Defining material parameters 
      E     =  PROPS(1) !young modulus
      nu    =  PROPS(2) !poisson modulus 
!     other material properties: alpha,beta...	
      alpha =  PROPS(3) !alpha variable of the exponential
      beta  =  PROPS(4) !beta variable of the exponential	
      Mf    =  PROPS(5) !Mf
      cc    =  PROPS(6) !c for the Lode angle
      m     =  PROPS(7) !m value for pt
      rhos  =  PROPS(8) !rhos
      ksis  =  PROPS(9) !ksis
      rhom  =  PROPS(10)!rhom 
      ksim  =  PROPS(11)!ksim 
      pc0   =  PROPS(12)!dummy variable to make units fit

      ps    = param(1)
      pm    = param(2)
      pc    = param(1)+param(2)
      pt    = m*param(2)  
!     initializing derivatives
      derve = 0.0
      dervp = 0.0
      dervs = 0.0
      derinv = 0.0

      call ELMOD1(PROPS(1),PROPS(2), ELMOD) 
      call pq(STRESS,p,q)
!     Evaluating plastic strain increment invariants
      strain(1,1)=xx(7)
      strain(2,2)=xx(8)
      strain(3,3)=xx(9) 
      strain(1,2)=xx(10)
      strain(1,3)=xx(11)
      strain(2,3)=xx(12) 
      strain(2,1)=xx(10)
      strain(3,1)=xx(11)
      strain(3,2)=xx(12)  
      call pq(strain,ep,eq) 
      ep = 3.0*ep  
      eq = 2.0*eq/3.0 
      Iun = 1.0
      call J23s(STRESS,J2s,J3s)
      call Idev1(Idev) 
      STRESSDEV=double3333x33(Idev,STRESS)
	  STRAINDEV=double3333x33(Idev,strain)
!     Lode angle 
      TrialTh  = cos3theta1(STRESS)  
      TrialGth = 2.0*cc/((1.0+cc)-(1.0-cc)*TrialTh)
!     auxilliary variable for the yield surface
      uu = exp(-((p-pt)/(pc-pt)-alpha)**2/beta)
      ff= q**2.0d0+Mf**2*TrialGth**2*uu*(p-pt)*(p-pc)
!     evaluating the derivatives with respect to the invariants
!     p
      derinv(1) = Mf**2*TrialGth**2*uu*(2*p-pt-pc)
     & -2.0*(ff-q**2)*((p-pt)/(pc-pt)-alpha)/(beta*(pc-pt))
!     q
      derinv(2) = 2.0
!     J2s
      derinv(3) = -(ff-q**2)*(9*sqrt(3.0)*(1.0-cc)*J3s)/
     &(2.0*(1.0+cc)*J2s**(2.5)-sqrt(27.0)*(1.0-cc)*J2s*J3s)
!     J3s 
      derinv(4) =  (ff-q**2)*(6.0*sqrt(3.0)*(1.0-cc))/
     &(2.0*(1.0+cc)*J2s**(1.5)-sqrt(27.0)*(1.0-cc)*J3s)
!     pt
      derinv(5) =-Mf**2*TrialGth**2*uu*(p-pc)
     &-2.0*(ff-q**2)*(p-pc)*((p-pt)/(pc-pt)-alpha)/(beta*(pc-pt)**2)
!     pc
      derinv(6) =-Mf**2*TrialGth**2*uu*(p-pt)
     &+2.0*(ff-q**2)*(p-pt)*((p-pt)/(pc-pt)-alpha)/(beta*(pc-pt)**2)

!     correcting for units
      derinv = derinv/PROPS(12)

!     evaluating the derivatives with respect to the plastic strains
      dervp = 0.0
      dfdps =  derinv(6)*ps*rhos
      dfdpm = (derinv(6)+m*derinv(5))*pm*rhom
!     with respect to ps - softening and hardening
!     volumetric strain
!      dervp(1) = dfdps*(rhos*ps/3.0)derinv(5)*0.0 +derinv(6)/3.0
      dervp(1) =-dfdps
      dervp(2) =-dfdps
      dervp(3) =-dfdps
      dervp(4) = 0.0
      dervp(5) = 0.0
      dervp(6) = 0.0
!     deviatoric strain
      if (eq>TINY) then
         temp2 = double33x3333(STRAINDEV,Idev)
!		 temp2(1,2) = temp2(1,2)+temp2(2,1)
!		 temp2(1,3) = temp2(1,3)+temp2(3,1)
!		 temp2(2,3) = temp2(3,2)+temp2(2,3)
         dervp(1) = dervp(1) + 2.0/3.0*dfdps*ksis/eq*temp2(1,1)
         dervp(2) = dervp(2) + 2.0/3.0*dfdps*ksis/eq*temp2(2,2)
         dervp(3) = dervp(3) + 2.0/3.0*dfdps*ksis/eq*temp2(3,3)
         dervp(4) = dervp(4) + 2.0/3.0*dfdps*ksis/eq*temp2(1,2)
         dervp(5) = dervp(5) + 2.0/3.0*dfdps*ksis/eq*temp2(1,3)
         dervp(6) = dervp(6) + 2.0/3.0*dfdps*ksis/eq*temp2(2,3)
      endif
!     with respect to pm - only softening
!     sign(a1,a2) = (sign of a2) * a1
!      dervp(1) = dervp(1) - dfdpm*sign(1.0D0,ep)
!      dervp(2) = dervp(1) - dfdpm*sign(1.0D0,ep)
!      dervp(3) = dervp(1) - dfdpm*sign(1.0D0,ep)

      dervp(1) = dervp(1) - dfdpm*sign(1.0D0,ep)
      dervp(2) = dervp(2) - dfdpm*sign(1.0D0,ep)
      dervp(3) = dervp(3) - dfdpm*sign(1.0D0,ep)  
      dervp(4) = dervp(4) + 0.0
      dervp(5) = dervp(5) + 0.0
      dervp(6) = dervp(6) + 0.0
!     deviatoric strain
      if (eq>TINY) then
         temp2 = double33x3333(STRAINDEV,Idev)
!		 temp2(1,2) = temp2(1,2)+temp2(2,1)
!		 temp2(1,3) = temp2(1,3)+temp2(3,1)
!		 temp2(2,3) = temp2(3,2)+temp2(2,3)		 
         dervp(1) = dervp(1) + 2.0/3.0*dfdpm*ksim/eq*temp2(1,1)
         dervp(2) = dervp(2) + 2.0/3.0*dfdpm*ksim/eq*temp2(2,2)
         dervp(3) = dervp(3) + 2.0/3.0*dfdpm*ksim/eq*temp2(3,3)
         dervp(4) = dervp(4) + 2.0/3.0*dfdpm*ksim/eq*temp2(1,2)
         dervp(5) = dervp(5) + 2.0/3.0*dfdpm*ksim/eq*temp2(1,3)
         dervp(6) = dervp(6) + 2.0/3.0*dfdpm*ksim/eq*temp2(2,3)
      endif
!     evaluating the derivatives with respect to the stresses
      temp2 = double33x3333(STRESSDEV,Idev)
!	  temp2(1,2) = temp2(1,2)+temp2(2,1)
!	  temp2(1,3) = temp2(1,3)+temp2(3,1)
!	  temp2(2,3) = temp2(3,2)+temp2(2,3)	  
!      call Iunit1(Iunit)
     
      tempder = delta*derinv(1)/3.0+3./2.*derinv(2)*temp2
!      tempder = delta*derinv(1)/3.0+3./2.*derinv(2)*STRESSDEV  

      dervs(1) = tempder(1,1)
      dervs(2) = tempder(2,2)
      dervs(3) = tempder(3,3)
      dervs(4) = tempder(1,2)
      dervs(5) = tempder(1,3)
      dervs(6) = tempder(2,3)

      temp2    = double33x3333(tempder,ELMOD)
!	  temp2(1,2) = temp2(1,2)+temp2(2,1)
!	  temp2(1,3) = temp2(1,3)+temp2(3,1)
!	  temp2(2,3) = temp2(3,2)+temp2(2,3)	  

      derve(1) = temp2(1,1)
      derve(2) = temp2(2,2)
      derve(3) = temp2(3,3)
      derve(4) = temp2(1,2)
      derve(5) = temp2(1,3)
      derve(6) = temp2(2,3)
	  
      return
C ---------------------------------------------------------------           
      END SUBROUTINE DerivF 
	  
	  
	  
	  
	  
      SUBROUTINE ludcmp(a,n,np,indx,d)
	  implicit none
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
!     Largest expected n, and a small number.
      PARAMETER (NMAX=500,TINY=1.0e-20) 
!     Given a matrix a(1:n,1:n), with physical dimension np by np, 
!     this routine replaces it by the LU decomposition of a rowwise 
!     permutation of itself. a and n are input. a is output,
!     arranged as in equation (2.3.14) above; indx(1:n) is an output vector 
!     that records the row permutation effected by the partial pivoting; d is 
!     output as ±1 depending on whether the number of row interchanges was 
!     even or odd, respectively. This routine is used in combination 
!     with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
!     vv stores the implicit scaling of each row.
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX) 
      d=1. !No row interchanges yet.
      do 12 i=1,n !Loop over rows to get the implicit scaling information
         aamax=0.
         do 11 j=1,n
             if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      enddo
!        No nonzero largest element.
         if (aamax.lt.TINY) stop 'singular matrix in ludcmp'
         vv(i)=1./aamax !Save the scaling.
 12    enddo
       do 19 j=1,n !This is the loop over columns of Crout’s method.
!      This is equation (2.3.12) except for i = j.
         do 14 i=1,j-1 
             sum=a(i,j)
             do 13 k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
 13          enddo 
             a(i,j)=sum
 14      enddo 
!        Initialize for the search for largest pivot element.
         aamax=0. 
         do 16 i=j,n
!        This is i = j of equation (2.3.12) and i = j+1. . .N
!        of equation (2.3.13).
             sum=a(i,j)  
             do 15 k=1,j-1
                 sum=sum-a(i,k)*a(k,j)
 15          enddo 
             a(i,j)=sum
             dum=vv(i)*abs(sum) !Figure of merit for the pivot.
             if (dum.ge.aamax) then !Is it better than the best so far?
                 imax=i
                 aamax=dum
             endif
 16      enddo 
         if (j.ne.imax)then !Do we need to interchange rows?
             do 17 k=1,n !Yes, do so...
                 dum=a(imax,k)
                 a(imax,k)=a(j,k)
                 a(j,k)=dum
 17          enddo 
             d=-d !...and change the parity of d.
             vv(imax)=vv(j) !Also interchange the scale factor.
         endif
         indx(j)=imax
         if(a(j,j).eq.0.) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to
!        the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to
!        substitute TINY for zero.
         if(j.ne.n)then !Now, finally, divide by the pivot element.
             dum=1./a(j,j)
             do 18 i=j+1,n
                 a(i,j)=a(i,j)*dum
 18          enddo 
         endif
 19   enddo  !Go back for the next column in the reduction.
      return
      END

		 
		 
		 
		 
		 
		 
		 
		 
		 
      SUBROUTINE lubksb(a,n,np,indx,b)
	  implicit none
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
!     Solves the set of n linear equations A · X = B. Here a is input, not as
!     the matrix A but rather as its LU decomposition, determined by the 
!     routine ludcmp. indx is input as the permutation vector returned by 
!     ludcmp. b(1:n) is input as the right-hand side vector B,
!     and returns with the solution vector X. a, n, np, and indx are not 
!     modified by this routine and can be left in place for successive calls 
!     with different right-hand sides b. This routine takes into account the 
!     possibility that b will begin with many zero elements, so it is 
!     efficient for use in matrix inversion.
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0 !When ii is set to a positive value, it will become the index
!     of the first nonvanishing element of b. We now do
!     the forward substitution, equation (2.3.6). The only new
!     wrinkle is to unscramble the permutation as we go.
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
             do 11 j=ii,i-1
                 sum=sum-a(i,j)*b(j)
 11          enddo 
         else if (sum.ne.0.) then
             ii=i !A nonzero element was encountered, so from now on we will
         endif !have to do the sums in the loop above.
         b(i)=sum
 12   enddo 
      do 14 i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
         sum=b(i)
         do 13 j=i+1,n
             sum=sum-a(i,j)*b(j)
 13      enddo
         b(i)=sum/a(i,i) !Store a component of the solution vector X.
 14   enddo 
      return !All done!
      END
	  
	  
	  
	  
	  
	  
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
C

       implicit none
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      integer nstatv,nfield,NOEL,NPT,LAYER,JRCD,
     1 KSPT,KSTEP,KINC,NDI,NSHR,MATLAYO,LACCFLA
      double precision FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2),ARRAY(15),COORD(*),PNEWDT,CELENT,DTIME
      integer JARRAY(15),JMAC(*),JMATYP(*)
      double precision n1,rho1,L1,T1,T1E,PHI1,P11,P12,P22,PI	
      double precision n2,rho2,L2,T2,T2E,PHI2,RA1,RA2	  
      double precision POREP, STR(4) ,STRN,STRT,E,ANU,RES,RCOOR
	   

  14  RETURN
      END	  
