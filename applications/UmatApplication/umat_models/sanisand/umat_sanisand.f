c------------------------------------------------------------------------------
c complete file suite for Dafalias & Manzari (2004) SANISAND model for sand
c------------------------------------------------------------------------------
      subroutine sanisand_umat(stress,statev,ddsdde,sse,spd,scd,
     &  rpl,ddsddt,drplde,drpldt,
     &  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c------------------------------------------------------------------------------
c user subroutine for Abaqus 6.3
c------------------------------------------------------------------------------
c
c Implemented constitutive law:
c -----------------------------
c Dafalias & Manzari(2004) SANISAND model for sand
c
c Ref: Dafalias, Y. F. and Manzari, M. T.  
c      Simple plasticity sand model accounting for fabric change effects
c      J. Engng. Mechanics, ASCE (2004), 130(6):622-634.
c
c ----------------------------------------------------------------------------
c The string for the material name may contain 9 characters.
c ----------------------------------------------------------------------------
c Material constants:
c   	
c     ---------------------------------------------------------------------
c     props(j)      
c     ---------------------------------------------------------------------
c        1     p_a        Atmospheric pressure 
c        2     e0         Void ratio on CSL at p = 0  
c        3     lambda     CSL parameter (e:p plane)
c        4     xi         CSL parameter (e:p plane)
c        5     M_c        Slope of CSL in q:p plane, TX compression
c        6     M_e        Slope of CSL in q:p plane, TX extension
c        7     mm         opening of yield surface cone
c        8     G0         Shear modulus constant
c        9     nu         Poisson's ratio
c        10    h0         Plastic modulus constant
c        11    c_h        Plastic modulus constant
c        12    n_b        Plastic modulus constant
c        13    A0         Dilatancy constant
c        14    n_d        Dilatancy constant
c        15    z_max      Fabric index constant
c        16    c_z        Fabric index constant
c        17    bulk_w     Pore water bulk modulus (undrained conditions)
c	 18    p_tmult    shift of mean stress pt=p_tmult*p_a	
c        19    initial value of void ratio
c     ----------------------------------------------------------------------
c
c     Solution dependent state variables (statev):
c     definition via sdvini
c
c     group 1: internal variables (14 variables)
c
c        1 ... alpha_11	  back stress, orientation of yield surface cone
c        2 ... alpha_22
c        3 ... alpha_33
c        4 ... alpha_12
c        5 ... alpha_13
c        6 ... alpha_23
c
c        7 ... void       void ratio
c
c        8 ... Fab_11     fabric tensor z
c        9 ... Fab_22
c       10 ... Fab_33
c       11 ... Fab_12
c       12 ... Fab_13
c       13 ... Fab_23
c
c       14 ... not used        
c
c     group 2: memory variables for shear reversal (SR) and other purposes
c
c       15 ... alpha_sr_11	 alpha value at stress reversal points (discrete update)
c       16 ... alpha_sr_22   
c       17 ... alpha_sr_33   
c       18 ... alpha_sr_12
c       19 ... alpha_sr_13
c       20 ... alpha_sr_23
c
c       21 ... not used
c       22 ... not used
c       23 ... not used
c       24 ... not used
c       25 ... not used
c       26 ... not used
c       27 ... not used
c
c       28 ... not used 
c
c     group 3: variables saved for post processing or other purposes
c
c       29 ... pore	    excess pore pressure (undrained case)
c       30 ... p	    mean effective stress
c       31 ... q	    deviator stress
c       32 ... z	    Lode parameter (cos(3theta))
c       33 ... dtsub	suggested size of first time substep
c       34 ... nfev	    number of function evaluation
c       35 ... not used
c       36 ... not used
c
c Authors: 
c     C. Tamagnini (tamag@unipg.it)
c     Dipartimento di Ingegneria Civile e Ambientale 
c     Università degli Studi di Perugia, Italy
c
c     M. Martinelli
c     Dipartimento di Ingegneria Strutturale e Geotecnica
c     Università di Roma "La Sapienza", Italy
c
c     C. Miriano
c     Dipartimento di Ingegneria Strutturale e Geotecnica
c     Università di Roma "La Sapienza", Italy
c
c Modifications: D. Masin, 2015
c
c NOTES: 
c     - sign convention for stress and strain: tension and extension positive
c     - stress and strain sign convention changed upon entering the SP algorithm
c     - tangent stiffness operator evaluated according to two alternative options
c       selected setting the logical flag "cons_lin":
c       cons_lin.eq.0  -> numerical linearization via
c                            direct perturbation of dstran
c       cons_lin.eq.1 -> continuum tangent stiffness
c                            (not optimal for full N-R iterative solver)
c       cons_lin.eq.2 -> elastic tangent stiffness
c                            (even less optimal for full N-R iterative solver, but sometimes more stable)
c
c
c Last change: 4/2013  
c
c----------------------------------------------------------------------------
c
      implicit none
c
      character*80 cmname
c
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc
c
      double precision stress(ntens), statev(nstatv),
     &  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),
     &  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     &  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp, 
     &  dtemp, pnewdt, celent
c
c ... 1. nasvdim    = maximum number of additional state variables
c     2. tolintT    = prescribed error tolerance for the adaptive 
c                     substepping scheme
c     3. maxnint    = maximum number of time substeps allowed.
c                     If the limit is exceeded abaqus is forced to reduce 
c                     the overall time step size (cut-back) 
c     4. DTmin      = minimum substeps size allowed.
c                     If the limit is exceeded abaqus is forced to reduce 
c                     the overall time step size (cut-back)
c     5. perturb    = perturbation par. for computation of Jacobian matrices
c     6. nfasv      = number of first additional state variable in statev field 
c     7. prsw       = switch for printing information
c
c ... declaration of local variables
c
      integer prsw,elprsw,cons_lin,abaqus,chiara,check_ff,drcor
c
      integer i,error,maxnint,nfev,mario_DT_test,inittension
      integer nparms,nasvdim,nfasv,nydim,nzdim
      integer nasvy,nasvz,nyact,nzact,plastic,testing
c
      double precision dot_vect
c	
      double precision parms(nprops),theta,tolintT,dtsub,DTmin,perturb
      double precision sig_n(6),sig_np1(6),DDtan(6,6),pore,PI
      double precision deps_np1(6),depsv_np1,norm_D2,norm_D,tolintTtest
      double precision eps_n(6),epsv_n,alphayield(6)
      double precision norm_deps2,norm_deps,pp,qq,cos3t,ddum
      double precision zero,tol_f,fact_thres,p_thres,stran_lim,eps_debug
      double precision p_atm,ptshift,phimob,tol_f_test,youngel,nuel
      double precision avoid,apsi,aec,yf_DM,fyield,psi_void_DM,Mb
      double precision dummy,sdev(6),I1,alpha(6),cM,tau(6),gth,etanorm
      double precision sinphinorm
c
      parameter (nasvdim  = 36)
      parameter (nydim    = 6+14)
      parameter (nzdim    = 14)
      parameter (tolintT  = 1.00d-3)
      parameter (tolintTtest = 1.0d-2) 
c
      parameter (maxnint  = 50000)
      parameter (DTmin    = 1.0d-18)
      parameter (perturb  = 1.0d-4)
      parameter (nfasv    = 1)
      parameter (prsw     = 0)
      parameter (cons_lin = 1)
	parameter (abaqus = 0)
c chiara
 	parameter (eps_debug = 0.9d-3)
c
      parameter (zero = 0.0d0)
      parameter (PI = 3.14159265358979323846264338327950288)
      parameter (fact_thres=0.000000001d0)
c
c ... additional state variables
c
      double precision  asv1(nydim-6),asv2(nzdim)
c
c ... solution vector (stresses, additional state variables)
c
      double precision  y(nydim),y_n(nydim),z(nzdim),z_n(nzdim)
c
c      common /z_nct_errcode/error
c	common /z_tolerance/tol_f
c	common /z_check_yeld/check_ff
c	common /z_drift_correction/drcor
c	common /z_threshold_pressure/p_thres
	
c
        tol_f=1.0d-6
        tol_f_test=1.0d-6
	check_ff=0
	drcor=1
	plastic=0
	phimob=0.0d0
	ptshift=0.0d0
c
c
c ... Error Management:
c     ----------------
c     error =  0 ... no problem in time integration
c     error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
c                    stress state), reduce time integration substeps
c     error =  3 ... problems in time integration, reduce abaqus load increment 
c                    (cut-back)
c     error=10 ... severe error, terminate calculation
c
      error=0
c
c ... check problem dimensions
c
      if (ndi.ne.3) then
c
        write(6,*) 'ERROR: this UMAT can be used only for elements'
        write(6,*) '       with 3 direct stress/strain components'
        write(6,*) 'noel = ',noel
        error=10
c
      endif
c      open(unit=6,position='Append',file=
c     .	'C:/users/david/data/Zhejiang/plaxis-sani/UMATdebug.txt')
c
c ... check material parameters and move them to array parms
c
      nparms=nprops
      call check_parms_DM(props,parms,nparms)
c
c ... print informations about time integration, useful when problems occur
c
	p_atm=parms(1)
	p_thres=fact_thres*p_atm
c
      elprsw = 0
      if (prsw .ne. 0) then
c
c ... print only in some defined elements
c
c	 if ((noel.eq.101).and.(npt.eq.1)) elprsw = 1
c
      endif
c
c ... define number of additional state variables
c
      call define(nasvy,nasvz)
      nyact = 6 + nasvy
      nzact = nasvz
      if (nyact.gt.nydim) then
        write(6,*) 'ERROR:nydim too small; program terminated'
        error=10
      elseif (nzact.gt.nzdim) then
        write(6,*) 'ERROR:nzdim too small; program terminated'
        error=10
      endif
c
c ... suggested time substep size, and initial excess pore pressure
c
      pore = statev(29)
c
c ... changes sign conventions for stresses and strains, compute 
c     current effective stress tensor and volumetric strain
c
      ptshift=parms(18)*parms(1)
      do i=1,3
          stress(i) = stress(i)-ptshift
      enddo

      call move_sig(stress,ntens,-1*ptshift,sig_n)
	
      call move_sig(stress,ntens,pore,sig_n)
      call move_eps(dstran,ntens,deps_np1,depsv_np1)
      call move_eps(stran,ntens,eps_n,epsv_n)
      
      norm_D2=dot_vect(2,deps_np1,deps_np1,6)
      norm_D=sqrt(norm_D2)
c chiara
	if (eps_n(1).gt.eps_debug) then
		chiara=1
	end if
c
c ... initialise void ratio and yield surface inclination
c
      if (statev(7) .lt. 0.001) then
            do i=1,6        
               alphayield(i)=zero
            end do
            call deviator(sig_n,alphayield,ddum,pp)
      	    avoid=0
            if(parms(19) .le. 5.0) then 
                   avoid=parms(19)
            else if(parms(19) .gt. 5.0) then
        	   apsi=parms(19)-10.0d0
        	   aec=parms(2)-parms(3)*(pp/parms(1))**parms(4)
        	   avoid=aec+apsi
            endif
            statev(7)=avoid
            do i=1,6        
               statev(i)=alphayield(i)/pp
               statev(i+14)=alphayield(i)/pp
            end do
      end if
c
c ... move vector of additional state variables into asv1(nydim) and asv2(nzdim)
c
      do i=1,nasvy
        asv1(i) = statev(i-1+nfasv)
      enddo
c
      do i=1,nasvz
        asv2(i) = statev(i-1+nfasv+nasvy)
      enddo
c
c --------------------
c ... Time integration
c --------------------
c
      call iniyz(y,nydim,z,nzdim,asv1,nasvy,asv2,nasvz,sig_n,ntens)
      call push(y,y_n,nydim)
      call push(z,z_n,nzdim)
c
      if (elprsw.ne.0) then
        write(6,*) '================================================='
        write(6,*) '         Call of UMAT - DM SANISAND model:       '
        write(6,*) '================================================='
        call wrista(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &              parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
      endif
c
c ... local integration using adaptive RKF-23 method with error control
c
      dtsub = dtime
      if((dtsub.le.zero).or.(dtsub.gt.dtime)) then
        dtsub = dtime
      end if
      
      testing=0
c     For use in PLAXIS, activate the following line
c      if(kstep.eq.1 .AND. kinc.eq.1) testing=1
c     For use in ABAQUS, the line above should be inactive
	
      if(norm_D.eq.0) testing=2
c     FEM asking for ddsdde only

      nfev = 0 ! initialisation

      if(testing.eq.1) then
        call rkf23_upd_DM(y,z,nyact,nasvy,nasvz,tolintTtest,maxnint,
     &         DTmin,deps_np1,parms,nparms,nfev,elprsw,
     &	       mario_DT_test,
     &         error,tol_f_test,check_ff,drcor,p_thres,plastic)
c ... give original state if the model fails without substepping
          if(error.ne.0) then
            do i=1,nyact        
               y(i)=y_n(i)
            end do
            error=0
          end if
      else if(testing.eq.2) then
            do i=1,nyact        
                  y(i)=y_n(i)
            end do
c ... Normal RKF23 integration
      else   !testing.eq.0
        call rkf23_upd_DM(y,z,nyact,nasvy,nasvz,tolintT,maxnint,
     &         DTmin,deps_np1,parms,nparms,nfev,elprsw,
     &	       mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
      end if
c
c ... error conditions (if any)
c
	if(mario_DT_test.eq.1) then
	call wrista(4,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &            parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
	endif
c
      if(error.eq.3) then
c
c ... reduce abaqus load increment
c
        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &            parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
        write(6,*) 'subroutine UMAT: reduce step size in ABAQUS'
        write(6,*) 'error 3 activated'

        if(abaqus.ne.0) then
            pnewdt = 0.25d0
        else
c ... write a message and return the original state
           do i=1,nyact        
                  y(i)=y_n(i)
           end do     
        endif
c			
c        if(abaqus.eq.0) then
c                write(6,*) 'analysis ended because number of time ' 
c		  write(6,*) 'substeps exceeded maximum number allowed'
c		  write(6,*) ' (maxnint)'
c	      call xit_DM
c	  endif
  
	return
c
      elseif(error.eq.10) then
      	write(6,*) 'error 10 activated'
c
        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &              parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
        call xit_DM
c
      endif
c
c ... update dtsub and nfev
c
      if(dtsub.le.0.0d0) then 
      	dtsub = 0
      else if(dtsub.ge.dtime) then 
      	dtsub = dtime
      end if
      statev(33)=dtsub
      statev(34)=dfloat(nfev)
c
c ... computation of Jacobian matrix
c
      error=0
      if(cons_lin.eq.0) then
c
c ... parameter of the numerical differentiation
c     double precision
c
        norm_deps2=dot_vect(2,deps_np1,deps_np1,ntens)
        norm_deps=dsqrt(norm_deps2)
        theta=perturb*max(norm_deps,1.0d-6)
c
c ... compute consistent tangent via numerical perturbation
c
        call pert_DM(y_n,y,z,nyact,nasvy,nasvz,tolintT,maxnint,DTmin,
     &       deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
      else
c        
        call tang_stiff(y,z,nyact,nasvy,nasvz,parms,nparms,
     &         DDtan,cons_lin,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
      endif
c
	if(error.ne.0) then
	 write(6,*) 'some error in pert or tang_stiff, value=',error
c
c        call wrista(2,y,nydim,deps_np1,dtime,coords,statev,nstatv,
c     &              parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c        call xit_DM
c
      endif
c
c ... convert solution (stress + cons. tangent) to abaqus format
c     update pore pressure and compute total stresses

      inittension=0
c     just checks for NAN
      call check_RKF_DM(inittension,y,nyact,nasvy,parms,nparms)
      if (inittension.ne.0) then
           do i=1,nyact        
                  y(i)=y_n(i)
           end do        
      end if
	
c
      call solout(stress,ntens,asv1,nasvy,asv2,nasvz,ddsdde,
     +            y,nydim,z,pore,depsv_np1,parms,nparms,DDtan)
c
c ... updated vector of additional state variables to abaqus statev vector
c
      do i=1,nasvy
        statev(i-1+nfasv) = asv1(i) 
      end do
c
      do i=1,nasvz
        statev(i-1+nfasv+nasvy) = asv2(i)
      enddo
c
c ... transfer additional information to statev vector
c
      do i=1,6
        sig_np1(i)=y(i)
      end do
      call inv_sig(sig_np1,pp,qq,cos3t)
c
      statev(29) = pore 
      statev(30) = pp
      statev(31) = qq
      statev(32) = cos3t
     
      cM=parms(6)/parms(5)
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
      call deviator(sig_np1,sdev,I1,pp)
      do i=1,6
        tau(i)=sdev(i)-pp*alpha(i)
      end do
      call lode_DM(tau,cM,cos3t,gth,dummy)
      etanorm=gth*qq/pp
      sinphinorm=3*etanorm/(6+etanorm)
      statev(33) = asin(sinphinorm)*180/PI
      statev(34) = nfev
      
c      check that bounding surtface is not violated 
c      if (noel.eq.324 .and. npt.eq.3) then
c      	    fyield=yf_DM(y,nyact,parms,nparms)  
c      	    apsi=psi_void_DM(statev(7),pp,parms,nparms)
c      	    Mb=parms(5)*dexp(-parms(12)*apsi)
c      	    write(6,*) 'fyield=',fyield
c      	    write(6,*) 'psi=',apsi
c      	    write(6,*) 'Mb=',Mb
c      	    write(6,*) 'qq/pp*gth=',qq/pp*gth
c      	    write(6,*) '---------------------'
c      end if
      
      do i=1,3
          stress(i) = stress(i)+ptshift
      enddo
c      close(6)
c
c -----------------------
c End of time integration
c -----------------------
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine alpha_th_DM(flag,n,gth,psi,parms,nparms,alpha)
c-----------------------------------------------------------------------------
c	calculate tensors:
c
c     alpha_c => flag = 1
c     alpha_b => flag = 2
c     alpha_d => flag = 3
c
c   Dafalias & Manzari (2004) SANISAND model for sands
c
c   variables allocated in parms
c
c   1     p_a        Atmospheric pressure 
c   2     e0         Void ratio on CSL at p = 0  
c   3     lambda     CSL parameter (e:p plane)
c   4     xi         CSL parameter (e:p plane)
c   5     M_c        Slope of CSL in q:p plane, TX compression
c   6     M_e        Slope of CSL in q:p plane, TX extension
c   7     mm         opening of yield surface cone
c   8     G0         Shear modulus constant
c   9     nu         Poisson's ratio
c   10    h0         Plastic modulus constant
c   11    c_h        Plastic modulus constant
c   12    n_b        Plastic modulus constant
c   13    A0         Dilatancy constant
c   14    n_d        Dilatancy constant
c   15    z_max      Fabric index constant
c   16    c_z        Fabric index constant
c   17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c	written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c
      integer flag,nparms,i
c
      double precision n(6),gth,psi,parms(nparms),alpha(6)
      double precision M_c,mm,n_b,n_d
c
      double precision M,alpha_th
      double precision two,three,sqrt23
c
      data two,three/2.0d0,3.0d0/
c
      sqrt23=dsqrt(two/three)
c
c ... recover material parameters
c
      M_c=parms(5)
	mm=parms(7)
	n_b=parms(12)
	n_d=parms(14)
c
c ... select which alpha tensor to evaluate
c
      if(flag.eq.1) then
c
c ... critical state cone
c
        M=M_c
c
      elseif(flag.eq.2) then 
c
c ... bounding surface cone
c
        M=M_c*dexp(-n_b*psi)
c
      else
c
c ... dilatancy cone
c
        M=M_c*dexp(n_d*psi)
c
      endif
c
c ... tensor alpha_ij
c
      alpha_th=M*gth-mm 
c
      do i=1,6
        alpha(i)=sqrt23*alpha_th*n(i)      
      end do
c	
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine check_crossing(y,y_tr,n,parms,nparms,prod)
c-----------------------------------------------------------------------------
c
c  computes 
c
c  prod := dsig_tr*grad(f)_k 
c
c  useful for checking if crossing of yield locus occurs whenever 
c  f_k = 0 and f_tr > 0 
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c-----------------------------------------------------------------------------
      implicit none
c
      integer i,n,nparms
c
      double precision dot_vect
c
      double precision y(n),y_tr(n),parms(nparms)
      double precision P(6),P1(6),dsig_tr(6)
      double precision prod
c
c ... gradient of yield surface at state y_k
c
      call grad_f_DM(y,n,parms,nparms,P,P1)
c
      do i=1,6
        dsig_tr(i)=y_tr(i)-y(i)
      end do ! i
c		  
      prod=dot_vect(1,P,dsig_tr,6)  
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine check_parms_DM(props,parms,nprops)
c-----------------------------------------------------------------------------
c     checks input material parameters for Dafalias & Manzari (2004)  
c     SANISAND model for sand
c
c     Material constants:
c   	
c     ---------------------------------------------------------------------
c     props(j)      
c     ---------------------------------------------------------------------
c        1     p_a        Atmospheric pressure 
c        2     e0         Void ratio on CSL at p = 0  
c        3     lambda     CSL parameter (e:p plane)
c        4     xi         CSL parameter (e:p plane)
c        5     M_c        Slope of CSL in q:p plane, TX compression
c	   6     M_e        Slope of CSL in q:p plane, TX extension
c        7     mm         opening of yield surface cone
c        8     G0         Shear modulus constant
c        9     nu         Poisson's ratio
c        10    h0         Plastic modulus constant
c        11    c_h        Plastic modulus constant
c        12    n_b        Plastic modulus constant
c        13    A0         Dilatancy constant
c        14    n_d        Dilatancy constant
c	   15    z_max      Fabric index constant
c        16    c_z        Fabric index constant
c        17    bulk_w     Pore water bulk modulus (undrained conditions)
c     ---------------------------------------------------------------------
c
c     Solution dependent state variables (statev):
c     definition via sdvini
c
c     group 1: internal variables (14 variables)
c
c        1 ... alpha_11	  back stress, orientation of yield surface cone
c        2 ... alpha_22
c        3 ... alpha_33
c        4 ... alpha_12
c        5 ... alpha_13
c        6 ... alpha_23
c
c        7 ... void       void ratio
c
c        8 ... Fab_11     fabric tensor z
c        9 ... Fab_22
c       10 ... Fab_33
c       11 ... Fab_12
c       12 ... Fab_13
c       13 ... Fab_23
c
c       14 ... not used        
c
c     group 2: memory variables for shear reversal (SR) and other purposes
c
c       15 ... alpha_sr_11	 alpha value at stress reversal points (discrete update)
c       16 ... alpha_sr_22   
c       17 ... alpha_sr_33   
c       18 ... alpha_sr_12
c       19 ... alpha_sr_13
c       20 ... alpha_sr_23
c
c       21 ... not used
c       22 ... not used
c       23 ... not used
c       24 ... not used
c       25 ... not used
c       26 ... not used
c       27 ... not used
c
c       28 ... not used 
c
c     group 3: variables saved for post processing or other purposes
c
c       29 ... pore	    excess pore pressure (undrained case)
c       30 ... p	    mean effective stress
c       31 ... q	    deviator stress
c       32 ... z	    Lode parameter (cos(3theta))
c       33 ... dtsub	suggested size of first time substep
c       34 ... nfev	    number of function evaluation
c       35 ... not used
c       36 ... not used
c
c-----------------------------------------------------------------------------
      implicit none
c
      integer nprops
c
      double precision props(nprops),parms(nprops)

      double precision p_a,e0,lambda,xi,M_c,M_e,mm
      double precision G0,nu,h0,c_h,n_b,A0,n_d,z_max
      double precision c_z,bulk_w,sinphi,PI,sinphiext
c
	  double precision zero
c
      parameter(zero=0.0d0)
      parameter(PI=3.14159265358979323846264338327950288)
c
c ... recover material parameters and initial state info
c
      p_a=props(1)
      e0=props(2)
      lambda=props(3)
	xi=props(4)
      M_c=props(5)
      M_e=props(6)
      mm=props(7)
      G0=props(8)
      nu=props(9) 
      h0=props(10)
      c_h=props(11)
      n_b=props(12)
      A0=props(13)
      n_d=props(14)
      z_max=props(15)
      c_z=props(16)
      bulk_w=props(17)
c
c ... move vector props into local vector parms
c
      call push(props,parms,nprops)

      if(parms(5) .gt. 5) then
      	      sinphi=sin(parms(5)/180*PI)
      	      parms(5)=6*sinphi/(3-sinphi)
      else
      	      sinphi=3*parms(5)/(6+parms(5))
      end if
      if(parms(6) .gt. 5) then
      	      sinphiext=sin(parms(6)/180*PI)
      	      parms(6)=6*sinphiext/(3+sinphiext)
      else if ((parms(6) .le. 5) .and. (parms(6) .gt. 0.01)) then
      	      sinphiext=3*parms(6)/(6-parms(6))
      else
      	      parms(6)=parms(5)*(3-sinphi)/(3+sinphi)
      end if

      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine define(nasvy,nasvz)
c-----------------------------------------------------------------------------
      implicit none 
      integer nasvy,nasvz
c
c number of additional state variables stored in vectors y and z
c must 14 each (otherwise change nasvdim in umat)
c
c        components of ASV(i) stored in y (14 variables)
c
c        1 ... alpha_11	  back stress, orientation of yield surface cone
c        2 ... alpha_22
c        3 ... alpha_33
c        4 ... alpha_12
c        5 ... alpha_13
c        6 ... alpha_23
c        7 ... void       void ratio
c        8 ... Fab_11	  fabric tensor
c        9 ... Fab_22
c       10 ... Fab_33
c       11 ... Fab_12
c       12 ... Fab_13
c       13 ... Fab_23
c       14 ... not used        
c
c       components of ASV(i) stored in z (14 variables)
c
c       15 ... alpha_sr_11	alpha value at stress reversal points (discrete update)
c       16 ... alpha_sr_22  
c       17 ... alpha_sr_33
c       18 ... alpha_sr_12
c       19 ... alpha_sr_13
c       20 ... alpha_sr_23
c       21 ... not used     
c       22 ... not used
c       23 ... not used
c       24 ... not used
c       25 ... not used
c       26 ... not used
c       27 ... not used
c       28 ... not used 
c
      nasvy = 14
      nasvz = 14
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c------------------------------------------------------------------------------
      subroutine deviator(t,s,trace,mean)
c------------------------------------------------------------------------------
c calculate deviator and trace of 2nd order tensor t(6)
c
c NOTE: Voigt notation is used with the following index conversion
c
c       11 -> 1
c       22 -> 2
c       33 -> 3
c       12 -> 4
c       13 -> 5
c       23 -> 6
c
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision t(6),s(6),trace,mean
      double precision one,three,onethird
c
      data one,three/1.0d0,3.0d0/
c
c ... some constants
c
      onethird=one/three
c
c ... trace and mean value
c
      trace=t(1)+t(2)+t(3)
      mean=onethird*trace
c
c ... deviator stress
c
      s(1)=t(1)-mean
      s(2)=t(2)-mean
      s(3)=t(3)-mean
      s(4)=t(4)
      s(5)=t(5)
      s(6)=t(6)
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c------------------------------------------------------------------------------
      double precision function distance(alpha_k,alpha,n)
c------------------------------------------------------------------------------
c computes distance function
c
c     d = (alpha^k_{ij}-alpha_{ij})n_{ij}   (k=sr,b,d)
c
c Dafalias & Manzari (2004) SANISAND model for sands
c
c------------------------------------------------------------------------------
      implicit none
c
      integer i
c
      double precision dot_vect
c
      double precision alpha_k(6),alpha(6),n(6),delta(6)
c
      do i=1,6
        delta(i)=alpha_k(i)-alpha(i)
      end do
c
      distance=dot_vect(1,delta,n,6)
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c------------------------------------------------------------------------------
      double precision function dot_vect(flag,a,b,n)
c------------------------------------------------------------------------------
c dot product of a 2nd order tensor, stored in Voigt notation
c
c flag = 1 -> vectors are stresses in Voigt notation
c flag = 2 -> vectors are strains in Voigt notation
c flag = 3 -> ordinary dot product between R^n vectors
c------------------------------------------------------------------------------
      implicit none
      integer i,n,flag
      double precision a(n),b(n)
      double precision zero,half,one,two,coeff
c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c
      if(flag.eq.1) then
c
c ... stress tensor (or the like)
c
        coeff=two
c
      elseif(flag.eq.2) then
c
c ... strain tensor (or the like)
c
        coeff=half
c
      else
c
c ... standard vectors
c
        coeff=one
c	
      end if
c
      dot_vect=zero
c
      do i=1,n
        if(i.le.3) then
          dot_vect = dot_vect+a(i)*b(i)
        else
          dot_vect = dot_vect+coeff*a(i)*b(i)
        end if
      end do
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine drift_corr_DM(y,n,z,nasvz,parms,nparms,tol,switch2,
     & mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c  performs consistent drift correction (see Sloan et al. 2001)
c  Dafalias & Manzari(2004) SANISAND model for sand
c
c  written 8/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c
      double precision dot_vect
c
	integer switch2,mario_DT_test
c
      external matmul
      double precision yf_DM
c 
      integer n,nasvz,nparms,i,n_drift,max_ndrift,switch
	integer iter, itermax
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

c
      double precision y(n),y0(n),y1(n), z(nasvz),parms(nparms)
      double precision gradf(6),gradf1(6),gradg(6),gradg1(6)
      double precision DDe(6,6),UU(6),VV(6),h_alpha(6),Kpm1,p1,pp1
      double precision f0,tol,zero,one,denom,fnm1,p,three,onethird,f0_p
	double precision factor,f1,p_atm
c
      parameter(zero=0.0d0,one=1.0d0,three=3.0d0)
      parameter(max_ndrift=10000, itermax=1000)
c
c      common /z_nct_errcode/error	
c
c ... initialize constants and vectors
c
      call push(y,y0,n)
c
c ... check if current state is inside the elastic nucleus
c
c Chiara Miriano 15 maggio 2009
c	onethird=one/three
c	p=(y0(1)+y0(2)+y0(3))*onethird
c	if(p.lt.zero) then
c		do i=1,3
c			y0(i)=y0(i)-p
c		end do
c	p=(y0(1)+y0(2)+y0(3))*onethird
c	end if
c Chiara Miriano  15 maggio 2009
c
      f0=yf_DM(y0,n,parms,nparms)
	onethird=one/three
	p=(y0(1)+y0(2)+y0(3))*onethird
	
      n_drift=0
	switch=0
	f0_p=f0/p
c	p_atm=parms(1)
	
c	if(p.lt.(p_atm/100)) f0_p=f0_p/1000000
c	f0_p=f0
	switch2=0
c
      do while(f0_p.gt.tol)
cccc		do while(f0.gt.tol)
c
        fnm1=f0
c
c ... current state outside yield surface, correct it until f0 < ftol
c
        n_drift=n_drift+1
c
c ... elastic stiffness and gradients of f and g
c
        call el_stiff_DM(y0,n,parms,nparms,DDe,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
        call grad_f_DM(y0,n,parms,nparms,gradf,gradf1)
        call grad_g_DM(y0,n,parms,nparms,gradg,gradg1)
c
c ... vectors UU=DDe*gradg and VV=DDe*gradf		
c
	  call matmul(DDe,gradg1,UU,6,6,1)
	  call matmul(DDe,gradf1,VV,6,6,1)
c
c ... hardening function h_alpha and plastic modulus (1/Kp)
c
       call plast_mod_DM(y0,n,z,nasvz,parms,nparms,h_alpha,Kpm1,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
	if (switch2.gt.zero) return
c
        if(one/Kpm1.le.zero) then 
c          write(6,*) 'ERROR: subroutine DRIFT_CORR:'
          write(6,*) 'subcritical softening condition'
c          write(6,*) 'Kp = ',one/Kpm1
c          error=10
	   error=3
c		switch2=1
         return
cc		switch=1
        end if
c
c ... correction for stress (y(1):y(6))
c	
	  if(switch.eq.0) then
		do i=1,6
			y1(i)=y0(i)-Kpm1*f0*UU(i)
		end do
c
c ... correction for hardening variable alpha (y(7):y(12))
c
		do i=1,6
			y1(6+i)=y0(6+i)+Kpm1*f0*h_alpha(i)
		end do
		do i=13,n
			y1(i)=y0(i)
		end do	
c
c ... recompute drift at the new state
c
		f0=yf_DM(y1,n,parms,nparms)
		if(f0.gt.fnm1) then
			switch=1
c	switch2=1
	p1=(y1(1)+y1(2)+y1(3))*onethird
c	write(*,*)'switch2=',switch2
c	write(*,*)'p=',p
c	write(*,*)'p1=',p1
c	return
		else
			call push(y1,y0,n)
		end if		
	  else
c
c ... normal correction in place of consistent correction
c
c		call push(y,y0,n)
		call push(y0,y1,n)
		f0=yf_DM(y0,n,parms,nparms)
c		f1=yf_DM(y1,n,parms,nparms)
		denom=dot_vect(1,gradf,gradf,6)
	factor=one
	f1=f0
ccc	iter=0
ccc	do while(f1.ge.f0)
ccc		iter=iter+1
ccc		if (iter.gt.itermax) then
ccc		error=10
ccc		endif
		do i=1,6
			y1(i)=y0(i)-f0*gradf(i)/denom/factor
		end do
		do i=13,n
			y1(i)=y0(i)
		end do
		pp1=(y1(1)+y1(2)+y1(3))*onethird
ccc		factor=factor*2
			if(pp1.lt.zero)then
c			write(*,*)'pp1<0=',pp1
				switch2=1
				return
			endif

		f1=yf_DM(y1,n,parms,nparms)
ccc	enddo


c		write(*,*) 'drift_corr: normal correction'
		call push(y1,y0,n)

c		pause
	  end if
c
c ... recompute drift at the new state
c
	  f0=yf_DM(y0,n,parms,nparms)
	  p=(y0(1)+y0(2)+y0(3))*onethird
	  f0_p=f0/p
c	  if(p.lt.(p_atm/100)) f0_p=f0_p/1000000
c	  f0_p=f0	  
c
        if(n_drift.gt.max_ndrift) then
          write(6,*) 'ERROR: subroutine DRIFT_CORR:'
          write(6,*) 'too many iterations, increase tolerance'
          write(6,*) 'n_drift = ',n_drift
          write(6,*) 'drift = ',f0_p
c          error=10
c	   error=3
c          return
	  f0_p=0
        end if		

	
c ... bottom of while loop
c
c	  write(*,*) 'drift_corr:      f0 = ',f0
c	  write(*,*) 'drift_corr: n_drift = ',n_drift
	end do
c
c ... return corrected stress and q into vector y
c
      call push(y0,y,n)      
c	
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine el_stiff_DM(y,n,parms,nparms,DDe,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c------------------------------------------------------------------------------
c subroutine to compute elastic stiffness 
c Dafalias &Manzari SANISAND model (2004)
c
c material parameters
c
c  1     p_a        Atmospheric pressure 
c  2     e0         Void ratio on CSL at p = 0  
c  3     lambda     CSL parameter (e:p plane)
c  4     xi         CSL parameter (e:p plane)
c  5     M_c        Slope of CSL in q:p plane, TX compression
c  6     M_e        Slope of CSL in q:p plane, TX extension
c  7     mm         opening of yield surface cone
c  8     G0         Shear modulus constant
c  9     nu         Poisson's ratio
c  10    h0         Plastic modulus constant
c  11    c_h        Plastic modulus constant
c  12    n_b        Plastic modulus constant
c  13    A0         Dilatancy constant
c  14    n_d        Dilatancy constant
c  15    z_max      Fabric index constant
c  16    c_z        Fabric index constant
c  17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  NOTE: soil mechanics convention (compression positive)
c        all stress and strain vectors are 6-dimensional
c------------------------------------------------------------------------------
      implicit none
c
      integer i,j,n,nparms
c
      double precision y(n),parms(nparms)
      double precision p_a,G0,nu,ratio
      double precision sig1,sig2,sig3,p,void
      double precision coeff1,coeff2
      double precision Kt,Gt,fe
      double precision Id(6,6),IxI(6,6),DDe(6,6)
      double precision zero,half,one,two,three 
      double precision pp,p_thres_E,tenm3
c
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

      parameter(zero=1.0d0,half=0.5d0)
      parameter(one=1.0d0,two=2.0d0,three=3.0d0)
	parameter(p_thres_E=0.001d0)
cc      parameter(tenm3=0.001d0)
c
c	common /z_threshold_pressure/p_thres
c
c ... initialize matrices
c
      call pzero(Id,36)
      call pzero(IxI,36)
      call pzero(DDe,36)
c
      Id(1,1)=one
      Id(2,2)=one
      Id(3,3)=one
      Id(4,4)=half
      Id(5,5)=half
      Id(6,6)=half
c
      IxI(1,1)=one
      IxI(2,1)=one
      IxI(3,1)=one
      IxI(1,2)=one
      IxI(2,2)=one
      IxI(3,2)=one
      IxI(1,3)=one
      IxI(2,3)=one
      IxI(3,3)=one
c
c ... recover material parameters
c
      p_a=parms(1)
      G0=parms(8)
      nu=parms(9) 
c
c ... recover state variables
c
      sig1=y(1)
      sig2=y(2)
      sig3=y(3)
c
      void=y(13)
c
c ... mean stress
c
      p=(sig1+sig2+sig3)/three
c
	pp=p
ccc	if(p.lt.p_thres_E)then
ccc		pp=p_thres_E
ccc	end if
	if(p.lt.p_thres)then
		pp=p_thres
	end if
c
c ... max. shear modulus, tangent shear modulus Gt, tangent bulk modulus Kt
c
      ratio=three*(one-two*nu)/(two*(one+nu))
      fe=(2.97d0-void)*(2.97d0-void)/(one+void)
      Gt=G0*p_a*fe*dsqrt(pp/p_a)
      Kt=Gt/ratio
c
c ... elastic stiffness, stored in matrix DDe(6,6)
c
      coeff1=Kt-two*Gt/three
      coeff2=two*Gt
c
      do i=1,6
        do j=1,6
          DDe(i,j)=coeff1*IxI(i,j)+coeff2*Id(i,j)
        end do
      end do
c
      return
      end 
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
	  subroutine f_hypoelas_DM(y,n,parms,nparms,deps,F,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c
c ... computes the function F(y) for (hypo)elastic processes
c     Dafalias & Manzari SANISAND Model (2004)
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c  written by Tamagnini 10/2008
c
c-----------------------------------------------------------------------------
c
      implicit none
c
      external matmul
c
      integer n,m,nparms
c
      double precision y(n),parms(nparms),deps(6)
      double precision depsv,void
      double precision F(n),De(6,6),dsig_e(6)
      double precision one
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      data one/1.0d0/
c
      call pzero(F,n)
c
c ... void ratio and volum. strain increment
c
      void = y(13)
      depsv=deps(1)+deps(2)+deps(3)
c
c ... elastic stiffness matrix
c
      call el_stiff_DM(y,n,parms,nparms,De,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
      call matmul(De,deps,dsig_e,6,6,1)
c
      F(1)=dsig_e(1)
      F(2)=dsig_e(2)
      F(3)=dsig_e(3)
      F(4)=dsig_e(4)
      F(5)=dsig_e(5)
      F(6)=dsig_e(6)
c
      F(13)=-(one+void)*depsv
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine f_plas_DM(y,n,nasvy,z,nz,parms,nparms,deps,kRK,nfev,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c calculate coefficient kRK from current state (stored in y and z) and
c strain increment deps
c 
c Dafalias & Manzari (2004) SANISAND model for sands
c
c variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c
      integer n,nz,nasvy,nparms,i,nfev
c
	integer switch2,mario_DT_test
c
      double precision y(n),z(nz),kRK(n),parms(nparms),deps(6)
      double precision F_sig(6),F_q(nasvy)
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
	double precision zero
c
	parameter(zero=0.0d0)
c
c      common /z_nct_errcode/error	
c
c ... update counter for the number of function f(y) evaluations
c
      nfev=nfev+1
c
c ... initialize kRK
c
      call pzero(kRK,n)
c	
c ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
c
      call get_F_sig_q(y,n,nasvy,z,nz,parms,nparms,deps,F_sig,F_q,
     & switch2,mario_DT_test,error)
		
	if(switch2.gt.zero) return

      if(error.eq.10) return
c
      do i=1,6
        kRK(i)=F_sig(i)
      end do			 
c	
      do i=1,nasvy
        kRK(6+i)=F_q(i)
      end do			 
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine get_F_sig_q(y,n,nasvy,z,nz,parms,nparms,deps,F_sig,F_q,
     & switch2,mario_DT_test,error)
c-----------------------------------------------------------------------------
c
c  computes vectors F_sigma and F_q in F(y)
c  Dafalias & Manzari (2004) SANISAND model for sands
c
c  variables allocated in vector y(n):
c
c variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
      external matmul
c		
	integer switch2,mario_DT_test
c 
      integer nparms,n,nasvy,nz
c
      double precision y(n),z(nz),parms(nparms),deps(6)
      double precision Dep(6,6),HH(nasvy,6),F_sig(6),F_q(nasvy)
      double precision zero
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

	parameter(zero=0.0d0)
      p_thres = 0.001d0
c
c ... compute tangent operators
c
      call get_tan_DM(y,n,nasvy,z,nz,parms,nparms,Dep,HH,switch2,
     & mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
	if(switch2.gt.zero) then
c		write(*,*) 'get_tan - switch2>0'
		return
	endif
c
c ... compute F_sig=Dep*deps
c
      call matmul(Dep,deps,F_sig,6,6,1)
c
c ... compute F_q=HH*deps
c
      call matmul(HH,deps,F_q,nasvy,6,1)
c	
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine get_tan_DM(y,ny,nasvy,z,nz,parms,nparms,Dep,Hep,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c computes matrices Dep and Hep
c Dafalias & Manzari (2004) SANISAND model for sands
c
c variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c material properties allocated in vector parms(nparms):
c
c   1     p_a        Atmospheric pressure 
c   2     e0         Void ratio on CSL at p = 0  
c   3     lambda     CSL parameter (e:p plane)
c   4     xi         CSL parameter (e:p plane)
c   5     M_c        Slope of CSL in q:p plane, TX compression
c   6     M_e        Slope of CSL in q:p plane, TX extension
c   7     mm         opening of yield surface cone
c   8     G0         Shear modulus constant
c   9     nu         Poisson's ratio
c   10    h0         Plastic modulus constant
c   11    c_h        Plastic modulus constant
c   12    n_b        Plastic modulus constant
c   13    A0         Dilatancy constant
c   14    n_d        Dilatancy constant
c   15    z_max      Fabric index constant
c   16    c_z        Fabric index constant
c   17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c  NOTE: stress and strain convention: compression positive
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
      external matmul
c 
      integer nparms,ny,nz,nasvy,i,j,switch2,iter,iter_max,switch4
	integer mario_DT_test
c
      double precision dot_vect,distance,psi_void_DM
c
      double precision y(ny),z(nz),parms(nparms)
      double precision De(6,6),Dep(6,6),Hep(nasvy,6),m(6)
      double precision LL(6),LL1(6),RR(6),RR1(6),U(6),V(6)
c
      double precision p_a,e0,lambda,xi,M_c,M_e,cM,mm
      double precision G0,nu,h0,c_h,n_b
      double precision A0,n_d,z_max,c_z,bulk_w
c
      double precision sig(6),alpha(6),void,Fab(6)
      double precision alpha_sr(6),alpha_b(6)
      double precision s(6),tau(6),n(6)
	double precision norm2,norm,I1,p,psi,cos3t,gth,dgdth
	double precision b0,d_sr,hh,db
      double precision Hplas,LDeR,Kp,Kpm1
      double precision mtrR,brack_mtrR,tol_ff,tol_dil,Hvs
      double precision h_alpha(6),h_fab(6),HH_alpha(6,6),HH_fab(6,6)
	double precision yf_DM,ff0,chvoid
c
      double precision zero,tiny,half,one,two,three,large,kappa
      double precision onethird,twothird
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(tiny=1.0d-15,large=1.0e15)
c Heaviside function parameter
	parameter(kappa=3.0d2)
c
c      common /z_nct_errcode/error	
c
      data m/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
c
	switch2=zero
	switch4=zero
	iter=0
	iter_max=1e3
c
c ... initialize constants and vectors
c
      onethird=one/three
      twothird=two/three
      half=one/two
c
      call pzero(Dep,36)
      call pzero(Hep,6*nasvy)
c
c ... recover material parameters
c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
      M_c=parms(5)
      M_e=parms(6)
      mm=parms(7)
      G0=parms(8)
      nu=parms(9) 
      h0=parms(10)
      c_h=parms(11)
      n_b=parms(12)
      A0=parms(13)
      n_d=parms(14)
      z_max=parms(15)
      c_z=parms(16)
      bulk_w=parms(17)
c	  
	cM=M_e/M_c
c
c ... recover state variables
c
      do i=1,6
        sig(i)=y(i)	
      end do !i
c
      do i=1,6
        alpha(i)=y(6+i)	
      end do !i
c
      void=y(13)
c
      do i=1,6
        Fab(i)=y(13+i)	
      end do !i
c
      do i=1,6
        alpha_sr(i)=z(i)	
      end do !i	  
c
c ... deviator stress and mean pressure
c
      call deviator(sig,s,I1,p)
c
c ... stress ratio tensor and unit vector n
c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
      if(norm.lt.tiny) then
        norm=tiny
      endif
      do i=1,6
        n(i)=tau(i)/norm
c		if(n(i).lt.tiny)then
c		n(i)=zero
c		endif
      end do
c
c ... elastic stiffness
c
      call el_stiff_DM(y,ny,parms,nparms,De,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
c ... gradient of yield function L and tensor V=De*L
c
      call grad_f_DM(y,ny,parms,nparms,LL,LL1)
      call matmul(De,LL1,V,6,6,1)
c
c ... gradient of plastic potential R and tensor U=De*R
c
      call grad_g_DM(y,ny,parms,nparms,RR,RR1)
      call matmul(De,RR1,U,6,6,1)
c
c ... plastic modulus functions b0 and hh
c	  
	  if (dabs(p).gt.zero) then
	    chvoid=c_h*void
	    if(chvoid.ge.1) then
c	     error=3
             chvoid=0.99999
	    end if
	    b0=G0*h0*(one-chvoid)/dsqrt(p/p_a)
	  else
	    b0=large
	  end if
c
      d_sr=distance(alpha,alpha_sr,n)
	  if (d_sr.lt.zero) then
		call push(alpha,alpha_sr,6)
c		write(*,*) 'alpha updated'
	  end if
c
	  if (d_sr.lt.tiny) then
		d_sr=tiny
	  end if
c
	hh=b0/d_sr

c
c ... hardening function h_alpha	  
c
      psi=psi_void_DM(void,p,parms,nparms)
      call lode_DM(tau,cM,cos3t,gth,dgdth)
      call alpha_th_DM(2,n,gth,psi,parms,nparms,alpha_b)
      db=distance(alpha_b,alpha,n)
c
      do i=1,6
        h_alpha(i)=twothird*hh*(alpha_b(i)-alpha(i))
      end do
c
c ... hardening function h_fab
c
 	mtrR=-RR(1)-RR(2)-RR(3)
      brack_mtrR=half*(mtrR+dabs(mtrR))
c
c chiara heaviside function
c
cc	tol_dil=tol_ff*A0
c
c	if((-tol_dil.lt.mtrR).and.(tol_dil.gt.mtrR)) then
c		Hvs=one/(1+exp(two*mtrR*kappa))
c		bracK_mtrR=Hvs*mtrR
c	endif
c
      do i=1,6
	  h_fab(i)=-c_z*brack_mtrR*(z_max*n(i)+Fab(i))
      end do
c
c ... plastic moduli Hplas and Kp
c
      Hplas=twothird*hh*p*db

	if(Hplas.gt.1e+15) then
c		write(*,*)'Hplas'

	endif
c
      LDeR=dot_vect(1,LL1,U,6)
c
      Kp=LDeR+Hplas

	ff0=yf_DM(y,ny,parms,nparms)
c .......................................................................
	
	if(mario_DT_test.eq.zero) then

		if(LDeR.lt.zero) then
		switch2=1
c		write(*,*)'LDeR < 0'
		return
		endif

	

		if(Kp.lt.zero) then
		switch2=1
c		write(*,*)'function get_tan: Kp < zero'
		return
		endif
	
	else
		

		if(LDeR.le.zero) then
		switch2=1
c		error=3
c		write(*,*)'subroutine get_tan_DM: LDeR < 0'
		return
		endif
	endif



	if(Kp.lt.zero)then
c		write(6,*)'function get_tan: Kp < 0'
	error=3
	return
	endif
	call push(alpha_sr,z,6)
      Kpm1=one/Kp
c	if (Kpm1.lt.tiny) then
c		Kpm1=zero
c	  end if
c
c ... elastoplastic stiffness matrix
c
      do i=1,6
        do j=1,6
          Dep(i,j)=De(i,j)-Kpm1*U(i)*V(j)
        end do !j
      end do !i
c
c ... hardening tensor H_alpha
      do i=1,6
        do j=1,6
          HH_alpha(i,j)=Kpm1*h_alpha(i)*V(j)
        end do !j
      end do !i
c
c ... hardening tensor H_fab
c		
      do i=1,6
        do j=1,6
          HH_fab(i,j)=Kpm1*h_fab(i)*V(j)
        end do !j
      end do !i
c
c ... Build tangent evolution matrix Hep(nasv,6) row-wise
c
      do j=1,6
c
        Hep(1,j) =HH_alpha(1,j)			! alpha(1)
        Hep(2,j) =HH_alpha(2,j)			! alpha(2)
        Hep(3,j) =HH_alpha(3,j)			! alpha(3)
        Hep(4,j) =HH_alpha(4,j)			! alpha(4)
        Hep(5,j) =HH_alpha(5,j)			! alpha(5)
        Hep(6,j) =HH_alpha(6,j)			! alpha(6)
        Hep(7,j) =-(one+void)*m(j)		! void
        Hep(8,j) =HH_fab(1,j)			! Fab(1)
        Hep(9,j) =HH_fab(2,j)			! Fab(2)
        Hep(10,j)=HH_fab(3,j)			! Fab(3)
        Hep(11,j)=HH_fab(4,j)			! Fab(4)
        Hep(12,j)=HH_fab(5,j)			! Fab(5)
        Hep(13,j)=HH_fab(6,j)			! Fab(6)
c
      end do !j
c	
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine grad_f_DM(y,ny,parms,nparms,gradf,gradf1)
c------------------------------------------------------------------------------
c ... subroutine to compute the gradient of yield function at state y 
c     Dafalias & Manzari (2004) SNAISAND model for sands
c
c     P(i)  = grad(f) stress-like vector in Voigt notation
c     P1(i) = grad(f) strain-like vector in Voigt notation
c
c ... variables allocated in parms
c
c  1     p_a        Atmospheric pressure 
c  2     e0         Void ratio on CSL at p = 0  
c  3     lambda     CSL parameter (e:p plane)
c  4     xi         CSL parameter (e:p plane)
c  5     M_c        Slope of CSL in q:p plane, TX compression
c  6     M_e        Slope of CSL in q:p plane, TX extension
c  7     mm         opening of yield surface cone
c  8     G0         Shear modulus constant
c  9     nu         Poisson's ratio
c  10    h0         Plastic modulus constant
c  11    c_h        Plastic modulus constant
c  12    n_b        Plastic modulus constant
c  13    A0         Dilatancy constant
c  14    n_d        Dilatancy constant
c  15    z_max      Fabric index constant
c  16    c_z        Fabric index constant
c  17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c  variables allocated in vector y(n):
c
c       y(1) = sig(1)
c       y(2) = sig(2)
c       y(3) = sig(3)
c       y(4) = sig(4)
c       y(5) = sig(5)
c       y(6) = sig(6)
c       y(7) = alpha(1)
c       y(8) = alpha(2)
c       y(9) = alpha(3)
c      y(10) = alpha(4)
c      y(11) = alpha(5)
c      y(12) = alpha(6)
c      y(13) = void
c      y(14) = Fab(1)
c      y(15) = Fab(2)
c      y(16) = Fab(3)
c      y(17) = Fab(4)
c      y(18) = Fab(5)
c      y(19) = Fab(6)
c      y(20) = not used
c
c  NOTE: soil mechanics convention (compression positive)
c        all stress and strain vectors are 6-dimensional
c------------------------------------------------------------------------------
      implicit none
c
      double precision dot_vect
c
      integer ny,nparms,i
c
      double precision parms(nparms),y(ny),gradf(6),gradf1(6),del(6)
      double precision mm,sig(6),s(6),r(6),I1,p
      double precision alpha(6),tau(6),n(6)
      double precision norm,norm2,v,vv
      double precision one,two,three,sqrt23,onethird,small
	double precision n1,n2
c
      parameter(one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(small=1.0d-10)
	parameter(n1=0.816496580927739,n2=-0.40824829046385)
c
      data del/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
c
      sqrt23=dsqrt(two/three)
      onethird=one/three
c
      call pzero(n,6)
c
c ... recover material parameters
c
      mm=parms(7)
c
c ... recover state variables
c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
c
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
c
c ... deviator stress and mean pressure
c
      call deviator(sig,s,I1,p)
c
c ... reduced stress tensor and unit vector
c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
c
	if(norm.lt.small) then
	norm=small
	endif
c
      do i=1,6
        n(i)=tau(i)/norm
      enddo
c
c	norm_n=dot_vect(1,n,n,6)
c
c	if(norm2.lt.small) then
c	    n(1)=n1
c		n(2)=n2
c		n(3)=n2
c      endif
c
c
c ... coefficient V
c
      if(dabs(p).lt.small) then
        do i=1,6
          r(i)=s(i)/small
        enddo
      else
        do i=1,6
          r(i)=s(i)/p
        enddo
      endif
      v=dot_vect(1,r,n,6)
      vv=-onethird*v
c
c ... gradient of f
c
      do i=1,6
        gradf(i)=n(i)+vv*del(i)
        if(i.le.3) then
          gradf1(i)=gradf(i)
        else
          gradf1(i)=two*gradf(i)
        endif
      enddo
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c-----------------------------------------------------------------------------
      subroutine grad_g_DM(y,ny,parms,nparms,gradg,gradg1)
c------------------------------------------------------------------------------
c ... subroutine to compute the gradient of plastic potential at state y 
c     Dafalias & Manzari (2004) SANISAND model for sands
c
c     gradg(i)  = grad(g) stress-like vector in Voigt notation
c     gradg1(i) = grad(g) strain-like vector in Voigt notation
c
c ... variables allocated in parms
c
c  1     p_a        Atmospheric pressure 
c  2     e0         Void ratio on CSL at p = 0  
c  3     lambda     CSL parameter (e:p plane)
c  4     xi         CSL parameter (e:p plane)
c  5     M_c        Slope of CSL in q:p plane, TX compression
c  6     M_e        Slope of CSL in q:p plane, TX extension
c  7     mm         opening of yield surface cone
c  8     G0         Shear modulus constant
c  9     nu         Poisson's ratio
c  10    h0         Plastic modulus constant
c  11    c_h        Plastic modulus constant
c  12    n_b        Plastic modulus constant
c  13    A0         Dilatancy constant
c  14    n_d        Dilatancy constant
c  15    z_max      Fabric index constant
c  16    c_z        Fabric index constant
c  17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c  variables allocated in vector y(n):
c
c       y(1) = sig(1)
c       y(2) = sig(2)
c       y(3) = sig(3)
c       y(4) = sig(4)
c       y(5) = sig(5)
c       y(6) = sig(6)
c       y(7) = alpha(1)
c       y(8) = alpha(2)
c       y(9) = alpha(3)
c      y(10) = alpha(4)
c      y(11) = alpha(5)
c      y(12) = alpha(6)
c      y(13) = void
c      y(14) = Fab(1)  (stress--like)
c      y(15) = Fab(2)
c      y(16) = Fab(3)
c      y(17) = Fab(4)
c      y(18) = Fab(5)
c      y(19) = Fab(6)
c      y(20) = not used
c
c  NOTE: soil mechanics convention (compression positive)
c        all stress and strain vectors are 6-dimensional
c------------------------------------------------------------------------------
      implicit none
c
      double precision dot_vect,distance,psi_void,psi_void_DM
c
      integer ny,nparms,i
c
      double precision M_c,M_e,cM,A0
c
      double precision parms(nparms),y(ny),gradg(6),gradg1(6)
      double precision sig(6),s(6),alpha(6),Fab(6),I1,p
      double precision n(6),n2(6),tau(6),Rdev(6)
      double precision Ad,alpha_d(6),dd	  
	double precision cos3t,gth,dgdth
      double precision void,psi,dil,dil3
c
      double precision temp1,temp2,temp3,temp4
      double precision norm,norm2
      double precision zero,one,two,three,six
	double precision half,sqrt6,onethird,small,del(6)
c
	integer chiara
c
      parameter(half=0.5d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      parameter(zero=0.0d0,small=1.0d-10)
c
      data del/1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
c
      sqrt6=dsqrt(six)
      onethird=one/three
c
      call pzero(n,6)
c
c ... recover material parameters
c
      M_c=parms(5)
      M_e=parms(6)
      A0=parms(13)
c
      cM=M_e/M_c
c
c ... recover state variables
c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
c
      alpha(1)=y(7)
      alpha(2)=y(8) 
      alpha(3)=y(9) 
      alpha(4)=y(10) 
      alpha(5)=y(11) 
      alpha(6)=y(12)
c
      void=y(13)
c
      Fab(1)=y(14)
      Fab(2)=y(15)
      Fab(3)=y(16)
      Fab(4)=y(17)
      Fab(5)=y(18)
      Fab(6)=y(19)
c
c ... deviator stress and mean pressure
c
      call deviator(sig,s,I1,p)
c
c ... stress ratio tensor, unit tensor n and tensor n^2 
c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
c
      if(norm.lt.small) then
	  norm=small
      endif
      do i=1,6
          n(i)=tau(i)/norm
      enddo
c
      n2(1)=n(1)*n(1)+n(4)*n(4)+n(5)*n(5)
      n2(2)=n(4)*n(4)+n(2)*n(2)+n(6)*n(6)
      n2(3)=n(6)*n(6)+n(5)*n(5)+n(3)*n(3)
      n2(4)=n(1)*n(4)+n(4)*n(2)+n(6)*n(5)
      n2(5)=n(5)*n(1)+n(6)*n(4)+n(3)*n(5)
      n2(6)=n(4)*n(5)+n(2)*n(6)+n(6)*n(3)
c
c ... state parameter psi
c
      psi=psi_void_DM(void,p,parms,nparms)
c
c ... Lode angle; functions g(theta) and (1/g)dg/dtheta
c
      call lode_DM(tau,cM,cos3t,gth,dgdth)
c
c ... vector Rdev
c	C&M
      temp1=one+three*cos3t*dgdth
      temp2=-three*sqrt6*dgdth
      do i=1,6
          Rdev(i)=temp1*n(i)+temp2*(n2(i)-onethird*del(i))
      enddo
c	  	  
c ... dilatancy function
c
	temp3=dot_vect(1,Fab,n,6)
	temp4=half*(temp3+dabs(temp3))
      Ad=A0*(one+temp4)	  
c	  
      call alpha_th_DM(3,n,gth,psi,parms,nparms,alpha_d)
      dd = distance(alpha_d,alpha,n)
c Chiara
	if((psi.gt.zero).and.(dd.lt.zero)) then
		dd=zero
	endif
c end Chiara
c
      dil=Ad*dd
      dil3=onethird*dil
c
      do i=1,6
        gradg(i)=Rdev(i)+dil3*del(i)
        if(i.le.3) then
          gradg1(i)=gradg(i)
        else
          gradg1(i)=two*gradg(i)
        endif
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine iniyz(y,nydim,z,nzdim,qq1,nasvy,qq2,nasvz,sig,ntens)
c-----------------------------------------------------------------------------
c  initializes the vectors of state variables
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c-----------------------------------------------------------------------------
      implicit none
c
      integer i,nydim,nzdim,nasvy,nasvz,ntens
c
      double precision y(nydim),z(nzdim)
      double precision qq1(nasvy),qq2(nasvz),sig(ntens)
c
      call pzero(y,nydim)
      call pzero(z,nzdim)
c
      do i=1,ntens
        y(i) = sig(i)
      enddo
c
      do i=1,nasvy
        y(6+i) = qq1(i)
      enddo
c
      do i=1,nasvz
        z(i) = qq2(i)
      enddo
c
      return
      end
c-----------------------------------------------------------------------------
	subroutine intersect_DM(y0,y1,y_star,n,parms,nparms,tol_ff,
     &    xi,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c
c ... finds the intersection point between the stress path 
c     and the yield surface (Papadimitriou & Bouckovalas model)
c     using Newton method
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c ... variables allocated in parms
c
c   1     p_a        Atmospheric pressure 
c   2     e0         Void ratio on CSL at p = 0  
c   3     lambda     CSL parameter (e:p plane)
c   4     xi         CSL parameter (e:p plane)
c   5     M_c        Slope of CSL in q:p plane, TX compression
c   6     M_e        Slope of CSL in q:p plane, TX extension
c   7     mm         opening of yield surface cone
c   8     G0         Shear modulus constant
c   9     nu         Poisson's ratio
c   10    h0         Plastic modulus constant
c   11    c_h        Plastic modulus constant
c   12    n_b        Plastic modulus constant
c   13    A0         Dilatancy constant
c   14    n_d        Dilatancy constant
c	15    z_max      Fabric index constant
c   16    c_z        Fabric index constant
c   17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c  Tamagnini 10/2008
c
c-----------------------------------------------------------------------------
c
      implicit none
c
      integer n,nparms,maxiter,kiter,i,kiter_bis,bisect
c
      double precision yf_DM,dot_vect
c
      double precision parms(nparms),y0(n),y1(n),y_star(n),y05(n)
      double precision tol_ff,fy_star,err,dfdxi,dfdxi_m1,xi,fy05
	double precision dxi, xip1
      double precision sig0(6),sig1(6),dsig(6),P_star(6),P1_star(6)
      double precision zero,one,half,three,onethird
	double precision pp_star,low, fy11, fy00, xi_max, xi_i, pp05
	double precision y00(n),y11(n)
	
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0,three=3.0d0)
	parameter(low=1.0d-10)
c
c      common /z_nct_errcode/error	
c
      xi=one
      maxiter=5000
      kiter=0
	bisect=0
	kiter_bis=0
c
      do i=1,6
        sig0(i)=y0(i)
        sig1(i)=y1(i)
        dsig(i)=sig1(i)-sig0(i)
      end do !i
c
      call push(y1,y_star,n)
c
      fy_star=yf_DM(y_star,n,parms,nparms)
	onethird=one/three
	pp_star=(y_star(1)+y_star(2)+y_star(3))*onethird
      err=dabs(fy_star/pp_star)
	if(pp_star.gt.one) err=dabs(fy_star)	
c	
c
c	
	if(bisect.eq.0) then
c
c
c ... start Newton iteration
c
      do while ((err.gt.tol_ff).and.(bisect.eq.0))
c
        kiter=kiter+1
c
        call grad_f_DM(y_star,n,parms,nparms,P_star,P1_star)
c
        dfdxi=dot_vect(1,P_star,dsig,6)
		if (dfdxi.lt.low) then
ccc			dfdxi=low
ccc			write(*,*)'dfdxi.lt.zero'
			bisect=1
		endif
        dfdxi_m1=one/dfdxi
c
c ... search direction
c
        dxi=-dfdxi_m1*fy_star
        xip1=xi+dxi
c
c	Line search
c
	  do while ((xip1.lt.zero).or.(xip1.gt.one))
		dxi=half*dxi
	    xip1=xi+dxi
	  end do
c
	  xi=xip1
c
c	End Line search
c
c ... find new intersection point and yield function value
c
        do i=1,n
          y_star(i)=y0(i)+xi*(y1(i)-y0(i))
        end do !i
c
        fy_star=yf_DM(y_star,n,parms,nparms)
		if (fy_star.lt.zero) then
c			write(*,*)'fy_star.lt.zero'
			bisect=1
		else
		onethird=one/three
		pp_star=(y_star(1)+y_star(2)+y_star(3))*onethird
		err=dabs(fy_star/pp_star)
		if(pp_star.gt.one) err=dabs(fy_star)	
cccc		err=dabs(fy_star)
		endif
c
	if (kiter.gt.maxiter+1) then
          write(6,*) 'ERROR: max no. of iterations exceeded'
          write(6,*) 'Subroutine INTERSECT_DM'
          write(6,*) 'err = ',err
c          error=10
c	   error=3
c          return 
	err=0
        end if

	
      end do ! bottom of Newton iteration
c
c ... check that 0 < xi < 1 (intersection point between initial and final states)
c
      if((xi.lt.zero).and.(xi.gt.one)) then 
c
        write(6,*) 'ERROR: the intersection point found lies'
        write(6,*) '       outside the line connecting initial'
        write(6,*) '       and trial stress states'
        write(6,*) 'Subroutine INTERSECT_DM'
        write(6,*) 'xi = ',xi
        xi = zero
c	    error=10
c	    error=3
        return 
c	
	endif
	endif

c
	if(bisect.eq.1) then
c ... start bisection method
ccc	write(*,*) 'bisection method'
c
c		find f((a+b)/2)
	
	do i=1,n
		y00(i)=y0(i)
		y11(i)=y1(i)
	enddo
	fy00 =yf_DM(y00,n,parms,nparms)
	fy11 =yf_DM(y11,n,parms,nparms)
	do i=1,n
		y05(i)=y0(i)
	enddo
	pp05=(y05(1)+y05(2)+y05(3))*onethird	
	fy05 =yf_DM(y05,n,parms,nparms)
cccc	err = dabs(fy05)
	err=abs(fy05/pp05)
	if(pp05.gt.one) err=dabs(fy05)	

	do while(err.gt.tol_ff)
		kiter_bis=kiter_bis+1
c
		do i=1,6
			y05(i)=half*(y00(i)+y11(i))
		enddo	
	
		fy05 =yf_DM(y05,n,parms,nparms)
		pp05=(y05(1)+y05(2)+y05(3))*onethird
cccc		err = dabs(fy05)
		err=abs(fy05/pp05)
		if(pp05.gt.one) err=dabs(fy05)	
	
		if(fy05.lt.zero) then
			call push(y05,y00,n)
		else
			call push(y05,y11,n)
		endif
		
		if (kiter_bis.gt.maxiter+1) then
			write(6,*) 'ERROR: max no. of iterations exceeded'
			write(6,*) 'Subroutine INTERSECT_DM - bisection'
			write(6,*) 'err = ',err
			err=0
c          		error=10
c	                error=3
c          		return 
	  	endif
	enddo

	do i=1,n
		y_star(i)=y05(i)
	enddo
		
c	xi= (y05(1)-y0(1))/(y1(1)-y0(1))

	xi_max=zero
	do i=1,6
		if((y1(i)-y0(i)).ne.zero) then
		  xi_i= (y05(i)-y0(i))/(y1(i)-y0(i))
		  if(xi_i.gt.xi_max) then
			xi_max = xi_i
		  endif
		endif
	enddo
	xi = xi_max
	
c ... end bisection method	
c
	
	endif

c 
      return
      end
c
c------------------------------------------------------------------------------
      subroutine inv_sig(sig,pp,qq,cos3t)
c------------------------------------------------------------------------------
c calculate invariants of stress tensor
c Dafalias &Manzari SANISAND model (2004)
c
c NOTE: Voigt notation is used with the following index conversion
c
c       11 -> 1
c       22 -> 2
c       33 -> 3
c       12 -> 4
c       13 -> 5
c       23 -> 6
c
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision sig(6),sdev(6),s2(6)
      double precision I1,J2bar,J2bar_sq,J3bar,trs2,trs3
      double precision pp,qq,cos3t,numer,denom
c
      double precision zero,one,two,three
      double precision onethird,half,onept5,sqrt3,tiny
c
      double precision dot_vect
c
      data zero,one,two,three/0.0d0,1.0d0,2.0d0,3.0d0/
      data tiny/1.0d-15/
c
c ... some constants
c
      onethird=one/three
      half=one/two
      onept5=three/two
      sqrt3=dsqrt(three)
c
c ... trace and mean stress
c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
c
c ... deviator stress
c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
c
c ... second invariants
c
      trs2=dot_vect(1,sdev,sdev,6)
      J2bar=half*trs2
      qq=dsqrt(onept5*trs2)
c
c ... components of (sdev_ij)(sdev_jk) (stress-like Voigt vector)
c
      s2(1)=sdev(1)*sdev(1)+sdev(4)*sdev(4)+sdev(5)*sdev(5)
      s2(2)=sdev(4)*sdev(4)+sdev(2)*sdev(2)+sdev(6)*sdev(6)
      s2(3)=sdev(6)*sdev(6)+sdev(5)*sdev(5)+sdev(3)*sdev(3)
      s2(4)=sdev(1)*sdev(4)+sdev(4)*sdev(2)+sdev(6)*sdev(5)
      s2(5)=sdev(5)*sdev(1)+sdev(6)*sdev(4)+sdev(3)*sdev(5)
      s2(6)=sdev(4)*sdev(5)+sdev(2)*sdev(6)+sdev(6)*sdev(3)
c	     
c ... Lode angle
c
      if(trs2.lt.tiny) then 
c
        cos3t=one
c		
      else
c
        trs3=dot_vect(1,sdev,s2,6)
c
        J3bar=onethird*trs3
        J2bar_sq=dsqrt(J2bar)
        numer=three*sqrt3*J3bar
        denom=two*(J2bar_sq**3)
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
          cos3t=cos3t/dabs(cos3t)
        end if
c
      end if 
c
      return
      end
c
c------------------------------------------------------------------------------
      subroutine lode_DM(r,cM,cos3t,gth,dgdth)
c------------------------------------------------------------------------------
c calculate cos(3*theta) from deviatoric stress ratio tensor 
c tau(i) = s(i)-p*alpha(i) stored in vector r(6)
c
c computes functions g(theta) and (1/g)dg/dtheta from:
c a) Argyris function (Argyris=1) (as in the original paper)
c b) Van Eekelen function (Argyris=0) (more appropriate for high friction angles)
c
c       gth   = g(theta)
c       dgdth = (1/g)dg/dtheta
c
c NOTE: Voigt notation is used with the following index conversion (ABAQUS)
c
c       11 -> 1
c       22 -> 2
c       33 -> 3
c       12 -> 4
c       13 -> 5
c       23 -> 6
c       
c SM stress convention: compression positive
c
c------------------------------------------------------------------------------
c
      implicit none
c
      integer Argyris
c
      double precision r(6),r2(6)
      double precision trr2,trr3,J2bar,J3bar,J2bar_sq
      double precision cM,n_VE,n_VEm1,numer,denom,cos3t
c
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
	double precision alpha,beta,gth,dgdth
      double precision one,two,three
      double precision onethird,half,sqrt3,tiny
c
      double precision dot_vect
c
      data one,two,three/1.0d0,2.0d0,3.0d0/
      data tiny,n_VE/1.0d-15,-0.25d0/
	data Argyris/0/
c
c ... some constants
c
      onethird=one/three
      half=one/two
      sqrt3=dsqrt(three)
c
c ... second invariant
c
      trr2=dot_vect(1,r,r,6)
      J2bar=half*trr2
c
c ... components of (r_ij)(r_jk) (stress-like Voigt vector)
c
      r2(1)=r(1)*r(1)+r(4)*r(4)+r(5)*r(5)
      r2(2)=r(4)*r(4)+r(2)*r(2)+r(6)*r(6)
      r2(3)=r(6)*r(6)+r(5)*r(5)+r(3)*r(3)
      r2(4)=r(1)*r(4)+r(4)*r(2)+r(6)*r(5)
      r2(5)=r(5)*r(1)+r(6)*r(4)+r(3)*r(5)
      r2(6)=r(4)*r(5)+r(2)*r(6)+r(6)*r(3)
c	     
c ... Lode angle
c
      if(trr2.lt.tiny) then 
c
        cos3t=one
c		
      else
c
        trr3=dot_vect(1,r,r2,6)
c
        J3bar=onethird*trr3
        J2bar_sq=dsqrt(J2bar)
        numer=three*sqrt3*J3bar
        denom=two*(J2bar_sq**3)
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
          cos3t=cos3t/dabs(cos3t)
        end if
c
      end if 
c	     
c ... g function and its derivative
c
	  if (Argyris.ne.0) then
c
c ... Argyris function
c
	  gth=two*cM/((one+cM)-(one-cM)*cos3t)
	  dgdth=(1-cM)*gth/(two*cM) 
c
        else
c
c ... Van Eekelen function
c
        n_VEm1=one/n_VE
c
        tmp1=one/(two**n_VE)
	  tmp2=cM**n_VEm1
        tmp3=one+tmp2
        tmp4=one-tmp2
c
	  alpha=tmp1*(tmp3**n_VE)
        beta=tmp4/tmp3
c
        tmp5=(one+beta*cos3t)**n_VE
        tmp6=one+beta*cos3t
c
        gth=alpha*tmp5
        dgdth=n_VE*beta/tmp6
c
        end if
c
      return
      end
c
c------------------------------------------------------------------------------
      subroutine matmul(a,b,c,l,m,n)
c------------------------------------------------------------------------------
c matrix multiplication
c------------------------------------------------------------------------------
      implicit none
c
      integer i,j,k,l,m,n
c
      double precision a(l,m),b(m,n),c(l,n)
c
      do i=1,l
        do j=1,n
          c(i,j) = 0.0d0
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine move_eps(dstran,ntens,deps,depsv)
c-----------------------------------------------------------------------------
c Move strain/strain increment stran/dstran into eps/deps, 
c computes volumetric strain/strain increment and switches 
c sign convention from solid to soil mechanics
c
c NOTE: 
c   stran  = strain tensor (extension positive)
c	eps    = strain tensor (compression positive)
c	epsv   = vol. strain (compression positive)
c   dstran = strain increment tensor (extension positive)
c	deps   = strain increment tensor (compression positive)
c	depsv  = vol. strain increment (compression positive)
c
c   eps/deps has always 6 components
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision deps(6),dstran(ntens),depsv
c
      call pzero(deps,6)
c
      do i=1,ntens
        deps(i) = -dstran(i)
      enddo
c
      depsv=deps(1)+deps(2)+deps(3)
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine move_sig(stress,ntens,pore,sig)
c-----------------------------------------------------------------------------
c Computes effective stress from total stress (stress) and pore pressure (pore)
c and switches sign convention from solid to soil mechanics
c
c NOTE: stress = total stress tensor (tension positive)
c       pore   = exc. pore pressure (undrained conds., compression positive)
c       sig    = effective stress (compression positive)
c
c       sig has always 6 components
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision sig(6),stress(ntens),pore
c
      call pzero(sig,6)
c
      do i=1,ntens
        if(i.le.3) then
          sig(i) = -stress(i)-pore
        else
          sig(i) = -stress(i)
        end if
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine norm_res_DM(y_til,y_hat,ny,norm_R)
c-----------------------------------------------------------------------------
c  evaluate norm of residual vector Res=||y_hat-y_til||
c  Dafalias & Manzari(2004) SANISAND model for sand
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
	  implicit none
c 
      integer ny,i
c
      double precision y_til(ny),y_hat(ny)
      double precision err(ny),norm_R2,norm_R
      double precision sig_hat(6),sig_til(6),del_sig(6)
      double precision alpha_hat(6),alpha_til(6),del_alpha(6)
	double precision Fab_hat(6),Fab_til(6),del_Fab(6)
      double precision void_hat,void_til,del_void
      double precision norm_sig2,norm_alpha2,norm_Fab2
      double precision norm_sig,norm_alp,norm_Fab
      double precision dot_vect,zero
c
      parameter(zero=0.0d0)
c
      call pzero(err,ny)
c
c ... recover stress tensor and internal variables
c
      do i=1,6
        sig_hat(i)=y_hat(i)
        sig_til(i)=y_til(i)
        del_sig(i)=dabs(sig_hat(i)-sig_til(i))
      end do
c
      do i=1,6
        alpha_hat(i)=y_hat(6+i)
        alpha_til(i)=y_til(6+i)
        del_alpha(i)=dabs(alpha_hat(i)-alpha_til(i))
      end do
c
      void_hat=y_hat(13)
      void_til=y_til(13)
      del_void=dabs(void_hat-void_til)
c
      do i=1,6
        Fab_hat(i)=y_hat(13+i)
        Fab_til(i)=y_til(13+i)
        del_Fab(i)=dabs(Fab_hat(i)-Fab_til(i))
      end do
c
c ... relative error norms
c
      norm_sig2=dot_vect(1,sig_hat,sig_hat,6)
      norm_alpha2=dot_vect(1,alpha_hat,alpha_hat,6)
	norm_Fab2=dot_vect(1,Fab_hat,Fab_hat,6)
      norm_sig=dsqrt(norm_sig2)
      norm_alp=dsqrt(norm_alpha2)
	norm_Fab=dsqrt(norm_Fab2)
c
      if(norm_sig.gt.zero) then
        do i=1,6
          err(i)=del_sig(i)/norm_sig
        end do
      end if
c
      if(norm_alp.gt.zero) then
        do i=1,6
          err(6+i)=del_alpha(i)/norm_alp
        end do
      end if
c
      err(13)=del_void/void_hat
c
c      if(norm_Fab.gt.zero) then
c        do i=1,6
c          err(13+i)=del_Fab(i)/norm_Fab
c        end do
c      end if
c
c chiara 4 maggio 2010
c
      do i=1,6
		if((Fab_til(i).ne.zero).and.(norm_Fab.gt.zero)) then
			err(13+i)=del_Fab(i)/norm_Fab
		end if
      end do    
c
c ... global relative error norm
c
	norm_R2=dot_vect(3,err,err,ny)
      norm_R=dsqrt(norm_R2)
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine pert_DM(y_n,y_np1,z,n,nasvy,nasvz,err_tol,
     &                   maxnint,DTmin,deps_np1,parms,
     &                   nparms,nfev,elprsw,theta,ntens,DD,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c  compute numerically consistent tangent stiffness 
c  Dafalias & Manzari (2004) SANISAND model for sand
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c 
      integer elprsw
c
      integer ntens,jj,kk
      integer n,nasvy,nasvz,nparms,nfev
      integer maxnint,mario_DT_test
c
      double precision y_n(n),y_np1(n),y_star(n),z(nasvz),parms(nparms)
      double precision err_tol
      double precision theta,DTmin
      double precision deps_np1(6),deps_star(6)
      double precision dsig(6),DD(6,6)
      double precision zero,three
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      parameter(zero=0.0d0,three=3.0d0)
c
c      common /z_nct_errcode/error
c      common /z_plastic_flag/plastic
c
c ... initialize DD and y_star
c 
      call pzero(DD,36)
      call pzero(y_star,n)
c
      if(plastic.eq.0) then
c
c ... elastic process, DD = De (explicit)
c
        call el_stiff_DM(y_np1,n,parms,nparms,DD,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
      else
c
c ... plastic process, DD computed using numerical perturbation
c     loop over strain basis vectors
c
	do jj=1,ntens
        call push(y_n,y_star,n)
        call push(deps_np1,deps_star,ntens)
c
c        do jj=1,ntens
c
c ... perturbed strain increment
c
	      deps_star(jj)=deps_star(jj)+theta
c
c ... perturbed final state, stored in y_star
c
          if(error.ne.10) then
            call rkf23_upd_DM(y_star,z,n,nasvy,nasvz,err_tol,maxnint,
     &           DTmin,deps_star,parms,nparms,nfev,elprsw,
     &		 mario_DT_test,
     &           error,tol_f,check_ff,drcor,p_thres,plastic)
          end if
c
c ... stiffness components in column jj (kk = row index, jj = column index)
c
          do kk=1,ntens
            dsig(kk)=y_star(kk)-y_np1(kk)
            DD(kk,jj)=dsig(kk)/theta
          end do !kk
c
        end do !jj
c
      end if	
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine plast_mod_DM(y,ny,z,nz,parms,nparms,h_alpha,Kpm1,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c  computes vector h_hard and plastic modulus Kp
c  Dafalias & Manzari(2004) SANISAND model for sand
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)
c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c  material properties allocated in vector parms(nparms):
c
c   1     p_a        Atmospheric pressure 
c   2     e0         Void ratio on CSL at p = 0  
c   3     lambda     CSL parameter (e:p plane)
c   4     xi         CSL parameter (e:p plane)
c   5     M_c        Slope of CSL in q:p plane, TX compression
c   6     M_e        Slope of CSL in q:p plane, TX extension
c   7     mm         opening of yield surface cone
c   8     G0         Shear modulus constant
c   9     nu         Poisson's ratio
c   10    h0         Plastic modulus constant
c   11    c_h        Plastic modulus constant
c   12    n_b        Plastic modulus constant
c   13    A0         Dilatancy constant
c   14    n_d        Dilatancy constant
c   15    z_max      Fabric index constant
c   16    c_z        Fabric index constant
c   17    bulk_w     Pore water bulk modulus (undrained conditions)
c
c
c  NOTE: stress and strain convention: compression positive
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
      external matmul
c 
      integer nparms,ny,nz,i,iter,switch2,mario_DT_test
c
      double precision dot_vect,distance,ref_db,psi_void,psi_void_DM
c
      double precision y(ny),z(nz),parms(nparms)
      double precision De(6,6),h_alpha(6)
      double precision LL(6),LL1(6),RR(6),RR1(6),U(6),V(6)
c
      double precision p_a,e0,lambda,xi,M_c,M_e,cM,mm
      double precision G0,nu,h0,c_h,n_b
      double precision A0,n_d,z_max,c_z,bulk_w
c
      double precision sig(6),alpha(6),void,Fab(6)
      double precision alpha_sr(6),alpha_b(6)
      double precision s(6),tau(6),n(6),I1,p,psi,cos3t
	double precision b0,d_sr,hh,db,gth,dgdth
      double precision HHp,LDeR,Kp,Kpm1,norm2,norm,chvoid
c
      double precision zero,tiny,one,two,three,large
      double precision twothird
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(tiny=1.0d-15,large=1.0e15)
c
c      common /z_nct_errcode/error	
c
c ... initialize constants and vectors
c
      twothird=two/three
c
c ... recover material parameters
c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
      M_c=parms(5)
      M_e=parms(6)
      mm=parms(7)
      G0=parms(8)
      nu=parms(9) 
      h0=parms(10)
      c_h=parms(11)
      n_b=parms(12)
      A0=parms(13)
      n_d=parms(14)
      z_max=parms(15)
      c_z=parms(16)
      bulk_w=parms(17)
c	  
	  cM=M_e/M_c
c
	switch2=zero
	iter=0
c
c ... recover state variables
c
      do i=1,6
        sig(i)=y(i)	
      end do !i
c
      do i=1,6
        alpha(i)=y(6+i)	
      end do !i
c
      void=y(13)
c
      do i=1,6
        Fab(i)=y(13+i)	
      end do !i
c
      do i=1,6
        alpha_sr(i)=z(i)	
      end do !i	  
c
c ... deviator stress and mean pressure
c
      call deviator(sig,s,I1,p)
c
c ... stress ratio tensor and unit vector n
c
      do i=1,6
        tau(i)=s(i)-p*alpha(i)
      end do ! i
c
      norm2=dot_vect(1,tau,tau,6)
      norm=dsqrt(norm2)
      if(norm.lt.tiny) then
        norm=tiny
      endif
      do i=1,6
        n(i)=tau(i)/norm
c		if(n(i).lt.tiny)then
c		n(i)=zero
c		endif
      end do
c
c ... elastic stiffness
c
      call el_stiff_DM(y,ny,parms,nparms,De,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
c ... gradient of yield function L and tensor V=De*L
c
      call grad_f_DM(y,ny,parms,nparms,LL,LL1)
      call matmul(De,LL1,V,6,6,1)
c
c ... gradient of plastic potential R and tensor U=De*R
c
      call grad_g_DM(y,ny,parms,nparms,RR,RR1)
      call matmul(De,RR1,U,6,6,1)
c
c ... plastic modulus functions b0 and hh
c	  
	  if (dabs(p).gt.zero) then
	    chvoid=c_h*void
	    if(chvoid.ge.1) then
c	     error=3
             chvoid=0.99999
	    end if
	    b0=G0*h0*(one-chvoid)/dsqrt(p/p_a)
	  else
	    b0=large
	  end if
c
      d_sr=distance(alpha,alpha_sr,n)
	  if (d_sr.lt.zero) then
		call push(alpha,alpha_sr,6)
c		write(*,*) 'alpha updated'
	  end if
c
	  if (d_sr.lt.tiny) then
		d_sr=tiny
	  end if
c
	hh=b0/d_sr

c
c ... hardening function h_alpha	  
c
      psi=psi_void_DM(void,p,parms,nparms)
      call lode_DM(tau,cM,cos3t,gth,dgdth)
      call alpha_th_DM(2,n,gth,psi,parms,nparms,alpha_b)
      db=distance(alpha_b,alpha,n)
c
      do i=1,6
        h_alpha(i)=twothird*hh*(alpha_b(i)-alpha(i))
      end do
c
c ... plastic moduli HHp and Kp
c
      HHp=twothird*hh*p*db

	if(HHp.gt.1e+15) then
c		write(*,*)'Hplas'

	endif
c
      LDeR=dot_vect(1,LL1,U,6)
c
      Kp=LDeR+HHp

c .......................................................................
	
	if(mario_DT_test.eq.zero) then

		if(LDeR.lt.zero) then
		switch2=1
c		write(*,*)'LDeR < 0'
		return
		endif

	

		if(Kp.lt.zero) then
		switch2=1
c		write(*,*)'function plast_mod: Kp < zero'
		return
		endif
	
	else
		

		if(LDeR.le.zero) then
			switch2=1
c		error=3
c		write(*,*)'subroutine plast_mod_DM: LDeR < 0'
		return
		endif
	endif



	if(Kp.lt.zero)then
c		write(6,*)'function plast_mod: Kp < zero'
	error=3
	return
	endif






	call push(alpha_sr,z,6)



      Kpm1=one/Kp
c	
      return
      end
c
c------------------------------------------------------------------------------
      double precision function psi_void_DM(void,p,parms,nparms)
c------------------------------------------------------------------------------
c computes state parameter psi (pyknotropy factor)
c Dafalias & Manzari (2004) SANISAND model for sands
c------------------------------------------------------------------------------
      implicit none
c
      integer nparms
c
      double precision void,p,parms(nparms)
      double precision p_a,e0,lambda,xi,ec
c
c ... recover material parameters
c
      p_a=parms(1)
      e0=parms(2)
      lambda=parms(3)
	xi=parms(4)
c
      ec=e0-lambda*(p/p_a)**xi
      psi_void_DM=void-ec
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine push(a,b,n)
c-----------------------------------------------------------------------------
c push vector a into vector b
c-----------------------------------------------------------------------------
      implicit none
      integer i,n
      double precision a(n),b(n) 
c
      do i=1,n
        b(i)=a(i)
      enddo
c
      return
      end
c
      subroutine pzero(v,nn)
c
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Zero real array of data

c      Inputs:
c         nn     - Length of array

c      Outputs:
c         v(*)   - Array with zero values
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer n,nn
      double precision v(nn)

      do n = 1,nn
        v(n) = 0.0d0
      end do ! n

      end
c
c-----------------------------------------------------------------------------
      subroutine rkf23_upd_DM(y,z,n,nasvy,nasvz,err_tol,maxnint,DTmin,
     &           deps_np1,parms,nparms,nfev,elprsw,
     &		 mario_DT_test,
     &           error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c  Dafalias & Manzari(2004) SANISAND model for sand
c
c  numerical solution of y'=f(y)
c  explicit, adapive RKF23 scheme with local time step extrapolation
c
c  variables allocated in vector y(n):
c
c	y(1)	= sig(1)
c	y(2)	= sig(2)
c	y(3)	= sig(3)
c	y(4)	= sig(4)
c	y(5)	= sig(5)
c	y(6)	= sig(6)
c	y(7)	= alpha(1)
c	y(8)	= alpha(2)
c	y(9)	= alpha(3)
c	y(10)	= alpha(4)
c	y(11)	= alpha(5)
c	y(12)	= alpha(6)
c   y(13)   = void
c   y(14)   = Fab(1)
c   y(15)   = Fab(2)
c   y(16)   = Fab(3)
c   y(17)   = Fab(4)
c   y(18)   = Fab(5)
c   y(19)   = Fab(6)
c   y(20)   = not used
c
c  variables allocated in vector z(nasvz):
c
c	z(1)	= alpha_sr(1)
c	z(2)	= alpha_sr(2)
c	z(3)	= alpha_sr(3)
c	z(4)	= alpha_sr(4)

c	z(5)	= alpha_sr(5)
c	z(6)	= alpha_sr(6)
c	z(7)	= not used
c	z(8)	= not used
c	z(9)	= not used
c	z(10)	= not used
c	z(11)	= not used
c	z(12)	= not used
c	z(13)	= not used
c	z(14)	= not used
c
c-----------------------------------------------------------------------------
      implicit none
c
      integer elprsw,mario,switch2,switch3,mario2
	integer mario_DT, mario_DT_test
c
      integer n,nasvy,nasvz,nparms,i,ksubst,kreject,nfev
      integer maxnint,attempt,maxnint_1
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      double precision y(n),z(nasvz),deps_np1(6)
      double precision parms(nparms),DTmin,err_tol,err_tol_1, err_tol_n
      double precision zero,half,one,two,three,four,six
      double precision ptnine,one6,one3,two3,temp,prod,pt1
	double precision z1(nasvz),deps_np1_star(6), z_k(nasvz)
c
      double precision y_k(n),y_tr(n),y_star(n),yf_DM,y_k1(n)
      double precision y_2(n),y_3(n),y_til(n),y_hat(n)
      double precision p_atm,tol_ff,ff_tr,ff_k
      double precision T_k,DT_k,xi
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision norm_R,S_hull
      double precision Fab(6),dev_fab(5),I1,f_p,absfp2
      double precision pp,onethird,ptone,p_thres2,tol_ff1,pp_k,pp_tr
	double precision ff_k_pp_k,ff_tr_pp_tr,pp_3,pp_2,pp_hat,ten,min_y_tr
	double precision iter, pp_kk
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
      parameter(pt1=1.0d-3,ptone=0.1d0,ten=1.0d1)
c      parameter(drcor=1)
c
c      common /z_nct_errcode/error
c      common /z_plastic_flag/plastic
c	common /z_tolerance/tol_f
c	common /z_check_yield/check_ff
c	common /z_drift_correction/drcor
c	common /z_threshold_pressure/p_thres
c
c
c ... initialize y_k vector and other variables
c
      call pzero(y_k,n)
c
      one6=one/six
      one3=one/three
      two3=two/three
c
      plastic=0
	mario = 0
	mario_DT=0
	mario_DT_test=0

c
	iter=iter+1
ccccc	write(*,*)'iter =', iter
ccccc	if(iter.gt.145) then
ccccc		write(*,*)'stop'
ccccc	endif
c
c ... start of update process
c
      call push(y,y_k,n)
	call push(z,z_k,nasvz)
c
c ... set tolerance for yield function
c
      p_atm=parms(1)
      tol_ff=tol_f*p_atm
c
c ... yield function at current state, T = 0.0
c
c	if (check_ff) then
c
	ff_k=yf_DM(y_k,n,parms,nparms)
	onethird=one/three
	pp_k=(y_k(1)+y_k(2)+y_k(3))*onethird
	ff_k_pp_k=ff_k/pp_k
	if(pp_k.gt.one) ff_k_pp_k=ff_k	

c
c ... abort execution if initial state is outside the YS
c
c	if (check_ff) then
c
	if (ff_k_pp_k.gt.tol_ff) then 
cccc	if (ff_k.gt.tol_ff) then

c		write(6,*) 'ERROR: initial state is outside the YS'
c		write(6,*) 'Subroutine RKF23_UPDATE'
c		write(*,*) 'f = ',ff_k_pp_k
c		write(*,*) 'f = ',ff_k
c 		error=10
c	        error=3
c 		return 
c
                call drift_corr_DM(y_k,n,z1,nasvz,parms,nparms,tol_ff,
     &            switch2,mario_DT_test,
     &            error,tol_f,check_ff,drcor,p_thres,plastic)

	end if
c	endif
c
c ... compute trial solution (single step) and update z
c
	deps_np1_star=deps_np1
	call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
	pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
	if((pp_k.gt.(p_thres+p_thres))) then
		do while(pp_tr.le.p_thres)
			call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
			pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird 
			deps_np1_star=deps_np1_star*half
c			write(*,*) 'pp_k>p_thres'
		end do
	elseif((pp_k.le.(p_thres+p_thres))
     &.and.(pp_tr.gt.(p_thres+p_thres))) then
		deps_np1_star=deps_np1
		call trial_state(y_k,n,parms,nparms,deps_np1_star,y_tr,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
		pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
c		write(*,*) 'pp_k<=p_thres and p_tr>p_thres'
	elseif((pp_k.le.(p_thres+p_thres))
     &.and.(pp_tr.le.(p_thres+p_thres))) then
		 call push(y_k,y_tr,n)
		write(*,*) 'pp_k<=p_thres and pp_tr<=p_thres'
	endif

		
		
c
c ... yield function at trial state
c        
      ff_tr=yf_DM(y_tr,n,parms,nparms)
	pp_tr=(y_tr(1)+y_tr(2)+y_tr(3))*onethird
	ff_tr_pp_tr=ff_tr/pp_tr
	if(pp_tr.gt.one) ff_tr_pp_tr=ff_tr
c
c ... compute scalar product of dsig_tr and grad(f) 
c     to check if crossing of yl occurs
c
      call check_crossing(y_k,y_tr,n,parms,nparms,prod)
c        
c ... check whether plastic loading, elastic unloading or mixed ep loading
c
      if (ff_tr_pp_tr.lt.tol_ff) then
cccc		if (ff_tr.lt.tol_ff) then
c
c ... Case 1: Elastic unloading inside the yield function: 
c             trial state is the final state
c
        call push(y_tr,y_k,n)
c        
      else
c
c ... Case 2: Some plastic loading occurs during the step
c
c ................................................................................................................
c
	if(pp_tr.lt.p_thres) then
c	(situation not admissible)
c
ccc		p_thres2=ptone*p_thres
c
		write(*,*) 'Low mean pressure, p=',pp_k
c
c
c........................................................................................................
c
c
	else
c
        	if ((ff_k_pp_k.lt.(-tol_ff)).or.(prod.lt.zero)) then
cccc			if ((ff_k.lt.(-tol_ff)).or.(prod.lt.zero)) then
c
c ... Case 2a: the initial part of the stress path is inside the YL;
c              find the intersection point and update current state at the YL
		call intersect_DM(y_k,y_tr,y_star,n,parms,nparms,tol_ff,xi,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
          		call push(y_star,y_k,n)
c
c		Chiara 27.08.09 beginning
c
			plastic=1
c
c		Chiara 27.08.09 end    
c        
        	else
c
c Case 2b: all the stress path lies outside the YL
c
        	xi=zero
c
c		Chiara 27.08.09 beginning
c
		plastic=1
c
c		Chiara 27.08.09 end    
c
        	end if
c
c
c initialize normalized time and normalized time step
c
c
        T_k=xi     
        DT_k=(one-xi)
        ksubst=0
        kreject=0
        nfev=0 
	  attempt=1
	  maxnint_1=maxnint
	  err_tol_1=err_tol 
	  err_tol_n=err_tol  
	  switch3=0           
c
c ... start substepping 
c
ccc	mario = zero
c
        do while((T_k.lt.one).and.(mario.eq.zero)
     &.and.(mario_DT.eq.zero)) !**********************************
c
ccc.and.(attempt.ne.3)
          ksubst=ksubst+1
c
c ... write substepping info
c
cccc		if(iter.gt.145)then
c			if(ksubst.gt.82830) then
c          write(*,*) ksubst,T_k,DT_k,pp_hat
c	endif
cccc	endif
c
c ... check for maximum number of substeps
c
c 
          if((ksubst.gt.maxnint_1).or.(switch3.eq.1)) then
c
			if(attempt.eq.1) then
					maxnint_1=2.0*maxnint
					err_tol_1=1000.0*err_tol
		write(*,*) 'I had to increase tolintT,', 'T_k=', T_k
c					call push(y_hat,y_k,n)
c					call push(z1,z,nasvz)
c					call push(y_k,y,n)
     	   				attempt=2
					DT_k=1-T_k
c			return
c				
			elseif(attempt.eq.2) then
ccc					write(*,*) 'number of substeps ',ksubst
ccc          		 	write(*,*) 'is too big, step rejected'
	      				write(*,*) 'attempt number ',attempt
ccc					write(*,*) 'T_k =',T_k		
ccccc					error=3
ccccc					write(*,*) 'mario1 = one'
					mario=one
					call push(z1,z,nasvz)
					call push(y_k,y,n)
ccccc					attempt=3
                    return
			endif
		endif
c
  
				
c
c ... build RK functions
c
		call push(z_k,z1,nasvz)
cc		call push(y_k,y_k1,nasvy)
		pp_kk=(y_k(1)+y_k(2)+y_k(3))*onethird

c
          call f_plas_DM(y_k,n,nasvy,z1,nasvz,parms,nparms,
     +                   deps_np1,kRK_1,nfev,switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
          if(error.eq.10) return

		if (switch2.eq.zero) then
c
c ... find y_2
c
			temp=half*DT_k
c
			do i=1,n
			y_2(i)=y_k(i)+temp*kRK_1(i)
			end do
c
			pp_2=(y_2(1)+y_2(2)+y_2(3))*onethird
c
			if(pp_2.gt.zero)then
cc	if((y_2(1).gt.zero).and.(y_2(2).gt.zero).and.(y_2(3).gt.zero))then
c		
			call f_plas_DM(y_2,n,nasvy,z1,nasvz,parms,nparms,
     +		               deps_np1,kRK_2,nfev,switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
			if(error.eq.10) return

			if (switch2.eq.zero) then
c					
c ... find y_3
c
				do i=1,n
				y_3(i)=y_k(i)-DT_k*kRK_1(i)+two*DT_k*kRK_2(i)
				end do
c
				pp_3=(y_3(1)+y_3(2)+y_3(3))*onethird
c
				if(pp_3.gt.zero)then 
cc		if((y_3(1).gt.zero).and.(y_3(2).gt.zero).and.(y_3(3).gt.zero))then
c
					call f_plas_DM(y_3,n,nasvy,z1,nasvz,parms,nparms,
     +                   deps_np1,kRK_3,nfev,switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
					if(error.eq.10) return
					if (switch2.eq.zero) then
c		 		
c ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
c
						do i=1,n	
						y_til(i)=y_k(i)+DT_k*kRK_2(i)
						y_hat(i)=y_k(i)+DT_k*
     &               (one6*kRK_1(i)+two3*kRK_2(i)+one6*kRK_3(i))
						end do
c
c ... local error estimate
c
						call norm_res_DM(y_til,y_hat,n,norm_R)
c
c ... time step size estimator according to Hull
c	
						S_hull=ptnine*DT_k*(err_tol/norm_R)**one3
c
						if(norm_R.eq.zero) then
ccc							write(*,*)'RKF23: norm_R=0'
ccc							error=3
						endif
c
c
c
c
      if ((norm_R.lt.err_tol).and.(attempt.ne.2).and.
     &      (attempt.ne.3)) then 				
c
c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_2 = one'
cccccc				error =3
				mario=one
			endif
c
            
c ... correct drift from yield surface using Sloan algorithm
c
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,
     &              switch2,mario_DT_test,
     &              error,tol_f,check_ff,drcor,p_thres,plastic)
				end if
				pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif	
c
c
			end if
c
c
		if ((norm_R.lt.err_tol_1).and.(attempt.eq.2)
     &   .and.(switch2.eq.zero)) then 				
c
c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_3 = one'
				mario=one
cccccc				error =3
			endif
c
            
c ... correct drift from yield surface using Sloan algorithm
c
cc			if (switch2.eq.zero) then
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
				end if

			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif !switch2.eq.zero
cc			endif !switch2.eq.zero	
c
c
		end if !(norm_R.lt.err_tol_1).and.(attempt.eq.2).and.(switch2.eq.zero)
c
c
c
			if((norm_R.gt.err_tol).and.(attempt.ne.3)
     & .and.(switch2.eq.zero)) then
c
c ... substep is not accepted, recompute with new (smaller) substep size DT
c
c		  write(*,*) 'substep size ',DT_k,'  norm_R =',norm_R
c
            	DT_k=max(DT_k/four,S_hull)
c
c ... check for minimum step size
c
            		if(DT_k.lt.DTmin) then
             			write(*,*) 'substep size ',DT_k,
     &		         		' is too small, step rejected'
						DT_k= one - T_k
						mario2=1
cc						err_tol_n=100.0*err_tol
						write(*,*) 'err_tol_n=100*err_tol=',err_tol_n
						switch3=1
c						write(*,*) 'mario_4 = one'
c						error =3
            		end if
c					
          	end if
c
		if((norm_R.lt.err_tol_n).and.(attempt.ne.3)
     &	.and.(switch2.eq.zero).and.(mario2.eq.one)) then
c
c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
c
            

			pp_hat=(y_hat(1)+y_hat(2)+y_hat(1))*onethird
			if (pp_hat.lt.p_thres) then
				write(*,*) 'mario_31 = one; pp_hat<p_thres'
				mario=one
cccccc				error=3
			endif
c
            
c ... correct drift from yield surface using Sloan algorithm
c
cc			if (switch2.eq.zero) then
				if(drcor.ne.0) then
	       call drift_corr_DM(y_hat,n,z1,nasvz,parms,nparms,tol_ff,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
				end if

			if (switch2.eq.zero) then
				call push(y_hat,y_k,n)
				call push(z1,z_k,nasvz)

				T_k=T_k+DT_k;
				DT_k=min(four*DT_k,S_hull)
				DT_k=min((one-T_k),DT_k)
			endif !switch2.eq.zero
cc			endif !switch2.eq.zero
		mario2=zero	
cc					
         end if
c
c
			if (attempt.eq.3) then
c
            		write(*,*) 'Third attempt: solution accepted'
		  		
c
c ... correct drift from yield surface using Sloan algorithm
c
		  	if(drcor.ne.0) then
				call drift_corr_DM(y_k,n,z_k,nasvz,parms,nparms,tol_ff,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
	      		endif

			call push(y_hat,y_k,n)
c
            		T_k=T_k+DT_k;
c
          	end if


c
c
c ... if switch2=1 (if switch2=1 occurs in accepted solution part)
			if (switch2.ne.zero) then
				DT_k=DT_k/four
	            		if(DT_k.lt.DTmin) then
						DT_k= one - T_k
ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - switch2>0, mario5 = one'
ccccccccc						error =3
ccccccccc						return
cccccc						attempt=3
            		end if	
			endif
			
	else
c ... if switch2=1 (in kRk_3)
		DT_k=DT_k/four
           		if(DT_k.lt.DTmin) then
							DT_k= one - T_k
ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - in kRk_3, mario6 = one, pp_3=',pp_3
ccccccccc						error =3
ccccccccc						return
cccccc						attempt=3
            		end if	
	endif
		
		else								
cc ... if pp_3 lower than zero
            	DT_k=DT_k/four
				if(DT_k.lt.DTmin) then
					DT_k= one - T_k
ccccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - pp_3<0, mario7=one, pp_2=',pp_2
ccccccccc						error =3
ccccccccc						return
cccccc						attempt=3
            			end if
		endif

	else
c ... if switch2=1 (in kRk_2)
		DT_k=DT_k/four
				if(DT_k.lt.DTmin) then
					DT_k= one - T_k
cccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - kRk_2, mario8=one, pp_2=',pp_2
cccccccc						error =3
cccccccc						return
cccccc						attempt=3
            			end if
	endif


	else
c ... if pp_2 lower than zero
		DT_k=DT_k/four
			if(DT_k.lt.DTmin) then
					DT_k= one - T_k
ccccccc							mario_DT=1
						mario_DT_test=1
c				write(*,*) 'mario_DT - pp_2<0 - mario9=one, pp_kk=',pp_kk
cccccccc						error =3
cccccccc						return
cccccc						attempt=3
            			end if
	endif
c
	else
c ... if switch2=1 (in kRk_1)
		DT_k=DT_k/four
			if(DT_k.lt.DTmin) then
					DT_k= one - T_k
cccccccc						mario_DT=1
						mario_DT_test=1
				write(*,*) 'mario_DT - kRk_1, mario10=one, pp_kk=',pp_kk
c
c
ccccccccc						error =3
ccccccccc							return
cccccc						attempt=3
            			end if
	endif

c
c
c
c ... bottom of while loop
c
        end do !*****************************************************
c
       end if 
      end if
c
	if(mario.eq.1) then !stop solution, keep current configuration
cccccc					DT_k=1-T_k
					write(*,*) 'attempt number ',attempt
      				write(*,*) 'accept solution'
					call push(z_k,z,nasvz)
					call push(y_k,y,n)
c	error=3
ccc	attempt=1
	else if(mario_DT.eq.1) then	 !abort solution, keep previous configuration
					write(*,*) 'mario_DT'
c					call push(z1,z,nasvz)
c					call push(y_k,y,n)
	else		
c
      call push(y_k,y,n)
	call push(z_k,z,nasvz)
c
	endif
c
		if(drcor.ne.0) then
	call drift_corr_DM(y,n,z,nasvz,parms,nparms,tol_ff,
     & switch2,mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
	    endif
		if (switch2.ne.zero) then
			write(*,*) 'RKF23 - drift_corr - switch2>0'
cccccc			error=3
			return
		endif 
c
      return
c
 2000 format('Substep no.',i4,'- T_k=',d12.4,'- DT_k=',d12.4,
     &' -pp_hat=',d12.4)
c
      end
c
!-----------------------------------------------------------------------------
      subroutine solout(stress,ntens,asv1,nasvy,asv2,nasvz,ddsdde,
     +                  y,nydim,z,pore,depsv_np1,parms,nparms,DD)
!-----------------------------------------------------------------------------
! copy the vector of state variables to umat output
! Dafalias &Manzari SANISAND model (2004)
!
! modified 4/2008 (Tamagnini)
!
! NOTE: solid mechanics convention for stress and strain components
!       pore is always positive in compression
!
! depsv_np1 = vol. strain increment, compression positive 
! y(1:6)    = effective stress, compression positive 
! pore      = excess pore pressure, compression positive 
! stress    = total stress, tension positive
!
!-----------------------------------------------------------------------------
      implicit none
!
      integer nydim,nasvy,nasvz,nparms,ntens,i,j
!
      double precision y(nydim),z(nasvz),asv1(nasvy),asv2(nasvz)
      double precision stress(ntens),ddsdde(ntens,ntens),DD(6,6)
      double precision parms(nparms),bulk_w,pore,depsv_np1 
!
      bulk_w=parms(17)
!
! ... update excess pore pressure (if undrained cond.), compression positive
!
      pore=pore+bulk_w*depsv_np1
!
! updated total stresses (effective stresses stored in y(1:6))
!
      do i=1,ntens
		if (i.le.3) then
          stress(i) = -y(i)-pore
		else
          stress(i) = -y(i)
        end if
      enddo
!
! additional state variables
!
      do i=1,nasvy
		asv1(i) = y(6+i)
      enddo
!
      do i=1,nasvz
		asv2(i) = z(i)
      enddo
!
! consistent tangent stiffness
!
      do j=1,ntens
        do i=1,ntens
		if((i.le.3).and.(j.le.3)) then
          ddsdde(i,j) = DD(i,j)+bulk_w
		else
          ddsdde(i,j) = DD(i,j)
		end if        
        end do
      enddo
!
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine tang_stiff(y,z,n,nasvy,nasvz,parms,nparms,DD,cons_lin,
     &           error,tol_f,check_ff,drcor,p_thres,plastic)
c-----------------------------------------------------------------------------
c  compute continuum tangent stiffness at the end of the step
c  Dafalias & Manzari (2004) SANISAND model for sand
c
c  written 10/2008 (Tamagnini)
c-----------------------------------------------------------------------------
      implicit none
c 
      integer switch2,mario_DT_test
c
      integer n,nasvy,nasvz,nparms,cons_lin
c
      double precision y(n),z(nasvz),parms(nparms)
      double precision DD(6,6),HH(nasvy,6)
      double precision zero,three
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres

      parameter(zero=0.0d0,three=3.0d0)
      mario_DT_test = 1
c
c      common /z_nct_errcode/error
c      common /z_plastic_flag/plastic
c
c ... initialize DD
c 
      call pzero(DD,36)
c
      if(plastic.eq.1 .and. cons_lin.eq.1) then
c
c ... plastic process
c
	  call get_tan_DM(y,n,nasvy,z,nasvz,parms,nparms,DD,HH,switch2,
     &    mario_DT_test,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
     
       else       
             call el_stiff_DM(y,n,parms,nparms,DD,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)

c	
c ... if switch2=0 (LDeR<0) try Elastic tangent matrix 
c
! 	if(switch2.ne.zero) then
c        error=3
! 		write(*,*) 'tang_stiff - switch2=1 - elastoplastic_matrix=0'
! 	endif
c
      end if	
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine trial_state(y_k,n,parms,nparms,deps,y_tr,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)               
c-----------------------------------------------------------------------------
c
c ... computes the trial stress state (freezing plastic flow) 
c
c     eps  = strain at the beginning of the step
c     deps = strain increment for this step
c
c	NOTE:
c
c	mode = 1	single step Runge-Kutta 3rd order
c	mode = 2	single step forward Euler
c
c-----------------------------------------------------------------------------
c
      implicit none
c
      integer n,m,nparms,mode,i
c
      double precision y_k(n),parms(nparms),deps(6)
      double precision y_2(n),y_3(n),y_tr(n)
      double precision one,two,three,six
      double precision kRK_1(n),kRK_2(n),kRK_3(n)
      double precision DT_k,DTk05,DTk2,DTk6,DTk23
      
      integer error,check_ff,drcor,plastic
      double precision tol_f,p_thres
c
      parameter(mode=1)
c
      data one,two,three,six/1.0d0,2.0d0,3.0d0,6.0d0/
c
      DT_k=one
c
c ... compute F_el
c
      call f_hypoelas_DM(y_k,n,parms,nparms,deps,kRK_1,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
      if (mode.eq.1) then
c
c ... 3rd order RK - build F function
c
        DTk05=DT_k/two
        DTk2=two*DT_k
        DTk6=DT_k/six
        DTk23=two*DT_k/three
c
        do i=1,n
          y_2(i)=y_k(i)+DTk05*kRK_1(i)
        end do ! i
c
        call f_hypoelas_DM(y_2,n,parms,nparms,deps,kRK_2,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
        do i=1,n
          y_3(i)=y_k(i)-DT_k*kRK_1(i)+DTk2*kRK_2(i)
        end do ! i
c
        call f_hypoelas_DM(y_3,n,parms,nparms,deps,kRK_3,
     &         error,tol_f,check_ff,drcor,p_thres,plastic)
c
        do i=1,n
          y_tr(i)=y_k(i)+DTk6*kRK_1(i)+DTk23*kRK_2(i)+DTk6*kRK_3(i)
        end do ! i
c
      else
c
c ... forward Euler
c
        do i=1,n
          y_tr(i)=y_k(i)+DT_k*kRK_1(i)
        end do ! i
c
      end if
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine wrista(mode,y,nydim,deps_np1,dtime,coords,statev,
     &           nstatv,parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c-----------------------------------------------------------------------------
c ... subroutine for managing output messages
c
c     mode
c
c     all = writes: kstep, kinc, noel, npt
c     2   = writes also: error message,coords(3),parms(nparms),ndi,
c           nshr,stress(nstress),deps(nstress),dtime,statev(nstatv)
c     3   = writes also: stress(nstress),deps(nstress),dtime,statev(nstatv)
c
c-----------------------------------------------------------------------------
      implicit none
c
      integer mode,nydim,nstatv,nparms,noel,npt,ndi,nshr,kstep,kinc,i    
c
      double precision y(nydim),statev(nstatv),parms(nparms)
      double precision deps_np1(6),coords(3),dtime
c
c ... writes for mode = 2
c
      if (mode.eq.2) then
        write(6,*) '==================================================='
        write(6,*) '   ERROR: abaqus job failed during call of UMAT'
        write(6,*) '==================================================='
        write(6,*) ' state dump: '
        write(6,*) 
      endif
c
c ... writes for all mode values
c
	if(mode.ne.4) then
      write(6,*) 'Step: ',kstep, 'increment: ',kinc,
     & 'element: ', noel, 'Integration point: ',npt
      write(6,*) 
	endif
c
c ... writes for mode = 2
c
      if (mode.eq.2) then
        write(6,*) 'Co-ordinates of material point:'
c        write(6,*) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',
c     &    coords(3)
        write(6,*) 
        write(6,*) 'Material parameters:'
        write(6,*) 
        do i=1,nparms
          write(6,105) 'prop(',i,') = ',parms(i)
        enddo 
        write(6,*)
        write(6,102) 'No. of mean components:  ',ndi
        write(6,102) 'No. of shear components: ',nshr
        write(6,*)
      endif
c
c ... writes for mode = 2 or 3
c
      if ((mode.eq.2).or.(mode.eq.3)) then
        write(6,*) 'Stresses:'
        write(6,*) 
        write(6,101) 'sigma(1) = ',y(1)
        write(6,101) 'sigma(2) = ',y(2)
        write(6,101) 'sigma(3) = ',y(3)
        write(6,101) 'sigma(4) = ',y(4)
        write(6,101) 'sigma(5) = ',y(5)
        write(6,101) 'sigma(6) = ',y(6)
        write(6,*) 
        write(6,*) 'Strain increment:'
        write(6,*) 
        write(6,101) 'deps_np1(1) = ',deps_np1(1)
        write(6,101) 'deps_np1(2) = ',deps_np1(2)
        write(6,101) 'deps_np1(3) = ',deps_np1(3)
        write(6,101) 'deps_np1(4) = ',deps_np1(4)
        write(6,101) 'deps_np1(5) = ',deps_np1(5)
        write(6,101) 'deps_np1(6) = ',deps_np1(6)
        write(6,*) 
        write(6,*) 'Time increment:'
        write(6,*) 
        write(6,108) 'dtime = ',dtime
        write(6,*) 
        write(6,*) 'Internal state variables:'
        write(6,*) 
        write(6,109) 'alpha_11     = ',statev(1)
        write(6,109) 'alpha_22     = ',statev(2)
        write(6,109) 'alpha_33     = ',statev(3)
        write(6,109) 'alpha_12     = ',statev(4)
        write(6,109) 'alpha_13     = ',statev(5)
        write(6,109) 'alpha_23     = ',statev(6)
        write(6,109) 'void         = ',statev(7)
        write(6,109) 'Fab_11       = ',statev(8)
        write(6,109) 'Fab_22       = ',statev(9)
        write(6,109) 'Fab_33       = ',statev(10)
        write(6,109) 'Fab_12       = ',statev(11)
        write(6,109) 'Fab_13       = ',statev(12)
        write(6,109) 'Fab_23       = ',statev(13)
        write(6,109) 'dummy        = ',statev(14)
        write(6,109) 'alpha_sr_11  = ',statev(15)
        write(6,109) 'alpha_sr_22  = ',statev(16)
        write(6,109) 'alpha_sr_33  = ',statev(17)
        write(6,109) 'alpha_sr_12  = ',statev(18)
        write(6,109) 'alpha_sr_13  = ',statev(19)
        write(6,109) 'alpha_sr_23  = ',statev(20)
        write(6,109) 'dummy        = ',statev(21)
        write(6,109) 'dummy        = ',statev(22)
        write(6,109) 'dummy        = ',statev(23)
        write(6,109) 'dummy        = ',statev(24)
        write(6,109) 'dummy        = ',statev(25)
        write(6,109) 'dummy        = ',statev(26)
        write(6,109) 'dummy        = ',statev(27)
        write(6,109) 'dummy        = ',statev(28)
        write(6,*) 
        write(6,*) '==================================================='
c
      endif

	if (mode.eq.4) then
	      write(6,*) 'Step: ',kstep, 'increment: ',kinc,
     & 'element: ', noel, 'Integration point: ',npt
c        write(6,*) 'DT<DTmin:x1=',coords(1),' x2 = ',
c     &    coords(2),' x3 = ',coords(3)
c	      write(6,*) 
	endif
c
101   format(1X,a10,f50.44)
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a12,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a10,f50.44)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
c       
      return
      end
c
c-----------------------------------------------------------------------------
      double precision function yf_DM(y,ny,parms,nparms)
c-----------------------------------------------------------------------------
c     compute yield function of Dafalias & Manzari (2004) 
c     SANISAND model for sands
c-----------------------------------------------------------------------------
c
      implicit none
c
      integer ny,nparms
c
      double precision dot_vect
      double precision y(ny),parms(nparms)
      double precision mm,zero,one,two,three,sqrt23,norm2
      double precision sig(6),s(6),trace,p,alpha(6),sbar(6)
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
c
c ... some constants and material parameters
c
      sqrt23=dsqrt(two/three)
      mm=parms(7)
c   	
c ... recover state variables
c
      sig(1)=y(1)
      sig(2)=y(2)
      sig(3)=y(3)
      sig(4)=y(4)
      sig(5)=y(5)
      sig(6)=y(6)
c
      alpha(1)=y(7)
      alpha(2)=y(8)
      alpha(3)=y(9)
      alpha(4)=y(10)
      alpha(5)=y(11)
      alpha(6)=y(12)
c
c ... mean stress and deviator stress
c
      call deviator(sig,s,trace,p)
c
      sbar(1)=s(1)-p*alpha(1)
      sbar(2)=s(2)-p*alpha(2)
      sbar(3)=s(3)-p*alpha(3)
      sbar(4)=s(4)-p*alpha(4)
      sbar(5)=s(5)-p*alpha(5)
      sbar(6)=s(6)-p*alpha(6)
c
c ... compute yield function
c
      norm2=dot_vect(1,sbar,sbar,6)
      yf_DM=dsqrt(norm2)-sqrt23*mm*p
c
      return
      end
c
c.....MODIFICATO
c
	subroutine xit_DM()
	stop
	return
	end
c
      
c------------------------------------------------------------------------------
      subroutine inv_sig_full(sig,pp,qq,cos3t,I1,I2,I3)
c------------------------------------------------------------------------------
c calculate invariants of stress tensor
c
c NOTE: Voigt notation is used with the following index conversion
c
c       11 -> 1
c       22 -> 2
c    33 -> 3
c       12 -> 4
c       13 -> 5
c       23 -> 6
c
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision sig(6),sdev(6)
      double precision eta(6),eta_d(6),eta_d2(6)
      double precision xmin1,xmin2,xmin3
      double precision tretadev3,pp,qq,cos3t,I1,I2,I3
      double precision norm2,norm2sig,norm2eta,numer,denom
c
      double precision half,one,two,three,six
      double precision onethird,threehalves,sqrt6,tiny
c
      double precision dot_vect
c
      data half,one/0.5d0,1.0d0/
      data two,three,six/2.0d0,3.0d0,6.0d0/
      data tiny/1.0d-18/
c
c ... some constants
c
      onethird=one/three
      threehalves=three/two
      sqrt6=dsqrt(six)
c
c ... trace and mean stress
c
      I1=sig(1)+sig(2)+sig(3)
      pp=onethird*I1
c
c ... deviator stress
c
      sdev(1)=sig(1)-pp
      sdev(2)=sig(2)-pp
      sdev(3)=sig(3)-pp
      sdev(4)=sig(4)
      sdev(5)=sig(5)
      sdev(6)=sig(6)
c
c ... normalized stress and dev. normalized stress
c

      if(I1.ne.0) then
         eta(1)=sig(1)/I1
         eta(2)=sig(2)/I1
         eta(3)=sig(3)/I1
         eta(4)=sig(4)/I1
         eta(5)=sig(5)/I1
        eta(6)=sig(6)/I1
      else
        eta(1)=sig(1)/tiny
        eta(2)=sig(2)/tiny
        eta(3)=sig(3)/tiny
        eta(4)=sig(4)/tiny
        eta(5)=sig(5)/tiny
        eta(6)=sig(6)/tiny        
      end if
c
      eta_d(1)=eta(1)-onethird
      eta_d(2)=eta(2)-onethird
      eta_d(3)=eta(3)-onethird
      eta_d(4)=eta(4)
      eta_d(5)=eta(5)
      eta_d(6)=eta(6)
c
c ... second invariants
c
      norm2=dot_vect(1,sdev,sdev,6)
      norm2sig=dot_vect(1,sig,sig,6)
      norm2eta=dot_vect(1,eta_d,eta_d,6)
c
      qq=dsqrt(threehalves*norm2)
      I2=half*(norm2sig-I1*I1)
c
c ... components of (eta_d_ij)(eta_d_jk)
c
      eta_d2(1)=eta_d(1)*eta_d(1)+eta_d(4)*eta_d(4)+eta_d(5)*eta_d(5)
      eta_d2(2)=eta_d(4)*eta_d(4)+eta_d(2)*eta_d(2)+eta_d(6)*eta_d(6)
      eta_d2(3)=eta_d(6)*eta_d(6)+eta_d(5)*eta_d(5)+eta_d(3)*eta_d(3)
      eta_d2(4)=eta_d(1)*eta_d(4)+eta_d(4)*eta_d(2)+eta_d(6)*eta_d(5)
      eta_d2(5)=eta_d(5)*eta_d(1)+eta_d(6)*eta_d(4)+eta_d(3)*eta_d(5)
      eta_d2(6)=eta_d(4)*eta_d(5)+eta_d(2)*eta_d(6)+eta_d(6)*eta_d(3)
c           
c ... Lode angle
c
      if(norm2eta.lt.tiny) then 
c
        cos3t=-one
c               
      else
c
        tretadev3=dot_vect(1,eta_d,eta_d2,6)
c
        numer=-sqrt6*tretadev3
        denom=(dsqrt(norm2eta))**3
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
             cos3t=cos3t/dabs(cos3t)
        end if
c
      end if 
c
c ... determinant
c
      xmin1=sig(2)*sig(3)-sig(6)*sig(6)
      xmin2=sig(4)*sig(3)-sig(6)*sig(5)
      xmin3=sig(4)*sig(6)-sig(5)*sig(2)
c
      I3=sig(1)*xmin1-sig(4)*xmin2+sig(5)*xmin3

c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine check_RKF_DM(error_RKF,y,ny,nasv,parms,nparms)
c-----------------------------------------------------------------------------
c Checks is RKF23 solout vector y is OK for hypoplasticity
c-----------------------------------------------------------------------------
      implicit none
c
        integer error_RKF,ny,nasv,i,nparms,testnan,iopt
c
        double precision y(ny),parms(nparms)
        double precision sig(6),pmean,sig_star(6)
        double precision I1,I2,I3,pp,qq,cos3t
        double precision ptshift,minstress,sin2phim,tolerance
        double precision OCR,omega,fSBS,sensit,cos2phic
        double precision coparam,sin2phicco
c
c	check for NAN
 	testnan=0
        do i=1,ny
       	  call umatisnan_DM(y(i),testnan)
        end do
        if(testnan.eq.1) error_RKF=1
c
      return
      end
c
c-----------------------------------------------------------------------------
      subroutine umatisnan_DM(chcknum,testnan)
c-----------------------------------------------------------------------------
c
c  checks whether number is NaN
c
c-----------------------------------------------------------------------------
        double precision chcknum
        integer testnan

	    if (.not.(chcknum .ge. 0. .OR. chcknum .lt. 0.)) testnan=1        
	    if (chcknum .gt. 1.d30) testnan=1        
	    if (chcknum .lt. -1.d30) testnan=1        
 	    if (chcknum .ne. chcknum) testnan=1        
       
        return
        end         

