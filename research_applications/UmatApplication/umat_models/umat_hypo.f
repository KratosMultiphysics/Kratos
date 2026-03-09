! Copyright (C)  2009  C. Tamagnini, E. Sellari, D. Masin, P.A. von Wolffersdorff
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
!  USA.

c------------------------------------------------------------------------------
      subroutine umat_hypo(stress,statev,ddsdde,sse,spd,scd,
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
c MASIN HYPO - Masin hypoplastic model with intergranular strains
c D. Masin (2005) A hypoplastic constitutive model for clays. IJNAMG, 29:311-336
c
c Implementation based on:
c Fellin, W. and Ostermann, A. (2002):
c Consistent tangent operators for constitutive rate equations.
c International Journal for Numerical and Analytical Methods in Geomechanics
c
c ----------------------------------------------------------------------------
c The string for the material name may contain 9 characters.
c ----------------------------------------------------------------------------
c Material constants:
c       
c     ---------------------------------------------------------------------
c     props(j)      
c     ---------------------------------------------------------------------
c        1      phi 
c        2      p_t 
c        3      h_s
c        4      n
c        5      e_d0
c        6      e_c0
c        7      e_i0
c        8      alpha
c        9      beta
c        10      m_R 
c        11      m_T
c        12      RR
c        13     beta_r 
c        14     chi
c        15     bulk_w
c        16		e0         
c
c     ----------------------------------------------------------------------
c
c Solution dependent state variables (statev):
c definition via sdvini
c
c        1 ... del_11  intergranular strain component
c        2 ... del_22  intergranular strain component
c        3 ... del_33  intergranular strain component
c        4 ... del_12  intergranular strain component
c        5 ... del_13  intergranular strain component
c        6 ... del_23  intergranular strain component
c        7 ... void    void ratio
c        8 ... pore    excess pore pressure (undrained conditions, bulk_w>>0)
c        9 ... p       mean stress (o)
c       10 ... nfev    number of function evaluation
c       11 ... phi_mob phi_mob in degrees
c       12 ... stiff   ratio actual stiffness/hypoplastic stiffness
c       13 ... rho     normalised length of intergr. strain rho
c
c       For undrained analyses with penalty approach:
c       bulk_w: bulk modulus of water phase
c
c Authors: 
c     C. Tamagnini (tamag@unipg.it)
c     E. Sellari
c     Dipartimento di Ingegneria Civile e Ambientale 
c     Universit? degli Studi di Perugia, Italy
c
c NOTE: sign convention for stress and strain: tension and extension positive
c         version with numerical linearization via direct perturbation of dstran
c
c
c Last change: 12/2005  
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
c     5. perturb    = perturbation parameter for numerical computation of Jacobian matrices
c     6. nfasv      = number of first additional state variable in statev field 
c     7. prsw       = switch for printing information
c
c ... declaration of local variables
c
        logical prsw,elprsw
c
      integer i,error,maxnint,nfev,testnan,maxninttest
        integer nparms,nasvdim,nfasv,nydim,nasv,nyact,testing
c
        double precision dot_vect_h
c       
      double precision parms(nprops),theta,tolintT,dtsub,DTmin,perturb
      double precision sig_n(6),sig_np1(6),DDtan(6,6),pore
        double precision deps_np1(6),depsv_np1,norm_deps,tolintTtest
        double precision norm_deps2,pp,qq,cos3t,I1,I2,I3,norm_D2,norm_D
c
      parameter (nasvdim = 13)
      parameter (nydim = 6+nasvdim)
c       parameter (tolintT = 1.0d-3) ...orig value...
        parameter (tolintT = 1.0d-3) 
        parameter (tolintTtest = 1.0d-1) 
c
c       parameter (maxnint = 1000) ...orig value...
        parameter (maxnint = 10000)
        parameter (maxninttest = 1000)
        parameter (DTmin = 1.0d-17)
        parameter (perturb = 1.0d-5)
        parameter (nfasv = 1)
        parameter (prsw = .true.)

c
c ... additional state variables
c
      double precision  asv(nasvdim)
c
c ... solution vector (stresses, additional state variables)
c
      double precision  y(nydim),y_n(nydim),dy(nydim)
c
        common /z_nct_errcode/error
c
c ... Error Management:
c     ----------------
c     error =  0 ... no problem in time integration
c     error =  1 ... problems in evaluation of the time rate, (e.g. undefined 
c                    stress state), reduce time integration substeps
c     error =  3 ... problems in time integration, reduce abaqus load increment 
c                    (cut-back)
c     error = 10 ... severe error, terminate calculation
c
      error=0
c
c ... check problem dimensions
c
      write(6,*) 'in umat_hypo, line 178'
                
      if (ndi.ne.3) then
c
                write(6,*) 'ERROR: this UMAT can be used only for elm.'
                write(6,*) 'with 3 direct stress/strain components'
                write(6,*) 'noel = ',noel
                error=10
c
      endif
c
c ... check material parameters and move them to array parms(nparms)
c
      call check_parms_h(props,nprops,parms,nparms)
      write(6,*) 'in umat_hypo, line 192'
c
c ... print informations about time integration, useful when problems occur
c
      elprsw = .false.
      write(6,*) 'in umat_hypo, line 197'
      if (prsw) then
c
c ... print only in some defined elements
c
      write(6,*) 'in umat_hypo, line 202'

                if ((noel.eq.101).and.(npt.eq.1)) elprsw = .false.
      endif
      write(6,*) 'in umat_hypo, line 204'

c
c ... define number of additional state variables
c
      write(6,*) 'in umat_hypo, line 206'
      call define_h(nasv)
      nyact = 6 + nasv
      if (nyact.gt.nydim) then
          write(6,*) 'ERROR: nasvdim too small in UMAT'
          error=10
      endif
      write(6,*) 'in umat_hypo, line 213'
c
c ... suggested time substep size, and initial excess pore pressure
c
      write(6,*) 'in umat_hypo, line 215'
      dtsub = statev(13)
        pore = -statev(8)
      write(6,*) 'in umat_hypo, line 217'

c
c ... vector of additional state variables
c
      do i=1,nasv
        asv(i) = statev(i-1+nfasv)
      enddo
c
c ... compute volume strain increment and current effective stress tensor
c
       do i=1,6        
            sig_n(i)=0
            deps_np1(i)=0
       end do
        call move_sig_h(stress,ntens,pore,sig_n)
        call move_eps_h(dstran,ntens,deps_np1,depsv_np1)

	  norm_D2=dot_vect_h(2,deps_np1,deps_np1,6)
	  norm_D=sqrt(norm_D2)

c ... check whether the strain rate from the ABAQUS is not NAN	  
	  testnan=0
	  call umatisnan_h(norm_D,testnan)
	  if (testnan .eq. 1) then 
	     call wrista_h(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &              parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
	     write(6,*) 'Error in integration, noel ',noel
	     write(6,*) 'Try to decrease the global step size'
	     call xit_h
	  end if
      write(6,*) 'in umat_hypo, line 249'

c
c --------------------
c ... Time integration
c --------------------
c
        call iniy_h(y,nydim,nasv,ntens,sig_n,asv)
        call push_h(y,y_n,nydim)
c
      if (elprsw) then
        write(6,*) '==================================================='
        write(6,*) 'Call of umat:'
        write(6,*) '==================================================='
        call wrista_h(3,y,nydim,deps_np1,dtime,coords,statev,nstatv,
     &              parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
      endif
c
c ... parameter of the numerical differentiation: sqrt(macheps)*||deps||
c     double precision
c
      norm_deps2=dot_vect_h(2,deps_np1,deps_np1,ntens)
        norm_deps=dsqrt(norm_deps2)
      theta=-perturb*max(norm_deps,1.0d-6)  ! negative sign for compression
c
c ... local integration using adaptive RKF-23 method, consistent Jacobian and error estimation
c
      if((dtsub.le.0.0d0).or.(dtsub.gt.dtime)) then
        dtsub = dtime
      endif
c
c	  if testing==1 PLAXIS is testing for the initial strain increment.
	  testing=0
c	  For use in ABAQUS, comment out the following line
	  if(kstep.eq.1 .AND. kinc.eq.1) testing=1
	  if(norm_D.eq.0) testing=2
c 	  FEM asking for ddsdde only
      if(testing.eq.2) then
            do i=1,nyact        
                  y(i)=y_n(i)
            end do      
	  else if(testing.eq.1) then
          call rkf23_update_h(y,nyact,nasv,dtsub,tolintTtest,
     &                      maxninttest,DTmin,
     &                      deps_np1,parms,nparms,nfev,elprsw,dtime)
c	      give original state if the model fails without substepping
     	  if(error.eq.3) then
            do i=1,nyact        
               y(i)=y_n(i)
            end do
            error=0
          end if
c 	  Normal RKF23 integration
      else if(error.ne.10) then
                call rkf23_update_h(y,nyact,nasv,dtsub,tolintT,
     &                      maxnint,DTmin,
     &                      deps_np1,parms,nparms,nfev,elprsw,dtime)
      end if
      write(6,*) 'in umat_hypo, line 307'


c
c ... error conditions (if any)
c
      if (error.eq.3) then
c
c ... reduce abaqus load increment
c
                pnewdt = 0.25d0
c
                write(6,*) 'subroutine UMAT: reduce step size in ABAQUS'
                call wrista_h(1,y,nydim,deps_np1,dtime,
     &                coords,statev,nstatv,
     &                parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c               call xit_h
c               return
     
c ... do not do anything, we are the most likely close to the tensile region
            do i=1,nyact        
                  y(i)=y_n(i)
            end do
c
      elseif (error.eq.10) then
c
                call wrista_h(2,y,nydim,deps_np1,dtime,
     &                coords,statev,nstatv,
     &                parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
                call xit_h
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
      statev(13)=dtsub
        statev(10)=dfloat(nfev)
c
c ... compute consistent tangent via numerical perturbation
c
      write(6,*) 'in umat_hypo, line 352'

        call perturbate_h(y_n,y,nyact,nasv,dtsub,tolintT,maxnint,DTmin,
     &      deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan,dtime)
     
c
c ... convert solution (stress + cons. tangent) to abaqus format
c     update pore pressure and compute total stresses 
c
      call solout_h(stress,ntens,asv,nasv,ddsdde,
     +            y,nydim,pore,depsv_np1,parms,nparms,DDtan)
     
c
c ... updated vector of additional state variables to abaqus statev vector
c
      do i=1,nasv
        statev(i-1+nfasv) = asv(i) 
      end do
c
c ... transfer additional information to statev vector
c
      do i=1,6
                sig_np1(i)=y(i)
        end do
      pp=-(sig_np1(1)+sig_np1(2)+sig_np1(3))/3
c
      statev(8) = -pore 
      statev(9) = pp

      call calc_statev_h(sig_np1,statev,parms,nparms,nasv,
     & 	nasvdim,deps_np1)
      write(6,*) 'in umat_hypo, line 382'

c
c -----------------------
c End of time integration
c -----------------------
c
      return
      end
c------------------------------------------------------------------------------
c-----------------------------------------------------------------------------
      subroutine check_parms_h(props,nprops,parms,nparms)
c-----------------------------------------------------------------------------
c checks input material parameters 
c
c written 10/2004 (Tamagnini & Sellari)
c-----------------------------------------------------------------------------
      implicit none
c
      integer nprops,nparms,i,error
c
      double precision props(nprops),parms(nprops)
        double precision zero,one,four,pi,pi_deg
        double precision phi_deg,phi,hs,en,ed0,ec0,ei0,alpha,beta
        double precision m_R,m_T,r_uc,beta_r,chi,bulk_w,p_t
c
        parameter(zero=0.0d0,one=1.0d0,four=4.0d0,pi_deg=180.0d0)
c
        common /z_nct_errcode/error
c
        nparms=nprops
c
      do i=1,nprops
                parms(i)=props(i)
      enddo
c
c ... recover material parameters
c
        phi_deg=parms(1)
        hs    =parms(3)
        en    =parms(4)
        ed0   =parms(5)
        ec0   =parms(6)
        ei0   =parms(7)
        alpha =parms(8)
        beta  =parms(9)
        m_R=parms(10) 
        m_T=parms(11)
        r_uc=parms(12)
        beta_r=parms(13)
        chi=parms(14)
        bulk_w=parms(15)        
        p_t=parms(2)
        
c
      pi=four*datan(one)
        phi=phi_deg*pi/pi_deg
        parms(1)=phi
c
        if(phi.le.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'phi = ',phi
                error = 10
                return 
c
        end if
c
        if(m_R.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'm_R = ',m_R
                error = 10 
                return 
c
        end if
c
        if(m_T.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'm_T = ',m_T
                error = 10 
                return 
c
        end if
c
        if(r_uc.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'r_uc = ',r_uc
                error = 10 
                return 
c
        end if
c
        if(beta_r.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'beta_r = ',beta_r
                error = 10 
                return 
c
        end if
c
        if(chi.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'chi = ',chi
                error = 10 
                return 
c
        end if
c 
        if(bulk_w.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'bulk_w = ',bulk_w
                error = 10 
                return 
c
        end if
c 
        if(p_t.lt.zero) then
c       
                write(6,*) 'ERROR: subroutine CHECK_PARMS:'
                write(6,*) 'p_t = ',p_t
                error = 10 
                return 
c
        end if
c 
      return
      end
c-----------------------------------------------------------------------------
      subroutine define_h(nasv)
c-----------------------------------------------------------------------------
      implicit none 
      integer nasv
c
c number of additional state variables 
c must be less than  18 (otherwise change nasvdim in umat)
c
c    nasv(1) ... del_11  intergranular strain component
c    nasv(2) ... del_22  intergranular strain component
c    nasv(3) ... del_33  intergranular strain component
c    nasv(4) ... del_12  intergranular strain component
c    nasv(5) ... del_13  intergranular strain component
c    nasv(6) ... del_23  intergranular strain component
c    nasv(7) ... void    void ratio
c
c modified 6/2005 (Tamagnini, Sellari & Miriano)
c
      nasv = 7
      return
      end
c------------------------------------------------------------------------------
      double precision function dot_vect_h(flag,a,b,n)
c------------------------------------------------------------------------------
c dot product of a 2nd order tensor, stored in Voigt notation
c created 10/2004 (Tamagnini & Sellari)
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
        dot_vect_h=zero
c
        do i=1,n
                if(i.le.3) then
                      dot_vect_h = dot_vect_h+a(i)*b(i)
                else
                      dot_vect_h = dot_vect_h+coeff*a(i)*b(i)
                end if
        end do
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine get_F_sig_q_h(sig,q,nasv,parms,nparms,deps,F_sig,F_q)
c-----------------------------------------------------------------------------
c
c  finds vectors F_sigma and F_q in F(y)
c
c  written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
        implicit none
        double precision dot_vect_h
        
c 
      integer nparms,nasv,ii
c
        double precision sig(6),q(nasv),parms(nparms),deps(6)
        double precision MM(6,6),HH(nasv,6),F_sig(6),F_q(nasv)
        double precision LL(6,6),NN(6),norm_D,norm_D2
        integer istrain
c
c ... compute tangent operators
c
		if(parms(10) .le. 0.5) then
			istrain=0 
		else 
			istrain=1
		end if

        call get_tan_h(deps,sig,q,nasv,parms,nparms,MM,
     .        HH,LL,NN,istrain)
c
c ... compute F_sig=MM*deps
c

		if (istrain .eq. 1) then
        		call matmul_h(MM,deps,F_sig,6,6,1)
        else 
        		call matmul_h(LL,deps,F_sig,6,6,1)
		        norm_D2=dot_vect_h(2,deps,deps,6)
		        norm_D=sqrt(norm_D2)
                do ii=1,6
                     F_sig(ii)=F_sig(ii)+NN(ii)*norm_D
                end do
        end if

c
c ... compute F_q=HH*deps
c
        call matmul_h(HH,deps,F_q,nasv,6,1)
c       
        return
        end
c-----------------------------------------------------------------------------
      subroutine get_tan_h(deps,sig,q,nasv,parms,nparms,MM,HH,
     .		 LL,NN,istrain)
c-----------------------------------------------------------------------------
c  computes matrices M and H for Masin hypoplastic model for clays
c  version with intergranular strains
c
c  NOTE: stress and strain convention: tension and extension positive
c
c  written 6/2005 (Tamagnini & Sellari)
c-----------------------------------------------------------------------------
        implicit none
c 
      integer nparms,nasv,i,j,error
c
        double precision dot_vect_h
c
        double precision sig(6),q(nasv),parms(nparms),deps(6)
        double precision eta(6),eta_dev(6),del(6),void,sig_star(6)
        double precision eta_del(6),eta_delta(6),eta_eps(6)
        double precision norm_del,norm_del2,norm_deps,norm_deps2,eta_dn2
        double precision pp,qq,cos3t,I1,I2,I3,tanpsi
        double precision a,a2,FF,fd,fs
        double precision num,den,aF,Fa2,eta_n2,norm_m,norm_m2
        double precision II(6,6),IU(6,6)
        double precision MM(6,6),HH(nasv,6),LL(6,6),NN(6),AA(6,6),m(6)
        integer istrain
        double precision m_dir(6),m_dir1(6),Leta(6),H_del(6,6),H_e(6)
        double precision load,rho
        double precision zero,tiny,half,one,two,three,six,eight,nine
        double precision onethird,sqrt3,twosqrt2,sqrt2,oneeight,ln2m1
        double precision temp1,temp2,temp3,temp4
        double precision phi,hs,en,ed0,ec0,ei0,alpha,beta,r_uc
        double precision m_R,m_T,beta_r,chi,bulk_w,p_t,sinphi,sinphi2
        double precision ec,ed,ei,bauer,fb,fe,sq2,sq3,sq6,az
c
        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
        parameter(tiny=1.0d-17,half=0.5d0,eight=8.0d0,nine=9.0d0)
        parameter(sq2=1.4142135623730951455d0,
     &          sq3=1.7320508075688771931d0,
     &          sq6=2.4494897427831778813d0)

c
        common /z_nct_errcode/error     
c
c ... initialize constants and vectors
c
        onethird=one/three
        sqrt3=dsqrt(three)
        twosqrt2=two*dsqrt(two)
        sqrt2=dsqrt(two)
        oneeight=one/eight
        onethird=one/three
        ln2m1=one/dlog(two)
c
        do i=1,6
                do j=1,6
                        MM(i,j)=zero
                        LL(i,j)=zero
                        II(i,j)=zero
                        IU(i,j)=zero
                        H_del(i,j)=zero
                end do
                eta_del(i)=zero
                eta_delta(i)=zero
                eta_eps(i)=zero
        end do
c
        do i=1,nasv
                do j=1,6
                        HH(i,j)=zero
                end do
        end do
c
c ... fourth order identity tensors in Voigt notation
c
        II(1,1)=one
        II(2,2)=one
        II(3,3)=one
        II(4,4)=half
        II(5,5)=half
        II(6,6)=half
c
        IU(1,1)=one
        IU(2,2)=one
        IU(3,3)=one
        IU(4,4)=one
        IU(5,5)=one
        IU(6,6)=one
c
c ... recover material parameters
c
        phi	  =parms(1)
        hs    =parms(3)
        en    =parms(4)
        ed0   =parms(5)
        ec0   =parms(6)
        ei0   =parms(7)
        alpha =parms(8)
        beta  =parms(9)
        m_R=parms(10) 
        m_T=parms(11)
        r_uc=parms(12)
        beta_r=parms(13)
        chi=parms(14)
        bulk_w=parms(15)
        p_t=parms(2)
c
        sinphi=dsin(phi)
        sinphi2=sinphi*sinphi

c
c ... recover internal state variables
c
        del(1)=q(1)
        del(2)=q(2)
        del(3)=q(3)
        del(4)=q(4)
        del(5)=q(5)
        del(6)=q(6)
        void=q(7)
c
c ... axis translation due to cohesion (p_t>0)
c
        sig_star(1)=sig(1)-p_t
        sig_star(2)=sig(2)-p_t
        sig_star(3)=sig(3)-p_t
        sig_star(4)=sig(4)
        sig_star(5)=sig(5)
        sig_star(6)=sig(6)
c
c ... strain increment and intergranular strain directions
c
        norm_deps2=dot_vect_h(2,deps,deps,6)
        norm_del2=dot_vect_h(2,del,del,6)
        norm_deps=dsqrt(norm_deps2)
        norm_del=dsqrt(norm_del2)
c
        if(norm_del.ge.tiny) then
c
                do i=1,6
                        eta_del(i)=del(i)/norm_del
                end do
c
        end if
c
        eta_delta(1)=eta_del(1)
        eta_delta(2)=eta_del(2)
        eta_delta(3)=eta_del(3)
        eta_delta(4)=half*eta_del(4)
        eta_delta(5)=half*eta_del(5)
        eta_delta(6)=half*eta_del(6)
c
        if(norm_deps.ge.tiny) then
c
                do i=1,6
                        eta_eps(i)=deps(i)/norm_deps
                end do
c
        end if
c
c ... auxiliary stress tensors
c
        call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
c
c        if (pp.gt.tiny) then
c
c ... if mean stress is negative, return with MM = 0, HH = 0 and error = 10 (severe)
c
c                write(6,*) 'ERROR: subroutine GET_TAN:'
c                write(6,*) 'Mean stress is positive (tension): p = ',pp
c                error = 10
c                return 
c
c        end if
c
        
        eta(1)=sig_star(1)/I1
        eta(2)=sig_star(2)/I1
        eta(3)=sig_star(3)/I1
        eta(4)=sig_star(4)/I1
        eta(5)=sig_star(5)/I1
        eta(6)=sig_star(6)/I1   
c
        eta_dev(1)=eta(1)-onethird
        eta_dev(2)=eta(2)-onethird
        eta_dev(3)=eta(3)-onethird
        eta_dev(4)=eta(4)
        eta_dev(5)=eta(5)
        eta_dev(6)=eta(6)
c
c ... functions a and F
c
        eta_dn2=dot_vect_h(1,eta_dev,eta_dev,6)
        tanpsi=sqrt3*dsqrt(eta_dn2)
        temp1=oneeight*tanpsi*tanpsi+
     &    (two-tanpsi*tanpsi)/(two+sqrt2*tanpsi*cos3t)
        temp2=tanpsi/twosqrt2
c
        a=sqrt3*(three-sin(phi))/(twosqrt2*sin(phi))
        a2=a*a
        FF=dsqrt(temp1)-temp2
c
c ... barotropy and pyknotropy functions
c
	    bauer=dexp(-(-I1/hs)**en)
      	ed = ed0*bauer
      	ec = ec0*bauer
      	ei = ei0*bauer

      	temp1=three+a*a-a*sq3*((ei0-ed0)/(ec0-ed0))**alpha
      	if(temp1.lt.zero) stop 'factor fb not defined'
	    fb=hs/en/temp1*(one+ei)/ei*(ei0/ec0)**beta*(-I1/hs)**(one-en)
    	fe=(ec/void)**beta
            	
        fs=fb*fe
c
		if(void.ge.ed) then
        	fd=((void-ed)/(ec-ed))**alpha
        else
        	fd=0
        end if
c
c
c ... tensor L
c
        eta_n2=dot_vect_h(1,eta,eta,6)
        do i = 1,6
                do j=1,6
                        LL(i,j)=(II(i,j)*FF*FF+
     &	                  a2*eta(i)*eta(j))/eta_n2
                end do
        end do

c
c ... tensor NN
c

       do i=1,6
         NN(i) = FF*a*(eta(i)+eta_dev(i))/eta_n2
       enddo
        
c
c ... BEGIN INTERGR. STRAIN
c

        if(istrain .eq. 1) then
c
c ... loading function
c
        load=dot_vect_h(2,eta_del,eta_eps,6)
c
c ... intergranular strain--related tensors
c
        rho=norm_del/r_uc
c
        if (rho.gt.one) then
                rho=one
        end if
c
        call matmul_h(LL,eta_del,Leta,6,6,1)
c
c ... tangent stiffness M(sig,q,eta_eps)
c
        temp1=((rho**chi)*m_T+(one-rho**chi)*m_R)*fs
c
        if (load.gt.zero) then
c    
                temp2=(rho**chi)*(one-m_T)*fs
                temp3=(rho**chi)*fs*fd
c
                do i=1,6
                  do j=1,6
                    AA(i,j)=temp2*Leta(i)*eta_delta(j)
     &                      +temp3*NN(i)*eta_delta(j)
                    MM(i,j)=temp1*LL(i,j)+AA(i,j)
                  end do
                end do
c
        else
c
                temp4=(rho**chi)*(m_R-m_T)*fs
c
                do i=1,6
                  do j=1,6
                        AA(i,j)=temp4*Leta(i)*eta_delta(j)
                        MM(i,j)=temp1*LL(i,j)+AA(i,j)
                  end do
                end do
c
        end if
c
c ... intergranular strain evolution function
c     NOTE: H_del transforms a strain-like vector into a strain-like vector
c           eta_del(i) instead of eta_delta(i)
c           I = 6x6 unit matrix
c
        if (load.gt.zero) then
c
                do i=1,6
                  do j=1,6
                H_del(i,j)=IU(i,j)-(rho**beta_r)*eta_del(i)*eta_delta(j)
                  end do
                end do
c
        else
c
                do i=1,6
              H_del(i,i)=one
                end do
c
        end if
c
c ... void ratio evolution function (tension positive)
c
        do i=1,6 
                if (i.le.3) then
                  H_e(i)=one+void
                else
              H_e(i)=zero
                end if
        end do
c
c ... assemble hardening matrix
c
        do i=1,nasv
                if (i.le.6) then
                        do j=1,6
                                HH(i,j)=H_del(i,j)
                        end do
                else
                        do j=1,6
                                HH(i,j)=H_e(j)
                        end do
                end if
        end do
c       
c ... end istrain
        else if (istrain .eq. 0) then
        
        do i=1,6 
                if (i.le.3) then
                  H_e(i)=one+void
                else
              H_e(i)=zero
                end if
        end do        
        do i=1,nasv
                if (i.le.6) then
                        do j=1,6
                                HH(i,j)=0
                        end do
                else
                        do j=1,6
                                HH(i,j)=H_e(j)
                        end do
                end if
        end do        
c ... end istrain/noistrain switch        
        end if

        do i=1,6
           do j=1,6
                LL(i,j)=LL(i,j)*fs
           end do
           NN(i)=NN(i)*fs*fd
        end do        

        return
        end
c-----------------------------------------------------------------------------
      subroutine iniy_h(y,nydim,nasv,ntens,sig,qq)
c-----------------------------------------------------------------------------
c initializes the vector of state variables
c-----------------------------------------------------------------------------
      implicit none
c
      integer i,nydim,nasv,ntens
c
      double precision y(nydim),qq(nasv),sig(ntens)
c
      do i=1,nydim
        y(i) = 0
      enddo
c
      do i=1,ntens
        y(i) = sig(i)
      enddo
c
c additional state variables
c
      do i=1,nasv
        y(6+i) = qq(i)
      enddo
c
      return
      end
c------------------------------------------------------------------------------
      subroutine inv_eps_h(eps,eps_v,eps_s,sin3t)
c------------------------------------------------------------------------------
c calculate invariants of strain tensor
c------------------------------------------------------------------------------
c
      implicit none
c
      integer i
c
      double precision eps(6),edev(6),edev2(6),ev3
        double precision tredev3,eps_v,eps_s,sin3t
        double precision norm2,numer,denom
c
      double precision zero,one,two,three,six
      double precision onethird,twothirds,sqrt6
c
      data zero,one,two,three,six/0.0d0,1.0d0,2.0d0,3.0d0,6.0d0/
c
c ... some constants
c
        onethird=one/three
        twothirds=two/three
        sqrt6=dsqrt(six)
c
c ... volumetric strain
c
      eps_v=eps(1)+eps(2)+eps(3)
c
      ev3=onethird*eps_v
c
c ... deviator strain
c
        edev(1)=eps(1)-ev3
        edev(2)=eps(2)-ev3
        edev(3)=eps(3)-ev3
        edev(4)=eps(4)/two
        edev(5)=eps(5)/two
        edev(6)=eps(6)/two
c
c ... second invariant
c
        norm2=edev(1)*edev(1)+edev(2)*edev(2)+edev(3)*edev(3)+
     &      two*(edev(4)*edev(4)+edev(5)*edev(5)+edev(6)*edev(6))
c
        eps_s=dsqrt(twothirds*norm2)
c
c ... components of (edev_ij)(edev_jk)
c
        edev2(1)=edev(1)*edev(1)+edev(4)*edev(4)+edev(5)*edev(5)
        edev2(2)=edev(4)*edev(4)+edev(2)*edev(2)+edev(6)*edev(6)
        edev2(3)=edev(6)*edev(6)+edev(5)*edev(5)+edev(3)*edev(3)
        edev2(4)=two*(edev(1)*edev(4)+edev(4)*edev(2)+edev(6)*edev(5))
        edev2(5)=two*(edev(5)*edev(1)+edev(6)*edev(4)+edev(3)*edev(5))
        edev2(6)=two*(edev(4)*edev(5)+edev(2)*edev(6)+edev(6)*edev(3))
c            
c ... Lode angle
c
        if(eps_s.eq.zero) then 
c
                sin3t=-one
c               
        else
c
                tredev3=zero
                do i=1,6
                        tredev3=tredev3+edev(i)*edev2(i)
                end do
c
                numer=sqrt6*tredev3
                denom=(dsqrt(norm2))**3
                sin3t=numer/denom
                if(dabs(sin3t).gt.one) then
                        sin3t=sin3t/dabs(sin3t)
                end if
c
        end if 
c
      return
      end
c------------------------------------------------------------------------------
      subroutine inv_sig_h(sig,pp,qq,cos3t,I1,I2,I3)
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
        double precision dot_vect_h
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
        eta(1)=sig(1)/I1
        eta(2)=sig(2)/I1
        eta(3)=sig(3)/I1
        eta(4)=sig(4)/I1
        eta(5)=sig(5)/I1
        eta(6)=sig(6)/I1
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
        norm2=dot_vect_h(1,sdev,sdev,6)
        norm2sig=dot_vect_h(1,sig,sig,6)
        norm2eta=dot_vect_h(1,eta_d,eta_d,6)
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
                tretadev3=dot_vect_h(1,eta_d,eta_d2,6)
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
c------------------------------------------------------------------------------
      subroutine matmul_h(a,b,c,l,m,n)
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
c-----------------------------------------------------------------------------
      subroutine move_asv_h(asv,nasv,qq_n)
c-----------------------------------------------------------------------------
c move internal variables in vector qq_n and changes intergranular strain 
c from continuum to soil mechanics convention
c
c NOTE: del has always 6 components
c
c written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer nasv,i
      double precision asv(nasv),qq_n(nasv),zero 
c
        parameter(zero=0.0d0)
c
      do i=1,nasv
                qq_n(i)=zero
      enddo
c
c ... intergranular strain tensor stored in qq_n(1:6)
c
      do i=1,6
                qq_n(i) = -asv(i)
      enddo
c
c ... void ratio stored in qq_n(7)
c
        qq_n(7) = asv(7) 
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine move_eps_h(dstran,ntens,deps,depsv)
c-----------------------------------------------------------------------------
c Move strain increment dstran into deps and computes 
c volumetric strain increment
c
c NOTE: all strains negative in compression; deps has always 6 components
c
c written 7/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision deps(6),dstran(ntens),depsv
c
      do i=1,ntens
                deps(i) = dstran(i)
      enddo
c
        depsv=deps(1)+deps(2)+deps(3)
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine move_sig_h(stress,ntens,pore,sig)
c-----------------------------------------------------------------------------
c computes effective stress from total stress (stress) and pore pressure (pore)
c
c NOTE: stress = total stress tensor (tension positive)
c         pore   = exc. pore pressure (undrained conds., compression positive)
c         sig    = effective stress (tension positive)
c
c       sig has always 6 components
c
c written 7/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
      implicit none
      integer ntens,i
      double precision sig(6),stress(ntens),pore,zero 
c
        parameter(zero=0.0d0)
c
      do i=1,6
                sig(i)=zero
      enddo
c
      do i=1,ntens
                if(i.le.3) then
                        sig(i) = stress(i)+pore
                else
                        sig(i) = stress(i)
                end if
      enddo
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine norm_res_h(y_til,y_hat,ny,nasv,norm_R)
c-----------------------------------------------------------------------------
c  evaluate norm of residual vector Res=||y_hat-y_til||
c
c  written 6/2005 (Tamagnini, Sellari & Miriano)
c-----------------------------------------------------------------------------
        implicit none
c 
      integer ny,nasv,ng,k,i,testnan
c
        double precision y_til(ny),y_hat(ny),void_til,void_hat,del_void
        double precision err(ny),norm_R2,norm_R
        double precision norm_sig2,norm_q2,norm_sig,norm_q
        double precision sig_hat(6),sig_til(6),del_sig(6)
        double precision q_hat(nasv),q_til(nasv),del_q(nasv)
        double precision dot_vect_h,zero
c
        parameter(zero=0.0d0)
c
        ng=6*nasv
        k=42+nasv
c
        do i=1,ny
                err(i)=zero
        end do
c
c ... recover stress tensor and internal variables
c
        do i=1,6
                sig_hat(i)=y_hat(i)
                sig_til(i)=y_til(i)
                del_sig(i)=dabs(sig_hat(i)-sig_til(i))
        end do
c
        do i=1,nasv-1
                q_hat(i)=y_hat(6+i)
                q_til(i)=y_til(6+i)
                del_q(i)=dabs(q_hat(i)-q_til(i))
        end do
c
        void_hat=y_hat(6+nasv)
        void_til=y_til(6+nasv)
        del_void=dabs(void_hat-void_til)
c
c ... relative error norms
c
        norm_sig2=dot_vect_h(1,sig_hat,sig_hat,6)
        norm_q2=dot_vect_h(2,q_hat,q_hat,6)
        norm_sig=dsqrt(norm_sig2)
        norm_q=dsqrt(norm_q2)
c
        if(norm_sig.gt.zero) then
                do i=1,6
                        err(i)=del_sig(i)/norm_sig
                end do
        end if
c
        if(norm_q.gt.zero) then
                do i=1,nasv-1
                err(6+i)=del_q(i)/norm_q
                end do
        end if
c
        err(6+nasv)=del_void/void_hat
c
c ... global relative error norm
c
        norm_R2=dot_vect_h(3,err,err,ny)
        norm_R=dsqrt(norm_R2)
c
 	  	testnan=0
  	    call umatisnan_h(norm_sig,testnan)
  	    call umatisnan_h(norm_q,testnan)
  	    call umatisnan_h(void_hat,testnan)
		if(testnan.eq.1) then
			norm_R=1.d20
		end if

        return
        end

c-----------------------------------------------------------------------------
      subroutine perturbate_h(y_n,y_np1,n,nasv,dtsub,err_tol,maxnint,
     &    DTmin,deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DD, dtime)
c-----------------------------------------------------------------------------
c
c  compute numerically consistent tangent stiffness
c
c  written 12/2005 (Tamagnini)
c-----------------------------------------------------------------------------
        implicit none
c 
        logical elprsw
c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev
        integer maxnint,error
c
      double precision y_n(n),y_np1(n),y_star(n),parms(nparms)
        double precision dtsub,err_tol,DTmin, dtime
        double precision theta,sig(6),q(nasv)
        double precision deps_np1(6),deps_star(6)
        double precision dsig(6),DD(6,6),HHtmp(nasv,6)
        double precision LL(6,6),NN(6)
        integer istrain
        double precision zero
c
      parameter(zero=0.0d0)
c
        common /z_nct_errcode/error
c
c ... initialize DD and y_star
c 
		if(parms(10) .le. 0.5) then
			istrain=0 
		else 
			istrain=1
		end if

        do kk=1,6
                do jj=1,6
                        DD(kk,jj)=zero
                end do
        end do
        do i=1,6
                sig(i)=y_n(i)
        end do
        do i=1,nasv
                q(i)=y_n(6+i)
        end do
        
        call push_h(y_n,y_star,n)

        if(error.ne.10) then
          call get_tan_h(deps_np1,sig,q,nasv,parms,nparms,
     .          	DD,HHtmp,LL,NN,istrain)                
        end if
        if(istrain .eq. 0) then
          do kk=1,6
                do jj=1,6
                        DD(kk,jj)=LL(kk,jj)
                end do
          end do
        end if

        return
        end        
        
c-----------------------------------------------------------------------------
      subroutine push_h(a,b,n)
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
c-----------------------------------------------------------------------------
      subroutine rhs_h(y,ny,nasv,parms,nparms,deps,kRK,nfev)
c-----------------------------------------------------------------------------
c calculate coefficient kRK from current state y and strain increment deps
c Masin hypoplastic model for clays with intergranular strains
c
c written 12/2005 (Tamagnini & Sellari)
c-----------------------------------------------------------------------------
      implicit none
c
        integer error,ny,nparms,nasv,i,nfev
c
      double precision zero,one,two,four 
        double precision y(ny),kRK(ny),parms(nparms),deps(6)
        double precision sig(6),q(nasv)
        double precision F_sig(6),F_q(nasv)
c
        common /z_nct_errcode/error     
c
        parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0)
c
c ... update counter for the number of function f(y) evaluations
c
        nfev=nfev+1
c
c ... initialize kRK
c
        do i=1,ny
                kRK(i)=zero
        end do
c
c ... recover current state variables (sig,q)                   
c
        do i=1,6
                sig(i)=y(i)
        end do
c
      do i=1,nasv
                q(i)=y(6+i)
        end do
c       
c ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
c
        call get_F_sig_q_h(sig,q,nasv,parms,nparms,deps,F_sig,F_q)
        if(error.eq.10) return
c
        do i=1,6
c
                kRK(i)=F_sig(i)
c
        end do                   
c       
        do i=1,nasv
c
                kRK(6+i)=F_q(i)
c
        end do                   
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine rkf23_update_h(y,n,nasv,dtsub,err_tol,maxnint,DTmin,
     &                        deps_np1,parms,nparms,nfev,elprsw,dtime)
c-----------------------------------------------------------------------------
c
c  numerical solution of y'=f(y)
c  explicit, adapive RKF23 scheme with local time step extrapolation
c
c  Tamagnini, Sellari & Miriano 6/2005
c
c-----------------------------------------------------------------------------
        implicit none
c
        logical elprsw
c
      integer n,nasv,nparms,i,ksubst,kreject,nfev
        integer maxnint,error,error_RKF
c
      double precision y(n),parms(nparms),dtsub,err_tol,DTmin
        double precision zero,half,one,two,three,four,six
        double precision ptnine,onesixth,onethird,twothirds,temp
c
        double precision deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
        double precision y_hat(n)
        double precision T_k,DT_k,dtime
        double precision kRK_1(n),kRK_2(n),kRK_3(n)
        double precision norm_R,S_hull
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0,ptnine=0.9d0)
c
        common /z_nct_errcode/error     
c
c ... initialize y_k vector and other variables
c
        do i=1,n
                y_k(i)=zero
        end do
c
        onesixth=one/six
        onethird=one/three
        twothirds=two/three
c
c ... start of update process
c
        T_k=zero      
        DT_k=dtsub/dtime
        ksubst=0
        kreject=0
        nfev=0
c
        do i=1,n
                y_k(i)=y(i)
        end do
c
c ... start substepping 
c
        do while(T_k.lt.one) 
c
                ksubst=ksubst+1
c
c ... write substepping info
c
c               write(*,1234) ksubst,T_k,DT_k
c1234           format('Substep no.',i4,' -- T_k = ',d12.4,' -- DT_k = ',d12.4)
c
c ... check for maximum number of substeps
c
                if(ksubst.gt.maxnint) then
                       	write(*,*) 'number of substeps ',ksubst,
     &                             ' is too big, step rejected'
                        error=3
                        return
                end if          
c
c ... build RK functions
c

                call rhs_h(y_k,n,nasv,parms,nparms,deps_np1,kRK_1,nfev)
                if(error.eq.10) return
c
c ... find y_2
c
                temp=half*DT_k
c
                 do i=1,n
                        y_2(i)=y_k(i)+temp*kRK_1(i)
                end do

c               
                call rhs_h(y_2,n,nasv,parms,nparms,deps_np1,kRK_2,nfev)
                if(error.eq.10) return
c                                       
c ... find y_3
c
                do i=1,n
                        y_3(i)=y_k(i)-DT_k*kRK_1(i)+two*DT_k*kRK_2(i)
                end do
c

                call rhs_h(y_3,n,nasv,parms,nparms,deps_np1,kRK_3,nfev)
                
                if(error.eq.10) return
c                               
c ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
c
                do i=1,n        

                        y_til(i)=y_k(i)+DT_k*kRK_2(i)
                        y_hat(i)=y_k(i)+DT_k*
     &          (onesixth*kRK_1(i)+twothirds*kRK_2(i)+onesixth*kRK_3(i))
                end do
c
c ... local error estimate
c

                call norm_res_h(y_til,y_hat,n,nasv,norm_R)
c				check if output y_hat can be used as an input into the next step
				error_RKF=0
                call check_RKF_h(error_RKF,y_hat,n,nasv,parms,nparms)

                if (error_RKF.ne.0) then
                	norm_R=1.d20
                	error_RKF=0
                end if
c
c ... time step size estimator according to Hull
c       	
				if(norm_R .ne. 0) then
                	S_hull=ptnine*DT_k*(err_tol/norm_R)**onethird
                else
                	S_hull=1
                end if
c

        if (norm_R.lt.err_tol) then                             
c
c ... substep is accepted, update y_k and T_k and estimate new substep size DT_k
c
                        do i=1,n        
                                y_k(i)=y_hat(i)
                        end do
c
                        T_k=T_k+DT_k
                        DT_k=min(four*DT_k,S_hull)
						dtsub=DT_k*dtime
                        DT_k=min((one-T_k),DT_k)        
c
                else
c
c ... substep is not accepted, recompute with new (smaller) substep size DT
c
                        DT_k=max(DT_k/four,S_hull)
c
c ... check for minimum step size
c
                        if(DT_k.lt.DTmin) then
                              write(6,*) 'substep size ',DT_k,
     &                             ' is too small, step rejected'
                              error=3
                              return
                        end if          
c                                       
                end if                                                  
c
c ... bottom of while loop
c
        end do
        
c
c ... recover final state
c
        do i=1,n
                y(i)=y_k(i)
        end do
c

        return
      end
c

c-----------------------------------------------------------------------------
      subroutine check_RKF_h(error_RKF,y,ny,nasv,parms,nparms)
c-----------------------------------------------------------------------------
c Checks is RKF23 solout vector y is OK for hypoplasticity
c-----------------------------------------------------------------------------
      implicit none
c
        integer error_RKF,ny,nasv,i,nparms,testnan,iopt
c
        double precision y(ny),parms(nparms)
        double precision sig(6),pmean,sig_star(6)
        double precision xN1(3),xN2(3),xN3(3),S(3),P,Q,tmin
        double precision p_t
c
        p_t    =parms(2)
        do i=1,6
                sig(i)=y(i)
        end do

        sig_star(1)=sig(1)-p_t
        sig_star(2)=sig(2)-p_t
        sig_star(3)=sig(3)-p_t
        sig_star(4)=sig(4)
        sig_star(5)=sig(5)
        sig_star(6)=sig(6)
                
    	pmean=-(sig_star(1)+sig_star(2)+sig_star(3))/3
    	
c		check for positive mean stress
        if(pmean .le. 0) then
        	error_RKF=1
        end if
c
c		calculate minimum principal stress
c
		iopt=0
        Call PrnSig_h(iopt, sig_star, xN1, xN2, xN3,
     &		S(1),S(2),S(3), P, Q)
        tmin     = 1.0d+20
        do i=1,3
                if(tmin .ge. -S(i)) then
                	 tmin=-S(i)
                endif	 
        enddo 
c
c		check for tension
c
        if(tmin .le. 0) then
                error_RKF=1
        end if
        
c		check for NAN
 	  	testnan=0
        do i=1,ny
       	  call umatisnan_h(y(i),testnan)
        end do
        if(testnan.eq.1) error_RKF=1
        
c
c
      return
      end
c

c-----------------------------------------------------------------------------
      subroutine solout_h(stress,ntens,asv,nasv,ddsdde,y,nydim,
     +                  pore,depsv_np1,parms,nparms,DD)
c-----------------------------------------------------------------------------
c copy the vector of state variables to umat output
c modified 7/2005 (Tamagnini, Sellari)
c
c NOTE: solid mechanics convention for stress and strain components
c       pore is always positive in compression
c-----------------------------------------------------------------------------
      implicit none
c
      integer nydim,nasv,nparms,ntens,i,j
c
      double precision y(nydim),asv(nasv),stress(ntens)
        double precision ddsdde(ntens,ntens),DD(6,6)
        double precision parms(nparms),bulk_w,pore,depsv_np1 
c
        bulk_w=parms(15)
c
c ... update excess pore pressure (if undrained conditions), compression positive
c
        pore=pore-bulk_w*depsv_np1
c
c updated total stresses (effective stresses stored in y(1:6))
c
      do i=1,ntens
                if (i.le.3) then
                        stress(i) = y(i)-pore
                else
                        stress(i) = y(i)
                end if
        enddo
c
c additional state variables (first 6 components are intergranular strains)
c
      do i=1,nasv
                asv(i) = y(6+i)
      enddo
c
c consistent tangent stiffness
c
      do j=1,ntens
        do i=1,ntens
          ddsdde(i,j) = DD(i,j)      
        enddo
      enddo
c
      do j=1,3
        do i=1,3
          ddsdde(i,j) = ddsdde(i,j)+bulk_w        
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------------
      subroutine wrista_h(mode,y,nydim,deps_np1,dtime,coords,statev,
     &           nstatv,parms,nparms,noel,npt,ndi,nshr,kstep,kinc)
c-----------------------------------------------------------------------------
c ... subroutine for managing output messages
c
c     mode
c
c     all = writes:             kstep, kinc, noel, npt
c       2   = writes also:      error message,coords(3),parms(nparms),ndi,nshr,stress(nstress)
c                                               deps(nstress),dtime,statev(nstatv)
c     3   = writes also:        stress(nstress),deps(nstress),dtime,statev(nstatv)
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
        write(6,*) 'ERROR: abaqus job failed during call of UMAT'
        write(6,*) '==================================================='
        write(6,*) 'state dump:'
        write(6,*) 
      endif
c
c ... writes for all mode values
c
      write(6,111) 'Step: ',kstep, 'increment: ',kinc,
     & 'element: ', noel, 'Integration point: ',npt
      write(6,*) 
c
c ... writes for mode = 2
c
      if (mode.eq.2) then
        write(6,*) 'Co-ordinates of material point:'
        write(6,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',
     &    coords(3)
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
        write(6,*) 'Internal variables:'
        write(6,*) 
        write(6,109) 'del(1) = ',statev(1)
        write(6,109) 'del(2) = ',statev(2)
        write(6,109) 'del(3) = ',statev(3)
        write(6,109) 'del(4) = ',statev(4)
        write(6,109) 'del(5) = ',statev(5)
        write(6,109) 'del(6) = ',statev(6)
        write(6,109) 'void   = ',statev(7)
        write(6,*) 
        write(6,*) '==================================================='
c
      endif
c
101   format(1X,a15,e10.4)
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,f12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,f12.4)
109   format(1X,a6,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
c       
      return
      end

      
c-----------------------------------------------------------------------------
      subroutine calc_statev_h(stress,statev,parms,nparms,nasv,
     & nasvdim,deps)
c-----------------------------------------------------------------------------
c
c  computes additional state variables for postprocessing
c
c-----------------------------------------------------------------------------
        implicit none
c 
        logical elprsw
c
      integer ntens,jj,kk,i
      integer n,nasv,nparms,nfev,nasvdim
        integer maxnint,error
c
      double precision parms(nparms),dot_vect_h
        double precision stress(6),statev(nasvdim)
        double precision deps(6),tmax,tmin
        double precision MM(6,6),HHtmp(nasv,6)
        double precision LL(6,6),NN(6)
        integer istrain
        double precision zero,two,four,iopt,three
        double precision I1,I2,I3,cos3t,pp,qq
        double precision sin2phi,sinphi,sig_star(6),p_t
        double precision norm_del,norm_del2,del(6)
c
      parameter(zero=0.0d0,two=2.0d0,four=4.0d0,three=3.0d0)
c

c ... calc phimob (statev 11) from Matsuoka-Nakai YS

      p_t    =parms(2)
      do i=1,3
              sig_star(i)=stress(i)-p_t
      end do
      do i=4,6
              sig_star(i)=stress(i)
      end do

      call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
	  if(I3 .ne. 0) then
        sin2phi=(9.d0+I1*I2/I3)/(1.d0+I1*I2/I3)
      else 
      	sin2phi=0
      end if
	  if(sin2phi .lt. 0) then
        sin2phi=0
      end if 
	  if(sin2phi .gt. 1) then
        sin2phi=1
      end if 
      sinphi=sqrt(sin2phi)
      
      statev(11)= asin(sinphi)*
     .				   180.0d0/3.141592d0
 		

c ... calc norm. length of intergr. strain rho (statev 12)
		if(parms(10) .le. 0.5) then
			istrain=0 
		else 
			istrain=1
		end if

		if(istrain .eq. 1) then
        
        	do i=1,6
				del(i)=statev(i)
        	enddo       
        
	        norm_del2=dot_vect_h(2,del,del,6)
    	    norm_del=dsqrt(norm_del2)
        	statev(12)=norm_del/parms(12)
        	
        else
        	statev(12)=0
        end if

        return
        end        
            
c-----------------------------------------------------------------------------
      subroutine umatisnan_h(chcknum,testnan)
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
      
c-----------------------------------------------------------------------------
        subroutine xit_h
c-----------------------------------------------------------------------------
        stop
c
        return
        end

C***********************************************************************
      Subroutine PrnSig_h(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3_h(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a_h(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
C***********************************************************************
      Subroutine Eig_3_h(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),
     *          xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            v1p=c*v(1,ip)-s*v(1,iq)
            v2p=c*v(2,ip)-s*v(2,iq)
            v3p=c*v(3,ip)-s*v(3,iq)
            v(1,iq)=s*v(1,ip)+c*v(1,iq)
            v(2,iq)=s*v(2,ip)+c*v(2,iq)
            v(3,iq)=s*v(3,ip)+c*v(3,iq)
            v(1,ip)=v1p
            v(2,ip)=v2p
            v(3,ip)=v3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
      End Do
      Return
      End ! Eig_3

      Subroutine Eig_3a_h(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      Return
      End ! Eig_3a
      
