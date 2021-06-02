//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Montanino
//

#include "custom_elements/compressible_biphase_dg_navier_stokes_explicit.h"

namespace Kratos {

unsigned int	conta_dg(
    unsigned int i, 
    unsigned int j, 
    unsigned int k, 
    unsigned int l, 
    unsigned int maxi, 
    unsigned int maxj, 
    unsigned int maxk, 
    unsigned int maxl)
{
	
	unsigned int s;
	
	s = l + maxl*k + maxl*maxk*j + maxl*maxk*maxj*i;
	
	return s;
}
/*
void ShockCapturing(const double mu,
               const double lambda,
               const double c_v,
			   const double ros,
               const double h,
			   array_1d<double,2>& dt_diff,
			   array_1d<double,2>& ds_diff,
               array_1d<double,4>& tau,
               array_1d<double,2>& q,
               double ro_el,
               array_1d<double,10>& gradU, // 5*2
               array_1d<double,5>& Residual,
			   array_1d<double,5>& U,
			   array_1d<double,2>& a_el,
			   double norm_u,
			   double SpeedSound)
{
    const int SpaceDimension  = 2;
	
    double alpha  	= 2.0;
	double beta  	= 0.5;
    const double tol  = 1e-32;
    const double tol2 = 1e-32;                               

    unsigned int i;

    double a_ds = 0.0;
	double a_dt = 0.0;
	double v_sc = 0.0;
    double k_sc = 0.0;
    
	array_1d<double,SpaceDimension>  res_m;
	array_1d<double,SpaceDimension>  t_el;
    double ds_ref 	= U(1);
	double ds_ref2 	= U(1)*U(1);
	double dt_ref 	= U(0);
	double dt_ref2 	= U(0)*U(0);
	double res_dt 	= Residual(0);
	double res_ds 	= Residual(1);
	double res_e 	= Residual(4);
    double norm_res_ds;
	double norm_res_dt;
	double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;
	    
	t_el(0) = -a_el(1);
	t_el(1) =  a_el(0);

    res_m(0) = Residual(2); 
    res_m(1) = Residual(3);

	norm_res_ds = sqrt(res_ds*res_ds);
	norm_res_dt = sqrt(res_dt*res_dt);
    norm_res_m 	= sqrt(res_m(0)*res_m(0) + res_m(1)*res_m(1));
    norm_res_e 	= sqrt(res_e*res_e);

	double normnormdt = 0.0;              // Frobenius norm of dust density
	double normnormds = 0.0;              // Frobenius norm of dust density
    
	double ZsuY = sqrt(res_ds*res_ds/ds_ref2);
	double ZsuYt = sqrt(res_dt*res_dt/dt_ref2);

	for (i = 0; i < 2; i++)      normnormdt += gradU(i)*gradU(i)/dt_ref2;
	for (i = 2; i < 4; i++)      normnormds += gradU(i)*gradU(i)/ds_ref2;
	
	for (i = 4; i < 8; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    
    
	double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 8; i < 10; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
	if (norm_grade > tol)         k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);
	
	if (fabs(ds_ref) < tol2)	  
		a_ds = 0.0;		
	else a_ds = ZsuY*pow(normnormds,0.5*beta - 1.0)*pow(0.5*h,beta);
	
	a_dt = ZsuYt*pow(normnormdt,0.5*beta - 1.0)*pow(0.5*h,beta);
	
	for (i = 0; i < 2; i++)       dt_diff[i] = 0.0; //a_dt*gradU[i];
	for (i = 0; i < 2; i++)       ds_diff[i] = a_ds*gradU[i+2];

	for (i = 0; i < 4; i++)       tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 2; i++)       q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);

}
*/

void ShockCapturing2dg(const double mu,
               const double lambda,
               const double c_v,
			   const double ros,
               const double h,
			   array_1d<double,2>& dt_diff,
			   array_1d<double,2>& ds_diff,
               array_1d<double,4>& tau,
               array_1d<double,2>& q,
               double ro_el,
               array_1d<double,10>& gradU, // 5*2
               array_1d<double,5>& Residual,
			   array_1d<double,5>& U,
			   array_1d<double,2>& a_el,
			   double norm_u,
			   double SpeedSound)
{
    const int SpaceDimension  = 2;
	
    double alpha  	  = 0.5;
	double alpha_dc   = 1.0;
	double alpha_de	  = 0.5;
	const double tol  = 1e-32;
    const double tol2 = 1e-32;                               

    unsigned int i;

    double a_ds = 0.0;
	double a_dt = 0.0;
	double v_sc = 0.0;
    double k_sc = 0.0;
    
	array_1d<double,SpaceDimension>  res_m;
	array_1d<double,SpaceDimension>  t_el;
    double dt_ref 	= U(0);
	double ds_ref 	= U(1);
	double res_dt 	= Residual(0);
	double res_ds 	= Residual(1);
	double res_e 	= Residual(4);
    double norm_res_ds;
	double norm_res_dt;
	double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;

	double normgraddt = sqrt(gradU(0)*gradU(0) + gradU(1)*gradU(1));
	double normgradds = sqrt(gradU(2)*gradU(2) + gradU(3)*gradU(3));
	    
	t_el(0) = -a_el(1);
	t_el(1) =  a_el(0);

    res_m(0) = Residual(2); 
    res_m(1) = Residual(3);

	norm_res_ds = sqrt(res_ds*res_ds);
	norm_res_dt = sqrt(res_dt*res_dt);
    norm_res_m 	= sqrt(res_m(0)*res_m(0) + res_m(1)*res_m(1));
    norm_res_e 	= sqrt(res_e*res_e);

	for (i = 4; i < 8; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    
    
	double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 8; i < 10; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
	if (norm_grade > tol)         k_sc = 0.5*h*alpha_de*(norm_res_e/norm_grade);
	
	if (normgraddt > tol) 		  a_dt = 0.5*h*alpha_dc*norm_res_dt/normgraddt; 
	if (normgradds > tol) 		  a_ds = 0.5*h*alpha_dc*norm_res_ds/normgradds;

	
//	Isotropic DC density

	for (i = 0; i < 2; i++)       dt_diff[i] = a_dt*gradU[i];
	for (i = 0; i < 2; i++)       ds_diff[i] = a_ds*gradU[i+2];


//	Crosswind DC density

/*
	const double ads = a_el[0]*gradU[2] + a_el[1]*gradU[3];
	const double adt = a_el[0]*gradU[0] + a_el[1]*gradU[1];

	const double tds = t_el[0]*gradU[2] + t_el[1]*gradU[3];
	const double tdt = t_el[0]*gradU[0] + t_el[1]*gradU[1];

	const double k_dt_stream = 0.0;
	const double k_ds_stream = 1.0;

	const double k_dt_cross = 0.0;
	const double k_ds_cross = 0.5;

	for (i = 0; i < 2; i++){
		dt_diff[i] = k_dt_stream*a_dt*adt*a_el[i] + k_dt_cross*a_dt*tdt*t_el[i];
		ds_diff[i] = k_ds_stream*a_ds*ads*a_el[i] + k_ds_cross*a_ds*tds*t_el[i];
	}
*/

	for (i = 0; i < 4; i++)       tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 2; i++)       q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);


}


/*
void LocalMassMatrix(array_1d<double,15>& LumpedMassMatrix,const array_1d<double,3>& N, const unsigned int nodesElement)
{

    
}
*/

template<>
void CompressibleBiphaseDGNavierStokesExplicit<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}


template<>
void CompressibleBiphaseDGNavierStokesExplicit<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

template<>
void CompressibleBiphaseDGNavierStokesExplicit<2>::ComputeGaussPointLHSContribution(BoundedMatrix<double,15,15>& lhs, const ElementDataStruct& data)
{   
    
}



template<>
void CompressibleBiphaseDGNavierStokesExplicit<2>::ComputeGaussPointRHSContribution(array_1d<double,15>& rhs, const ElementDataStruct& data)
{
    const unsigned  int nodesElement = 3;
    const unsigned  int SpaceDimension = 2;
    const unsigned  int nScalarVariables = SpaceDimension + 3;
    const unsigned  int nNodalVariables = nScalarVariables*nodesElement;
    const double h = data.h;

    unsigned int i, j, k, l, m, s, t, tt, p, pp;
    
    const unsigned int size3 = nScalarVariables * SpaceDimension * nScalarVariables;
    const unsigned int size4 = nScalarVariables * SpaceDimension * nScalarVariables * nScalarVariables;
    const unsigned int sizeK = SpaceDimension * SpaceDimension * nScalarVariables * SpaceDimension;
    const unsigned int sizeKT = nScalarVariables;


    array_1d<double, size3>     A;
    array_1d<double, size4>     dAdU;
    array_1d<double, sizeK>     K;
    array_1d<double, sizeKT>    KT;

    array_1d<double,nNodalVariables>                     LumpedMassMatrix;
    array_1d<double,nScalarVariables>                    U_gauss;
    array_1d<double,nNodalVariables>                     U;
    array_1d<double,nNodalVariables>                     Un;
    array_1d<double,nNodalVariables>                     up;
    array_1d<double,nScalarVariables>                    L;
    array_1d<double,nScalarVariables>                    Residual;
    array_1d<double,nNodalVariables*nScalarVariables>    Lstar;
    array_1d<double,nScalarVariables*SpaceDimension>     G;
    array_1d<double,nScalarVariables*nScalarVariables>   S;
    array_1d<double,nScalarVariables*nScalarVariables>   B;
    array_1d<double,SpaceDimension>                      dt_diff;
	array_1d<double,SpaceDimension>                      ds_diff;
	array_1d<double,SpaceDimension*SpaceDimension>       tau;
    array_1d<double,SpaceDimension>                      q;
	array_1d<double,SpaceDimension>                      a_el;
    array_1d<double,nScalarVariables*nNodalVariables>    NN;
    array_1d<double,nScalarVariables*SpaceDimension>     gradU;
    array_1d<double,(nScalarVariables*SpaceDimension)*nNodalVariables>   gradV;
    array_1d<double,nScalarVariables>                    invtauStab;
	array_1d<double,nScalarVariables>                    switchStab;

    array_1d<double,nNodalVariables>     FConv;
    array_1d<double,nNodalVariables>     FStab;
    array_1d<double,nNodalVariables>     FDiff;
    array_1d<double,nNodalVariables>     F;
    
	
    
    const double& ctau = 0.5;   // This coefficient multiplies the divergence of the velocity in the calculation of tau. 
                                // In 3d would be 0.66667
    
    const double& dt = data.dt;

    // In this implementation this function returns only nodal forces, rgardless of the time integration scheme used.
    // const double& bdf0 = data.bdf0;
    // const double& bdf1 = data.bdf1;
    // const double& bdf2 = data.bdf2;

    const BoundedMatrix<double,nodesElement,nScalarVariables>& UU = data.U;			
    const BoundedMatrix<double,nodesElement,nScalarVariables>& UUn = data.Un;
    const BoundedMatrix<double,nodesElement,nScalarVariables>& Up = data.Up;    // Useful for the stabilizing part.
    BoundedMatrix<double,nodesElement,nScalarVariables> UUp = Up;

    const BoundedMatrix<double,nodesElement,SpaceDimension>& f_ext = data.f_ext;			
    const array_1d<double,nodesElement>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double Cv = data.c_v;
    const double gamma = data.gamma;
    const double Cp = Cv * gamma;
    const double R = Cp - Cv;
    const double ros = data.ros;
    const double Cs = data.c_s;

	const double sw_conv = 1.0;
    const double sw_diff = 1.0;
    const double sw_stab = 1.0;


    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

	switchStab[0] = 1.0;
	switchStab[1] = 1.0;
	switchStab[2] = 1.0;
	switchStab[3] = 1.0;
	switchStab[4] = 1.0;    

    // Get shape function values
    const array_1d<double,nodesElement>& N = data.N;					       
    const BoundedMatrix<double,nodesElement,SpaceDimension>& DN = data.DN_DX;	

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,SpaceDimension> f_gauss = prod(trans(f_ext), N);      
    const double r_gauss = N(0)*r(0) + N(1)*r(1) + N(2)*r(2);
    

    // Define U and Udot
    for (i = 0; i < nodesElement; i++){
        for (j = 0; j < nScalarVariables; j++){
            
    //        UUp(i,j) = (UU(i,j) - UUn(i,j))/dt;   // This could be done better, since we have an estimation of the udot from the previous substep
            
        }
    }

    for (i = 0; i < nScalarVariables; i++){
        U_gauss(i) = 0;
        for (j = 0; j < nodesElement; j++){
            U_gauss(i) += N(j)*UU(j,i);

        }
    }

    for (i = 0; i < nScalarVariables*SpaceDimension; i++)       gradU(i) = 0.0;
// Compute the gradient of the variables /dro/dx,dro/dy, dm1/dx dm1/dy dm2/dx dm2/dy de/dx de/dy)
    for (i = 0; i < nScalarVariables; i++){         // Controllare
        for (k = 0; k < nodesElement; k++){

            gradU(i*SpaceDimension + 0) += DN(k,0)*UU(k,i);
            gradU(i*SpaceDimension + 1) += DN(k,1)*UU(k,i);

        }
    }

    for (i = 0; i < nNodalVariables*nScalarVariables; i++ )  NN(i) = 0.0;

    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file
   
    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file    
    for ( i = 0; i < nScalarVariables * SpaceDimension * nNodalVariables; i++)  gradV[i] = 0.0;

	for (i = 0; i < nodesElement; i++ ){
        
        gradV[0*nNodalVariables + nScalarVariables*i    ] = DN(i,0);
        gradV[1*nNodalVariables + nScalarVariables*i    ] = DN(i,1);

        gradV[2*nNodalVariables + nScalarVariables*i + 1] = DN(i,0);
        gradV[3*nNodalVariables + nScalarVariables*i + 1] = DN(i,1);

        gradV[4*nNodalVariables + nScalarVariables*i + 2] = DN(i,0);
        gradV[5*nNodalVariables + nScalarVariables*i + 2] = DN(i,1);

        gradV[6*nNodalVariables + nScalarVariables*i + 3] = DN(i,0);
        gradV[7*nNodalVariables + nScalarVariables*i + 3] = DN(i,1);

        gradV[8*nNodalVariables + nScalarVariables*i + 4] = DN(i,0);  // Verify these lines
        gradV[9*nNodalVariables + nScalarVariables*i + 4] = DN(i,1);

	}


    const double DG_el = U_gauss(0);
    const double DS_el = U_gauss(1);
    const double m1_el = U_gauss(2);
    const double m2_el = U_gauss(3);
    const double etot_el = U_gauss(4);

    const double DTOT_el = DG_el + DS_el; 

    const double  kk = DS_el/DG_el;

    // HERE DEFINE MIXTURE COEFFICIENTS CV_M, ECC.

    const double mu_mixture = mu/(1 + kk);        // Controllare questo spaccimmo di termine (trovare su Neri et al.)

	const double Cp_mixture = (Cp + kk*Cs)/(1.0 + kk);
	const double Cv_mixture = (Cv + kk*Cs)/(1.0 + kk);

	const double Prandtl = 0.5;

//	const double lambda_mixture = Cp_mixture*mu_mixture/Prandtl;
	const double lambda_mixture = lambda;

	// Fine variazione

    const double u1_el = m1_el/DTOT_el;
	const double u2_el = m2_el/DTOT_el;

	const double norm2u = u1_el*u1_el + u2_el*u2_el;
	const double norm_u = sqrt(norm2u);
	const double norm2m = m1_el*m1_el + m2_el*m2_el;

	a_el(0) = u1_el/norm_u;
	a_el(1) = u2_el/norm_u;

	if (norm_u == 0){
		a_el(0) = 0.0;
		a_el(1) = 0.0;	 
	}	



    const double	DTOT_el2 	= DTOT_el*DTOT_el;
	const double	DTOT_el3 	= DTOT_el*DTOT_el2;
	
	const double	coeff 	= DTOT_el2*ros;
	
	const double	Dsmros 	= DS_el - ros;
	const double	Dsmros2 = Dsmros*Dsmros;
	const double	Dsmros3 = Dsmros*Dsmros2;
	
	const double	Cmixed  = (Cv*DG_el + Cs*DS_el);
	const double	Cmixed2 = Cmixed*Cmixed;
	const double	Cmixed3 = Cmixed2*Cmixed;

	const double	Cpmixed = Cs*DS_el + (Cv + R)*DG_el;
	
	const double	Cv2Dg3 	= Cv*Cv*DG_el*DG_el*DG_el;
	const double	Cs2Ds3 	= Cs*Cs*DS_el*DS_el*DS_el;
	const double	Ds3	   	= DS_el*DS_el*DS_el;
	
	
	
	
	const double	p_el   	= (DG_el*(2*DTOT_el*etot_el - norm2m)*R)/(2*DTOT_el*Cmixed);
		
	const double	pdg 	=  Cs*DS_el*etot_el*R/Cmixed2 + ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*norm2m*R)/(2*DTOT_el2*Cmixed2);
	const double	pds 	= -Cs*DG_el*etot_el*R/Cmixed2 + (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*norm2m*R)/(2*DTOT_el2*Cmixed2);
	const double	pm1 	= -(DG_el*m1_el*R)/(DTOT_el*Cmixed);
	const double	pm2 	= -(DG_el*m2_el*R)/(DTOT_el*Cmixed);
	const double	pet 	= DG_el*R/Cmixed;
	
	const double	pdgdg 	=  -(2*Cs*Cv*DS_el*etot_el*R)/Cmixed3 - ((Cv2Dg3 - Cs2Ds3 - Cs*Cv*DS_el*DS_el*(3*DG_el + DS_el))*norm2m*R)/(DTOT_el3*Cmixed3);
	const double	pdgds 	=  (Cs*(Cv*DG_el - Cs*DS_el)*etot_el*R)/Cmixed3 - ((Cv2Dg3 - Cs*Cv*DS_el*DS_el*(3*DG_el + DS_el))*norm2m*R)/(DTOT_el3*Cmixed3);
	const double	pdgm1 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m1_el*R)/(DTOT_el2*Cmixed2);
	const double	pdgm2 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m2_el*R)/(DTOT_el2*Cmixed2);
	const double	pdget 	=  Cs*DS_el*R/Cmixed2;
	
	const double	pdsds 	= (2*Cs*Cs*DG_el*etot_el*R)/Cmixed3 - (DG_el*(Cv*Cv*DG_el*DG_el + Cs*Cv*DG_el*(DG_el + 3*DS_el) + 
			Cs*Cs*(DG_el*DG_el + 3*DG_el*DS_el + 3*DS_el*DS_el))*norm2m*R)/(DTOT_el3*Cmixed3);
	const double	pdsm1 	= (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*m1_el*R)/(DTOT_el2*Cmixed2);
	const double	pdsm2 	= (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*m2_el*R)/(DTOT_el2*Cmixed2);
	const double	pdset 	= -Cs*DG_el*R/Cmixed2;
	
	const double	pm1m1 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const double	pm2m2 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const double    Temperature = p_el/(DG_el*R);

	const double gas_concentration = 1 - DS_el/ros;
	const double gas_density = DG_el/gas_concentration;

	// printf("eps_g = %.3e\n", gas_concentration);

    double SpeedSound2  = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);  // verificare sound_speed
	double SpeedSound;

	SpeedSound = sqrt(SpeedSound2);
		 

	// printf("gamma = %.3e - SpeedSound = %.3e\n", gamma, SpeedSound); 

	if (dt > h/(SpeedSound + norm_u))   printf("dt = %.3e dt_crit = %.3e\n", dt, h/(SpeedSound + norm_u));

	if (DTOT_el < 0)    				printf("dtot_el = %.3e \n", DTOT_el);

    for (i = 0; i < size3; i++)     A(i) = 0.0;
	for (i = 0; i < size4; i++)     dAdU(i) = 0.0;
			
	// Build A

	A[conta_dg(0,0,0,0,5,2,5,1)] =  DS_el*m1_el/DTOT_el2; 
	A[conta_dg(0,0,1,0,5,2,5,1)] = -DG_el*m1_el/DTOT_el2;
	A[conta_dg(0,0,2,0,5,2,5,1)] =  DG_el/DTOT_el;
	
	A[conta_dg(0,1,0,0,5,2,5,1)] =  DS_el*m2_el/DTOT_el2; 
	A[conta_dg(0,1,1,0,5,2,5,1)] = -DG_el*m2_el/DTOT_el2;
	A[conta_dg(0,1,3,0,5,2,5,1)] =  DG_el/DTOT_el;
	
	A[conta_dg(1,0,0,0,5,2,5,1)] = -DS_el*m1_el/DTOT_el2; 
	A[conta_dg(1,0,1,0,5,2,5,1)] =  DG_el*m1_el/DTOT_el2;
	A[conta_dg(1,0,2,0,5,2,5,1)] =  DS_el/DTOT_el;
	
	A[conta_dg(1,1,0,0,5,2,5,1)] = -DS_el*m2_el/DTOT_el2; 
	A[conta_dg(1,1,1,0,5,2,5,1)] =  DG_el*m2_el/DTOT_el2;
	A[conta_dg(1,1,3,0,5,2,5,1)] =  DS_el/DTOT_el;
	
	A[conta_dg(2,0,0,0,5,2,5,1)] =  -m1_el*m1_el/DTOT_el2 + pdg;
	A[conta_dg(2,0,1,0,5,2,5,1)] =  -m1_el*m1_el/DTOT_el2 + pds;
	A[conta_dg(2,0,2,0,5,2,5,1)] =  2*m1_el/DTOT_el + pm1;
	A[conta_dg(2,0,3,0,5,2,5,1)] =  pm2;
	A[conta_dg(2,0,4,0,5,2,5,1)] =  pet;
	
	A[conta_dg(2,1,0,0,5,2,5,1)] =  -m1_el*m2_el/DTOT_el2;
	A[conta_dg(2,1,1,0,5,2,5,1)] =  -m1_el*m2_el/DTOT_el2;
	A[conta_dg(2,1,2,0,5,2,5,1)] =  m2_el/DTOT_el;
	A[conta_dg(2,1,3,0,5,2,5,1)] =  m1_el/DTOT_el;
	
	A[conta_dg(3,0,0,0,5,2,5,1)] =  -m1_el*m2_el/DTOT_el2;
	A[conta_dg(3,0,1,0,5,2,5,1)] =  -m1_el*m2_el/DTOT_el2;
	A[conta_dg(3,0,2,0,5,2,5,1)] =  m2_el/DTOT_el;
	A[conta_dg(3,0,3,0,5,2,5,1)] =  m1_el/DTOT_el;
		
	A[conta_dg(3,1,0,0,5,2,5,1)] =  -m2_el*m2_el/DTOT_el2 + pdg;
	A[conta_dg(3,1,1,0,5,2,5,1)] =  -m2_el*m2_el/DTOT_el2 + pds;
	A[conta_dg(3,1,2,0,5,2,5,1)] =  pm1;
	A[conta_dg(3,1,3,0,5,2,5,1)] =  2*m2_el/DTOT_el + pm2;
	A[conta_dg(3,1,4,0,5,2,5,1)] =  pet;
	
	A[conta_dg(4,0,0,0,5,2,5,1)] =  -(m1_el*(etot_el - DTOT_el*pdg + p_el))/DTOT_el2;
	A[conta_dg(4,0,1,0,5,2,5,1)] =  -(m1_el*(etot_el - DTOT_el*pds + p_el))/DTOT_el2;
	A[conta_dg(4,0,2,0,5,2,5,1)] =  (etot_el + m1_el*pm1 + p_el)/DTOT_el;
	A[conta_dg(4,0,3,0,5,2,5,1)] =  m1_el*pm2/DTOT_el;
	A[conta_dg(4,0,4,0,5,2,5,1)] =	m1_el*(1.0 + pet)/DTOT_el;
	
	A[conta_dg(4,1,0,0,5,2,5,1)] =  -(m2_el*(etot_el - DTOT_el*pdg + p_el))/DTOT_el2;
	A[conta_dg(4,1,1,0,5,2,5,1)] =  -(m2_el*(etot_el - DTOT_el*pds + p_el))/DTOT_el2;
	A[conta_dg(4,1,2,0,5,2,5,1)] =	m2_el*pm1/DTOT_el;  
	A[conta_dg(4,1,3,0,5,2,5,1)] =  (etot_el + m2_el*pm2 + p_el)/DTOT_el;
	A[conta_dg(4,1,4,0,5,2,5,1)] =  m2_el*(1.0 + pet)/DTOT_el;



    //	Build dAdU

dAdU[conta_dg(0,0,0,0,5,2,5,5)] = -2*DS_el*m1_el/DTOT_el3;	 
	dAdU[conta_dg(0,0,0,1,5,2,5,5)] =  (DG_el - DS_el)*m1_el/DTOT_el3;
	dAdU[conta_dg(0,0,0,2,5,2,5,5)] = DS_el/DTOT_el2;
		
	dAdU[conta_dg(0,0,1,0,5,2,5,5)] = (DG_el - DS_el)*m1_el/DTOT_el3; 
	dAdU[conta_dg(0,0,1,1,5,2,5,5)] = 2*DG_el*m1_el/DTOT_el3;
	dAdU[conta_dg(0,0,1,2,5,2,5,5)] = -DG_el/DTOT_el2;
	 
	dAdU[conta_dg(0,0,2,0,5,2,5,5)] = DS_el/DTOT_el2; 
	dAdU[conta_dg(0,0,2,1,5,2,5,5)] = -DG_el/DTOT_el2;
			
	/////////////////////////////////////////////////////////////
	
	dAdU[conta_dg(0,1,0,0,5,2,5,5)] = -2*DS_el*m2_el/DTOT_el3; 
	dAdU[conta_dg(0,1,0,1,5,2,5,5)] =  (DG_el - DS_el)*m2_el/DTOT_el3;
	dAdU[conta_dg(0,1,0,3,5,2,5,5)] = DS_el/DTOT_el2;
		
	dAdU[conta_dg(0,1,1,0,5,2,5,5)] =  (DG_el - DS_el)*m2_el/DTOT_el3;
	dAdU[conta_dg(0,1,1,1,5,2,5,5)] =  2*DG_el*m2_el/DTOT_el3;
	dAdU[conta_dg(0,1,1,3,5,2,5,5)] =  -DG_el/DTOT_el2; 
		 
	dAdU[conta_dg(0,1,3,0,5,2,5,5)] = DS_el/DTOT_el2; 
	dAdU[conta_dg(0,1,3,1,5,2,5,5)] = -DG_el/DTOT_el2;
	
// 1  ////////////////////////////////////////////////////////
			
	dAdU[conta_dg(1,0,0,0,5,2,5,5)] = 2*DS_el*m1_el/DTOT_el3; 
	dAdU[conta_dg(1,0,0,1,5,2,5,5)] = (-DG_el + DS_el)*m1_el/DTOT_el3;
	dAdU[conta_dg(1,0,0,2,5,2,5,5)] = -DS_el/DTOT_el2;
	
	dAdU[conta_dg(1,0,1,0,5,2,5,5)] = (-DG_el + DS_el)*m1_el/DTOT_el3;  
	dAdU[conta_dg(1,0,1,1,5,2,5,5)] = -2*DG_el*m1_el/DTOT_el3;
	dAdU[conta_dg(1,0,1,2,5,2,5,5)] = DG_el/DTOT_el2;
	 
	dAdU[conta_dg(1,0,2,0,5,2,5,5)] = -DS_el/DTOT_el2; 
	dAdU[conta_dg(1,0,2,1,5,2,5,5)] =  DG_el/DTOT_el2;
	
/////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_dg(1,1,0,0,5,2,5,5)] =  2*DS_el*m2_el/DTOT_el3;
	dAdU[conta_dg(1,1,0,1,5,2,5,5)] =  (-DG_el + DS_el)*m2_el/DTOT_el3;
	dAdU[conta_dg(1,1,0,3,5,2,5,5)] = -DS_el/DTOT_el2; 
	
	dAdU[conta_dg(1,1,1,0,5,2,5,5)] =  (-DG_el + DS_el)*m2_el/DTOT_el3;
	dAdU[conta_dg(1,1,1,1,5,2,5,5)] = -2*DG_el*m2_el/DTOT_el3;
	dAdU[conta_dg(1,1,1,3,5,2,5,5)] = DG_el/DTOT_el2; 
	 
	dAdU[conta_dg(1,1,3,0,5,2,5,5)] = -DS_el/DTOT_el2;
	dAdU[conta_dg(1,1,3,1,5,2,5,5)] =  DG_el/DTOT_el2;
	
			
// 2  /////////////////////////////////////////////////////////////
			
	dAdU[conta_dg(2,0,0,0,5,2,5,5)] = -2*m1_el*m1_el/DTOT_el3 + pdgdg;	 
	dAdU[conta_dg(2,0,0,1,5,2,5,5)] = -2*m1_el*m1_el/DTOT_el3 + pdgds;
	dAdU[conta_dg(2,0,0,2,5,2,5,5)] = -2*m1_el/DTOT_el2 + pdgm1;
	dAdU[conta_dg(2,0,0,3,5,2,5,5)] =  pdgm2;
	dAdU[conta_dg(2,0,0,4,5,2,5,5)] =  pdget;
	
	dAdU[conta_dg(2,0,1,0,5,2,5,5)] = 2*m1_el*m1_el/DTOT_el3 + pdgds;
	dAdU[conta_dg(2,0,1,1,5,2,5,5)] = 2*m1_el*m1_el/DTOT_el3 + pdsds;
	dAdU[conta_dg(2,0,1,2,5,2,5,5)] = -2*m1_el/DTOT_el2 + pdsm1;
	dAdU[conta_dg(2,0,1,3,5,2,5,5)] = pdsm2;
	dAdU[conta_dg(2,0,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta_dg(2,0,2,0,5,2,5,5)] = -2*m1_el/DTOT_el2 + pdgm1;
	dAdU[conta_dg(2,0,2,1,5,2,5,5)] = -2*m1_el/DTOT_el2 + pdsm1;
	dAdU[conta_dg(2,0,2,2,5,2,5,5)] = 2.0/DTOT_el + pm1m1;
	 
	dAdU[conta_dg(2,0,3,0,5,2,5,5)] = pdgm2; 
	dAdU[conta_dg(2,0,3,1,5,2,5,5)] = pdsm2;
	dAdU[conta_dg(2,0,3,3,5,2,5,5)] = pm2m2; 
				
	dAdU[conta_dg(2,0,4,0,5,2,5,5)] = pdget;
	dAdU[conta_dg(2,0,4,1,5,2,5,5)] = pdset;
	 
	dAdU[conta_dg(2,1,0,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;	 
	dAdU[conta_dg(2,1,0,1,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;
	dAdU[conta_dg(2,1,0,2,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(2,1,0,3,5,2,5,5)] = -m1_el/DTOT_el2;
	
	dAdU[conta_dg(2,1,1,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3; 
	dAdU[conta_dg(2,1,1,1,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;
	dAdU[conta_dg(2,1,1,2,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(2,1,1,3,5,2,5,5)] = -m1_el/DTOT_el2;
	 
	dAdU[conta_dg(2,1,2,0,5,2,5,5)] = -m2_el/DTOT_el2; 
	dAdU[conta_dg(2,1,2,1,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(2,1,2,3,5,2,5,5)] = 1.0/DTOT_el; 
			
	dAdU[conta_dg(2,1,3,0,5,2,5,5)] = -m1_el/DTOT_el2; 
	dAdU[conta_dg(2,1,3,1,5,2,5,5)] = -m1_el/DTOT_el2;
	dAdU[conta_dg(2,1,3,2,5,2,5,5)] = 1.0/DTOT_el;
	
			
// 3 /////////////////////////////////////////////////////////////////////////////
			
	dAdU[conta_dg(3,0,0,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;	 
	dAdU[conta_dg(3,0,0,1,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;
	dAdU[conta_dg(3,0,0,2,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(3,0,0,3,5,2,5,5)] = -m1_el/DTOT_el2;
		
	dAdU[conta_dg(3,0,1,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3; 
	dAdU[conta_dg(3,0,1,1,5,2,5,5)] = 2*m1_el*m2_el/DTOT_el3;
	dAdU[conta_dg(3,0,1,2,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(3,0,1,3,5,2,5,5)] = -m1_el/DTOT_el2;
	 
	dAdU[conta_dg(3,0,2,0,5,2,5,5)] = -m2_el/DTOT_el2; 
	dAdU[conta_dg(3,0,2,1,5,2,5,5)] = -m2_el/DTOT_el2;
	dAdU[conta_dg(3,0,2,3,5,2,5,5)] = 1.0/DTOT_el; 
			
	dAdU[conta_dg(3,0,3,0,5,2,5,5)] = -m1_el/DTOT_el2; 
	dAdU[conta_dg(3,0,3,1,5,2,5,5)] = -m1_el/DTOT_el2;
	dAdU[conta_dg(3,0,3,2,5,2,5,5)] = 1.0/DTOT_el;
	
//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_dg(3,1,0,0,5,2,5,5)] = 2*m2_el*m2_el/DTOT_el3 + pdgdg;	 
	dAdU[conta_dg(3,1,0,1,5,2,5,5)] = 2*m2_el*m2_el/DTOT_el3 + pdgds;
	dAdU[conta_dg(3,1,0,2,5,2,5,5)] = pdgm1;
	dAdU[conta_dg(3,1,0,3,5,2,5,5)] = -2*m2_el/DTOT_el2 + pdgm2;
	dAdU[conta_dg(3,1,0,4,5,2,5,5)] = pdget;
	
	dAdU[conta_dg(3,1,1,0,5,2,5,5)] = 2*m2_el*m2_el/DTOT_el3 + pdgds;
	dAdU[conta_dg(3,1,1,1,5,2,5,5)] = 2*m2_el*m2_el/DTOT_el3 + pdsds;
	dAdU[conta_dg(3,1,1,2,5,2,5,5)] = pdsm1;
	dAdU[conta_dg(3,1,1,3,5,2,5,5)] = -2*m2_el/DTOT_el2 + pdsm2;
	dAdU[conta_dg(3,1,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta_dg(3,1,2,0,5,2,5,5)] = pdgm1;
	dAdU[conta_dg(3,1,2,1,5,2,5,5)] = pdsm1;
	dAdU[conta_dg(3,1,2,2,5,2,5,5)] = pm1m1;
	
	dAdU[conta_dg(3,1,3,0,5,2,5,5)] = -2*m2_el/DTOT_el2 + pdgm2;
	dAdU[conta_dg(3,1,3,1,5,2,5,5)] = -2*m2_el/DTOT_el2 + pdsm2;
	dAdU[conta_dg(3,1,3,3,5,2,5,5)] = 2.0/DTOT_el + pm2m2; 
	
	dAdU[conta_dg(3,1,4,0,5,2,5,5)] = pdget;
	dAdU[conta_dg(3,1,4,1,5,2,5,5)] = pdset;
			
// 4 //////////////////////////////////////////////////////////////////////
			
	dAdU[conta_dg(4,0,0,0,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - 2*pdg*DTOT_el + pdgdg*DTOT_el2)/DTOT_el3;
	dAdU[conta_dg(4,0,0,1,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - pdg*DTOT_el - pds*DTOT_el + pdgds*DTOT_el2)/DTOT_el3;
	dAdU[conta_dg(4,0,0,2,5,2,5,5)] = (-etot_el + DS_el*(pdg + m1_el*pdgm1) - m1_el*pm1 - p_el + (pdg + m1_el*pdgm1)*DG_el)/DTOT_el2;
	dAdU[conta_dg(4,0,0,3,5,2,5,5)] = (m1_el*(DS_el*pdgm2 - pm2 + pdgm2*DG_el))/DTOT_el2;
	dAdU[conta_dg(4,0,0,4,5,2,5,5)] = (m1_el*(-1 + DS_el*pdget - pet + pdget*DG_el))/DTOT_el2;
	
	dAdU[conta_dg(4,0,1,0,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - pdg*DTOT_el - pds*DTOT_el + pdgds*DTOT_el2)/DTOT_el3;
	dAdU[conta_dg(4,0,1,1,5,2,5,5)] = (m1_el*(2*(etot_el + p_el) - 2*pds*DTOT_el + pdsds*DTOT_el2))/DTOT_el3;
	dAdU[conta_dg(4,0,1,2,5,2,5,5)] = (-etot_el + DS_el*(pds + m1_el*pdsm1) - m1_el*pm1 - p_el + (pds + m1_el*pdsm1)*DG_el)/DTOT_el2;
	dAdU[conta_dg(4,0,1,3,5,2,5,5)] = m1_el*(pdsm2/DTOT_el - pm2/DTOT_el2);
	dAdU[conta_dg(4,0,1,4,5,2,5,5)] = (m1_el*(-1 + DTOT_el*pdset - pet))/DTOT_el2;
	 
	dAdU[conta_dg(4,0,2,0,5,2,5,5)] = (-etot_el + DS_el*(pdg + m1_el*pdgm1) - m1_el*pm1 - p_el + (pdg + m1_el*pdgm1)*DG_el)/DTOT_el2;
	dAdU[conta_dg(4,0,2,1,5,2,5,5)] = (-etot_el + DS_el*(pds + m1_el*pdsm1) - m1_el*pm1 - p_el + (pds + m1_el*pdsm1)*DG_el)/DTOT_el2;
	dAdU[conta_dg(4,0,2,2,5,2,5,5)] = (m1_el*pm1m1 + 2*pm1)/DTOT_el;
	dAdU[conta_dg(4,0,2,3,5,2,5,5)] = pm2/DTOT_el;
	dAdU[conta_dg(4,0,2,4,5,2,5,5)] = (1 + pet)/DTOT_el;
			
	dAdU[conta_dg(4,0,3,0,5,2,5,5)] = (m1_el*(DS_el*pdgm2 - pm2 + pdgm2*DG_el))/DTOT_el2;
	dAdU[conta_dg(4,0,3,1,5,2,5,5)] = m1_el*(pdsm2/DTOT_el - pm2/DTOT_el2);
	dAdU[conta_dg(4,0,3,2,5,2,5,5)] = pm2/DTOT_el;
	dAdU[conta_dg(4,0,3,3,5,2,5,5)] = m1_el*pm2m2/DTOT_el;
			
	dAdU[conta_dg(4,0,4,0,5,2,5,5)] = (m1_el*(-1 + DS_el*pdget - pet + pdget*DG_el))/DTOT_el2;
	dAdU[conta_dg(4,0,4,1,5,2,5,5)] = (m1_el*(-1 + DTOT_el*pdset - pet))/DTOT_el2;
	dAdU[conta_dg(4,0,4,2,5,2,5,5)] = (1 + pet)/DTOT_el;
	 
	dAdU[conta_dg(4,1,0,0,5,2,5,5)] = m2_el*(2*(etot_el + p_el) - 2*pdg*DTOT_el + pdgdg*DTOT_el2)/DTOT_el3;	 
	dAdU[conta_dg(4,1,0,1,5,2,5,5)] = (m2_el*pdgds)/DTOT_el - (m2_el*(pdg + pds))/DTOT_el2 + (2*m2_el*(etot_el + p_el))/DTOT_el3;
	dAdU[conta_dg(4,1,0,2,5,2,5,5)] = (DTOT_el*m2_el*pdgm1 - m2_el*pm1)/DTOT_el2; 
	dAdU[conta_dg(4,1,0,3,5,2,5,5)] = (-etot_el - p_el + DTOT_el*(pdg + m2_el*pdgm2) - m2_el*pm2)/DTOT_el2;
	dAdU[conta_dg(4,1,0,4,5,2,5,5)] = (m2_el*(-1.0 + DTOT_el*pdget - pet))/DTOT_el2;
	
	dAdU[conta_dg(4,1,1,0,5,2,5,5)] = (m2_el*pdgds)/DTOT_el - (m2_el*(pdg + pds))/DTOT_el2 + (2*m2_el*(etot_el + p_el))/DTOT_el3;
	dAdU[conta_dg(4,1,1,1,5,2,5,5)] = m2_el*(-2*pds/DTOT_el2 + pdsds/DTOT_el + 2*(etot_el + p_el)/DTOT_el3);
	dAdU[conta_dg(4,1,1,2,5,2,5,5)] = m2_el*(DTOT_el*pdsm1 - pm1)/DTOT_el2;
	dAdU[conta_dg(4,1,1,3,5,2,5,5)] = (-etot_el + DTOT_el*(pds + m2_el*pdsm2) - m2_el*pm2 - p_el)/DTOT_el2;
	dAdU[conta_dg(4,1,1,4,5,2,5,5)] = (m2_el*(-1 + DTOT_el*pdset - pet))/DTOT_el2;
	 
	dAdU[conta_dg(4,1,2,0,5,2,5,5)] = (DTOT_el*m2_el*pdgm1 - m2_el*pm1)/DTOT_el2;
	dAdU[conta_dg(4,1,2,1,5,2,5,5)] = m2_el*(DTOT_el*pdsm1 - pm1)/DTOT_el2;
	dAdU[conta_dg(4,1,2,2,5,2,5,5)] = m2_el*pm1m1/DTOT_el;
	dAdU[conta_dg(4,1,2,3,5,2,5,5)] = pm1/DTOT_el;
			
	dAdU[conta_dg(4,1,3,0,5,2,5,5)] = (-etot_el - p_el + DTOT_el*(pdg + m2_el*pdgm2) - m2_el*pm2)/DTOT_el2;
	dAdU[conta_dg(4,1,3,1,5,2,5,5)] = (-etot_el + DTOT_el*(pds + m2_el*pdsm2) - m2_el*pm2 - p_el)/DTOT_el2;  
	dAdU[conta_dg(4,1,3,2,5,2,5,5)] = pm1/DTOT_el;
	dAdU[conta_dg(4,1,3,3,5,2,5,5)] = (m2_el*pm2m2 + 2*pm2)/DTOT_el;
	dAdU[conta_dg(4,1,3,4,5,2,5,5)] = (1 + pet)/DTOT_el;
			
	dAdU[conta_dg(4,1,4,0,5,2,5,5)] = (m2_el*(-1.0 + DTOT_el*pdget - pet))/DTOT_el2;
	dAdU[conta_dg(4,1,4,1,5,2,5,5)] = (m2_el*(-1 + DTOT_el*pdset - pet))/DTOT_el2;
	dAdU[conta_dg(4,1,4,3,5,2,5,5)] = (1 + pet)/DTOT_el;

    for (i = 0; i < nScalarVariables*nScalarVariables; i++)   S(i) = 0.0;
																				// CONTROLLARE S					
    S(2*nScalarVariables + 0) = f_gauss(0);
	S(2*nScalarVariables + 1) = f_gauss(0);
    S(3*nScalarVariables + 0) = f_gauss(1);
	S(3*nScalarVariables + 1) = f_gauss(1);
    S(4*nScalarVariables + 0) = r_gauss;
	S(4*nScalarVariables + 1) = r_gauss;
    S(4*nScalarVariables + 2) = f_gauss(0);
    S(4*nScalarVariables + 3) = f_gauss(1);

    for (i = 0; i < nodesElement*nScalarVariables*nScalarVariables; i++)     Lstar[i] = 0.0;

    for (i = 0; i < nodesElement; i++){
			
		pp = i*nScalarVariables*nScalarVariables;
		
		Lstar[pp + 2*nScalarVariables + 0] = N[i]*S[2*nScalarVariables + 0];
		Lstar[pp + 3*nScalarVariables + 0] = N[i]*S[3*nScalarVariables + 0];

		Lstar[pp + 2*nScalarVariables + 1] = N[i]*S[2*nScalarVariables + 1];	// Controllare
		Lstar[pp + 3*nScalarVariables + 1] = N[i]*S[3*nScalarVariables + 1];	// Controllare
				
		for (k = 0; k < nScalarVariables - 1; k++){
			Lstar[pp + 4*nScalarVariables + k] = N[i]*S[4*nScalarVariables + k];
		}
    }

    for (k = 0; k < nScalarVariables; k++){
		for ( s = 0; s < nNodalVariables; s++){
			
			p = s*nScalarVariables + k;
			
			for (i = 0; i < nScalarVariables; i++){
				
				pp = i*nNodalVariables + s;

				for (j = 0; j < SpaceDimension; j++){

					t = conta_dg(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

					Lstar[p] += (-A[t]*gradV[(SpaceDimension*i+j)*nNodalVariables + s]);
                    
	
				}
			}
		}
	}

    for (i = 0; i < nScalarVariables*nScalarVariables; i++) B[i] = 0.0;

    for (i = 0; i < nScalarVariables; i++){
        for (j = 0; j < SpaceDimension; j++){
            for (k = 0; k < nScalarVariables; k++){
                pp = i * nScalarVariables + k;
                for (m = 0; m < nScalarVariables; m++){
                    tt = conta_dg(i,j,k,m,nScalarVariables,SpaceDimension,nScalarVariables,nScalarVariables);

                    B[pp]  += dAdU[tt]*gradU[m*SpaceDimension + j];

                }
            }
        }
	}

    for (j = 0; j < nodesElement; j++){
        pp = j*nScalarVariables*nScalarVariables;
        for (i = 0; i < nScalarVariables; i++){
            p = i*nScalarVariables;
            for (k = 0; k < nScalarVariables; k++){
                Lstar[pp + p + k] -= N[j]*B[p + k];
            }
        }
    }


	invtauStab[0] =	stab_c2*(norm_u + SpeedSound)/h; 
    invtauStab[1] =	stab_c2*(norm_u + SpeedSound)/h;
	invtauStab[2] =	stab_c1*mu_mixture/(DTOT_el*h*h) + invtauStab[0];
	invtauStab[3] =	invtauStab[2];
	invtauStab[4] = stab_c1*lambda_mixture/(DTOT_el*Cp_mixture*h*h) + invtauStab[0];  

	
// controllare L
    L[0] = 0.0;
    L[1] = 0.0;
    L[2] = -S[2*nScalarVariables + 0]*U_gauss[0] - S[2*nScalarVariables + 1]*U_gauss[1];
    L[3] = -S[3*nScalarVariables + 0]*U_gauss[0] - S[3*nScalarVariables + 1]*U_gauss[1];
    L[4] = 0.0;
    

    for (k = 0; k < nScalarVariables - 1 ; k++){
        L[4] -= S[4*nScalarVariables + k]*U_gauss[k];
    }

    for (i = 0; i < nScalarVariables; i++ ){
		for (k = 0; k < nScalarVariables; k++){
			for (j = 0; j < SpaceDimension; j++){

				s = conta_dg(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

				L[i] += A[s]*gradU[k*SpaceDimension + j];
			}
		}

		Residual[i] = -L[i];

		for (k = 0; k < nodesElement; k++){
            Residual[i] -= N[k]*UUp(k,i);
        }
	}

	for (i = 0; i < nScalarVariables; i++){
        for (k = 0; k < nodesElement; k++){
			FConv[i + k*nScalarVariables] = N[k]*L[i];
		}
    }

    // Build diffusive term: stress tensor and thermal diffusion

    
    // HERE FiND TAU AND q

    for (i = 0; i < sizeK;  i++)    K[i]  = 0.0;
    for (i = 0; i < sizeKT; i++)    KT[i] = 0.0;
			
	K[conta_dg(0,0,0,0,2,2,5,2)] =  (-2.0 + ctau)*m1_el/DTOT_el2;
	K[conta_dg(0,0,0,1,2,2,5,2)] =  ctau*m2_el/DTOT_el2;
	K[conta_dg(0,0,1,0,2,2,5,2)] =  (-2.0 + ctau)*m1_el/DTOT_el2;
	K[conta_dg(0,0,1,1,2,2,5,2)] =  ctau*m2_el/DTOT_el2;
	K[conta_dg(0,0,2,0,2,2,5,2)] =  (2.0 - ctau)/DTOT_el;
	K[conta_dg(0,0,3,1,2,2,5,2)] = -ctau/DTOT_el;  
	
	K[conta_dg(0,1,0,0,2,2,5,2)] = -m2_el/DTOT_el2;
	K[conta_dg(0,1,0,1,2,2,5,2)] = -m1_el/DTOT_el2;
	K[conta_dg(0,1,1,0,2,2,5,2)] = -m2_el/DTOT_el2;
	K[conta_dg(0,1,1,1,2,2,5,2)] = -m1_el/DTOT_el2;  
	K[conta_dg(0,1,2,1,2,2,5,2)] =  1.0/DTOT_el; 
	K[conta_dg(0,1,3,0,2,2,5,2)] =  1.0/DTOT_el;
	
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	K[conta_dg(1,0,0,0,2,2,5,2)] = -m2_el/DTOT_el2;
	K[conta_dg(1,0,0,1,2,2,5,2)] = -m1_el/DTOT_el2;
	K[conta_dg(1,0,1,0,2,2,5,2)] = -m2_el/DTOT_el2;
	K[conta_dg(1,0,1,1,2,2,5,2)] = -m1_el/DTOT_el2;
	K[conta_dg(1,0,2,1,2,2,5,2)] =  1.0/DTOT_el;
	K[conta_dg(1,0,3,0,2,2,5,2)] =  1.0/DTOT_el; 
	
	K[conta_dg(1,1,0,0,2,2,5,2)] =  ctau*m1_el/DTOT_el2;
	K[conta_dg(1,1,0,1,2,2,5,2)] =  (ctau - 2.0)*m2_el/DTOT_el2;
	K[conta_dg(1,1,1,0,2,2,5,2)] =  ctau*m1_el/DTOT_el2;
	K[conta_dg(1,1,1,1,2,2,5,2)] =  (ctau - 2.0)*m2_el/DTOT_el2;
	K[conta_dg(1,1,2,0,2,2,5,2)] = -ctau/DTOT_el;
	K[conta_dg(1,1,3,1,2,2,5,2)] =  (2.0 - ctau)/DTOT_el;


	
	KT[0] =  pdg/(DG_el*R) - p_el/(DG_el*DG_el*R);
	KT[1] =  pds/(DG_el*R);
	KT[2] =  pm1/(DG_el*R);
	KT[3] =  pm2/(DG_el*R);
	KT[4] =  pet/(DG_el*R);


    for (i = 0; i < SpaceDimension; i++){
        
        q(i) = 0.0;
		dt_diff(i) = 0.0;
		ds_diff(i) = 0.0;

        for (j = 0; j < SpaceDimension; j++){

            tau(i*SpaceDimension + j) = 0.0;
            
            for (l = 0; l < nScalarVariables; l++){
                for (m = 0; m < SpaceDimension; m++){

                    tau(i*SpaceDimension + j) += mu_mixture*K[conta_dg(i,j,l,m,2,2,5,2)]*gradU[l*SpaceDimension + m];

                }
            }
        }
        for (l = 0; l < nScalarVariables; l++){
            q(i) -= lambda_mixture*KT(l)*gradU(l*SpaceDimension + i);
        }
    }

    ShockCapturing2dg(mu_mixture, lambda_mixture, Cv_mixture, ros, h, dt_diff, ds_diff,tau, q, DTOT_el, gradU, Residual,U_gauss,a_el,norm_u,SpeedSound);

    // Build diffusive term: Diffusion tensor

	for ( i = 0; i < nScalarVariables*SpaceDimension; i++ )    G[i] = 0.0;
		
	for (i = 0; i < SpaceDimension; i++){

		for (j = 0; j < SpaceDimension; j++)
			G[(i + 2)*SpaceDimension + j] = -tau[i*SpaceDimension + j];
    }

	for (j = 0; j < SpaceDimension; j++){

		G[0*SpaceDimension + j] = -dt_diff[j];
		G[1*SpaceDimension + j] = -ds_diff[j];   // Numerical diffusivity for dust
		G[4*SpaceDimension + j] = q[j];

		for (i = 0; i < SpaceDimension; i++)
            G[4*SpaceDimension + j] += (-U_gauss[i + 2]/DTOT_el*tau[i*SpaceDimension + j]);

	}


    // Build diffusive term: Diffusion force

	for (s = 0; s < nNodalVariables; s++){

		FDiff[s] = 0.0;

		for (i = 0; i < nScalarVariables; i++){
			for ( j = 0; j < SpaceDimension; j++){
				FDiff[s] -= gradV[(SpaceDimension*i + j)*nNodalVariables + s]*G[i*SpaceDimension + j];
            }
        }
 	}

    // Stabilizing residual part

	for (s = 0; s < nNodalVariables; s++){

		FStab[s] = 0.0;

		for (k = 0; k < nScalarVariables; k++){
			FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k]*switchStab[k];
		}
	}

    // Force contribuution at the Gauss Point   

	int check = 1;

    for (i = 0; i < nNodalVariables; i++){

		rhs[i] = - sw_conv*FConv[i] - sw_diff*FDiff[i] - sw_stab*FStab[i];

		if (std::isnan(FConv[i]) == 1 || std::isnan(FDiff[i]) == 1 || std::isnan(FStab[i]) == 1){
			printf("%d %.3e %.3e %.3e\n", i, FConv[i], FDiff[i], FStab[i]);
			check = 0;
		}

	}
	if (check == 0)	{
		printf("%.3e %.3e %.3e %.3e %.3e \n", U_gauss(0), U_gauss(1), U_gauss(2), U_gauss(3), U_gauss(4));

		for (s = 0; s < nNodalVariables; s++){

			FStab[s] = 0.0;

			for (k = 0; k < nScalarVariables; k++){
				FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k]*switchStab[k];

				printf("%d %d %.3e %.3e %.3e\n", s, k, Lstar[s*nScalarVariables + k], Residual[k], invtauStab[k]);
			}
		}

		printf("stab_c2 = %.3e - norm_u = %.3e - SpeedSound = %.3e - stab_c1 = %.3e - mu_mixture = %.3e - lambda_mixture = %.3e"
			"Cp_mixture = %.3e \n", stab_c2, norm_u, SpeedSound, stab_c1, mu_mixture, lambda_mixture, Cp_mixture);

		double sps2 = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);
		
		printf("spsound2 = %.3e\n",sps2);

		printf("DG_el = %.3e - Cpmixed = %.3e - en-u2 = %.3e - DTOT_el = %.3e - Cmixed2 = %.3e\n",
		DG_el, Cpmixed, 2*etot_el - DTOT_el*norm2u, DTOT_el, Cmixed2);
		
		printf("%.3e %.3e \n",2*etot_el, DTOT_el*norm2u);

		printf("cos cos cos \n\n");
		


		abort();
	}

	
	
}


}
