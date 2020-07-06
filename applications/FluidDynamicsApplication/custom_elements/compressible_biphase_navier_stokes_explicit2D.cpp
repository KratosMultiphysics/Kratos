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

#include "custom_elements/compressible_biphase_navier_stokes_explicit.h"

namespace Kratos {

unsigned int	conta(
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

void ShockCapturing(const double mu,
               const double lambda,
               const double c_v,
               const double h,
               array_1d<double,4>& tau,
               array_1d<double,2>& q,
               double ro_el,
               array_1d<double,10>& gradU, // 5*2
               array_1d<double,5>& Residual)
{
    const int SpaceDimension = 2;
    
   const double alpha = 5.0;                               // Algorithm constant
   const double tol = 1e-16;                               

    unsigned int i;

    double v_sc = 0.0;
    double k_sc = 0.0;
                                          //Shock capturing viscosity
    array_1d<double,SpaceDimension> res_m;
    double res_e = Residual(4);
    double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;
    

    res_m(0) = Residual(2); 
    res_m(1) = Residual(3);

    norm_res_m = sqrt(res_m(0) * res_m(0) + res_m(1) * res_m(1));
    
    for (i = 4; i < 8; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    norm_res_e = sqrt(res_e*res_e);

    double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 8; i < 10; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
 
    if (norm_grade > tol)         k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);

    for (i = 0; i < 4; i++)       tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 2; i++)       q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);

//    printf("%.3e %.3e %.3e %.3e %.3e \n",res_m(0),res_m(1),res_e,v_sc,k_sc);
    
}

void LocalMassMatrix(array_1d<double,15>& LumpedMassMatrix,const array_1d<double,3>& N, const unsigned int nodesElement)
{

    
}


template<>
void CompressibleBiphaseNavierStokesExplicit<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}


template<>
void CompressibleBiphaseNavierStokesExplicit<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

template<>
void CompressibleBiphaseNavierStokesExplicit<2>::ComputeGaussPointLHSContribution(BoundedMatrix<double,15,15>& lhs, const ElementDataStruct& data)
{   
    
}



template<>
void CompressibleBiphaseNavierStokesExplicit<2>::ComputeGaussPointRHSContribution(array_1d<double,15>& rhs, const ElementDataStruct& data)
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
    array_1d<double,SpaceDimension*SpaceDimension>       tau;
    array_1d<double,SpaceDimension>                      q;
    array_1d<double,nScalarVariables*nNodalVariables>    NN;
    array_1d<double,nScalarVariables*SpaceDimension>     gradU;
    array_1d<double,(nScalarVariables*SpaceDimension)*nNodalVariables>   gradV;
    array_1d<double,nScalarVariables>                    invtauStab;

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
    BoundedMatrix<double,nodesElement,nScalarVariables> UUp;

    const BoundedMatrix<double,nodesElement,SpaceDimension>& f_ext = data.f_ext;			
    const array_1d<double,nodesElement>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double Cv = data.c_v;
    const double gamma = data.gamma;
    const double Cp = Cv*gamma;
    const double R = Cp - Cv;
    const double ros = data.ros;
    const double Cs = data.c_s;

	const double sw_conv = 1.0;
    const double sw_diff = 1.0;
    const double sw_stab = 1.0;


    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;    

    // Get shape function values
    const array_1d<double,nodesElement>& N = data.N;					       
    const BoundedMatrix<double,nodesElement,SpaceDimension>& DN = data.DN_DX;	

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,SpaceDimension> f_gauss = prod(trans(f_ext), N);      
    const double r_gauss = N(0)*r(0) + N(1)*r(1) + N(2)*r(2);
    

    // Define U and Udot
    for (i = 0; i < nodesElement; i++){
        for (j = 0; j < nScalarVariables; j++){
            
            UUp(i,j) = (UU(i,j) - UUn(i,j))/dt;   // This could be done better, since we have an estimation of the udot from the previous substep
            
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

    const double rom = DG_el + DS_el; 

    const double  kk = DS_el/DG_el;

    // HERE DEFINE MIXTURE COEFFICIENTS CV_M, ECC.

    const double mu_mixture = mu/(1 + kk);        // Controllare questo spaccimmo di termine (trovare su Neri et al.)

	const double Cp_mixture = (Cp + kk*Cs)/(1.0 + kk);
	const double Cv_mixture = (Cv + kk*Cs)/(1.0 + kk);

	const double Prandtl = 0.5;

	const double lambda_mixture = lambda; // Cp_mixture*mu_mixture/Prandtl;

	// Fine variazione

    const double u1_el = m1_el/rom;
	const double u2_el = m2_el/rom;

    const double norm2u = u1_el*u1_el + u2_el*u2_el;
	const double norm_u = sqrt(norm2u);
	const double norm2m = m1_el*m1_el + m2_el*m2_el;

    const double rom2 	= rom*rom;
	const double rom3 	= rom*rom2;
	
	const double Cmixed  = (Cv*DG_el + Cs*DS_el);
	const double Cmixed2 = Cmixed*Cmixed;
	const double Cmixed3 = Cmixed*Cmixed2;
	
	const double Cv2Dg3 = Cv*Cv*DG_el*DG_el*DG_el;
	const double Cs2Ds3 = Cs*Cs*DS_el*DS_el*DS_el;
		
	const double p_el   = (DG_el*(2*rom*etot_el - norm2m)*R)/(2*rom*Cmixed);
		
	const double pdg 	=  Cs*DS_el*etot_el*R/Cmixed2 + ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*norm2m*R)/(2*rom2*Cmixed2);
	const double pds 	= -Cs*DG_el*etot_el*R/Cmixed2 + (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*norm2m*R)/(2*rom2*Cmixed2);
	const double pm1 	= -(DG_el*m1_el*R)/(rom*Cmixed);
	const double pm2 	= -(DG_el*m2_el*R)/(rom*Cmixed);
	const double pet 	= DG_el*R/Cmixed;
	
	const double pdgdg 	=  -(2*Cs*Cv*DS_el*etot_el*R)/Cmixed3 - ((Cv2Dg3 - Cs2Ds3 - Cs*Cv*DS_el*DS_el*(3*DG_el + DS_el))*norm2m*R)/(rom3*Cmixed3);
	const double pdgds 	=  (Cs*(Cv*DG_el - Cs*DS_el)*etot_el*R)/Cmixed3 - ((Cv2Dg3 - Cs*Cv*DS_el*DS_el*(3*DG_el + DS_el))*norm2m*R)/(rom3*Cmixed3);
	const double pdgm1 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m1_el*R)/(rom2*Cmixed2);
	const double pdgm2 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m2_el*R)/(rom2*Cmixed2);
	const double pdget 	=  Cs*DS_el*R/Cmixed2;
	
	const double pdsds 	= (2*Cs*Cs*DG_el*etot_el*R)/Cmixed3 - (DG_el*(Cv*Cv*DG_el*DG_el + Cs*Cv*DG_el*(DG_el + 3*DS_el) + 
			Cs*Cs*(DG_el*DG_el + 3*DG_el*DS_el + 3*DS_el*DS_el))*norm2m*R)/(rom3*Cmixed3);
	const double pdsm1 	= (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*m1_el*R)/(rom2*Cmixed2);
	const double pdsm2 	= (DG_el*((Cs + Cv)*DG_el + 2*Cs*DS_el)*m2_el*R)/(rom2*Cmixed2);
	const double pdset 	= -Cs*DG_el*R/Cmixed2;
	
	const double pm1m1 	= -DG_el*R/(rom*Cmixed);
	
	const double pm2m2 	= -DG_el*R/(rom*Cmixed);
	
	const double Temperature = (2*rom*etot_el - norm2m)/(2*rom*Cmixed);

	const double gas_concentration = 1 - DS_el/ros;
	const double gas_density = DG_el/gas_concentration;

	// printf("eps_g = %.3e\n", gas_concentration);

    const double SpeedSound  = sqrt(gamma*R*Temperature)*sqrt(gas_density/(gas_concentration * rom));  // verificare sound_speed

	// printf("gamma = %.3e - SpeedSound = %.3e\n", gamma, SpeedSound); 

    for (i = 0; i < size3; i++)     A(i) = 0.0;
	for (i = 0; i < size4; i++)     dAdU(i) = 0.0;
			
	// Build A

	A[conta(0,0,0,0,5,2,5,1)] =  DS_el*m1_el/rom2; 
	A[conta(0,0,1,0,5,2,5,1)] = -DG_el*m1_el/rom2;
	A[conta(0,0,2,0,5,2,5,1)] =  DG_el/rom;
	
	A[conta(0,1,0,0,5,2,5,1)] =  DS_el*m2_el/rom2; 
	A[conta(0,1,1,0,5,2,5,1)] = -DG_el*m2_el/rom2;
	A[conta(0,1,3,0,5,2,5,1)] =  DG_el/rom;
	
	A[conta(1,0,0,0,5,2,5,1)] = -DS_el*m1_el/rom2; 
	A[conta(1,0,1,0,5,2,5,1)] =  DG_el*m1_el/rom2;
	A[conta(1,0,2,0,5,2,5,1)] =  DS_el/rom;
	
	A[conta(1,1,0,0,5,2,5,1)] = -DS_el*m2_el/rom2; 
	A[conta(1,1,1,0,5,2,5,1)] =  DG_el*m2_el/rom2;
	A[conta(1,1,3,0,5,2,5,1)] =  DS_el/rom;
	
	A[conta(2,0,0,0,5,2,5,1)] =  -m1_el*m1_el/rom2 + pdg;
	A[conta(2,0,1,0,5,2,5,1)] =  -m1_el*m1_el/rom2 + pds;
	A[conta(2,0,2,0,5,2,5,1)] =  2*m1_el/rom + pm1;
	A[conta(2,0,3,0,5,2,5,1)] =  pm2;
	A[conta(2,0,4,0,5,2,5,1)] =  pet;
	
	A[conta(2,1,0,0,5,2,5,1)] =  -m1_el*m2_el/rom2;
	A[conta(2,1,1,0,5,2,5,1)] =  -m1_el*m2_el/rom2;
	A[conta(2,1,2,0,5,2,5,1)] =  m2_el/rom;
	A[conta(2,1,3,0,5,2,5,1)] =  m1_el/rom;
	
	A[conta(3,0,0,0,5,2,5,1)] =  -m1_el*m2_el/rom2;
	A[conta(3,0,1,0,5,2,5,1)] =  -m1_el*m2_el/rom2;
	A[conta(3,0,2,0,5,2,5,1)] =  m2_el/rom;
	A[conta(3,0,3,0,5,2,5,1)] =  m1_el/rom;
		
	A[conta(3,1,0,0,5,2,5,1)] =  -m2_el*m2_el/rom2 + pdg;
	A[conta(3,1,1,0,5,2,5,1)] =  -m2_el*m2_el/rom2 + pds;
	A[conta(3,1,2,0,5,2,5,1)] =  pm1;
	A[conta(3,1,3,0,5,2,5,1)] =  2*m2_el/rom + pm2;
	A[conta(3,1,4,0,5,2,5,1)] =  pet;
	
	A[conta(4,0,0,0,5,2,5,1)] =  -(m1_el*(etot_el - rom*pdg + p_el))/rom2;
	A[conta(4,0,1,0,5,2,5,1)] =  -(m1_el*(etot_el - rom*pds + p_el))/rom2;
	A[conta(4,0,2,0,5,2,5,1)] =  (etot_el + m1_el*pm1 + p_el)/rom;
	A[conta(4,0,3,0,5,2,5,1)] =  m1_el*pm2/rom;
	A[conta(4,0,4,0,5,2,5,1)] =	m1_el*(1.0 + pet)/rom;
	
	A[conta(4,1,0,0,5,2,5,1)] =  -(m2_el*(etot_el - rom*pdg + p_el))/rom2;
	A[conta(4,1,1,0,5,2,5,1)] =  -(m2_el*(etot_el - rom*pds + p_el))/rom2;
	A[conta(4,1,2,0,5,2,5,1)] =	m2_el*pm1/rom;  
	A[conta(4,1,3,0,5,2,5,1)] =  (etot_el + m2_el*pm2 + p_el)/rom;
	A[conta(4,1,4,0,5,2,5,1)] =  m2_el*(1.0 + pet)/rom;


    //	Build dAdU

	dAdU[conta(0,0,0,0,5,2,5,5)] = -2*DS_el*m1_el/rom3;	 
	dAdU[conta(0,0,0,1,5,2,5,5)] =  (DG_el - DS_el)*m1_el/rom3;
	dAdU[conta(0,0,0,2,5,2,5,5)] = DS_el/rom2;
		
	dAdU[conta(0,0,1,0,5,2,5,5)] = (DG_el - DS_el)*m1_el/rom3; 
	dAdU[conta(0,0,1,1,5,2,5,5)] = 2*DG_el*m1_el/rom3;
	dAdU[conta(0,0,1,2,5,2,5,5)] = -DG_el/rom2;
	 
	dAdU[conta(0,0,2,0,5,2,5,5)] = DS_el/rom2; 
	dAdU[conta(0,0,2,1,5,2,5,5)] = -DG_el/rom2;
			
	
	dAdU[conta(0,1,0,0,5,2,5,5)] = -2*DS_el*m2_el/rom3; 
	dAdU[conta(0,1,0,1,5,2,5,5)] =  (DG_el - DS_el)*m2_el/rom3;
	dAdU[conta(0,1,0,3,5,2,5,5)] = DS_el/rom2;
		
	dAdU[conta(0,1,1,0,5,2,5,5)] =  (DG_el - DS_el)*m2_el/rom3;
	dAdU[conta(0,1,1,1,5,2,5,5)] =  2*DG_el*m2_el/rom3;
	dAdU[conta(0,1,1,3,5,2,5,5)] =  -DG_el/rom2; 
		 
	dAdU[conta(0,1,3,0,5,2,5,5)] = DS_el/rom2; 
	dAdU[conta(0,1,3,1,5,2,5,5)] = -DG_el/rom2;
	
			
	dAdU[conta(1,0,0,0,5,2,5,5)] = 2*DS_el*m1_el/rom3; 
	dAdU[conta(1,0,0,1,5,2,5,5)] = (-DG_el + DS_el)*m1_el/rom3;
	dAdU[conta(1,0,0,2,5,2,5,5)] = -DS_el/rom2;
	
	dAdU[conta(1,0,1,0,5,2,5,5)] = (-DG_el + DS_el)*m1_el/rom3;  
	dAdU[conta(1,0,1,1,5,2,5,5)] = -2*DG_el*m1_el/rom3;
	dAdU[conta(1,0,1,2,5,2,5,5)] = DG_el/rom2;
	 
	dAdU[conta(1,0,2,0,5,2,5,5)] = -DS_el/rom2; 
	dAdU[conta(1,0,2,1,5,2,5,5)] =  DG_el/rom2;
	
	dAdU[conta(1,1,0,0,5,2,5,5)] =  2*DS_el*m2_el/rom3;
	dAdU[conta(1,1,0,1,5,2,5,5)] =  (-DG_el + DS_el)*m2_el/rom3;
	dAdU[conta(1,1,0,3,5,2,5,5)] = -DS_el/rom2; 
	
	dAdU[conta(1,1,1,0,5,2,5,5)] =  (-DG_el + DS_el)*m2_el/rom3;
	dAdU[conta(1,1,1,1,5,2,5,5)] = -2*DG_el*m2_el/rom3;
	dAdU[conta(1,1,1,3,5,2,5,5)] = DG_el/rom2; 
	 
	dAdU[conta(1,1,3,0,5,2,5,5)] = -DS_el/rom2;
	dAdU[conta(1,1,3,1,5,2,5,5)] =  DG_el/rom2;
	
			
	dAdU[conta(2,0,0,0,5,2,5,5)] = -2*m1_el*m1_el/rom3 + pdgdg;	 
	dAdU[conta(2,0,0,1,5,2,5,5)] = -2*m1_el*m1_el/rom3 + pdgds;
	dAdU[conta(2,0,0,2,5,2,5,5)] = -2*m1_el/rom2 + pdgm1;
	dAdU[conta(2,0,0,3,5,2,5,5)] =  pdgm2;
	dAdU[conta(2,0,0,4,5,2,5,5)] =  pdget;
	
	dAdU[conta(2,0,1,0,5,2,5,5)] = 2*m1_el*m1_el/rom3 + pdgds;
	dAdU[conta(2,0,1,1,5,2,5,5)] = 2*m1_el*m1_el/rom3 + pdsds;
	dAdU[conta(2,0,1,2,5,2,5,5)] = -2*m1_el/rom2 + pdsm1;
	dAdU[conta(2,0,1,3,5,2,5,5)] = pdsm2;
	dAdU[conta(2,0,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta(2,0,2,0,5,2,5,5)] = -2*m1_el/rom2 + pdgm1;
	dAdU[conta(2,0,2,1,5,2,5,5)] = -2*m1_el/rom2 + pdsm1;
	dAdU[conta(2,0,2,2,5,2,5,5)] = 2.0/rom + pm1m1;
	 
	dAdU[conta(2,0,3,0,5,2,5,5)] = pdgm2; 
	dAdU[conta(2,0,3,1,5,2,5,5)] = pdsm2;
	dAdU[conta(2,0,3,3,5,2,5,5)] = pm2m2; 
				
	dAdU[conta(2,0,4,0,5,2,5,5)] = pdget;
	dAdU[conta(2,0,4,1,5,2,5,5)] = pdset;
	 
	dAdU[conta(2,1,0,0,5,2,5,5)] = 2*m1_el*m2_el/rom3;	 
	dAdU[conta(2,1,0,1,5,2,5,5)] = 2*m1_el*m2_el/rom3;
	dAdU[conta(2,1,0,2,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(2,1,0,3,5,2,5,5)] = -m1_el/rom2;
	
	dAdU[conta(2,1,1,0,5,2,5,5)] = 2*m1_el*m2_el/rom3; 
	dAdU[conta(2,1,1,1,5,2,5,5)] = 2*m1_el*m2_el/rom3;
	dAdU[conta(2,1,1,2,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(2,1,1,3,5,2,5,5)] = -m1_el/rom2;
	 
	dAdU[conta(2,1,2,0,5,2,5,5)] = -m2_el/rom2; 
	dAdU[conta(2,1,2,1,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(2,1,2,3,5,2,5,5)] = 1.0/rom; 
			
	dAdU[conta(2,1,3,0,5,2,5,5)] = -m1_el/rom2; 
	dAdU[conta(2,1,3,1,5,2,5,5)] = -m1_el/rom2;
	dAdU[conta(2,1,3,2,5,2,5,5)] = 1.0/rom;
	
		
	dAdU[conta(3,0,0,0,5,2,5,5)] = 2*m1_el*m2_el/rom3;	 
	dAdU[conta(3,0,0,1,5,2,5,5)] = 2*m1_el*m2_el/rom3;
	dAdU[conta(3,0,0,2,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(3,0,0,3,5,2,5,5)] = -m1_el/rom2;
		
	dAdU[conta(3,0,1,0,5,2,5,5)] = 2*m1_el*m2_el/rom3; 
	dAdU[conta(3,0,1,1,5,2,5,5)] = 2*m1_el*m2_el/rom3;
	dAdU[conta(3,0,1,2,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(3,0,1,3,5,2,5,5)] = -m1_el/rom2;
	 
	dAdU[conta(3,0,2,0,5,2,5,5)] = -m2_el/rom2; 
	dAdU[conta(3,0,2,1,5,2,5,5)] = -m2_el/rom2;
	dAdU[conta(3,0,2,3,5,2,5,5)] = 1.0/rom; 
			
	dAdU[conta(3,0,3,0,5,2,5,5)] = -m1_el/rom2; 
	dAdU[conta(3,0,3,1,5,2,5,5)] = -m1_el/rom2;
	dAdU[conta(3,0,3,2,5,2,5,5)] = 1.0/rom;
	
	 
	dAdU[conta(3,1,0,0,5,2,5,5)] = 2*m2_el*m2_el/rom3 + pdgdg;	 
	dAdU[conta(3,1,0,1,5,2,5,5)] = 2*m2_el*m2_el/rom3 + pdgds;
	dAdU[conta(3,1,0,2,5,2,5,5)] = pdgm1;
	dAdU[conta(3,1,0,3,5,2,5,5)] = -2*m2_el/rom2 + pdgm2;
	dAdU[conta(3,1,0,4,5,2,5,5)] = pdget;
	
	dAdU[conta(3,1,1,0,5,2,5,5)] = 2*m2_el*m2_el/rom3 + pdgds;
	dAdU[conta(3,1,1,1,5,2,5,5)] = 2*m2_el*m2_el/rom3 + pdsds;
	dAdU[conta(3,1,1,2,5,2,5,5)] = pdsm1;
	dAdU[conta(3,1,1,3,5,2,5,5)] = -2*m2_el/rom2 + pdsm2;
	dAdU[conta(3,1,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta(3,1,2,0,5,2,5,5)] = pdgm1;
	dAdU[conta(3,1,2,1,5,2,5,5)] = pdsm1;
	dAdU[conta(3,1,2,2,5,2,5,5)] = pm1m1;
	
	dAdU[conta(3,1,3,0,5,2,5,5)] = -2*m2_el/rom2 + pdgm2;
	dAdU[conta(3,1,3,1,5,2,5,5)] = -2*m2_el/rom2 + pdsm2;
	dAdU[conta(3,1,3,3,5,2,5,5)] = 2.0/rom + pm2m2; 
	
	dAdU[conta(3,1,4,0,5,2,5,5)] = pdget;
	dAdU[conta(3,1,4,1,5,2,5,5)] = pdset;
			
			
	dAdU[conta(4,0,0,0,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - 2*pdg*rom + pdgdg*rom2)/rom3;
	dAdU[conta(4,0,0,1,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - pdg*rom - pds*rom + pdgds*rom2)/rom3;
	dAdU[conta(4,0,0,2,5,2,5,5)] = (-etot_el + DS_el*(pdg + m1_el*pdgm1) - m1_el*pm1 - p_el + (pdg + m1_el*pdgm1)*DG_el)/rom2;
	dAdU[conta(4,0,0,3,5,2,5,5)] = (m1_el*(DS_el*pdgm2 - pm2 + pdgm2*DG_el))/rom2;
	dAdU[conta(4,0,0,4,5,2,5,5)] = (m1_el*(-1 + DS_el*pdget - pet + pdget*DG_el))/rom2;
	
	dAdU[conta(4,0,1,0,5,2,5,5)] = m1_el*(2*(etot_el + p_el) - pdg*rom - pds*rom + pdgds*rom2)/rom3;
	dAdU[conta(4,0,1,1,5,2,5,5)] = (m1_el*(2*(etot_el + p_el) - 2*pds*rom + pdsds*rom2))/rom3;
	dAdU[conta(4,0,1,2,5,2,5,5)] = (-etot_el + DS_el*(pds + m1_el*pdsm1) - m1_el*pm1 - p_el + (pds + m1_el*pdsm1)*DG_el)/rom2;
	dAdU[conta(4,0,1,3,5,2,5,5)] = m1_el*(pdsm2/rom - pm2/rom2);
	dAdU[conta(4,0,1,4,5,2,5,5)] = (m1_el*(-1 + rom*pdset - pet))/rom2;
	 
	dAdU[conta(4,0,2,0,5,2,5,5)] = (-etot_el + DS_el*(pdg + m1_el*pdgm1) - m1_el*pm1 - p_el + (pdg + m1_el*pdgm1)*DG_el)/rom2;
	dAdU[conta(4,0,2,1,5,2,5,5)] = (-etot_el + DS_el*(pds + m1_el*pdsm1) - m1_el*pm1 - p_el + (pds + m1_el*pdsm1)*DG_el)/rom2;
	dAdU[conta(4,0,2,2,5,2,5,5)] = (m1_el*pm1m1 + 2*pm1)/rom;
	dAdU[conta(4,0,2,3,5,2,5,5)] = pm2/rom;
	dAdU[conta(4,0,2,4,5,2,5,5)] = (1 + pet)/rom;
			
	dAdU[conta(4,0,3,0,5,2,5,5)] = (m1_el*(DS_el*pdgm2 - pm2 + pdgm2*DG_el))/rom2;
	dAdU[conta(4,0,3,1,5,2,5,5)] = m1_el*(pdsm2/rom - pm2/rom2);
	dAdU[conta(4,0,3,2,5,2,5,5)] = pm2/rom;
	dAdU[conta(4,0,3,3,5,2,5,5)] = m1_el*pm2m2/rom;
			
	dAdU[conta(4,0,4,0,5,2,5,5)] = (m1_el*(-1 + DS_el*pdget - pet + pdget*DG_el))/rom2;
	dAdU[conta(4,0,4,1,5,2,5,5)] = (m1_el*(-1 + rom*pdset - pet))/rom2;
	dAdU[conta(4,0,4,2,5,2,5,5)] = (1 + pet)/rom;
	 
	dAdU[conta(4,1,0,0,5,2,5,5)] = m2_el*(2*(etot_el + p_el) - 2*pdg*rom + pdgdg*rom2)/rom3;	 
	dAdU[conta(4,1,0,1,5,2,5,5)] = (m2_el*pdgds)/rom - (m2_el*(pdg + pds))/rom2 + (2*m2_el*(etot_el + p_el))/rom3;
	dAdU[conta(4,1,0,2,5,2,5,5)] = (rom*m2_el*pdgm1 - m2_el*pm1)/rom2; 
	dAdU[conta(4,1,0,3,5,2,5,5)] = (-etot_el - p_el + rom*(pdg + m2_el*pdgm2) - m2_el*pm2)/rom2;
	dAdU[conta(4,1,0,4,5,2,5,5)] = (m2_el*(-1.0 + rom*pdget - pet))/rom2;
	
	dAdU[conta(4,1,1,0,5,2,5,5)] = (m2_el*pdgds)/rom - (m2_el*(pdg + pds))/rom2 + (2*m2_el*(etot_el + p_el))/rom3;
	dAdU[conta(4,1,1,1,5,2,5,5)] = m2_el*(-2*pds/rom2 + pdsds/rom + 2*(etot_el + p_el)/rom3);
	dAdU[conta(4,1,1,2,5,2,5,5)] = m2_el*(rom*pdsm1 - pm1)/rom2;
	dAdU[conta(4,1,1,3,5,2,5,5)] = (-etot_el + rom*(pds + m2_el*pdsm2) - m2_el*pm2 - p_el)/rom2;
	dAdU[conta(4,1,1,4,5,2,5,5)] = (m2_el*(-1 + rom*pdset - pet))/rom2;
	 
	dAdU[conta(4,1,2,0,5,2,5,5)] = (rom*m2_el*pdgm1 - m2_el*pm1)/rom2;
	dAdU[conta(4,1,2,1,5,2,5,5)] = m2_el*(rom*pdsm1 - pm1)/rom2;
	dAdU[conta(4,1,2,2,5,2,5,5)] = m2_el*pm1m1/rom;
	dAdU[conta(4,1,2,3,5,2,5,5)] = pm1/rom;
			
	dAdU[conta(4,1,3,0,5,2,5,5)] = (-etot_el - p_el + rom*(pdg + m2_el*pdgm2) - m2_el*pm2)/rom2;
	dAdU[conta(4,1,3,1,5,2,5,5)] = (-etot_el + rom*(pds + m2_el*pdsm2) - m2_el*pm2 - p_el)/rom2;  
	dAdU[conta(4,1,3,2,5,2,5,5)] = pm1/rom;
	dAdU[conta(4,1,3,3,5,2,5,5)] = (m2_el*pm2m2 + 2*pm2)/rom;
	dAdU[conta(4,1,3,4,5,2,5,5)] = (1 + pet)/rom;
			
	dAdU[conta(4,1,4,0,5,2,5,5)] = (m2_el*(-1.0 + rom*pdget - pet))/rom2;
	dAdU[conta(4,1,4,1,5,2,5,5)] = (m2_el*(-1 + rom*pdset - pet))/rom2;
	dAdU[conta(4,1,4,3,5,2,5,5)] = (1 + pet)/rom;

    for (i = 0; i < nScalarVariables*nScalarVariables; i++)   S(i) = 0.0;

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
		Lstar[pp + 2*nScalarVariables + 1] = N[i]*S[2*nScalarVariables + 1];
		Lstar[pp + 3*nScalarVariables + 0] = N[i]*S[3*nScalarVariables + 0];
		Lstar[pp + 3*nScalarVariables + 1] = N[i]*S[3*nScalarVariables + 1];
		
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

					t = conta(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

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
                    tt = conta(i,j,k,m,nScalarVariables,SpaceDimension,nScalarVariables,nScalarVariables);

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
	invtauStab[2] =	stab_c1*mu_mixture/(rom*h*h) + invtauStab[0];
	invtauStab[3] =	invtauStab[2];
	invtauStab[4] = stab_c1*lambda_mixture/(rom*Cp_mixture*h*h) + invtauStab[0];  

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

				s = conta(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

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
			
	K[conta(0,0,0,0,2,2,5,2)] =  (-2.0 + ctau)*m1_el/rom2;
	K[conta(0,0,0,1,2,2,5,2)] =  ctau*m2_el/rom2;
	K[conta(0,0,1,0,2,2,5,2)] =  (-2.0 + ctau)*m1_el/rom2;
	K[conta(0,0,1,1,2,2,5,2)] =  ctau*m2_el/rom2;
	K[conta(0,0,2,0,2,2,5,2)] =  (2.0 - ctau)/rom;
	K[conta(0,0,3,1,2,2,5,2)] = -ctau/rom;  
	
	K[conta(0,1,0,0,2,2,5,2)] = -m2_el/rom2;
	K[conta(0,1,0,1,2,2,5,2)] = -m1_el/rom2;
	K[conta(0,1,1,0,2,2,5,2)] = -m2_el/rom2;
	K[conta(0,1,1,1,2,2,5,2)] = -m1_el/rom2;  
	K[conta(0,1,2,1,2,2,5,2)] =  1.0/rom; 
	K[conta(0,1,3,0,2,2,5,2)] =  1.0/rom;
	
	
	K[conta(1,0,0,0,2,2,5,2)] = -m2_el/rom2;
	K[conta(1,0,0,1,2,2,5,2)] = -m1_el/rom2;
	K[conta(1,0,1,0,2,2,5,2)] = -m2_el/rom2;
	K[conta(1,0,1,1,2,2,5,2)] = -m1_el/rom2;
	K[conta(1,0,2,1,2,2,5,2)] =  1.0/rom;
	K[conta(1,0,3,0,2,2,5,2)] =  1.0/rom; 
	
	K[conta(1,1,0,0,2,2,5,2)] =  ctau*m1_el/rom2;
	K[conta(1,1,0,1,2,2,5,2)] =  (ctau - 2.0)*m2_el/rom2;
	K[conta(1,1,1,0,2,2,5,2)] =  ctau*m1_el/rom2;
	K[conta(1,1,1,1,2,2,5,2)] =  (ctau - 2.0)*m2_el/rom2;
	K[conta(1,1,2,0,2,2,5,2)] = -ctau/rom;
	K[conta(1,1,3,1,2,2,5,2)] =  (2.0 - ctau)/rom;

	
	KT[0] =  pdg/(DG_el*R) - p_el/(DG_el*DG_el*R);
	KT[1] =  pds/(DG_el*R);
	KT[2] =  pm1/(DG_el*R);
	KT[3] =  pm2/(DG_el*R);
	KT[4] =  pet/(DG_el*R);

    for (i = 0; i < SpaceDimension; i++){
        
        q(i) = 0.0;

        for (j = 0; j < SpaceDimension; j++){

            tau(i*SpaceDimension + j) = 0.0;
            
            for (l = 0; l < nScalarVariables; l++){
                for (m = 0; m < SpaceDimension; m++){

                    tau(i*SpaceDimension + j) += mu_mixture*K[conta(i,j,l,m,2,2,5,2)]*gradU[l*SpaceDimension + m];

                }
            }
        }
        for (l = 0; l < nScalarVariables; l++){
            q(i) -= lambda_mixture*KT(l)*gradU(l*SpaceDimension + i);
        }
    }

    ShockCapturing(mu_mixture, lambda_mixture, Cv_mixture, h, tau, q, rom, gradU, Residual);

    // Build diffusive term: Diffusion tensor

	for ( i = 0; i < nScalarVariables*SpaceDimension; i++ )    G[i] = 0.0;
		
	for (i = 0; i < SpaceDimension; i++){
		for (j = 0; j < SpaceDimension; j++)
			G[(i + 2)*SpaceDimension + j] = -tau[i*SpaceDimension + j];
    }

	for (j = 0; j < SpaceDimension; j++){

		G[4*SpaceDimension + j] = q[j];

		for (i = 0; i < SpaceDimension; i++)
            G[4*SpaceDimension + j] += (-U_gauss[i + 2]/rom*tau[i*SpaceDimension + j]);

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
			FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k];
		}
	}

    // Force contribuution at the Gauss Point   

    for (i = 0; i < nNodalVariables; i++){

		rhs[i] = - sw_conv*FConv[i] - sw_diff*FDiff[i] - sw_stab*FStab[i];

	}
}


}
