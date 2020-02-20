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

#include "custom_elements/compressible_navier_stokes_explicit.h"

namespace Kratos {

unsigned int	cont(
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
               array_1d<double,8>& gradU, // 4*2
               array_1d<double,4>& R)
{
    const int SpaceDimension = 2;
    
   const double alpha = 0.8;                               // Algorithm constant
   const double tol = 1e-3;                               

    unsigned int i;

    double v_sc = 0.0;
    double k_sc = 0.0;
                                          //Shock capturing viscosity
    array_1d<double,SpaceDimension> res_m;
    double res_e = R(3);
    double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;
    

    res_m(0) = R(1); 
    res_m(1) = R(2);

    norm_res_m = sqrt(res_m(0) * res_m(0) + res_m(1) * res_m(1));
    
    for (i = 2; i < 6; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    norm_res_e = sqrt(res_e*res_e);

    double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 6; i < 8; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
 
    if (norm_grade>tol)        k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);

    for (i = 0; i < 4; i++)     tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 2; i++)     q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);

//    printf("%.3e %.3e %.3e %.3e %.3e \n",res_m(0),res_m(1),res_e,v_sc,k_sc);
    
}

void LocalMassMatrix(array_1d<double,12>& LumpedMassMatrix,const array_1d<double,3>& N, const unsigned int nodesElement)
{

    unsigned int t;

    for (unsigned int i = 0; i < nodesElement; i++){

        t = 3*i;

        LumpedMassMatrix(t)     = 0.0;
        LumpedMassMatrix(t + 1) = 0.0;
        LumpedMassMatrix(t + 2) = 0.0;
        LumpedMassMatrix(t + 3) = 0.0;

        for (unsigned int j = 0; j < nodesElement; j++){
            LumpedMassMatrix(t)     += N(i)*N(j);
            LumpedMassMatrix(t + 1) += N(i)*N(j);
            LumpedMassMatrix(t + 2) += N(i)*N(j);
            LumpedMassMatrix(t + 3) += N(i)*N(j);
        }

    }
}


template<>
void CompressibleNavierStokesExplicit<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int SpaceDimension = 2;
    unsigned int nScalarVariables = SpaceDimension + 2;
    unsigned int nodesElement = 3;
    unsigned int nNodalVariables  = nodesElement*nScalarVariables;

    if (rResult.size() != nNodalVariables)
        rResult.resize(nNodalVariables, false);         // A che serve sta roba?

    for(unsigned int i=0; i<nodesElement; i++)
    {
       rResult[i*(nScalarVariables)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();	// Collect current value of nodal variable into a local vector U?
       rResult[i*(nScalarVariables)+1]  =  this->GetGeometry()[i].GetDof(MOMENTUM_X).EquationId();
       rResult[i*(nScalarVariables)+2]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Y).EquationId();
       rResult[i*(nScalarVariables)+3]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void CompressibleNavierStokesExplicit<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int SpaceDimension = 2;
    unsigned int nScalarVariables = SpaceDimension + 2;
    unsigned int nodesElement = 3;
    unsigned int nNodalVariables  = nodesElement*nScalarVariables;

    if (ElementalDofList.size() != nNodalVariables)     // A che serve sta roba ?
        ElementalDofList.resize(nNodalVariables);

    for(unsigned int i=0; i<nodesElement; i++)
    {
        ElementalDofList[i*(nScalarVariables)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);		// Difference with the operation made before? 
        ElementalDofList[i*(nScalarVariables)+1]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_X);
        ElementalDofList[i*(nScalarVariables)+2]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Y);
        ElementalDofList[i*(nScalarVariables)+3]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}

template<>
void CompressibleNavierStokesExplicit<2>::ComputeGaussPointLHSContribution(BoundedMatrix<double,12,12>& lhs, const ElementDataStruct& data)
{   
    const unsigned int nodesElement = 3;
    const unsigned int SpaceDimension = 2;
    const unsigned int nScalarVariables = SpaceDimension + 2;
    const unsigned int nNodalVariables = nScalarVariables*nodesElement;

    array_1d<double,nNodalVariables>      LumpedMassMatrix;

    const array_1d<double,nodesElement>&   N = data.N;

    LocalMassMatrix(LumpedMassMatrix, N, nodesElement);

    for (unsigned int i = 0; i < nNodalVariables; i++){
        for (unsigned int j = 0; j < nNodalVariables; j++){
            lhs(i,j) = 0.0;
            if (i == j)     lhs(i,j) = LumpedMassMatrix(i);
        }
    }
}



template<>
void CompressibleNavierStokesExplicit<2>::ComputeGaussPointRHSContribution(array_1d<double,12>& rhs, const ElementDataStruct& data)
{
    const int nodesElement = 3;
    const int SpaceDimension = 2;
    const int nScalarVariables = SpaceDimension + 2;
    const int nNodalVariables = nScalarVariables*nodesElement;
    const double h = data.h;

    unsigned int i, j, k, m, s, t, tt, p, pp;
    
    const unsigned int size3 = nScalarVariables * SpaceDimension * nScalarVariables;
    const unsigned int size4 = nScalarVariables * SpaceDimension * nScalarVariables * nScalarVariables;

    array_1d<double, size3> A;
    array_1d<double, size4> dAdU;

    array_1d<double,nNodalVariables>                     LumpedMassMatrix;
    array_1d<double,nScalarVariables>                    U_gauss;
    array_1d<double,nNodalVariables>                     U;
    array_1d<double,nNodalVariables>                     Un;
    array_1d<double,nNodalVariables>                     up;
    array_1d<double,nScalarVariables>                    L;
    array_1d<double,nScalarVariables>                    R;
    array_1d<double,nNodalVariables*nScalarVariables>    Lstar;
    array_1d<double,nScalarVariables*SpaceDimension>     G;
    array_1d<double,nScalarVariables*nScalarVariables>   S;
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
    
	
    // Replace with right coefficient according to the time scheme used
    const double& ctau = 0.5;
    
    const double& dt = data.dt;

    // const double& bdf0 = data.bdf0;
    // const double& bdf1 = data.bdf1;
    // const double& bdf2 = data.bdf2;

    const BoundedMatrix<double,nodesElement,nScalarVariables>& UU = data.U;			// Lo mismo de antes // Da tenere come sono e riordinare
    const BoundedMatrix<double,nodesElement,nScalarVariables>& UUn = data.Un;
    const BoundedMatrix<double,nodesElement,nScalarVariables>& Up = data.Up;

    const BoundedMatrix<double,nodesElement,SpaceDimension>& f_ext = data.f_ext;			
    const array_1d<double,nodesElement>& r = data.r;
    const double mu = data.mu;
    // const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.c_v;
    const double gamma = data.gamma;
    const double cp = cv*gamma;
  
    const double sw_conv = 1.0;
    const double sw_diff = 1.0;
    const double sw_stab = 1.0;


    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;    

    // Get shape function values
    const array_1d<double,nodesElement>& N = data.N;					        // N (in which point?)
    const BoundedMatrix<double,nodesElement,SpaceDimension>& DN = data.DN_DX;	// DN (in which point ?)

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,SpaceDimension> f_gauss = prod(trans(f_ext), N);      // Da tenere così
    const double r_gauss = N(0)*r(0) + N(1)*r(1) + N(2)*r(2);
    //    const BoundedMatrix<double,dim,BlockSize> grad_U = prod(trans(DN), U); /*TO DO: WRONG MIGHT BE REMOVED*/   // Appost
    //    const array_1d<double,nScalarVariables> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N); 	// OK


    // Define U and Udot
    for (i = 0; i < nodesElement; i++){
        for (j = 0; j < nScalarVariables; j++){
            
            p = i*nScalarVariables + j; 

            U(p) = UU(i,j);
            Un(p) = UUn(i,j);

            up(p) = (U(p) - Un(p))/dt;
            //up(p) = Up(i,j);
            
        }
    }

    for (i = 0; i < nScalarVariables; i++){
        U_gauss(i) = 0;
        for (j = 0; j < nodesElement; j++){
            U_gauss(i) += N(j)*U(i + j*nScalarVariables);
        }
    }

    for (i = 0; i < nScalarVariables*SpaceDimension; i++)       gradU(i) = 0.0;
// Compute the gradient of the variables /dro/dx,dro/dy, dm1/dx dm1/dy dm2/dx dm2/dy de/dx de/dy)
    for (i = 0; i < nScalarVariables; i++){         // Controllare
        for (k = 0; k < nodesElement; k++){

            gradU(i*SpaceDimension + 0) += DN(k,0)*U(i + k*nScalarVariables);
            gradU(i*SpaceDimension + 1) += DN(k,1)*U(i + k*nScalarVariables);

        }
    }

    for (i = 0; i < nNodalVariables*nScalarVariables; i++ )  NN(i) = 0.0;

    
    for (j = 0; j < nodesElement; j++){

        NN(0 * nNodalVariables + j*nScalarVariables    ) = N(j);
        NN(1 * nNodalVariables + j*nScalarVariables + 1) = N(j);
        NN(2 * nNodalVariables + j*nScalarVariables + 2) = N(j);
        NN(3 * nNodalVariables + j*nScalarVariables + 3) = N(j);

    }
    
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

	}


    const double ro_el = U_gauss(0);
    const double m1_el = U_gauss(1);
    const double m2_el = U_gauss(2);
    const double etot_el = U_gauss(3);

    const double u1_el = m1_el/ro_el;
	const double u2_el = m2_el/ro_el;

	const double norm2u = u1_el*u1_el + u2_el*u2_el;
	const double norm_u = sqrt(norm2u);
	const double norm2m = m1_el*m1_el + m2_el*m2_el;

	const double p_el = (gamma - 1.0)*(etot_el - 0.5 * norm2m/ro_el);
    

    const double SpeedSound = sqrt( gamma * (gamma - 1.0) * (etot_el / ro_el  - 0.5 * norm2u) );


    for (i = 0; i < size3; i++)     A(i) = 0.0;
	
	for (i = 0; i < size4; i++)      dAdU(i) = 0.0;
			
	const double    dpdro = 0.5*(gamma - 1)*norm2u;
	const double    dpdm1 = -(gamma - 1)*u1_el;
	const double    dpdm2 = -(gamma - 1)*u2_el;
	const double    dpde  =  (gamma - 1);

	const double    d2pdro2 	= -(gamma - 1)*norm2u/ro_el;
	const double    d2pdrodm1 	=  (gamma - 1)*u1_el/ro_el;
	const double    d2pdrodm2 	=  (gamma - 1)*u2_el/ro_el;

	const double    d2pdm12 	= -(gamma - 1)/ro_el;
	const double    d2pdm22 	= -(gamma - 1)/ro_el;

	// Build A

	A[cont(0,0,1,0,4,2,4,1)] = 1.0;
	A[cont(0,1,2,0,4,2,4,1)] = 1.0;

	A[cont(1,0,0,0,4,2,4,1)] = dpdro - u1_el*u1_el;
	A[cont(1,0,1,0,4,2,4,1)] = dpdm1 + 2*u1_el;
	A[cont(1,0,2,0,4,2,4,1)] = dpdm2;
	A[cont(1,0,3,0,4,2,4,1)] = dpde;

	A[cont(1,1,0,0,4,2,4,1)] = -u1_el*u2_el;
	A[cont(1,1,1,0,4,2,4,1)] = u2_el;
	A[cont(1,1,2,0,4,2,4,1)] = u1_el;

	A[cont(2,0,0,0,4,2,4,1)] = -u1_el*u2_el;
	A[cont(2,0,1,0,4,2,4,1)] = u2_el;
	A[cont(2,0,2,0,4,2,4,1)] = u1_el;

	A[cont(2,1,0,0,4,2,4,1)] = dpdro - u2_el*u2_el;
	A[cont(2,1,1,0,4,2,4,1)] = dpdm1;
	A[cont(2,1,2,0,4,2,4,1)] = dpdm2 + 2*u2_el;
	A[cont(2,1,3,0,4,2,4,1)] = dpde;

	A[cont(3,0,0,0,4,2,4,1)] = dpdro*u1_el - (etot_el + p_el)*m1_el/(ro_el*ro_el);
	A[cont(3,0,1,0,4,2,4,1)] = (etot_el + p_el)/ro_el + u1_el * dpdm1;
	A[cont(3,0,2,0,4,2,4,1)] = u1_el * dpdm2;
	A[cont(3,0,3,0,4,2,4,1)] = u1_el * (1 + dpde);

	A[cont(3,1,0,0,4,2,4,1)] = dpdro*u2_el - (etot_el + p_el)*m2_el/(ro_el*ro_el);
	A[cont(3,1,1,0,4,2,4,1)] = u2_el * dpdm1;
	A[cont(3,1,2,0,4,2,4,1)] = (etot_el + p_el)/ro_el + u2_el * dpdm2;
	A[cont(3,1,3,0,4,2,4,1)] = u2_el * (1 + dpde);

    //	Build dAdU

	dAdU[cont(1,0,0,0,4,2,4,4)] = d2pdro2 + 2*u1_el*u1_el/ro_el;
	dAdU[cont(1,0,0,1,4,2,4,4)] = d2pdrodm1 - 2*u1_el/ro_el;
	dAdU[cont(1,0,0,2,4,2,4,4)] = d2pdrodm2;

	dAdU[cont(1,0,1,0,4,2,4,4)] = d2pdrodm1 - 2.0 * u1_el/ro_el;
	dAdU[cont(1,0,1,1,4,2,4,4)] = d2pdm12 + 2.0/ro_el;

	dAdU[cont(1,0,2,0,4,2,4,4)] = d2pdrodm2;
	dAdU[cont(1,0,2,2,4,2,4,4)] = d2pdm22;

	/// ------

	dAdU[cont(1,1,0,0,4,2,4,4)] = 2*u1_el*u2_el/ro_el;
	dAdU[cont(1,1,0,1,4,2,4,4)] = -u2_el/ro_el;
	dAdU[cont(1,1,0,2,4,2,4,4)] = -u1_el/ro_el;

	dAdU[cont(1,1,1,0,4,2,4,4)] = -u2_el/ro_el;
	dAdU[cont(1,1,1,2,4,2,4,4)] = 1.0/ro_el;

	dAdU[cont(1,1,2,0,4,2,4,4)] = -u1_el/ro_el;
	dAdU[cont(1,1,2,1,4,2,4,4)] = 1.0/ro_el;

	/// ------

	dAdU[cont(2,0,0,0,4,2,4,4)] = 2*u1_el*u2_el/ro_el;
	dAdU[cont(2,0,0,1,4,2,4,4)] = -u2_el/ro_el;
	dAdU[cont(2,0,0,2,4,2,4,4)] = -u1_el/ro_el;

	dAdU[cont(2,0,1,0,4,2,4,4)] = -u2_el/ro_el;
	dAdU[cont(2,0,1,2,4,2,4,4)] = 1.0/ro_el;

	dAdU[cont(2,0,2,0,4,2,4,4)] = -u1_el/ro_el;
	dAdU[cont(2,0,2,1,4,2,4,4)] = 1.0/ro_el;



	dAdU[cont(2,1,0,0,4,2,4,4)] = d2pdro2 + 2*u2_el*u2_el/ro_el;
	dAdU[cont(2,1,0,1,4,2,4,4)] = d2pdrodm1;
	dAdU[cont(2,1,0,2,4,2,4,4)] = d2pdrodm2 - 2*u2_el/ro_el;

	dAdU[cont(2,1,1,0,4,2,4,4)] = d2pdrodm1;
	dAdU[cont(2,1,1,1,4,2,4,4)] = d2pdm12;

	dAdU[cont(2,1,2,0,4,2,4,4)] = d2pdrodm2 - 2*u2_el/ro_el;
	dAdU[cont(2,1,2,2,4,2,4,4)] = d2pdm22 + 2.0/ro_el;

	/// ------

	dAdU[cont(3,0,0,0,4,2,4,4)] = (m1_el * (2*etot_el + 2*p_el + ro_el*(-2*dpdro + d2pdro2*ro_el)))/(ro_el*ro_el*ro_el);
	dAdU[cont(3,0,0,1,4,2,4,4)] = -(etot_el + dpdm1*m1_el + p_el - dpdro*ro_el - d2pdrodm1*m1_el*ro_el)/(ro_el*ro_el);
	dAdU[cont(3,0,0,2,4,2,4,4)] = (m1_el*(-dpdm2 + d2pdrodm2*ro_el))/(ro_el*ro_el);
	dAdU[cont(3,0,0,3,4,2,4,4)] = -(1 + dpde) * m1_el/(ro_el*ro_el);

	dAdU[cont(3,0,1,0,4,2,4,4)] = -(etot_el + dpdm1*m1_el + p_el - dpdro*ro_el - d2pdrodm1*m1_el*ro_el)/(ro_el*ro_el);
	dAdU[cont(3,0,1,1,4,2,4,4)] = (2*dpdm1 + d2pdm12*m1_el)/ro_el;
	dAdU[cont(3,0,1,2,4,2,4,4)] = dpdm2/ro_el;
	dAdU[cont(3,0,1,3,4,2,4,4)] = (1 + dpde)/ro_el;

	dAdU[cont(3,0,2,0,4,2,4,4)] = (m1_el*(-dpdm2 + d2pdrodm2*ro_el))/(ro_el * ro_el);
	dAdU[cont(3,0,2,1,4,2,4,4)] = dpdm2/ro_el;
	dAdU[cont(3,0,2,2,4,2,4,4)] = (d2pdm22*m1_el)/ro_el;

	dAdU[cont(3,0,3,0,4,2,4,4)] = -(1 + dpde)*m1_el/(ro_el*ro_el);
	dAdU[cont(3,0,3,1,4,2,4,4)] = (1 + dpde)/ro_el;

	dAdU[cont(3,1,0,0,4,2,4,4)] = (m2_el*(2*etot_el + 2*p_el - 2*dpdro*ro_el + d2pdro2*ro_el*ro_el))/(ro_el*ro_el*ro_el);
	dAdU[cont(3,1,0,1,4,2,4,4)] = (m2_el*(-dpdm1 + d2pdrodm1*ro_el))/(ro_el*ro_el);
	dAdU[cont(3,1,0,2,4,2,4,4)] = -(etot_el + p_el - dpdro*ro_el + m2_el*(dpdm2 - d2pdrodm2*ro_el))/(ro_el*ro_el);
	dAdU[cont(3,1,0,3,4,2,4,4)] = (-1 - dpde)*m2_el/(ro_el*ro_el);

	dAdU[cont(3,1,1,0,4,2,4,4)] = (m2_el*(-dpdm1 + d2pdrodm1*ro_el))/(ro_el*ro_el);
	dAdU[cont(3,1,1,1,4,2,4,4)] = (d2pdm12*m2_el)/ro_el;
	dAdU[cont(3,1,1,2,4,2,4,4)] = dpdm1/ro_el;

	dAdU[cont(3,1,2,0,4,2,4,4)] = -(etot_el + p_el - dpdro*ro_el + m2_el*(dpdm2 - d2pdrodm2*ro_el))/(ro_el*ro_el);
	dAdU[cont(3,1,2,1,4,2,4,4)] = dpdm1/ro_el;
	dAdU[cont(3,1,2,2,4,2,4,4)] = (2*dpdm2 + d2pdm22*m2_el)/ro_el;
	dAdU[cont(3,1,2,3,4,2,4,4)] = (1 + dpde)/ro_el;

	dAdU[cont(3,1,3,0,4,2,4,4)] = (-1.0 - dpde)*m2_el/(ro_el*ro_el);
	dAdU[cont(3,1,3,2,4,2,4,4)] = (1.0 + dpde)/ro_el;

    for (i = 0; i < nScalarVariables*nScalarVariables; i++)   S(i) = 0.0;

    S(1*nScalarVariables + 0) = f_gauss(0);
    S(2*nScalarVariables + 0) = f_gauss(1);
    S(3*nScalarVariables + 0) = r_gauss;
    S(3*nScalarVariables + 1) = f_gauss(0);
    S(3*nScalarVariables + 2) = f_gauss(1);



    for (k = 0; k < nScalarVariables; k++){
		for ( s = 0; s < nNodalVariables; s++){
			
			p = s*nScalarVariables + k;
			
			Lstar[p] = 0.0;

			for (i = 0; i < nScalarVariables; i++){
				
				pp = i*nNodalVariables + s;

				Lstar[p] += (-NN[pp]*S[i*nScalarVariables + k]);

				for (j = 0; j < SpaceDimension; j++){

					t = cont(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

					Lstar[p] += (-A[t]*gradV[(SpaceDimension*i+j)*nNodalVariables + s]);
                    
					for (m = 0; m < nScalarVariables; m++){

						tt = cont(i,j,k,m,nScalarVariables,SpaceDimension,nScalarVariables,nScalarVariables);

						Lstar[p] += (-NN[pp]*dAdU[tt]*gradU[m*SpaceDimension + j]);

					}
				}

			}
		}
	}

	invtauStab[0] =	stab_c2*(norm_u + SpeedSound)/h; 
	invtauStab[1] =	stab_c1*mu/(ro_el*h*h) + invtauStab[0]; 
	invtauStab[2] =	invtauStab[1];
	invtauStab[3] = stab_c1*lambda/(ro_el*cp*h*h) + invtauStab[0];  

    L[0] = 0.0;
    L[1] = -S[1*nScalarVariables + 0]*U_gauss[0];
    L[2] = -S[2*nScalarVariables + 0];
    L[3] = 0.0;
    
    for (k = 0; k < SpaceDimension - 1 ; k++){
        L[3] -= S[3*nScalarVariables + k]*U_gauss[k];
    }

    for (i = 0; i < nScalarVariables; i++ ){

		R[i] = 0.0;

		for (k = 0; k < nScalarVariables; k++){

//			L[i] -= S[i*nScalarVariables + k]*U_gauss[k];

			for (j = 0; j < SpaceDimension; j++){

				s = cont(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

				L[i] += A[s]*gradU[k*SpaceDimension + j];
			}
		}

		R[i] = -L[i];

		// for (k = 0; k < nNodalVariables; k++){

		// 	R[i] -= NN[i*nNodalVariables + k]*up[k];

		// }
        
        for (k = 0; k < nodesElement; k++){
            R[i] -= N[k]*up[i + nScalarVariables*k];
        }


	}

	// for (s = 0; s < nNodalVariables; s++){

	// 	FConv[s] = 0.0;

	// 	for (k = 0; k < nodesElement; k++){
	// 		FConv[s] += NN[i*nNodalVariables + s]*L[i];
	// 	}
	// }

    for (i = 0; i < nScalarVariables; i++){
        for (k = 0; k < nodesElement; k++){
			FConv[i + k*nScalarVariables] = N[k]*L[i];
		}
    }


    // Build diffusive term: stress tensor and thermal diffusion

    double  divrom, divm;

    divm   = gradU[1*SpaceDimension + 0] + gradU[2*SpaceDimension + 1];
	divrom = U_gauss[1]*gradU[0*SpaceDimension + 0] + U_gauss[2]*gradU[0*SpaceDimension + 1];

	for (i = 0; i < SpaceDimension; i++){

		q[i] = lambda*etot_el/(cv*ro_el*ro_el)*gradU[0*SpaceDimension + i] - lambda/(ro_el*cv)*gradU[3*SpaceDimension + i];

		for (j = 0; j < SpaceDimension; j++){

			tau[i*SpaceDimension + j] = mu/ro_el * (gradU[(i+1)*SpaceDimension + j] + gradU[(j+1)*SpaceDimension + i])
					  - mu/(ro_el * ro_el) * (U_gauss[i+1] * gradU[0*SpaceDimension + j] + U_gauss[j+1] * gradU[0*SpaceDimension + i]);

			q[i] 	 += (-lambda * U_gauss[j+1] * U_gauss[j+1]/(ro_el*ro_el*ro_el*cv) * gradU[0*SpaceDimension + i] +
						  lambda * U_gauss[j+1] / (ro_el * ro_el * cv) * gradU[(j+1)*SpaceDimension + i]);

		}

		tau[i*SpaceDimension + i] += (-ctau*mu/ro_el*divm + ctau*mu/(ro_el*ro_el)*divrom);	// Controllare coefficienti di tau in 2d

	}



    ShockCapturing(mu, lambda, cv, h, tau, q, ro_el, gradU, R);



    // Build diffusive term: Diffusion tensor

	for ( i = 0; i < nScalarVariables*SpaceDimension; i++ )    G[i] = 0.0;
		
	for (i = 0; i < SpaceDimension; i++){
		for (j = 0; j < SpaceDimension; j++)
			G[(i+1)*SpaceDimension + j] = -tau[i*SpaceDimension + j];
    }

	for (j = 0; j < SpaceDimension; j++){

		G[3*SpaceDimension + j] = q[j];

		for (i = 0; i < SpaceDimension; i++)
            G[3*SpaceDimension + j] += (-U_gauss[i+1]/ro_el*tau[i*SpaceDimension + j]);

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
			FStab[s] += Lstar[s*nScalarVariables + k]*R[k]/invtauStab[k];
		}
	}

    // Force contribuution at the Gauss Point   

    for (i = 0; i < nNodalVariables; i++){

		rhs[i] = - sw_conv*FConv[i] - sw_diff*FDiff[i] - sw_stab*FStab[i];

	}


}


}
