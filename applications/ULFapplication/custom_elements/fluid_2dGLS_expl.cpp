//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Pavel Ryzhakov and Julio Marti 
//



//#define GRADPN_FORM

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/fluid_2dGLS_expl.h"
#include "utilities/math_utils.h"
#include "ULF_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
  
  //THIS IS A COMPRESSIBLE FLUID ELEMENT, WITH GLS STABILIZATION, RUNGE-KUTTA Momentum Time integration, FRACTIONAL STEP
  //************************************************************************************
  //************************************************************************************
  Fluid2DGLS_expl::Fluid2DGLS_expl(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!
  }
  //************************************************************************************
  //************************************************************************************
  Fluid2DGLS_expl::Fluid2DGLS_expl(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {


  }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::CalculateLumpedMass()
  {
    //note that for the compressible case, rho will be also a variable

    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

    //double Area;
    //GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
    double Area = GeometryUtils::CalculateVolume2D(GetGeometry());
    double lumped_mass_fac = Area * 0.33333333333333333;
    //filling in the diagonal of the lumped mass matrix,  (later I can change it to vector...)
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho0;
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho1;
    GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho2;
  }
  //************************************************************************************
  //************************************************************************************
  Element::Pointer Fluid2DGLS_expl::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    KRATOS_TRY

      return Element::Pointer(new Fluid2DGLS_expl(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
  }

  Fluid2DGLS_expl::~Fluid2DGLS_expl()
  {
  }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
  }

  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  {
    //KRATOS_WATCH("Empty function for this element")
    //KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
  }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::CalculateGalerkinMomentumResidual(VectorType& GalerkinRHS)
  {
    KRATOS_TRY
      ///////////////////////NECESSARY LOCALS///////////////////////////////////////////
      boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
    boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
    array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
    boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
    boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
    boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
    array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
    array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
    ///////////////////////////////////////////////////////////////////////////////////


    //first we compute  the force term and pressure gradient terms:
    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);
    //if (Area<0.0000000001)  KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    //getting the velocity on the nodes and other necessary variabless
    const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);//
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);//
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);


    //====================================================================
    //calculating viscosity and density
    double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
    double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

    //VISCOUS CONTRIBUTION
    // += Laplacian * nu; --> ONE GAUSS POINT
    //msWorkMatrix is used now to store the element laplacian 3x3
 
    /*   
	 noalias(msWorkMatrix) = Area*density*nu * prod(msDN_DX,trans(msDN_DX));

	 //x comp
	 GalerkinRHS[0]=-1.0*(msWorkMatrix(0,0)*vel0[0]+msWorkMatrix(0,1)*vel1[0]+msWorkMatrix(0,2)*vel2[0]);
	 //y comp
	 GalerkinRHS[1]=-1.0*(msWorkMatrix(0,0)*vel0[1]+msWorkMatrix(0,1)*vel1[1]+msWorkMatrix(0,2)*vel2[1]);

	 //x comp
	 GalerkinRHS[2]=-1.0*(msWorkMatrix(1,0)*vel0[0]+msWorkMatrix(1,1)*vel1[0]+msWorkMatrix(1,2)*vel2[0]);
	 //y comp
	 GalerkinRHS[3]=-1.0*(msWorkMatrix(1,0)*vel0[1]+msWorkMatrix(1,1)*vel1[1]+msWorkMatrix(1,2)*vel2[1]);

	 //x comp
	 GalerkinRHS[4]=-1.0*(msWorkMatrix(2,0)*vel0[0]+msWorkMatrix(2,1)*vel1[0]+msWorkMatrix(2,2)*vel2[0]);
	 //y comp
	 GalerkinRHS[5]=-1.0*(msWorkMatrix(2,0)*vel0[1]+msWorkMatrix(2,1)*vel1[1]+msWorkMatrix(2,2)*vel2[1]);
    
    */
    /// ANOTHER FORM OF VISCOUS MATRIX CONTRIBUTION
    ////////////////////////////////////////

    boost::numeric::ublas::bounded_matrix<double,3,6> msB;
    boost::numeric::ublas::bounded_matrix<double,3,3> ms_constitutive_matrix;
    boost::numeric::ublas::bounded_matrix<double,3,6> ms_temp;
    boost::numeric::ublas::bounded_matrix<double,6,6> DampingMatrix;



    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    
    //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
    //filling matrix B
    for (unsigned int i=0; i<NumberOfNodes; i++)
      {
        unsigned int index = dim*i;
        msB(0,index+0)=msDN_DX(i,0);
        msB(0,index+1)= 0.0;
        msB(1,index+0)=0.0;
        msB(1,index+1)= msDN_DX(i,1);
        msB(2,index+0)= msDN_DX(i,1);
        msB(2,index+1)= msDN_DX(i,0);
      }

    //constitutive tensor
    ms_constitutive_matrix(0,0) = (4.0/3.0)*nu*density;
    ms_constitutive_matrix(0,1) = -2.0/3.0*nu*density;
    ms_constitutive_matrix(0,2) = 0.0;
    ms_constitutive_matrix(1,0) = -2.0/3.0*nu*density;
    ms_constitutive_matrix(1,1) = 4.0/3.0*nu*density;
    ms_constitutive_matrix(1,2) = 0.0;
    ms_constitutive_matrix(2,0) = 0.0;
    ms_constitutive_matrix(2,1) = 0.0;
    ms_constitutive_matrix(2,2) = nu*density;

    //calculating viscous contributions
    ms_temp = prod( ms_constitutive_matrix , msB);
    noalias(DampingMatrix) = prod( trans(msB) , ms_temp);

    DampingMatrix *= Area;

    msAuxVec[0]=vel0[0];
    msAuxVec[1]=vel0[1];
    msAuxVec[2]=vel1[0];
    msAuxVec[3]=vel1[1];
    msAuxVec[4]=vel2[0];
    msAuxVec[5]=vel2[1];

    noalias(GalerkinRHS)=-prod(DampingMatrix, msAuxVec);

   	

    //and now we add the pressure gradient and the force term
    //external forces (component)

    const array_1d<double,3> body_force = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
							     GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
							     GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
    unsigned int number_of_nodes=3;
    for(unsigned int i = 0; i<number_of_nodes; i++)
      {
        //f=A*N_I*b, N_I=0.33333333 for 1 Gauss point
        GalerkinRHS[i*2] += body_force[0]* density * Area * 0.3333333333333;
        GalerkinRHS[i*2+1] += body_force[1] * density * Area * 0.3333333333333;

      }
    	

    //Now we shall add the Gp term(integrated by parts)
    double p_avg = msN[0]* p_n0 + msN[1] * p_n1 + msN[2] * p_n2;
    p_avg *= Area;
    //p_avg *= 0.0;
    GalerkinRHS[0] += msDN_DX(0, 0) * p_avg;
    GalerkinRHS[1] += msDN_DX(0, 1) * p_avg;
    GalerkinRHS[2] += msDN_DX(1, 0) * p_avg;
    GalerkinRHS[3] += msDN_DX(1, 1) * p_avg;
    GalerkinRHS[4] += msDN_DX(2, 0) * p_avg;
    GalerkinRHS[5] += msDN_DX(2, 1) * p_avg;
	

    array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_OLD_OLD);
    rhs0[0] += GalerkinRHS[0];
    rhs0[1] += GalerkinRHS[1];
    rhs0[2] = 0.0;

    array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_OLD_OLD);
    rhs1[0] += GalerkinRHS[2];
    rhs1[1] += GalerkinRHS[3];
    rhs1[2] = 0.0;

    array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_OLD_OLD);
    rhs2[0] += GalerkinRHS[4];
    rhs2[1] += GalerkinRHS[5];
    rhs2[2] = 0.0;


    KRATOS_CATCH("")

      }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::CalculateRHSVector(VectorType& Galerkin_RHS, double& dt)
  {
    KRATOS_TRY

      KRATOS_CATCH("")

      }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      ///////////////////////NECESSARY LOCALS///////////////////////////////////////////
      boost::numeric::ublas::bounded_matrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
    boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
    array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
    boost::numeric::ublas::bounded_matrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
    boost::numeric::ublas::bounded_matrix<double,2,6> msConvOp = ZeroMatrix(2,6);
    boost::numeric::ublas::bounded_matrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
    array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
    array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,3> ms_temp_vec_np = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
    ///////////////////////////////////////////////////////////////////////////////////

    if(rRightHandSideVector.size() != 3)
      {
        rLeftHandSideMatrix.resize(3,3,false);
        rRightHandSideVector.resize(3,false);
      }

    double dt = rCurrentProcessInfo[DELTA_TIME];

    //fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
    //so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2
    const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
    //old iteration can be used if we want to iterate between 1st and 2nd fractional steps

    const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

    double one_sixth = 0.166666666666667; 
    //in msAuxVec we store the velocity, (not the frac step vel, but u_n, the one that enters the stabilization)
    msAuxVec[0]=fv0[0];
    msAuxVec[1]=fv0[1];
    msAuxVec[2]=fv1[0];
    msAuxVec[3]=fv1[1];
    msAuxVec[4]=fv2[0];
    msAuxVec[5]=fv2[1];

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

    //calculating average density and viscosity
    double nu = 0.33333333333333*(nu0 + nu1 + nu2 );
    double density = 0.33333333333333*(rho0 + rho1 + rho2 );

    //ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
    //but with one integration N=0.333333333
    ms_vel_gauss[0] =  0.33333333333333*(fv0[0]+fv1[0]+fv2[0]);
    ms_vel_gauss[1] =  0.33333333333333*(fv0[1]+fv1[1]+fv2[1]);

    //calculating parameter tau (saved internally to each element)
    double h = sqrt(2.00*Area);
    double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
    norm_u = sqrt(norm_u);
    double tau = 1.00 / ( 1.0 * 4.00*nu/(h*h) + 1.0/dt);
    
    //AND NOW WE ADD THE RESPECTIVE CONTRIBUTIONS TO THE RHS AND LHS of THE SECOND FRAC STEP
    //we use Backward Euler for this step, therefore stab. contribution no RHS +=Tau1*(gradQ, residual)
    //								   and LHS +=Tau1*(gradQ, gradP)
    //laplacian term	       L = Dt * gradN * trans(gradN);
    //stabilization term       Spp = tau * gradN * trans(gradN);
    //WATCH OUT for DIVISION with RHO - check if it changes or not in case of Momentum being the primary Variable
    //
    //	msWorkMatrix stores the element laplacian
    //
    noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
    noalias(rLeftHandSideMatrix) = (one_sixth*dt + tau) * Area*msWorkMatrix;

    //////////////////////////////////////////////////////////
    ////////////		AND NOW RHS	//////////////////
    //////////////////////////////////////////////////////////

    //Dirichlet contribution  (that is: LHS*p_new)
    ms_temp_vec_np[0] = p0;
    ms_temp_vec_np[1] = p1;
    ms_temp_vec_np[2] = p2;
    //LHS is already multiplied by AREA
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);

    //NOW RHS-=dt L p_old
    //changing the meaning of temp_vec_np

    ms_temp_vec_np[0] = p0_old;
    ms_temp_vec_np[1] = p1_old;
    ms_temp_vec_np[2] = p2_old;

    noalias(rRightHandSideVector) += one_sixth*Area*dt* (prod(msWorkMatrix,ms_temp_vec_np)) ;

    //here we have the Du_tila term
    double Gaux;
    Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
    Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
    Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
    
    //RHS+=-Dv
    rRightHandSideVector[0] -= density*Area*Gaux * msN[0];
    rRightHandSideVector[1] -= density*Area*Gaux * msN[1];
    rRightHandSideVector[2] -= density*Area*Gaux * msN[2];

    //RHS = +tau*nablaN*f, we reuse aux
    //ms_aux0 stores ff_gauss;

    ms_aux0=0.33333333333333333*(ff0+ff1+ff2);
    //ms_aux1 - is the product of: (nabla q, f)
    ms_aux1[0]=msDN_DX(0,0)*ms_aux0[0]+msDN_DX(0,1)*ms_aux0[1];
    ms_aux1[1]=msDN_DX(1,0)*ms_aux0[0]+msDN_DX(1,1)*ms_aux0[1];
    ms_aux1[2]=msDN_DX(2,0)*ms_aux0[0]+msDN_DX(2,1)*ms_aux0[1];
    rRightHandSideVector += tau*density*Area*ms_aux1;


    //RHS += -tau*nablaN*du_gausspoint/dt
    //we reuse ms_vel_gauss to store the accelerations( (u_n - u_n-1)/dt)

    ms_vel_gauss[0]=0.33333333333*(fv0[0]+fv1[0]+fv2[0]-fv0_old[0]-fv1_old[0]-fv2_old[0])/dt;
    ms_vel_gauss[1]=0.33333333333*(fv0[1]+fv1[1]+fv2[1]-fv0_old[1]-fv1_old[1]-fv2_old[1])/dt;

    //and now we reuse ms_aux1

    ms_aux1=prod(msDN_DX,ms_vel_gauss);
    	
    noalias(rRightHandSideVector) -= tau*density*Area*ms_aux1;


    KRATOS_CATCH("")
      }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::FinalFractionalStep(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      KRATOS_THROW_ERROR(std::logic_error,  "METHOD NOT IMPL inside the element Final Fractional Step is done within the low_mach strategy.. " , "");

    KRATOS_CATCH("")
      }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
  {
    //if the VAr is NODAL_MASS, we calculate the lumped mass
    if(rVariable == NODAL_MASS)
      {
        CalculateLumpedMass();
      }
    else
      KRATOS_THROW_ERROR(std::logic_error,  "You are doing something wrong  FCT calculate... of nodal_mass with wring parameters.. " , "");
  }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::Calculate(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
  {

    if(rVariable == VELOCITY && Output[0]==1.0)
      //we use "Output" as a switch between 1st Frac Step and last Frac Step(Correction Step)
      {
        //here the residual will be temporarily written
        Vector TmpRhs(6);
        // first we write the Galerkin contributions to the momentum residual
        CalculateGalerkinMomentumResidual(TmpRhs);
        //and now the stabilization terms added
        double dt = rCurrentProcessInfo[DELTA_TIME];
        CalculateRHSVector(TmpRhs,  dt);

      }
    else if(rVariable == VELOCITY && Output[0]==2.0)
      {
        FinalFractionalStep(rCurrentProcessInfo);
      }
    else
      {
        KRATOS_THROW_ERROR(std::logic_error,  "You are doing something wrong in ur fractional step.... " , "");
      }
  }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().PointsNumber();


    if(rResult.size() != number_of_nodes)
      rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
      {
        rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
      }
    KRATOS_CATCH("")

      }
  //************************************************************************************
  //************************************************************************************
  void Fluid2DGLS_expl::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(ElementalDofList.size() != number_of_nodes)
      ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
      {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
      }
    KRATOS_CATCH("");
  }



} // Namespace Kratos


