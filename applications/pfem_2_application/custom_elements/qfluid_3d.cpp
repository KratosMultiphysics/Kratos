//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//


// System includes
#define QCOMP

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/qfluid_3d.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

  QFluid3D::QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)

  {
    //DO NOT ADD DOFS HERE!!!
  }


  QFluid3D::QFluid3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
    //InitializeAuxiliaries();
  }

  Element::Pointer QFluid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    KRATOS_TRY
      return Element::Pointer(new QFluid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
  }

  QFluid3D::~QFluid3D()
  {
  }


  void QFluid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber == 1) //first step of the fractional step solution
      {
        Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
      }
    else if (FractionalStepNumber == 4)//second step of the fractional step solution
      {
        Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
      }

    KRATOS_CATCH("")
      }


  void QFluid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_ERROR<< "method not implemented";
  }


  //calculation by component of the fractional step VELOCITY corresponding to the first stage
  void QFluid3D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    const unsigned int number_of_points = 4;
    const unsigned int dim = 3;
    unsigned int matsize = number_of_points *dim;

    //getting data for the given geometry
    BoundedMatrix<double,4,3> msDN_DX;
    BoundedMatrix<double,12,12> msMass= ZeroMatrix(12,12);
    //msMass=ZeroMatrix(12,12);

    array_1d<double,4> msN; //dimension = number of nodes
    BoundedMatrix<double,12,12> rDampMatrix = ZeroMatrix(12,12);

    BoundedMatrix<double,6,12> msB = ZeroMatrix(6,12);
    BoundedMatrix<double,6,6> ms_constitutive_matrix;
    array_1d<double,4> temp_vec_np;

    BoundedMatrix<double,6,12> ms_temp;
    double Volume;


    if(rLeftHandSideMatrix.size1() != matsize)
      rLeftHandSideMatrix.resize(matsize,matsize,false);

    if(rRightHandSideVector.size() != matsize)
      rRightHandSideVector.resize(matsize,false);

    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);

    const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);
    const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

    //calculating viscosity
    double nu = 0.25*(nu0 + nu1 + nu2 + nu3);


    double density = 0.25*(rho0 + rho1 + rho2 + rho3);

#if defined(QCOMP)
    //KRATOS_WATCH(density);
    //density=1000.0;
#else
    density = this->GetValue(DENSITY);
    nu=0.0;
#endif

    double dt = rCurrentProcessInfo[DELTA_TIME];

    msMass(0,0) = 0.25;
    msMass(1,1) = 0.25;
    msMass(2,2) = 0.25;
    msMass(3,3) = 0.25;
    msMass(4,4) = 0.25;
    msMass(5,5) = 0.25;
    msMass(6,6) = 0.25;
    msMass(7,7) = 0.25;
    msMass(8,8) = 0.25;
    msMass(9,9) = 0.25;
    msMass(10,10) = 0.25;
    msMass(11,11) = 0.25;

    CalculateViscousMatrix(rLeftHandSideMatrix, msDN_DX, nu, dt);

    noalias(rLeftHandSideMatrix) += msMass * density / dt;

    rLeftHandSideMatrix *= Volume;

    const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double,3>& force3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);

    array_1d<double,4> rhs_aux;
    for( unsigned int component_index = 0; component_index < dim; component_index++)
      {
        //external forces (component)
        double force_component = 0.25*(force0[component_index] + force1[component_index] + force2[component_index] + force3[component_index]);
        noalias(rhs_aux) = (force_component ) * msN * density;

        //adding pressure gradient (integrated by parts)
        double p_avg = p0old + p1old + p2old + p3old;
        p_avg *= 0.25;

#if defined(QCOMP)
        rhs_aux[0] += msDN_DX(0,component_index)*p_avg;
        rhs_aux[1] += msDN_DX(1,component_index)*p_avg;
        rhs_aux[2] += msDN_DX(2,component_index)*p_avg;
        rhs_aux[3] += msDN_DX(3,component_index)*p_avg;

#else
	//nada
#endif
	noalias(temp_vec_np) = ZeroVector(4);

	for(unsigned int iii = 0; iii<number_of_points; iii++)
	  {
	    const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
	    temp_vec_np[iii] += v[component_index] * density / dt;
	  }

        noalias(rhs_aux) += prod(msMass,temp_vec_np) ;

	//writing the rhs_aux in its place
        for( unsigned int i = 0; i < number_of_points; i++)
	  {
            rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
	  }
      }

    //multiplying by area
    rRightHandSideVector *= Volume;

    //LHS dirichlet contribution
    Vector fvvect(12);
    for( unsigned int component_index = 0; component_index < dim; component_index++)
      {
        fvvect[0 + component_index] = fv0[component_index];
        fvvect[3 + component_index] = fv1[component_index];
        fvvect[6 + component_index] = fv2[component_index];
        fvvect[9 + component_index] = fv3[component_index];
      }
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);


    KRATOS_CATCH("");
  }

  void QFluid3D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;


    unsigned int number_of_points = 4;

    if (rLeftHandSideMatrix.size1() != number_of_points)
      rLeftHandSideMatrix.resize(number_of_points, number_of_points);

    if (rRightHandSideVector.size() != number_of_points)
      rRightHandSideVector.resize(number_of_points);

    array_1d<double, 3 > vel_gauss;
    array_1d<double, 3 > aux;
    array_1d<double, 4 > temp_vec_np;

    //getting data for the given geometry
    double Volume;
    array_1d<double, 4 > N;
    BoundedMatrix<double, 4, 3 > DN_DX;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

    const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX);
    //const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

    const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX);
    //const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

    const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX);
    //const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

    const array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    const double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
    const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSUREAUX);
    //const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

    //calculating avergage density and viscosity
    //double nu = 0.25 * (nu0 + nu1 + nu2 + nu3);
    double density = 0.25 * (rho0 + rho1 + rho2 + rho3);
    //double nu=0.0;
    //#if defined(QCOMP)
    //nu=0.0;

    //#else
    density = this->GetValue(DENSITY);
    //#endif


    double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1] + vel_gauss[2] * vel_gauss[2];
    norm_u = sqrt(norm_u);
    //double tau = 0.0;

    double dt = rCurrentProcessInfo[DELTA_TIME];
    //dt=0.005;

    noalias(rLeftHandSideMatrix) = (dt / density) * prod(DN_DX, trans(DN_DX));


    //calculation of the RHS
    // RHS = -G*vfrac
    double Gaux;
    Gaux = DN_DX(0, 0) * fv0[0] + DN_DX(0, 1) * fv0[1] + DN_DX(0, 2) * fv0[2];
    Gaux += DN_DX(1, 0) * fv1[0] + DN_DX(1, 1) * fv1[1] + DN_DX(1, 2) * fv1[2];
    Gaux += DN_DX(2, 0) * fv2[0] + DN_DX(2, 1) * fv2[1] + DN_DX(2, 2) * fv2[2];
    Gaux += DN_DX(3, 0) * fv3[0] + DN_DX(3, 1) * fv3[1] + DN_DX(3, 2) * fv3[2];



    //
    vel_gauss[0] =  0.25*(fv0[0]+fv1[0]+fv2[0]+fv3[0] );
    vel_gauss[1] =  0.25*(fv0[1]+fv1[1]+fv2[1]+fv3[1]);
    vel_gauss[2] =  0.25*(fv0[2]+fv1[2]+fv2[2]+fv3[2]);

    rRightHandSideVector[0] +=(DN_DX(0,0) * vel_gauss[0] + DN_DX(0,1) * vel_gauss[1]+ DN_DX(0,2) * vel_gauss[2]);
    rRightHandSideVector[1] +=(DN_DX(1,0) * vel_gauss[0] + DN_DX(1,1) * vel_gauss[1] + DN_DX(1,2) * vel_gauss[2]);
    rRightHandSideVector[2] +=(DN_DX(2,0) * vel_gauss[0] + DN_DX(2,1) * vel_gauss[1] + DN_DX(2,2) * vel_gauss[2]);
    rRightHandSideVector[3] +=(DN_DX(3,0) * vel_gauss[0] + DN_DX(3,1) * vel_gauss[1] + DN_DX(3,2) * vel_gauss[2]);

    //dirichlet contribution
    temp_vec_np[0] = p0;
    temp_vec_np[1] = p1;
    temp_vec_np[2] = p2;
    temp_vec_np[3] = p3;
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);

    // RHS += dt * L * pold
    temp_vec_np[0] = p0old;
    temp_vec_np[1] = p1old;
    temp_vec_np[2] = p2old;
    temp_vec_np[3] = p3old;

    noalias(vel_gauss) = prod(trans(DN_DX), temp_vec_np);
    noalias(rRightHandSideVector) += (dt / density) * prod(DN_DX, vel_gauss);

    //multiplicating by the Volume
    rLeftHandSideMatrix *= Volume;
    rRightHandSideVector *= Volume;

    //adding contributions to nodal Volumes following the corresponding lumping term
    double nodal_contrib = 0.25 * Volume * density;

    double& m0 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
    m0 += nodal_contrib;

    double& m1 = GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
    m1 += nodal_contrib;

    double& m2 = GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
    m2 += nodal_contrib;

    double& m3 = GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
    m3 += nodal_contrib;


    KRATOS_CATCH("");
  }

  void QFluid3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY
      int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];


    if(FractionalStepNumber  == 5) //calculation of stabilization terms
      {

	//const unsigned int number_of_points = 4;
	//const unsigned int dim = 3;
	//unsigned int matsize = number_of_points *dim;

	//getting data for the given geometry
	BoundedMatrix<double,4,3> msDN_DX;
	//BoundedMatrix<double,12,12> msMass= ZeroMatrix(12,12);
	//msMass=ZeroMatrix(12,12);

	array_1d<double,12> GalerkinRHS = ZeroVector(12); //dimension = number of nodes

	array_1d<double,4> msN; //dimension = number of nodes

	double Volume;

	GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume);

	const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);


	const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

	const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);


	const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);


	double density = 0.25*(rho0 + rho1 + rho2 + rho3);

#if defined(QCOMP)

#else
	density = this->GetValue(DENSITY);
	nu=0.0;
#endif

	//double dt = CurrentProcessInfo[DELTA_TIME];

	const array_1d<double,3> body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +	GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE) + GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));
        unsigned int number_of_nodes=4;

        for(unsigned int i = 0; i<number_of_nodes; i++)
	  {
            GalerkinRHS[i*3] += body_force[0] * Volume * 0.25 * density;
            GalerkinRHS[i*3+1] += body_force[1] * Volume * 0.25 * density;
            GalerkinRHS[i*3+2] += body_force[2] * Volume * 0.25 * density;
	  }

	GetGeometry()[0].SetLock();
        array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
        rhs0[0] += GalerkinRHS[0] ;
        rhs0[1] += GalerkinRHS[1] ;
        rhs0[2] += GalerkinRHS[2] ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
        rhs1[0] += GalerkinRHS[3] ;
        rhs1[1] += GalerkinRHS[4] ;
        rhs1[2] += GalerkinRHS[5] ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
        rhs2[0] += GalerkinRHS[6] ;
        rhs2[1] += GalerkinRHS[7] ;
        rhs2[2] += GalerkinRHS[8] ;
        GetGeometry()[2].UnSetLock();

	GetGeometry()[3].SetLock();
        array_1d<double,3>& rhs3 = GetGeometry()[3].FastGetSolutionStepValue(FORCE);
        rhs3[0] += GalerkinRHS[9] ;
        rhs3[1] += GalerkinRHS[10] ;
        rhs3[2] += GalerkinRHS[11] ;
        GetGeometry()[3].UnSetLock();

	double nodal_contrib = 0.25 * Volume;

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[2].UnSetLock();

        GetGeometry()[3].SetLock();
        GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[3].UnSetLock();
      }
    else if(FractionalStepNumber == 6) //calculation of velocities
      {

    	double Area;
    	array_1d<double, 4 > msN;
    	BoundedMatrix<double, 4, 3 > DN_DX;

        //double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,msN,Area);

	const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    	const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
	const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    	const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

	double density = 0.25*(rho0 + rho1 + rho2 + rho3);
	density=1000.0;

	density = this->GetValue(DENSITY);
	KRATOS_ERROR<<  "method not implemented" ;

        double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
        double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX);

        double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
        double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX);

        double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
        double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX);

	double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
        double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSUREAUX);

	double p_grad =  DN_DX(0, 0) * (p0 - p0old) + DN_DX(1, 0) * (p1 - p1old) + DN_DX(2, 0) * (p2 - p2old) + DN_DX(3, 0) * (p3 - p3old);
       	double p_grad1 = DN_DX(0, 1) * (p0 - p0old) + DN_DX(1, 1) * (p1 - p1old) + DN_DX(2, 1) * (p2 - p2old) + DN_DX(3, 1) * (p3 - p3old);
       	double p_grad2 = DN_DX(0, 2) * (p0 - p0old) + DN_DX(1, 2) * (p1 - p1old) + DN_DX(2, 2) * (p2 - p2old) + DN_DX(3, 2) * (p3 - p3old);

	//geom[0].SetLock();
	GetGeometry()[0].FastGetSolutionStepValue(FORCE_X) -= p_grad * 0.25 * Area;
	GetGeometry()[0].FastGetSolutionStepValue(FORCE_Y) -= p_grad1 * 0.25 * Area;
	GetGeometry()[0].FastGetSolutionStepValue(FORCE_Z) -= p_grad2 * 0.25 * Area;
	//geom[0].UnSetLock();

	//geom[1].SetLock();
	GetGeometry()[1].FastGetSolutionStepValue(FORCE_X) -= p_grad * 0.25 * Area;
	GetGeometry()[1].FastGetSolutionStepValue(FORCE_Y) -= p_grad1 * 0.25 * Area;
	GetGeometry()[1].FastGetSolutionStepValue(FORCE_Z) -= p_grad2 * 0.25 * Area;
	//geom[1].UnSetLock();

	//geom[2].SetLock();
	GetGeometry()[2].FastGetSolutionStepValue(FORCE_X) -= p_grad * 0.25 * Area;
	GetGeometry()[2].FastGetSolutionStepValue(FORCE_Y) -= p_grad1 * 0.25 * Area;
	GetGeometry()[2].FastGetSolutionStepValue(FORCE_Z) -= p_grad2 * 0.25 * Area;
	//geom[2].UnSetLock();


	//geom[3].SetLock();
	GetGeometry()[3].FastGetSolutionStepValue(FORCE_X) -= p_grad * 0.25 * Area;
	GetGeometry()[3].FastGetSolutionStepValue(FORCE_Y) -= p_grad1 * 0.25 * Area;
	GetGeometry()[3].FastGetSolutionStepValue(FORCE_Z) -= p_grad2 * 0.25 * Area;
	//geom[3].UnSetLock();

        double nodal_contrib = 0.25 * Area;

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[2].UnSetLock();

        GetGeometry()[3].SetLock();
        GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[3].UnSetLock();

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
        GetGeometry()[2].UnSetLock();

        GetGeometry()[3].SetLock();
        GetGeometry()[3].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib;
        GetGeometry()[3].UnSetLock();

      }
    else if (FractionalStepNumber == 7)
      {
    	array_1d<double,4> temp_vec_np; //dimension = number of nodes
        array_1d<double,4> GalerkinRHS = ZeroVector(4); //dimension = number of nodes
        BoundedMatrix<double,4,3> msDN_DX; // = ZeroMatrix(3,2);
        array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes


        double dt = CurrentProcessInfo[DELTA_TIME];

        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

	Matrix Mass(4,4);
        Mass=ZeroMatrix(4,4);
        Mass(0,0) = 0.25;
        Mass(1,1) = 0.25;
        Mass(2,2) = 0.25;
        Mass(3,3) = 0.25;

	double nodal_contrib = 0.25 * Area;

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[2].UnSetLock();

	GetGeometry()[3].SetLock();
        GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[3].UnSetLock();

        const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
        //const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        //const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
	const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);

        const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
        double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
        //const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        //const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
	const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);

        const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
        double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
        //const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        //const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
	const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);

	const array_1d<double,3>& vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
        double p_n3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);
        //const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
        //const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
	const double k3 = GetGeometry()[3].FastGetSolutionStepValue(BULK_MODULUS);

	//double density = 0.25*(rho0 + rho1 + rho2 + rho3);

	double bulk_modulus = 0.25*(k0 + k1 + k2 + k3) * dt;

	bulk_modulus=  -1.0 * dt * 5000.0;

	temp_vec_np[0]= p_n0;
	temp_vec_np[1]= p_n1;
	temp_vec_np[2]= p_n2;
	temp_vec_np[3]= p_n3;

        noalias(GalerkinRHS) = prod(Mass, temp_vec_np);
	//KRATOS_WATCH(GalerkinRHS);

        //double E_over_R = 24466.81;
    	//double C = 1.4e10;

	GalerkinRHS *= Area ;

        double Gaux;
	Gaux =  msDN_DX(0,0) * vel0[0] + msDN_DX(0,1) * vel0[1] + msDN_DX(0,2) * vel0[2];
	Gaux += msDN_DX(1,0) * vel1[0] + msDN_DX(1,1) * vel1[1] + msDN_DX(1,2) * vel1[2];
	Gaux += msDN_DX(2,0) * vel2[0] + msDN_DX(2,1) * vel2[1] + msDN_DX(2,2) * vel2[2];
	Gaux += msDN_DX(3,0) * vel3[0] + msDN_DX(3,1) * vel3[1] + msDN_DX(3,2) * vel3[2];


        double t1 = GetGeometry()[0].FastGetSolutionStepValue(YCH4);

        double t2 = GetGeometry()[1].FastGetSolutionStepValue(YCH4);

        double t3 = GetGeometry()[2].FastGetSolutionStepValue(YCH4);

	double t4 = GetGeometry()[3].FastGetSolutionStepValue(YCH4);

	double temp=t1 + t2 + t3 + t4;
	temp *= 0.25;
	if(temp>1000.0) temp=1000.0;

	//double E_over_R = 28961.49;
    	//double C = 1.19e15;

	constexpr double E_over_R = 24400.0;//28961.49;
    	constexpr double C = 2.18e12;//1.19e15;
        double aux_var= C * exp(-E_over_R/(temp));
	//////////////
        GalerkinRHS[0] += bulk_modulus * Area * Gaux * 0.25 + 1.0 * bulk_modulus * Area * aux_var * 0.25 ;

        GalerkinRHS[1] += bulk_modulus * Area * Gaux * 0.25 + 1.0 * bulk_modulus * Area * aux_var * 0.25 ;

        GalerkinRHS[2] += bulk_modulus * Area * Gaux * 0.25 + 1.0 * bulk_modulus * Area * aux_var * 0.25 ;

        GalerkinRHS[3] += bulk_modulus * Area * Gaux * 0.25 + 1.0 * bulk_modulus * Area * aux_var * 0.25 ;

	GetGeometry()[0].SetLock();
        double & rhs0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSUREAUX);
        rhs0 += GalerkinRHS[0] ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        double & rhs1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSUREAUX);
        rhs1 += GalerkinRHS[1] ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        double & rhs2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSUREAUX);
        rhs2 += GalerkinRHS[2] ;
        GetGeometry()[2].UnSetLock();

	GetGeometry()[3].SetLock();
        double & rhs3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSUREAUX);
        rhs3 += GalerkinRHS[3] ;
        GetGeometry()[3].UnSetLock();
	//KRATOS_WATCH(GalerkinRHS);
      }

    KRATOS_CATCH("");
  }

  void QFluid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dim = 3;

    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber == 1) //step 1
      {
        if(rResult.size() != number_of_nodes*dim)
	  rResult.resize(number_of_nodes*dim,false);

        for (unsigned int i=0; i<number_of_nodes; i++)
	  {
	    rResult[i*dim] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
	    rResult[i*dim+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
	    rResult[i*dim+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
	  }
      }
    else if(FractionalStepNumber == 4) // pressure correction step
      {
        if(rResult.size() != number_of_nodes)
	  rResult.resize(number_of_nodes,false);
        for (unsigned int i=0; i<number_of_nodes; i++)
	  rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
      }

  }

  void QFluid3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 3;

    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber == 1) //step 1
      {
        if(ElementalDofList.size() != number_of_nodes*dim)
	  ElementalDofList.resize(number_of_nodes*dim);

        for (unsigned int i=0; i<number_of_nodes; i++)
	  {
            ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(VELOCITY_X);
            ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
            ElementalDofList[i*dim+2] = GetGeometry()[i].pGetDof(VELOCITY_Z);
	  }
      }
    else if(FractionalStepNumber == 4) // pressure correction step
      {
        if(ElementalDofList.size() != number_of_nodes)
	  ElementalDofList.resize(number_of_nodes);

        for (unsigned int i=0; i<number_of_nodes; i++)
	  ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
      }
  }

  void QFluid3D::CalculateViscousMatrix(MatrixType& K, const BoundedMatrix<double,4,3>& DN_DX, const double& nu,const double& dtt )
  {

    double viscosity=nu;
    double dt=dtt;
    double k=25000.0 * dt;

#if defined(QCOMP)

  k=5000.0 * dt;
#else
    k=0.0;
#endif


    //k=0.0;
    const unsigned int number_of_nodes = 4;
    const unsigned int dim = 3;

    BoundedMatrix<double,6,12> msB = ZeroMatrix(6,12);
    BoundedMatrix<double,6,6> ms_constitutive_matrix;
    BoundedMatrix<double,6,6> msCapx;
    BoundedMatrix<double,6,12> ms_temp;
    BoundedMatrix<double,12,12> rDampMatrix;
    BoundedMatrix<double,6,12> B;


    //unsigned int start;

    for(unsigned int i = 0; i<number_of_nodes; i++)
      {
	unsigned int start = dim*i;

	msB(0,start) =	DN_DX(i,0);
	msB(1,start+1)=	DN_DX(i,1);
	msB(2,start+2)= DN_DX(i,2);
	msB(3,start) =	DN_DX(i,1);
	msB(3,start+1) = DN_DX(i,0);
	msB(4,start) =	DN_DX(i,2);
	msB(4,start+2) = DN_DX(i,0);
	msB(5,start+1)= DN_DX(i,2);
	msB(5,start+2) = DN_DX(i,1);
      }


    //const double& a = nu;
    //constitutive tensor
    ms_constitutive_matrix(0,0) = (4.0/3.0)*viscosity;
    ms_constitutive_matrix(0,1) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(0,2) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(0,3) = 0.0;
    ms_constitutive_matrix(0,4) = 0.0;
    ms_constitutive_matrix(0,5) = 0.0;

    ms_constitutive_matrix(1,0) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(1,1) = 4.0/3.0*viscosity;
    ms_constitutive_matrix(1,2) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(1,3) = 0.0;
    ms_constitutive_matrix(1,4) = 0.0;
    ms_constitutive_matrix(1,5) = 0.0;

    ms_constitutive_matrix(2,0) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(2,1) = -2.0/3.0*viscosity;
    ms_constitutive_matrix(2,2) = 4.0/3.0*viscosity;
    ms_constitutive_matrix(2,3) = 0.0;
    ms_constitutive_matrix(2,4) = 0.0;
    ms_constitutive_matrix(2,5) = 0.0;

    ms_constitutive_matrix(3,0) = 0.0;
    ms_constitutive_matrix(3,1) = 0.0;
    ms_constitutive_matrix(3,2) = 0.0;
    ms_constitutive_matrix(3,3) = viscosity;
    ms_constitutive_matrix(3,4) = 0.0;
    ms_constitutive_matrix(3,5) = 0.0;

    ms_constitutive_matrix(4,0) = 0.0;
    ms_constitutive_matrix(4,1) = 0.0;
    ms_constitutive_matrix(4,2) = 0.0;
    ms_constitutive_matrix(4,3) = 0.0;
    ms_constitutive_matrix(4,4) = viscosity;
    ms_constitutive_matrix(4,5) = 0.0;

    ms_constitutive_matrix(5,0) = 0.0;
    ms_constitutive_matrix(5,1) = 0.0;
    ms_constitutive_matrix(5,2) = 0.0;
    ms_constitutive_matrix(5,3) = 0.0;
    ms_constitutive_matrix(5,4) = 0.0;
    ms_constitutive_matrix(5,5) = viscosity;

    ms_temp = prod( ms_constitutive_matrix , msB);

    noalias(K) = prod( trans(msB) , ms_temp);

    //KRATOS_WATCH(K);
    //other matrix
    msCapx(0,0) = 1.0 * k;
    msCapx(0,1) = 1.0 * k;
    msCapx(0,2) = 1.0 * k;
    msCapx(0,3) = 0.0;
    msCapx(0,4) = 0.0;
    msCapx(0,5) = 0.0;

    msCapx(1,0) = 1.0 * k;
    msCapx(1,1) = 1.0 * k;
    msCapx(1,2) = 1.0 * k;
    msCapx(1,3) = 0.0;
    msCapx(1,4) = 0.0;
    msCapx(1,5) = 0.0;

    msCapx(2,0) = 1.0 * k;
    msCapx(2,1) = 1.0 * k;
    msCapx(2,2) = 1.0 * k;
    msCapx(2,3) = 0.0;
    msCapx(2,4) = 0.0;
    msCapx(2,5) = 0.0;

    msCapx(3,0) = 0.0;
    msCapx(3,1) = 0.0;
    msCapx(3,2) = 0.0;
    msCapx(3,3) = 0.0;
    msCapx(3,4) = 0.0;
    msCapx(3,5) = 0.0;

    msCapx(4,0) = 0.0;
    msCapx(4,1) = 0.0;
    msCapx(4,2) = 0.0;
    msCapx(4,3) = 0.0;
    msCapx(4,4) = 0.0;
    msCapx(4,5) = 0.0;

    msCapx(5,0) = 0.0;
    msCapx(5,1) = 0.0;
    msCapx(5,2) = 0.0;
    msCapx(5,3) = 0.0;
    msCapx(5,4) = 0.0;
    msCapx(5,5) = 0.0;
    ms_temp = prod( msCapx , msB);

    noalias(K) += prod( trans(msB) , ms_temp);

}

  //performs the Kroneker product of the Reduced Matrix with the identity matrix of
  //size "dimension" ADDING to the destination matrix
  inline void  QFluid3D::ExpandAndAddReducedMatrix(
						   MatrixType& Destination,
						   BoundedMatrix<double,4,4>& ReducedMatrix,
						   const unsigned int dimension)
  {
    KRATOS_TRY
      for (unsigned int i=0; i<ReducedMatrix.size2(); i++)
	{
	  int rowindex = i*dimension;
	  for (unsigned int j=0; j<ReducedMatrix.size2(); j++)
	    {
	      unsigned int colindex = j*dimension;
	      for(unsigned int ii=0; ii<dimension; ii++)
		Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
	    }
	}
    KRATOS_CATCH("")
      }

  inline double QFluid3D::CalculateH(double Volume)
  {
    double h = pow(6.00*Volume,0.3333333);


    return h;
  }

} // Namespace Kratos
