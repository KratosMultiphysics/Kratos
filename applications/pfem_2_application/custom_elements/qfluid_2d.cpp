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


//#define FIRST

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/qfluid_2d.h"
#include "utilities/math_utils.h"
#include "pfem_2_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


QFluid2D::QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry)
  : Element(NewId, pGeometry)
{
  //DO NOT ADD DOFS HERE!!!
}


  QFluid2D::QFluid2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
  }

  Element::Pointer QFluid2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    KRATOS_TRY

      return Element::Pointer(new QFluid2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
  }

  QFluid2D::~QFluid2D()
  {
  }


  void QFluid2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber ==1 ) //first step of the fractional step solution
      {
        Stage1(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);//,ComponentIndex);
      }
    else if (FractionalStepNumber == 4)//second step of the fractional step solution
      {
        //Stage2(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
      }


void QFluid2D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;

  unsigned int number_of_points = 3;
  unsigned int dim = 2;
  unsigned int matsize = number_of_points *dim;


  if(rLeftHandSideMatrix.size1() != matsize)
    rLeftHandSideMatrix.resize(matsize,matsize,false);

  if(rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false);
  //KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");

  BoundedMatrix<double, 3, 3 > WorkMatrix;

  array_1d<double, 3 > N;
  array_1d<double, 2 > vel_gauss;
  array_1d<double, 3 > temp_vec_np;
  array_1d<double, 3 > u_DN;
  array_1d<double, 3 > aux0;
  array_1d<double, 3 > aux1;
  array_1d<double, 3 > aux2;

  //BoundedMatrix<double,3,3> msaux_matrix;
  BoundedMatrix<double,6,6> msMassFactors;
  BoundedMatrix<double,3,2> msDN_DX;
  array_1d<double,3> msN; //dimension = number of nodes
  array_1d<double,2> ms_vel_gauss; //dimesion coincides with space dimension
  //array_1d<double,3> temp_vec_np; //dimension = number of nodes
  array_1d<double,3> ms_u_DN; //dimension = number of nodes
  //getting data for the given geometry

  if(rLeftHandSideMatrix.size1() != matsize)
    rLeftHandSideMatrix.resize(matsize,matsize,false);

  if(rRightHandSideVector.size() != matsize)
    rRightHandSideVector.resize(matsize,false);

  double dt = rCurrentProcessInfo[DELTA_TIME];

  //getting data for the given geometry
  double Area;
  //CalculateGeometryData(msDN_DX,msN,Area);
  GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

  //getting the velocity vector on the nodes

  //getting the velocity on the nodes
  const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
  //double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
  double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
  const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
  const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
  const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);

  const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
  //double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
  double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
  const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
  const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
  const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);

  const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
  //double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
  double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
  const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
  const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
  const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);


  //calculating viscosity
  double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
  double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
  double kk= 0.3333333333333333333333*(k0 + k1 + k2 );

  //getting the BDF2 coefficients (not fixed to allow variable time step)
  //the coefficients INCLUDE the time step
  //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

  //double h = sqrt(2.00*Area);
  //double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1];
  //norm_u = sqrt(norm_u);
  //double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );

  //Matrix msMassFactors(6,6);
  //    Matrix msaux_matrix(3,3);
  BoundedMatrix<double,3,3> msaux_matrix;

  msMassFactors(0,0) = 1.00/3.00;
  msMassFactors(0,1) = 0.00;
  msMassFactors(0,2) = 0.00;
  msMassFactors(0,3) = 0.00;
  msMassFactors(0,4) = 0.00;
  msMassFactors(0,5) = 0.00;

  msMassFactors(1,0) = 0.0;
  msMassFactors(1,1) = 1.00/3.00;
  msMassFactors(1,2) = 0.00;
  msMassFactors(1,3) = 0.00;
  msMassFactors(1,4) = 0.00;
  msMassFactors(1,5) = 0.00;

  msMassFactors(2,0) = 0.00;
  msMassFactors(2,1) = 0.00;
  msMassFactors(2,2) = 1.00/3.00;
  msMassFactors(2,3) = 0.00;
  msMassFactors(2,4) = 0.00;
  msMassFactors(2,5) = 0.00;

  msMassFactors(3,0) = 0.00;
  msMassFactors(3,1) = 0.00;
  msMassFactors(3,2) = 0.00;
  msMassFactors(3,3) = 1.00/3.00;
  msMassFactors(3,4) = 0.00;
  msMassFactors(3,5) = 0.00;

  msMassFactors(4,0) = 0.00;
  msMassFactors(4,1) = 0.00;
  msMassFactors(4,2) = 0.00;
  msMassFactors(4,3) = 0.00;
  msMassFactors(4,4) = 1.00/3.00;
  msMassFactors(4,5) = 0.00;

  msMassFactors(5,0) = 0.00;
  msMassFactors(5,1) = 0.00;
  msMassFactors(5,2) = 0.00;
  msMassFactors(5,3) = 0.00;
  msMassFactors(5,4) = 0.00;
  msMassFactors(5,5) = 1.00/3.00;

  density=905.0;
  kk *= dt;
  CalculateViscousMatrix(rLeftHandSideMatrix, msDN_DX, nu, dt, kk);

  //KRATOS_WATCH(rLeftHandSideMatrix);
  //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX

  //INERTIA CONTRIBUTION
  noalias(rLeftHandSideMatrix) += msMassFactors / dt * density;
  //KRATOS_WATCH(rLeftHandSideMatrix);

  //multiplication by the area
  rLeftHandSideMatrix *= Area;

  //CALCULATION OF THE RHS
  const array_1d<double,3>& force0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
  const array_1d<double,3>& force1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
  const array_1d<double,3>& force2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


  array_1d<double,3> rhs_aux;

  for( unsigned int component_index = 0; component_index < dim; component_index++)
    {
      //external forces (component)
      double force_component = 0.3333333333333333*(force0[component_index] + force1[component_index] + force2[component_index]);
      noalias(rhs_aux) = (force_component ) * msN * density;

      //adding pressure gradient (integrated by parts)
      double p_avg = p0old + p1old + p2old;
      p_avg *= 0.3333333333333333 ;
      /*rhs_aux[0] += msDN_DX(0,component_index)*p_avg * 1.0 ;
        rhs_aux[1] += msDN_DX(1,component_index)*p_avg * 1.0 ;
        rhs_aux[2] += msDN_DX(2,component_index)*p_avg * 1.0 ;*/

      double p_grad = msDN_DX(0, component_index) * 1.0 * p0old + msDN_DX(1, component_index) * 1.0*  p1old + msDN_DX(2, component_index) *1.0 * p2old;
      //p_grad /= density;
      // RHS = Fext - grad*pn

      rhs_aux[0] -= p_grad * 0.3333333333333333;//msN[0] ;
      rhs_aux[1] -= p_grad * 0.3333333333333333 ; //msN[1] ;
      rhs_aux[2] -= p_grad * 0.3333333333333333; //msN[2] ;

      for(unsigned int iii = 0; iii<number_of_points; iii++)
        {
	  const array_1d<double,3>& v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1) );
	  temp_vec_np[iii] = v[component_index] / dt * density;
        }
      noalias(rhs_aux) += prod(msMassFactors,temp_vec_np) ;

      //writing the rhs_aux in its place
      for( unsigned int i = 0; i < number_of_points; i++)
        {
	  rRightHandSideVector[i*dim + component_index] = rhs_aux[i];
        }
    }
  rRightHandSideVector *= Area;

  //KRATOS_WATCH(rRightHandSideVector);
  //LHS stabilization contribution
  Vector fvvect(6);
  for( unsigned int component_index = 0; component_index < dim; component_index++)
    {
      fvvect[0 + component_index] = fv0[component_index];
      fvvect[2 + component_index] = fv1[component_index];
      fvvect[4 + component_index] = fv2[component_index];
    }
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,fvvect);


  KRATOS_CATCH("");
}



  void QFluid2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

      ///////////////////////NECESSARY LOCALS///////////////////////////////////////////
      BoundedMatrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
    BoundedMatrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
    array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
    BoundedMatrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
    BoundedMatrix<double,2,6> msConvOp = ZeroMatrix(2,6);
    BoundedMatrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
    array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
    array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,3> temp_vec_np = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
    BoundedMatrix<double,6,6> msAuxMat1 = ZeroMatrix(6,6);
    BoundedMatrix<double,6,3> msAuxMat2 = ZeroMatrix(6,3);
    BoundedMatrix<double,2,2> msGrad_ug = ZeroMatrix(2,2);
    array_1d<double,6> msStabMomRes = ZeroVector(6); //dimension = number of nodes
    BoundedMatrix<double,6,3> msGradOp = ZeroMatrix(6,3);

    BoundedMatrix<double,3,3> msWorkMatrix1 = ZeroMatrix(3,3);
    ///////////////////////////////////////////////////////////////////////////////////

    if(rRightHandSideVector.size() != 3)
    {
      rRightHandSideVector.resize(3,false);
    }

    double dt = rCurrentProcessInfo[DELTA_TIME];

    //fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
    //so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2
    const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& fv0_n_1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    //double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& fv1_n_1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    //double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& fv2_n_1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
    //old iteration can be used if we want to iterate between 1st and 2nd fractional steps
    //double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);


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
    double tau = 1.00 / ( 4.00*nu/(h*h) +2.00*norm_u/h);// + 1.0/dt);
    //tau=0.0;
    noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
    noalias(msWorkMatrix1) = (dt/2.0 + tau) * Area*msWorkMatrix;


    //////////////////////////////////////////////////////////
    ////////////		AND NOW RHS	//////////////////
    //////////////////////////////////////////////////////////

    //Dirichlet contribution  (that is: LHS*p_new)
    temp_vec_np[0] = p0;
    temp_vec_np[1] = p1;
    temp_vec_np[2] = p2;
    //LHS is already multiplied by AREA
    noalias(rRightHandSideVector) = -prod(msWorkMatrix1,temp_vec_np);

    //NOW RHS-=dt L p_old
    //changing the meaning of temp_vec_np
    temp_vec_np[0] = p0_old;
    temp_vec_np[1] = p1_old;
    temp_vec_np[2] = p2_old;

    noalias(rRightHandSideVector) += Area* dt/2.0 * (prod(msWorkMatrix,temp_vec_np)) ;

    //***************************************************************************

    //here we have the Du_tila term
    double Gaux;
    Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1];
    Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1];
    Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1];
    rRightHandSideVector[0] -= density*Area*Gaux * msN[0];
    rRightHandSideVector[1] -= density*Area*Gaux * msN[1];
    rRightHandSideVector[2] -= density*Area*Gaux * msN[2];

    ms_aux0=0.33333333333333333*(ff0+ff1+ff2);
    //ms_aux1 - is the product of: (nabla q, f)
    ms_aux1[0]=msDN_DX(0,0)*ms_aux0[0]+msDN_DX(0,1)*ms_aux0[1];
    ms_aux1[1]=msDN_DX(1,0)*ms_aux0[0]+msDN_DX(1,1)*ms_aux0[1];
    ms_aux1[2]=msDN_DX(2,0)*ms_aux0[0]+msDN_DX(2,1)*ms_aux0[1];

    ms_vel_gauss[0]=0.33333333333*(fv0[0]+fv1[0]+fv2[0]-fv0_old[0]-fv1_old[0]-fv2_old[0])/dt;
    ms_vel_gauss[1]=0.33333333333*(fv0[1]+fv1[1]+fv2[1]-fv0_old[1]-fv1_old[1]-fv2_old[1])/dt;

    //and now we reuse ms_aux1
    ms_aux1=prod(msDN_DX,ms_vel_gauss);

    //contains d_ug/dx
    //x-component of u derived with resp to x, remember msAuxVec stores velocity
    msGrad_ug(0,0)=msDN_DX(0,0)*msAuxVec[0]+msDN_DX(1,0)*msAuxVec[2]+msDN_DX(2,0)*msAuxVec[4];
    //x-component of u derived with resp to y
    msGrad_ug(0,1)=msDN_DX(0,1)*msAuxVec[0]+msDN_DX(1,1)*msAuxVec[2]+msDN_DX(2,1)*msAuxVec[4];
    //y-component of u derived with resp to x
    msGrad_ug(1,0)=msDN_DX(0,0)*msAuxVec[1]+msDN_DX(1,0)*msAuxVec[3]+msDN_DX(2,0)*msAuxVec[5];
    //y-component of u derived with resp to y
    msGrad_ug(1,1)=msDN_DX(0,1)*msAuxVec[1]+msDN_DX(1,1)*msAuxVec[3]+msDN_DX(2,1)*msAuxVec[5];

    array_1d<double,2> a;
    a[0]=0.33333333333333*(msAuxVec[0]+msAuxVec[2]+msAuxVec[4])*msGrad_ug(0,0)+0.33333333333333*(msAuxVec[1]+msAuxVec[3]+msAuxVec[5])*msGrad_ug(0,1);
    a[1]=0.33333333333333*(msAuxVec[0]+msAuxVec[2]+msAuxVec[4])*msGrad_ug(1,0)+0.33333333333333*(msAuxVec[1]+msAuxVec[3]+msAuxVec[5])*msGrad_ug(1,1);

   KRATOS_CATCH("")
     }
  void QFluid2D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      ///////////////////////NECESSARY LOCALS///////////////////////////////////////////
      BoundedMatrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
    BoundedMatrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
    array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
    BoundedMatrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
    BoundedMatrix<double,2,6> msConvOp = ZeroMatrix(2,6);
    BoundedMatrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
    array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
    array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
    array_1d<double,3> temp_vec_np = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux0 = ZeroVector(3); //dimension = number of nodes
    array_1d<double,3> ms_aux1 = ZeroVector(3); //dimension = number of nodes
    BoundedMatrix<double,6,6> msAuxMat1 = ZeroMatrix(6,6);
    BoundedMatrix<double,6,3> msAuxMat2 = ZeroMatrix(6,3);
    BoundedMatrix<double,2,2> msGrad_ug = ZeroMatrix(2,2);
    array_1d<double,6> msStabMomRes = ZeroVector(6); //dimension = number of nodes
    BoundedMatrix<double,6,3> msGradOp = ZeroMatrix(6,3);
    array_1d<double, 2 > vel_gauss;
    BoundedMatrix<double,3,3> Tres = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> msResta = IdentityMatrix(3,3);
    BoundedMatrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
    array_1d<double,3> ms_aux2 = ZeroVector(3); //dimension = number of nodes

    BoundedMatrix<double,3,3> aux_matrix = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> aux_matrix1 = ZeroMatrix(3,3);

    mThisIntegrationMethod= GeometryData::GI_GAUSS_1;

    if(rRightHandSideVector.size() != 3)
      {
        rLeftHandSideMatrix.resize(3,3,false);
        rRightHandSideVector.resize(3,false);
      }
    KRATOS_CATCH("")
      }

  void QFluid2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
  {
    KRATOS_TRY

      int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
    BoundedMatrix<double, 3, 3 > MassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
    BoundedMatrix<double, 3, 3 > WorkMatrix;
    BoundedMatrix<double, 3, 2 > DN_DX;

    BoundedMatrix<double, 6, 6 > rDampMatrix;


    array_1d<double, 3 > N;
    array_1d<double, 2 > vel_gauss;
    array_1d<double, 3 > temp_vec_np;
    array_1d<double, 3 > u_DN;
    array_1d<double, 3 > aux0;
    array_1d<double, 3 > aux1;
    array_1d<double, 3 > aux2;

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
    //KRATOS_ERROR(std::logic_error, "method not implemented", "");
    if (FractionalStepNumber == 5) //calculation of stabilization terms
      {
        //KRATOS_ERROR(std::logic_error, "method not implemented", "");
        ///////////////////////NECESSARY LOCALS///////////////////////////////////////////
        BoundedMatrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
        array_1d<double,6> GalerkinRHS = ZeroVector(6); //dimension = number of nodes
        BoundedMatrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
        array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
        BoundedMatrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
        BoundedMatrix<double,2,6> msConvOp = ZeroMatrix(2,6);
        BoundedMatrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
        array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
	BoundedMatrix<double,3,6> msB;
        BoundedMatrix<double,3,3> ms_constitutive_matrix;
        BoundedMatrix<double,3,6> ms_temp;
	//double dt = CurrentProcessInfo[DELTA_TIME];
	double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

        //getting the velocity on the nodes and other necessary variabless
        const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
        //double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);

        const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
        //const double k0 = GetGeometry()[0].FastGetSolutionStepValue(BULK_MODULUS);

        const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
        //double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
        const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
        //const double k1 = GetGeometry()[1].FastGetSolutionStepValue(BULK_MODULUS);

	const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
        //double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);

        const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
        //const double k2 = GetGeometry()[2].FastGetSolutionStepValue(BULK_MODULUS);

        //====================================================================
        //calculating viscosity and density
        double nu = 0.333333333333333333333333*(nu0 + nu1 + nu2 );
        double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );
        //double bulk_modulus = 0.3333333333333333333333*(k0 + k1 + k2 ) * dt;
        nu=0.0;
        density=1000.0;

        noalias(msWorkMatrix) = Area*  nu  * prod(msDN_DX,trans(msDN_DX));

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

        ms_adv_vel[0] = msN[0]*(vel0[0])+msN[1]*(vel1[0])+msN[2]*(vel2[0]);
        ms_adv_vel[1] = msN[0]*(vel0[1])+msN[1]*(vel1[1])+msN[2]*(vel2[1]);

        const array_1d<double,3> body_force = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+ GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +	GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
        unsigned int number_of_nodes=3;

        for(unsigned int i = 0; i<number_of_nodes; i++)
	  {
            GalerkinRHS[i*2] += body_force[0] * Area * 0.3333333333333 * density;
            GalerkinRHS[i*2+1] += body_force[1] * Area * 0.3333333333333 * density;
	  }

	double p_avg=0.333333333333*(p_n0+p_n1+p_n2)*Area;

        GalerkinRHS[0]+= 0.0 * msDN_DX(0,0)*p_avg ;
        GalerkinRHS[1]+= 0.0 * msDN_DX(0,1)*p_avg ;

        GalerkinRHS[2]+= 0.0 * msDN_DX(1,0)*p_avg ;
        GalerkinRHS[3]+= 0.0 * msDN_DX(1,1)*p_avg ;

        GalerkinRHS[4]+= 0.0 * msDN_DX(2,0)*p_avg ;
        GalerkinRHS[5]+= 0.0 * msDN_DX(2,1)*p_avg ;

 	double p_grad = DN_DX(0, 0) * p_n0 + DN_DX(1, 0) * p_n1 + DN_DX(2, 0) * p_n2;
        double p_grad1 = DN_DX(0, 1) * p_n0 + DN_DX(1, 1) * p_n1 + DN_DX(2, 1) * p_n2;

	p_avg = N[0] * p_n0 + N[1] * p_n1 + N[2] * p_n2;
#if defined(FIRST)
        GalerkinRHS[0] -=  0.0;
        GalerkinRHS[1] -=  0.0;

        GalerkinRHS[2] -=  0.0;
        GalerkinRHS[3] -=  0.0;

        GalerkinRHS[4] -=  0.0;
        GalerkinRHS[5] -=  0.0;
#else
        GalerkinRHS[0] +=  msN[0] * p_grad * Area ;  //ojo con el signo
        GalerkinRHS[1] +=  msN[0] * p_grad1 * Area;

        GalerkinRHS[2] +=  msN[1] * p_grad * Area;
        GalerkinRHS[3] +=  msN[1] * p_grad1 * Area;

        GalerkinRHS[4] +=  msN[2] * p_grad * Area;
        GalerkinRHS[5] +=  msN[2] * p_grad1 * Area;

#endif

        GetGeometry()[0].SetLock();
        array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(FORCE);
        rhs0[0] += GalerkinRHS[0] ;
        rhs0[1] += GalerkinRHS[1] ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(FORCE);
        rhs1[0] += GalerkinRHS[2] ;
        rhs1[1] += GalerkinRHS[3] ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(FORCE);
        rhs2[0] += GalerkinRHS[4] ;
        rhs2[1] += GalerkinRHS[5] ;
        GetGeometry()[2].UnSetLock();

        double nodal_contrib = 0.333333333333333333333333333 * Area;

	//KRATOS_WATCH(nodal_contrib);
	//KRATOS_WATCH(density);
        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[2].UnSetLock();
      }
    else if (FractionalStepNumber == 6) //calculation of velocities
      {
	///////////////////////NECESSARY LOCALS///////////////////////////////////////////
        BoundedMatrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
        array_1d<double,6> GalerkinRHS = ZeroVector(6); //dimension = number of nodes
        BoundedMatrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
        array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
        BoundedMatrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
        BoundedMatrix<double,2,6> msConvOp = ZeroMatrix(2,6);
        BoundedMatrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
        array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

	//double density = this->GetValue(DENSITY);
	//double density = 1000.0;

        const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
        const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
        const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

	double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

        //density = this->GetValue(DENSITY);
        //density=1000.0;
        double nodal_contrib = 0.333333333333333333333333333 * Area;

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib * density;
        GetGeometry()[2].UnSetLock();


        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib * density ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib * density ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA) += nodal_contrib * density ;
        GetGeometry()[2].UnSetLock();
      }

    else if (FractionalStepNumber == 7)
      {
	///////////////////////NECESSARY LOCALS///////////////////////////////////////////
        BoundedMatrix<double,3,3> msWorkMatrix = ZeroMatrix(3,3);
        array_1d<double,6> GalerkinRHS = ZeroVector(3); //dimension = number of nodes
        BoundedMatrix<double,3,2> msDN_DX = ZeroMatrix(3,2);
        array_1d<double,3> msN = ZeroVector(3); //dimension = number of nodes
        BoundedMatrix<double,6,2> msShapeFunc = ZeroMatrix(6,2);
        BoundedMatrix<double,2,6> msConvOp = ZeroMatrix(2,6);
        BoundedMatrix<double,6,6> msAuxMat = ZeroMatrix(6,6);
        array_1d<double,6> msAuxVec = ZeroVector(6); //dimension = number of nodes
        array_1d<double,2> ms_adv_vel = ZeroVector(2); //dimesion coincides with space dimension
        array_1d<double,2> ms_vel_gauss = ZeroVector(2); //dimesion coincides with space dimension
	array_1d<double,3> temp_vec_np = ZeroVector(3); //dimension = number of nodes

	double dt = CurrentProcessInfo[DELTA_TIME];

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Area);

    	Matrix Mass(3,3);
    	Mass(0,0)=2.0;
    	Mass(1,1)=2.0;
    	Mass(2,2)=2.0;
    	Mass(0,1)=1.0;
    	Mass(0,2)=1.0;
    	Mass(1,0)=1.0;
    	Mass(1,2)=1.0;
    	Mass(2,0)=1.0;
    	Mass(2,1)=1.0;
    	Mass/=12.0;

        Mass *=Area;

        double nodal_contrib = 0.333333333333333333333333333 * Area;
        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib ;
        GetGeometry()[2].UnSetLock();

        //getting the velocity on the nodes and other necessary variabless
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
        double bulk_modulus = -0.3333333333333333333333*(k0 + k1 + k2 ) * dt;
        //double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

        temp_vec_np[0] = p_n0;
        temp_vec_np[1] = p_n1;
        temp_vec_np[2] = p_n2;

        noalias(GalerkinRHS) = prod(MassFactors, temp_vec_np);

	GalerkinRHS *= Area * 1.0 ;

        double Gaux;
        Gaux =  msDN_DX(0,0) * vel0[0] + msDN_DX(0,1) * vel0[1];
        Gaux += msDN_DX(1,0) * vel1[0] + msDN_DX(1,1) * vel1[1];
        Gaux += msDN_DX(2,0) * vel2[0] + msDN_DX(2,1) * vel2[1];

	constexpr double E_over_R = 28961.49;//24466.81;
    	constexpr double C = 1.19e15;

        double t1 = GetGeometry()[0].FastGetSolutionStepValue(YCH4);
        //double Arr1 = C * exp(-E_over_R/(t1));
        double t2 = GetGeometry()[1].FastGetSolutionStepValue(YCH4);
        //double Arr2 = C * exp(-E_over_R/(t2));
        double t3 = GetGeometry()[2].FastGetSolutionStepValue(YCH4);
        //double Arr3 = C * exp(-E_over_R/(t3));

        //double Aver = 0.33333333333333*(Arr1 + Arr2 + Arr3 );

	double temp=t1 + t2 + t3;
	temp *= 0.33333333333333;
	if(temp>1000.0) temp=1000.0;

        double aux_var= C * exp(-E_over_R/(temp));
        GalerkinRHS[0] += bulk_modulus * Area * Gaux * 0.33333333333333 + 0.0 * bulk_modulus *  aux_var * 0.33333333333333 * Area ;
        GalerkinRHS[1] += bulk_modulus * Area * Gaux * 0.33333333333333 + 0.0 * bulk_modulus *  aux_var * 0.33333333333333 * Area ;
        GalerkinRHS[2] += bulk_modulus * Area * Gaux * 0.33333333333333 + 0.0 * bulk_modulus *  aux_var * 0.33333333333333 * Area ;

	//double ddd=C * exp(-E_over_R/(temp));
	//KRATOS_WATCH(ddd);

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
    }
    KRATOS_CATCH("");
}

void QFluid2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;

    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber == 1) //step 1
      {
        if(rResult.size() != number_of_nodes*dim)
	  rResult.resize(number_of_nodes*dim,false);

        for (unsigned int i=0; i<number_of_nodes; i++)
	  {
            rResult[i*dim] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
            rResult[i*dim+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();

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

  void QFluid2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
  {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;

    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    if(FractionalStepNumber == 1) //step 1
      {
	if(ElementalDofList.size() != number_of_nodes*dim)
	  ElementalDofList.resize(number_of_nodes*dim);

        for (unsigned int i=0; i<number_of_nodes; i++)
	  {
            ElementalDofList[i*dim] = GetGeometry()[i].pGetDof(VELOCITY_X);
            ElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
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


  double QFluid2D::ComputeSmagorinskyViscosity(const BoundedMatrix<double, 3, 2 > & DN_DX,
					       const double& h,
					       const double& C,
					       const double nu
					       )
  {
    BoundedMatrix<double, 2, 2 > dv_dx = ZeroMatrix(2, 2);

    // Compute Symmetric Grad(u). Note that only the lower half of the matrix is filled
    for (unsigned int k = 0; k < 3; ++k)
      {
        const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(FRACT_VEL);
        for (unsigned int i = 0; i < 2; ++i)
	  {
            for (unsigned int j = 0; j < i; ++j) // Off-diagonal
	      dv_dx(i, j) += 0.5 * (DN_DX(k, j) * rNodeVel[i] + DN_DX(k, i) * rNodeVel[j]);
            dv_dx(i, i) += DN_DX(k, i) * rNodeVel[i]; // Diagonal
	  }
      }

    // Norm[ Grad(u) ]
    double NormS(0.0);
    for (unsigned int i = 0; i < 2; ++i)
      {
        for (unsigned int j = 0; j < i; ++j)
	  NormS += 2.0 * dv_dx(i, j) * dv_dx(i, j); // Using symmetry, lower half terms of the matrix are added twice
        NormS += dv_dx(i, i) * dv_dx(i, i); // Diagonal terms
      }

    NormS = sqrt(NormS);

    // Total Viscosity
    return 2.0 * C * C * h * h * NormS;
  }

  void QFluid2D::CalculateViscousMatrix( MatrixType& K, const BoundedMatrix<double,3,2>& DN_DX, const double& nu, const double& dt, const double& kk)//,const double& bulk )
{

  //double deltat = dt;
  double viscosity=nu;

  double k=kk;

  //k=0.0;
  K(0,0) = 2.0 * pow(DN_DX(0,0), 2) * viscosity + viscosity * pow(DN_DX(0,1), 2)           + (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(0,0);
  K(0,1) = DN_DX(0,1) * viscosity * DN_DX(0,0)                                    + (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(0,1);
  K(0,2) = 2.0 * DN_DX(0,0) * viscosity * DN_DX(1,0) + DN_DX(0,1) * viscosity * DN_DX(1,1) + (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(1,0);
  K(0,3) = DN_DX(0,1) * viscosity * DN_DX(1,0)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(1,1);
  K(0,4) = 2.0 * DN_DX(0,0) * viscosity * DN_DX(2,0) + DN_DX(0,1) * viscosity * DN_DX(2,1)+ (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(2,0);
  K(0,5) = DN_DX(0,1) * viscosity * DN_DX(2,0)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(0,0) *  DN_DX(2,1);
  K(1,0) = DN_DX(0,1) * viscosity * DN_DX(0,0)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(0,0);
  K(1,1) = 2.0 * viscosity * pow(DN_DX(0,1), 2) + pow(DN_DX(0,0), 2) * viscosity         + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(0,1);
  K(1,2) = DN_DX(0,0) * viscosity * DN_DX(1,1)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(1,0);
  K(1,3) = 2.0 * DN_DX(0,1) * viscosity * DN_DX(1,1) + DN_DX(0,0) * viscosity * DN_DX(1,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(1,1);
  K(1,4) = DN_DX(0,0) * viscosity * DN_DX(2,1)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(2,0);
  K(1,5) = 2.0 * DN_DX(0,1) * viscosity * DN_DX(2,1) + DN_DX(0,0) * viscosity * DN_DX(2,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(0,1) *  DN_DX(2,1);
  K(2,0) = 2.0 * DN_DX(0,0) * viscosity * DN_DX(1,0) + DN_DX(0,1) * viscosity * DN_DX(1,1) + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(0,0);
  K(2,1) = DN_DX(0,0) * viscosity * DN_DX(1,1)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(0,1);
  K(2,2) = 2.0 * pow(DN_DX(1,0), 2) * viscosity + viscosity * pow(DN_DX(1,1), 2)           + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(1,0);
  K(2,3) = DN_DX(1,1) * viscosity * DN_DX(1,0)                                    + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(1,1);
  K(2,4) = 2.0 * DN_DX(1,0) * viscosity * DN_DX(2,0) + DN_DX(1,1) * viscosity * DN_DX(2,1) + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(2,0);
  K(2,5) = DN_DX(1,1) * viscosity * DN_DX(2,0)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(1,0) *  DN_DX(2,1);
  K(3,0) = DN_DX(0,1) * viscosity * DN_DX(1,0)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(0,0);
  K(3,1) = 2.0 * DN_DX(0,1) * viscosity * DN_DX(1,1) + DN_DX(0,0) * viscosity * DN_DX(1,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(0,1);
  K(3,2) = DN_DX(1,1) * viscosity * DN_DX(1,0)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(1,0);
  K(3,3) = 2.0 * viscosity * pow(DN_DX(1,1), 2) + pow(DN_DX(1,0), 2) * viscosity          + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(1,1);
  K(3,4) = DN_DX(1,0) * viscosity * DN_DX(2,1)                                   + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(2,0);
  K(3,5) = 2.0 * DN_DX(1,1) * viscosity * DN_DX(2,1) + DN_DX(1,0) * viscosity * DN_DX(2,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(1,1) *  DN_DX(2,1);
  K(4,0) = 2.0 * DN_DX(0,0) * viscosity * DN_DX(2,0) + DN_DX(0,1) * viscosity * DN_DX(2,1) + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(0,0);
  K(4,1) = DN_DX(0,0) * viscosity * DN_DX(2,1)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(0,1);
  K(4,2) = 2.0 * DN_DX(1,0) * viscosity * DN_DX(2,0) + DN_DX(1,1) * viscosity * DN_DX(2,1) + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(1,0);
  K(4,3) = DN_DX(1,0) * viscosity * DN_DX(2,1)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(1,1);
  K(4,4) = 2.0 * pow(DN_DX(2,0), 2) * viscosity + viscosity * pow(DN_DX(2,1), 2)           + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(2,0);
  K(4,5) = DN_DX(2,1) * viscosity * DN_DX(2,0)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(2,0) *  DN_DX(2,1);
  K(5,0) = DN_DX(0,1) * viscosity * DN_DX(2,0)                                      + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(0,0);
  K(5,1) = 2.0 * DN_DX(0,1) * viscosity * DN_DX(2,1) + DN_DX(0,0) * viscosity * DN_DX(2,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(0,1);
  K(5,2) = DN_DX(1,1) * viscosity * DN_DX(2,0)                                    + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(1,0);
  K(5,3) = 2.0 * DN_DX(1,1) * viscosity * DN_DX(2,1) + DN_DX(1,0) * viscosity * DN_DX(2,0) + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(1,1);
  K(5,4) = DN_DX(2,1) * viscosity * DN_DX(2,0)                                     + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(2,0);
  K(5,5) = 2.0 * viscosity * pow(DN_DX(2,1), 2) + pow(DN_DX(2,0), 2) * viscosity         + (k-2.0 / 3.0 * viscosity) * DN_DX(2,1) *  DN_DX(2,1);


  //filling the symmetric part
  /*for(unsigned int i = 1; i<K.size1(); i++)
    for(unsigned int j = 0; j<i; j++)
    K(i,j) = K(j,i);*/
}

  //************************************************************************************

  int QFluid2D::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      return 0;


    KRATOS_CATCH("");
  }


} // Namespace Kratos
