//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_UTILITIES_H_INCLUDED )
#define  KRATOS_RIGID_BODY_UTILITIES_H_INCLUDED

#ifdef FIND_MAX
#undef FIND_MAX
#endif

#define FIND_MAX(a, b) ((a)>(b)?(a):(b))

// System includes
#include <cmath>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

// External includes
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.

  /** Computes Volume, Weight and Inertia Tensor of a given Mesh
   *  Mesh elements must have Properties with the variable THICKNESS in 2D in order to compute volume
   *  Mesh elements must have Properties with the variable DENSITY in order to compute mass
   *  Mesh nodes    must have a the nodal variable VOLUME_ACCELERATION in order to compute weight
   */

  class KRATOS_API(CONTACT_MECHANICS_APPLICATION) RigidBodyUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::MeshType::GeometryType                          GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidBodyUtilities(){ mEchoLevel = 0;  mParallel = true; };

    RigidBodyUtilities(bool Parallel){ mEchoLevel = 0;  mParallel = Parallel; };


    /// Destructor.
    ~RigidBodyUtilities(){};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //**************************************************************************
    //**************************************************************************

    Vector GetVolumeAcceleration(ModelPart& rModelPart)
    {
      KRATOS_TRY

      Vector VolumeAcceleration = ZeroVector(3);

      ModelPart::NodeType& rNode = rModelPart.Nodes().front();

      if( rNode.SolutionStepsDataHas(VOLUME_ACCELERATION) ){
	array_1d<double,3>& rVolumeAcceleration = rNode.FastGetSolutionStepValue(VOLUME_ACCELERATION);
	for(unsigned int i=0; i<3; i++)
	  VolumeAcceleration[i] = rVolumeAcceleration[i];
      }

      //alternative:

      // double num_nodes = 0;

      // for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); in++)
      // 	{

      // 	  if( in->SolutionStepsDataHas(VOLUME_ACCELERATION) ){ //temporary, will be checked once at the beginning only
      // 	    VolumeAcceleration += in->FastGetSolutionStepValue(VOLUME_ACCELERATION);
      // 	    num_nodes+=1;
      // 	  }

      // 	}

      // VolumeAcceleration/=num_nodes;

      //std::cout<<" VolumeAcceleration "<<VolumeAcceleration<<std::endl;

      return VolumeAcceleration;


      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    double GetElasticModulus(ModelPart& rModelPart)
    {
      KRATOS_TRY

      double ElasticModulus = 0;

      ModelPart::ElementType& rElement = rModelPart.Elements().front();

      ElasticModulus = rElement.GetProperties()[YOUNG_MODULUS];

      return ElasticModulus;


      KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    double VolumeCalculation(ModelPart& rModelPart)
    {
      KRATOS_TRY

      double ModelVolume = 0;

      bool mAxisymmetric = false;

      if( mAxisymmetric ){

	//std::cout<<"  AXISYMMETRIC MODEL "<<std::endl;

	double radius = 0;
	double two_pi = 6.28318530717958647693;


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rModelPart.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelVolumePartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension > 2 )
		std::cout<<" Axisymmetric problem with dimension: "<<dimension<<std::endl;

	      radius = 0;
	      for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		radius += ie->GetGeometry()[i].X();

	      radius/=double(ie->GetGeometry().size());

	      ModelVolumePartition[k] += ie->GetGeometry().Area() * two_pi * radius ;


	    }

	}

	for(int i=0; i<number_of_threads; i++)
	  ModelVolume += ModelVolumePartition[i];

      }
      else{


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rModelPart.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelVolumePartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {
	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();
	      if( dimension == 2){
		ModelVolumePartition[k] += ie->GetGeometry().Area() *  ie->GetProperties()[THICKNESS];
	      }
	      else if( dimension == 3){
		ModelVolumePartition[k] += ie->GetGeometry().Volume();
	      }
	      else{
		//do nothing.
	      }

	    }

	}

	for(int i=0; i<number_of_threads; i++)
	  ModelVolume += ModelVolumePartition[i];


      }

      return ModelVolume;


      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    double MassCalculation(ModelPart& rModelPart)
    {
      KRATOS_TRY

      double ModelMass = 0;

      bool mAxisymmetric = false;

      if( mAxisymmetric ){

	//std::cout<<"  AXISYMMETRIC MODEL "<<std::endl;

	double radius = 0;
	double two_pi = 6.28318530717958647693;

	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rModelPart.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension > 2 )
		std::cout<<" Axisymmetric problem with dimension: "<<dimension<<std::endl;

	      radius = 0;
	      for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		radius += ie->GetGeometry()[i].X();

	      radius/=double(ie->GetGeometry().size());

	      ModelMassPartition[k] += ie->GetGeometry().Area() * two_pi * radius * ie->GetProperties()[DENSITY];


	    }
	}

	for(int i=0; i<number_of_threads; i++)
	  ModelMass += ModelMassPartition[i];

      }
      else{


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rModelPart.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {
	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();
	      if( dimension == 2){
		ModelMassPartition[k] += ie->GetGeometry().Area() *  ie->GetProperties()[THICKNESS] * ie->GetProperties()[DENSITY];
	      }
	      else if( dimension == 3){
		ModelMassPartition[k] += ie->GetGeometry().Volume() * ie->GetProperties()[DENSITY];
	      }
	      else{
		//do nothing.
	      }

	    }
	}

	for(int i=0; i<number_of_threads; i++)
	  ModelMass += ModelMassPartition[i];


      }

      return ModelMass;


      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    Vector CalculateCenterOfMass(ModelPart& rModelPart)
    {
      KRATOS_TRY


      ModelPart::MeshType& rMesh = rModelPart.GetMesh();

      return this->CenterOfMassCalculation(rMesh);


      KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************

    Vector CenterOfMassCalculation(ModelPart::MeshType& rMesh)
    {

      KRATOS_TRY

      double ModelMass = 0;
      Vector CenterOfMass = ZeroVector(3);

      double Thickness = 1;

      bool mAxisymmetric = false;

      if( mAxisymmetric ){

	//std::cout<<"  AXISYMMETRIC MODEL "<<std::endl;

	double radius = 0;
	double two_pi = 6.28318530717958647693;


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rMesh.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension > 2 )
		std::cout<<" Axisymmetric problem with dimension: "<<dimension<<std::endl;

	      radius = 0;
	      for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		radius += ie->GetGeometry()[i].X();

	      radius/=double(ie->GetGeometry().size());

	      ModelMassPartition[k] += ie->GetGeometry().Area() * two_pi * radius * ie->GetProperties()[DENSITY];
	    }
	}

	for(int i=0; i<number_of_threads; i++)
	  ModelMass += ModelMassPartition[i];

      }
      else{


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rMesh.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);

	std::vector<Vector> CenterOfMassPartition(number_of_threads);
	for(int i=0; i<number_of_threads; i++)
	  CenterOfMassPartition[i] = ZeroVector(3);


        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension == 2){
		Thickness = ie->GetProperties()[THICKNESS];
		ModelMassPartition[k] += ie->GetGeometry().Area() * Thickness * ie->GetProperties()[DENSITY];

	      }
	      else if( dimension == 3){
		ModelMassPartition[k] += ie->GetGeometry().Volume() * ie->GetProperties()[DENSITY];
	      }
	      else{
		//do nothing.
	      }


	      const unsigned int number_of_nodes  = ie->GetGeometry().PointsNumber();

	      GeometryType::IntegrationMethod IntegrationMethod = ie->GetGeometry().GetDefaultIntegrationMethod();

	      GeometryType::JacobiansType J;
	      J = ie->GetGeometry().Jacobian( J, IntegrationMethod );

	      Matrix mInvJ;
	      double mDetJ = 0;

	      const Matrix& Ncontainer = ie->GetGeometry().ShapeFunctionsValues( IntegrationMethod );

	      //reading integration points
	      const GeometryType::IntegrationPointsArrayType& integration_points = ie->GetGeometry().IntegrationPoints( IntegrationMethod );

	      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{

		  MathUtils<double>::InvertMatrix( J[PointNumber], mInvJ, mDetJ );

		  double IntegrationWeight = mDetJ * integration_points[PointNumber].Weight() * Thickness;

		  Vector N = row( Ncontainer, PointNumber);

		  for ( unsigned int i = 0; i < number_of_nodes; i++ )
		    {

		      array_1d<double, 3> &CurrentPosition  = ie->GetGeometry()[i].Coordinates();

		      for ( unsigned int j = 0; j < dimension; j++ )
			{
			  CenterOfMassPartition[k][j] += IntegrationWeight * N[i] * ie->GetProperties()[DENSITY] * CurrentPosition[j];
			}
		    }
		}

	    }

	}

	for(int i=0; i<number_of_threads; i++){
	  ModelMass += ModelMassPartition[i];
	  CenterOfMass += CenterOfMassPartition[i];
	}

      }

      if(ModelMass!=0)
	CenterOfMass /= ModelMass;
      else
	std::cout<<"   [ WARNING: RIGID BODY WITH NO MASS ] "<<std::endl;

      return CenterOfMass;


      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************


    Matrix CalculateInertiaTensor(ModelPart& rModelPart)
    {
      KRATOS_TRY

      ModelPart::MeshType& rMesh = rModelPart.GetMesh();

      return this->InertiaTensorCalculation(rMesh);

      KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************


    Matrix InertiaTensorCalculation(ModelPart::MeshType& rMesh)
    {
      KRATOS_TRY

      double ModelMass = 0;
      Vector Center = CenterOfMassCalculation(rMesh);

      array_1d<double, 3> CenterOfMass;
      for ( unsigned int k = 0; k < 3; k++ )
	{
	  CenterOfMass[k] = Center[k];
	}

      Matrix InertiaTensor = ZeroMatrix(3,3);

      double Thickness = 1;

      bool mAxisymmetric = false;

      if( mAxisymmetric ){

	//std::cout<<"  AXISYMMETRIC MODEL "<<std::endl;

	double radius = 0;
	double two_pi = 6.28318530717958647693;

		//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rMesh.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);

        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension > 2 )
		std::cout<<" Axisymmetric problem with dimension: "<<dimension<<std::endl;

	      radius = 0;
	      for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
		radius += ie->GetGeometry()[i].X();

	      radius/=double(ie->GetGeometry().size());

	      ModelMassPartition[k] += ie->GetGeometry().Area() * two_pi * radius * ie->GetProperties()[DENSITY];
	    }
	}

	for(int i=0; i<number_of_threads; i++)
	  ModelMass += ModelMassPartition[i];


      }
      else{


	//Search and calculation in parallel
	ModelPart::ElementsContainerType& pElements = rMesh.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

	if( mParallel == false )
	  number_of_threads = 1;

	vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

	Vector ModelMassPartition = ZeroVector(number_of_threads);
	std::vector<Matrix> InertiaTensorPartition(number_of_threads);
	for(int i=0; i<number_of_threads; i++)
	  InertiaTensorPartition[i] = ZeroMatrix(3,3);


        #pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];


	  for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
	    {

	      const unsigned int dimension        = ie->GetGeometry().WorkingSpaceDimension();

	      if( dimension == 2){
		Thickness = ie->GetProperties()[THICKNESS];
		ModelMassPartition[k] += ie->GetGeometry().Area() * Thickness * ie->GetProperties()[DENSITY];

	      }
	      else if( dimension == 3){
		ModelMassPartition[k] += ie->GetGeometry().Volume() * ie->GetProperties()[DENSITY];
	      }
	      else{
		//do nothing.
	      }

	      const unsigned int number_of_nodes  = ie->GetGeometry().PointsNumber();

	      GeometryType::IntegrationMethod IntegrationMethod = ie->GetGeometry().GetDefaultIntegrationMethod();

	      GeometryType::JacobiansType J;
	      J = ie->GetGeometry().Jacobian( J, IntegrationMethod );

	      Matrix mInvJ;
	      double mDetJ = 0;

	      const Matrix& Ncontainer = ie->GetGeometry().ShapeFunctionsValues( IntegrationMethod );

	      //reading integration points
	      const GeometryType::IntegrationPointsArrayType& integration_points = ie->GetGeometry().IntegrationPoints( IntegrationMethod );

	      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{

		  MathUtils<double>::InvertMatrix( J[PointNumber], mInvJ, mDetJ );

		  double IntegrationWeight = mDetJ * integration_points[PointNumber].Weight() * Thickness;

		  Vector N = row( Ncontainer, PointNumber);

		  double MassAB = 0;
		  Matrix OuterAB = ZeroMatrix(3,3);
		  double InnerAB = 0;
		  Matrix Identity = IdentityMatrix(3);

		  for ( unsigned int i = 0; i < number_of_nodes; i++ )
		    {
		      array_1d<double, 3> CurrentPositionA  = ie->GetGeometry()[i].Coordinates();

		      CurrentPositionA -= CenterOfMass;

		      for ( unsigned int j = 0; j < number_of_nodes; j++ )
			{
			  array_1d<double, 3> CurrentPositionB  = ie->GetGeometry()[j].Coordinates();

			  CurrentPositionB -= CenterOfMass;

			  MassAB = IntegrationWeight * N[i] * N[j] * ie->GetProperties()[DENSITY];

			  OuterAB = outer_prod( CurrentPositionA, CurrentPositionB );
			  InnerAB = inner_prod( CurrentPositionA, CurrentPositionB );

			  InertiaTensorPartition[k] += MassAB * ( InnerAB * Identity - OuterAB );

			}

		    }
		}

	    }
	}

	for(int i=0; i<number_of_threads; i++){
	  ModelMass += ModelMassPartition[i];
	  InertiaTensor += InertiaTensorPartition[i];
	}
      }

      // test to clean negligible values:
      //clear_numerical_errors(InertiaTensor);

      return InertiaTensor;


      KRATOS_CATCH( "" )
    }



    //**************************************************************************
    //**************************************************************************


    Matrix InertiaTensorMainAxesCalculation(ModelPart& rModelPart,Matrix& MainAxes)
    {
      KRATOS_TRY

      Matrix InertiaTensor = ZeroMatrix(3,3);
      InertiaTensor = this->CalculateInertiaTensor(rModelPart);

      this->InertiaTensorToMainAxes(InertiaTensor, MainAxes);

      return InertiaTensor;


      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************


    void InertiaTensorToMainAxes(Matrix& InertiaTensor, Matrix& MainAxes)
    {
      KRATOS_TRY

      MainAxes = ZeroMatrix(3,3);
      Vector MainInertias = ZeroVector(3);
      EigenVectors3x3(InertiaTensor, MainAxes ,MainInertias);

      InertiaTensor = ZeroMatrix(3,3);
      for(unsigned int i=0; i<3; i++)
	InertiaTensor(i,i) = MainInertias[i];


      KRATOS_CATCH( "" )
    }



    void TransferRigidBodySkinToModelPart(ModelPart& rSourceModelPart, ModelPart& rDestinationModelPart, const int mesh_id){

        ElementsContainerType&  mesh_elements = rSourceModelPart.GetMesh(mesh_id).Elements();

        for(unsigned int i=0; i<mesh_elements.size(); i++) {
            ElementsContainerType::ptr_iterator pTubeElement = mesh_elements.ptr_begin()+i;
            (rDestinationModelPart.Elements()).push_back(*pTubeElement);
        }

        NodesContainerType&  mesh_nodes = rSourceModelPart.GetMesh(mesh_id).Nodes();

        for(unsigned int i=0; i<mesh_nodes.size(); i++) {
            NodesContainerType::ptr_iterator pTubeNode = mesh_nodes.ptr_begin()+i;
            (rDestinationModelPart.Nodes()).push_back(*pTubeNode);
        }

    }


    //**************************************************************************
    //**************************************************************************



    /**
     * level of echo for the error calculation
     */
    virtual void SetEchoLevel(int Level)
    {
      mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
      return mEchoLevel;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}


  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    int mEchoLevel;

    bool mParallel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    inline void clear_numerical_errors(Matrix& rMatrix)
    {
      double Maximum = 0;

      //get max value
      for(unsigned int i=0; i<rMatrix.size1(); i++)
	{
	for(unsigned int j=0; j<rMatrix.size2(); j++)
	  {
	    if(Maximum<fabs(rMatrix(i,j)))
	      Maximum = fabs(rMatrix(i,j));
	  }
	}

      //set zero value if value is < Maximum *1e-5
      double Critical = Maximum*1e-5;

      for(unsigned int i=0; i<rMatrix.size1(); i++)
	{
	for(unsigned int j=0; j<rMatrix.size2(); j++)
	  {
	    if(fabs(rMatrix(i,j))<Critical)
	      rMatrix(i,j) = 0;
	  }
	}

    }


    /**
     * calculates Eigenvectors and the EigenValues of given square matrix A(3x3)
     * The QR Algorithm with shifts is used
     * @param A the given symmetric 3x3 real matrix
     * @param V matrix to store the eigenvectors in rows
     * @param d matrix to store the eigenvalues
     */

    static inline void EigenVectors3x3(const Matrix& A, Matrix&V, Vector& d)
    {

      if(A.size1()!=3 || A.size2()!=3)
	std::cout<<" GIVEN MATRIX IS NOT 3x3  eigenvectors calculation "<<std::endl;

      Vector e = ZeroVector(3);
      V = ZeroMatrix(3,3);

      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  V(i,j) = A(i,j);
	}
      }

      // *******************//
      //Symmetric Housholder reduction to tridiagonal form

      for (int j = 0; j < 3; j++) {
	d[j] = V(2,j);
      }

      // Householder reduction to tridiagonal form.

      for (int i = 2; i > 0; i--) {

	// Scale to avoid under/overflow.

	double scale = 0.0;
	double h = 0.0;
	for (int k = 0; k < i; k++) {
	  scale = scale + fabs(d[k]);
	}
	if (scale == 0.0) {
	  e[i] = d[i-1];
	  for (int j = 0; j < i; j++) {
	    d[j] = V(i-1,j);
	    V(i,j) = 0.0;
	    V(j,i) = 0.0;
	  }
	} else {

	  // Generate Householder vector.

	  for (int k = 0; k < i; k++) {
	    d[k] /= scale;
	    h += d[k] * d[k];
	  }
	  double f = d[i-1];
	  double g = sqrt(h);
	  if (f > 0) {
	    g = -g;
	  }
	  e[i] = scale * g;
	  h = h - f * g;
	  d[i-1] = f - g;
	  for (int j = 0; j < i; j++) {
	    e[j] = 0.0;
	  }

	  // Apply similarity transformation to remaining columns.

	  for (int j = 0; j < i; j++) {
	    f = d[j];
	    V(j,i) = f;
	    g = e[j] + V(j,j) * f;
	    for (int k = j+1; k <= i-1; k++) {
	      g += V(k,j) * d[k];
	      e[k] += V(k,j) * f;
	    }
	    e[j] = g;
	  }
	  f = 0.0;
	  for (int j = 0; j < i; j++) {
	    e[j] /= h;
	    f += e[j] * d[j];
	  }
	  double hh = f / (h + h);
	  for (int j = 0; j < i; j++) {
	    e[j] -= hh * d[j];
	  }
	  for (int j = 0; j < i; j++) {
	    f = d[j];
	    g = e[j];
	    for (int k = j; k <= i-1; k++) {
	      V(k,j) -= (f * e[k] + g * d[k]);
	    }
	    d[j] = V(i-1,j);
	    V(i,j) = 0.0;
	  }
	}
	d[i] = h;
      }

      // Accumulate transformations.

      for (int i = 0; i < 2; i++) {
	V(2,i) = V(i,i);
	V(i,i) = 1.0;
	double h = d[i+1];
	if (h != 0.0) {
	  for (int k = 0; k <= i; k++) {
	    d[k] = V(k,i+1) / h;
	  }
	  for (int j = 0; j <= i; j++) {
	    double g = 0.0;
	    for (int k = 0; k <= i; k++) {
	      g += V(k,i+1) * V(k,j);
	    }
	    for (int k = 0; k <= i; k++) {
	      V(k,j) -= g * d[k];
	    }
	  }
	}
	for (int k = 0; k <= i; k++) {
	  V(k,i+1) = 0.0;
	}
      }
      for (int j = 0; j < 3; j++) {
	d[j] = V(2,j);
	V(2,j) = 0.0;
      }
      V(2,2) = 1.0;
      e[0] = 0.0;

      // *******************//

      // Symmetric tridiagonal QL algorithm.

      for (int i = 1; i < 3; i++) {
	e[i-1] = e[i];
      }
      e[2] = 0.0;

      double f = 0.0;
      double tst1 = 0.0;
      double eps = std::pow(2.0,-52.0);
      for (int l = 0; l < 3; l++) {

	// Find small subdiagonal element

	tst1 = FIND_MAX(tst1,fabs(d[l]) + fabs(e[l]));
	int m = l;
	while (m < 3) {
	  if (fabs(e[m]) <= eps*tst1) {
	    break;
	  }
	  m++;
	}

	// If m == l, d[l] is an eigenvalue,
	// otherwise, iterate.

	if (m > l) {
	  int iter = 0;
	  do {
	    iter = iter + 1;  // (Could check iteration count here.)

	    // Compute implicit shift

	    double g = d[l];
	    double p = (d[l+1] - g) / (2.0 * e[l]);
	    double r = hypot2(p,1.0);
	    if (p < 0) {
	      r = -r;
	    }
	    d[l] = e[l] / (p + r);
	    d[l+1] = e[l] * (p + r);
	    double dl1 = d[l+1];
	    double h = g - d[l];
	    for (int i = l+2; i < 3; i++) {
	      d[i] -= h;
	    }
	    f = f + h;

	    // Implicit QL transformation.

	    p = d[m];
	    double c = 1.0;
	    double c2 = c;
	    double c3 = c;
	    double el1 = e[l+1];
	    double s = 0.0;
	    double s2 = 0.0;
	    for (int i = m-1; i >= l; i--) {
	      c3 = c2;
	      c2 = c;
	      s2 = s;
	      g = c * e[i];
	      h = c * p;
	      r = hypot2(p,e[i]);
	      e[i+1] = s * r;
	      s = e[i] / r;
	      c = p / r;
	      p = c * d[i] - s * g;
	      d[i+1] = h + s * (c * g + s * d[i]);

	      // Accumulate transformation.

	      for (int k = 0; k < 3; k++) {
		h = V(k,i+1);
		V(k,i+1) = s * V(k,i) + c * h;
		V(k,i) = c * V(k,i) - s * h;
	      }
	    }
	    p = -s * s2 * c3 * el1 * e[l] / dl1;
	    e[l] = s * p;
	    d[l] = c * p;

	    // Check for convergence.

	  } while (fabs(e[l]) > eps*tst1);
	}
	d[l] = d[l] + f;
	e[l] = 0.0;
      }

      // Sort eigenvalues and corresponding vectors.

      for (int i = 0; i < 2; i++) {
	int k = i;
	double p = d[i];
	for (int j = i+1; j < 3; j++) {
	  if (d[j] < p) {
	    k = j;
	    p = d[j];
	  }
	}
	if (k != i) {
	  d[k] = d[i];
	  d[i] = p;
	  for (int j = 0; j < 3; j++) {
	    p = V(j,i);
	    V(j,i) = V(j,k);
	    V(j,k) = p;
	  }
	}
      }

      V = trans(V);

      // *******************//


    }


    static double hypot2(double x, double y) {
      return std::sqrt(x*x+y*y);
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}


  }; // Class RigidBodyUtilities

  ///@}
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

} // namespace Kratos.

#endif // KRATOS_RIGID_BODY_UTILITIES_H_INCLUDED defined
