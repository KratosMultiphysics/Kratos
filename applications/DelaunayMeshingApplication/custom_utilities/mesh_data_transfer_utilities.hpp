//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_MESH_DATA_TRANSFER_UTILITIES_H_INCLUDED)
#define KRATOS_MESH_DATA_TRANSFER_UTILITIES_H_INCLUDED

// System includes
#include <stdlib.h>

// Project includes
#include "includes/variables.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/variables_list_data_value_container.h"

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
  /** Detail class definition.
   */
  class KRATOS_API(DELAUNAY_MESHING_APPLICATION) MeshDataTransferUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of data transfer
    KRATOS_CLASS_POINTER_DEFINITION( MeshDataTransferUtilities );

    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;
    typedef std::vector<Node<3>::Pointer >                       PointPointerVector;

    /**
     * Flags related to the meshing parameters
     */

    KRATOS_DEFINE_LOCAL_FLAG( NODE_TO_ELEMENT );
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENT_TO_NODE );
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENT_TO_ELEMENT );
    KRATOS_DEFINE_LOCAL_FLAG( INITIALIZE_MASTER_CONDITION );
    KRATOS_DEFINE_LOCAL_FLAG( MASTER_ELEMENT_TO_MASTER_CONDITION );

    struct TransferParameters
    {

      KRATOS_CLASS_POINTER_DEFINITION(TransferParameters);

      Flags    Options;

      std::vector<const Variable<double>* >               DoubleVariables;
      std::vector<const Variable<array_1d<double,3> >* > Array1DVariables;
      std::vector<const Variable<Vector>* >               VectorVariables;
      std::vector<const Variable<Matrix>* >               MatrixVariables;

      bool VariablesSetFlag;

      // setting refining variables (generally for python interface)
      void Set(Flags ThisFlag)
      {
	Options.Set(ThisFlag);
      };

      void Reset(Flags ThisFlag)
      {
	Options.Reset(ThisFlag);
      };

      void SetOptions(const Flags&  rOptions)
      {
	Options=rOptions;
      };


      void SetVariable(const Variable<double>& pVariable)
      {
	VariablesSetFlag = true;
	DoubleVariables.push_back(&pVariable);
      }

      void SetVariable(const Variable<array_1d<double,3> >& pVariable)
      {
	VariablesSetFlag = true;
	Array1DVariables.push_back(&pVariable);
      }


      void SetVariable(const Variable<Vector>& pVariable)
      {
	VariablesSetFlag = true;
	VectorVariables.push_back(&pVariable);
      }


      void SetVariable(const Variable<Matrix>& pVariable)
      {
	VariablesSetFlag = true;
	MatrixVariables.push_back(&pVariable);
      }

      Flags GetOptions()
      {
	return Options;
      };

      void Initialize ()
      {
	VariablesSetFlag = false;
      };

    };


    struct BoundaryVariables
    {

      double DoubleVariable;
      array_1d<double, 3> Array1DVariable;
      Vector VectorVariable;
      Matrix MatrixVariable;

      void Initialize(const unsigned int& dimension, const unsigned int& voigt_size)
      {
	DoubleVariable = 0;
	Array1DVariable.clear();
	VectorVariable.resize(voigt_size);
	noalias(VectorVariable) = ZeroVector(voigt_size);
	MatrixVariable.resize(dimension, dimension);
	noalias(MatrixVariable) = IdentityMatrix(dimension);
      }

    };


    struct BoundaryVariableArrays
    {
      unsigned int array_size;
      std::vector<double> DoubleVariableArray;
      std::vector<array_1d<double,3> > Array1DVariableArray;
      std::vector<Vector> VectorVariableArray;
      std::vector<Matrix> MatrixVariableArray;

      void Initialize(const unsigned int& size)
      {
	array_size = size;
	DoubleVariableArray.resize(size);
	Array1DVariableArray.resize(size);
	VectorVariableArray.resize(size);
	MatrixVariableArray.resize(size);
      }

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshDataTransferUtilities()
    {
    } //


    /// Copy constructor.
    MeshDataTransferUtilities(MeshDataTransferUtilities const& rOther)
    {
    } //


    /// Destructor.
    virtual ~MeshDataTransferUtilities() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************

    void TransferData(ModelPart& rModelPart,
		      const Element & rReferenceElement,
		      PointPointerVector &list_of_new_centers,
		      std::vector<Geometry<Node<3> > >& list_of_new_vertices,
		      Flags Options);


    //*******************************************************************************************
    //*******************************************************************************************

    void InitializeBoundaryData(Condition::Pointer rCurrentCondition,
				const TransferParameters& rTransferVariables,
				const ProcessInfo& rCurrentProcessInfo);



    //*******************************************************************************************
    //*******************************************************************************************

    void TransferInitialBoundaryData(Condition::Pointer rCurrentCondition,
				     const TransferParameters& rTransferVariables,
				     BoundaryVariables& rVariables);

    //*******************************************************************************************
    //*******************************************************************************************

    void TransferCurrentBoundaryData(Element::Pointer rCurrentElement,
				     Condition::Pointer rCurrentCondition,
				     const TransferParameters& rTransferVariables,
				     BoundaryVariables& rVariables,
				     BoundaryVariableArrays& rVariableArrays,
				     const ProcessInfo& rCurrentProcessInfo);


    //*******************************************************************************************
    //*******************************************************************************************

    void TransferBoundaryData(Condition::Pointer rCurrentCondition,
			      Condition::Pointer rReferenceCondition,
			      const TransferParameters& rTransferVariables);


    //*******************************************************************************************
    //*******************************************************************************************

    void TransferBoundaryData(Element::Pointer rCurrentElement,
			      Condition::Pointer rCurrentCondition,
			      const TransferParameters& rTransferVariables,
			      const ProcessInfo& rCurrentProcessInfo);


    //*******************************************************************************************
    //*******************************************************************************************

    void TransferBoundaryData(const TransferParameters& rTransferVariables,
			      ModelPart& rModelPart);

    //*******************************************************************************************
    //*******************************************************************************************


    void TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
				       ModelPart& rModelPart);



    //*******************************************************************************************
    //*******************************************************************************************


    void TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
				       const Variable<double>& rCriticalVariable,
				       const double& CriticalValue,
				       ModelPart& rModelPart);



    //*******************************************************************************************
    //*******************************************************************************************
    void TransferElementalValuesToNodes( const TransferParameters& rTransferVariables,
					 ModelPart& rModelPart);



    //*******************************************************************************************
    //*******************************************************************************************
    void TransferNodalValuesToElements(ModelPart& rModelPart,
				       const Element & rReferenceElement,
				       PointPointerVector &list_of_new_centers,
				       std::vector<Geometry<Node<3> > >& list_of_new_vertices);


    //*******************************************************************************************
    //*******************************************************************************************
    void TransferElementalValuesToNodes(ModelPart& rModelPart,
					const Element & rReferenceElement,
					PointPointerVector &list_of_new_centers,
					std::vector<Geometry<Node<3> > >& list_of_new_vertices);


    //*******************************************************************************************
    //*******************************************************************************************

    void TransferElementalValuesToElements(ModelPart& rModelPart,
					   const Element & rReferenceElement,
					   PointPointerVector &list_of_new_centers,
					   std::vector<Geometry<Node<3> > >& list_of_new_vertices);



    //*******************************************************************************************
    //*******************************************************************************************

    inline void CalculateCenterAndSearchRadius(const std::vector<std::vector<double> >& rPointCoordinates,
					       std::vector<double>& rCenter, double& rRadius)
    {

      if( rPointCoordinates.size() == 3 ){

	CalculateCenterAndSearchRadius( rPointCoordinates[0][0], rPointCoordinates[0][1],
					rPointCoordinates[1][0], rPointCoordinates[1][1],
					rPointCoordinates[2][0], rPointCoordinates[2][1],
					rCenter[0], rCenter[1], rRadius);
      }
      else if( rPointCoordinates.size() == 4 ){

	CalculateCenterAndSearchRadius( rPointCoordinates[0][0], rPointCoordinates[0][1], rPointCoordinates[0][2],
					rPointCoordinates[1][0], rPointCoordinates[1][1], rPointCoordinates[1][2],
					rPointCoordinates[2][0], rPointCoordinates[2][1], rPointCoordinates[2][2],
					rPointCoordinates[3][0], rPointCoordinates[3][1], rPointCoordinates[3][2],
					rCenter[0], rCenter[1], rCenter[2], rRadius);
      }
      else{
	 KRATOS_THROW_ERROR( std::logic_error,"Number of points supplied out of range ERROR", "" )
      }


    }

    //*******************************************************************************************
    //*******************************************************************************************


    inline void CalculateCenterAndSearchRadius(const double x0, const double y0,
					       const double x1, const double y1,
					       double& xc, double& yc, double& R)
    {
      xc = 0.5*(x0+x1);
      yc = 0.5*(y0+y1);

      R = sqrt( (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) );
    }


    inline void CalculateCenterAndSearchRadius(const double x0, const double y0, const double z0,
					       const double x1, const double y1, const double z1,
					       const double x2, const double y2, const double z2,
					       double& xc, double& yc, double& zc, double& R)
    {
      xc = 0.3333333333333333333*(x0+x1+x2);
      yc = 0.3333333333333333333*(y0+y1+y2);
      zc = 0.3333333333333333333*(z0+z1+z2);

      double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
      double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
      double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);

      R = R1;
      if(R2 > R) R = R2;
      if(R3 > R) R = R3;

      R = sqrt(R);
    }


    inline void CalculateCenterAndSearchRadius(const double x0, const double y0,
					       const double x1, const double y1,
					       const double x2, const double y2,
					       double& xc, double& yc, double& R)
    {
      xc = 0.3333333333333333333*(x0+x1+x2);
      yc = 0.3333333333333333333*(y0+y1+y2);

      double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
      double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
      double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

      R = R1;
      if(R2 > R) R = R2;
      if(R3 > R) R = R3;

      R = sqrt(R);
    }

    inline void CalculateCenterAndSearchRadius(const double x0, const double y0, const double z0,
					       const double x1, const double y1, const double z1,
					       const double x2, const double y2, const double z2,
					       const double x3, const double y3, const double z3,
					       double& xc, double& yc, double& zc, double& R)
    {
      xc = 0.25*(x0+x1+x2+x3);
      yc = 0.25*(y0+y1+y2+y3);
      zc = 0.25*(z0+z1+z2+z3);

      double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
      double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
      double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);
      double R4 = (xc-x3)*(xc-x3) + (yc-y3)*(yc-y3) + (zc-z3)*(zc-z3);

      R = R1;
      if(R2 > R) R = R2;
      if(R3 > R) R = R3;
      if(R4 > R) R = R4;

      R = sqrt(R);
    }


    //*******************************************************************************************
    //*******************************************************************************************


    void FillVectorData( VariablesList& rVariablesList,
			 Node<3>::Pointer pnode);

    void Interpolate2Nodes( Geometry<Node<3> > &geom,
			    const std::vector<double>& N,
			    VariablesList& rVariablesList,
			    Node<3>& pnode);

    void Interpolate( Geometry<Node<3> >& geom,
		      const std::vector<double>& N,
		      VariablesList& rVariablesList,
		      Node<3>::Pointer pnode,
		      double& alpha);

    VariablesListDataValueContainer InterpolateVariables( Geometry<Node<3> >& geom,
							  const std::vector<double>& N,
							  VariablesList& rVariablesList,
							  Node<3>::Pointer pnode,
							  double& alpha);

    void InterpolateData( Geometry<Node<3> >& geom,
			  const std::vector<double>& N,
			  unsigned int step_data_size,
			  Node<3>::Pointer pnode,
			  double& alpha);


    VariablesListDataValueContainer InterpolateVariablesData( Geometry<Node<3> >& geom,
							      const std::vector<double>& N,
							      unsigned int step_data_size,
							      Node<3>::Pointer pnode,
							      double& alpha);


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
      return "";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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
    int mEchoLevel;

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


    ///@}
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    MeshDataTransferUtilities& operator=(MeshDataTransferUtilities const& rOther);

    ///@}
    ///@name Private Operations
    ///@{


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

  }; // Class MeshDataTransferUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MeshDataTransferUtilities& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshDataTransferUtilities& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif //  KRATOS_MESH_DATA_TRANSFER_UTILITITES_H_INCLUDED defined
