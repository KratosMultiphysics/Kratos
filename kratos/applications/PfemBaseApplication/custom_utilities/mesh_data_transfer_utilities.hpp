//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MESH_DATA_TRANSFER_UTILITIES_H_INCLUDED)
#define  KRATOS_MESH_DATA_TRANSFER_UTILITITES_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
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
  class MeshDataTransferUtilities
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
    KRATOS_DEFINE_LOCAL_FLAG( MASTER_ELEMENT_TO_NODE );
    KRATOS_DEFINE_LOCAL_FLAG( INITIALIZATION );

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


      void Initialize ()
      {
	VariablesSetFlag = false;
      };

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
		      Flags Options,
		      ModelPart::IndexType MeshId=0);
	                  
      
    //*******************************************************************************************
    //*******************************************************************************************

    void InitializeBoundaryData(Condition::Pointer rCondition, 
				const TransferParameters& rTransferVariables);


    //*******************************************************************************************
    //*******************************************************************************************

    void TransferBoundaryData(Condition::Pointer rCurrentCondition, 
			      Condition::Pointer rReferenceCondition, 
			      const TransferParameters& rTransferVariables);



    //*******************************************************************************************
    //*******************************************************************************************

    void TransferBoundaryData(const TransferParameters& rTransferVariables,
			      ModelPart& rModelPart,			      
			      ModelPart::IndexType MeshId=0);
	
    //*******************************************************************************************
    //*******************************************************************************************

	
    void TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
				       ModelPart& rModelPart,
				       ModelPart::IndexType MeshId=0);
	
      

    //*******************************************************************************************
    //*******************************************************************************************

	
    void TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
				       const Variable<double>& rCriticalVariable,
				       const double& CriticalValue,
				       ModelPart& rModelPart,
				       ModelPart::IndexType MeshId=0);
	


    //*******************************************************************************************
    //*******************************************************************************************
    void TransferElementalValuesToNodes( const TransferParameters& rTransferVariables,
					 ModelPart& rModelPart,
					 ModelPart::IndexType MeshId=0);
	


    //*******************************************************************************************
    //*******************************************************************************************
    void TransferNodalValuesToElements(ModelPart& rModelPart,
				       const Element & rReferenceElement,
				       PointPointerVector &list_of_new_centers,
				       std::vector<Geometry<Node<3> > >& list_of_new_vertices,
				       ModelPart::IndexType MeshId=0);
	

    //*******************************************************************************************
    //*******************************************************************************************
    void TransferElementalValuesToNodes(ModelPart& rModelPart,
					const Element & rReferenceElement,
					PointPointerVector &list_of_new_centers,
					std::vector<Geometry<Node<3> > >& list_of_new_vertices,
					ModelPart::IndexType MeshId=0);
	

    //*******************************************************************************************
    //*******************************************************************************************

    void TransferElementalValuesToElements(ModelPart& rModelPart,
					   const Element & rReferenceElement,
					   PointPointerVector &list_of_new_centers,
					   std::vector<Geometry<Node<3> > >& list_of_new_vertices,				       
					   ModelPart::IndexType MeshId=0);
	



    //*******************************************************************************************
    //*******************************************************************************************

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
	

    inline double CalculateVol(const double x0, const double y0,
			       const double x1, const double y1,
			       const double x2, const double y2)
    {
      return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }

    inline bool CalculatePosition(const double x0, const double y0,
				  const double x1, const double y1,
				  const double x2, const double y2,
				  const double xc, const double yc,
				  array_1d<double,3>& N)
    {
      double area = CalculateVol(x0,y0,x1,y1,x2,y2);

      if(area < 1e-20)
	{
	  KRATOS_THROW_ERROR( std::logic_error,"element with zero area found", "" )
	    }

      N[0] = CalculateVol(x1,y1,x2,y2,xc,yc)  / area;
      N[1] = CalculateVol(x2,y2,x0,y0,xc,yc)  / area;
      N[2] = CalculateVol(x0,y0,x1,y1,xc,yc)  / area;

      double tol = 1e-4;
      double upper_limit = 1.0+tol;
      double lower_limit = -tol;

      if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
	return true;

      return false;
    }

    //*******************************************************************************************
    //*******************************************************************************************
    void Interpolate( Geometry<Node<3> >& geom,
		      const array_1d<double,3>& N,
		      unsigned int step_data_size,
		      Node<3>::Pointer pnode);
	

    VariablesListDataValueContainer InterpolateVariables( Triangle2D3<Node<3> >& geom,
							  const array_1d<double,3>& N,
							  unsigned int step_data_size,
							  Node<3>::Pointer pnode);


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

#endif // KRATOS_MESH_DATA_TRANSFER_UTILITIES_H_INCLUDED  defined 


