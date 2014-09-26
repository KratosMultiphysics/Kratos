//
//   Project Name:                       MachiningApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                     May 2014 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MESH_MODELER_H_INCLUDED )
#define  KRATOS_MESH_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "processes/find_nodal_h_process.h"

#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/boundary_skin_build_process.hpp"

#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/change_tip_elements_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"


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
class MeshModeler
{
protected:
    

    struct RefineVariables
    {      
    private:

      //Pointer variables
      const Variable< double >* mpDissipationVariable; 
      const Variable< double >* mpErrorVariable; 

    public: 
      int      NumberOfElements;     
      double   SizeFactor;          //nodal h  size factor
      double   CriticalRadius;      //critical area   size
      double   CriticalSide;        //critical length size

      
      double   CriticalDissipation; //critical dissipative variable value
      double   ReferenceError;      //critical error percentage
     

      void SetDissipationVariable (const Variable<double>& rVariable)
      {
	mpDissipationVariable = &rVariable;
      };
      void SetErrorVariable       (const Variable<double>& rVariable)
      {
	mpErrorVariable = &rVariable;
      };

      const Variable<double>& GetDissipationVariable()
      {
	return *mpDissipationVariable;
      };

      const Variable<double>& GetErrorVariable()
      {
	return *mpErrorVariable;
      };

      void Initialize (){
	
	NumberOfElements    = 0;
	SizeFactor          = 0;  
	CriticalRadius      = 0;  
	CriticalSide        = 0;  
	CriticalDissipation = 0;
	ReferenceError      = 0;

      };

    };


    struct InfoVariables
    {      
    public: 

      int   NumberOfElements;
      int   NumberOfNodes;
      int   NumberOfConditions;

      int   CriticalElements;
      
      int   InsertedNodes;
      int   RemovedNodes;
      int   InsertedBoundaryNodes;
      int   InsertedConditions;

      int   NumberOfNewElements;
      int   NumberOfNewNodes;
      int   NumberOfNewConditions;

      bool  GeometricalSmoothingRequired;
      bool  MechanicalSmoothingRequired;

      void Initialize (){
	
	NumberOfElements   = 0;
	NumberOfNodes      = 0;
	NumberOfConditions = 0;

	CriticalElements = 0;
      
	InsertedNodes = 0;
	RemovedNodes  = 0;
	
	InsertedBoundaryNodes = 0;
	InsertedConditions = 0;
	
	NumberOfNewElements = 0;
	NumberOfNewNodes = 0;
	NumberOfNewConditions = 0;

	GeometricalSmoothingRequired = false;
	MechanicalSmoothingRequired  = false;

      };
      
      void CheckGeometricalSmoothing(){

	if( InsertedNodes * 100 > NumberOfNodes || RemovedNodes * 100 > NumberOfNodes ){
	  GeometricalSmoothingRequired = true;
	}
	else if( (InsertedNodes + RemovedNodes) * 200> NumberOfNodes ){
	  GeometricalSmoothingRequired = true;
	}
	else{
	  GeometricalSmoothingRequired = false;
	}

      }


      void CheckMechanicalSmoothing(){

	if( CriticalElements > 0 || NumberOfNewElements != 0 )
	  MechanicalSmoothingRequired = true;
	else
	  MechanicalSmoothingRequired = false;

      }



    };


    struct BoundingBoxVariables
    {
    public:
      bool    IsSetFlag;
      double  Radius;
      Vector  Center;
      Vector  Velocity;

      void Initialize (){
	IsSetFlag = false;
	Radius = 0;
	Center = ZeroVector(3);
	Velocity = ZeroVector(3);
      };

    };


    struct MeshingVariables
    {
    protected:
      //Pointer variables
      const Element   *mpReferenceElement;
      const Condition *mpReferenceCondition;

    public:

      void SetReferenceElement   (const Element   & rElement)
      {
	mpReferenceElement=&rElement;
      };
      void SetReferenceCondition (const Condition & rCondition)
      {
	mpReferenceCondition=&rCondition;
      };

      Element const&    GetReferenceElement   ()
      {
	return *mpReferenceElement;
      };
      Condition const&  GetReferenceCondition ()
      {
	return *mpReferenceCondition;
      };

      //Other variables (general struct access)
      Flags   MeshingOptions;
      Flags   RefiningOptions;
      
      //General meshing flags
      bool    RemeshFlag;
      bool    RefineFlag;
      bool    ConstrainedFlag;
      bool    MeshSmoothingFlag;
      bool    JacobiSmoothingFlag;
      bool    AvoidTipElementsFlag;
      bool    RigidWallSetFlag;

      bool    NodalIdsSetFlag;
      double  AlphaParameter;
      double  OffsetFactor;

      std::vector<int> NodalPreIds;
      std::vector<int> NodalNewIds;

      std::vector<int> PreservedElements;
      std::vector<std::vector<int> > NeighbourList; 

      RefineVariables  Refine;  
       
      InfoVariables RemeshInfo;

      std::vector<RigidWallBoundingBox::Pointer> RigidWalls;

      BoundingBoxVariables BoundingBox;
 
      void Initialize (){
  
	RemeshFlag = false;
	RefineFlag = false;
	ConstrainedFlag = false;
	MeshSmoothingFlag = false;
	JacobiSmoothingFlag = false;
	AvoidTipElementsFlag = false;
	RigidWallSetFlag = false;
	NodalIdsSetFlag = false;
	
	AlphaParameter = 0;
	OffsetFactor  = 0;

	Refine.Initialize();
	RemeshInfo.Initialize();
	BoundingBox.Initialize();
	
      };
    };

  

public:

    /**
     * Flags related to the meshing parameters
     */

    //meshing options
    KRATOS_DEFINE_LOCAL_FLAG ( REMESH );
    KRATOS_DEFINE_LOCAL_FLAG ( RECONNECT );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_MESH );

    KRATOS_DEFINE_LOCAL_FLAG ( CONSTRAINED_MESH );
    KRATOS_DEFINE_LOCAL_FLAG ( BOUNDARIES_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( NEIGHBOURS_SEARCH );

    KRATOS_DEFINE_LOCAL_FLAG ( SET_DOF );

    KRATOS_DEFINE_LOCAL_FLAG ( CONTACT_SEARCH );

    //refining options
    KRATOS_DEFINE_LOCAL_FLAG ( SELECT_ELEMENTS );
    KRATOS_DEFINE_LOCAL_FLAG ( PASS_ALPHA_SHAPE );

    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_INSERT_NODES );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ADD_NODES );

    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ELEMENTS );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_BOUNDARY );

    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_NODES );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_ON_BOUNDARY );

    KRATOS_DEFINE_LOCAL_FLAG ( CRITERION_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG ( CRITERION_ENERGY );
    KRATOS_DEFINE_LOCAL_FLAG ( CRITERION_DISTANCE );

    KRATOS_DEFINE_LOCAL_FLAG ( ENGAGED_NODES );

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshModeler
    KRATOS_CLASS_POINTER_DEFINITION( MeshModeler );

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    //typedef ModelPart::GeometricalDataContainerType GeometricalDataContainerType;

    //typedef GeometricalDataContainerType::GeometricalDataType GeometricalDataType;

    //typedef ModelPart::GeometryType GeometryType;

    //typedef ModelPart::GeometriesContainerType GeometriesContainerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshModeler() {}

    /// Destructor.
    virtual ~MeshModeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //************* STARTING METHODS

    /**
     * Called to initialize the modeler
     * Must be called before any calculation is done
     */
    void Initialize (int NumberOfDomains);
  
    
    /**
     * level of echo for the mesh modeler
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
	mDataTransferUtilities.SetEchoLevel(Level);
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }

    /**
     * Remesh information is given to the modeler
     */
    void SetRemeshData (Element   const& rReferenceElement,
			Condition const& rReferenceCondition,
			bool RemeshFlag           = false,
			bool ConstrainedFlag      = false,
			bool MeshSmoothingFlag    = false,
			bool JacobiSmoothingFlag  = false,
			bool AvoidTipElementsFlag = false,
			double AlphaParameter     = 1.4,
			double OffsetFactor       = 1.0,
			int    MeshId             = 0);

    /**
     * Refine information is given to the modeler
     */
    void SetRefineData (bool   RefineFlag  = false,
			double SizeFactor  = 0.5,
			double Dissipation = 40,
			double Radius      = 0.00004,
			double Error       = 2,
			int    MeshId      = 0);

    /**
     * Walls of the domain are given to the modeler (refining properties)
     */
    void SetRigidWall (RigidWallBoundingBox::Pointer pRigidWall);


    /**
     * Boxes to refine are given to the modeler
     */
    void SetRefiningBox (double Radius,
			 Vector Center,
			 Vector Velocity);


    void SetRefiningBox (double Radius,
			 Vector Center,
			 Vector Velocity,
			 int MeshId);


    //*******************************************************************************************
    //*******************************************************************************************

    /**
     * Mesh Modeler :: Remesh all ModelPart
     */
    virtual void GenerateMesh(ModelPart& rModelPart);
   

      
    //*******************************************************************************************
    //*******************************************************************************************

    virtual void GenerateMesh(ModelPart& ThisModelPart, Element const& rReferenceElement, Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR( std::logic_error, "This modeler CAN NOT be used for mesh generation.", "" )
    }

    virtual void GenerateNodes(ModelPart& ThisModelPart)
    {
        KRATOS_ERROR( std::logic_error, "This modeler CAN NOT be used for node generation.", "" )
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MeshModeler";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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

    std::vector<MeshingVariables> mMeshingVariables;
    
    ModelerUtilities              mModelerUtilities;
    MeshDataTransferUtilities     mDataTransferUtilities;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Assignment operator.
    MeshModeler& operator=(MeshModeler const& rOther);

    /// Copy constructor.
    MeshModeler(MeshModeler const& rOther){};


    ///@}
    ///@name Protected Operations
    ///@{

    

   /**
     * Mesh Modeler :: Variables Transfer without remeshing
     */
    virtual void PerformTransferOnly(ModelPart& rModelPart,
				     MeshingVariables& rMeshingVariables,
				     ModelPart::IndexType MeshId=0){};

    /**
     * Mesh Modeler :: Delaunay Tessellation
     */
    virtual void GenerateDT (ModelPart& rModelPart,
			     MeshingVariables& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};
    
    /**
     * Mesh Modeler :: Constrained Delaunay Tessellation
     */
    virtual void GenerateCDT(ModelPart& rModelPart,
			     MeshingVariables& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};
    
    /**
     * Mesh Modeler :: Refined Delaunay Tessellation
     */
    virtual void GenerateRDT(ModelPart& rModelPart,
			     MeshingVariables& rMeshingVariables,
			     ModelPart::IndexType MeshId=0){};

    /**
     * Mesh Modeler :: Refined Constrained Delaunay Tessellation
     */
    virtual void GenerateRCDT(ModelPart& rModelPart,
			      MeshingVariables& rMeshingVariables,
			      ModelPart::IndexType MeshId=0){};
        

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
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class MeshModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MeshModeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MeshModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MESH_MODELER_H_INCLUDED  defined 


