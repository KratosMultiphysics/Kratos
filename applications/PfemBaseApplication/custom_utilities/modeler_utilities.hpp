//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MODELER_UTILITIES_H_INCLUDED )
#define  KRATOS_MODELER_UTILITIES_H_INCLUDED

// External includes

// System includes

// Project includes
#include "processes/process.h"
#include "containers/flags.h"
#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"
//#include "includes/kratos_flags.h"

#include "pfem_base_application_variables.h"

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
class ModelerUtilities
{
public:
  

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( ModelerUtilities );


    typedef Node<3>                                                       PointType;
    typedef Node<3>::Pointer                                       PointPointerType;
    typedef Geometry< Node<3> >                                        GeometryType;
    typedef std::vector<PointPointerType>                        PointPointerVector;
    typedef ModelPart::PropertiesType                                PropertiesType;
    typedef ModelPart::PropertiesContainerType              PropertiesContainerType;
    typedef ModelPart::MeshType                                            MeshType;
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;
    typedef MeshDataTransferUtilities::TransferParameters    TransferParametersType;


    /**
     * Flags related to the meshing parameters
     */

    //meshing options

    //(configuration)
    KRATOS_DEFINE_LOCAL_FLAG ( REMESH );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE );
    KRATOS_DEFINE_LOCAL_FLAG ( RECONNECT );
    KRATOS_DEFINE_LOCAL_FLAG ( CONSTRAINED );
    KRATOS_DEFINE_LOCAL_FLAG ( CONTACT_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( MESH_SMOOTHING );
    KRATOS_DEFINE_LOCAL_FLAG ( VARIABLES_SMOOTHING );
     
    //execution options (tessellation)
    KRATOS_DEFINE_LOCAL_FLAG ( NEIGHBOURS_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( BOUNDARIES_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( SET_DOF );

    //removing options

    //(configuration)
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_NODES );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_NODES_ON_DISTANCE );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_NODES_ON_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_NODES_ON_THRESHOLD );
  
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_BOUNDARY_NODES );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_BOUNDARY_NODES_ON_DISTANCE );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_BOUNDARY_NODES_ON_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG ( REMOVE_BOUNDARY_NODES_ON_THRESHOLD );

    //refining options

    //(configuration)
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ADD_NODES );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_INSERT_NODES );

    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ELEMENTS );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ELEMENTS_ON_DISTANCE );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ELEMENTS_ON_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_ELEMENTS_ON_THRESHOLD );

    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_BOUNDARY );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_BOUNDARY_ON_DISTANCE );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_BOUNDARY_ON_ERROR );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE_BOUNDARY_ON_THRESHOLD );

    //execution options      
    //(select)
    KRATOS_DEFINE_LOCAL_FLAG ( SELECT_ELEMENTS );
    KRATOS_DEFINE_LOCAL_FLAG ( PASS_ALPHA_SHAPE );
    KRATOS_DEFINE_LOCAL_FLAG ( ENGAGED_NODES );

   

    struct InfoParameters
    { 
   
     KRATOS_CLASS_POINTER_DEFINITION(InfoParameters);

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
	
	if( InsertedNodes > NumberOfNodes * 0.002 || RemovedNodes > NumberOfNodes * 0.002 ){
	  GeometricalSmoothingRequired = true;
	}
	else if( (InsertedNodes + RemovedNodes) > NumberOfNodes * 0.004 ){
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


    struct RefiningParameters
    {      

    KRATOS_CLASS_POINTER_DEFINITION(RefiningParameters);

    private:

      //Pointer variables
      const Variable< double >* mpThresholdVariable; 
      const Variable< double >* mpErrorVariable; 

    public: 

      Flags    RefiningOptions;     //configuration refining options
      Flags    RemovingOptions;     //configuration removing options
      
      int      NumberOfElements;     

      double   Alpha;               //critical alpha parameter

      double   CriticalRadius;      //critical area   size
      double   CriticalSide;        //critical length size

      
      double   ReferenceThreshold;  //critical variable threshold value
      double   ReferenceError;      //critical error percentage
     

      SpatialBoundingBox::Pointer  RefiningBox;
      bool     RefiningBoxSetFlag;

      
      // setting refining variables (generally for python interface)

      void SetRefiningOptions(const Flags&  rOptions)
      {
	RefiningOptions=rOptions;
      };

      void SetRemovingOptions(const Flags&  rOptions)
      {
	RemovingOptions=rOptions;
      };

      void SetAlphaParameter( const double rAlpha)
      {
	Alpha = rAlpha;
      };

      void SetCriticalRadius( const double rCriticalRadius )
      {
	CriticalRadius = rCriticalRadius;
      };

      void SetCriticalSide( const double rCriticalSide )
      {
	CriticalSide = rCriticalSide;
      };

      void SetParameters (const double rAlpha, const double rCriticalRadius, const double rCriticalSide)
      {
	Alpha   = rAlpha;
	CriticalRadius = rCriticalRadius;
	CriticalSide = rCriticalSide;
      }

      void SetRefiningBox ( SpatialBoundingBox::Pointer pRefiningBox )
      {
	RefiningBoxSetFlag =true;
	RefiningBox = pRefiningBox;
      }

      void SetReferenceThreshold( const double rReferenceThreshold )
      {
	ReferenceThreshold = rReferenceThreshold;
      };

      void SetReferenceError( const double rReferenceError )
      {
	ReferenceError = rReferenceError;
      };


      void SetThresholdVariable (const Variable<double>& rVariable)
      {
	mpThresholdVariable = &rVariable;
      };

      void SetErrorVariable       (const Variable<double>& rVariable)
      {
	mpErrorVariable = &rVariable;
      };

      const Variable<double>& GetThresholdVariable()
      {
	return *mpThresholdVariable;
      };

      const Variable<double>& GetErrorVariable()
      {
	return *mpErrorVariable;
      };

      void Initialize (){
	
	NumberOfElements    = 0;
	Alpha               = 0;
	CriticalRadius      = 0;  
	CriticalSide        = 0;  
	ReferenceThreshold  = 0;
	ReferenceError      = 0;

      };

    };


    struct MeshingParameters
    {

    KRATOS_CLASS_POINTER_DEFINITION(MeshingParameters);
  
    protected:

      //Pointer variables
      const Element   *mpReferenceElement;
      const Condition *mpReferenceCondition;

    public:

      //General configuration flags
      Flags   Options;
      
      //Local execution flags
      Flags   ExecutionOptions;

      //General configuration variables
      double  AlphaParameter;
      double  OffsetFactor;

      SpatialBoundingBox::Pointer  MeshingBox;
      bool MeshingBoxSetFlag;

      TransferParametersType TransferVariables;
      bool TransferVariablesSetFlag;

      //Local execution variables
      bool    NodalIdsSetFlag;

      std::vector<int> NodalPreIds;
      std::vector<int> NodalNewIds;

      std::vector<int> PreservedElements;
      std::vector<std::vector<int> > NeighbourList; 

      InfoParameters::Pointer        Info;
      RefiningParameters::Pointer  Refine;
     
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

      void SetAlphaParameter( const double rAlpha)
      {
	AlphaParameter = rAlpha;
      };

      void SetOffsetFactor(const double rOffsetFactor)
      {
	OffsetFactor=rOffsetFactor;
      };

      void SetInfoParameters(InfoParameters::Pointer rInfo)
      {
	Info=rInfo;
      };

      void SetRefiningParameters(RefiningParameters::Pointer rRefine)
      {
	Refine=rRefine;
      };

      void SetMeshingBox(SpatialBoundingBox::Pointer pMeshingBox)
      {
	MeshingBoxSetFlag =true;
	MeshingBox = pMeshingBox;
      }


      void SetTransferParameters(TransferParametersType& rTransferVariables)
      {
	TransferVariablesSetFlag =true;
	TransferVariables = rTransferVariables;
      }


      void SetTransferVariable(const Variable<double>& rTransferVariable)
      {
	TransferVariablesSetFlag =true;
	TransferVariables.SetVariable(rTransferVariable);
      }

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

 
      void Initialize (){

	AlphaParameter = 0;

	OffsetFactor = 0;

	MeshingBoxSetFlag = false;

	TransferVariablesSetFlag = false;
  
	NodalIdsSetFlag = false;

	// RemeshInfo.Initialize();
	// Refine.Initialize();	
      };
    };


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelerUtilities() {} //

    /// Destructor.
    virtual ~ModelerUtilities() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
   
    void SetDomainLabels   (ModelPart& rModelPart);

    void BuildTotalMesh    (ModelPart& rModelPart, int EchoLevel);
    
    void CleanMeshFlags    (ModelPart& rModelPart, ModelPart::IndexType MeshId=0);

    void CleanRemovedNodes (ModelPart& rModelPart, ModelPart::IndexType MeshId=0);


    //*******************************************************************************************
    //*******************************************************************************************

    bool CheckSubdomain     (Geometry<Node<3> >& rGeometry);
    
    bool CheckInnerCentre   (Geometry<Node<3> >& rGeometry);
    
    bool CheckOuterCentre   (Geometry<Node<3> >& rGeometry, double& rOffsetFactor, bool& rSelfContact);

    bool CheckGeometryShape (Geometry<Node<3> >& rGeometry, double& rShape);         

    //*******************************************************************************************
    //*******************************************************************************************

      
    //computes geometry radius
    double& ComputeRadius   (double& rRadius, double& rVolume, std::vector<Vector >& rVertices,const unsigned int& dimension);

    //returns false if it should be removed
    bool AlphaShape         (double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension);

    //returns false if it should be removed
    bool ShrankAlphaShape   (double AlphaParameter, Geometry<Node<3> >& rGeometry, double& rOffsetFactor, const unsigned int dimension);

    //returns the nodal h relative to a single boundary node
    double FindBoundaryH    (Node<3>& BoundaryPoint);
    
    //writes a list of particles telling if they are set as boundary or not
    void CheckParticles     (ModelPart& rModelPart,ModelPart::IndexType MeshId=0);
    
    //*******************************************************************************************
    //*******************************************************************************************
   
    static inline double CalculateSideLength(PointType& P1,PointType& P2)
    {
      return sqrt( (P1.X()-P2.X())*(P1.X()-P2.X()) + (P1.Y()-P2.Y())*(P1.Y()-P2.Y()) );
    };

    static inline double CalculateTriangleRadius(Geometry< Node<3> >& rGeometry)
    {

      double L1 = CalculateSideLength (rGeometry[0],rGeometry[1]);
      double L2 = CalculateSideLength (rGeometry[1],rGeometry[2]);
      double L3 = CalculateSideLength (rGeometry[2],rGeometry[0]);
      
      
      double Area = rGeometry.Area();
      
      
      double Rcrit = Area*2/(L1+L2+L3);
      
      return Rcrit;
      
    };
  

    static inline double CalculateTriangleRadius(const double x0, const double y0,
						 const double x1, const double y1,
						 const double x2, const double y2,
						 double& Area)
    {
      
      double L1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
      double L2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
      double L3 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));     
      
      Area = fabs( 0.5 * ( (x0*y1) - (x0*y2) - (x1*y0) + (x1*y2) + (x2*y0) - (x2*y1) ) );
      
      // std::cout<< " Area "<<Area<<" L1 "<<L1<<" L2 "<<L2<<" L3 "<<L3<<std::endl;
      
      double Rcrit = Area*2/(L1+L2+L3);
      
      return Rcrit;
    }

    static inline double CalculateAverageSideLength(const double x0, const double y0,
						    const double x1, const double y1,
						    const double x2, const double y2) 
    { 
      double length_0 = sqrt( x0*x0 + y0*y0 );
      double length_1 = sqrt( x1*x1 + y1*y1 );
      double length_2 = sqrt( x2*x2 + y2*y2 );
      
      return 0.5*( length_0 + length_1 + length_2 );
    };

    //*******************************************************************************************
    //*******************************************************************************************				           

    bool CheckConditionInBox (Condition::Pointer& pCondition, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo);
    
    bool CheckElementInBox   (Element::Pointer& pElement, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo);
    
    bool CheckVerticesInBox  (Geometry<Node<3> >& rGeometry, SpatialBoundingBox& rRefiningBox, ProcessInfo& rCurrentProcessInfo);
    

    //*******************************************************************************************
    //*******************************************************************************************		

    Condition::Pointer FindMasterCondition(Condition::Pointer& pCondition, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found);

    Condition::Pointer FindMasterCondition(Condition::Pointer& pCondition, PointType& pSlaveNode, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found);
    

    //*******************************************************************************************
    //*******************************************************************************************		

    bool CheckContactActive    (GeometryType& rConditionGeometry, bool& rSemiActiveContact, std::vector<bool>& rSemiActiveNodes);
    
    bool CheckNodeCloseWallTip (std::vector<SpatialBoundingBox::Pointer>& rRigidWalls, PointType& rNode, ProcessInfo& rCurrentProcessInfo, double& rFactor);

    double CheckCriticalRadius (ModelPart& rModelPart, double& rCriticalRadius, unsigned int MeshId = 0);

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
    ModelerUtilities& operator=(ModelerUtilities const& rOther);

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

}; // Class ModelerUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      ModelerUtilities& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const ModelerUtilities& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_MODELER_UTILITIES_H_INCLUDED  defined 

