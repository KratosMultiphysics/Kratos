//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_MESHER_UTILITIES_H_INCLUDED )
#define  KRATOS_MESHER_UTILITIES_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/spatial_bounding_box.hpp"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "delaunay_meshing_application_variables.h"

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
class KRATOS_API(DELAUNAY_MESHING_APPLICATION) MesherUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MesherUtilities
    KRATOS_CLASS_POINTER_DEFINITION( MesherUtilities );


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


    enum ContactElementType //(contact domain definition)
    {
      NonContact,   //this is not a contact element
      PointToFace,  //2D/3D classical well defined node to face contact
      EdgeToEdge,   //3D edge to edge complex contact element
      PointToPoint, //2D/3D not compatible element belong to multiple bodies
      Undefined     //to be defined later
    };


    /**
     * Flags related to the meshing parameters
     */

    //meshing options

    //(configuration)
    KRATOS_DEFINE_LOCAL_FLAG ( REMESH );
    KRATOS_DEFINE_LOCAL_FLAG ( REFINE );
    KRATOS_DEFINE_LOCAL_FLAG ( RECONNECT );
    KRATOS_DEFINE_LOCAL_FLAG ( TRANSFER );
    KRATOS_DEFINE_LOCAL_FLAG ( CONSTRAINED );
    KRATOS_DEFINE_LOCAL_FLAG ( CONTACT_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( MESH_SMOOTHING );
    KRATOS_DEFINE_LOCAL_FLAG ( VARIABLES_SMOOTHING );

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
    KRATOS_DEFINE_LOCAL_FLAG ( INITIALIZE_MESHER_INPUT );
    KRATOS_DEFINE_LOCAL_FLAG ( FINALIZE_MESHER_INPUT );

    KRATOS_DEFINE_LOCAL_FLAG ( TRANSFER_KRATOS_NODES_TO_MESHER );
    KRATOS_DEFINE_LOCAL_FLAG ( TRANSFER_KRATOS_ELEMENTS_TO_MESHER );
    KRATOS_DEFINE_LOCAL_FLAG ( TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER );
    KRATOS_DEFINE_LOCAL_FLAG ( TRANSFER_KRATOS_FACES_TO_MESHER );

    KRATOS_DEFINE_LOCAL_FLAG ( SELECT_TESSELLATION_ELEMENTS );
    KRATOS_DEFINE_LOCAL_FLAG ( KEEP_ISOLATED_NODES );

    //execution options (tessellation) //not needed any more, just the strings definition...to set.
    KRATOS_DEFINE_LOCAL_FLAG ( NEIGHBOURS_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( BOUNDARIES_SEARCH );
    KRATOS_DEFINE_LOCAL_FLAG ( SET_DOF );

    KRATOS_DEFINE_LOCAL_FLAG ( PASS_ALPHA_SHAPE );


    struct MeshContainer
    {

    KRATOS_CLASS_POINTER_DEFINITION(MeshContainer);

    protected:

      double* mpPointList;
      int*    mpElementList;

      double* mpElementSizeList;
      int*    mpElementNeighbourList;

      int     mNumberOfPoints;
      int     mNumberOfElements;

    public:


      //flags to set when the pointers are created (true) or deleted (false)
      bool PointListFlag;
      bool ElementListFlag;
      bool ElementSizeListFlag;
      bool ElementNeighbourListFlag;

      void SetPointList(double* &rPointList) { mpPointList = rPointList; }
      void SetElementList(int* &rElementList) { mpElementList = rElementList; };
      void SetElementSizeList(double* &rElementSizeList) { mpElementSizeList = rElementSizeList; };
      void SetElementNeighbourList(int* &rElementNeighbourList) { mpElementNeighbourList = rElementNeighbourList; };
      void SetNumberOfPoints(int& rNumberOfPoints) {mNumberOfPoints = rNumberOfPoints;};
      void SetNumberOfElements(int& rNumberOfElements) {mNumberOfElements = rNumberOfElements;};

      double* GetPointList() { return mpPointList; };
      int*    GetElementList() { return mpElementList; };
      double* GetElementSizeList() { return mpElementSizeList; };
      int*    GetElementNeighbourList() { return mpElementNeighbourList; };

      int&    GetNumberOfPoints()   {return mNumberOfPoints;};
      int&    GetNumberOfElements() {return mNumberOfElements;};

      void CreatePointList(const unsigned int NumberOfPoints, const unsigned int Dimension)
      {
	if( mpPointList ){
	  delete [] mpPointList;
	}
	mNumberOfPoints = NumberOfPoints;
	mpPointList     = new double[NumberOfPoints * Dimension];
	PointListFlag   = true;
      }

      void CreateElementList(const unsigned int NumberOfElements, const unsigned int NumberOfVertices)
      {
	if( mpElementList ){
	  delete [] mpElementList;
	}
	mNumberOfElements = NumberOfElements;
	mpElementList     = new int[NumberOfElements * NumberOfVertices];
	ElementListFlag   = true;
      }

      void CreateElementSizeList(const unsigned int NumberOfElements)
      {
	if( mpElementSizeList ){
	  delete [] mpElementSizeList;
	}
	mpElementSizeList     = new double[NumberOfElements];
	ElementSizeListFlag   = true;
      }

      void CreateElementNeighbourList(const unsigned int NumberOfElements, const unsigned int NumberOfFaces)
      {
	if( mpElementNeighbourList ){
	  delete [] mpElementNeighbourList;
	}
	mpElementNeighbourList     = new int[NumberOfElements * NumberOfFaces];
	ElementNeighbourListFlag   = true;
      }

      void Initialize()
      {
	mpPointList            = (double*) NULL;
	mpElementList          = (int*)    NULL;
	mpElementSizeList      = (double*) NULL;
	mpElementNeighbourList = (int*)    NULL;
	mNumberOfPoints        = 0;
	mNumberOfElements      = 0;

	PointListFlag            = false;
	ElementListFlag          = false;
	ElementSizeListFlag      = false;
	ElementNeighbourListFlag = false;
      }


      void Finalize()
      {
	if( mpPointList && PointListFlag ){
	  delete [] mpPointList;
	}

	if( mpElementList && ElementListFlag ){
	  delete [] mpElementList;
	}

	if( mpElementSizeList && ElementSizeListFlag ){
	  delete [] mpElementSizeList;
	}

	if( mpElementNeighbourList && ElementNeighbourListFlag ){
	  delete [] mpElementNeighbourList;
	}

	Initialize();
      }

    };


    struct MeshingInfoParameters
    {

     KRATOS_CLASS_POINTER_DEFINITION(MeshingInfoParameters);

    public:

      //initial and total
      unsigned int   NumberOfElements;
      unsigned int   NumberOfNodes;
      unsigned int   NumberOfConditions;

      //added
      unsigned int   NumberOfNewElements;
      unsigned int   NumberOfNewNodes;
      unsigned int   NumberOfNewConditions;

      //total for all refining boxes
      unsigned int   InsertedNodes;
      unsigned int   RemovedNodes;
      unsigned int   InsertedBoundaryNodes;
      unsigned int   InsertedBoundaryConditions;

      unsigned int   CriticalElements;

      bool  GeometricalSmoothingRequired;
      bool  MechanicalSmoothingRequired;

      double InitialMeshVolume;

      void Initialize (){

	NumberOfElements   = 0;
	NumberOfNodes      = 0;
	NumberOfConditions = 0;

	NumberOfNewElements = 0;
	NumberOfNewNodes = 0;
	NumberOfNewConditions = 0;

	InsertedNodes = 0;
	RemovedNodes  = 0;

	InsertedBoundaryNodes = 0;
	InsertedBoundaryConditions = 0;

	CriticalElements = 0;

	GeometricalSmoothingRequired = false;
	MechanicalSmoothingRequired  = false;

      };


      void SetNumberOfNodes(unsigned int NumberNodes)
      {
	NumberOfNodes = NumberNodes;
      }

      void SetNumberOfElements(unsigned int NumberElements)
      {
	NumberOfElements = NumberElements;
      }

      void SetNumberOfConditions(unsigned int NumberConditions)
      {
	NumberOfConditions = NumberConditions;
      }

      unsigned int GetNumberOfNodes()
      {
	return NumberOfNodes;
      }

      unsigned int GetNumberOfElements()
      {
	return NumberOfElements;
      }

      unsigned int GetNumberOfConditions()
      {
	return NumberOfConditions;
      }

      void SetNumberOfNewNodes(unsigned int NumberNodes)
      {
	NumberOfNewNodes = NumberNodes;
      }

      void SetNumberOfNewElements(unsigned int NumberElements)
      {
	NumberOfNewElements = NumberElements;
      }

      void SetNumberOfNewConditions(unsigned int NumberConditions)
      {
	NumberOfNewConditions = NumberConditions;
      }

      void SetInitialMeshVolume(double Volume)
      {
	InitialMeshVolume = Volume;
      }

      double GetInitialMeshVolume()
      {
	return InitialMeshVolume;
      }

      bool CheckGeometricalSmoothing(){

	//std::cout<<" Inserted Nodes "<<InsertedNodes<<" NumNodes "<<NumberOfNodes<<" Removed "<<RemovedNodes<<std::endl;

	if( InsertedNodes > NumberOfNodes * 0.001 || RemovedNodes > NumberOfNodes * 0.001 ){
	  GeometricalSmoothingRequired = true;
	}
	else if( (InsertedNodes + RemovedNodes) > NumberOfNodes * 0.002 ){
	  GeometricalSmoothingRequired = true;
	}
	else{
	  GeometricalSmoothingRequired = false;
	}

	return GeometricalSmoothingRequired;
      }

      bool CheckMechanicalSmoothing(){

	if( CriticalElements > 0 || NumberOfNewElements != 0 )
	  MechanicalSmoothingRequired = true;
	else
	  MechanicalSmoothingRequired = false;

	return MechanicalSmoothingRequired;
      }

      int GetInsertedNodes(){
	return InsertedNodes;
      }


    };

    struct RefineBodyInfo
    {
      //number of entities refined/removed due to criterion:
      unsigned int on_distance;
      unsigned int on_threshold;
      unsigned int on_error;

      void Initialize(){
	on_distance=0;
	on_threshold=0;
	on_error=0;
      }

      void EchoStats(){
	std::cout<<"   BodyInfo [on_distance:"<<on_distance<<", on_threshold:"<<on_threshold<<", on_error:"<<on_error<<"]"<<std::endl;
      }

    };

    struct RefineBoundaryInfo
    {
      //number of entities refined/removed due to criterion:
      unsigned int on_distance;
      unsigned int on_threshold;
      unsigned int on_error;

      unsigned int in_contact;
      unsigned int in_concave_boundary;
      unsigned int in_convex_boundary;

      void Initialize(){
	on_distance=0;
	on_threshold=0;
	on_error=0;
	in_contact=0;
	in_concave_boundary=0;
	in_convex_boundary=0;
      }

      void EchoStats(){
	std::cout<<"   BoundaryInfo [on_distance:"<<on_distance<<", on_threshold:"<<on_threshold<<", on_error:"<<on_error<<"]"<<std::endl;
	std::cout<<"                [in_contact:"<<in_contact<<", in_concave_boundary:"<<in_concave_boundary<<", in_convex_boundary:"<<in_convex_boundary<<"]"<<std::endl;
      }

    };


    struct RefiningInfoParameters
    {
      //refine/remove body
      RefineBodyInfo BodyElementsRefined; //refined elements
      RefineBodyInfo BodyNodesRemoved; //removed nodes

      //refine/remove boundary
      RefineBoundaryInfo BoundaryConditionsRefined; //refined boundary conditions
      RefineBoundaryInfo BoundaryNodesRemoved; //removed boundary nodes
    };


    struct RefiningParameters
    {

    KRATOS_CLASS_POINTER_DEFINITION(RefiningParameters);

    private:

      //Pointer variables
      const Variable< double >* mpThresholdVariable;
      const Variable< double >* mpErrorVariable;

    public:

      //reference sizes
      Flags    RefiningOptions;     //configuration refining options
      Flags    RemovingOptions;     //configuration removing options

      double   Alpha;               //critical alpha parameter

      double   CriticalRadius;      //critical area   size
      double   CriticalSide;        //critical length size

      double   ReferenceThreshold;  //critical variable threshold value
      double   ReferenceError;      //critical error percentage

      //computed sizes
      double   InitialRadius;       //initial mesh radius/nodal-h
      double   MeanVolume;          //mean element area/volume

      //info parameters
      RefiningInfoParameters Info;

      //applied in the spatial box
      bool     RefiningBoxSetFlag;
      SpatialBoundingBox::Pointer  RefiningBox;


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

      double GetMeanVolume()
      {
	return MeanVolume;
      };

      void SetMeanVolume( const double rMeanVolume )
      {
	MeanVolume = rMeanVolume;
      };

      void SetCriticalRadius( const double rCriticalRadius )
      {
	CriticalRadius = rCriticalRadius;
      };

      void SetInitialRadius(const double rInitialRadius)
      {
      	KRATOS_TRY
	  InitialRadius=rInitialRadius;
      	KRATOS_CATCH(" ")

      	  }

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

      Flags GetRefiningOptions()
      {
	return RefiningOptions;
      };

      Flags GetRemovingOptions()
      {
	return RemovingOptions;
      };


      double GetReferenceThreshold()
      {
	return ReferenceThreshold;
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

	Alpha               = 0;

	CriticalRadius      = 0;
	CriticalSide        = 0;

	ReferenceThreshold  = 0;
	ReferenceError      = 0;

	InitialRadius       = 0;
	MeanVolume          = 0;

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

      //SubModelPart Name
      std::string SubModelPartName;

      //General configuration flags
      Flags   Options;

      //Local execution flags
      Flags   ExecutionOptions;    //configuration meshing options

      //General configuration variables
      double  AlphaParameter;
      double  OffsetFactor;

      //Options for the mesher
      std::string TessellationFlags;
      std::string TessellationInfo;

      bool TransferVariablesSetFlag;

      //Local execution variablesg
      bool InputInitializedFlag;
      bool MeshElementsSelectedFlag;

      unsigned int ConditionMaxId;
      unsigned int ElementMaxId;
      unsigned int NodeMaxId;

      std::vector<int> NodalPreIds;
      std::vector<int> PreservedElements;

      std::vector<BoundedVector<double, 3> > Holes;

      //Mesher pointers to the mesh structures
      MeshContainer       InMesh;
      MeshContainer      OutMesh;
      MeshContainer      MidMesh;

      std::vector<std::vector<int> > NeighbourList;

      //Global Meshing info
      MeshingInfoParameters::Pointer      Info;

      //Global Remining parameters
      RefiningParameters::Pointer       Refine;

      //some local bbx-based refining parameters
      //can be defined here: std::vector<RefiningParameters::Pointer> LocalRefineVector;

      //Global Tranfer parameters
      TransferParametersType::Pointer Transfer;

      PropertiesType::Pointer       Properties;

      //Global Meshing box
      bool MeshingBoxSetFlag;
      SpatialBoundingBox::Pointer  MeshingBox;


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

      void SetSubModelPartName(std::string const& rSubModelPartName)
      {
	SubModelPartName = rSubModelPartName;
      };

      void SetExecutionOptions(const Flags&  rOptions)
      {
	ExecutionOptions=rOptions;
      };

      void SetTessellationFlags(std::string rFlags)
      {
	TessellationFlags=rFlags;
      };

      void SetTessellationInfo(std::string rInfo)
      {
	TessellationInfo=rInfo;
      };


      void SetAlphaParameter( const double rAlpha)
      {
	AlphaParameter = rAlpha;
      };

      void SetOffsetFactor(const double rOffsetFactor)
      {
	OffsetFactor=rOffsetFactor;
      };

      void SetInfoParameters(MeshingInfoParameters::Pointer rInfo)
      {
	Info=rInfo;
      };

      void SetRefiningParameters(RefiningParameters::Pointer rRefine)
      {
	Refine=rRefine;
      };


      void SetProperties(PropertiesType::Pointer rProperties)
      {
	Properties=rProperties;
      };


      void SetMeshingBox(SpatialBoundingBox::Pointer pMeshingBox)
      {
	MeshingBoxSetFlag =true;
	MeshingBox = pMeshingBox;
      }


      void SetTransferParameters(TransferParametersType::Pointer rTransferVariables)
      {
	TransferVariablesSetFlag =true;
	Transfer = rTransferVariables;
      }


      void SetTransferVariable(const Variable<double>& rTransferVariable)
      {
	TransferVariablesSetFlag =true;
	Transfer->SetVariable(rTransferVariable);
      }

      void SetReferenceElement   (const Element   & rElement)
      {
	mpReferenceElement=&rElement;
      };

      void SetReferenceCondition (const Condition & rCondition)
      {
	mpReferenceCondition=&rCondition;
      };

      void SetHoles(std::vector<BoundedVector<double, 3> >& rHoles)
      {
	Holes = rHoles;
      }

      std::vector<BoundedVector<double, 3> >& GetHoles()
      {
	return Holes;
      }

      std::string GetSubModelPartName()
      {
	return SubModelPartName;
      };

      Flags GetOptions()
      {
	return Options;
      };

      MeshingInfoParameters::Pointer GetInfoParameters()
      {
	return Info;
      };

      TransferParametersType::Pointer GetTransferParameters()
      {
	return Transfer;
      };

      RefiningParameters::Pointer GetRefiningParameters()
      {
	return Refine;
      };

      PropertiesType::Pointer GetProperties()
      {
	return Properties;
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

	OffsetFactor   = 0;

	MeshingBoxSetFlag        = false;

	TransferVariablesSetFlag = false;

	InputInitializedFlag     = false;
	MeshElementsSelectedFlag = false;

	ConditionMaxId = 0;
	ElementMaxId = 0;
	NodeMaxId = 0;

	InMesh.Initialize();
	OutMesh.Initialize();
	MidMesh.Initialize();

	// RemeshInfo.Initialize();
	// Refine.Initialize();
      };

      void InitializeMeshing(){
	MeshElementsSelectedFlag   = false;
        PreservedElements.clear();
	PreservedElements.resize(0);
      };

      void FinalizeMeshing(){
	MeshElementsSelectedFlag   = false;
        PreservedElements.clear();
	PreservedElements.resize(0);
      };

    };



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MesherUtilities() {} //

    /// Destructor.
    virtual ~MesherUtilities() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void SetModelPartNameToConditions (ModelPart& rModelPart);

    void SetModelPartNameToNodes (ModelPart& rModelPart);

    double ComputeModelPartVolume (ModelPart& rModelPart);

    //*******************************************************************************************
    //*******************************************************************************************

    bool CheckSubdomain     (Geometry<Node<3> >& rGeometry);

    bool CheckInnerCentre   (Geometry<Node<3> >& rGeometry);

    bool CheckOuterCentre   (Geometry<Node<3> >& rGeometry, double& rOffsetFactor, bool& rSelfContact);

    bool CheckSliver        (Geometry<Node<3> >& rGeometry);

    ContactElementType  CheckContactElement  (Geometry<Node<3> >& rGeometry, std::vector<int>& rSlaveVertices);

    double GetAndCompareSideLenghts  (Geometry<Node<3> >& rGeometry, double& rMaximumSideLength, double& rMinimumSideLength);

    bool CheckGeometryShape (Geometry<Node<3> >& rGeometry, int& rShape);

    //*******************************************************************************************
    //*******************************************************************************************


    //computes geometry radius
    double& ComputeRadius   (double& rRadius, double& rVolume, std::vector<Vector >& rVertices,const unsigned int& dimension);

    //returns false if it should be removed
    bool AlphaShape         (double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension);
    bool AlphaShape         (double AlphaParameter, Geometry<Node<3> >& rGeometry, const unsigned int dimension, const double MeanMeshSize);

    //returns false if it should be removed
    bool ShrankAlphaShape   (double AlphaParameter, Geometry<Node<3> >& rGeometry, double& rOffsetFactor, const unsigned int dimension);

    //returns the nodal h relative to a single boundary node
    double FindBoundaryH    (Node<3>& BoundaryPoint);

    //writes a list of particles telling if they are set as boundary or not
    void CheckParticles     (ModelPart& rModelPart);

    //computes velocity norms of the geometry
    bool CheckRelativeVelocities (Geometry<Node<3> >& rGeometry, const double& rRelativeFactor);

    //computes prediction of volume decrease of the geometry
    bool CheckVolumeDecrease(GeometryType& rVertices, const unsigned int& rDimension,const double& rTolerance,double& VolumeChange);

    //computes prediction of volume after a projection of the displacement geometry
    double GetMovedVolume(GeometryType& rVertices, const unsigned int& rDimension, double MovementFactor);

    //computes deformation gradient determinant
    double GetDeformationGradientDeterminant(GeometryType& rVertices, const unsigned int& rDimension);

    //*******************************************************************************************
    //*******************************************************************************************

    static inline double CalculateSideLength(PointType& P1,PointType& P2)
    {
      return sqrt( (P1.X()-P2.X())*(P1.X()-P2.X()) + (P1.Y()-P2.Y())*(P1.Y()-P2.Y()) + (P1.Z()-P2.Z())*(P1.Z()-P2.Z()));
    };

    static inline double CalculateBoundarySize(Geometry< Node<3> >& rGeometry)
    {

      if( rGeometry.size() == 2 ){
	return CalculateSideLength (rGeometry[0],rGeometry[1]);
      }
      else if( rGeometry.size() == 3 ){
	return 2*CalculateTriangleRadius(rGeometry);
      }
      else{
	std::cout<<"BOUNDARY SIZE NOT COMPUTED :: geometry not correct"<<std::endl;
	return 0;
      }

    };

    static inline double CalculateTriangleRadius(Geometry< Node<3> >& rGeometry)
    {

      double L1 = CalculateSideLength (rGeometry[0],rGeometry[1]);
      double L2 = CalculateSideLength (rGeometry[1],rGeometry[2]);
      double L3 = CalculateSideLength (rGeometry[2],rGeometry[0]);


      double Area = rGeometry.Area();

      //inradius
      double Rcrit = Area*2/(L1+L2+L3);

      return Rcrit;

    };

    static inline double CalculateTetrahedronRadius(Geometry< Node<3> >& rGeometry)
    {

      //edges
      double L1 = CalculateSideLength (rGeometry[0],rGeometry[1]);
      double L2 = CalculateSideLength (rGeometry[1],rGeometry[2]);
      double L3 = CalculateSideLength (rGeometry[2],rGeometry[3]);
      double L4 = CalculateSideLength (rGeometry[3],rGeometry[0]);
      double L5 = CalculateSideLength (rGeometry[3],rGeometry[1]);
      double L6 = CalculateSideLength (rGeometry[2],rGeometry[0]);

      //inradius
      double S   = 0.5*(L1+L4+L5);  //semiperimeter
      double R1  = sqrt( S*(S-L1)*(S-L4)*(S-L5) ) / S;

      S   = 0.5*(L2+L3+L5);  //semiperimeter
      double R2  = sqrt( S*(S-L2)*(S-L3)*(S-L5) ) / S;

      S   = 0.5*(L3+L4+L6);  //semiperimeter
      double R3  = sqrt( S*(S-L3)*(S-L4)*(S-L6) ) / S;

      S   = 0.5*(L1+L2+L6);  //semiperimeter
      double R4  = sqrt( S*(S-L1)*(S-L2)*(S-L6) ) / S;

      S = 1.0/(R1*R1) + 1.0/(R2*R2) + 1.0/(R3*R3) + 1.0/(R4*R4);

      double Rcrit = sqrt(2.0/S);  //this is always bigger than the inradius

      return Rcrit;

    };

    static inline double CalculateElementRadius(Geometry< Node<3> >& rGeometry)
    {

      if(  rGeometry.size() == 3 )
	return CalculateTriangleRadius( rGeometry );
      else
	return CalculateTetrahedronRadius( rGeometry );

    };

    static inline double CalculateElementRadius(Geometry< Node<3> >& rGeometry, double& rDomainSize)
    {

      if( rGeometry.size() == 3 )
	return CalculateTriangleRadius( rGeometry[0].X(), rGeometry[0].Y(),
					rGeometry[1].X(), rGeometry[1].Y(),
					rGeometry[2].X(), rGeometry[2].Y(),
					rDomainSize );
      else
	return CalculateTetrahedronRadius( rGeometry[0].X(), rGeometry[0].Y(), rGeometry[0].Z(),
					   rGeometry[1].X(), rGeometry[1].Y(), rGeometry[1].Z(),
					   rGeometry[2].X(), rGeometry[2].Y(), rGeometry[2].Z(),
					   rGeometry[3].X(), rGeometry[3].Y(), rGeometry[3].Z(),
					   rDomainSize );

    };

    static inline double CalculateTriangleArea(const double x0, const double y0,
					       const double x1, const double y1,
					       const double x2, const double y2)
    {
      return 0.5*( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) );
    }


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

      //inradius
      double Rcrit = Area*2/(L1+L2+L3);

      return Rcrit;
    }



    static inline double CalculateTetrahedronVolume(const double x0, const double y0, const double z0,
						    const double x1, const double y1, const double z1,
						    const double x2, const double y2, const double z2,
						    const double x3, const double y3, const double z3)
    {
      //volume
      double Volume = 0;

      Volume  = CalculateDeterminant(x1,y1,z1,x2,y2,z2,x3,y3,z3);
      Volume -= CalculateDeterminant(x0,y0,z0,x2,y2,z2,x3,y3,z3);
      Volume += CalculateDeterminant(x0,y0,z0,x1,y1,z1,x3,y3,z3);
      Volume -= CalculateDeterminant(x0,y0,z0,x1,y1,z1,x2,y2,z2);

      Volume *= (1.0/6.0);

      return Volume;
    }


    static inline double CalculateTetrahedronRadius(const double x0, const double y0, const double z0,
						    const double x1, const double y1, const double z1,
						    const double x2, const double y2, const double z2,
						    const double x3, const double y3, const double z3,
						    double& Volume)
    {

      //edges
      double L1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));
      double L2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
      double L3 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3));
      double L4 = sqrt((x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) + (z3-z0)*(z3-z0));
      double L5 = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1));
      double L6 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));

      //volume
      Volume = CalculateTetrahedronVolume(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

      //inradius
      double S   = 0.5*(L1+L4+L5);  //semiperimeter
      double R1  = sqrt( S*(S-L1)*(S-L4)*(S-L5) ) / S;

      S   = 0.5*(L2+L3+L5);  //semiperimeter
      double R2  = sqrt( S*(S-L2)*(S-L3)*(S-L5) ) / S;

      S   = 0.5*(L3+L4+L6);  //semiperimeter
      double R3  = sqrt( S*(S-L3)*(S-L4)*(S-L6) ) / S;

      S   = 0.5*(L1+L2+L6);  //semiperimeter
      double R4  = sqrt( S*(S-L1)*(S-L2)*(S-L6) ) / S;

      S = 1.0/(R1*R1) + 1.0/(R2*R2) + 1.0/(R3*R3) + 1.0/(R4*R4);

      double Rcrit = sqrt(2.0/S);  //this is always bigger than the inradius

      return Rcrit;
    }


    static inline double CalculateDeterminant(const double x0, const double y0, const double z0,
					      const double x1, const double y1, const double z1,
					      const double x2, const double y2, const double z2)
    {
      return (x0*y1*z2-x0*y2*z1-x1*y0*z2+x1*y2*z0+x2*y0*z1-x2*y1*z0);
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
    //*******************************************************************************************
    //*******************************************************************************************
    static inline bool CalculatePosition(const std::vector<std::vector<double> >& rPointCoordinates,
					 const std::vector<double>& rCenter, std::vector<double>& rShapeFunctionsN)
    {

      if( rPointCoordinates.size() == 3 ){

	return CalculatePosition( rPointCoordinates[0][0], rPointCoordinates[0][1],
				  rPointCoordinates[1][0], rPointCoordinates[1][1],
				  rPointCoordinates[2][0], rPointCoordinates[2][1],
				  rCenter[0], rCenter[1], rShapeFunctionsN );
      }
      else if( rPointCoordinates.size() == 4 ){

	return CalculatePosition( rPointCoordinates[0][0], rPointCoordinates[0][1], rPointCoordinates[0][2],
				   rPointCoordinates[1][0], rPointCoordinates[1][1], rPointCoordinates[1][2],
				   rPointCoordinates[2][0], rPointCoordinates[2][1], rPointCoordinates[2][2],
				   rPointCoordinates[3][0], rPointCoordinates[3][1], rPointCoordinates[3][2],
				   rCenter[0], rCenter[1], rCenter[2], rShapeFunctionsN );
      }
      else{
	 KRATOS_THROW_ERROR( std::logic_error,"Number of points supplied out of range ERROR", "" )

      }

      return false;

    }

    //*******************************************************************************************
    //*******************************************************************************************

    static inline bool CalculatePosition(const double& x0, const double& y0, const double& z0,
					 const double& x1, const double& y1, const double& z1,
					 const double& x2, const double& y2, const double& z2,
					 const double& x3, const double& y3, const double& z3,
					 const double& xc, const double& yc, const double& zc,
					 std::vector<double>& rShapeFunctionsN)
    {
      double volume = CalculateTetrahedronVolume(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

      if(volume < 1e-25)
	{
	  std::cout<<" ERROR LS: tetrahedral element with zero area found: "<<volume<<" position ("<<x0<<", "<<y0<<", "<<z0<<") ("<<x1<<", "<<y1<<", "<<z1<<") ("<<x2<<", "<<y2<<", "<<z2<<") ("<<x3<<", "<<y3<<", "<<z3<<") "<<std::endl;
	}

      if( rShapeFunctionsN.size() != 4 ){
	rShapeFunctionsN.resize(4);
	//std::fill( rShapeFunctionsN.begin(), rShapeFunctionsN.end(), 0 );
      }

      rShapeFunctionsN[0] = CalculateTetrahedronVolume(xc,yc,zc,x1,y1,z1,x2,y2,z2,x3,y3,z3) / volume;
      rShapeFunctionsN[1] = CalculateTetrahedronVolume(x0,y0,z0,xc,yc,zc,x2,y2,z2,x3,y3,z3) / volume;
      rShapeFunctionsN[2] = CalculateTetrahedronVolume(x0,y0,z0,x1,y1,z1,xc,yc,zc,x3,y3,z3) / volume;
      rShapeFunctionsN[3] = CalculateTetrahedronVolume(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) / volume;

      //std::cout<<" N "<<rShapeFunctionsN[0]<<" "<<rShapeFunctionsN[1]<<" "<<rShapeFunctionsN[2]<<" "<<rShapeFunctionsN[3]<<" "<<std::endl;

      double tol = 1e-5;
      double upper_limit = 1.0+tol;
      double lower_limit = -tol;

      if(rShapeFunctionsN[0] >= lower_limit && rShapeFunctionsN[1] >= lower_limit && rShapeFunctionsN[2] >= lower_limit && rShapeFunctionsN[3] >= lower_limit &&
	 rShapeFunctionsN[0] <= upper_limit && rShapeFunctionsN[1] <= upper_limit && rShapeFunctionsN[2] <= upper_limit && rShapeFunctionsN[3] <= upper_limit ) //if the xc yc zc is inside the tetrahedron
	return true;
      return false;
    }


    //*******************************************************************************************
    //*******************************************************************************************


    static inline bool CalculatePosition(const double& x0, const double& y0,
					 const double& x1, const double& y1,
					 const double& x2, const double& y2,
					 const double& xc, const double& yc,
					 std::vector<double>& rShapeFunctionsN)
    {
      double area = CalculateTriangleArea(x0,y0,x1,y1,x2,y2);

      if(area < 1e-15)
	{
	  //KRATOS_THROW_ERROR( std::logic_error,"element with zero area found", "" );
	  std::cout<<" ERROR LS: triangle element with zero area found: "<<area<<" position ("<<x0<<", "<<y0<<") ("<<x1<<", "<<y1<<") ("<<x2<<", "<<y2<<") "<<std::endl;
	}

      if( rShapeFunctionsN.size() != 3 ){
	rShapeFunctionsN.resize(3);
	//std::fill( rShapeFunctionsN.begin(), rShapeFunctionsN.end(), 0 );
      }

      rShapeFunctionsN[0] = CalculateTriangleArea(x1,y1,x2,y2,xc,yc) / area;
      rShapeFunctionsN[1] = CalculateTriangleArea(x2,y2,x0,y0,xc,yc) / area;
      rShapeFunctionsN[2] = CalculateTriangleArea(x0,y0,x1,y1,xc,yc) / area;

      double tol = 1e-5;
      double upper_limit = 1.0+tol;
      double lower_limit = -tol;

      if(rShapeFunctionsN[0] >= lower_limit && rShapeFunctionsN[1] >= lower_limit && rShapeFunctionsN[2] >= lower_limit && rShapeFunctionsN[0] <= upper_limit && rShapeFunctionsN[1] <= upper_limit && rShapeFunctionsN[2] <= upper_limit) //if the xc yc is inside the triangle
	return true;

      return false;
    }

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

    bool CheckContactActive(GeometryType& rConditionGeometry, bool& rSemiActiveContact, std::vector<bool>& rSemiActiveNodes);

    bool CheckContactCurvature(GeometryType& rConditionGeometry, std::vector<array_1d<double, 3> >& rContactNormals);

    double CheckCriticalRadius(ModelPart& rModelPart, double rCriticalRadius);

    double GetMeanRadius(ModelPart& rModelPart, double& rCriticalRadius);

    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxNodeId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      unsigned int max_id = rModelPart.Nodes().back().Id();

      for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node!= rModelPart.NodesEnd(); ++i_node)
	{
	  if(i_node->Id() > max_id)
	    max_id = i_node->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }


   //*******************************************************************************************
   //*******************************************************************************************

    static inline double CalculateModelPartVolume(ModelPart& rModelPart)
    {
      KRATOS_TRY

	const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
      double ModelPartVolume = 0;
      if( dimension == 2 ){

	for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
	  {
            if( i_elem->GetGeometry().LocalSpaceDimension() != dimension )
              ModelPartVolume += i_elem->GetGeometry().Area();
	  }
      }
      else if( dimension == 3 ){ //dimension == 3

	for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
	  {
            if( i_elem->GetGeometry().LocalSpaceDimension() != dimension )
                ModelPartVolume += i_elem->GetGeometry().Volume();
	  }
      }

      return ModelPartVolume;

      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxConditionId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      unsigned int max_id = rModelPart.Conditions().back().Id();

      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); ++i_cond)
	{
	  if(i_cond->Id() > max_id)
	    max_id = i_cond->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    static inline unsigned int GetMaxElementId(ModelPart& rModelPart)
    {
      KRATOS_TRY

      unsigned int max_id = rModelPart.Elements().back().Id();

      for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(); i_elem!= rModelPart.ElementsEnd(); ++i_elem)
	{
	  if(i_elem->Id() > max_id)
	    max_id = i_elem->Id();
	}

      return max_id;

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

    bool FindCondition(Geometry< Node<3> >& rConditionGeometry ,Geometry< Node<3> >& rGeometry, DenseMatrix<unsigned int>& lpofa, DenseVector<unsigned int>& lnofa, unsigned int& iface);

    //*******************************************************************************************
    //*******************************************************************************************

    /**
     *  Set Nodes to mesh
     */
    void SetNodes(ModelPart& rModelPart,
		  MeshingParameters& rMeshingVariables);

    /**
     * Set Elements to mesh
     */
    void SetElements(ModelPart& rModelPart,
		     MeshingParameters& rMeshingVariables);

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
        std::stringstream buffer;
        buffer << "MesherUtilities";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MesherUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "MesherUtilities Data";
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
    MesherUtilities& operator=(MesherUtilities const& rOther);

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

}; // Class MesherUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				      MesherUtilities& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MesherUtilities& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESHER_UTILITIES_H_INCLUDED  defined
