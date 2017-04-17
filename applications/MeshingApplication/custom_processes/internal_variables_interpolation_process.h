// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_INTERNAL_VARIABLES_INTERPOLATION_PROCESS )
#define  KRATOS_INTERNAL_VARIABLES_INTERPOLATION_PROCESS

// System includes

// External includes

// Project includes
#include "meshing_application.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product

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
    
    #if !defined(INTERPOLATION_TYPES)
    #define INTERPOLATION_TYPES
        enum InterpolationTypes {CPT = 0, LST = 1, SFT = 2};
    #endif
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** @brief Custom Gauss Point container to be used by the search
 */
class GaussPointItem: public Point<3>
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of GaussPointItem
    KRATOS_CLASS_POINTER_DEFINITION( GaussPointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    GaussPointItem():
        Point<3>()
    {
    }

    GaussPointItem(const array_1d<double, 3> Coords):
        Point<3>(Coords)
    {}
    
    GaussPointItem(
        const array_1d<double, 3> Coords,
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        double Weight
    ):
        Point<3>(Coords),
        mpConstitutiveLaw(pConstitutiveLaw),
        mWeight(Weight)
    {}

    ///Copy constructor  (not really required)
    GaussPointItem(const GaussPointItem& rhs):
        Point<3>(rhs),
        mpOriginCond(rhs.mpOriginCond)
    {
    }

    /// Destructor.
    ~GaussPointItem(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Returns the point
     * @return The point
     */
    
    Point<3> GetPoint()
    {
        Point<3> Point(this->Coordinates());
        
        return Point;
    }
    
    /**
     * Set the point
     * @param The point
     */
    
    void SetPoint(const Point<3> Point)
    {
        this->Coordinates() = Point.Coordinates();
    }

    /**
     * Sets the Constitutive Law associated to the point
     * @param pConstitutiveLaw: The pointer to the Constitutive Law
     */

    void SetConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    {
        mpConstitutiveLaw = pConstitutiveLaw;
    }
    
    /**
     * Returns the Constitutive Law associated to the point
     * @return mpConstitutiveLaw: The pointer to the Constitutive Law associated to the point
     */

    ConstitutiveLaw::Pointer GetConstitutiveLaw()
    {
        return mpConstitutiveLaw;
    }
    
    /**
     * Sets the Constitutive Law associated to the point
     * @param pConstitutiveLaw: The pointer to the Constitutive Law
     */

    void SetConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    {
        mpConstitutiveLaw = pConstitutiveLaw;
    }
    
    /**
     * Returns the integration weigth associated to the point
     * @return mWeight: The pointer to the Constitutive Law associated to the point
     */

    double GetWeight()
    {
        return mWeight;
    }
    
    /**
     * Sets the integration weigth associated to the point
     * @param Weight: The integration weight
     */

    void SetWeight(double Weight)
    {
        mWeight = Weight;
    }
    
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

    ConstitutiveLaw::Pointer mpConstitutiveLaw; // The constitutive law pointer
    double mWeight;                             // The integration weight of the GP

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class GaussPointItem 

/** \brief InternalVariablesInterpolationProcess
 * This utilitiy has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * The process employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 */

class InternalVariablesInterpolationProcess 
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    
    // Type definitions for the tree
    typedef GaussGaussPointItem                                        PointType;
    typedef PointType::Pointer                             PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of InternalVariablesInterpolationProcess
    KRATOS_CLASS_POINTER_DEFINITION( InternalVariablesInterpolationProcess );
      
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    
    /**
     * The constructor of the search utility uses the following inputs:
     * @param rOriginMainModelPart: The model part from where interpolate values
     * @param rDestinationMainModelPart: The model part where we want to interpolate the values
     * @param ThisParameters: The parameters containing all the information needed
     */
    
    InternalVariablesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters = Parameters("{'allocation_size': 1000, 'bucket_size': 4, 'search_factor': 2, 'interpolation_type': 'LST', 'internal_variable_interpolation_list':[]}")
        )
    :mrOriginMainModelPart(rOriginMainModelPart),
     mrDestinationMainModelPart(rDestinationMainModelPart),
     mDimension(rMainModelPart.GetProcessInfo()[DOMAIN_SIZE]),
     mAllocationSize(ThisParameters["allocation_size"].GetInt()),
     mBucketSize(ThisParameters["bucket_size"].GetInt()),
     mSearchFactor(ThisParameters["search_factor"].GetDouble()),
     mThisInterpolationType(ConvertInter(ThisParameters["interpolation_type"].GetString()))
    {        
        // TODO: Add somethig if necessary
    }
    
    virtual ~InternalVariablesInterpolationProcess(){};

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * We execute the search relative to the old and new model part
     */
    
    virtual void Execute()
    {
    }
        
    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param TypeSearch: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
    
    void SearchGaussPoints()
    {  
        /** NOTE: There are mainly two ways to interpolate the internal variables (there are three, but just two are behave correctly)
         * CPT: Closest point transfer. It transfer the values from the closest GP
         * LST: Least-square projection transfer. It transfers from the closest GP from the old mesh
         * SFT: It transfer GP values to the nodes in the old mesh and then interpolate to the new mesh using the sahpe functions all the time (NOTE: THIS DOESN'T WORK, AND REQUIRES EXTRA STORE)
         */ 
        
        // We update the list of points
        UpdateGaussPointList();
        
        // Initialize values
        PointVector PointsFound(mAllocationSize);
        unsigned int NumberPointsFound = 0;    
        
        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTree TreePoints(mPointListDestination.begin(), mPointListDestination.end(), mBucketSize); 
        
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrMainModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            if (itCond->Is(SLAVE) == true)
            {
                Point<3> Center;
                const double SearchRadius = mSearchFactor * ContactUtilities::CenterAndRadius((*itCond.base()), Center);

                NumberPointsFound = TreePoints.SearchInRadius(Center, SearchRadius, PointsFound.begin(), mAllocationSize);
                
                if (NumberPointsFound > 0)
                {   
                    
                    for(unsigned int i = 0; i < NumberPointsFound; i++)
                    {   
                        Condition::Pointer pCondOrigin = PointsFound[i]->GetCondition();
                        
                        // TODO: Do something
                    }
                }
            }
        }
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

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    virtual std::string Info() const
    {
        return "InternalVariablesInterpolationProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
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
    
    ModelPart& rOriginMainModelPart;             // The origin model part
    ModelPart& rDestinationMainModelPart;        // The destination model part
    const unsigned int mDimension;               // Dimension size of the space
    const unsigned int mAllocationSize;          // Allocation size for the vectors and max number of potential results
    const unsigned int mBucketSize;              // Bucket size for kd-tree
    const double mSearchFactor;                  // The search factor to be considered
    std::vector<Variable<double>> mVariableList; // The list of variables to interpolate
    PointVector mPointListDestination;           // A list that contents the all the points (from nodes) from the modelpart 
    InterpolationTypes mThisInterpolationType;   // The interpolation type considered

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This converts the interpolation string to an enum
     * @param str: The string that you want to comvert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */
        
    InterpolationTypes ConvertInter(const std::string& str)
    {
        if(str == "CPT") 
        {
            return CPT;
        }
        else if(str == "LST") 
        {
            return LST;
        }
        else if(str == "SFT") 
        {
            return SFT;
        }
        else
        {
            return LST;
        }
    }
    
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

}; // Class InternalVariablesInterpolationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  InternalVariablesInterpolationProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InternalVariablesInterpolationProcess& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_INTERNAL_VARIABLES_INTERPOLATION_PROCESS  defined 
