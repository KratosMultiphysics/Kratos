#include <omp.h>
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
#include <omp.h>

// Project includes
#include "meshing_application.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"
// Include the trees
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 

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
class GaussPointItem 
    : public Point<3>
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
    {
    }
    
    GaussPointItem(
        const array_1d<double, 3> Coords,
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        const double Weight
        ):Point<3>(Coords),
          mpConstitutiveLaw(pConstitutiveLaw),
          mWeight(Weight)
    {
    }

    ///Copy constructor  (not really required)
    GaussPointItem(const GaussPointItem& GP):
        Point<3>(GP),
        mpConstitutiveLaw(GP.mpConstitutiveLaw),
        mWeight(GP.mWeight)
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
    typedef ModelPart::ElementsContainerType              ElementsArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    
    // Type definitions for the tree
    typedef GaussPointItem                                        PointType;
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
        Parameters ThisParameters =  Parameters(R"({})")
        )
    :mrOriginMainModelPart(rOriginMainModelPart),
     mrDestinationMainModelPart(rDestinationMainModelPart),
     mDimension(rDestinationMainModelPart.GetProcessInfo()[DOMAIN_SIZE])
     {
        Parameters DefaultParameters = Parameters(R"(
            {
                "allocation_size"                      : 1000, 
                "bucket_size"                          : 4, 
                "search_factor"                        : 2, 
                "interpolation_type"                   : "LST", 
                "internal_variable_interpolation_list" :[]
            })" );
        
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
        mAllocationSize = ThisParameters["allocation_size"].GetInt();
        mBucketSize = ThisParameters["bucket_size"].GetInt();
        mSearchFactor = ThisParameters["search_factor"].GetDouble();
        mThisInterpolationType = ConvertInter(ThisParameters["interpolation_type"].GetString());
         
        if (ThisParameters["internal_variable_interpolation_list"].IsArray() == true)
        {
            auto VariableArrayList = ThisParameters["internal_variable_interpolation_list"];
            
            for (unsigned int iVar = 0; iVar < VariableArrayList.size(); iVar++)
            {
                mInternalVariableList.push_back(KratosComponents<Variable<double>>::Get(VariableArrayList[iVar].GetString()));
            }
        }
        else
        {
            std::cout << "WARNING:: No variables to interpolate, look that internal_variable_interpolation_list is correctly defined in your parameters" << std::endl;
            mInternalVariableList.clear();
        }
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
        /** NOTE: There are mainly two ways to interpolate the internal variables (there are three, but just two are behave correctly)
         * CPT: Closest point transfer. It transfer the values from the closest GP
         * LST: Least-square projection transfer. It transfers from the closest GP from the old mesh
         * SFT: It transfer GP values to the nodes in the old mesh and then interpolate to the new mesh using the sahpe functions all the time (NOTE: THIS DOESN"T WORK, AND REQUIRES EXTRA STORE)
         */ 
        
        if (mThisInterpolationType == CPT && mInternalVariableList.size() > 0)
        {
            InterpolateGaussPointsCPT();
        }
        else if (mThisInterpolationType == LST && mInternalVariableList.size() > 0)
        {                        
            InterpolateGaussPointsLST();
        }
        else if (mThisInterpolationType == SFT && mInternalVariableList.size() > 0)
        {
            InterpolateGaussPointsSFT();
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
    
    // The model parts
    ModelPart& mrOriginMainModelPart;                    // The origin model part
    ModelPart& mrDestinationMainModelPart;               // The destination model part
    const unsigned int mDimension;                       // Dimension size of the space
    
    // The allocation parameters
    unsigned int mAllocationSize;                  // Allocation size for the vectors and max number of potential results
    unsigned int mBucketSize;                      // Bucket size for kd-tree
    
    // The seatch variables 
    double mSearchFactor;                          // The search factor to be considered
    PointVector mPointListOrigin;                        // A list that contents the all the gauss points from the origin modelpart 
    
    // Variables to interpolate
    std::vector<Variable<double>> mInternalVariableList; // The list of variables to interpolate
    InterpolationTypes mThisInterpolationType;           // The interpolation type considered

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This function creates a lists of gauss points ready for the search
     * @param ThisModelPart: The model part to consider
     */
    
    PointVector CreateGaussPointList(ModelPart& ThisModelPart)
    {
        PointVector ThisPointVector;
        
        GeometryData::IntegrationMethod ThisIntegrationMethod;
        
        // Iterate in the elements
        ElementsArrayType& pElements = ThisModelPart.Elements();
        auto numElements = pElements.end() - pElements.begin();
        
        const ProcessInfo& CurrentProcessInfo = ThisModelPart.GetProcessInfo();

        // Creating a buffer for parallel vector fill
        const unsigned int NumThreads = omp_get_max_threads();
        std::vector<PointVector> PointsBuffers(NumThreads);
        
        #pragma omp parallel
        {
            const unsigned int Id = omp_get_thread_num();
            
            #pragma omp for
            for(unsigned int i = 0; i < numElements; i++) 
            {
                auto itElem = pElements.begin() + i;
                
                // Getting the geometry
                Element::GeometryType& rThisGeometry = itElem->GetGeometry();
                
                // Getting the integration points
                ThisIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rThisGeometry.IntegrationPoints(ThisIntegrationMethod);
                const unsigned int IntegrationPointsNumber = IntegrationPoints.size();
                
                // Computing the Jacobian
                Vector VectorDetJ(IntegrationPointsNumber);
                rThisGeometry.DeterminantOfJacobian(VectorDetJ,ThisIntegrationMethod);
                
                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(IntegrationPointsNumber);
                itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,CurrentProcessInfo);
                    
                for (unsigned int iGaussPoint = 0; iGaussPoint < IntegrationPointsNumber; iGaussPoint++ )
                {
                    const array_1d<double, 3> LocalCoordinates = IntegrationPoints[iGaussPoint].Coordinates();
                    
                    // We compute the corresponding weight
                    const double Weight = VectorDetJ[iGaussPoint] * IntegrationPoints[iGaussPoint].Weight();
                    
                    // We compute the global coordinates
                    array_1d<double, 3> GlobalCoordinates;
                    GlobalCoordinates = rThisGeometry.GlobalCoordinates( GlobalCoordinates, LocalCoordinates );
                    
                    // We create the respective GP
                    PointTypePointer pPoint = PointTypePointer(new PointType(GlobalCoordinates, ConstitutiveLawVector[iGaussPoint], Weight));
                    (PointsBuffers[Id]).push_back(pPoint);
                }
            }
            
            // Combine buffers together
            #pragma omp single
            {
                for( auto& PointsBuffer : PointsBuffers) 
                {
                    std::move(PointsBuffer.begin(),PointsBuffer.end(),back_inserter(ThisPointVector));
                }
            }
        }
        
        return ThisPointVector;
    }

    /**
     * This method interpolate the values of the GP using the CPT method
     */
    
    void InterpolateGaussPointsCPT()
    {
        // We Initialize the process info
        const ProcessInfo& CurrentProcessInfo = mrDestinationMainModelPart.GetProcessInfo();
        
        // We update the list of points
        mPointListOrigin.clear();     
        mPointListOrigin = CreateGaussPointList(mrOriginMainModelPart);

        #pragma omp parallel firstprivate(mPointListOrigin)
        {
            // We initialize the intergration method
            GeometryData::IntegrationMethod ThisIntegrationMethod;
            
            // Create a tree
            // It will use a copy of mNodeList (a std::vector which contains pointers)
            // Copying the list is required because the tree will reorder it for efficiency
            KDTree TreePoints(mPointListOrigin.begin(), mPointListOrigin.end(), mBucketSize); 
            
            // Iterate over the destination elements
            ElementsArrayType& pElements = mrDestinationMainModelPart.Elements();
            auto numElements = pElements.end() - pElements.begin();
            
            #pragma omp for
            for(unsigned int i = 0; i < numElements; i++) 
            {
                auto itElem = pElements.begin() + i;
                
                // Getting the geometry
                Element::GeometryType& rThisGeometry = itElem->GetGeometry();
                
                // Getting the integration points
                ThisIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rThisGeometry.IntegrationPoints(ThisIntegrationMethod);
                const unsigned int IntegrationPointsNumber = IntegrationPoints.size();
                
                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(IntegrationPointsNumber);
                itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,CurrentProcessInfo);
            
                for (unsigned int iGaussPoint = 0; iGaussPoint < IntegrationPointsNumber; iGaussPoint++ )
                {
                    // We compute the global coordinates
                    const array_1d<double, 3> LocalCoordinates = IntegrationPoints[iGaussPoint].Coordinates();
                    array_1d<double, 3> GlobalCoordinates;
                    GlobalCoordinates = rThisGeometry.GlobalCoordinates( GlobalCoordinates, LocalCoordinates );
                    
                    PointTypePointer pGPOrigin = TreePoints.SearchNearestPoint(GlobalCoordinates);
                    
                    for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                    {
                        Variable<double> ThisVar = mInternalVariableList[iVar];
                        
                        double OriginValue;
                        OriginValue = (pGPOrigin->GetConstitutiveLaw())->GetValue(ThisVar, OriginValue);
                        
                        (ConstitutiveLawVector[iGaussPoint])->SetValue(ThisVar, OriginValue, CurrentProcessInfo);
                    }
                }
            }
        }
    }
    
    /**
     * This method interpolate the values of the GP using the LST method
     */
        
    void InterpolateGaussPointsLST()
    {
        // We Initialize the process info
        const ProcessInfo& CurrentProcessInfo = mrDestinationMainModelPart.GetProcessInfo();
        
        // We update the list of points
        mPointListOrigin.clear();     
        mPointListOrigin = CreateGaussPointList(mrOriginMainModelPart);
        
        #pragma omp parallel firstprivate(mPointListOrigin)
        {
            // We initialize the intergration method
            GeometryData::IntegrationMethod ThisIntegrationMethod;
            
            // Initialize values
            PointVector PointsFound(mAllocationSize);
            std::vector<double> PointsDistances(mAllocationSize);
            unsigned int NumberPointsFound = 0;
            
            // Create a tree
            // It will use a copy of mNodeList (a std::vector which contains pointers)
            // Copying the list is required because the tree will reorder it for efficiency
            KDTree TreePoints(mPointListOrigin.begin(), mPointListOrigin.end(), mBucketSize); 
            
            // Iterate over the destination elements
            ElementsArrayType& pElements = mrDestinationMainModelPart.Elements();
            auto numElements = pElements.end() - pElements.begin();
            
            #pragma omp for
            for(unsigned int i = 0; i < numElements; i++) 
            {
                auto itElem = pElements.begin() + i;
                
                // Getting the geometry
                Element::GeometryType& rThisGeometry = itElem->GetGeometry();
                
                // Getting the integration points
                ThisIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rThisGeometry.IntegrationPoints(ThisIntegrationMethod);
                const unsigned int IntegrationPointsNumber = IntegrationPoints.size();
                
                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(IntegrationPointsNumber);
                itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,CurrentProcessInfo);
            
                // Computing the radius
                const double Radius = mSearchFactor *  (mDimension == 2 ? std::sqrt(rThisGeometry.Area()) : std::cbrt(rThisGeometry.Volume()));
            
                // We get the NODAL_H vector
                Vector NodalHVector(rThisGeometry.size());
                for (unsigned int iNode = 0; iNode < rThisGeometry.size(); iNode++)
                {
                    if ( rThisGeometry[iNode].SolutionStepsDataHas( NODAL_H ) == false )
                    {
                        KRATOS_ERROR << "NODAL_H is not defined in the node ID: " << rThisGeometry[iNode].Id() << std::endl;
                    }
                    
                    NodalHVector[iNode] = rThisGeometry[iNode].FastGetSolutionStepValue(NODAL_H);
                }
                
                for (unsigned int iGaussPoint = 0; iGaussPoint < IntegrationPointsNumber; iGaussPoint++ )
                {
                    // We compute the global coordinates
                    const array_1d<double, 3> LocalCoordinates = IntegrationPoints[iGaussPoint].Coordinates();
                    array_1d<double, 3> GlobalCoordinates;
                    GlobalCoordinates = rThisGeometry.GlobalCoordinates( GlobalCoordinates, LocalCoordinates );
                    
                    // We compute the pondered characteristic length
                    Vector N( rThisGeometry.size() );
                    rThisGeometry.ShapeFunctionsValues( N, LocalCoordinates );
                    const double CharacteristicLength = inner_prod(N, NodalHVector);
                    
                    NumberPointsFound = TreePoints.SearchInRadius(GlobalCoordinates, Radius, PointsFound.begin(), PointsDistances.begin(), mAllocationSize);
                    
                    if (NumberPointsFound > 0)
                    {    
                        for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                        {
                            Variable<double> ThisVar = mInternalVariableList[iVar];
                            
                            double WeightingFunctionNumerator   = 0.0;
                            double WeightingFunctionDenominator = 0.0;
                            double OriginValue;
                            
                            for (unsigned int iPointFound = 0; iPointFound < NumberPointsFound; iPointFound++)
                            {
                                PointTypePointer pGPOrigin = PointsFound[iPointFound];
                                
                                const double Distance = PointsDistances[iPointFound];
                                
                                OriginValue = (pGPOrigin->GetConstitutiveLaw())->GetValue(ThisVar, OriginValue);
                                
                                const double PonderatedWeight = pGPOrigin->GetWeight() * std::exp( -4.0 * Distance * Distance /(CharacteristicLength * CharacteristicLength));
                                
                                WeightingFunctionNumerator   += PonderatedWeight * OriginValue;
                                WeightingFunctionDenominator += PonderatedWeight;
                            }
                            
                            const double DestinationValue = WeightingFunctionNumerator/WeightingFunctionDenominator;
                            
                            (ConstitutiveLawVector[iGaussPoint])->SetValue(ThisVar, DestinationValue, CurrentProcessInfo);
                        }
                    }
                    else
                    {
                        std::cout << "WARNING:: It wasn't impossible to find any Gauss Point from where interpolate the internal variables" << std::endl;
                    }
                }
            }
        }
    }
    
    /**
     * This method interpolate the values of the GP using the SFT method
     */
    
    void InterpolateGaussPointsSFT()
    {
        // Initialize some values
        GeometryData::IntegrationMethod ThisIntegrationMethod;
                
        // Iterate in the nodes to initialize the values
        NodesArrayType& pNode = mrOriginMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        /* Nodes */
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
            {
                Variable<double> ThisVar = mInternalVariableList[iVar];
                
                itNode->SetValue(ThisVar, 0.0);
            }
        }
        
        // Iterate in the elements to ponderate the values
        ElementsArrayType& pElementsOrigin = mrOriginMainModelPart.Elements();
        auto numElements = pElementsOrigin.end() - pElementsOrigin.begin();
        
        const ProcessInfo& OriginProcessInfo = mrOriginMainModelPart.GetProcessInfo();
        
        /* Elements */
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElementsOrigin.begin() + i;
            
            // Getting the geometry
            Element::GeometryType& rThisGeometry = itElem->GetGeometry();
            
            // Getting the integration points
            ThisIntegrationMethod = itElem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rThisGeometry.IntegrationPoints(ThisIntegrationMethod);
            const unsigned int IntegrationPointsNumber = IntegrationPoints.size();
            
            // Computing the Jacobian
            Vector VectorDetJ(IntegrationPointsNumber);
            rThisGeometry.DeterminantOfJacobian(VectorDetJ,ThisIntegrationMethod);
            
            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(IntegrationPointsNumber);
            itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,OriginProcessInfo);
            
            // We initialize the total weigth
            double TotalWeight = 0.0;
            
            for (unsigned int iGaussPoint = 0; iGaussPoint < IntegrationPointsNumber; iGaussPoint++ )
            {
                array_1d<double, 3> LocalCoordinates = IntegrationPoints[iGaussPoint].Coordinates();
                
                // We compute the corresponding weight
                const double Weight = VectorDetJ[iGaussPoint] * IntegrationPoints[iGaussPoint].Weight();
                TotalWeight += Weight;
                
                // We compute the pondered characteristic length
                Vector N( rThisGeometry.size() );
                rThisGeometry.ShapeFunctionsValues( N, LocalCoordinates );
                
                // We compute the global coordinates
                array_1d<double, 3> GlobalCoordinates;
                GlobalCoordinates = rThisGeometry.GlobalCoordinates( GlobalCoordinates, LocalCoordinates );
                
                for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                {
                    Variable<double> ThisVar = mInternalVariableList[iVar];
                    
                    double OriginValue;
                    OriginValue = ConstitutiveLawVector[iGaussPoint]->GetValue(ThisVar, OriginValue);
                    
                    // We sum all the contributions
                    for (unsigned int iNode = 0; iNode < rThisGeometry.size(); iNode++)
                    {
                        #pragma omp atomic
                        rThisGeometry[iNode].GetValue(ThisVar) += N[iNode] * OriginValue * Weight;
                    }
                }
            }
            
            // We divide by the total weight
            for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
            {
                Variable<double> ThisVar = mInternalVariableList[iVar];
                
                for (unsigned int iNode = 0; iNode < rThisGeometry.size(); iNode++)
                {
                    #pragma omp critical
                    rThisGeometry[iNode].GetValue(ThisVar) /= TotalWeight;
                }
            }
        }
        
        // We interpolate to the new nodes
        if (mDimension == 2)
        {
            // We create the locator
            BinBasedFastPointLocator<2> PointLocator = BinBasedFastPointLocator<2>(mrOriginMainModelPart);
            PointLocator.UpdateSearchDatabase();
            
            // Iterate in the nodes
            NodesArrayType& pNode = mrDestinationMainModelPart.Nodes();
            auto numNodes = pNode.end() - pNode.begin();
            
            /* Nodes */
            #pragma omp parallel for 
            for(unsigned int i = 0; i < numNodes; i++) 
            {
                auto itNode = pNode.begin() + i;
                
                Vector N;
                Element::Pointer pElement;
                
                const bool found = PointLocator.FindPointOnMeshSimplified(itNode->Coordinates(), N, pElement, mAllocationSize);
                
                if (found == false)
                {
                    std::cout << "WARNING: GP not found (interpolation not posible)" << std::endl;
                    std::cout << "\t X:"<< itNode->X() << "\t Y:"<< itNode->Y() << std::endl;
                }
                else
                {
                    for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                    {
                        Variable<double> ThisVar = mInternalVariableList[iVar];
                        
                        Vector Values(pElement->GetGeometry().size());
                        
                        for (unsigned int iNode = 0; iNode < pElement->GetGeometry().size(); iNode++)
                        {
                            Values[iNode] = pElement->GetGeometry()[iNode].GetValue(ThisVar);
                        }
                        
                        itNode->GetValue(ThisVar) = inner_prod(Values, N);
                    }
                }
            }
        }
        else
        {
            // We create the locator
            BinBasedFastPointLocator<3> PointLocator = BinBasedFastPointLocator<3>(mrOriginMainModelPart);
            PointLocator.UpdateSearchDatabase();
            
            // Iterate in the nodes
            NodesArrayType& pNode = mrDestinationMainModelPart.Nodes();
            auto numNodes = pNode.end() - pNode.begin();
            
            /* Nodes */
            #pragma omp parallel for 
            for(unsigned int i = 0; i < numNodes; i++) 
            {
                auto itNode = pNode.begin() + i;
                
                Vector N;
                Element::Pointer pElement;
                
                const bool found = PointLocator.FindPointOnMeshSimplified(itNode->Coordinates(), N, pElement, mAllocationSize);
                
                if (found == false)
                {
                    std::cout << "WARNING: Node "<< itNode->Id() << " not found (interpolation not posible)" << std::endl;
                    std::cout << "\t X:"<< itNode->X() << "\t Y:"<< itNode->Y() << "\t Z:"<< itNode->Z() << std::endl;
                }
                else
                {
                    for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                    {
                        Variable<double> ThisVar = mInternalVariableList[iVar];
                        
                        Vector Values(pElement->GetGeometry().size());
                        
                        for (unsigned int iNode = 0; iNode < pElement->GetGeometry().size(); iNode++)
                        {
                            Values[iNode] = pElement->GetGeometry()[iNode].GetValue(ThisVar);
                        }
                        
                        itNode->GetValue(ThisVar) = inner_prod(Values, N);
                    }
                }
            }
        }
        
        // Finally we interpolate to the new GP
        ElementsArrayType& pElementsDestination = mrDestinationMainModelPart.Elements();
        numElements = pElementsDestination.end() - pElementsDestination.begin();
        
        const ProcessInfo& DestinationProcessInfo = mrOriginMainModelPart.GetProcessInfo();
        
        /* Elements */
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElementsDestination.begin() + i;
            
            // Getting the geometry
            Element::GeometryType& rThisGeometry = itElem->GetGeometry();
            
            // Getting the integration points
            ThisIntegrationMethod = itElem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rThisGeometry.IntegrationPoints(ThisIntegrationMethod);
            const unsigned int IntegrationPointsNumber = IntegrationPoints.size();
            
            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(IntegrationPointsNumber);
            itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,DestinationProcessInfo);
            
            for (unsigned int iGaussPoint = 0; iGaussPoint < IntegrationPointsNumber; iGaussPoint++ )
            {
                array_1d<double, 3> LocalCoordinates = IntegrationPoints[iGaussPoint].Coordinates();
                
                // We compute the pondered characteristic length
                Vector N( rThisGeometry.size() );
                rThisGeometry.ShapeFunctionsValues( N, LocalCoordinates );
                
                // We compute the global coordinates
                array_1d<double, 3> GlobalCoordinates;
                GlobalCoordinates = rThisGeometry.GlobalCoordinates( GlobalCoordinates, LocalCoordinates );
                
                Vector Values(rThisGeometry.size() );
                
                for (unsigned int iVar = 0; iVar < mInternalVariableList.size(); iVar++)
                {
                    Variable<double> ThisVar = mInternalVariableList[iVar];
                    
                    for (unsigned int iNode = 0; iNode < rThisGeometry.size(); iNode++)
                    {
                        Values[iNode] = rThisGeometry[iNode].GetValue(ThisVar);
                    }
                    
                    const double DestinationValue = inner_prod(Values, N);
                    
                    ConstitutiveLawVector[iGaussPoint]->SetValue(ThisVar, DestinationValue, DestinationProcessInfo);
                }
            }
        }
    }
    
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
