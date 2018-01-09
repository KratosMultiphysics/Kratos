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
#include "utilities/openmp_utils.h"
#include "meshing_application.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "custom_includes/gauss_point_item.h"
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
        );

    ~InternalVariablesInterpolationProcess() override= default;

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

    void Execute() override;

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

    std::string Info() const override
    {
        return "InternalVariablesInterpolationProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
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
    unsigned int mAllocationSize;                        // Allocation size for the vectors and max number of potential results
    unsigned int mBucketSize;                            // Bucket size for kd-tree

    // The seatch variables
    double mSearchFactor;                                // The search factor to be considered
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

    PointVector CreateGaussPointList(ModelPart& ThisModelPart);

    /**
     * This method interpolate the values of the GP using the CPT method
     */

    void InterpolateGaussPointsCPT();

    /**
     * This method interpolate the values of the GP using the LST method
     */

    void InterpolateGaussPointsLST();

    /**
     * This method interpolate the values of the GP using the SFT method
     */

    void InterpolateGaussPointsSFT();

        // Iterate in the nodes to initialize the values
        NodesArrayType& pNode = mrOriginMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();

        /* Nodes */
        #pragma omp parallel for
        for(int i = 0; i < numNodes; i++)
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
        for(int i = 0; i < numElements; i++)
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
#pragma omp critical
                    {
                      rThisGeometry[iNode].GetValue(ThisVar) += N[iNode] * OriginValue * Weight;
                    }
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
            for(int i = 0; i < numNodes; i++)
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
            for(int i = 0; i < numNodes; i++)
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
        for(int i = 0; i < numElements; i++)
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

    InterpolationTypes ConvertInter(const std::string& str);

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
