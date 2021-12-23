// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_INTERNAL_VARIABLES_INTERPOLATION_PROCESS )
#define  KRATOS_INTERNAL_VARIABLES_INTERPOLATION_PROCESS

// System includes

// External includes

// Project includes
#include "meshing_application_variables.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "custom_includes/gauss_point_item.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"
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

    /// Definition of node type
    typedef Node<3> NodeType;

    /// Definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    /// Type definitions for the tree
    typedef GaussPointItem                                        PointType;
    typedef PointType::Pointer                             PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class InternalVariablesInterpolationProcess
 * @ingroup MeshingApplication
 * @brief This utilitiy has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * @details The process employs the projection.h from MeshingApplication, which works internally using a kd-tree
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(MESHING_APPLICATION) InternalVariablesInterpolationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Definitions for the variables
    typedef Variable<double>             DoubleVarType;
    typedef Variable<array_1d<double, 3>> ArrayVarType;
    typedef Variable<Vector>             VectorVarType;
    typedef Variable<Matrix>             MatrixVarType;

    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ElementsContainerType              ElementsArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;

    /// Pointer definition of InternalVariablesInterpolationProcess
    KRATOS_CLASS_POINTER_DEFINITION( InternalVariablesInterpolationProcess );

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enum it used to list the different types of interpolations available
     */
    enum class InterpolationTypes {
        CLOSEST_POINT_TRANSFER = 0, /// Closest Point Transfer. It transfer the values from the closest GP
        LEAST_SQUARE_TRANSFER = 1, /// Least-Square projection Transfer. It transfers from the closest GP from the old mesh
        SHAPE_FUNCTION_TRANSFER = 2  /// Shape Function Transfer. It transfer GP values to the nodes in the old mesh and then interpolate to the new mesh using the shape functions all the time
        };

    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rOriginMainModelPart The model part from where interpolate values
     * @param rDestinationMainModelPart The model part where we want to interpolate the values
     * @param ThisParameters The parameters containing all the information needed
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
     * @brief We execute the search relative to the old and new model part
     * @details There are mainly two ways to interpolate the internal variables (there are three, but just two are behave correctly)
     * - CLOSEST_POINT_TRANSFER: Closest Point Transfer. It transfer the values from the closest GP
     * - LEAST_SQUARE_TRANSFER: Least-Square projection Transfer. It transfers from the closest GP from the old mesh
     * - SHAPE_FUNCTION_TRANSFER: Shape Function Transfer. It transfer GP values to the nodes in the old mesh and then interpolate to the new mesh using the shape functions all the time
     * @note SFT THIS DOESN'T WORK, AND REQUIRES EXTRA STORE
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
    ///@name Private Type Definitions
    ///@{

    // Auxiliar search struct
    template<std::size_t TDim>
    struct auxiliar_search
    {
        auxiliar_search(ModelPart& rModelPart)
            : point_locator(rModelPart)
        {
            // We create the locator
            point_locator.UpdateSearchDatabase();
        }

        BinBasedFastPointLocator<TDim> point_locator;
        Vector N;
        Element::Pointer p_element;
    };

    ///@}
    ///@name Private Static Member Variables
    ///@{

    ///@}
    ///@name Private Member Variables
    ///@{

    // The model parts
    ModelPart& mrOriginMainModelPart;                    /// The origin model part
    ModelPart& mrDestinationMainModelPart;               /// The destination model part
    const std::size_t mDimension;                        /// Dimension size of the space

    // The allocation parameters
    std::size_t mAllocationSize;                    /// Allocation size for the vectors and max number of potential results
    std::size_t mBucketSize;                        /// Bucket size for kd-tree

    // The seatch variables
    double mSearchFactor;                           /// The search factor to be considered
    PointVector mPointListOrigin;                   /// A list that contents the all the gauss points from the origin modelpart

    // Variables to interpolate
    std::vector<std::string> mInternalVariableList; /// The list of internal variables to interpolate
    InterpolationTypes mThisInterpolationType;      /// The interpolation type considered

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function creates a lists of gauss points ready for the search
     * @param ThisModelPart The model part to consider
     */

    PointVector CreateGaussPointList(ModelPart& ThisModelPart);

    /**
     * @brief This method interpolate the values of the GP using the Closest Point Transfer method
     */

    void InterpolateGaussPointsClosestPointTransfer();

    /**
     * @brief This method interpolate the values of the GP using the Least-Square projection Transfer method
     */

    void InterpolateGaussPointsLeastSquareTransfer();

    /**
     * @brief This method interpolate the values of the GP using the Shape Function Transfer method
     */

    void InterpolateGaussPointsShapeFunctionTransfer();

    /**
     * @brief This method computes the total number of variables to been interpolated
     * @return The total number of variables to be interpolated
     */
    std::size_t ComputeTotalNumberOfVariables();

    /**
     * @brief This method saves the values on the gauss point object
     * @param rThisVar The variable to transfer
     * @param pPointOrigin The pointer to the current GP
     * @param itElemOrigin The origin element iterator to save on the auxiliar point
     * @param GaussPointId The index of te current GaussPoint computed
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void SaveValuesOnGaussPoint(
        const Variable<TVarType>& rThisVar,
        PointTypePointer pPointOrigin,
        ElementsArrayType::iterator itElemOrigin,
        const IndexType GaussPointId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        std::vector<TVarType> values;
        itElemOrigin->CalculateOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
        pPointOrigin->SetValue(rThisVar, values[GaussPointId]);
    }

    /**
     * @brief Simply gets a origin value from a CL and it sets on the destination CL
     * @param rThisVar The variable to transfer
     * @param pOriginConstitutiveLaw The CL on the original mesh
     * @param pDestinationConstitutiveLaw The Cl on the destination mesh
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void GetAndSetDirectVariableOnConstitutiveLaw(
        const Variable<TVarType>& rThisVar,
        ConstitutiveLaw::Pointer pOriginConstitutiveLaw,
        ConstitutiveLaw::Pointer pDestinationConstitutiveLaw,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        TVarType origin_value;
        origin_value = pOriginConstitutiveLaw->GetValue(rThisVar, origin_value);

        pDestinationConstitutiveLaw->SetValue(rThisVar, origin_value, rCurrentProcessInfo);
    }

    /**
     * @brief This method sets the value directly on the elementusing the value from the closest gauss point from the old mesh
     * @param rThisVar The variable to transfer
     * @param pPointOrigin The pointer to the current GP
     * @param itElemDestination The destination element iterato where to set the values
     * @param GaussPointId The index of te current GaussPoint computed
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void GetAndSetDirectVariableOnElements(
        const Variable<TVarType>& rThisVar,
        PointTypePointer pPointOrigin,
        ElementsArrayType::iterator itElemDestination,
        const IndexType GaussPointId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        std::vector<TVarType> values;
        itElemDestination->CalculateOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
        TVarType aux_value;
        values[GaussPointId] = pPointOrigin->GetValue(rThisVar, aux_value);
        itElemDestination->SetValuesOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
    }

    /**
     * @brief Gets a origin value from near points and it sets on the destination CL using a weighted proportion
     * @param rThisVar The variable to transfer
     * @param NumberOfPointsFound The number of points found during the search
     * @param PointsFound The list of points found
     * @param PointDistances The distances of the points found
     * @param CharacteristicLenght The characteristic length of the problem
     * @param pDestinationConstitutiveLaw The Cl on the destination mesh
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void GetAndSetWeightedVariableOnConstitutiveLaw(
        const Variable<TVarType>& rThisVar,
        const std::size_t NumberOfPointsFound,
        PointVector& PointsFound,
        const std::vector<double>& PointDistances,
        const double CharacteristicLenght,
        ConstitutiveLaw::Pointer pDestinationConstitutiveLaw,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        TVarType weighting_function_numerator = rThisVar.Zero();
        double weighting_function_denominator = 0.0;
        TVarType origin_value;

        for (std::size_t i_point_found = 0; i_point_found < NumberOfPointsFound; ++i_point_found) {
            PointTypePointer p_gp_origin = PointsFound[i_point_found];

            const double distance = PointDistances[i_point_found];

            origin_value = (p_gp_origin->GetConstitutiveLaw())->GetValue(rThisVar, origin_value);

            const double ponderated_weight = p_gp_origin->GetWeight() * std::exp( -4.0 * distance * distance /std::pow(CharacteristicLenght, 2));

            weighting_function_numerator   += ponderated_weight * origin_value;
            weighting_function_denominator += ponderated_weight;
        }

        const TVarType destination_value = weighting_function_numerator/weighting_function_denominator;

        pDestinationConstitutiveLaw->SetValue(rThisVar, destination_value, rCurrentProcessInfo);
    }

    /**
     * @brief Gets a origin value from near points and it sets on the destination CL using a weighted proportion
     * @param rThisVar The variable to transfer
     * @param NumberOfPointsFound The number of points found during the search
     * @param PointsFound The list of points found
     * @param PointDistances The distances of the points found
     * @param CharacteristicLenght The characteristic length of the problem
     * @param itElemDestination The destination element iterato where to set the values
     * @param GaussPointId The index of te current GaussPoint computed
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void GetAndSetWeightedVariableOnElements(
        const Variable<TVarType>& rThisVar,
        const std::size_t NumberOfPointsFound,
        PointVector& PointsFound,
        const std::vector<double>& PointDistances,
        const double CharacteristicLenght,
        ElementsArrayType::iterator itElemDestination,
        const IndexType GaussPointId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        TVarType weighting_function_numerator = rThisVar.Zero();
        double weighting_function_denominator = 0.0;
        TVarType origin_value;

        for (std::size_t i_point_found = 0; i_point_found < NumberOfPointsFound; ++i_point_found) {
            PointTypePointer p_gp_origin = PointsFound[i_point_found];

            const double distance = PointDistances[i_point_found];

            origin_value = p_gp_origin->GetValue(rThisVar, origin_value);

            const double ponderated_weight = p_gp_origin->GetWeight() * std::exp( -4.0 * distance * distance /std::pow(CharacteristicLenght, 2));

            weighting_function_numerator   += ponderated_weight * origin_value;
            weighting_function_denominator += ponderated_weight;
        }

        const TVarType destination_value = weighting_function_numerator/weighting_function_denominator;

        std::vector<TVarType> values;
        itElemDestination->CalculateOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
        values[GaussPointId] = destination_value;
        itElemDestination->SetValuesOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
    }

    /**
     * @brief This method interpolates and add values from  the CL using shape functions
     * @param rThisGeometry The geometry of the element
     * @param rThisVar The variable to transfer
     * @param N The shape function used
     * @param pConstitutiveLaw The CL on the original mesh
     * @param Weight The integration weight
     */
    template<class TVarType>
    static inline void InterpolateAddVariableOnConstitutiveLaw(
        GeometryType& rThisGeometry,
        const Variable<TVarType>& rThisVar,
        const Vector& N,
        ConstitutiveLaw::Pointer& pConstitutiveLaw,
        const double Weight
        );

    /**
     * @brief This method interpolates and add values from the element using shape functions
     * @param rThisGeometry The geometry of the element
     * @param rThisVar The variable to transfer
     * @param N The shape function used
     * @param rElement The origin element to save on the auxiliar point
     * @param GaussPointId The index of te current GaussPoint computed
     * @param Weight The integration weight
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void InterpolateAddVariableOnElement(
        GeometryType& rThisGeometry,
        const Variable<TVarType>& rThisVar,
        const Vector& N,
        Element& rElement,
        const IndexType GaussPointId,
        const double Weight,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method ponderates a value by the total integration weight
     * @param rThisGeometry The geometry of the element
     * @param rThisVar The variable to transfer
     * @param TotalWeight The total integration weight
     */
    template<class TVarType>
    static inline void PonderateVariable(
        GeometryType& rThisGeometry,
        const Variable<TVarType>& rThisVar,
        const double TotalWeight
        );

    /**
     * @brief This method interpolates using shape functions and the values from the elemental nodes
     * @param rThisVar The variable to transfer
     * @param N The shape function used
     * @param rNode The current node
     * @param pElement The pointer to the current element
     */
    template<class TVarType>
    static inline void InterpolateToNode(
        const Variable<TVarType>& rThisVar,
        const Vector& N,
        NodeType& rNode,
        Element::Pointer pElement
        )
    {
        // An auxiliar value
        TVarType aux_value = rThisVar.Zero();

        // Interpolate with shape function
        const std::size_t number_nodes = pElement->GetGeometry().size();
        for (std::size_t i_node = 0; i_node < number_nodes; ++i_node)
            aux_value += N[i_node] * pElement->GetGeometry()[i_node].GetValue(rThisVar);

        rNode.SetValue(rThisVar, aux_value);
    }

    /**
     * @brief Gets a origin value from near points and it sets on the destination CL using a weighted proportion
     * @param rThisGeometry The geometry of the element
     * @param rThisVar The variable to transfer
     * @param N The shape function used
     * @param pDestinationConstitutiveLaw The Cl on the destination mesh
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void SetInterpolatedValueOnConstitutiveLaw(
        GeometryType& rThisGeometry,
        const Variable<TVarType>& rThisVar,
        const Vector& N,
        ConstitutiveLaw::Pointer pDestinationConstitutiveLaw,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // An auxiliar value
        TVarType destination_value = rThisVar.Zero();

        // Interpolate with shape function
        const std::size_t number_nodes = rThisGeometry.size();
        for (std::size_t i_node = 0; i_node < number_nodes; ++i_node)
            destination_value += N[i_node] * rThisGeometry[i_node].GetValue(rThisVar);

        pDestinationConstitutiveLaw->SetValue(rThisVar, destination_value, rCurrentProcessInfo);
    }

    /**
     * @brief Gets a origin value from near points and it sets on the destination CL using a weighted proportion
     * @param rThisGeometry The geometry of the element
     * @param rThisVar The variable to transfer
     * @param N The shape function used
     * @param rElement The destination element where to set the values
     * @param GaussPointId The index of te current GaussPoint computed
     * @param rCurrentProcessInfo The process info
     */
    template<class TVarType>
    static inline void SetInterpolatedValueOnElement(
        GeometryType& rThisGeometry,
        const Variable<TVarType>& rThisVar,
        const Vector& N,
        Element& rElement,
        const IndexType GaussPointId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // An auxiliar value
        TVarType destination_value = rThisVar.Zero();

        // Interpolate with shape function
        const std::size_t number_nodes = rThisGeometry.size();
        for (std::size_t i_node = 0; i_node < number_nodes; ++i_node)
            destination_value += N[i_node] * rThisGeometry[i_node].GetValue(rThisVar);

        std::vector<TVarType> values;
        rElement.CalculateOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
        values[GaussPointId] = destination_value;
        rElement.SetValuesOnIntegrationPoints(rThisVar, values, rCurrentProcessInfo);
    }

    /**
     * @brief This method interpolates values to all the nodes from the old to the new mesh
     * @tparam TDim The dimension
     */
    template<std::size_t TDim>
    void InterpolateToNodes()
    {
        // Iterate over nodes
        NodesArrayType& r_nodes_array = mrDestinationMainModelPart.Nodes();

        /* Nodes */
        block_for_each(r_nodes_array, auxiliar_search<TDim>(mrOriginMainModelPart),
        [this](Node<3>& rNode, auxiliar_search<TDim>& aux) {

            const bool old_entity = rNode.IsDefined(OLD_ENTITY) ? rNode.Is(OLD_ENTITY) : false;
            if (!old_entity) {
                const bool found = aux.point_locator.FindPointOnMeshSimplified(rNode.Coordinates(), aux.N, aux.p_element, mAllocationSize);

                if (!found) {
                    KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING: Node "<< rNode.Id() << " not found (interpolation not posible)" <<  "\t X:"<< rNode.X() << "\t Y:"<< rNode.Y() << "\t Z:"<< rNode.Z() << std::endl;
                } else {
                    for (auto& variable_name : mInternalVariableList) {
                        if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                            const auto& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                            InterpolateToNode(r_variable, aux.N, rNode, aux.p_element);
                        } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                            const auto& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                            InterpolateToNode(r_variable, aux.N, rNode, aux.p_element);
                        } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                            const auto& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                            InterpolateToNode(r_variable, aux.N, rNode, aux.p_element);
                        } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                            const auto& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                            InterpolateToNode(r_variable, aux.N, rNode, aux.p_element);
                        } else {
                            KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                        }
                    }
                }
            }
        });
    }

    /**
     * @brief This converts the interpolation string to an enum
     * @param Str The string that you want to comvert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */
    InterpolationTypes ConvertInter(const std::string& Str);

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

template<> 
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<double>(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    );
template<> 
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<array_1d<double, 3>>(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    );
template<> 
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<Vector>(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    );
template<> 
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<Matrix>(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    );

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    );

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    );

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    );

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    );

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const double TotalWeight
    );

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const double TotalWeight
    );

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const double TotalWeight
    );

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const double TotalWeight
    );

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
