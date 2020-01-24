//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Saransh Saxena
//

#if !defined(KRATOS_SIMPLE_ERROR_CALCULATOR_PROCESS)
#define KRATOS_SIMPLE_ERROR_CALCULATOR_PROCESS

// System includes
#include <map>

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/ublas_interface.h"


namespace Kratos{
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

/**
 * @class MetricsTemperatureGradientProcess
 * @ingroup ConvectionDiffusionApplication
 * @brief This is an error estimator based on comparison of smoothed gauss-point gradiant integration on nodal quantities and 
 * the fe based nodal quantities.
 * @author Saransh Saxena
 */

template<std::size_t TDim>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) MetricsTemperatureGradientProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Containers definition
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    /// The definition of the node, element and geometry type
    typedef Node <3>                                                                NodeType;
    typedef Geometry<Node<3>>                                                   GeometryType;

    /// Definition of the indextype
    typedef std::size_t                                                            IndexType;

    /// Matrix and integration type definition
    typedef BoundedMatrix<double, TDim, TDim>                                     MatrixType;
    typedef GeometryType::ShapeFunctionsGradientsType      ShapeFunctionDerivativesArrayType;

    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Pointer definition of MetricsTemperatureGradientProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetricsTemperatureGradientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    MetricsTemperatureGradientProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters);

    /// Destructor.
    ~MetricsTemperatureGradientProcess() override = default;

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
     * @brief We initialize the Metrics of the MMG solver using the Simple Error Calculator approach
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MetricsTemperatureGradientProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetricsTemperatureGradientProcess";
    }

    /// Print object"s data.
    void PrintData(std::ostream& rOStream) const override
    {
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
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ModelPart& mrThisModelPart; /// The model part to compute

    std::map<int,int> mNodeMap;      /// Map of local indexing to global node index
    std::map<int,int> mElemMap;      /// Map of local indexing to global element index

    double mMinSize;                 /// The minimal size of the elements
    double mMaxSize;                 /// The maximal size of the elements
    Vector mNodalArea;               /// Vector to store Nodal Area
    Vector mElementError;            /// Vector to store Element Error
    bool mHistoricalResults;

    std::size_t mEchoLevel;             /// The echo level
    bool mNodalAveragingH;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    void CreateMap();

    /**
     * @brief This method computes the Nodal Area for the nodes in the Model Part
     */
    void CalculateNodalArea(Vector& nodal_area) const;

    /**
     * @brief This method computes the Geometry Data such as Shape Functions, Gauss Weights and Shape Function Derivatives
     */
    void CalculateGeomData(GeometryType& r_geom, Matrix& ShapeFunctions, ShapeFunctionDerivativesArrayType& ShapeDerivatives,Vector& DetJ, Vector& GaussWeights, unsigned int& NumGPoints);

    /**
     * @brief This method computes the temperature gradient at nodes using the ShapeFunctions
     */
    void CalculateNodalTempGradient(Vector& nodal_area);
    
    /**
     * @brief This method computes the Nodal error between the C1 and C0 continuous definitions
     */
    void CalculateElementError(Vector& nodal_area);

    /**
     * @brief In this final step the metric is computed for MMGProcess
     */
    void CalculateMetric();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MetricsTemperatureGradientProcess& operator=(MetricsTemperatureGradientProcess const& rOther);

    /// Copy constructor.
    //MetricsTemperatureGradientProcess(MetricsTemperatureGradientProcess const& rOther);

    ///@}
};// class MetricsTemperatureGradientProcess
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, class TVarType>
inline std::istream& operator >> (std::istream& rIStream,
                                  MetricsTemperatureGradientProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MetricsTemperatureGradientProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
};// namespace Kratos.
#endif /* KRATOS_SIMPLE_ERROR_CALCULATOR_PROCESS defined */
