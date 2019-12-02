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

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/element.h"
#include "processes/process.h"

namespace Kratos{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    typedef std::size_t SizeType;

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
 * @class SimpleErrorCalculatorProcess
 * @ingroup ConvectionDiffusionApplication
 * @brief This is an error estimator based on comparison of smoothed gauss-point gradiant integration on nodal quantities and 
 * the fe based nodal quantities.
 * @author Saransh Saxena
 */

template<SizeType TDim>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) SimpleErrorCalculatorProcess
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
    typedef Node <TDim>                                                             NodeType;
    typedef Geometry<Node<TDim>>                                                GeometryType;

    /// Definition of the iterators
    typedef WeakPointerVector< Element >::iterator                         WeakElementItType;
    typedef NodesArrayType::iterator                                              NodeItType;
    typedef ElementsArrayType::iterator                                        ElementItType;

    /// Definition of the indextype
    typedef std::size_t                                                            IndexType;

    /// Matrix and integration type definition
    typedef BoundedMatrix<double, TDim, TDim>                                     MatrixType;
    typedef GeometryType::ShapeFunctionsGradientsType      ShapeFunctionDerivativesArrayType;

    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Pointer definition of SimpleErrorCalculatorProcess
    KRATOS_CLASS_POINTER_DEFINITION(SimpleErrorCalculatorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    SimpleErrorCalculatorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~SimpleErrorCalculatorProcess() override = default;

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
        return "SimpleErrorCalculatorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleErrorCalculatorProcess";
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

    double mMinSize;                 /// The minimal size of the elements
    double mMaxSize;                 /// The maximal size of the elements

    SizeType mEchoLevel;             /// The echo level
    std::string mReferenceVariable;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the Error Ratio and Writes it to the nodes/elements
     */
    void ErrorEstimatorImplementation();

    // /**
    //  * @brief This method computes the Gauss Point Quantites for the element
    //  */
    // void CalculateGaussPointData(Vector& rGaussWeights,
    //                              Matrix& rNContainer,
    //                              ShapeFunctionDerivativesArrayType& rDN_DX);

    /**
     * @brief This method computes the temperature gradient at nodes using the ShapeFunctions
     */
    void CalculateNodalTempGradient();

    /**
     * @brief In this final step the metric is computed for MMGProcess
     */
    void CalculateMetricScalar();

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
    SimpleErrorCalculatorProcess& operator=(SimpleErrorCalculatorProcess const& rOther);

    /// Copy constructor.
    //SimpleErrorCalculatorProcess(SimpleErrorCalculatorProcess const& rOther);

    ///@}
};// class SimpleErrorCalculatorProcess
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, class TVarType>
inline std::istream& operator >> (std::istream& rIStream,
                                  SimpleErrorCalculatorProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SimpleErrorCalculatorProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
};// namespace Kratos.
#endif /* KRATOS_SIMPLE_ERROR_CALCULATOR_PROCESS defined */
