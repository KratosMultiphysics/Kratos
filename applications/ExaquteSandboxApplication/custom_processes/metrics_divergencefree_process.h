//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//                   Brendan Keith
//

#if !defined(KRATOS_DIVERGENCEFREE_METRICS_PROCESS)
#define KRATOS_DIVERGENCEFREE_METRICS_PROCESS

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
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
 * @class MetricDivergenceFreeProcess
 * @ingroup MeshingApplication
 * @brief This class is can be used to compute the metrics of the model part with different strategies, exploting the divergence.
 * The developed strategies are:
 * mean distribution strategy: this strategy builds a metric which tries to uniformly spread the refinement indicator in the whole domain
 * maximum strategy: this strategy build a metric which aims at refining only if the refinement indicator belong to the interval [coefficient*max(refinement indicator),max(refinement indicator)]
 * @author Riccardo Tosi
 */
template<SizeType TDim>
class KRATOS_API(MESHING_APPLICATION) MetricDivergenceFreeProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Containers definition
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    /// The definition of the node type
    typedef Node                                                                NodeType;

    /// Definition of the iterators
    typedef WeakPointerVector< Element >::iterator                         WeakElementItType;
    typedef NodesArrayType::iterator                                              NodeItType;
    typedef ElementsArrayType::iterator                                        ElementItType;

    /// Definition of the indextype
    typedef std::size_t                                                            IndexType;

    /// Matrix type definition
    typedef BoundedMatrix<double, TDim, TDim> MatrixType;

    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Pointer definition of MetricDivergenceFreeProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetricDivergenceFreeProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    MetricDivergenceFreeProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~MetricDivergenceFreeProcess() override = default;

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
     * @brief We initialize the Metrics of the MMG solver using the Divergence Free Metric matrix approach
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
        return "MetricDivergenceFreeProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetricDivergenceFreeProcess";
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
    enum RefinementStrategies {MaximumStrategy, MeanDistributionStrategy, GlobalToleranceStrategy};
    RefinementStrategies mRefinementStrategy; /// Refinement strategy
    SizeType mEchoLevel;             /// The echo level
    std::string mReferenceVariable;

    // Mean dsitribution strategy
    std::string mMeanStrategyReferenceNorm;
    double mMeanStrategyTargetRefinementCoefficient;
    double mMeanStrategyRefinementBound;
    double mMeanStrategyDivergenceFreeOverAllDomain;

    // Maximum strategy
    double mMaxStrategyTargetRefinementCoefficient;
    double mMaxStrategyRefinementCoefficient;
    double mDivergenceFreeMaxValue;

    // Interpolation error strategy
    double mGlobalErrorStrategyGlobalTolerance;
    double mGlobalErrorStrategyMeshConstant;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes each refinement strategy
     */
    void InitializeRefinementStrategy();

    /**
     * @brief In this final step the metric is computed
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
    MetricDivergenceFreeProcess& operator=(MetricDivergenceFreeProcess const& rOther);

    /// Copy constructor.
    //MetricDivergenceFreeProcess(MetricDivergenceFreeProcess const& rOther);

    ///@}
};// class MetricDivergenceFreeProcess
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, class TVarType>
inline std::istream& operator >> (std::istream& rIStream,
                                  MetricDivergenceFreeProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MetricDivergenceFreeProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_DIVERGENCEFREE_METRICS_PROCESS defined */
