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

#if !defined(KRATOS_ERROR_METRICS_PROCESS)
#define KRATOS_ERROR_METRICS_PROCESS

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
 * @class MetricErrorProcess
 * @ingroup MeshingApplication
 * @brief This class is can be used to compute the metrics of the model part with a error already computed
 * @author Vicente Mataix Ferrandiz
 */
template<SizeType TDim>
class KRATOS_API(MESHING_APPLICATION) MetricErrorProcess
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
    typedef Node <3>                                                                NodeType;

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

    /// Pointer definition of MetricErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetricErrorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    MetricErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~MetricErrorProcess() override = default;

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
     * @brief We initialize the Metrics of the MMG sol using the Hessian Metric matrix approach
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
        return "MetricErrorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetricErrorProcess";
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

    double mMinSize;           /// The minimal size of the elements
    double mMaxSize;           /// The maximal size of the elements

    bool mSetElementNumber;    /// Determines if a target number of elements for the new mesh is set
    SizeType mElementNumber;   /// The target number of elements for the new mesh
    double mTargetError;       /// The overall target error for the new mesh
    bool mAverageNodalH;       /// Determines if the nodal h is averaged from the surrounding elements or if the lowest value is taken

    SizeType mEchoLevel;       /// The echo level

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method estimates the new element size
     */
    void CalculateElementSize();

    /**
     * @brief In this final step the metric is computed
     */
    void CalculateMetric();

    /**
     * @brief This computes the element size depending of the geometry and it assigns to the ELEMENT_H variable
     * @param itElement The element iterator
     */
    void ComputeElementSize(ElementItType itElement);

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
    MetricErrorProcess& operator=(MetricErrorProcess const& rOther);

    /// Copy constructor.
    //MetricErrorProcess(MetricErrorProcess const& rOther);

    ///@}
};// class MetricErrorProcess
///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, class TVarType>
inline std::istream& operator >> (std::istream& rIStream,
                                  MetricErrorProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MetricErrorProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_ERROR_METRICS_PROCESS defined */
