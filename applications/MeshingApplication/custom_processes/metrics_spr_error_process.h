// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Anna Rehr
//  Co-author   :    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SPR_ERROR_METRICS_PROCESS)
#define KRATOS_SPR_ERROR_METRICS_PROCESS

// System includes
#include <omp.h>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "meshing_application.h"
#include "processes/find_nodal_neighbours_process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
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

    /// Definition of the size type
    typedef std::size_t                                                             SizeType;

    /// Definition of the spaces
    typedef UblasSpace<double, CompressedMatrix, Vector>                     SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>                                LocalSpaceType;

    /// The definition of linear solvers
    typedef LinearSolver<SparseSpaceType, LocalSpaceType>                   LinearSolverType;
    
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
 * @class SPRMetricProcess
 * @ingroup MeshingApplication
 * @brief This class is can be used to compute the metrics of the model part with a superconvergent patch recovery approach
 * @tparam TDim The dimension to be computed
 * @author Anna Rehr
 */
template<SizeType TDim>
class SPRMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of SPRMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(SPRMetricProcess);
    
    /// The Voigt notation size
    static constexpr SizeType SigmaSize = (TDim == 2) ? 3 : 6;

    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    
    SPRMetricProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})"),
        LinearSolverType::Pointer pLinearSolver = nullptr
        );
    
    /// Destructor.
    virtual ~SPRMetricProcess() {}
    
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
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
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
    virtual std::string Info() const override
    {
        return "SPRMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SPRMetricProcess";
    }

    /// Print object"s data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    
    ModelPart& mThisModelPart;                /// The model part to compute

    double mMinSize;                          /// The minimal size of the elements
    double mMaxSize;                          /// The maximal size of the elements

    double mPenaltyNormal;                    /// The normal penalty
    double mPenaltyTangent;                   /// The tangent penalty

    SizeType mEchoLevel;                      /// The echo level

    bool mSetElementNumber;                   /// Determines if a target number of elements for the new mesh is set
    SizeType mElementNumber;                  /// The target number of elements for the new mesh
    double mTargetError;                      /// The overall target error for the new mesh
    bool mAverageNodalH;                      /// Determines if the nodal h is averaged from the surrounding elements or if the lowest value is taken

    LinearSolverType::Pointer mpLinearSolver; /// The linear solver considered
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * 
     */
    double SuperconvergentPatchRecovery();

    /**
     * 
     */
    void CalculatePatch(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

    /** Calculates the recovered stress at a node in the case of a standard patch without contact BC
    * @param itNode the node for which the recovered stress should be calculated
    * @param itPatchNode the center node of the patch
    */
    void CalculatePatchStandard(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

    /**
     * It calculates the recovered stress at a node where contact BCs are regarded
     * @param itNode the node for which the recovered stress should be calculated
     * @param itPatchNode the center node of the patch
     */
    void CalculatePatchContact(
        NodeItType itNode,
        NodeItType itPatchNode,
        SizeType NeighbourSize,
        Vector& rSigmaRecovered
        );

    /**
     * Sets the element size
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
    SPRMetricProcess& operator=(SPRMetricProcess const& rOther)
    {
        return *this;
    };

    /// Copy constructor.
    //SPRMetricProcess(SPRMetricProcess const& rOther);

};// class SPRMetricProcess

};// namespace Kratos.
#endif /* KRATOS_SPR_ERROR_METRICS_PROCESS defined */
