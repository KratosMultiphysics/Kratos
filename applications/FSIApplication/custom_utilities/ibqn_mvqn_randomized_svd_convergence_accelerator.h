//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR)
#define  KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR

// System includes

// External includes

// Project includes
#include "utilities/dense_svd_decomposition.h"

// Application includes
#include "ibqn_mvqn_convergence_accelerator.h"

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

/// Forward declaration of standard MVQN
/// This is required since the include ward avoids the inclusion of the standard MVQN
template<class TSparseSpace, class TDenseSpace>
class MVQNRandomizedSVDConvergenceAccelerator;

/** @brief Interface Block Newton convergence accelerator
 * Interface Block Newton equations convergence accelerator
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class IBQNMVQNRandomizedSVDConvergenceAccelerator: public IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( IBQNMVQNRandomizedSVDConvergenceAccelerator );

    typedef IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>      BaseType;
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::DenseVectorType                           VectorType;
    typedef typename BaseType::DenseVectorPointerType             VectorPointerType;

    typedef typename BaseType::DenseMatrixType                           MatrixType;
    typedef typename BaseType::DenseMatrixPointerType             MatrixPointerType;

    typedef MVQNRandomizedSVDConvergenceAccelerator<TSparseSpace, TDenseSpace>        MVQNType;
    typedef typename BaseType::MVQNPointerType                                 MVQNPointerType;

    typedef typename DenseSingularValueDecomposition<TDenseSpace>::Pointer DenseSVDPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new IBQNMVQNRandomizedSVDConvergenceAccelerator object
     * MVQN convergence accelerator Json settings constructor
     * @param rConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit IBQNMVQNRandomizedSVDConvergenceAccelerator(
        DenseSVDPointerType pDenseSVD,
        Parameters ConvAcceleratorParameters)
    : BaseType()
    {
        ConvAcceleratorParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        BaseType::SetInitialRelaxationOmega(ConvAcceleratorParameters["w_0"].GetDouble());

        // Set the subdomains MVQN randomized SVD convergence accelerator pointers
        // Note that we call the simplified constructor with default zero relaxation omega and IBQN switch
        const double cut_off_tol = ConvAcceleratorParameters["cut_off_tol"].GetDouble();
        const unsigned int jacobian_modes = ConvAcceleratorParameters["jacobian_modes"].GetInt();
        MVQNPointerType p_convergence_accelerator_left = Kratos::make_unique<MVQNType>(pDenseSVD, jacobian_modes, cut_off_tol);
        MVQNPointerType p_convergence_accelerator_right = Kratos::make_unique<MVQNType>(pDenseSVD, jacobian_modes, cut_off_tol);
        BaseType::SetLeftConvergenceAcceleratorPointer(p_convergence_accelerator_left);
        BaseType::SetRightConvergenceAcceleratorPointer(p_convergence_accelerator_right);
    }

    /**
     * Copy Constructor.
     */
    IBQNMVQNRandomizedSVDConvergenceAccelerator(const IBQNMVQNRandomizedSVDConvergenceAccelerator& rOther) = delete;

    /**
     * Destructor.
     */
    virtual ~IBQNMVQNRandomizedSVDConvergenceAccelerator(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Parameters GetDefaultParameters() const override
    {
        Parameters ibqn_mvqn_default_parameters(R"({
            "solver_type"            : "IBQN_MVQN_randomized_SVD",
            "jacobian_modes"         : 10,
            "w_0"                    : 0.825,
            "abs_cut_off_tol"        : 1e-8
        })");

        return ibqn_mvqn_default_parameters;
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


    ///@}
    ///@name Friends
    ///@{


    ///@}
}; /* Class IBQNMVQNRandomizedSVDConvergenceAccelerator */

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR defined */
