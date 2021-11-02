//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "factories/factory.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

/** @brief Butcher tableau for Runge-Kutta method.
 * 
 * Contains all info necessary of a particular RK method.
 * It specifies the coefficients of the particular Runge-Kutta method.
 *
 * Child class must provide methods:
 * - static const MatrixType() GenerateRKMatrix
 * - static const VectorType() GenerateWeights
 * - static const VectorType() GenerateNodes
 * - std::string Name() const
 */
template<unsigned int TOrder, typename ChildClass>
class ButcherTableau
{
public:
    typedef BoundedMatrix<double, TOrder-1, TOrder-1> MatrixType;
    typedef array_1d<double, TOrder> VectorType;

    constexpr unsigned int Order() {return TOrder;}
    constexpr unsigned int Size() {return TOrder-1;}
    
    constexpr MatrixRow<const MatrixType> GetRKMatrixRow(const unsigned int SubStepIndex)
    {
        return row(mA, SubStepIndex-1);
    }

    constexpr const VectorType& GetWeights()
    {
        return mB;
    }

    constexpr double GetNode(const unsigned int SubStepIndex)
    {
        return mC(SubStepIndex);
    }

    virtual std::string Name() const = 0;

protected:
    const MatrixType mA = ChildClass::GenerateRKMatrix();   // Runge-Kutta matrix
    const VectorType mB = ChildClass::GenerateWeights();    // Weights vector
    const VectorType mC = ChildClass::GenerateNodes();      // Nodes vector
};


class ButcherTableauForwardEuler : public ButcherTableau<1, ButcherTableauForwardEuler>
{
public:
    static const MatrixType GenerateRKMatrix()
    {
        MatrixType A = ZeroMatrix(0, 0);
        return A;
    }

    static const VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0;
        return B;
    }

    static const VectorType GenerateNodes()
    {
        VectorType C;
        C(0) = 0.0;
        return C;
    }

    std::string Name() const override
    {
        return "ButcherTableauForwardEuler";
    }
};


class ButcherTableauMidPointMethod : public ButcherTableau<2, ButcherTableauMidPointMethod>
{
public:
    static const MatrixType GenerateRKMatrix()
    {
        MatrixType A = ZeroMatrix(1, 1);
        A(0,0) = 0.5;
        return A;
    }

    static const VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 0.0;
        B(1) = 1.0;
        return B;
    }

    static const VectorType GenerateNodes()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 0.5;
        return C;
    }

    std::string Name() const override
    {
        return "ButcherTableauMidPointMethod";
    }
};

/** @brief Explicit total variation diminishing 3rd order Runge-Kutta
 *
 * @details Implementation of Proposition 3.2 of:
 * Gottlieb, Sigal, and Chi-Wang Shu. "Total variation diminishing Runge-Kutta schemes."
 * Mathematics of computation 67.221 (1998): 73-85.
 */
class ButcherTableauRK3TVD : public ButcherTableau<3, ButcherTableauRK3TVD>
{
public:
    static const MatrixType GenerateRKMatrix()
    {
        MatrixType A = ZeroMatrix(2, 2);
        A(0,0) = 1;
        A(1,0) = 0.25;
        A(1,1) = 0.25;
        return A;
    }

    static const VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0 / 6.0;
        B(1) = 1.0 / 6.0;
        B(2) = 2.0 / 3.0;
        return B;
    }

    static const VectorType GenerateNodes()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 1.0;
        C(2) = 0.5;
        return C;
    }

    std::string Name() const override
    {
        return "ButcherTableauRK3TVD";
    }
};


class ButcherTableauRK4 : public ButcherTableau<4, ButcherTableauRK4>
{
public:
    static const MatrixType GenerateRKMatrix()
    {
        MatrixType A = ZeroMatrix(3, 3);
        A(0,0) = 0.5;
        A(1,1) = 0.5;
        A(2,2) = 1.0;
        return A;
    }

    static const VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0 / 6.0;
        B(1) = 1.0 / 3.0;
        B(2) = 1.0 / 3.0;
        B(3) = 1.0 / 6.0;
        return B;
    }

    static const VectorType GenerateNodes()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 0.5;
        C(2) = 0.5;
        C(3) = 1.0;
        return C;
    }

    std::string Name() const override
    {
        return "ButcherTableauRK4";
    }
};

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
 * 
 * Formulation: for i = 0...N substeps
 *
 * - The diferential equation is u_t = f(t, u)
 * 
 * The Runge-Kutta method is:
 * - k^(i) = f(t + c_i*dt, u^(i-1))             -> Where u^(0) := u^n
 * - u^(i) = u^n + \sum_{j=1}^i A_{ij} k^(j)    -> Intermediate steps. u^(N) is not needed, therefore neither is A[N, :]
 * - u^{n+1} = u^n + \sum_{i} B_{i} k(u^(i))    -> Solution
 *
 * 
 */
template <class TSparseSpace, class TDenseSpace, class TButcherTableau>
class ExplicitSolvingStrategyRungeKutta : public ExplicitSolvingStrategy<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    // The base solving strategy class definition
    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    // The base class definition
    typedef ExplicitSolvingStrategy<TSparseSpace, TDenseSpace> BaseType;

    /// The definition of the current class
    typedef ExplicitSolvingStrategyRungeKutta<TSparseSpace, TDenseSpace, TButcherTableau> ClassType;

    // The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The DOF type
    typedef typename BaseType::DofType DofType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategyRungeKutta);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit ExplicitSolvingStrategyRungeKutta()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilder, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit ExplicitSolvingStrategyRungeKutta(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategyRungeKutta(const ExplicitSolvingStrategyRungeKutta &Other) = delete;

    /** Destructor.
     */
    ~ExplicitSolvingStrategyRungeKutta() override = default;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "explicit_solving_strategy_runge_kutta"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "explicit_solving_strategy_runge_kutta";
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "ExplicitSolvingStrategyRungeKutta with tableau " << mButcherTableau.Name();
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    const TButcherTableau mButcherTableau;


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void SolveWithLumpedMassMatrix() override
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const unsigned int dof_size = p_explicit_bs->GetEquationSystemSize();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Set the auxiliary RK vectors
        LocalSystemVectorType u_n(dof_size); // TODO: THIS IS INEFICCIENT. CREATE A UNORDERED_SET WITH THE IDOF AND VALUE AS ENTRIES. THIS HAS TO BE OPTIONAL
        LocalSystemMatrixType rk_K(dof_size, mButcherTableau.Order());

        // Perform the RK 3 update
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;

        // Set the previous step solution in the current buffer position
        // Note that we set the 0 position of the buffer to be equal to the values in step n (not n+1)
        // Additionally, we save in an auxiliary vector the value of the fixed DOFs, which is also taken from the previous time step
        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                double& r_u_0 = it_dof->GetSolutionStepValue(0);
                const double& r_u_1 = it_dof->GetSolutionStepValue(1);
                if (it_dof->IsFixed()) {
                    u_n(i_dof) = r_u_0;
                }
                r_u_0 = r_u_1;
            }
        );

        // Calculate the RK3 intermediate sub steps
        for(unsigned int i=1; i<mButcherTableau.Order(); ++i)
        {
            PerformRungeKuttaIntermediateSubStep(i, u_n, rk_K);
        }
        PerformRungeKuttaLastSubStep(rk_K);

        // Do the final solution update
        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    const MatrixRow<LocalSystemMatrixType> substeps_k = row(rk_K, i_dof);
                    r_u = r_u_old + (dt / mass) * inner_prod(mButcherTableau.GetWeights(), substeps_k);

                    if(i_dof == 0)
                    {
                        std::cout << "K=";
                        for(unsigned int i=0; i<mButcherTableau.Order(); ++i)
                        {
                            std::cout << "  " << substeps_k[i];
                        }
                        std::cout << std::endl;
                    }
                } else {
                    r_u = u_n(i_dof);
                }
            }
        );
    }

    /**
     * @brief Initialize the Runge-Kutta intermediate substep
     * This method is intended to implement all the operations required before each Runge-Kutta intermediate substep
     */
    virtual void InitializeRungeKuttaIntermediateSubStep() {};

    /**
     * @brief Finalize the Runge-Kutta intermediate substep
     * This method is intended to implement all the operations required after each Runge-Kutta intermediate substep
     */
    virtual void FinalizeRungeKuttaIntermediateSubStep() {};

    /**
     * @brief Initialize the Runge-Kutta last substep
     * This method is intended to implement all the operations required before each Runge-Kutta last substep
     */
    virtual void InitializeRungeKuttaLastSubStep() {};

    /**
     * @brief Finalize the Runge-Kutta last substep
     * This method is intended to implement all the operations required after each Runge-Kutta last substep
     */
    virtual void FinalizeRungeKuttaLastSubStep() {};

    /**
     * @brief Performs an intermediate RK4 step
     * This functions performs all the operations required in an intermediate RK4 sub step
     * @param SubStepIndex The sub step index
     * @param SubStepCoefficients The sub step coefficients (these are saved as member variables)
     * @param rFixedDofsValues The vector containing the step n+1 values of the fixed DOFs
     * @param rIntermediateStepResidualVector The vector to store the intermediate sub step residual
     */
    virtual void PerformRungeKuttaIntermediateSubStep(
        const IndexType SubStepIndex,
        const LocalSystemVectorType& rFixedDofsValues,
        LocalSystemMatrixType& rIntermediateStepResidualVectors)
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Get model part and information
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = SubStepIndex;

        // Perform the intermidate sub step update
        InitializeRungeKuttaIntermediateSubStep();
        p_explicit_bs->BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rIntermediateStepResidualVectors(i_dof, SubStepIndex-1) = r_res;
                // Do the DOF update
                double& r_u = it_dof->GetSolutionStepValue(0);
                const double& r_u_old = it_dof->GetSolutionStepValue(1);
                if (!it_dof->IsFixed()) {
                    const double mass = r_lumped_mass_vector(i_dof);
                    const auto coeff = mButcherTableau.GetRKMatrixRow(SubStepIndex);
                    const auto previous_k = row(rIntermediateStepResidualVectors, i_dof);
                    r_u = r_u_old + (dt / mass) * inner_prod(coeff, previous_k);
                } else {
                    const double delta_u = rFixedDofsValues(i_dof) - r_u_old;
                    const auto coeff = mButcherTableau.GetNode(SubStepIndex);
                    r_u = r_u_old + coeff * delta_u;
                }
            }
        );

        FinalizeRungeKuttaIntermediateSubStep();
    }

    /**
     * @brief Performs the last RK4 step
     * This functions performs all the operations required in the last RK4 sub step
     * @param rLastStepResidualVector The vector to store the last sub step residual
     */
    virtual void PerformRungeKuttaLastSubStep(LocalSystemMatrixType& rLastStepResidualVector)
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();

        // Get model part
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();

        // Set the RUNGE_KUTTA_STEP value. This has to be done prior to the InitializeRungeKuttaStep()
        r_process_info.GetValue(RUNGE_KUTTA_STEP) = mButcherTableau.Size();

        // Perform the last sub step residual calculation
        InitializeRungeKuttaLastSubStep();
        p_explicit_bs->BuildRHS(r_model_part);

        IndexPartition<int>(r_dof_set.size()).for_each(
            [&](int i_dof){
                const auto it_dof = r_dof_set.begin() + i_dof;
                // Save current value in the corresponding vector
                const double& r_res = it_dof->GetSolutionStepReactionValue();
                rLastStepResidualVector(i_dof, mButcherTableau.Order() - 1) = r_res;
            }
        );

        FinalizeRungeKuttaLastSubStep();
    }

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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
};

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA  defined */
