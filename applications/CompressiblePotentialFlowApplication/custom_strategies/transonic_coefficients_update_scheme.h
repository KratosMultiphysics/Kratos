//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

#if !defined(KRATOS_TRANSONIC_COEFFICIENTS_UPDATE_SCHEME_H )
#define  KRATOS_TRANSONIC_COEFFICIENTS_UPDATE_SCHEME_H

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/entities_utilities.h"
#include "custom_processes/compute_nodal_value_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

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

/**
 * @class TransonicCoefficientsUpdateScheme
 * @ingroup CompressiblePotentialFlowApplication
 * @brief This class is a reimplementation of residualBasedIncrementalUpdateStaticScheme that is 
 * used in a transonic potential simulation to update both the CRITICAL_MACH and 
 * UPWIND_FACTOR_CONSTANT values, between nonlinear iterations, to a user-defined value when 
 * a user-defined residual tolerance is reached. This update is done only once.
 * @details The only operation done in this  scheme is the update of the database and values, no predict is done
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @see Scheme
 */
template<class TSparseSpace,class TDenseSpace>
class TransonicCoefficientsUpdateScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransonicCoefficientsUpdateScheme
    KRATOS_CLASS_POINTER_DEFINITION(TransonicCoefficientsUpdateScheme);

    typedef Scheme<TSparseSpace,TDenseSpace>                                 BaseSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;
    typedef TransonicCoefficientsUpdateScheme<TSparseSpace, TDenseSpace> ClassType;

    /// DoF array type definition
    typedef typename BaseType::DofsArrayType                                  DofsArrayType;

    /// Data type definition
    typedef typename BaseType::TDataType                                          TDataType;
    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType                          TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType                          TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemVectorType                  LocalSystemVectorType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemMatrixType                  LocalSystemMatrixType;

    /// Elements containers definition
    typedef ModelPart::ElementsContainerType                              ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType                          ConditionsArrayType;

    /// The definition of the vector containing the equation ids
    typedef Element::EquationIdVectorType                              EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Default onstructor.
    */
    explicit TransonicCoefficientsUpdateScheme()
        : BaseType()
    {  
    }

    /**
     * @brief Constructor. The pseudo static scheme (parameters)
     * @param ThisParameters Dummy parameters
     */
    explicit TransonicCoefficientsUpdateScheme(Parameters ThisParameters) : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }


    /** Destructor.
    */
    ~TransonicCoefficientsUpdateScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseSchemeType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief This is the place to initialize the Scheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        rModelPart.GetProcessInfo()[CRITICAL_MACH]          = rCriticalMach;
        rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = rUpwindFactorConstant;
        rModelPart.GetProcessInfo()[MACH_LIMIT]             = std::sqrt(rMachNumberSquaredLimit);

        BaseType::SetSchemeIsInitialized(true);
        KRATOS_CATCH("")
    }

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        // Update the upwind factor constant and critical mach
        if (rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] < rUpdateTransonicTolerance &&
            rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] > 1                       &&
            updated_values == false){
            
            double mach_1   = 0.0;
            double mach_2   = 0.0;

            if (rUpdateCriticalMach < 0.0)
            {
                // Compute nonal MACH values
                const std::vector<std::string> variable_array = {"MACH"};
                ComputeNodalValueProcess ComputeNodalValueProcess(rModelPart, variable_array);
                ComputeNodalValueProcess.Execute();

                // Select only the body model part
                auto& r_model = rModelPart.GetRootModelPart();
                auto& r_model_part = r_model.GetSubModelPart(rModelPartName);
                const int n_nodes = r_model_part.NumberOfNodes();

                // Search the biggest MACH value
                IndexPartition<std::size_t>(n_nodes).for_each([&](std::size_t index)
                {
                    const auto it_node = r_model_part.NodesBegin() + index;
                    const auto r_mach = it_node->GetValue(MACH);
                    if (mach_1 < r_mach) mach_1 = r_mach;
                });
                
                if (mach_1 > 0.99)
                {
                    mach_2 = std::sqrt(((1.4-1)*std::pow(mach_1,2)+2)/(2*1.4*std::pow(mach_1,2)-(1.4-1))) * 1.2;
                    if (mach_2 > 0.99) mach_2 = 0.99;    
                }
                
                rModelPart.GetProcessInfo()[CRITICAL_MACH] = mach_2;
            } else
            {
                rModelPart.GetProcessInfo()[CRITICAL_MACH] = rUpdateCriticalMach;
            }

            if (rUpdateUpwindFactorConstant > 0.0)
            {
                rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = rUpdateUpwindFactorConstant;
            }

            updated_values = true;
        }

        // Initialize non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "Please_entry_the reference_model_part",
            "critical_mach"                 : 0.92,
            "upwind_factor_constant"        : 2.0,
            "update_critical_mach"          : -1.0,
            "update_upwind_factor_constant" : -1.0,
            "update_transonic_tolerance"    : 1e-3,
            "mach_number_squared_limit"     : 3.0
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        rModelPartName              = ThisParameters["model_part_name"].GetString();
        rCriticalMach               = ThisParameters["critical_mach"].GetDouble();
        rUpdateCriticalMach         = ThisParameters["update_critical_mach"].GetDouble();
        rUpwindFactorConstant       = ThisParameters["upwind_factor_constant"].GetDouble();
        rUpdateUpwindFactorConstant = ThisParameters["update_upwind_factor_constant"].GetDouble();
        rUpdateTransonicTolerance   = ThisParameters["update_transonic_tolerance"].GetDouble();
        rMachNumberSquaredLimit     = ThisParameters["mach_number_squared_limit"].GetDouble();
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "transonic_coefficients_update_scheme";
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "TransonicCoefficientsUpdateScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    std::string rModelPartName;
    bool updated_values = false;
    int rEchoLevel                     =   0;
    double rCriticalMach               = 0.0;
    double rUpdateCriticalMach         = 0.0;
    double rUpwindFactorConstant       = 0.0;
    double rUpdateUpwindFactorConstant = 0.0;
    double rUpdateTransonicTolerance   = 0.0;
    double rMachNumberSquaredLimit     = 0.0;

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

}; // Class TransonicCoefficientsUpdateScheme
}  // namespace Kratos

#endif /* KRATOS_TRANSONIC_COEFFICIENTS_UPDATE_SCHEME_H  defined */
