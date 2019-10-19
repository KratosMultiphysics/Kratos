//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_FEMDEM_RESIDUAL_CRITERIA )
#define  KRATOS_FEMDEM_RESIDUAL_CRITERIA

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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
 * @class FemDemResidualCriteria
 * @ingroup KratosCore
 * @brief This is a convergence criteria that employes the residual as criteria
 * @details The reactions from the RHS are not computed in the residual
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace
         >
class FemDemResidualCriteria
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FemDemResidualCriteria );

    /// The definition of the base ConvergenceCriteria
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    /// The data type
    typedef typename BaseType::TDataType TDataType;

    /// The dofs array type
    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// The sparse matrix type
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    /// The dense vector type
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /// Definition of the IndexType
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    //* Constructor.
    explicit FemDemResidualCriteria(Kratos::Parameters Settings)
        : BaseType()
    {
        if (Settings.Has("residual_absolute_tolerance")) {
            mAlwaysConvergedNorm = Settings["residual_absolute_tolerance"].GetDouble();
        } else if (Settings.Has("absolute_tolerance")) {
            mAlwaysConvergedNorm = Settings["absolute_tolerance"].GetDouble();
        } else {
            KRATOS_WARNING("FemDemResidualCriteria") << "residual_absolute_tolerance or absolute_tolerance nor defined on settings. Using default 1.0e-9" << std::endl;
            mAlwaysConvergedNorm = 1.0e-9;
        }
        if (Settings.Has("residual_relative_tolerance")) {
            mRatioTolerance = Settings["residual_relative_tolerance"].GetDouble();
        } else if (Settings.Has("relative_tolerance")) {
            mRatioTolerance = Settings["relative_tolerance"].GetDouble();
        } else {
            KRATOS_WARNING("FemDemResidualCriteria") << "residual_relative_tolerance or relative_tolerance nor defined on settings. Using default 1.0e-4" << std::endl;
            mRatioTolerance = 1.0e-4;
        }
        if (Settings.Has("max_iteration")) {
            mMaxIterations = Settings["max_iteration"].GetInt();
        }

        this->mActualizeRHSIsNeeded = true;
    }

    //* Constructor.
    explicit FemDemResidualCriteria(
        TDataType NewRatioTolerance,
        TDataType AlwaysConvergedNorm)
        : BaseType(),
          mRatioTolerance(NewRatioTolerance),
          mAlwaysConvergedNorm(AlwaysConvergedNorm)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    //* Copy constructor.
    explicit FemDemResidualCriteria( FemDemResidualCriteria const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
      ,mReferenceDispNorm(rOther.mReferenceDispNorm)
    {
        this->mActualizeRHSIsNeeded = true;
    }

    //* Destructor.
    ~FemDemResidualCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Criterias that need to be called after getting the solution
     * @details Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        if (rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] == 1) {
            KRATOS_INFO("") << "___________________________________________________________________" << std::endl;
            KRATOS_INFO("") << "|    ITER     |     RATIO      |    ABS_NORM    |    CONVERGED    |" << std::endl;
        }
        const SizeType size_b = TSparseSpace::Size(rb);
        if (size_b != 0) { //if we are solving for something

            SizeType size_residual;
            CalculateResidualNorm(rModelPart, mCurrentResidualNorm, size_residual, rDofSet, rb);

            TDataType ratio = 0.0;
            if(mInitialResidualNorm < std::numeric_limits<TDataType>::epsilon()) {
                ratio = 0.0;
            } else {
                ratio = mCurrentResidualNorm/mInitialResidualNorm;
            }

            const TDataType float_size_residual = static_cast<TDataType>(size_residual);
            const TDataType absolute_norm = (mCurrentResidualNorm / float_size_residual);

            rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
            rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

            if (ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm) { // Converged
                if (rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] < 10) {
                    std::cout <<"|      " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << "      |  " 
                        << std::scientific << ratio << "  |  " 
                        << absolute_norm << "  |" << "      TRUE       |"<< std::endl;                
                } else {
                    std::cout <<"|      " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << "     |  "
                        << std::scientific << ratio << "  |  " 
                        << absolute_norm << "  |" << "      TRUE       |"<< std::endl;  
                }
                return true;
            } else {
                if (rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] < 10) {
                    std::cout <<"|      " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << "      |  " 
                        << std::scientific << ratio << "  |  " 
                        << absolute_norm << "  |" << "      FALSE      |"<< std::endl;                
                } else {
                    std::cout <<"|      " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << "     |  " 
                        << std::scientific << ratio << "  |  " 
                        << absolute_norm << "  |" << "      FALSE      |"<< std::endl;  
                }
                if (rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] == mMaxIterations) {
                    KRATOS_INFO("") << " ATTENTION! SOLUTION STEP NOT CONVERGED AFTER " <<  mMaxIterations << "ITERATIONS" << std::endl;
                }
                return false;
            }
        } else {
            return true;
        }
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);

        // Filling mActiveDofs when MPC exist
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            mActiveDofs.resize(rDofSet.size());

            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mActiveDofs.size()); ++i) {
                mActiveDofs[i] = true;
            }

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(rDofSet.size()); ++i) {
                const auto it_dof = rDofSet.begin() + i;
                if (it_dof->IsFixed()) {
                    mActiveDofs[it_dof->EquationId()] = false;
                }
            }

            for (const auto& r_mpc : rModelPart.MasterSlaveConstraints()) {
                for (const auto& r_dof : r_mpc.GetMasterDofsVector()) {
                    mActiveDofs[r_dof->EquationId()] = false;
                }
                for (const auto& r_dof : r_mpc.GetSlaveDofsVector()) {
                    mActiveDofs[r_dof->EquationId()] = false;
                }
            }
        }

        SizeType size_residual;
        CalculateResidualNorm(rModelPart, mInitialResidualNorm, size_residual, rDofSet, rb);
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual + reactions)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);
        KRATOS_INFO("") << "|_____________|________________|________________|_________________|" << std::endl;
        KRATOS_INFO("") << "" << std::endl;
    }

    ///@}
    ///@name Operations
    ///@{

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
        return "FemDemResidualCriteria";
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

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void CalculateResidualNorm(
        ModelPart& rModelPart,
        TDataType& rResidualSolutionNorm,
        SizeType& rDofNum,
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        TDataType residual_solution_norm = TDataType();
        SizeType dof_num = 0;

        // Auxiliar values
        TDataType residual_dof_value = 0.0;
        const auto it_dof_begin = rDofSet.begin();
        const int number_of_dof = static_cast<int>(rDofSet.size());

        // Loop over Dofs
        if (rModelPart.NumberOfMasterSlaveConstraints() > 0) {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:residual_solution_norm, dof_num)
            for (int i = 0; i < number_of_dof; i++) {
                auto it_dof = it_dof_begin + i;

                const IndexType dof_id = it_dof->EquationId();

                if (mActiveDofs[dof_id]) {
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);
                    residual_solution_norm += std::pow(residual_dof_value, 2);
                    dof_num++;
                }
            }
        } else {
            #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:residual_solution_norm, dof_num)
            for (int i = 0; i < number_of_dof; i++) {
                auto it_dof = it_dof_begin + i;

                if (!it_dof->IsFixed()) {
                    const IndexType dof_id = it_dof->EquationId();
                    residual_dof_value = TSparseSpace::GetValue(rb,dof_id);
                    residual_solution_norm += std::pow(residual_dof_value, 2);
                    dof_num++;
                }
            }
        }

        rDofNum = dof_num;
        rResidualSolutionNorm = std::sqrt(residual_solution_norm);
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

    TDataType mRatioTolerance;      /// The ratio threshold for the norm of the residual

    TDataType mInitialResidualNorm; /// The reference norm of the residual

    TDataType mCurrentResidualNorm; /// The current norm of the residual

    TDataType mAlwaysConvergedNorm; /// The absolute value threshold for the norm of the residual

    TDataType mReferenceDispNorm;   /// The norm at the beginning of the iterations

    std::vector<bool> mActiveDofs;  /// This vector contains the dofs that are active

    int mMaxIterations;

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

}; // Class ClassName

///@}

///@name Type Definitions
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_NEW_DISPLACEMENT_CRITERIA  defined
