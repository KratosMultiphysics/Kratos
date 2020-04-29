//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Carlos Roig and Ruben Zorrilla
//

#ifndef KRATOS_TRILINOS_MIXED_VECTOR_SCALAR_CRITERIA_H
#define	KRATOS_TRILINOS_MIXED_VECTOR_SCALAR_CRITERIA_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Application includes

namespace Kratos
{
///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for mixed vector-scalar problems.
/**
 This class implements a convergence control based on a nodal vector variable and
 a nodal scalar variable. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template< class TSparseSpace, class TDenseSpace >
class TrilinosMixedVectorScalarCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosMixedVectorScalarCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename DofUpdater<TSparseSpace >::UniquePointer DofUpdaterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param VectorRatioTolerance Relative tolerance for vector variable error
     * @param VectorAbsTolerance Absolute tolerance for vector variable error
     * @param ScalarRatioTolerance Relative tolerance for scalar variable error
     * @param ScalarAbsTolerance Absolute tolerance for scalar variable error
     * @param rVectorVariable Reference to the vector variable to check
     * @param rScalarVariable Reference to the scalar variable to check
     */
    TrilinosMixedVectorScalarCriteria(
        TDataType VectorRatioTolerance,
        TDataType VectorAbsTolerance,
        TDataType ScalarRatioTolerance,
        TDataType ScalarAbsTolerance,
        const Variable<array_1d<double,3>>& rVectorVariable,
        const Variable<double>& rScalarVariable)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
        mVectorRatioTolerance(VectorRatioTolerance),
        mVectorAbsTolerance(VectorAbsTolerance),
        mScalarRatioTolerance(ScalarRatioTolerance),
        mScalarAbsTolerance(ScalarAbsTolerance),
        mrVectorVariable(rVectorVariable),
        mrScalarVariable(rScalarVariable)
    {}

    /// Destructor.
    ~TrilinosMixedVectorScalarCriteria() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Compute relative and absoute error.
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b) override
    {
        // If required, initialize the Epetra vector import
        if(!mEpetraImportIsInitialized) {
            InitializeEpetraImport(rDofSet, Dx);
        }

        // Check if we are solving for something
        if (TSparseSpace::Size(Dx) != 0) {
            // Initialize values
            TDataType vector_solution_norm = 0.0;
            TDataType scalar_solution_norm = 0.0;
            TDataType vector_correction_norm = 0.0;
            TDataType scalar_correction_norm = 0.0;

            int n_dofs = rDofSet.size();
            unsigned int n_vector_dof = 0;
            unsigned int n_scalar_dof = 0;
            const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
            const int rank = r_data_comm.Rank();

            // Do the local Dx vector import
            Epetra_Vector local_dx( mpDofImport->TargetMap() );
            int i_err = local_dx.Import(Dx, *mpDofImport, Insert);
            KRATOS_ERROR_IF_NOT(i_err == 0) << "Local Dx import failed!" << std::endl;

            // Loop over Dofs
#pragma omp parallel reduction(+:vector_solution_norm,scalar_solution_norm,vector_correction_norm,scalar_correction_norm,n_vector_dof,n_scalar_dof)
            {
                int dof_id;
                TDataType dof_value;
                TDataType dof_dx;

#pragma omp for
                for (int i = 0; i < n_dofs; i++) {
                    auto it_dof = rDofSet.begin() + i;
                    if (it_dof->IsFree() && it_dof->GetSolutionStepValue(PARTITION_INDEX) == rank) {
                        dof_id = it_dof->EquationId();
                        dof_value = it_dof->GetSolutionStepValue();
                        dof_dx = local_dx[mpDofImport->TargetMap().LID(dof_id)];

                        // Accumulate the values in the corresponding variables
                        // We assume that if the DOF does not correspond to the scalar variable it is a component of the array one
                        const auto& r_dof_var = it_dof->GetVariable();
                        if (r_dof_var == mrScalarVariable) {
                            scalar_solution_norm += dof_value * dof_value;
                            scalar_correction_norm += dof_dx * dof_dx;
                            ++n_scalar_dof;
                        } else {
                            vector_solution_norm += dof_value * dof_value;
                            vector_correction_norm += dof_dx * dof_dx;
                            ++n_vector_dof;
                        }
                        // if ((r_dof_var == VELOCITY_X) || (r_dof_var == VELOCITY_Y) || (r_dof_var == VELOCITY_Z))
                        // {
                        //     vector_solution_norm += dof_value * dof_value;
                        //     vector_correction_norm += dof_dx * dof_dx;
                        //     ++n_vector_dof;
                        // }
                        // else
                        // {
                        //     scalar_solution_norm += dof_value * dof_value;
                        //     scalar_correction_norm += dof_dx * dof_dx;
                        //     ++n_scalar_dof;
                        // }
                    }
                }
            }

            // Accumulate values between ranks
            const unsigned int global_n_vector_dof = r_data_comm.SumAll(n_vector_dof);
            const unsigned int global_n_scalar_dof = r_data_comm.SumAll(n_scalar_dof);
            const TDataType global_vector_solution_norm = r_data_comm.SumAll(vector_solution_norm);
            const TDataType global_scalar_solution_norm = r_data_comm.SumAll(scalar_solution_norm);
            const TDataType global_vector_correction_norm = r_data_comm.SumAll(vector_correction_norm);
            const TDataType global_scalar_correction_norm = r_data_comm.SumAll(scalar_correction_norm);

            // Calculate convergence values
            const double zero_tol = 1.0e-12;
            // if(global_vector_solution_norm < zero_tol)
            //     global_vector_solution_norm = 1.0;
            // if(global_scalar_solution_norm < zero_tol)
            //     global_scalar_solution_norm = 1.0;

            const TDataType vector_ratio = std::sqrt(global_vector_correction_norm/(global_vector_solution_norm < zero_tol ? 1.0 : global_vector_solution_norm));
            const TDataType scalar_ratio = std::sqrt(global_scalar_correction_norm/(global_scalar_solution_norm < zero_tol ? 1.0 : global_scalar_solution_norm));
            const TDataType vector_abs = std::sqrt(global_vector_correction_norm)/ static_cast<TDataType>(global_n_vector_dof);
            const TDataType scalar_abs = std::sqrt(global_scalar_correction_norm)/ static_cast<TDataType>(global_n_scalar_dof);

            KRATOS_INFO_IF("", this->GetEchoLevel() > 0)
                << "CONVERGENCE CHECK:\n"
                << " VELOC.: ratio = " << vector_ratio << "; exp.ratio = " << mVectorRatioTolerance << " abs = " << vector_abs << " exp.abs = " << mVectorAbsTolerance << "\n"
                << " PRESS.: ratio = " << scalar_ratio << "; exp.ratio = " << mScalarRatioTolerance << " abs = " << scalar_abs << " exp.abs = " << mScalarAbsTolerance << std::endl;

            if ((vector_ratio <= mVectorRatioTolerance || vector_abs <= mVectorAbsTolerance) &&
                (scalar_ratio <= mScalarRatioTolerance || scalar_abs <= mScalarAbsTolerance)) {
                KRATOS_INFO_IF("", this->GetEchoLevel() > 0) << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
                return true;
            } else {
                return false;
            }
        } else {
            // Case in which all the DOFs are constrained!
            return true;
        }
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
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const TDataType mVectorRatioTolerance;
    const TDataType mVectorAbsTolerance;

    const TDataType mScalarRatioTolerance;
    const TDataType mScalarAbsTolerance;

    const Variable<array_1d<double,3>>& mrVectorVariable;
    const Variable<double>& mrScalarVariable;

    bool mEpetraImportIsInitialized = false;
    std::unique_ptr<Epetra_Import> mpDofImport = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Initialize the DofUpdater in preparation for a subsequent UpdateDofs call.
     *  The DofUpdater needs to be initialized only if the dofset changes.
     *  If the problem does not require creating/destroying nodes or changing the
     *  mesh graph, it is in general enough to intialize this tool once at the
     *  begining of the problem.
     *  If the dofset only changes under certain conditions (for example because
     *  the domain is remeshed every N iterations), it is enough to call the
     *  Clear method to let this class know that its auxiliary data has to be re-generated
     *  and Initialize will be called as part of the next UpdateDofs call.
     * @param rDofSet rDofSet The list of degrees of freedom.
     * @param rDx rDx The update vector.
     */
    void InitializeEpetraImport(
        const DofsArrayType &rDofSet,
        const TSystemVectorType &rDx)
    {
        int number_of_dofs = rDofSet.size();
        int system_size = TSparseSpace::Size(rDx);
        std::vector<int> index_array(number_of_dofs);

        // Filling the array with the global ids
        unsigned int counter = 0;
        for (typename DofsArrayType::const_iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof) {
            const int id = i_dof->EquationId();
            if (id < system_size) {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator new_end = std::unique(index_array.begin(), index_array.end());
        index_array.resize(new_end - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rDx.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        KRATOS_ERROR_IF(check_size < system_size)
            << "DOF count is not correct. There are less dofs then expected.\n"
            << "Expected number of active dofs: " << system_size << ", DOFs found: " << check_size << std::endl;

        // Defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(), &(*(index_array.begin())), 0, rDx.Comm());

        // Defining the import instance
        std::unique_ptr<Epetra_Import> p_dof_import(new Epetra_Import(dof_update_map, rDx.Map()));
        mpDofImport.swap(p_dof_import);

        mEpetraImportIsInitialized = true;
    }

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
///@} // Kratos classes

///@} // Application group
}

#endif // KRATOS_TRILINOS_MIXED_VECTOR_SCALAR_CRITERIA_H