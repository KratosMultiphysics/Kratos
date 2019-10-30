//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <unordered_map>
#include <cmath>

/* Project includes */
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/openmp_utils.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "Epetra_FEVector.h"
#include "trilinos_space.h"
#endif

// Include base h
#include "generic_convergence_criteria.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::GenericConvergenceCriteria(
    double rRatioTolerance, double rAbsTolerance)
    : ConvergenceCriteria<TSparseSpace, TDenseSpace>(),
      mRatioTolerance(rRatioTolerance),
      mAbsTolerance(rAbsTolerance)
{
    mSolutionVariables = "";
}

template <class TSparseSpace, class TDenseSpace>
GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::~GenericConvergenceCriteria()
{
}

template <class TSparseSpace, class TDenseSpace>
bool GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::PostCriteria(
    ModelPart& rModelPart,
    DofsArrayType& rDofSet,
    const TSystemMatrixType& A,
    const TSystemVectorType& Dx,
    const TSystemVectorType& b)
{
    if (SparseSpaceType::Size(Dx) != 0) // if we are solving for something
    {
        double solution_norm, increase_norm, dof_size;

        this->CalculateConvergenceCheckNorms(
            solution_norm, increase_norm, dof_size, rModelPart, rDofSet, A, Dx, b);

        if (solution_norm == 0.0)
            solution_norm = 1.0;

        const double ratio = increase_norm / solution_norm;
        const double ratio_abs = increase_norm / dof_size;

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int iteration = r_current_process_info[NL_ITERATION_NUMBER];

        if (this->GetEchoLevel() > 0)
        {
            std::stringstream msg;
            msg << "[" << iteration << "] CONVERGENCE CHECK: ";
            msg << mSolutionVariables;
            msg << ": ratio = " << std::scientific << ratio
                << "; exp.ratio = " << std::scientific << mRatioTolerance;
            msg << ": abs = " << std::scientific << ratio_abs
                << "; exp.abs = " << std::scientific << mAbsTolerance << std::endl;
            KRATOS_INFO(this->Info()) << msg.str();
        }

        if ((std::abs(ratio) > mRatioTolerance) && (std::abs(ratio_abs) > mAbsTolerance))
            return false;

        KRATOS_INFO_IF(this->Info(), this->GetEchoLevel() > 0)
            << "CONVERGENCE CHECK: " << mSolutionVariables
            << ": *** CONVERGENCE IS ACHIEVED ***\n";

        return true;
    }
    else // in this case all the displacements are imposed!
    {
        return true;
    }
}

template <class TSparseSpace, class TDenseSpace>
void GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::Initialize(ModelPart& rModelPart)
{
    BaseType::mConvergenceCriteriaIsInitialized = true;
}

template <class TSparseSpace, class TDenseSpace>
void GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::InitializeSolutionStep(
    ModelPart& rModelPart,
    DofsArrayType& rDofSet,
    const TSystemMatrixType& A,
    const TSystemVectorType& Dx,
    const TSystemVectorType& b)
{
    if (mSolutionVariables == "")
    {
        const int number_of_dofs = rDofSet.size();

        // create the variable list
        std::unordered_map<std::string, int> variables_map;
        for (int i_dof = 0; i_dof < number_of_dofs; ++i_dof)
        {
            variables_map.insert(std::make_pair(
                (rDofSet.begin() + i_dof)->GetVariable().Name(), 1));
        }

        for (const auto& pair : variables_map)
        {
            mSolutionVariables += pair.first + ", ";
        }

        // remove last two characters
        mSolutionVariables.pop_back();
        mSolutionVariables.pop_back();
    }
}

template <class TSparseSpace, class TDenseSpace>
std::string GenericConvergenceCriteria<TSparseSpace, TDenseSpace>::Info() const
{
    return "GenericConvergenceCriteria";
}

template <>
void GenericConvergenceCriteria<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>::CalculateConvergenceCheckNorms(
    double& rSolutionNorm,
    double& rIncreaseNorm,
    double& rDofSize,
    ModelPart& rModelPart,
    DofsArrayType& rDofSet,
    const TSystemMatrixType& A,
    const TSystemVectorType& Dx,
    const TSystemVectorType& b)
{
    int NumDofs = rDofSet.size();

    double solution_norm{0.0}, increase_norm{0.0};
    int dof_num{0};

    // Set a partition for OpenMP
    PartitionVector DofPartition;
    int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::DivideInPartitions(NumDofs, NumThreads, DofPartition);

    // Loop over Dofs
#pragma omp parallel reduction(+ : solution_norm, increase_norm, dof_num)
    {
        int k = OpenMPUtils::ThisThread();
        typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
        typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k + 1];

        std::size_t DofId;
        TDataType DofValue;
        TDataType DofIncr;

        for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
        {
            if (itDof->IsFree())
            {
                DofId = itDof->EquationId();
                DofValue = itDof->GetSolutionStepValue(0);
                DofIncr = Dx[DofId];

                solution_norm += DofValue * DofValue;
                increase_norm += DofIncr * DofIncr;
                dof_num += 1;
            }
        }
    }

    rSolutionNorm = std::sqrt(solution_norm);
    rIncreaseNorm = std::sqrt(increase_norm);
    rDofSize = static_cast<double>(dof_num);
}

// template instantiations

template class GenericConvergenceCriteria<UblasSpace<double, CompressedMatrix, Vector>,
                                          UblasSpace<double, Matrix, Vector>>;

#ifdef KRATOS_USING_MPI
template <>
void GenericConvergenceCriteria<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>::CalculateConvergenceCheckNorms(
    double& rSolutionNorm,
    double& rIncreaseNorm,
    double& rDofSize,
    ModelPart& rModelPart,
    DofsArrayType& rDofSet,
    const TSystemMatrixType& A,
    const TSystemVectorType& Dx,
    const TSystemVectorType& b)
{
    int NumDofs = rDofSet.size();

    double solution_norm{0.0}, increase_norm{0.0};
    int dof_num{0};

    // Set a partition for OpenMP
    PartitionVector DofPartition;
    int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::DivideInPartitions(NumDofs, NumThreads, DofPartition);

    const Communicator& r_communicator = rModelPart.GetCommunicator();
    const int my_pid = r_communicator.MyPID();

    // Loop over Dofs
#pragma omp parallel reduction(+ : solution_norm, increase_norm, dof_num)
    {
        int k = OpenMPUtils::ThisThread();
        typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
        typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k + 1];

        std::size_t DofId;
        TDataType DofValue;
        TDataType DofIncr;

        for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
        {
            if (itDof->IsFree() && itDof->GetSolutionStepValue(PARTITION_INDEX) == my_pid)
            {
                DofId = itDof->EquationId();
                DofValue = itDof->GetSolutionStepValue(0);
                DofIncr = SparseSpaceType::GetValue(Dx, DofId);

                solution_norm += DofValue * DofValue;
                increase_norm += DofIncr * DofIncr;
                dof_num += 1;
            }
        }
    }

    std::vector<double> residual_norms = {increase_norm, solution_norm,
                                          static_cast<double>(dof_num)};
    const std::vector<double>& total_residual_norms =
        r_communicator.GetDataCommunicator().SumAll(residual_norms);

    rSolutionNorm = std::sqrt(total_residual_norms[1]);
    rIncreaseNorm = std::sqrt(total_residual_norms[0]);
    rDofSize = total_residual_norms[2];
}

// template instantiations

template class GenericConvergenceCriteria<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
                                          UblasSpace<double, Matrix, Vector>>;

#endif

} // namespace Kratos
