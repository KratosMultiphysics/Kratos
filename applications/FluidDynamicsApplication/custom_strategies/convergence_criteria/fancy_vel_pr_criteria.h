//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#ifndef KRATOS_FANCY_VEL_PR_CRITERIA_H
#define	KRATOS_FANCY_VEL_PR_CRITERIA_H

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "utilities/bprinter_utility.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#if !defined(_WIN32)
	#include "utilities/color_utilities.h"
#endif
namespace Kratos
{
///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for fluid problems.
/**
 This class implements a convergence control based on nodal velocity and
 pressure values. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template<   class TSparseSpace,
            class TDenseSpace >
class FancyVelPrCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FancyVelPrCriteria );
    
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >  BaseType;

    typedef TSparseSpace                               SparseSpaceType;

    typedef typename BaseType::TDataType                     TDataType;

    typedef typename BaseType::DofsArrayType             DofsArrayType;

    typedef typename BaseType::TSystemMatrixType     TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType     TSystemVectorType;
    
    typedef OpenMPUtils::PartitionVector               PartitionVector;
    
    typedef std::size_t                                        KeyType;
    
    typedef boost::shared_ptr<BprinterUtility> TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param VelRatioTolerance Relative tolerance for velocity error
     * @param VelAbsTolerance Absolute tolerance for velocity error
     * @param PrsRatioTolerance Relative tolerance for presssure error
     * @param PrsAbsTolerance Absolute tolerance for presssure error
     */
    FancyVelPrCriteria(  
        TDataType VelRatioTolerance,
        TDataType VelAbsTolerance,
        TDataType PrsRatioTolerance,
        TDataType PrsAbsTolerance,
        TablePrinterPointerType pTable = nullptr,
        const bool StandaloneTable = true,
        const bool PrintingOutput = false 
        )
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mpTable(pTable),
          mStandaloneTable(StandaloneTable),
          mPrintingOutput(PrintingOutput),
          mTableIsInitialized(false)
    {
        mVelRatioTolerance = VelRatioTolerance;
        mVelAbsTolerance = VelAbsTolerance;

        mPrRatioTolerance = PrsRatioTolerance;
        mPrAbsTolerance = PrsAbsTolerance;
    }

    /// Destructor.
    ~FancyVelPrCriteria() override {}

    ///@}
    ///@name Operators
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
    bool PostCriteria(  ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        const TSystemMatrixType& A,
                        const TSystemVectorType& Dx,
                        const TSystemVectorType& b ) override
    {
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
        {
            if (mpTable != nullptr)
            {
                const unsigned int nl_iteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                mpTable->AddToRow<unsigned int>(nl_iteration);
            }
        }
        
        bool criterion_result;
        
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            // Initialize
            TDataType vel_solution_norm = 0.0;
            TDataType pr_solution_norm = 0.0;
            TDataType vel_increase_norm = 0.0;
            TDataType pr_increase_norm = 0.0;
            unsigned int vel_dof_num(0),pr_dof_num(0);

            // Set a partition for OpenMP
            const int num_dofs = rDofSet.size();
            PartitionVector DofPartition;
            const int num_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(num_dofs,num_threads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:vel_solution_norm,pr_solution_norm,vel_increase_norm,pr_increase_norm,vel_dof_num,pr_dof_num)
            {
                const int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t dof_id;
                TDataType dof_value;
                TDataType dof_incr;

                for (typename DofsArrayType::iterator it_dof = DofBegin; it_dof != DofEnd; ++it_dof)
                {
                    if (it_dof->IsFree())
                    {
                        dof_id = it_dof->EquationId();
                        dof_value = it_dof->GetSolutionStepValue(0);
                        dof_incr = Dx[dof_id];

                        KeyType curr_var = it_dof->GetVariable().Key();
                        if ((curr_var == VELOCITY_X) || (curr_var == VELOCITY_Y) || (curr_var == VELOCITY_Z))
                        {
                            vel_solution_norm += dof_value * dof_value;
                            vel_increase_norm += dof_incr * dof_incr;
                            ++vel_dof_num;
                        }
                        else
                        {
                            pr_solution_norm += dof_value * dof_value;
                            pr_increase_norm += dof_incr * dof_incr;
                            ++pr_dof_num;
                        }
                    }
                }
            }

            if(vel_solution_norm == 0.0) vel_solution_norm = 1.0;
            if(pr_solution_norm == 0.0) pr_solution_norm = 1.0;

            TDataType vel_ratio = sqrt(vel_increase_norm/vel_solution_norm);
            TDataType pr_ratio = sqrt(pr_increase_norm/pr_solution_norm);

            TDataType vel_abs = sqrt(vel_increase_norm)/ static_cast<TDataType>(vel_dof_num);
            TDataType pr_abs = sqrt(pr_increase_norm)/ static_cast<TDataType>(pr_dof_num);

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                if (mpTable != nullptr)
                {
                    std::cout.precision(4);
                    auto& Table = mpTable->GetTable();
                    Table  << vel_ratio  << mVelRatioTolerance  << vel_abs  << mVelAbsTolerance << pr_ratio  << mPrRatioTolerance  << pr_abs  << mPrAbsTolerance;
                }
                else
                {
                    std::cout.precision(4);
                    if (mPrintingOutput == false)
                    {
                    #if !defined(_WIN32)
                        std::cout << BOLDFONT("VELOCITY AND PRESSURE CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                        std::cout << BOLDFONT("\tVELOCITY: RATIO = ") << vel_ratio << BOLDFONT(" EXP.RATIO = ") << mVelRatioTolerance << BOLDFONT(" ABS = ") << vel_abs << BOLDFONT(" EXP.ABS = ") << mVelAbsTolerance << std::endl;
                        std::cout << BOLDFONT("\tPRESSURE: RATIO = ") << pr_ratio << BOLDFONT(" EXP.RATIO = ") << mPrRatioTolerance << BOLDFONT(" ABS = ") << pr_abs << BOLDFONT(" EXP.ABS = ") << mPrAbsTolerance << std::endl;
                    #else
                        std::cout << "VELOCITY AND PRESSURE CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                        std::cout << "\tVELOCITY: RATIO = " << vel_ratio << " EXP.RATIO = " << mVelRatioTolerance << " ABS = " << vel_abs << " EXP.ABS = " << mVelAbsTolerance << std::endl;
                        std::cout << "\tPRESSURE: RATIO = " << pr_ratio << " EXP.RATIO = " << mPrRatioTolerance << " ABS = " << pr_abs << " EXP.ABS = " << mPrAbsTolerance << std::endl;
                    #endif
                    }
                    else
                    {
                        std::cout << "VELOCITY AND PRESSURE CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tNL ITERATION: " << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] << std::endl;
                        std::cout << "\tVELOCITY: RATIO = " << vel_ratio << " EXP.RATIO = " << mVelRatioTolerance << " ABS = " << vel_abs << " EXP.ABS = " << mVelAbsTolerance << std::endl;
                        std::cout << "\tPRESSURE: RATIO = " << pr_ratio << " EXP.RATIO = " << mPrRatioTolerance << " ABS = " << pr_abs << " EXP.ABS = " << mPrAbsTolerance << std::endl;
                    }
                }
            }
            
            if (    (vel_ratio <= mVelRatioTolerance || vel_abs <= mVelAbsTolerance) &&
                    (pr_ratio <= mPrRatioTolerance || pr_abs <= mPrAbsTolerance) )
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    if (mpTable != nullptr)
                    {
                        auto& Table = mpTable->GetTable();
                        if (mPrintingOutput == false)
                        {
                        #if !defined(_WIN32)
                            Table << BOLDFONT(FGRN("       Achieved"));
                        #else
                            Table << "Achieved";
                        #endif
                        }
                        else
                        {
                            Table << "Achieved";
                        }
                    }
                    else
                    {
                        if (mPrintingOutput == false)
                        {
                        #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tVelocity and pressure") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        #else
                            std::cout << "\tVelocity and pressure convergence is achieved" << std::endl;
                        #endif
                        }
                        else
                        {
                            std::cout << "\tVelocity and pressure convergence is achieved" << std::endl;
                        }
                    }
                }
                
                criterion_result = true;
            }
            else
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    if (mpTable != nullptr)
                    {
                        auto& Table = mpTable->GetTable();
                        if (mPrintingOutput == false)
                        {
                        #if !defined(_WIN32)
                            Table << BOLDFONT(FRED("   Not achieved"));
                        #else
                            Table << "Not achieved";
                        #endif
                        }
                        else
                        {
                            Table << "Not achieved";
                        }
                    }
                }
                
                criterion_result = false;
            }
        }
        else //in this case all the displacements are imposed!
        {
            criterion_result = true;
        }
        
        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
        {
            if (mpTable != nullptr)
            {
                mpTable->PrintFooter();
            }
        }
        
        return criterion_result;
    }

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize( ModelPart& rModelPart	) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
        
        auto& Table = mpTable->GetTable();
        if (mStandaloneTable == true) Table.AddColumn("ITER", 4);
        Table.AddColumn("VEL. RAT", 10);
        Table.AddColumn("EXP. RAT", 10);
        Table.AddColumn("VEL. ABS", 10);
        Table.AddColumn("EXP. ABS", 10);
        Table.AddColumn("PRE. RAT", 10);
        Table.AddColumn("EXP. RAT", 10);
        Table.AddColumn("PRE. ABS", 10);
        Table.AddColumn("EXP. ABS", 10);
        Table.AddColumn("CONVERGENCE", 15);
        mTableIsInitialized = true;
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void InitializeSolutionStep(    
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b 
        ) override
    {
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0 && mStandaloneTable == true)
        {
            std::cout.precision(4);
            if (mPrintingOutput == false)
            {
            #if !defined(_WIN32)
                std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #else
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #endif
            }
            else
            {
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            }
                
            if (mpTable != nullptr)
            {
                mpTable->PrintHeader();
            }
        }
    }

    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeSolutionStep(  
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {}

    ///@} // Operations

private:

    TDataType mVelRatioTolerance;    // The ratio of tolerance of the velocity
    
    TDataType mVelAbsTolerance;      // The absolute tolerance of the velocity

    TDataType mPrRatioTolerance;     // The ratio of tolerance of the pressure
    
    TDataType mPrAbsTolerance;       // The absolute tolerance of the pressure
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    
    bool mStandaloneTable;           // If the table is not appended to any other table
    
    bool mPrintingOutput;            // If the colors and bold are printed
    
    bool mTableIsInitialized;        // If the table is already initialized 
};

///@} // Kratos classes

///@} // Application group
}

#endif	/* _VEL_PR_CRITERIA_H */
