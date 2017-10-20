// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "custom_utilities/contact_utilities.h"
#include "utilities/table_stream_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "utilities/color_utilities.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
///@{
    
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

/** @brief Custom convergence criteria for the mortar condition 
 */
template<class TSparseSpace, class TDenseSpace>
class ALMFrictionlessMortarConvergenceCriteria : public virtual  BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ALMFrictionlessMortarConvergenceCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaBaseType;
    
    typedef BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >          BaseType;

    typedef TSparseSpace                                                 SparseSpaceType;

    typedef typename BaseType::TDataType                                       TDataType;

    typedef typename BaseType::DofsArrayType                               DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;
    
    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;
    
    typedef ModelPart::NodesContainerType                                 NodesArrayType;
    
    typedef boost::shared_ptr<TableStreamUtility>                TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ALMFrictionlessMortarConvergenceCriteria(
        double Tolerance = std::numeric_limits<double>::epsilon(),
        TablePrinterPointerType pTable = nullptr,
        const bool PrintingOutput = false 
        ) : BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >(),
        mTolerance(Tolerance),
        mpTable(pTable),
        mPrintingOutput(PrintingOutput),
        mTableIsInitialized(false)
    {
    }

    ///Copy constructor 
    ALMFrictionlessMortarConvergenceCriteria( ALMFrictionlessMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
      ,mpTable(rOther.mpTable)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    /// Destructor
    ~ALMFrictionlessMortarConvergenceCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{
    
    /**
     * Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
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
        const TSystemVectorType& b
        ) override
    {
        // We call the base class
        BaseType::PostCriteria(rModelPart, rDofSet, A, Dx, b);
        
        // Defining the convergence
        unsigned int is_converged = 0;
        
//         const double& epsilon = rModelPart.GetProcessInfo()[INITIAL_PENALTY]; 
        const double& scale_factor = rModelPart.GetProcessInfo()[SCALE_FACTOR]; 
        
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        #pragma omp parallel for 
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            
            const double& epsilon = it_node->GetValue(INITIAL_PENALTY);
            
            // Check if the node is slave
            bool node_is_slave = true;
            if ((it_node)->IsDefined(SLAVE))
            {
                node_is_slave = (it_node)->Is(SLAVE);
            }
            
            if (node_is_slave == true)
            {
                const double augmented_normal_pressure = scale_factor * (it_node)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) + epsilon * (it_node)->FastGetSolutionStepValue(WEIGHTED_GAP);     
                    
                (it_node)->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (augmented_normal_pressure < mTolerance * scale_factor) // NOTE: This could be conflictive (< or <=)
                {
                    if ((it_node)->Is(ACTIVE) == false )
                    {
                        (it_node)->Set(ACTIVE, true);
                        #pragma omp atomic
                        is_converged += 1;
                    }
                }
                else
                {
                    if ((it_node)->Is(ACTIVE) == true )
                    {
                        (it_node)->Set(ACTIVE, false);
                        #pragma omp atomic
                        is_converged += 1;
                    }
                }
            }
        }
        
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                auto& table = mpTable->GetTable();
                if (is_converged == 0)
                {
                    if (mPrintingOutput == false)
                    {
                        table << BOLDFONT(FGRN("       Achieved"));
                    }
                    else
                    {
                        table << "Achieved";
                    }
                }
                else
                {
                    if (mPrintingOutput == false)
                    {
                        table << BOLDFONT(FRED("   Not achieved"));
                    }
                    else
                    {
                        table << "Not achieved";
                    }
                }
            }
            else
            {
                if (is_converged == 0)
                {
                    if (mPrintingOutput == false)
                    {
                        std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    }
                }
                else
                {
                    if (mPrintingOutput == false)
                    {
                        std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    }
                    else
                    {
                        std::cout << "\tActive set convergence is not achieved" << std::endl;
                    }
                }
            }
        }
        
        return (is_converged == 0 ? true : false);
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        ConvergenceCriteriaBaseType::mConvergenceCriteriaIsInitialized = true;
        
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            auto& table = mpTable->GetTable();
            table.AddColumn("ACTIVE SET CONV", 15);
            mTableIsInitialized = true;
        }
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
        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") );  
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
    {
    }

    /**
     * This function initializes the non linear iteration
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
    }
    
    /**
     * This function finalizes the non linear iteration
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
    }
    
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    
    double mTolerance;               // Tolerance considered in contact check
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    bool mPrintingOutput;            // If the colors and bold are printed
    bool mTableIsInitialized;        // If the table is already initialized
    
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
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class ALMFrictionlessMortarConvergenceCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H  defined */

