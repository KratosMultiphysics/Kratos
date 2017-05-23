// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/bprinter_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "custom_utilities/color_utilities.h"

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
    
    typedef boost::shared_ptr<BprinterUtility>                   TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    ALMFrictionlessMortarConvergenceCriteria(TablePrinterPointerType pTable = nullptr)
        : BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >(),
        mpTable(pTable)
    {
    }

    ///Copy constructor 
    ALMFrictionlessMortarConvergenceCriteria( ALMFrictionlessMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    virtual ~ALMFrictionlessMortarConvergenceCriteria() {}

    ///@}
    ///@name Operators
    ///@{
        
    /**
     * Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        // We update the normals if necessary
        if (rModelPart.GetProcessInfo()[CONSIDER_NORMAL_VARIATION] == true)
        {
            // Update normal of the conditions
            ContactUtilities::ComputeNodesMeanNormalModelPart( rModelPart.GetSubModelPart("Contact") ); 
        }
        
        return true;
    }
    
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
        BaseType::CalculateContactReactions(rModelPart, rDofSet, b);
        
        // We define the tolerance
        const double Tolerance = 1.0e-6;
//         const double Tolerance = std::numeric_limits<double>::epsilon();
        
        // Defining the convergence
        bool IsConverged = true;
        
        const double Epsilon = rModelPart.GetProcessInfo()[PENALTY_PARAMETER]; 
        const double ScaleFactor = rModelPart.GetProcessInfo()[SCALE_FACTOR]; 
        
        NodesArrayType& pNodes = rModelPart.GetSubModelPart("Contact").Nodes();
        auto numNodes = pNodes.end() - pNodes.begin();

        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNodes.begin() + i;
            
            // Check if the node is slave
            bool NodeIsSlave = true;
            if ((itNode)->IsDefined(SLAVE))
            {
                NodeIsSlave = (itNode)->Is(SLAVE);
            }
            
            if (NodeIsSlave == true)
            {
                const double AugmentedNormalPressure = ScaleFactor * (itNode)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) + Epsilon * (itNode)->FastGetSolutionStepValue(WEIGHTED_GAP);     
                    
                (itNode)->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, AugmentedNormalPressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (AugmentedNormalPressure < - Tolerance) // NOTE: This could be conflictive (< or <=)
                {
                    if ((itNode)->Is(ACTIVE) == false )
                    {
                        (itNode)->Set(ACTIVE, true);
                        IsConverged = false;
                    }
                }
                else
                {
                    (itNode)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = 0.0; // NOTE: To clear the value (can affect future iterations)
                    
                    if ((itNode)->Is(ACTIVE) == true )
                    {
                        (itNode)->Set(ACTIVE, false);
                        IsConverged = false;
                    }
                }
            }
        }
        
        // We update the pairs if necessary
        if (rModelPart.GetProcessInfo()[CONSIDER_PAIR_VARIATION] == true)
        {
            ConditionsArrayType& pConditions = rModelPart.GetSubModelPart("Contact").Conditions();
            auto numConditions = pConditions.end() - pConditions.begin();

            #pragma omp parallel for 
            for(unsigned int i = 0; i < numConditions; i++) 
            {
                auto itCond = pConditions.begin() + i;
                if ( (itCond)->Is(ACTIVE) == true )
                {
                    bool DeactivateCondition = true;
                    
                    for(unsigned int iNode = 0; iNode < itCond->GetGeometry().size(); iNode++)
                    {
                        if (itCond->GetGeometry()[iNode].Is(ACTIVE) == true)
                        {
                            DeactivateCondition = false;
                            break;
                        }
                    }
                    
                    // We deactivate the condition if necessary
                    if (DeactivateCondition == true)
                    {
                        itCond->Set(ACTIVE, false);
                    }
                }
            }
        }
        
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                auto& Table = mpTable->GetTable();
                if (IsConverged == true)
                {
                    Table << "Archieved";
                }
                else
                {
                    Table << "Not archieved";
                }
            }
            else
            {
                if (IsConverged == true)
                {
                    std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                }
                else
                {
                    std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                }
            }
        }
        
        return IsConverged;
    }
    
    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        ConvergenceCriteriaBaseType::mConvergenceCriteriaIsInitialized = true;
        
        if (mpTable != nullptr)
        {
            auto& Table = mpTable->GetTable();
            Table.AddColumn("ACTIVE SET CONV", 15);
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
    
    TablePrinterPointerType mpTable;
    
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

