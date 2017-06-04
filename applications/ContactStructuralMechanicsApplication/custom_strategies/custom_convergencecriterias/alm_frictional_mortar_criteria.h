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

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "custom_utilities/contact_utilities.h"
#include "custom_utilities/bprinter_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#if !defined(_WIN32)
	#include "custom_utilities/color_utilities.h"
#endif

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
class ALMFrictionalMortarConvergenceCriteria : public virtual  BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ALMFrictionalMortarConvergenceCriteria );

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
    ALMFrictionalMortarConvergenceCriteria(        
        double Tolerance = std::numeric_limits<double>::epsilon(),
        TablePrinterPointerType pTable = nullptr
        ) : BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >(),
        mTolerance(Tolerance),
        mpTable(pTable),
        mTableIsInitialized(false)
    {
    }

    ///Copy constructor 
    ALMFrictionalMortarConvergenceCriteria( ALMFrictionalMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    virtual ~ALMFrictionalMortarConvergenceCriteria() {}

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
        
        // Defining the convergence
        bool IsConvergedActive = true;
        bool IsConvergedSlip = true;
        
        const double Epsilon = rModelPart.GetProcessInfo()[PENALTY_PARAMETER]; 
        const double ScaleFactor = rModelPart.GetProcessInfo()[SCALE_FACTOR];
        const double TangentFactor = rModelPart.GetProcessInfo()[TANGENT_FACTOR];
        
        const array_1d<double,3> ZeroVector(0.0);
        
        NodesArrayType& NodesArray = rModelPart.GetSubModelPart("Contact").Nodes();
        int numNodes = static_cast<int>(NodesArray.size());

        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = NodesArray.begin() + i;
            
            // Check if the node is slave
            bool NodeIsSlave = true;
            if ((itNode)->IsDefined(SLAVE))
            {
                NodeIsSlave = (itNode)->Is(SLAVE);
            }
            
            if (NodeIsSlave == true)
            {
                const array_1d<double,3> LagrangeMultiplier = (itNode)->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                const array_1d<double,3> NodalNormal = (itNode)->GetValue(NORMAL);
                const double NormalLagrangeMultiplier = inner_prod(NodalNormal, LagrangeMultiplier);
                
                const double AugmentedNormalPressure = ScaleFactor * NormalLagrangeMultiplier + Epsilon * (itNode)->FastGetSolutionStepValue(WEIGHTED_GAP);     
                
                itNode->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, AugmentedNormalPressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)
                
                if (AugmentedNormalPressure < mTolerance * ScaleFactor) // NOTE: This could be conflictive (< or <=)
                {
                    if ((itNode)->Is(ACTIVE) == false )
                    {
                        (itNode)->Set(ACTIVE, true);
                        IsConvergedActive = false;
                    }
                    
                    // Computing the augmented tangent pressure
                    const array_1d<double,3> TangentLagrangeMultiplier = LagrangeMultiplier - NormalLagrangeMultiplier * NodalNormal;
                    const double LambdaTangent = norm_2(TangentLagrangeMultiplier); 
                    
                    // The friction coefficient
                    const double mu = (itNode)->GetValue(WEIGHTED_FRICTION);
                    
                    // Finally we compute the augmented tangent pressure
                    const double gt = (itNode)->FastGetSolutionStepValue(WEIGHTED_SLIP);
                    const double AugmentedTangentPressure = std::abs(ScaleFactor * LambdaTangent + TangentFactor * Epsilon * gt) + mu * AugmentedNormalPressure;
                    
                    (itNode)->SetValue(AUGMENTED_TANGENT_CONTACT_PRESSURE, AugmentedTangentPressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)
                    
                    if (AugmentedTangentPressure <= 0.0) // TODO: Check if it is minor equal or just minor
                    {
                        if ((itNode)->Is(SLIP) == true )
                        {
                            (itNode)->Set(SLIP, false);
                            IsConvergedSlip = false;
                        }
                    }
                    else
                    {
                        if ((itNode)->Is(SLIP) == false )
                        {
                            (itNode)->Set(SLIP, true);
                            IsConvergedSlip = false;
                        }
                    }   
                }
                else
                {
                    (itNode)->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER) = ZeroVector; // NOTE: To clear the value (can affect future iterations)
                    
                    if ((itNode)->Is(ACTIVE) == true )
                    {
                        (itNode)->Set(ACTIVE, false);
                        IsConvergedActive = false;
                    }
                }
            }
        }
        
        // We update the pairs if necessary
        if (rModelPart.GetProcessInfo()[CONSIDER_PAIR_VARIATION] == true)
        {
            ConditionsArrayType& ConditionsArray = rModelPart.GetSubModelPart("Contact").Conditions();
            int numConditions = static_cast<int>(ConditionsArray.size());

            #pragma omp parallel for 
            for(int i = 0; i < numConditions; i++) 
            {
                auto itCond = ConditionsArray.begin() + i;
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
                if (IsConvergedActive == true)
                {
                    #if !defined(_WIN32)
                            Table << BOLDFONT(FGRN("       Achieved"));
                    #else
                            Table << "Achieved";
                    #endif
                }
                else
                {
                    #if !defined(_WIN32)
                            Table << BOLDFONT(FRED("   Not achieved"));
                    #else
                            Table << "Not achieved";
                    #endif
                }
                if (IsConvergedSlip == true)
                {
                    #if !defined(_WIN32)
                            Table << BOLDFONT(FGRN("       Achieved"));
                    #else
                            Table << "Achieved";
                    #endif
                }
                else
                {
                    #if !defined(_WIN32)
                            Table << BOLDFONT(FRED("   Not achieved"));
                    #else
                            Table << "Not achieved";
                    #endif
                }
            }
            else
            {
                if (IsConvergedActive == true)
                {
                    #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    #else
                            std::cout << "\tActive set convergence is achieved" << std::endl;
                    #endif
                }
                else
                {
                    #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    #else
                            std::cout << "\tActive set convergence is not achieved" << std::endl;
                    #endif
                }
                
                if (IsConvergedSlip == true)
                {
                    #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    #else
                            std::cout << "\tSlip/stick set convergence is achieved" << std::endl;
                    #endif
                }
                else
                {
                    #if !defined(_WIN32)
                            std::cout << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    #else
                            std::cout << "\tSlip/stick set  convergence is not achieved" << std::endl;
                    #endif
                }
            }
        }
        
        return IsConvergedActive && IsConvergedSlip;
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
            mpTable->AddColumn("ACTIVE SET CONV.", 16);
            mpTable->AddColumn("SLIP/STICK CONV.", 16);
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

}; // Class ALMFrictionalMortarConvergenceCriteria 

///@name Explicit Specializations
///@{

}  // namespace Kratos 

#endif /* KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H  defined */

