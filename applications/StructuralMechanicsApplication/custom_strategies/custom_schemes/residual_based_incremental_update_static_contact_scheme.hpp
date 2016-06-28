// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

#if !defined(RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME)
#define RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME

/* System Includes */

/* External Includes */

/* Project includes */
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/point.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/contact_utilities.h"
#include "includes/condition.h"

namespace Kratos {

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
    
/** \brief  Short class definition.
This class provides the implementation of the condition contribution methods of the solution strategy
for a contact problem employing the standard lagrange multipliers constraining formulation 
*/

template<class TSparseSpace, class TDenseSpace >
class ResidualBasedIncrementalUpdateStaticContactScheme :
    public ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace >
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticContactScheme );

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    
    typedef typename BaseType::ElementsArrayType   ElementsArrayType;      
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    
    typedef Point<3> PointType;
    
    
    ResidualBasedIncrementalUpdateStaticContactScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>( )
    {}
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;
        
        ( rCurrentCondition )->InitializeNonLinearIteration( CurrentProcessInfo );
        BaseType::Condition_CalculateSystemContributions( rCurrentCondition, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );
        
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY;
        
        ( rCurrentCondition )->InitializeNonLinearIteration( CurrentProcessInfo );
        BaseType::Condition_Calculate_RHS_Contribution( rCurrentCondition, RHS_Contribution, EquationId, CurrentProcessInfo );
        
        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    void InitializeNonLinIteration(
      ModelPart& r_model_part,
      TSystemMatrixType& A,
      TSystemVectorType& Dx,
      TSystemVectorType& b
    )
    {
        KRATOS_TRY;
        
        // Update normal of the conditions
        ContactUtilities::ComputeNodesMeanNormalModelPart( r_model_part );
        
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        
        ElementsArrayType& pElements = r_model_part.Elements();
        for( typename ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it )
        {
            (it)->InitializeNonLinearIteration(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for( typename ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it )
        {
            if ( it->Is( ACTIVE ) )
            {
                (it)->InitializeNonLinearIteration(CurrentProcessInfo);
            }
        }
        
        KRATOS_CATCH( "" );
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont: The condition to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
      Condition::Pointer rCurrentCondition,
      ProcessInfo& CurrentProcessInfo
    )
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentConditiont: The element to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(
      Element::Pointer rCurrentElement,
      ProcessInfo& CurrentProcessInfo
    )
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }
    
};

}  // namespace Kratos

#endif /* RESIDUAL_BASED_INCREMENTAL_UPDATE_STATIC_CONTACT_SCHEME */
