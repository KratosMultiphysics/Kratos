//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Josep Maria Carbonell
//  comming from      SolidMechanicsApplication
//
//  Co-author:        Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_RESIDUAL_BASED_NEWMARK_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_NEWMARK_DISPLACEMENT_SCHEME 

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "residual_based_bossak_displacement_scheme.hpp"

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

/** @brief Newmark integration scheme (for dynamic problems)
 */

template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedNewmarkDisplacementScheme: public ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewmarkDisplacementScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                        BaseType;

    typedef typename BaseType::TDataType                                           TDataType;

    typedef typename BaseType::DofsArrayType                                   DofsArrayType;

    typedef typename Element::DofsVectorType                                  DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                           TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                           TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                   LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                   LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                               ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    
    typedef typename BaseType::Pointer                                       BaseTypePointer;
  
    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>  DerivedBaseType;
    
    typedef typename BaseType::LocalSystemComponents               LocalSystemComponentsType;

    ///@}
    ///@name Life Cycle
    ///@{
    ResidualBasedNewmarkDisplacementScheme()
      :DerivedBaseType(0.0)
    {
    }

    /** Copy Constructor.
     */
    ResidualBasedNewmarkDisplacementScheme(ResidualBasedNewmarkDisplacementScheme& rOther)
      :DerivedBaseType(rOther)
    {
    }

    /**
     * Clone 
     */
    virtual BaseTypePointer Clone()
    {
      return BaseTypePointer( new ResidualBasedNewmarkDisplacementScheme(*this) );
    }
    
    /** Destructor.
     */
    virtual ~ResidualBasedNewmarkDisplacementScheme
    () {}

    ///@}
    ///@name Operators
    ///@{

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
    
    ///@}
    ///@name Friends
    ///@{

protected:
    ///@}
    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Protected  Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations*
    ///@{
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
    ///@{
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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedNewmarkDisplacementScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_NEWMARK_DISPLACEMENT_SCHEME E defined */
