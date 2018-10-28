//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main authors:  Josep Maria Carbonell
//                 Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_RESIDUAL_BASED_NEWMARK_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_NEWMARK_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"

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
 * @class ResidualBasedNewmarkDisplacementScheme
 * @ingroup KratosCore
 * @brief Bossak integration scheme (for dynamic problems) for displacements
 * @details This is a dynamic implicit scheme based of the Newmark algorithm for displacements
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedNewmarkDisplacementScheme
    : public ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
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

    /**
     * @brief Constructor. The Newmark method (parameters)
     * @param ThisParameters Dummy parameters
     */
    explicit ResidualBasedNewmarkDisplacementScheme(Parameters ThisParameters)
      :DerivedBaseType(0.0)
    {
    }

    /**
     * @brief Default constructor. The Newmark method
     */
    explicit ResidualBasedNewmarkDisplacementScheme()
      :DerivedBaseType(0.0)
    {
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedNewmarkDisplacementScheme(ResidualBasedNewmarkDisplacementScheme& rOther)
      :DerivedBaseType(rOther)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedNewmarkDisplacementScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedNewmarkDisplacementScheme
    () override {}

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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedNewmarkDisplacementScheme";
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
