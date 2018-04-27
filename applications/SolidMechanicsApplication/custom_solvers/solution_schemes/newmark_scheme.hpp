//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NEWMARK_SCHEME_H_INCLUDED)
#define  KRATOS_NEWMARK_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/dynamic_scheme.hpp"

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
  class NewmarkScheme: public DynamicScheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;

    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef DynamicScheme<TSparseSpace,TDenseSpace>                           DerivedType;

    typedef typename DerivedType::NodeType                                       NodeType;

    typedef typename BaseType::IntegrationType                            IntegrationType;
    typedef typename BaseType::IntegrationPointerType              IntegrationPointerType;
    typedef typename BaseType::IntegrationMethodsVectorType  IntegrationMethodsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NewmarkScheme(IntegrationMethodsVectorType& rTimeIntegrationMethods, Flags& rOptions)
        :DerivedType(rTimeIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    NewmarkScheme(IntegrationMethodsVectorType& rTimeIntegrationMethods)
        :DerivedType(rTimeIntegrationMethods)
    {
    }

    /// Copy Constructor.
    NewmarkScheme(NewmarkScheme& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new NewmarkScheme(*this) );
    }

    /// Destructor.
    ~NewmarkScheme() override {}

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NewmarkScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NewmarkScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NewmarkScheme Data";
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

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class NewmarkScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEWMARK_SCHEME_H_INCLUDED defined
