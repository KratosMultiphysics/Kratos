//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                LHauser  $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NON_ASSOCIATIVE_YIELD_SURFACE_H_INCLUDED )
#define  KRATOS_NON_ASSOCIATIVE_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
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

  /// Short class definition.
  /** Detail class definition.
   */
  template<class THardeningRule>
  class NonAssociativeYieldSurface : 
     public YieldSurface<THardeningRule>
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef YieldSurface<THardeningRule>                                  BaseType;
    typedef typename YieldSurface<THardeningRule>::Pointer         BaseTypePointer;
    typedef THardeningRule                                       HardeningRuleType;
    typedef typename THardeningRule::PlasticDataType               PlasticDataType;
    typedef typename THardeningRule::InternalVariablesType   InternalVariablesType;


    /// Pointer definition of NonAssociativeYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( NonAssociativeYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonAssociativeYieldSurface() : BaseType() 
     {
        this->mpPlasticPotential = NULL;
     }

    NonAssociativeYieldSurface( BaseTypePointer const & rpPlasticPotential)
       : BaseType() {
          this->mpPlasticPotential = rpPlasticPotential;
       }

    /// Copy constructor.
    NonAssociativeYieldSurface(NonAssociativeYieldSurface const& rOther) : BaseType(rOther) {
       this->mpPlasticPotential = rOther.mpPlasticPotential;
    }

    /// Assignment operator.
    NonAssociativeYieldSurface& operator=(NonAssociativeYieldSurface const& rOther)
    {
       BaseType::operator=(rOther);
       mpPlasticPotential = rOther.mpPlasticPotential;

       return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new NonAssociativeYieldSurface(*this));
    }

    /// Destructor.
    virtual ~NonAssociativeYieldSurface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{



    /**
     * Calculate Plastic Potential Condition Stresses derivative
     */

    virtual VectorType& CalculateDeltaPlasticPotential(const PlasticDataType& rVariables, VectorType& rDeltaPlasticPotential) override
    {
      KRATOS_TRY

      if ( mpPlasticPotential) {
         rDeltaPlasticPotential = mpPlasticPotential->CalculateDeltaStressYieldCondition( rVariables, rDeltaPlasticPotential );
      }else {
         rDeltaPlasticPotential = this->CalculateDeltaStressYieldCondition( rVariables, rDeltaPlasticPotential );
      }
      return rDeltaPlasticPotential;

      KRATOS_CATCH(" ")
    }


    /**
     * Calculate Yield Condition Stresses derivative
     */

    virtual VectorType& CalculateDeltaStressInvPlasticPotential(const PlasticDataType& rVariables, VectorType& rDeltaStressInvPlasticPotential) override
    {
      KRATOS_TRY

      if ( mpPlasticPotential) {
         rDeltaStressInvPlasticPotential = mpPlasticPotential->CalculateDeltaStressInvYieldCondition( rVariables, rDeltaStressInvPlasticPotential );
      }else {
         rDeltaStressInvPlasticPotential = this->CalculateDeltaStressInvYieldCondition( rVariables, rDeltaStressInvPlasticPotential );
      }
      return rDeltaStressInvPlasticPotential;

      KRATOS_CATCH(" ")
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


    /// Turn back information as a string.
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "NonAssociativeYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "NonAssociativeYieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NonAssociativeYieldSurface Data";
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

    BaseTypePointer   mpPlasticPotential;

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
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NonAssociativeYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_YIELD_SURFACE_H_INCLUDED  defined
