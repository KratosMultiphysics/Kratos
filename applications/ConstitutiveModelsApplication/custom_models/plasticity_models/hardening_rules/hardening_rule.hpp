//
//   Project Name:        KratosConstitutiveModelApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                   July 2017 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_HARDENING_RULE_H_INCLUDED )
#define  KRATOS_HARDENING_RULE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"

#include "constitutive_models_application_variables.h"

#include "custom_models/constitutive_model_data.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HardeningRule
  {
  protected:

    //warning::this variable is going to be shadowed by they derived classes
    //if any problem is detected an alternative method must be used instead
    constexpr static std::size_t VarSize = 1;

  public:

    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                  MatrixType;
    typedef ConstitutiveModelData::VectorType                  VectorType;
    typedef ConstitutiveModelData::ModelData                ModelDataType;
    typedef ConstitutiveModelData::MaterialData          MaterialDataType;

    template<std::size_t TVarSize>
    struct InternalVariables
    {
      //internal variables
      array_1d<double, TVarSize> Variables;

      //default constructor (to initialize Variables)
      InternalVariables() { Variables.clear(); };

      const array_1d<double, TVarSize>& GetVariables() {return Variables;};

      unsigned int size() {return TVarSize;}

    private:

      friend class Serializer;

      void save(Serializer& rSerializer) const
      {
	rSerializer.save("Variables",Variables);
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("Variables",Variables);
      };

    };


    template<std::size_t TVarSize>
    struct PlasticModelData
    {
    private:

      Flags*               mpState;
      const ModelDataType* mpModelData;

    public:

      //flow rule internal variables
      double TrialStateFunction;
      double StressNorm;

      //hardening rule internal variables
      double RateFactor;

      //internal variables
      InternalVariables<TVarSize>      Internal;
      InternalVariables<TVarSize> DeltaInternal;

      //strain matrix
      MatrixType StrainMatrix; //wildcard strain (cauchy green tensors or infinitessimal tensor)

      //Set Data Pointers
      void SetState           (Flags& rState)                    {mpState = &rState;};
      void SetModelData       (const ModelDataType&  rModelData) {mpModelData = &rModelData;};

      //Get Data Pointers
      const ModelDataType&    GetModelData                () const {return *mpModelData;};
      const MaterialDataType& GetMaterialParameters       () const {return mpModelData->GetMaterialParameters();};

      //Get non const Data
      Flags& State                                        () {return *mpState;};

      //Get const Data
      const Flags&  GetState              () const {return *mpState;};
      const double& GetTrialStateFunction () const {return TrialStateFunction;};
      const double& GetStressNorm         () const {return StressNorm;};
      const double& GetRateFactor         () const {return RateFactor;};

      const InternalVariables<TVarSize>&      GetInternal () const {return Internal;};
      const InternalVariables<TVarSize>& GetDeltaInternal () const {return DeltaInternal;};

      const MatrixType&                   GetStrainMatrix () const {return StrainMatrix;};

      const array_1d<double,TVarSize>& GetInternalVariables       () const {return Internal.Variables;};
      const array_1d<double,TVarSize>& GetDeltaInternalVariables  () const {return DeltaInternal.Variables;};

    };

    typedef InternalVariables<VarSize>   InternalVariablesType;
    typedef PlasticModelData<VarSize>          PlasticDataType;

    /// Pointer definition of HardeningRule
    KRATOS_CLASS_POINTER_DEFINITION( HardeningRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HardeningRule() {}

    /// Copy constructor.
    HardeningRule(HardeningRule const& rOther) {}

    /// Assignment operator.
    HardeningRule& operator=(HardeningRule const& rOther)
      {
	return *this;
      }

    /// Clone.
    virtual HardeningRule::Pointer Clone() const
    {
      return Kratos::make_shared<HardeningRule>(*this);
    }

    /// Destructor.
    virtual ~HardeningRule() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Hardening functions
     */

    virtual double& CalculateHardening(const PlasticDataType& rVariables, double& rHardening)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the HardeningRule base class ... illegal operation" << std::endl;

      return rHardening;

      KRATOS_CATCH(" ")
    }


    /**
     * Calculate Hardening function derivatives
     */

    virtual double& CalculateDeltaHardening(const PlasticDataType& rVariables, double& rDeltaHardening)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the HardeningRule base class ... illegal operation" << std::endl;

      return rDeltaHardening;

      KRATOS_CATCH(" ")
    }


    virtual double& CalculateDeltaThermalHardening(const PlasticDataType& rVariables, double& rDeltaThermalHardening)
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the HardeningRule base class ... illegal operation" << std::endl;

      return rDeltaThermalHardening;

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
    virtual std::string Info() const
    {
      std::stringstream buffer;
      buffer << "HardeningRule" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "HardeningRule";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "HardeningRule Data";
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
    friend class Serializer;


    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class HardeningRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HARDENING_RULE_H_INCLUDED  defined
