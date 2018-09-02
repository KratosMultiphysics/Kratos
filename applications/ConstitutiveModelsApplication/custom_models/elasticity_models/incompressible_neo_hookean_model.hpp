//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED )
#define  KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_neo_hookean_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IncompressibleNeoHookeanModel : public IsochoricNeoHookeanModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleNeoHookeanModel
    KRATOS_CLASS_POINTER_DEFINITION( IncompressibleNeoHookeanModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IncompressibleNeoHookeanModel() : IsochoricNeoHookeanModel() {}

    /// Copy constructor.
    IncompressibleNeoHookeanModel(IncompressibleNeoHookeanModel const& rOther) : IsochoricNeoHookeanModel(rOther) {}

    /// Assignment operator.
    IncompressibleNeoHookeanModel& operator=(IncompressibleNeoHookeanModel const& rOther)
    {
      IsochoricNeoHookeanModel::operator=(rOther);
      return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IncompressibleNeoHookeanModel>(*this);
    }

    /// Destructor.
    ~IncompressibleNeoHookeanModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    // Simplyfied methods must be implemented for performance purposes
    /**
     * Calculate Stresses
     */

    /**
     * Calculate Constitutive Components
     */

    /**
     * Check
     */
    int Check(const Properties& rMaterialProperties,
		      const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      IsochoricNeoHookeanModel::Check(rMaterialProperties,rCurrentProcessInfo);

      return 0;

      KRATOS_CATCH(" ")
    };


    ///@}
    ///@name Access
    ///@{

    /**
     * method to ask the constitutive model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
                                std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
    {
      KRATOS_TRY

      HyperElasticModel::GetDomainVariablesList(rScalarVariables, rComponentVariables);

      rScalarVariables.push_back(PRESSURE);

      KRATOS_CATCH(" ")
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressibleHyperElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleHyperElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IncompressibleHyperElasticModel Data";
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


    double& AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						       const unsigned int& a, const unsigned int& b,
						       const unsigned int& c, const unsigned int& d) override
    {
      KRATOS_TRY

      return IsochoricMooneyRivlinModel::AddVolumetricConstitutiveComponent(rVariables,rCabcd,a,b,c,d);

      KRATOS_CATCH(" ")
    }

    //************// dW

    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //dU/dJ
    {
      KRATOS_TRY

      const ModelDataType&  rValues = rVariables.GetModelData();

      rDerivative = rValues.GetPressure();

      return rDerivative;

      KRATOS_CATCH(" ")
    };


    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override //ddU/dJdJ
    {
      KRATOS_TRY

      rDerivative = 0.0;

      return rDerivative;

      KRATOS_CATCH(" ")
    };


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


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsochoricNeoHookeanModel )
    }


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IncompressibleNeoHookeanModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_NEO_HOOKEAN_MODEL_H_INCLUDED  defined
