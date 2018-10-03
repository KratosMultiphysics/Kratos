//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_INCOMPRESSIBLE_HYPO_ELASTIC_MODEL_H_INCLUDED)
#define  KRATOS_INCOMPRESSIBLE_HYPO_ELASTIC_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/isochoric_hypo_elastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) IncompressibleHypoElasticModel : public IsochoricHypoElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleHypoElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( IncompressibleHypoElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IncompressibleHypoElasticModel() : IsochoricHypoElasticModel() {}

    /// Copy constructor.
    IncompressibleHypoElasticModel(IncompressibleHypoElasticModel const& rOther) : IsochoricHypoElasticModel(rOther) {}

    /// Assignment operator.
    IncompressibleHypoElasticModel& operator=(IncompressibleHypoElasticModel const& rOther)
    {
	HypoElasticModel::operator=(rOther);
	return *this;
    }

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return Kratos::make_shared<IncompressibleHypoElasticModel>(*this);
    }

    /// Destructor.
    ~IncompressibleHypoElasticModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    // Simplyfied methods must be implemented for performance purposes

    //************************************************************************************
    //************************************************************************************

    int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY

      if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value" << std::endl;

      if(POISSON_RATIO.Key() == 0){
        KRATOS_ERROR << "POISSON_RATIO has Key zero invalid value" << std::endl;
      }
      else{
        const double& nu = rMaterialProperties[POISSON_RATIO];
        if( nu < -0.999 && nu > -1.01 )
          KRATOS_ERROR << "POISSON_RATIO has an invalid value" << std::endl;
      }

      return 0;


      KRATOS_CATCH(" ")
    }


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

      HypoElasticModel::GetDomainVariablesList(rScalarVariables, rComponentVariables);

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
        buffer << "IncompressibleHypoElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleHypoElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "IncompressibleHypoElasticModel Data";
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

    void CalculateAndAddVolumetricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector) override
    {
      KRATOS_TRY

      const ModelDataType& rValues = rVariables.GetModelData();

      //volumetric stress
      const double& Pressure = rValues.GetPressure();
      for(unsigned int i = 0; i < 3; i++)
        rStressVector[i] += Pressure;

      KRATOS_CATCH(" ")
    }


    // set the default volumetric function for the incompressible case

    void CalculateAndAddVolumetricStrainEnergy(ElasticDataType& rVariables, double& rVolumetricDensityFunction) override
    {
      KRATOS_TRY

      KRATOS_ERROR << "calling the class function in IncompressibleHypoElasticModel ... illegal operation" << std::endl;

      KRATOS_CATCH(" ")
    }


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HypoElasticModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HypoElasticModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class IncompressibleHypoElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_HYPO_ELASTIC_MODEL_H_INCLUDED  defined
