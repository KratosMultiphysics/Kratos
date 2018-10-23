//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_OGDEN_MODEL_H_INCLUDED )
#define  KRATOS_OGDEN_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/elasticity_models/hyper_elastic_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) OgdenModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of OgdenModel
    KRATOS_CLASS_POINTER_DEFINITION( OgdenModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OgdenModel();

    /// Copy constructor.
    OgdenModel(OgdenModel const& rOther);

    /// Assignment operator.
    OgdenModel& operator=(OgdenModel const& rOther);

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    ~OgdenModel() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

   /**
     * Calculate Strain Energy Density Functions
     */
    void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override;

    /**
     * Calculate Constitutive Tensor
     */
    void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override;


    /**
     * Check
     */
    int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override;

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
        std::stringstream buffer;
        buffer << "OgdenModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "OgdenModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "OgdenModel Data";
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

    /**
     * Calculate Stresses
     */
    void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override;

    /**
     * Calculate Constitutive Tensor
     */
    void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix) override;

    virtual void CalculateAndAddConstitutiveTensorB(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    /**
     * Calculate Constitutive Components
     */
    double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
				     const array_1d<double,3>& rVectorDerivative,
				     const unsigned int& a, const unsigned int& b,
				     const unsigned int& c, const unsigned int& d);

    /**
     * Checks eigen values coincidence
     */
    void GetEigenCoincidence(const array_1d<double,3>& rStrainEigenValues,
			     array_1d<unsigned int,3>& Order,
			     unsigned int& option);

    /**
     * Calculate Derivative of a general isotropic tensor
     */
    double& CalculateIsotropicTensorDerivative(const MatrixType& rStrainMatrix,
					       const MatrixType& rStrainEigenVectors,
					       const array_1d<double,3>& rStrainEigenValues,
					       const MatrixType& rStressDerivatives,
					       const array_1d<double,3>& rStressEigenValues,
					       const array_1d<double,6>& rOptionFactors,
					       const unsigned int& rOption,
					       double &rCabcd,
					       const unsigned int& a, const unsigned int& b,
					       const unsigned int& c, const unsigned int& d);


    /**
     * Calculate Tensor Derivative Factors
     */
    void CalculateDerivativeFactors(array_1d<double,6>& rDerivativeFactors, const MatrixType& rStressDerivatives, const array_1d<double,3>& rStressEigenValues, const array_1d<double,3>& rStrainEigenValues, const array_1d<unsigned int,3>& rOrder);


    //************// Strain Data


    void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables) override;

    virtual void CalculateMainStresses(HyperElasticDataType& rVariables, array_1d<double,3>& rMainStresses);

    virtual void CalculateMainStressDerivatives(HyperElasticDataType& rVariables, MatrixType& rStressDerivatives);

    virtual double& CalculateStressDerivativesI(HyperElasticDataType& rVariables, double& rValue, const unsigned int& i, const unsigned int& j);
    virtual double& CalculateStressDerivativesII(HyperElasticDataType& rVariables, double& rValue, const unsigned int& i, const unsigned int& j);

    //************//W

    void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction) override;

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

    using HyperElasticModel::AddConstitutiveComponent;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class OgdenModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OGDEN_MODEL_H_INCLUDED  defined
