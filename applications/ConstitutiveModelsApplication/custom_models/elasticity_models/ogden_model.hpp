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
#include "custom_models/elasticity_models/hyperelastic_model.hpp"

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
    virtual ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    virtual ~OgdenModel();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Check
     */
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "OgdenModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "OgdenModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    virtual void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix) override;

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix) override;


    /**
     * Calculate Constitutive Components
     */

    virtual double& AddConstitutiveComponentA(HyperElasticDataType& rVariables, array_1d<double,6>& rScalingFactors,
					      MatrixType& rStressDerivatives, double &rCabcd,
					      const unsigned int& a, const unsigned int& b,
					      const unsigned int& c, const unsigned int& d);
    
    virtual double& AddConstitutiveComponentB(HyperElasticDataType& rVariables,
					      array_1d<double,6>& rScalingFactors, double &rCabcd,
					      const unsigned int& a, const unsigned int& b,
					      const unsigned int& c, const unsigned int& d);

    virtual double& AddConstitutiveComponentC(HyperElasticDataType& rVariables,
					      MatrixType& rStressDerivatives, double &rCabcd,
					      const unsigned int& a, const unsigned int& b,
					      const unsigned int& c, const unsigned int& d);

    //************// Strain Data


    virtual void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables) override;

    virtual void CalculateScalingFactors(array_1d<double,6>& rScalingFactors, MatrixType& rStressDerivatives, array_1d<double,3>& rStressEigenValues, array_1d<double,3>& rStrainEigenValues, array_1d<unsigned int,3>& rPermutation);



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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticModel )
    }

    virtual void load(Serializer& rSerializer) override
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
