//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MOONEY_RIVLIN_MODEL_H_INCLUDED )
#define  KRATOS_MOONEY_RIVLIN_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MooneyRivlinModel : public HyperElasticModel
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MooneyRivlinModel
    KRATOS_CLASS_POINTER_DEFINITION( MooneyRivlinModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MooneyRivlinModel();

    /// Copy constructor.
    MooneyRivlinModel(MooneyRivlinModel const& rOther);

    /// Assignment operator.
    MooneyRivlinModel& operator=(MooneyRivlinModel const& rOther);

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    ~MooneyRivlinModel() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Check
     */
    int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "MooneyRivlinModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MooneyRivlinModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "MooneyRivlinModel Data";
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
     * Calculate Constitutive Components
     */

    double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
					     const unsigned int& a, const unsigned int& b,
					     const unsigned int& c, const unsigned int& d) override;


    //************// Strain Data


    void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables) override;

    void CalculateScalingFactors(HyperElasticDataType& rVariables) override;


    //************// dW

    double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override; //dU/dJ

    double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative) override; //ddU/dJdJ


    virtual double& GetFunction1stI1Derivative(HyperElasticDataType& rVariables, double& rDerivative); //dW/dI1

    virtual double& GetFunction1stI2Derivative(HyperElasticDataType& rVariables, double& rDerivative); //dW/dI2

    virtual double& GetFunction1stI3Derivative(HyperElasticDataType& rVariables, double& rDerivative); //dW/dI3


    virtual double& GetFunction2ndI1Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI1dI1

    virtual double& GetFunction2ndI2Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI2dI2

    virtual double& GetFunction2ndI3Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI3dI3

    // the implementation of the crossed derivatives have to be added for a more general form (usually they are zero)
    // virtual double& GetFunction2ndI2I1Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI2dI1
    // virtual double& GetFunction2ndI3I1Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI3dI1
    // virtual double& GetFunction2ndI1I2Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI1dI2
    // virtual double& GetFunction2ndI3I2Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI3dI2
    // virtual double& GetFunction2ndI1I3Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI1dI3
    // virtual double& GetFunction2ndI2I3Derivative(HyperElasticDataType& rVariables, double& rDerivative); //ddW/dI2dI3

    //isochoric volumetric slit

    MatrixType& GetIsochoricRightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dC'/dC

    double& GetIsochoricRightCauchyGreenDerivative(const StrainData& rStrain,
						   double& rDerivative,
						   const double& a,
						   const double& b,
						   const double& c,
						   const double& d); //dC'/dC


    MatrixType& GetIsochoricLeftCauchyGreenDerivative(const StrainData& rStrain,
						      MatrixType& rDerivative); //db'/db

    double& GetIsochoricLeftCauchyGreenDerivative(const StrainData& rStrain,
						  double& rDerivative,
						  const double& a,
						  const double& b,
						  const double& c,
						  const double& d); //db'/db

    //************// right cauchy green: C

    MatrixType& GetI1RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dI1/dC

    MatrixType& GetI2RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dI2/dC

    MatrixType& GetI3RightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dI3/dC



    double& GetFourthOrderUnitTensor(double& rValue,
				     const double& a,
				     const double& b,
				     const double& c,
				     const double& d); //ddC/dCdC or ddb/dbdb


    double& GetInverseRightCauchyGreenDerivative(const StrainData& rStrain,
						 double& rDerivative,
						 const double& a,
						 const double& b,
						 const double& c,
						 const double& d); //dC^-1/dC


    //Invariants 1st derivatives by components
    double& GetI1RightCauchyGreen1stDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b); //dI1/dC


    double& GetI2RightCauchyGreen1stDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b); //dI2/dC


    double& GetI3RightCauchyGreen1stDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b); //dI3/dC





    //Invariants Square of the 1st derivatives by components
    double& GetI1RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						     double& rDerivative,
						     const double& a,
						     const double& b,
						     const double& c,
						     const double& d); //dI1/dC * dI2/dC


    double& GetI2RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						     double& rDerivative,
						     const double& a,
						     const double& b,
						     const double& c,
						     const double& d); //dI2/dC * dI3/dC


    double& GetI3RightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						     double& rDerivative,
						     const double& a,
						     const double& b,
						     const double& c,
						     const double& d); //dI3/dC * dI3/dC



    //Invariants 2nd derivatives by components
    double& GetI1RightCauchyGreen2ndDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b,
					       const double& c,
					       const double& d); //ddI1/dCdC


    double& GetI2RightCauchyGreen2ndDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b,
					       const double& c,
					       const double& d); //ddI2/dCdC


    double& GetI3RightCauchyGreen2ndDerivative(const StrainData& rStrain,
					       double& rDerivative,
					       const double& a,
					       const double& b,
					       const double& c,
					       const double& d); //ddI3/dCdC


    //************// left cauchy green : b

    MatrixType& GetI1LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dI1/db

    MatrixType& GetI2LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative);  //dI2/db

    MatrixType& GetI3LeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dI3/db


    //Invariants 1st derivatives by components
    double& GetI1LeftCauchyGreen1stDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b); //dI1/db


    double& GetI2LeftCauchyGreen1stDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b); //dI2/db


    double& GetI3LeftCauchyGreen1stDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b); //dI3/db



    //Invariants Square of the 1st derivatives by components
    double& GetI1LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						    double& rDerivative,
						    const double& a,
						    const double& b,
						    const double& c,
						    const double& d); //dI1/db * dI1/db


    double& GetI2LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						    double& rDerivative,
						    const double& a,
						    const double& b,
						    const double& c,
						    const double& d); //dI2/db * dI2/db


    double& GetI3LeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						    double& rDerivative,
						    const double& a,
						    const double& b,
						    const double& c,
						    const double& d); //dI3/db * dI3/db



    //Invariants 2nd derivatives by components
    double& GetI1LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b,
					      const double& c,
					      const double& d); //ddI1/dbdb


    double& GetI2LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b,
					      const double& c,
					      const double& d); //ddI2/dbdb


    double& GetI3LeftCauchyGreen2ndDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b,
					      const double& c,
					      const double& d); //ddI3/dbdb


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

  }; // Class MooneyRivlinModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOONEY_RIVLIN_MODEL_H_INCLUDED  defined


