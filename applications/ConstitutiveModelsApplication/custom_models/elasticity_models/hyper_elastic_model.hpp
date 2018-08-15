//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_HYPER_ELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_HYPER_ELASTIC_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_models/constitutive_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HyperElasticModel : public ConstitutiveModel
  {
  protected:


    struct StrainInvariants
    {
      double I1;
      double I2;
      double I3;

      double J;
      double J_13;

    };


    struct HyperElasticFactors
    {
      double Alpha1;  //1st derivative I1
      double Alpha2;  //1st derivative I2
      double Alpha3;  //1st derivative I3
      double Alpha4;  //1st derivative J

      double Beta1;   //2nd derivative I1
      double Beta2;   //2nd derivative I2
      double Beta3;   //2nd derivative I3
      double Beta4;   //2nd derivative J

      // the implementation of the crossed derivatives have to be added for a more general form (usually they are zero)
      // double Gamma21;  //2nd derivative ddW/dI2dI1
      // double Gamma31;  //2nd derivative ddW/dI3dI1
      // double Gamma12;  //2nd derivative ddW/dI1dI2
      // double Gamma32;  //2nd derivative ddW/dI3dI2
      // double Gamma13;  //2nd derivative ddW/dI1dI3
      // double Gamma23;  //2nd derivative ddW/dI2dI3

    };

    struct StrainEigenData
    {
       array_1d<double,3> Values;
       MatrixType        Vectors;
    };


    struct StrainData
    {

      StrainInvariants Invariants;
      StrainEigenData  Eigen;

      MatrixType       Matrix; //left(b) or right(C) cauchy green
      MatrixType       InverseMatrix; //insverse right(C) cauchy green

    };


    struct HyperElasticModelData
    {
    private:

      Flags*               mpState;
      const ModelDataType* mpModelData;

    public:

      HyperElasticFactors    Factors;
      StrainData              Strain;

      //Set Data Pointers
      void SetState           (Flags& rState)                    {mpState = &rState;};
      void SetModelData       (const ModelDataType&  rModelData) {mpModelData = &rModelData;};

      //Get Data Pointers
      const ModelDataType&    GetModelData                () const {return *mpModelData;};
      const MaterialDataType& GetMaterialParameters       () const {return mpModelData->GetMaterialParameters();};

      //Get non const Data
      Flags& State                                        () {return *mpState;};

      //Get const Data
      const Flags&  GetState                              () const {return *mpState;};


    };


  public:

    ///@name Type Definitions
    ///@{
    typedef HyperElasticModelData        HyperElasticDataType;

    /// Pointer definition of HyperElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HyperElasticModel();

    /// Copy constructor.
    HyperElasticModel(HyperElasticModel const& rOther);

    /// Assignment operator.
    HyperElasticModel& operator=(HyperElasticModel const& rOther);

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override;


    /// Destructor.
    ~HyperElasticModel() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize member data
     */
    void InitializeModel(ModelDataType& rValues) override;

    /**
     * Finalize member data
     */
    void FinalizeModel(ModelDataType& rValues) override;


    /**
     * Calculate Strain Energy Density Functions
     */
    void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override;


    /**
     * Calculate Stresses
     */
    void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;

    void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;

    void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;


    /**
     * Calculate Constitutive Tensor
     */
    void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override;

    void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override;

    void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override;


    /**
     * Calculate Stress and Constitutive Tensor
     */
    void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;

    void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;

    void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;


    /**
     * Check
     */
    int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    void SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue,
			  const ProcessInfo& rCurrentProcessInfo ) override
    {
      KRATOS_TRY

      // A method to compute the initial linear strain from the stress is needed
      //if(rThisVariable == INITIAL_STRESS_VECTOR)

      // A method to compute the initial linear strain from the stress is needed
      // if(rThisVariable == INITIAL_STRAIN_VECTOR){
      //   this->mHistoryVector = rValue;
      // }

      KRATOS_CATCH(" ")
    }


    void SetValue(const Variable<Matrix>& rThisVariable, const Matrix& rValue,
			  const ProcessInfo& rCurrentProcessInfo ) override
    {
      KRATOS_TRY

      // A method to compute the initial linear strain from the stress is needed
      //if(rThisVariable == INITIAL_STRESS_VECTOR)

      // A method to compute the initial linear strain from the stress is needed
      // if(rThisVariable == INITIAL_STRAIN_VECTOR){
      //   this->mHistoryVector = rValue;
      // }

      KRATOS_CATCH(" ")
    }

    /**
     * method to ask the plasticity model the list of variables (dofs)  needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
					std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
    {
      KRATOS_TRY

      rComponentVariables.push_back(DISPLACEMENT);

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
        buffer << "HyperElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HyperElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "HyperElasticModel Data";
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    static const MatrixType msIdentityMatrix;


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
    virtual void CalculateAndAddStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix);

    virtual void CalculateAndAddIsochoricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix);

    virtual void CalculateAndAddVolumetricStressTensor(HyperElasticDataType& rVariables, MatrixType& rStressMatrix);

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddIsochoricConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddVolumetricConstitutiveTensor(HyperElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    /**
     * Calculate Constitutive Components
     */

    virtual double& AddConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
					     const unsigned int& a, const unsigned int& b,
					     const unsigned int& c, const unsigned int& d);


    virtual double& AddIsochoricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						      const unsigned int& a, const unsigned int& b,
						      const unsigned int& c, const unsigned int& d);


    virtual double& AddVolumetricConstitutiveComponent(HyperElasticDataType& rVariables, double &rCabcd,
						       const unsigned int& a, const unsigned int& b,
						       const unsigned int& c, const unsigned int& d);


    //************// Strain Data


    virtual void CalculateStrainData(ModelDataType& rValues, HyperElasticDataType& rVariables);

    virtual void CalculateInvariants(HyperElasticDataType& rVariables);

    virtual void CalculateScalingFactors(HyperElasticDataType& rVariables);

    void CalculateStrainInvariants(const MatrixType& rStrainMatrix, double& rI1, double& rI2, double& rI3);


    //************//W

    virtual void CalculateAndAddIsochoricStrainEnergy(HyperElasticDataType& rVariables, double& rIsochoricDensityFunction);

    virtual void CalculateAndAddVolumetricStrainEnergy(HyperElasticDataType& rVariables, double& rVolumetricDensityFunction);


    //************// dW

    virtual double& GetVolumetricFunction1stJDerivative(HyperElasticDataType& rVariables, double& rDerivative); //dU/dJ

    virtual double& GetVolumetricFunction2ndJDerivative(HyperElasticDataType& rVariables, double& rDerivative); //ddU/dJdJ


    //************// right cauchy green: C
    MatrixType& GetJRightCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dJ/dC

    double& GetJRightCauchyGreen1stDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b); ///dJ/dC

    double& GetJRightCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						    double& rDerivative,
						    const double& a,
						    const double& b,
						    const double& c,
						    const double& d); //dJ/dC * dJ/dC

    double& GetJRightCauchyGreen2ndDerivative(const StrainData& rStrain,
					      double& rDerivative,
					      const double& a,
					      const double& b,
					      const double& c,
					      const double& d); //ddJ/dCdC

    //************// left cauchy green : b

    MatrixType& GetJLeftCauchyGreenDerivative(const StrainData& rStrain, MatrixType& rDerivative); //dJ/db

    double& GetJLeftCauchyGreen1stDerivative(const StrainData& rStrain,
					     double& rDerivative,
					     const double& a,
					     const double& b); //dJ/db

    double& GetJLeftCauchyGreenSquare1stDerivative(const StrainData& rStrain,
						   double& rDerivative,
						   const double& a,
						   const double& b,
						   const double& c,
						   const double& d); //dJ/db * dJ/db

    double& GetJLeftCauchyGreen2ndDerivative(const StrainData& rStrain,
					     double& rDerivative,
					     const double& a,
					     const double& b,
					     const double& c,
					     const double& d); //ddJ/dbdb


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveModel )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveModel )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class HyperElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HYPER_ELASTIC_MODEL_H_INCLUDED  defined
