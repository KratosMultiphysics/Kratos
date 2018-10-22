//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_HYPO_ELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_HYPO_ELASTIC_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HypoElasticModel : public ConstitutiveModel
  {
  protected:

    struct ElasticModelData
    {
    private:

      Flags*               mpState;
      const ModelDataType* mpModelData;

    public:

      BoundedMatrix<double,6,6> ConstitutiveTensor;
      MatrixType                StrainMatrix;

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
    typedef ElasticModelData              ElasticDataType;


    /// Pointer definition of HypoElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( HypoElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HypoElasticModel();

    /// Copy constructor.
    HypoElasticModel(HypoElasticModel const& rOther);

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override;

    /// Assignment operator.
    HypoElasticModel& operator=(HypoElasticModel const& rOther);

    /// Destructor.
    ~HypoElasticModel() override;


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
    int Check(const Properties& rProperties, const ProcessInfo& rCurrentProcessInfo) override;

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

      rComponentVariables.push_back(VELOCITY);

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
        buffer << "HypoElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HypoElasticModel";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "HypoElasticModel Data";
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
     * Stress update for a hypoelastic model
     */
    void AddHistoricalStress(ModelDataType& rValues, MatrixType& rStressMatrix);

    /**
     * Calculate Stresses
     */
    virtual void CalculateAndAddStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    virtual void CalculateAndAddIsochoricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    virtual void CalculateAndAddVolumetricStressTensor(ElasticDataType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables);

    virtual void CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables);


    virtual void CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables);

    //************//

    void InitializeElasticData(ModelDataType& rValues, ElasticDataType& rVariables);

    virtual void CalculateAndAddIsochoricStrainEnergy(ElasticDataType& rVariables, double& rIsochoricDensityFunction);

    virtual void CalculateAndAddVolumetricStrainEnergy(ElasticDataType& rVariables, double& rIsochoricDensityFunction);


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

  }; // Class HypoElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HYPO_ELASTIC_MODEL_H_INCLUDED  defined
