//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_LINEAR_ELASTIC_MODEL_H_INCLUDED )
#define  KRATOS_LINEAR_ELASTIC_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/constitutive_law_utilities.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) LinearElasticModel
  {
  protected:

    struct ElasticModelData
    {
    private:

      Flags*               mpState;
      const ModelDataType* mpModelData;
      
    public:
      
      bounded_matrix<double,6,6>   ConstitutiveTensor;
      
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

      
    /// Pointer definition of LinearElasticModel
    KRATOS_CLASS_POINTER_DEFINITION( LinearElasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.    
    LinearElasticModel() {}

    /// Copy constructor.
    LinearElasticModel(LinearElasticModel const& rOther) {}

    /// Clone.
    virtual LinearElasticModel::Pointer Clone() const
    {
      return ( LinearElasticModel::Pointer(new LinearElasticModel(*this)) );
    }

    /// Assignment operator.
    LinearElasticModel& operator=(LinearElasticModel const& rOther) {return *this;}


    /// Destructor.
    virtual ~LinearElasticModel(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Strain Energy Density Functions
     */
    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction) override;

      
    /**
     * Calculate Stresses
     */    
    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;

    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;


    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override; 
    
    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override; 

    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override; 
    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;

    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;
    
    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;

    
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "LinearElasticModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "LinearElasticModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "LinearElasticModel Data";
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
    virtual void CalculateAndAddStressTensor(ElasticVariablesType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    virtual void CalculateAndAddIsochoricStressTensor(ElasticVariablesType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    virtual void CalculateAndAddVolumetricStressTensor(ElasticVariablesType& rVariables, VectorType& rStrainVector, VectorType& rStressVector);

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddConstitutiveTensor(ElasticDataType& rVariables);

    virtual void CalculateAndAddIsochoricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    virtual void CalculateAndAddVolumetricConstitutiveTensor(ElasticDataType& rVariables, Matrix& rConstitutiveMatrix);

    //************//

    void InitializeElasticData(ModelDataType& rValues, HyperElasticDataType& rVariables);    

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticityModel )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticityModel )      
    }

     ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class LinearElasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LINEAR_ELASTIC_MODEL_H_INCLUDED  defined 


