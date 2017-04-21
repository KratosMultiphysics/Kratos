//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_ELASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_ELASTICITY_MODEL_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ElasticityModel
  {
  protected:

    using VoigtIndexType = const unsigned int(*)[2];

  public:

    ///@name Type Definitions
    ///@{  
    typedef ConstitutiveModelData::SizeType                    SizeType;
    typedef ConstitutiveModelData::VectorType                VectorType;
    typedef ConstitutiveModelData::MatrixType                MatrixType;
    typedef ConstitutiveModelData::ModelData              ModelDataType;
    typedef ConstitutiveModelData::MaterialData        MaterialDataType; 

    typedef ConstitutiveModelData::StrainMeasureType  StrainMeasureType;   
    typedef ConstitutiveModelData::StressMeasureType  StressMeasureType;   
    
    /// Pointer definition of ElasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( ElasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.    
    ElasticityModel() {}

    /// Copy constructor.
    ElasticityModel(ElasticityModel const& rOther) {}

    /// Clone.
    virtual ElasticityModel::Pointer Clone() const
    {
      return ( ElasticityModel::Pointer(new ElasticityModel(*this)) );
    }

    /// Assignment operator.
    ElasticityModel& operator=(ElasticityModel const& rOther) {return *this;}


    /// Destructor.
    virtual ~ElasticityModel(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    

    /**
     * Calculate Strain Energy Density Functions
     */
    virtual void CalculateStrainEnergy(ModelDataType& rValues, double& rDensityFunction)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
      
    /**
     * Calculate Stresses
     */    
    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }

    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }

    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
    }

    
    /**
     * Check
     */    
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_ERROR << "calling ElasticityModel base class " << std::endl;
      return 0;
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
        buffer << "ElasticityModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ElasticityModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "ElasticityModel Data";
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

  }; // Class ElasticityModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ELASTICITY_MODEL_H_INCLUDED  defined 


