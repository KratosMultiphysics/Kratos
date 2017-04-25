//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/constitutive_model_utilities.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ConstitutiveModel
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
    
    /// Pointer definition of ConstitutiveModel
    KRATOS_CLASS_POINTER_DEFINITION( ConstitutiveModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.    
    ConstitutiveModel() {}

    /// Copy constructor.
    ConstitutiveModel(ConstitutiveModel const& rOther) {}

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const
    {
      return ( ConstitutiveModel::Pointer(new ConstitutiveModel(*this)) );
    }

    /// Assignment operator.
    ConstitutiveModel& operator=(ConstitutiveModel const& rOther) {return *this;}


    /// Destructor.
    virtual ~ConstitutiveModel(){}


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
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
      
    /**
     * Calculate Stresses
     */    
    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }

    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutive) 
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }

    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutive)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
    }

    
    /**
     * Check
     */    
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_ERROR << "calling ConstitutiveModel base class " << std::endl;
      return 0;
    }
    
    ///@}
    ///@name Access
    ///@{
        
    /**
     * method to ask the constituitve model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    virtual void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
					std::vector<Variable<array_1d<double,3> > >& rComponentVariables)
    {
      KRATOS_TRY
	
      KRATOS_ERROR << "calling the Constitutive Model base class Variables List... illegal operation" << std::endl;
	
      KRATOS_CATCH(" ")
    }
    
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
        buffer << "ConstitutiveModel";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ConstitutiveModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "ConstitutiveModel Data";
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

  }; // Class ConstitutiveModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSTITUTIVE_MODEL_H_INCLUDED  defined 


