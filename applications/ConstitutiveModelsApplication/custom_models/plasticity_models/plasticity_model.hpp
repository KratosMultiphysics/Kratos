//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/elasticity_models/elasticity_model.hpp"
#include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp"

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
  template<class TElasticityModel, class TYieldCriterion>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) PlasticityModel
  {
  public:
    
    ///@name Type Definitions
    ///@{
    typedef ConstitutiveModelData::SizeType                               SizeType;
    typedef ConstitutiveModelData::VoigtIndexType                   VoigtIndexType;
    typedef ConstitutiveModelData::MatrixType                           MatrixType;
    typedef ConstitutiveModelData::VectorType                           VectorType;
    typedef ConstitutiveModelData::ModelData                         ModelDataType;

    typedef TElasticityModel                                   ElasticityModelType;
    typedef TYieldCriterion                                     YieldCriterionType;
    typedef typename TYieldCriterion::PlasticDataType              PlasticDataType;

    typedef typename TElasticityModel::Pointer          ElasticityModelTypePointer;
    typedef typename TYieldCriterion::Pointer            YieldCriterionTypePointer;

    typedef typename TYieldCriterion::InternalVariablesType  InternalVariablesType;
    
    /// Pointer definition of PlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( PlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlasticityModel() {}


    /// Constructor.
    PlasticityModel(ElasticityModelTypePointer pElasticityModel, YieldCriterionTypePointer pYieldCriterion) : mpElasticityModel(pElasticityModel), mpYieldCriterion(pYieldCriterion){}

    
    /// Copy constructor.
    PlasticityModel(PlasticityModel const& rOther) : mpElasticityModel(rOther.mpElasticityModel), mpYieldCriterion(rOther.mpYieldCriterion) {}

    /// Assignment operator.
    PlasticityModel& operator=(PlasticityModel const& rOther) {return *this;}

    /// Clone.
    PlasticityModel<TElasticityModel,TYieldCriterion>::Pointer Clone() const
    {
      return (PlasticityModel<TElasticityModel,TYieldCriterion>::Pointer(new PlasticityModel(*this)));
    }
    
    /// Destructor.
    virtual ~PlasticityModel() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Stresses
     */

    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_TRY
	
      KRATOS_THROW_ERROR( std::logic_error, "calling the PlasticityModel base class ... illegal operation!!", "" )
	
      KRATOS_CATCH(" ")
    }

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_TRY
	
      KRATOS_THROW_ERROR( std::logic_error, "calling the PlasticityModel base class ... illegal operation!!", "" )
	
      KRATOS_CATCH(" ")
    }

    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix)
    {
      KRATOS_TRY
	
      mpElasticityModel->CalculateVolumetricStress(rValues,rStressMatrix);
	
      KRATOS_CATCH(" ")
    }

    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) 
    {
      KRATOS_THROW_ERROR( std::logic_error, "calling PlasticityModel base class ..", "" )
    }
    
    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) 
    {
      KRATOS_THROW_ERROR( std::logic_error, "calling PlasticityModel base class ..", "" )
    }
    
    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) 
    {
      KRATOS_TRY
	
      mpElasticityModel->CalculateVolumetricConstitutiveTensor(rValues,rConstitutiveMatrix);
	
      KRATOS_CATCH(" ")
    }

    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
    {
      KRATOS_THROW_ERROR( std::logic_error, "calling PlasticityModel base class ..", "" )
    }
    
    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
    {
      KRATOS_THROW_ERROR( std::logic_error, "calling PlasticityModel base class ..", "" )
    }
    
    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix)
    {
      KRATOS_THROW_ERROR( std::logic_error, "calling PlasticityModel base class ..", "" )
    }



    ///@}
    ///@name Access
    ///@{

    ElasticityModel::Pointer pGetElasticityModel() {return mpElasticityModel;};
    
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
      buffer << "PlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PlasticityModel";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
    
    ElasticityModelType   mpElasticityModel;
    YieldCriterionType    mpYieldCriterion;
    
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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Serialization
    ///@{    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
      rSerializer.save("mpElasticityModel",mpElasticityModel);
      rSerializer.save("mpYieldCriterion",mpYieldCriterion);
    }

    virtual void load(Serializer& rSerializer)
    {
      rSerializer.load("mpElasticityModel",mpElasticityModel);
      rSerializer.load("mpYieldCriterion",mpYieldCriterion);
    }


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class PlasticityModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@} 
  ///@name Input and output 
  ///@{

  
  ///@}

  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_PLASTICITY_MODEL_H_INCLUDED  defined 


