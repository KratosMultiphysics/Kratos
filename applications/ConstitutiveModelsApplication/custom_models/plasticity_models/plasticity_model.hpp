//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/constitutive_model.hpp"
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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) PlasticityModel : public ConstitutiveModel
  {
  public:
    
    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef TElasticityModel                                   ElasticityModelType;
    typedef typename TElasticityModel::Pointer              ElasticityModelPointer;

    //yield criterion
    typedef TYieldCriterion                                     YieldCriterionType;
    typedef typename TYieldCriterion::Pointer                YieldCriterionPointer;

    //common types
    typedef ConstitutiveModelData::SizeType                               SizeType;
    typedef ConstitutiveModelData::VoigtIndexType                   VoigtIndexType;
    typedef ConstitutiveModelData::MatrixType                           MatrixType;
    typedef ConstitutiveModelData::VectorType                           VectorType;
    typedef ConstitutiveModelData::ModelData                         ModelDataType;
    typedef typename TYieldCriterion::PlasticDataType              PlasticDataType;
    typedef typename TYieldCriterion::InternalVariablesType  InternalVariablesType;

    
    /// Pointer definition of PlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( PlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlasticityModel() : ConstitutiveModel()
    {
      KRATOS_TRY
	
      mpElasticityModel = ElasticityModelPointer( new ElasticityModelType() );
      mpYieldCriterion  = YieldCriterionPointer( new YieldCriterionType() );
	
      KRATOS_CATCH(" ")
    }


    /// Constructor.
    PlasticityModel(ElasticityModelPointer pElasticityModel, YieldCriterionPointer pYieldCriterion) : ConstitutiveModel(), mpElasticityModel(pElasticityModel), mpYieldCriterion(pYieldCriterion) {}

    
    /// Copy constructor.
    PlasticityModel(PlasticityModel const& rOther) : ConstitutiveModel(rOther), mpElasticityModel(rOther.mpElasticityModel), mpYieldCriterion(rOther.mpYieldCriterion) {}

    /// Assignment operator.
    PlasticityModel& operator=(PlasticityModel const& rOther) {return *this;}

    /// Clone.
    ConstitutiveModel::Pointer Clone() const override
    {
      return ( PlasticityModel::Pointer(new PlasticityModel(*this)) );
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

    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY
	
      KRATOS_ERROR << "calling the PlasticityModel base class ... illegal operation" << std::endl;
	
      KRATOS_CATCH(" ")
    }

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY
	
      KRATOS_ERROR << "calling the PlasticityModel base class ... illegal operation" << std::endl;
	
      KRATOS_CATCH(" ")
    }

    virtual void CalculateVolumetricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {
      KRATOS_TRY
	
      mpElasticityModel->CalculateVolumetricStressTensor(rValues,rStressMatrix);
	
      KRATOS_CATCH(" ")
    }

    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_ERROR << "calling PlasticityModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_ERROR << "calling PlasticityModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY
	
      mpElasticityModel->CalculateVolumetricConstitutiveTensor(rValues,rConstitutiveMatrix);
	
      KRATOS_CATCH(" ")
    }

    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_ERROR << "calling PlasticityModel base class " << std::endl;
    }
    
    virtual void CalculateIsochoricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_ERROR << "calling PlasticityModel base class " << std::endl;
    }
    
    virtual void CalculateVolumetricStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_ERROR << "calling PlasticityModel base class " << std::endl;
    }


    /**
     * Check
     */    
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY
	
      return 0;

      KRATOS_CATCH(" ")
    }
    
    ///@}
    ///@name Access
    ///@{

    /**
     * Has Values
     */   
    virtual bool Has(const Variable<double>& rThisVariable) override {return false;}
    
    /**
     * Set Values
     */
    virtual void SetValue(const Variable<double>& rVariable,
                  const double& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override {}   
    /**
     * Get Values
     */
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) override { rValue=0; return rValue;}

    
    /**
     * method to ask the plasticity model the list of variables (dofs) needed from the domain
     * @param rScalarVariables : list of scalar dofs
     * @param rComponentVariables :  list of vector dofs
     */
    virtual void GetDomainVariablesList(std::vector<Variable<double> >& rScalarVariables,
					std::vector<Variable<array_1d<double,3> > >& rComponentVariables) override
    {
      KRATOS_TRY

      mpElasticityModel->GetDomainVariablesList(rScalarVariables, rComponentVariables);
 	
      KRATOS_CATCH(" ")
    }

    
    ConstitutiveModel::Pointer pGetElasticityModel() {return mpElasticityModel;};    
    
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
      buffer << "PlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "PlasticityModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "PlasticityModel Data";
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
    
    ElasticityModelPointer   mpElasticityModel;
    YieldCriterionPointer    mpYieldCriterion;
    
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

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveModel )

      rSerializer.save("mpElasticityModel",mpElasticityModel);
      rSerializer.save("mpYieldCriterion",mpYieldCriterion);
    }
    
    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveModel )

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


