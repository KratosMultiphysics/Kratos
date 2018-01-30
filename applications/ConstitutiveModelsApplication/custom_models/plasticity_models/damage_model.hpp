//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_DAMAGE_MODEL_H_INCLUDED )
#define  KRATOS_DAMAGE_MODEL_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/plasticity_model.hpp"

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
  template<class TElasticityModel, class TYieldSurface>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) DamageModel : public PlasticityModel<TElasticityModel,TYieldSurface>
  {
  public:
    
    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef TElasticityModel                               ElasticityModelType;

    //yield surface
    typedef TYieldSurface                                     YieldSurfaceType;
 
    //base type
    typedef PlasticityModel<ElasticityModelType,YieldSurfaceType>     BaseType;

    //common types
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::VectorType                           VectorType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::PlasticDataType                 PlasticDataType;
    typedef typename BaseType::InternalVariablesType     InternalVariablesType;

    typedef ConstitutiveModelData::StrainMeasureType         StrainMeasureType;   
    typedef ConstitutiveModelData::StressMeasureType         StressMeasureType;   

    
    /// Pointer definition of DamageModel
    KRATOS_CLASS_POINTER_DEFINITION( DamageModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DamageModel() : BaseType() {}
   
    /// Copy constructor.
    DamageModel(DamageModel const& rOther) :BaseType(rOther), mInternal(rOther.mInternal) {}

    /// Assignment operator.
    DamageModel& operator=(DamageModel const& rOther)
    {
      BaseType::operator=(rOther);
      mInternal = rOther.mInternal;
      return *this;
    }

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override
    {
      return ( DamageModel::Pointer(new DamageModel(*this)) );
    }

    /// Destructor.
    virtual ~DamageModel() {}


    ///@}
    ///@name Operators
    ///@{

   
    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize member data
     */    
    void InitializeModel(ModelDataType& rValues) override
    {
      KRATOS_TRY

      BaseType::InitializeModel(rValues);

      double& rDamageThreshold  = mInternal.Variables[0];


      //damage threshold properties
      const Properties& rMaterialProperties  =  rValues.GetMaterialProperties();
      rDamageThreshold =  rMaterialProperties[DAMAGE_THRESHOLD];


	
      KRATOS_CATCH(" ")
    }
    

    /**
     * Calculate Stresses
     */

    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override
    {      
      KRATOS_TRY
	
      PlasticDataType Variables;
      this->InitializeVariables(rValues,Variables);

      // calculate elastic stress
      this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);
      
      rValues.StressMatrix = rStressMatrix;  //store stress to ModelData StressMatrix

      // calculate damaged stress
      this->CalculateAndAddStressTensor(Variables,rStressMatrix);
      
      if( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES ) )
	this->UpdateInternalVariables(rValues, Variables, rStressMatrix);
     		
      KRATOS_CATCH(" ")	
    }
    
    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY

      //Initialize ConstitutiveMatrix
      rConstitutiveMatrix.clear();
      
      PlasticDataType Variables;
      this->InitializeVariables(rValues,Variables);

      //Calculate Stress Matrix

      MatrixType StressMatrix;
      // calculate elastic stress
      this->mElasticityModel.CalculateStressTensor(rValues,StressMatrix);
      
      rValues.StressMatrix = StressMatrix;  //store stress to ModelData StressMatrix

      // calculate damaged stress
      this->CalculateAndAddStressTensor(Variables,StressMatrix);

    
      //Calculate Constitutive Matrix

      // calculate elastic constitutive tensor
      this->mElasticityModel.CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);
    
      // calculate plastic constitutive tensor
      if( Variables.State().Is(ConstitutiveModelData::PLASTIC_REGION) )
      	this->CalculateAndAddPlasticConstitutiveTensor(Variables,rConstitutiveMatrix);
              
      Variables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
          
      KRATOS_CATCH(" ")
    }
    
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
    {
      KRATOS_TRY

      PlasticDataType Variables;
      this->InitializeVariables(rValues,Variables);

      //Calculate Stress Matrix
      
      // calculate elastic stress
      this->mElasticityModel.CalculateStressTensor(rValues,rStressMatrix);
      
      rValues.StressMatrix = rStressMatrix;  //store stress to ModelData StressMatrix

      // calculate damaged stress
      this->CalculateAndAddStressTensor(Variables,rStressMatrix);
            
      
      //Calculate Constitutive Matrix
    
      // calculate elastic constitutive tensor
      this->mElasticityModel.CalculateConstitutiveTensor(rValues,rConstitutiveMatrix);

      // calculate plastic constitutive tensor
      if( Variables.State().Is(ConstitutiveModelData::PLASTIC_REGION) ){
      	this->CalculateAndAddPlasticConstitutiveTensor(Variables,rConstitutiveMatrix);
      }

      Variables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
      
 
      if( rValues.State.Is(ConstitutiveModelData::UPDATE_INTERNAL_VARIABLES ) )
	this->UpdateInternalVariables(rValues, Variables, rStressMatrix);

      
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
      buffer << "DamageModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "DamageModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "DamageModel Data";	    
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
   
    // internal variables:
    InternalVariablesType  mInternal;
	
    ///@}
    ///@name Protected Operators
    ///@{

    
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculate Stresses
     */    
    virtual void CalculateAndAddStressTensor(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      double& rDamageThreshold = rVariables.Internal.Variables[0];

      //1.- Stress norm and Strengh factor for damage (StressNorm and RateFactor variables used)
      this->CalculateStressNorm(rVariables,rStressMatrix);
      
      //2.-Check yield condition
      rVariables.TrialStateFunction = this->mYieldSurface.CalculateYieldCondition(rVariables, rVariables.TrialStateFunction);

      if( rVariables.TrialStateFunction <= rDamageThreshold )
	{
	  rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,false);
	}
      else
	{
	  
	  //3.- Calculate the radial return
	  bool converged = this->CalculateRadialReturn(rVariables,rStressMatrix);
	  
	  if(!converged)
	    std::cout<<" ConstitutiveLaw did not converge "<<std::endl;
	  	  
	  //4.- Update back stress, plastic strain and stress
	  this->UpdateStressConfiguration(rVariables,rStressMatrix);
	  
	  
	  rVariables.State().Set(ConstitutiveModelData::PLASTIC_REGION,true);	    
	}
      
            
      rVariables.State().Set(ConstitutiveModelData::RETURN_MAPPING_COMPUTED,true);
   
      KRATOS_CATCH(" ")    
    }
    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddPlasticConstitutiveTensor(PlasticDataType& rVariables, Matrix& rConstitutiveMatrix)
    {
      KRATOS_TRY
    
      //Compute radial return check
      if( rVariables.State().IsNot(ConstitutiveModelData::RETURN_MAPPING_COMPUTED) )
	KRATOS_ERROR << "ReturnMapping has to be computed to perform the calculation" << std::endl;

      double& rDamage          = rVariables.DeltaInternal.Variables[0];

      //alternative way, compute damage if not computed
      //double rDamage = 0;
      //rDamage = this->mYieldSurface.CalculateStateFunction( rVariables, rDamage );
      
      double DeltaStateFunction = 0;
      DeltaStateFunction = this->mYieldSurface.CalculateDeltaStateFunction( rVariables, DeltaStateFunction );

      //tangent OPTION 1:
      
      // double& rDamageThreshold = rVariables.Internal.Variables[0];
      // const MatrixType& rStrainMatrix  = rVariables.GetStrainMatrix();

      // VectorType StrainVector;
      // StrainVector  = ConstitutiveModelUtilities::StrainTensorToVector(rStrainMatrix, StrainVector);

      // VectorType EffectiveStressVector;
      // noalias(EffectiveStressVector) = prod(rConstitutiveMatrix,StrainVector);

      
      // rConstitutiveMatrix *= (1-rDamage);
      // rConstitutiveMatrix += DeltaStateFunction * rDamageThreshold * outer_prod(EffectiveStressVector,EffectiveStressVector);
	
      //alternative tangent OPTION 2:
      
      const ModelDataType&  rModelData = rVariables.GetModelData();
      const MatrixType& rStressMatrix  = rModelData.GetStressMatrix();

      Vector EffectiveStressVector;
      EffectiveStressVector  = ConstitutiveModelUtilities::StressTensorToVector(rStressMatrix, EffectiveStressVector);

      VectorType EquivalentStrainVector;
      EquivalentStrainVector = this->CalculateEquivalentStrainDerivative(rVariables, rConstitutiveMatrix, EquivalentStrainVector);

      rConstitutiveMatrix *= (1-rDamage);
      rConstitutiveMatrix += DeltaStateFunction * outer_prod(EffectiveStressVector,EquivalentStrainVector);
      
      
      rVariables.State().Set(ConstitutiveModelData::CONSTITUTIVE_MATRIX_COMPUTED,true);
    
      KRATOS_CATCH(" ")
    }
    

    
    // calculate ratial return
    
    virtual bool CalculateRadialReturn(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      double& rDamageThreshold     = rVariables.Internal.Variables[0];
      double& rDamage              = rVariables.DeltaInternal.Variables[0];
	     
      double StateFunction         = rVariables.TrialStateFunction;
	
      if ( StateFunction >= rDamageThreshold )
	{
	  rDamageThreshold = StateFunction;
	       	
	  //Calculate State Function: (damage)
	  StateFunction = this->mYieldSurface.CalculateStateFunction( rVariables, StateFunction );

	  rDamage = StateFunction;

	  return true;

	}

      return false;

      KRATOS_CATCH(" ")
    }


    // auxiliar methods

    virtual void InitializeVariables(ModelDataType& rValues, PlasticDataType& rVariables)
    {
      KRATOS_TRY

      //set model data pointer
      rVariables.SetModelData(rValues);

      rValues.State.Set(ConstitutiveModelData::PLASTIC_REGION,false);

      rVariables.SetState(rValues.State);
    
      // RateFactor
      rVariables.RateFactor = 0;

      // Damage threshold variable
      rVariables.Internal = mInternal;

      // Damage variable
      rVariables.DeltaInternal.Variables.clear();

      // Flow Rule local variables
      rVariables.TrialStateFunction = 0;
      rVariables.StressNorm        = 0;

      // Set Strain
      rVariables.StrainMatrix = rValues.StrainMatrix;
      
      
      KRATOS_CATCH(" ")
    }
    
    virtual void UpdateStressConfiguration(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      double& rDamage = rVariables.DeltaInternal.Variables[0];
      
      //Stress Update: 
      rStressMatrix *= (1.0-rDamage);
       	      
      KRATOS_CATCH(" ")    
    }

    virtual void UpdateInternalVariables(ModelDataType& rValues, PlasticDataType& rVariables, const MatrixType& rStressMatrix)
    {
      KRATOS_TRY
      
      double& rDamageThreshold  = mInternal.Variables[0];

      //update mechanical variables
      rDamageThreshold = rVariables.Internal.Variables[0];     
      
      KRATOS_CATCH(" ")    
    }
    
    // calculate stress norm

    void CalculateStressNorm(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      // Compute strenght type parameter
      rVariables.RateFactor = 0.0;
        
      VectorType PrincipalStresses;
      noalias(PrincipalStresses) = ConstitutiveModelUtilities::EigenValuesDirectMethod(rStressMatrix);

      double Macaulay_PrincipalStress = 0.0;
      double Absolute_PrincipalStress = 0.0;
    
      for(unsigned int i=0; i<3; i++)
	{ 
	  if(PrincipalStresses[i] > 0.0)
	    {
	      Macaulay_PrincipalStress += PrincipalStresses[i];
	      Absolute_PrincipalStress += PrincipalStresses[i];
	    }
	  else
	    {
	      Absolute_PrincipalStress -= PrincipalStresses[i];
	    }
	}

      if(Absolute_PrincipalStress > 1.0e-20)
	{
	  rVariables.RateFactor = Macaulay_PrincipalStress/Absolute_PrincipalStress;
	}
      else
	{
	  rVariables.RateFactor = 0.5;
	}
    
      // Compute Equivalent Strain (rYieldCondition)
      const Matrix& rStrainMatrix = rVariables.GetStrainMatrix();
      MatrixType Auxiliar;
      noalias(Auxiliar) = prod(rStrainMatrix,rStressMatrix);
    
      rVariables.StressNorm = 0.0;
    
      for(unsigned int i=0; i<3; i++) 
	{
	  rVariables.StressNorm += Auxiliar(i,i);
	}
    
    
      KRATOS_CATCH(" ")   

    }

    
    // calculate equivalent strain derivative
    
    VectorType& CalculateEquivalentStrainDerivative(PlasticDataType& rVariables, const Matrix& rConstitutiveMatrix, VectorType& rEquivalentStrainDerivative)
    {
      KRATOS_TRY

      //The derivative of the equivalent strain with respect to the strain vector is obtained through the perturbation method

      VectorType StressVector;
      MatrixType StressMatrix;
      double EquivalentStrainForward  = 0.0;
      double EquivalentStrainBackward = 0.0;
      
      //Compute the strains perturbations in each direction of the vector
      const MatrixType& rStrainMatrix  = rVariables.GetStrainMatrix();
      VectorType StrainVector;
      StrainVector = ConstitutiveModelUtilities::StrainTensorToVector(rStrainMatrix, StrainVector);
      VectorType PerturbatedStrainVector;
      ConstitutiveModelUtilities::ComputePerturbationVector(PerturbatedStrainVector,StrainVector);

      for(unsigned int i = 0; i < StrainVector.size(); i++)
	{
	  //Forward perturbed equivalent strain
	  StrainVector[i] += PerturbatedStrainVector[i];
	  
	  rVariables.StrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(StrainVector,rVariables.StrainMatrix);
	  noalias(StressVector)   = prod(rConstitutiveMatrix, StrainVector);
	  StressMatrix            = ConstitutiveModelUtilities::StressVectorToTensor(StressVector,StressMatrix);

	  this->CalculateStressNorm(rVariables,StressMatrix);
	  EquivalentStrainForward = this->mYieldSurface.CalculateYieldCondition(rVariables, EquivalentStrainForward);

	  StrainVector[i] -= PerturbatedStrainVector[i];
	  
	  //Backward perturbed equivalent strain
	  StrainVector[i] -= PerturbatedStrainVector[i];
	  
	  rVariables.StrainMatrix = ConstitutiveModelUtilities::StrainVectorToTensor(StrainVector,rVariables.StrainMatrix);
	  noalias(StressVector)   = prod(rConstitutiveMatrix, StrainVector);
	  StressMatrix            = ConstitutiveModelUtilities::StressVectorToTensor(StressVector,StressMatrix);

	  this->CalculateStressNorm(rVariables,StressMatrix);
	  EquivalentStrainForward = this->mYieldSurface.CalculateYieldCondition(rVariables, EquivalentStrainForward);

	  StrainVector[i] += PerturbatedStrainVector[i];
        
	  rEquivalentStrainDerivative[i] = (EquivalentStrainForward - EquivalentStrainBackward) / (2.0 * PerturbatedStrainVector[i]);
	}

      return rEquivalentStrainDerivative;

      KRATOS_CATCH(" ")         
    }

    
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      rSerializer.save("InternalVariables",mInternal);
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("InternalVariables",mInternal);
    }
    
    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class DamageModel

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  ///@}
  
  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DAMAGE_MODEL_H_INCLUDED  defined 


