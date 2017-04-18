//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_MODEL_H_INCLUDED )
#define  KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_MODEL_H_INCLUDED


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
  template<class TElasticityModel, class TYieldCriterion>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonLinearAssociativePlasticModel : public PlasticityModel<TElasticityModel,TYieldCriterion>
  {
  public:
    
    ///@name Type Definitions
    ///@{
    typedef PlasticityModel<TElasticityModel,TYieldCriterion>         BaseType;
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::ElasticityModelType         ElasticityModelType;
    typedef typename BaseType::YieldCriterionType           YieldCriterionType;

    typedef typename BaseType::ElasticityModelTypePointer   ElasticityModelTypePointer;
    typedef typename BaseType::YieldCriterionTypePointer     YieldCriterionTypePointer;

    typedef typename BaseType::PlasticDataType                PlasticDataType;
    typedef typename BaseType::InternalVariablesType    InternalVariablesType;
    
  protected:

    struct ThermalVariables
    {
    public:
      double PlasticDissipation;
      double DeltaPlasticDissipation;

    private:

      friend class Serializer;

      // A private default constructor necessary for serialization
      
      void save(Serializer& rSerializer) const
      {
	rSerializer.save("PlasticDissipation",PlasticDissipation);
	rSerializer.save("DeltaPlasticDissipation",DeltaPlasticDissipation);
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("PlasticDissipation",PlasticDissipation);
	rSerializer.load("DeltaPlasticDissipation",DeltaPlasticDissipation);
      };
      
    };


    struct PlasticFactors
    {      
      double Beta0;
      double Beta1;
      double Beta2;
      double Beta3;
      double Beta4;	   

      MatrixType  Normal;
      MatrixType  Dev_Normal;
    };

     	
  public:

    
    /// Pointer definition of NonLinearAssociativePlasticModel
    KRATOS_CLASS_POINTER_DEFINITION( NonLinearAssociativePlasticModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonLinearAssociativePlasticModel();

    /// Constructor.
    NonLinearAssociativePlasticModel(ElasticityModelTypePointer pElasticityModel, YieldCriterionTypePointer pYieldCriterion);
    
    /// Copy constructor.
    NonLinearAssociativePlasticModel(NonLinearAssociativePlasticModel const& rOther);

    /// Assignment operator.
    NonLinearAssociativePlasticModel& operator=(NonLinearAssociativePlasticModel const& rOther);

    /// Clone.
    virtual BaseTypePointer Clone() const override;
    
    /// Destructor.
    virtual ~NonLinearAssociativePlasticModel();


    ///@}
    ///@name Operators
    ///@{

   
    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Stresses
     */

    virtual void CalculateStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;

    virtual void CalculateIsochoricStressTensor(ModelDataType& rValues, MatrixType& rStressMatrix) override;
    
    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateConstitutiveTensor(ModelDataType& rValues, Matrix& rConstitutiveMatrix) override; 
       
    /**
     * Calculate Stress and Constitutive Tensor
     */
    virtual void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override;
    

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
      buffer << "NonLinearAssociativePlasticModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "NonLinearAssociativePlasticModel";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

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
    InternalVariablesType  mPreviousInternal;

    ThermalVariables mThermalVariables;
	
    ///@}
    ///@name Protected Operators
    ///@{

    
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculate Stresses
     */
    virtual void CalculateAndAddIsochoricStressTensor(PlasticDataType& rVariables, MatrixType& rStressMatrix);

    /**
     * Calculate Constitutive Tensor
     */
    virtual void CalculateAndAddPlasticConstitutiveTensor(PlasticDataType& rVariables, Matrix& rConstitutiveMatrix);
    
    /**
     * Calculate Constitutive Components
     */    
    virtual double& AddPlasticConstitutiveComponent(PlasticDataType& rVariables,
						    PlasticFactors& rFactors, double& rCabcd,
						    const unsigned int& a, const unsigned int& b,
						    const unsigned int& c, const unsigned int& d);
    
    // calculate ratial return
    
    virtual bool CalculateRadialReturn        (PlasticDataType& rVariables, MatrixType& rStressMatrix);

    // implex protected methods
    virtual void CalculateImplexRadialReturn  (PlasticDataType& rVariables, MatrixType& rStressMatrix);


    // auxiliar methods

    virtual void InitializeVariables          (ModelDataType& rValues, PlasticDataType& rVariables);
     
    virtual void UpdateStressConfiguration    (PlasticDataType& rVariables, MatrixType& rStressMatrix);	 

    virtual void UpdateInternalVariables      (PlasticDataType& rVariables);
    
    virtual void CalculateScalingFactors      (PlasticDataType& rVariables, PlasticFactors& rScalingFactors);

    
    // energy calculation methods

    void CalculateThermalDissipation          (PlasticDataType& rVariables);

    // implex protected methods
    void CalculateImplexThermalDissipation    (PlasticDataType& rVariables);

    
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

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      rSerializer.save("InternalVariables",mInternal);
      rSerializer.save("PreviousInternalVariables",mPreviousInternal);
      rSerializer.save("ThermalVariables",mThermalVariables);
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("InternalVariables",mInternal);
      rSerializer.load("PreviousInternalVariables",mPreviousInternal);
      rSerializer.load("ThermalVariables",mThermalVariables);      
    }
    
    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NonLinearAssociativePlasticModel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}
  
  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_MODEL_H_INCLUDED  defined 


