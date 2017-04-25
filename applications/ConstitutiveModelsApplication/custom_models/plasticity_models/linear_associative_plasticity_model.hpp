//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_LINEAR_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_LINEAR_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_linear_associative_plasticity_model.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) LinearAssociativePlasticityModel : public NonLinearAssociativePlasticityModel<TElasticityModel,TYieldCriterion>
  {
  public:
    
    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef TElasticityModel                               ElasticityModelType;
    typedef typename  ElasticityModelType::Pointer      ElasticityModelPointer;

    //yield criterion
    typedef TYieldCriterion                                 YieldCriterionType;
    typedef typename YieldCriterionType::Pointer         YieldCriterionPointer;

    //derived type
    typedef NonLinearAssociativePlasticityModel<ElasticityModelType,YieldCriterionType>   DerivedType;
    
    //base type
    typedef PlasticityModel<ElasticityModelType,YieldCriterionType>   BaseType;

    //common types
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::PlasticDataType                 PlasticDataType;
    typedef typename BaseType::InternalVariablesType     InternalVariablesType;
    
    /// Pointer definition of LinearAssociativePlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( LinearAssociativePlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearAssociativePlasticityModel() : DerivedType() {}

    /// Constructor.
    LinearAssociativePlasticityModel(ElasticityModelPointer pElasticityModel, YieldCriterionPointer pYieldCriterion) : DerivedType(pElasticityModel, pYieldCriterion) {}
    
    /// Copy constructor.
    LinearAssociativePlasticityModel(LinearAssociativePlasticityModel const& rOther) : DerivedType(rOther) {}

    /// Assignment operator.
    LinearAssociativePlasticityModel& operator=(LinearAssociativePlasticityModel const& rOther)
    {
      DerivedType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual ConstitutiveModel::Pointer Clone() const override
    {
      return ( LinearAssociativePlasticityModel::Pointer(new LinearAssociativePlasticityModel(*this)) );
    }

    /// Destructor.
    virtual ~LinearAssociativePlasticityModel() {}


    ///@}
    ///@name Operators
    ///@{

   
    ///@}
    ///@name Operations
    ///@{


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
      buffer << "LinearAssociativePlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "LinearAssociativePlasticityModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "LinearAssociativePlasticityModel Data";	    
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

      
    // calculate ratial return
    virtual bool CalculateRadialReturn(PlasticDataType& rVariables, MatrixType& rStressMatrix)
    {
      KRATOS_TRY

      //start
      double DeltaStateFunction = 0;

      double& rEquivalentPlasticStrain     = rVariables.Internal.Variables[0];
      double& rDeltaGamma                  = rVariables.DeltaInternal.Variables[0];


      double StateFunction                = rVariables.TrialStateFunction;	

      //Calculate Delta State Function:
      DeltaStateFunction = this->mpYieldCriterion->CalculateDeltaStateFunction( rVariables, DeltaStateFunction );

      //Calculate DeltaGamma:
      rDeltaGamma  = StateFunction/DeltaStateFunction;
	       
      //Update Equivalent Plastic Strain:
      DeltaPlasticStrain        = sqrt(2.0/3.0) * rDeltaGamma;
      rEquivalentPlasticStrain += DeltaPlasticStrain;
	
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )  
    }
    
    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class LinearAssociativePlasticityModel

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  ///@}
  
  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_LINEAR_ASSOCIATIVE_PLASTICITY_MODEL_H_INCLUDED  defined 


