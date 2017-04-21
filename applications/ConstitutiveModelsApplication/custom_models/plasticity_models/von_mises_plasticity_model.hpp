//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_VON_MISES_PLASTICITY_MODEL_H_INCLUDED )
#define  KRATOS_VON_MISES_PLASTICITY_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/non_linear_associative_plasticity_model.hpp"
#include "custom_models/plasticity_models/yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_models/plasticity_models/hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"


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
  template<class TElasticityModel>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) VonMisesPlasticityModel : public NonLinearAssociativePlasticityModel<TElasticityModel, MisesHuberYieldCriterion<NonLinearIsotropicKinematicHardeningLaw> >
  {
  public:
    
    ///@name Type Definitions
    ///@{

    //elasticity model
    typedef TElasticityModel                               ElasticityModelType;
    typedef typename ElasticityModelType::Pointer       ElasticityModelPointer;

    //yield criterion
    typedef NonLinearIsotropicKinematicHardeningLaw                    HardeningLawType;
    typedef MisesHuberYieldCriterion<HardeningLawType>               YieldCriterionType;
    typedef typename YieldCriterionType::Pointer                  YieldCriterionPointer;

    //base type
    typedef NonLinearAssociativePlasticityModel<TElasticityModel,YieldCriterionType>  BaseType;

    //common types
    typedef typename BaseType::Pointer                         BaseTypePointer;
    typedef typename BaseType::SizeType                               SizeType;
    typedef typename BaseType::VoigtIndexType                   VoigtIndexType;
    typedef typename BaseType::MatrixType                           MatrixType;
    typedef typename BaseType::ModelDataType                     ModelDataType;
    typedef typename BaseType::MaterialDataType               MaterialDataType;
    typedef typename BaseType::PlasticDataType                 PlasticDataType;
    typedef typename BaseType::InternalVariablesType     InternalVariablesType;


    /// Pointer definition of VonMisesPlasticityModel
    KRATOS_CLASS_POINTER_DEFINITION( VonMisesPlasticityModel );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VonMisesPlasticityModel() : BaseType() {}

    /// Constructor.
    VonMisesPlasticityModel(ElasticityModelPointer pElasticityModel)
      :BaseType()
    {
      KRATOS_TRY

      this->mpElasticityModel = ElasticityModelPointer( new ElasticityModelType() );
      this->mpYieldCriterion  = YieldCriterionPointer( new YieldCriterionType() );

      std::cout<<" Von Mises Model "<<std::endl;
      this->mpElasticityModel->PrintInfo(std::cout);
     
      KRATOS_CATCH(" ")
    }
    
    /// Copy constructor.
    VonMisesPlasticityModel(VonMisesPlasticityModel const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    VonMisesPlasticityModel& operator=(VonMisesPlasticityModel const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    ElasticityModel::Pointer Clone() const override
    {
      return ( VonMisesPlasticityModel::Pointer(new VonMisesPlasticityModel(*this)) );
    }
    
    /// Destructor.
    virtual ~VonMisesPlasticityModel() {}


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
      buffer << "VonMisesPlasticityModel" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesPlasticityModel";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesPlasticityModel Data";
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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Serialization
    ///@{    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }
    
    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

  }; // Class VonMisesPlasticityModel

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

#endif // KRATOS_VON_MISES_PLASTICITY_MODEL_H_INCLUDED  defined 


