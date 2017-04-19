//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_VON_MISES_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_VON_MISES_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/plasticity_models/von_mises_plasticity_model.hpp"

namespace Kratos
{
  /**
   * Defines a hyperelastic-plastic isotropic constitutive law in 3D with von Mises plasticity model
   * the functionality is limited to large displacements plasticity.
   */
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) VonMisesHyperElasticPlastic3DLaw : public HyperElasticPlastic3DLaw
  {
  public:

    ///@name Type Definitions
    ///@{      

    typedef IsochoricNeoHookeanModel                                      HyperElasticModelType;
    typedef VonMisesPlasticityModel<HyperElasticModelType>                            ModelType;
    typedef typename ModelType::Pointer                                        ModelTypePointer;

    //base type
    typedef HyperElasticPlastic3DLaw                                                   BaseType;
      
    /// Pointer definition of VonMisesHyperElasticPlastic3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( VonMisesHyperElasticPlastic3DLaw );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VonMisesHyperElasticPlastic3DLaw() : BaseType()
    {
     KRATOS_TRY

     mpModel = ModelTypePointer( new ModelType() );
       
     KRATOS_CATCH(" ")	  
    }

    /// Constructor.
    //VonMisesHyperElasticPlastic3DLaw(ModelTypePointer pModel) :  BaseType(pModel) {} 

    /// Copy constructor.
    VonMisesHyperElasticPlastic3DLaw (const VonMisesHyperElasticPlastic3DLaw& rOther) : BaseType(rOther) {}

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override
    {
      return ( VonMisesHyperElasticPlastic3DLaw::Pointer(new VonMisesHyperElasticPlastic3DLaw(*this)) );
    }
    
    /// Destructor.
    virtual ~VonMisesHyperElasticPlastic3DLaw() {}

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
      buffer << "VonMisesHyperElasticPlastic3DLaw";
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesHyperElasticPlastic3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "VonMisesHyperElasticPlastic3DLaw Data";
      BaseType::PrintData(rOStream);
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

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )    	
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    VonMisesHyperElasticPlastic3DLaw& operator=(VonMisesHyperElasticPlastic3DLaw const& rOther){ return *this; }


    ///@}
  }; // Class VonMisesHyperElasticPlastic3DLaw
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.
#endif // KRATOS_VON_MISES_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
