//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_J2_ISOCHORIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED)
#define  KRATOS_J2_ISOCHORIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED

// System includes

// External includes 

// Project includes
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/plasticity_models/non_linear_associative_plastic_model.hpp"
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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) J2IsochoricNeoHookean3DLaw : public HyperElasticPlastic3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

      typedef IsochoricNeoHookeanModel                                                      ElasticityModelType;
      typedef NonLinearIsotropicKinematicHardeningLaw                                          HardeningLawType;
      typedef MisesHuberYieldCriterion<HardeningLawType>                                     YieldCriterionType;
      typedef NonLinearAssociativePlasticModel<ElasticityModelType, YieldCriterionType>        PlasticModelType;
      typedef typename PlasticModelType::Pointer                                            PlasticModelPointer;
      
      /// Pointer definition of J2IsochoricNeoHookean3DLaw
      KRATOS_CLASS_POINTER_DEFINITION(J2IsochoricNeoHookean3DLaw);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      J2IsochoricNeoHookean3DLaw() : HyperElasticPlastic3DLaw()
      {
	KRATOS_TRY    
	  
	mpModel = PlasticModelPointer( new PlasticModelType() );

	KRATOS_CATCH(" ")		
      }

      /// Constructor.
      //J2IsochoricNeoHookean3DLaw(ModelType::Pointer pModel) : HyperElasticPlastic3DLaw(pModel) {} 

      /// Copy constructor.
      J2IsochoricNeoHookean3DLaw(const J2IsochoricNeoHookean3DLaw& rOther) : HyperElasticPlastic3DLaw(rOther) {}

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
	return (J2IsochoricNeoHookean3DLaw::Pointer(new J2IsochoricNeoHookean3DLaw(*this)));
      }
      
      /// Destructor.
      virtual ~J2IsochoricNeoHookean3DLaw(){}
      

      ///@}
      ///@name Operators 
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 3; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 6; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures) override
      {
	KRATOS_TRY
	  
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );
	rFeatures.mOptions.Set( U_P_LAW );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
	
	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

	KRATOS_CATCH(" ")
      }
      
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
        buffer << "J2IsochoricNeoHookean3DLaw" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "J2IsochoricNeoHookean3DLaw";}

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
      }
      
      virtual void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
      }

      
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      J2IsochoricNeoHookean3DLaw& operator=(J2IsochoricNeoHookean3DLaw const& rOther){ return *this; }

        
      ///@}    
        
    }; // Class J2IsochoricNeoHookean3DLaw 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_J2_ISOCHORIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED  defined
