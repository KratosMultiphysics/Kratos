//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_HYPERELASTIC_PLASTIC_UP_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_PLASTIC_UP_3D_LAW_H_INCLUDED

// System includes

// External includes 

// Project includes
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HyperElasticPlasticUP3DLaw : public HyperElasticPlastic3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of HyperElasticPlasticUP3DLaw
      KRATOS_CLASS_POINTER_DEFINITION(HyperElasticPlasticUP3DLaw);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      HyperElasticPlasticUP3DLaw() : HyperElasticPlastic3DLaw() {}

      /// Constructor.
      HyperElasticPlasticUP3DLaw(ModelType::Pointer pModel) : HyperElasticPlastic3DLaw(pModel) {} 

      /// Copy constructor.
      HyperElasticPlasticUP3DLaw(const HyperElasticPlasticUP3DLaw& rOther) : HyperElasticPlastic3DLaw(rOther) {}

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
	return (HyperElasticPlasticUP3DLaw::Pointer(new HyperElasticPlasticUP3DLaw(*this)));
      }
      
      /// Destructor.
      virtual ~HyperElasticPlasticUP3DLaw(){}
      

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
        buffer << "HyperElasticPlasticUP3DLaw" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "HyperElasticPlasticUP3DLaw";}

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
      HyperElasticPlasticUP3DLaw& operator=(HyperElasticPlasticUP3DLaw const& rOther){ return *this; }

        
      ///@}    
        
    }; // Class HyperElasticPlasticUP3DLaw 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_HYPERELASTIC_PLASTIC_UP_3D_LAW_H_INCLUDED  defined
