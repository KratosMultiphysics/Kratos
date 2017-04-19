//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_HYPERELASTIC_AXISYMMETRIC_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_AXISYMMETRIC_2D_LAW_H_INCLUDED

// System includes

// External includes 

// Project includes
#include "custom_laws/hyperelastic_laws/hyperelastic_3D_law.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) HyperElasticAxisymmetric2DLaw : public HyperElastic3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{
     
      /// Pointer definition of HyperElasticAxisymmetric2DLaw
      KRATOS_CLASS_POINTER_DEFINITION(HyperElasticAxisymmetric2DLaw);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      HyperElasticAxisymmetric2DLaw() : HyperElastic3DLaw() {}

      /// Constructor.
      HyperElasticAxisymmetric2DLaw(ModelType::Pointer pModel) : HyperElastic3DLaw(pModel) {} 

      /// Copy constructor.
      HyperElasticAxisymmetric2DLaw(const HyperElasticAxisymmetric2DLaw& rOther) : HyperElastic3DLaw(rOther) {}

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
	return (HyperElasticAxisymmetric2DLaw::Pointer(new HyperElasticAxisymmetric2DLaw(*this)));
      }
      
      /// Destructor.
      virtual ~HyperElasticAxisymmetric2DLaw(){}
      

      ///@}
      ///@name Operators 
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 2; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 4; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures) override
      {
	KRATOS_TRY
	
	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

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
        buffer << "HyperElasticAxisymmetric2DLaw" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "HyperElasticAxisymmetric2DLaw";}

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
        
      /**
       * Get voigt index tensor:
       */
      VoigtIndexType GetVoigtIndexTensor() override 	
      {
	return this->msIndexVoigt2D4C;
      }
      
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )
      }
      
      virtual void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )
      }

      
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      HyperElasticAxisymmetric2DLaw& operator=(HyperElasticAxisymmetric2DLaw const& rOther){ return *this; }

        
      ///@}    
        
    }; // Class HyperElasticAxisymmetric2DLaw 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_HYPERELASTIC_AXISYMMETRIC_2D_LAW_H_INCLUDED  defined
