//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_3D_LAW_H_INCLUDED )
#define  KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_3D_LAW_H_INCLUDED

// System includes

// External includes 

// Project includes
#include "custom_laws/linear_elastic_laws/linear_elastic_3D_law.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) LinearElasticOrthotropic3DLaw : public LinearElastic3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

     /// Pointer definition of LinearElasticOrthotropic3DLaw
      KRATOS_CLASS_POINTER_DEFINITION(LinearElasticOrthotropic3DLaw);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      LinearElasticOrthotropic3DLaw() : LinearElastic3DLaw() {}

      /// Copy constructor.
      LinearElasticOrthotropic3DLaw(const LinearElasticOrthotropic3DLaw& rOther) : LinearElastic3DLaw(rOther) {}

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
	return (LinearElasticOrthotropic3DLaw::Pointer(new LinearElasticOrthotropic3DLaw(*this)));
      }
      
      /// Destructor.
      virtual ~LinearElasticOrthotropic3DLaw(){}
      

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
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ANISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
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

      
      /**
       * This function is designed to be called once to perform all the checks needed
       * on the input provided. Checks can be "expensive" as the function is designed
       * to catch user's errors.
       * @param rMaterialProperties
       * @param rElementGeometry
       * @param rCurrentProcessInfo
       * @return
       */
      int Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo) override
      {
  
	if(YOUNG_MODULUS_X.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_X))
	  KRATOS_ERROR << "YOUNG_MODULUS_X has Key zero or invalid value" << std::endl;

	if(YOUNG_MODULUS_Y.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Y))
	  KRATOS_ERROR << "YOUNG_MODULUS_Y has Key zero or invalid value" << std::endl;
	
	if(YOUNG_MODULUS_Z.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Z))
	  KRATOS_ERROR << "YOUNG_MODULUS_Z has Key zero or invalid value" << std::endl;

	if(POISSON_RATIO_XY.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XY))
	  KRATOS_ERROR << "POISSON_RATIO_XY has Key zero invalid value" << std::endl;

	if(POISSON_RATIO_YZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_YZ))
	  KRATOS_ERROR << "POISSON_RATIO_YZ has Key zero invalid value" << std::endl;

	if(POISSON_RATIO_XZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XZ))
	  KRATOS_ERROR << "POISSON_RATIO_XZ has Key zero invalid value" << std::endl;

        if(DENSITY.Key() == 0 || !rMaterialProperties.Has(DENSITY))
	  KRATOS_ERROR << "DENSITY has Key zero or invalid value" << std::endl;

        return 0;
	
      }
      
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const override
      {
	std::stringstream buffer;
        buffer << "LinearElasticOrthotropic3DLaw" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "LinearElasticOrthotropic3DLaw";}

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


      /**
       * calculates the linear elastic constitutive matrix in terms of Young's modulus and
       * Poisson ratio
       * @param E the Young's modulus
       * @param NU the Poisson ratio
       * @return the linear elastic constitutive matrix
       */
      void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
					 const Properties& rMaterialProperties) override
      {
	KRATOS_TRY
  
	// Orthotropic constitutive matrix
	double E1 = rMaterialProperties[YOUNG_MODULUS_X];
	double E2 = rMaterialProperties[YOUNG_MODULUS_Y];
	double E3 = rMaterialProperties[YOUNG_MODULUS_Z];

	double v12 = rMaterialProperties[POISSON_RATIO_XY];
	double v23 = rMaterialProperties[POISSON_RATIO_YZ];
	double v13 = rMaterialProperties[POISSON_RATIO_XZ];

	double P1 = 1.0/(E2*E2*v12*v12 + 2.0*E3*E2*v12*v13*v23 + E3*E2*v13*v13 - E1*E2 + E1*E3*v23*v23);
	double P2 = E1*E1;
	double P3 = E2*E2;
	double P4 = E1*v23 + E2*v12*v13;
	double P5 = E2*v12 + E3*v13*v23;
	double P6 = E3*E3;

	rConstitutiveMatrix(0, 0) = -P1*P2*(- E3*v23*v23 + E2);
	rConstitutiveMatrix(0, 1) = -E1*E2*P1*P5;
	rConstitutiveMatrix(0, 2) = -E2*E3*P1*(E1*v13 + E1*v12*v23);
	rConstitutiveMatrix(1, 0) = -E1*E2*P1*P5;
	rConstitutiveMatrix(1, 1) = -P1*P3*(- E3*v13*v13 + E1);
	rConstitutiveMatrix(1, 2) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 0) = -E1*E2*E3*P1*(v13 + v12*v23);
	rConstitutiveMatrix(2, 1) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 2) = -E2*E3*P1*(- E2*v12*v12 + E1);
	rConstitutiveMatrix(3, 3) = (E2*P2)/(P2 + v12*(P2 + P3) + E1*E2)/2.0;
	rConstitutiveMatrix(4, 4) = (E3*P3)/(P3 + v23*(P3 + P6) + E2*E3)/2.0;
	rConstitutiveMatrix(5, 5) = (E3*P2)/(P2 + v13*(P2 + P6) + E1*E3)/2.0;
    
	KRATOS_CATCH(" ")
      }    
      
        
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
      }
      
      virtual void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
      }

      
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      LinearElasticOrthotropic3DLaw& operator=(LinearElasticOrthotropic3DLaw const& rOther){ return *this; }

        
      ///@}    
        
    }; // Class LinearElasticOrthotropic3DLaw 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_LINEAR_ELASTIC_ORTHOTROPIC_3D_LAW_H_INCLUDED  defined
