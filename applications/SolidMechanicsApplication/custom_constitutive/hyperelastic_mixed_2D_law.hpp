//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_MIXED_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_MIXED_2D_LAW_H_INCLUDED

// System includes 

// External includes 

// Project includes
#include "custom_constitutive/hyperelastic_mixed_3D_law.hpp"


namespace Kratos
{
  /**
   * Defines a hyperelastic isotropic constitutive law in 3D Neohookean Model 
   * With stress split in an isochoric and volumetric parts
   * This material law is defined by the parameters:
   * 1) YOUNG MODULUS 
   * 2) POISSON RATIO
   * As there are no further parameters the functionality is limited 
   * to large displacements elasticity.
   */

  class HyperElasticMixed2DLaw : public HyperElasticMixed3DLaw
  {
  public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticMixed2DLaw
     */
    
    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticMixed2DLaw);
    
    /**
     * Life Cycle 
     */

    /**
     * Default constructor.
     */
    HyperElasticMixed2DLaw();
			
    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;
    
    /**
     * Copy constructor.
     */
    HyperElasticMixed2DLaw (const HyperElasticMixed2DLaw& rOther);
   

    /**
     * Assignment operator.
     */

    //HyperElasticMixed2DLaw& operator=(const HyperElasticMixed2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HyperElasticMixed2DLaw();
			
    /**
     * Operators 
     */
    
    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() { return 2; };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()         { return 3; };


     /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;
		
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


    /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
				       Vector& rStrainVector );


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
				 Vector& rStrainVector );
    

    /**
     * Calculates the isochoric constitutive matrix 
     * @param rMatrixIC can be the Identity or the inverse of the RightCauchyGreen tensor
     * @param rIsoStressVector the isochoric stress vector
     * @param rdetF the determinant of the deformation gradient
     * @param rTrace the trace of the RightCauchyGreen of the LeftCauchyGreen corresponding with the MatrixIC
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    void CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
					       const Vector & rIsoStressVector,
					       const double & rdetF0,
					       const double & rTrace,
					       const double & rLameLambda,
					       const double & rLameMu,
					       Matrix& rConstitutiveMatrix);

    /**
     * Calculates the volumetric constitutive matrix 
     * @param rMatrixIC can be the Identity or the inverse of the RightCauchyGreen tensor
     * @param rdetF the determinant of the deformation gradient
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */

    void CalculateVolumetricConstitutiveMatrix ( const Matrix & rMatrixIC,
						 const double & rdetF0,
						 const double & rLameLambda,
						 const double & rLameMu,
						 const GeometryType& rDomainGeometry,
						 const Vector & rShapeFunctions,
						 Matrix& rConstitutiveMatrix);
    
    /**
     * Calculates the isochoric constitutive matrix and makes a push-forward
     * @param rMatrixIC can be the Identity or the inverse of the RightCauchyGreen tensor
     * @param rIsoStressVector the isochoric stress vector
     * @param rF the total deformation gradient from the original to the current configuration
     * @param rdetF0 the determinant of the total deformation gradient
     * @param rTrace the trace of the RightCauchyGreen of the LeftCauchyGreen corresponding with the MatrixIC
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    void CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
					       const Vector & rIsoStressVector,
					       const Matrix & rF, 
					       const double & rdetF0,
					       const double & rTrace,
					       const double & rLameLambda,
					       const double & rLameMu,
					       Matrix& rConstitutiveMatrix);
    

    /**
     * Calculates the volumetric constitutive matrix and makes a push-forward
     * @param rMatrixIC can be the Identity or the inverse of the RightCauchyGreen tensor
     * @param rF the total deformation gradient from the original to the current configuration
     * @param rdetF0 the determinant of the total deformation gradient
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */

    void CalculateVolumetricConstitutiveMatrix ( const Matrix & rMatrixIC,
						 const Matrix & rF,			
						 const double & rdetF0,
						 const double & rLameLambda,
						 const double & rLameMu,
						 const GeometryType& rDomainGeometry,
						 const Vector & rShapeFunctions,
						 Matrix& rConstitutiveMatrix);  
		
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

 
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElasticMixed3DLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElasticMixed3DLaw);
    }



  }; // Class HyperElasticMixed2DLaw 
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_MIXED_2D_LAW_H_INCLUDED  defined 
