//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_MIXED_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_MIXED_3D_LAW_H_INCLUDED

// System includes 

// External includes 

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"


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

  class HyperElasticMixed3DLaw : public HyperElastic3DLaw
  {
  protected:

    /**
     * Parameters to be used in the volumetric and deviatoric split
     */
    struct VectorSplit
    {
      Vector  Isochoric;
      Vector  Volumetric;
    };

    struct MatrixSplit
    {
      Matrix  Isochoric;
      Matrix  Volumetric;
    };


  public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticMixed3DLaw
     */
    
    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticMixed3DLaw);
    
    /**
     * Life Cycle 
     */

    /**
     * Default constructor.
     */
    HyperElasticMixed3DLaw();
			
    
    /**
     * Copy constructor.
     */
    HyperElasticMixed3DLaw (const HyperElasticMixed3DLaw& rOther);
   

    /**
     * Assignment operator.
     */

    //HyperElasticMixed3DLaw& operator=(const HyperElasticMixed3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HyperElasticMixed3DLaw();
			
    /**
     * Operators 
     */
    
    /**
     * Operations needed by the base class:
     */

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues 
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);

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
     * Calculates the Pressure of the domain (element)
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rPressure the calculated pressure to be returned
     */
    double& CalculateDomainPressure (const GeometryType& rDomainGeometry,
				     const Vector & rShapeFunctions, 
				     double & rPressure);

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
    virtual void CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
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

    virtual void CalculateVolumetricConstitutiveMatrix ( const Matrix & rMatrixIC,
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
    virtual void CalculateIsochoricConstitutiveMatrix (const Matrix & rMatrixIC,
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

    virtual void CalculateVolumetricConstitutiveMatrix ( const Matrix & rMatrixIC,
							 const Matrix & rF,			
							 const double & rdetF0,
							 const double & rLameLambda,
							 const double & rLameMu,
							 const GeometryType& rDomainGeometry,
							 const Vector & rShapeFunctions,
							 Matrix& rConstitutiveMatrix);
  


    /**
     * Isochoric constitutive component
     */
    
    double& IsochoricConstitutiveComponent(double & rCabcd, 
					   const Matrix & rMatrixIC,
					   const double & rdetF0, 
					   const double & rTrace,
					   const double & rLameLambda, 
					   const double & rLameMu, 
					   const Matrix & rIsoStressMatrix,
					   const unsigned int& a, const unsigned int& b,
					   const unsigned int& c, const unsigned int& d);
    
    /**
     * Volumetric constitutive component
     */
    
    double& VolumetricConstitutiveComponent(double & rCabcd,
					    const Matrix & rMatrixIC, 
					    const double & rdetF0, 
					    const double & rLameLambda, 
					    const double & rLameMu,
					    const double & rPressure,
					    const unsigned int& a, const unsigned int& b, 
					    const unsigned int& c, const unsigned int& d);
    
    /**
     * Isochoric constitutive component push-forward
     */
    double& IsochoricConstitutiveComponent(double & rCabcd,
					   const Matrix & rMatrixIC, 
					   const Matrix & rF,
					   const double & rdetF0, 
					   const double & rTrace,
					   const double & rLameLambda, 
					   const double & rLameMu,
					   const Matrix & rIsoStressMatrix,
					   const unsigned int& a, const unsigned int& b, 
					   const unsigned int& c, const unsigned int& d);

   
    /**
     * Volumetric constitutive component push-forward
     */
    double& VolumetricConstitutiveComponent(double & rCabcd,
					    const Matrix & rMatrixIC, 
					    const Matrix & rF,
					    const double & rdetF0, 
					    const double & rLameLambda, 
					    const double & rLameMu,
					    const double & rPressure,
					    const unsigned int& a, const unsigned int& b, 
					    const unsigned int& c, const unsigned int& d);
    
 
    /**
     * Calculates the isochoric stress vector
     * @param rMatrixIC can be the inverse of the RightCauchyGreen or the LeftCauchy Green tensor
     * @param rIdentityMatrix can be the IdentityMatrix
     * @param rdetF0 the determinant of the total deformation gradient
     * @param rTrace the trace of the RightCauchyGreen of the LeftCauchyGreen corresponding with the MatrixIC
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculateIsochoricStress( const Matrix & rMatrixIC,
					   const Matrix & rIdentityMatrix,
					   const double & rdetF0,
					   const double & rTrace,
					   const double & rLameLambda, 
					   const double & rLameMu, 
					   StressMeasure rStressMeasure,
					   Vector& rIsoStressVector);

    /**
     * Calculates the volumetric stress vector
     * @param rMatrixIC can be the inverse of the RightCauchyGreen or the LeftCauchy Green tensor
     * @param rdetF0 the determinant of the total deformation gradient
     * @param rLameLambda lame paramenter lambda
     * @param rLameMu lame paramenter mu
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rVolStressVector vector where the stress result is stored
     */
    virtual void CalculateVolumetricStress( const Matrix & rMatrixIC,
					    const double & rdetF0,
					    const double & rLameLambda, 
					    const double & rLameMu, 
					    const GeometryType& rDomainGeometry,
					    const Vector & rShapeFunctions,
					    Vector& rVolStressVector );



		
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);
    }



  }; // Class HyperElasticMixed3DLaw 
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_MIXED_3D_LAW_H_INCLUDED  defined 
