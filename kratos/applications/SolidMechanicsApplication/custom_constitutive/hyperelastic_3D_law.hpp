//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 3D Neohookean Model
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO 
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) HyperElastic3DLaw : public ConstitutiveLaw
{
protected:


    struct MaterialResponseVariables
    {
        //general material properties
        double LameMu;
        double LameLambda;
      
        //general thermal properties
        double ThermalExpansionCoefficient;
        double ReferenceTemperature;

        //kinematic properties
        double J_pow13;
        double DeterminantF;
        double traceCG;               //LeftCauchyGreen or RightCauchyGreen
        Matrix CauchyGreenMatrix;     //LeftCauchyGreen or InverseRightCauchyGreen
        Matrix DeformationGradientF;  //Deformation Gradient Tensor in 3D
        Matrix Identity;

        //element properties
        const Vector*        mpShapeFunctionsValues;
        const GeometryType*  mpElementGeometry;

    public: 
      void SetShapeFunctionsValues (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};
      void SetElementGeometry      (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};
      const Vector& GetShapeFunctionsValues      () const {return *mpShapeFunctionsValues;};
      const GeometryType& GetElementGeometry     () const {return *mpElementGeometry;};
      

    };


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /**
     * Counted pointer of HyperElastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElastic3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElastic3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    HyperElastic3DLaw (const HyperElastic3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElastic3DLaw& operator=(const HyperElastic3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HyperElastic3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
        return 6;
    };


    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );


    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues );


    void InitializeSolutionStep( const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                 const Vector& rShapeFunctionsValues ,
                                 const ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep( const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry, //this is just to give the array of nodes
                               const Vector& rShapeFunctionsValues ,
                               const ProcessInfo& rCurrentProcessInfo);

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK1 (Parameters & rValues);

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters & rValues);


    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponsePK1 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponsePK2 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponseKirchhoff (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponseCauchy (Parameters & rValues);


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

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

    Matrix mInverseDeformationGradientF0;

    double mDeterminantF0;

    double mStrainEnergy;

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
    virtual void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
            Vector& rStrainVector );


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector );


    /**
     * Calculates the stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rStressVector vector where the stress result is stored
     */
    void CalculateStress( const MaterialResponseVariables& rElasticVariables,
                          StressMeasure rStressMeasure,
                          Vector& rStressVector);

   /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
                                           StressMeasure rStressMeasure,
					   Vector& rIsoStressVector);

    /**
     * Calculates the volumetric stress vector
     * @param rElasticResponseVariables the material variables
     * @param rVolStressVector vector where the stress result is stored
     */
    virtual void CalculateVolumetricStress( const MaterialResponseVariables & rElasticVariables,
                                            Vector& rVolStressVector );

    /**
     * Calculates the constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
					      Matrix& rConstitutiveMatrix);


    /**
     * Constitutive component
     */

    double& ConstitutiveComponent( double & rCabcd,
                                   const MaterialResponseVariables& rElasticVariables,
                                   const unsigned int& a, const unsigned int& b,
                                   const unsigned int& c, const unsigned int& d);

							      
    /**
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						       const Matrix & rIsoStressMatrix,
						       Matrix& rConstitutiveMatrix);


    /**
     * Constitutive isochoric component
     */
    double& IsochoricConstitutiveComponent( double & rCabcd,
                                            const MaterialResponseVariables& rElasticVariables,
                                            const Matrix & rIsoStressMatrix,
                                            const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);

    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
							Matrix& rConstitutiveMatrix);


    /**
     * Constitutive volumetric component
     */

    double& VolumetricConstitutiveComponent( double & rCabcd,
					     const MaterialResponseVariables& rElasticVariables,
					     const Vector& rFactors,
					     const unsigned int& a, const unsigned int& b,
					     const unsigned int& c, const unsigned int& d);


    /**
     * Calculates HyperElasticLaw Factor for the Neo-Hookean model
     * @param rElasticResponseVariables the material variables
     * @param rFactor the calculated factor to be returned
     */
    virtual double& CalculateVolumetricFactor (const MaterialResponseVariables & rElasticVariables,
					       double & rFactor);


    /**
     * Calculates the Pressure of the domain (element)
     * @param rElasticResponseVariables the material variables
     * @param rPressure the calculated pressure to be returned
     */
    virtual double& CalculateVolumetricPressure (const MaterialResponseVariables & rElasticVariables,
						 double & rPressure);


    /**
     * Calculates the Volumetric part factors
     * @param rElasticResponseVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    virtual Vector&  CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
							 Vector & rFactors);


    /**
     * Calculates the Temperature of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rTemperature the calculated temperature to be returned
     */
    virtual double& CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
						double & rTemperature);

    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     * if the matrix passed is 3D is does nothing
     * if the matrix passed is bigger or smaller throws an error
     * @param rMatrix : usually the DeformationGradientF
     */
    Matrix& Transform2DTo3D (Matrix& rMatrix);



    /**
      * Updates the material response:
      * Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void UpdateInternalVariables (Parameters & rValues);

    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    virtual bool CheckParameters(Parameters& rValues);


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


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
	rSerializer.save("mInverseDeformationGradientF0",mInverseDeformationGradientF0);
	rSerializer.save("mDeterminantF0",mDeterminantF0);
	rSerializer.save("mStrainEnergy",mStrainEnergy);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
	rSerializer.load("mInverseDeformationGradientF0",mInverseDeformationGradientF0);
	rSerializer.load("mDeterminantF0",mDeterminantF0);
	rSerializer.load("mStrainEnergy",mStrainEnergy);
    }


    ///@}

}; // Class HyperElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_3D_LAW_H_INCLUDED  defined 
