//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/hyperelastic_3D_law.hpp


#if !defined (KRATOS_HYPERELASTIC_3D_LAW_H_INCLUDED)
#define       KRATOS_HYPERELASTIC_3D_LAW_H_INCLUDED

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

class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) HyperElastic3DLaw
    : public ConstitutiveLaw
{
protected:


    struct MaterialResponseVariables
    {
        // general material properties
        double LameMu;
        double LameLambda;

        // general thermal properties
        double ThermalExpansionCoefficient;
        double ReferenceTemperature;

        // kinematic properties
        double DeterminantF;
        double traceCG;               // LeftCauchyGreen or RightCauchyGreen
        Matrix CauchyGreenMatrix;     // LeftCauchyGreen or InverseRightCauchyGreen
        Matrix DeformationGradientF;  // Deformation Gradient Tensor in 3D
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
    ConstitutiveLaw::Pointer Clone() const override;

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
    ~HyperElastic3DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };


    bool Has( const Variable<double>& rThisVariable ) override;
    bool Has( const Variable<Vector>& rThisVariable ) override;
    bool Has( const Variable<Matrix>& rThisVariable ) override;

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue ) override;


    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;
    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues ) override;


    void InitializeSolutionStep( const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                 const Vector& rShapeFunctionsValues ,
                                 const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep( const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry, //this is just to give the array of nodes
                               const Vector& rShapeFunctionsValues ,
                               const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;


    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //String Info() const override;
    /**
     * Print information about this object.
     */
    //void PrintInfo(std::ostream& rOStream) const override;
    /**
     * Print object's data.
     */
    //void PrintData(std::ostream& rOStream) const override;

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
     * Calculates HyperElasticLaw Factor for the Neo-Hookean model
     * @param rElasticResponseVariables the material variables
     * @param rFactor the calculated factor to be returned
     */
    virtual double& CalculateVolumetricFactor (const MaterialResponseVariables & rElasticVariables,
					       double & rFactor);


    /**
     * Calculates the Volumetric part factors
     * @param rElasticResponseVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    virtual Vector&  CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
							 Vector & rFactors);


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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
        rSerializer.save("mInverseDeformationGradientF0",mInverseDeformationGradientF0);
        rSerializer.save("mDeterminantF0",mDeterminantF0);
        rSerializer.save("mStrainEnergy",mStrainEnergy);
    }

    void load(Serializer& rSerializer) override
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
