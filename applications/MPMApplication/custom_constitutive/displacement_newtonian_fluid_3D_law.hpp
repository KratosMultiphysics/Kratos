//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Contri Alessandro
//
//  References:      This class is based on the file applications/MPMApplication/custom_constitutive/hyperelastic_3D_law.hpp


#if !defined (KRATOS_DISPLACEMENT_NEWTONIAN_3D_LAW_H_INCLUDED)
#define       KRATOS_DISPLACEMENT_NEWTONIAN_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"


namespace Kratos
{

/**
 * Defines a displacement-based newtonian fluid constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC VISCOSITY
 * 2) PRESSURE COEFFICIENT (Cole equation: p=BULK_MODULUS*((rho/rho_0)^(PRESSURE_COEFFICIENT)-1))
 * 1) BULK MODULUS
 * As there are no further parameters the functionality is limited
 * to large displacements.
 */

class KRATOS_API(MPM_APPLICATION) DispNewtonianFluid3DLaw
    : public ConstitutiveLaw
{
protected:

    struct MaterialResponseVariables
    {

	// general material properties
        double Mu; // shear modulus
        double BulkModulus;

	// kinematic properties
        double DeterminantF;
        double DeltaTime;
        Matrix DeformationGradientF;  // Deformation Gradient Tensor in 3D
        Matrix Identity;

        Matrix DeformationRate; // symmetric part of velocity gradient
	Matrix CauchyGreenMatrix;     // LeftCauchyGreen
        Matrix StressMatrix;

        const Vector* mpShapeFunctionsValues;
        const Matrix* mpShapeFunctionsDerivatives;
        const GeometryType* mpElementGeometry;

    public:
        void SetShapeFunctionsValues(const Vector& rShapeFunctionsValues)
        {
            mpShapeFunctionsValues = &rShapeFunctionsValues;
        };

        void SetElementGeometry(const GeometryType& rElementGeometry)
        {
            mpElementGeometry = &rElementGeometry;
        };

        const Vector& GetShapeFunctionsValues() const
        {
            return *mpShapeFunctionsValues;
        };

        const GeometryType& GetElementGeometry() const
        {
            return *mpElementGeometry;
        };
    };

public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of DispNewtonianFluid3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( DispNewtonianFluid3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    DispNewtonianFluid3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    DispNewtonianFluid3DLaw (const DispNewtonianFluid3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //DispNewtonianFluid3DLaw& operator=(const DispNewtonianFluid3DLaw& rOther);


    /**
     * Destructor.
     */
    ~DispNewtonianFluid3DLaw() override;

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
    SizeType GetStrainSize() const override
    {
        return 6;
    };

    bool Has(const Variable<double>& rThisVariable) override;
    bool Has(const Variable<Vector>& rThisVariable) override;
    bool Has(const Variable<Matrix>& rThisVariable) override;


    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Matrix>& rVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;


    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

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
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

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

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

   /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector );

    /**
     * Calculates the deformation rate
     * @param rViscousVariables
     */
    virtual void CalculateDeformationRate(MaterialResponseVariables& rViscousVariables);

    /**
     * Calculates the stress vector
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rStressVector vector where the stress result is stored
     */
    void CalculateStress(MaterialResponseVariables& rViscousVariables, StressMeasure rStressMeasure, Vector& rStressVector);

    /**
     * Calculates the volumetric stress vector
     * @param rViscousVariables the material variables
     * @param rVolStressVector vector where the stress result is stored
     */
    virtual void CalculateVolumetricStress(const MaterialResponseVariables& rViscousVariables, Vector& rVolStressVector);

    /**
     * Calculates the isochoric stress vector
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculateIsochoricStress(const MaterialResponseVariables& rViscousVariables, StressMeasure rStressMeasure, Vector& rIsoStressVector);


    /**
     * Calculates the constitutive matrix
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateConstitutiveMatrix ( const MaterialResponseVariables& rViscousVariables,
                                       Matrix& rConstitutiveMatrix);

    /**
     * Constitutive component
     */
    double& ConstitutiveComponent(double& rCabcd, const MaterialResponseVariables& rViscousVariables,
        const unsigned int& a, const unsigned int& b, const unsigned int& c, const unsigned int& d);

    /**
     * Calculates the volumetric constitutive matrix
     * @param rViscousVariables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateVolumetricConstitutiveMatrix(const MaterialResponseVariables& rViscousVariables, Matrix& rConstitutiveMatrix);

    /**
     * Constitutive volumetric component
     */
    double& VolumetricConstitutiveComponent(double& rCabcd, const MaterialResponseVariables& rViscousVariables, const Vector& rFactors, const unsigned int& a,
        const unsigned int& b, const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the isochoric constitutive matrix
     * @param rViscousVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateIsochoricConstitutiveMatrix(const MaterialResponseVariables& rViscousVariables, const Matrix& rIsoStressMatrix, Matrix& rConstitutiveMatrix);

    /**
     * Constitutive isochoric component
     */
    double& IsochoricConstitutiveComponent(double& rCabcd, const MaterialResponseVariables& rViscousVariables, const Matrix& rIsoStressMatrix, const unsigned int& a,
        const unsigned int& b, const unsigned int& c, const unsigned int& d);

    /**
     * Calculates the Pressure of the domain (element)
     * @param rViscousVariables the material variables
     * @param rPressure the calculated pressure to be returned
     */
    virtual double& CalculateVolumetricPressure(const MaterialResponseVariables& rViscousVariables,
        double& rPressure);

    /**
     * Calculates the Volumetric part factors
     * @param rViscousVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    virtual Vector& CalculateVolumetricPressureFactors(const MaterialResponseVariables& rViscousVariables, Vector& rFactors);

    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     * if the matrix passed is 3D is does nothing
     * if the matrix passed is bigger or smaller throws an error
     * @param rMatrix : usually the DeformationGradientF
     */
    Matrix& Transform2DTo3D(Matrix& rMatrix);

    /**
      * Updates the material response:
      * Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void UpdateInternalVariables(Parameters& rValues);

    /**
     * This function is designed to compute the deviatoric part of
     * a 3x3 matrix
     * @param rMatrix : input matrix to calculate the deviatoric part 
     * @param rDevMatrix : deviatoric part of rMatrix 
     */
    void CalculateDeviatoricPart(const Matrix& rMatrix, Matrix& rDevMatrix);

    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    virtual bool CheckParameters(Parameters& rValues);

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
            rSerializer.save("mInverseDeformationGradientF0", mInverseDeformationGradientF0);
        rSerializer.save("mDeterminantF0", mDeterminantF0);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
            rSerializer.load("mInverseDeformationGradientF0", mInverseDeformationGradientF0);
        rSerializer.load("mDeterminantF0", mDeterminantF0);
    }



}; // Class DispNewtonianFluid3DLaw
}  // namespace Kratos.
#endif // KRATOS_DISPLACEMENT_NEWTONIAN_3D_LAW_H_INCLUDED  defined
