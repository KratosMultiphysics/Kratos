//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


#if !defined (KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED)
#define KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"
#include "includes/ublas_interface.h"

namespace Kratos
{

/**
 * Defines a hencky-plastic isotropic constitutive law in 3D
 * The functionality is limited to large and small displacements
 */


class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) HenckyElasticPlastic3DLaw : public HyperElastic3DLaw
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
        Matrix  EigenValues;
        Matrix  EigenVectors;
        Matrix  Isochoric;
        Matrix  Volumetric;
        Matrix  Plastic;
    };

    struct PlasticMaterialResponseVariables
    {
        Matrix  TrialLeftStretchTensor;
        Matrix  InverseTrialLeftStretchTensor;
    };


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef MPMFlowRule::Pointer                MPMFlowRulePointer;
    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of HyperElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HenckyElasticPlastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyElasticPlastic3DLaw();


    HenckyElasticPlastic3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    HenckyElasticPlastic3DLaw (const HenckyElasticPlastic3DLaw& rOther)  ;


    /**
     * Assignment operator.
     */

    //HyperElasticPlastic3DLaw& operator=(const HyperElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~HenckyElasticPlastic3DLaw() override;

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


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    bool Has( const Variable<double>& rThisVariable ) override;
    bool Has( const Variable<Vector>& rThisVariable ) override;
    bool Has( const Variable<Matrix>& rThisVariable ) override;

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
    void InitializeMaterial( const Properties& rProps,
                             const GeometryType& rGeom,
                             const Vector& rShapeFunctionsValues ) override;


    void InitializeSolutionStep( const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                 const Vector& rShapeFunctionsValues,
                                 const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep( const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry, //this is just to give the array of nodes
                               const Vector& rShapeFunctionsValues,
                               const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    //virtual void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;



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
    unsigned int mPlasticRegion;
    Matrix mPlasticDeformationGradient;
    Matrix mElasticLeftCauchyGreen;

    MPMFlowRulePointer mpMPMFlowRule;

    YieldCriterionPointer mpYieldCriterion;

    HardeningLawPointer   mpHardeningLaw;

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS
        in a pullback fashion
    */
    void GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables);

    virtual Matrix SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix);

    virtual Vector SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix );

    virtual Matrix SetMatrixToAppropriateDimension(Matrix& rConstitutiveMatrix);

    virtual void CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables& rElasticVariables);

    virtual void CorrectKinematics(const PlasticMaterialResponseVariables& rPlasticVariables, Parameters & rValues, const MPMFlowRule::RadialReturnVariables rReturnMappingVariables, double& rDeterminantF, Matrix& rDeformationGradientF );

    virtual void CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen,const double& rAlpha, Matrix& rElastoPlasticMatrix, const MaterialResponseVariables& rElasticVariables);

    //virtual void CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
    //MPMFlowRule::RadialReturnVariables & rReturnMappingVariables,
    //Matrix& rConstitutiveMatrix);
    double& TensorComponent(double & rCabcd,
                            const Matrix& rMA, const Matrix& rMB,
                            const unsigned int& a, const unsigned int& b,
                            const unsigned int& c, const unsigned int& d);

    virtual void MyTensorProduct(const Matrix& rMA, const Matrix& rMB,
                                 Matrix& rEigenbasesProductMatrix);

    double& TensorComponent2(double & rCabcd,
                             const Matrix& rMA, const Matrix& rMB,
                             const unsigned int& a, const unsigned int& b,
                             const unsigned int& c, const unsigned int& d);

    virtual void MyTensorProduct2(const Matrix& rMA, const Matrix& rMB,
                                  Matrix& rEigenbasesProductMatrix);

    double& TensorComponent3(double & rCabcd,
                             const Matrix& rMA,
                             const unsigned int& a, const unsigned int& b,
                             const unsigned int& c, const unsigned int& d);

    virtual void MyTensorProduct3(const Matrix& rMA,
                                  Matrix& rEigenbasesProductMatrix);

    double& TensorComponent4(double & rCabcd,
                             const Matrix& rMA,
                             const unsigned int& a, const unsigned int& b,
                             const unsigned int& c, const unsigned int& d);

    virtual void MyTensorProduct4(const Matrix& rMA,
                                  Matrix& rEigenbasesProductMatrix);


    virtual Matrix CalculateEigenbases(const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rEigenbasesMatrix);

    //virtual void CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables,
    //const Matrix& rTrialElasticLeftCauchyGreen, const Matrix& rStressMatrix,
    //Matrix& rElastoPlasticTangentMatrix, Matrix& rConsistentMatrix );

    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS
        in the actual configuration
    */

    //double& PlasticConstitutiveComponent( double & rCabcd,
    //const MaterialResponseVariables& rElasticVariables,
    //const Matrix & rIsoStressMatrix,
    //const MPMFlowRule::PlasticFactors & rScalingFactors,
    //const unsigned int& a, const unsigned int& b,
    //const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressMatrix matrix where the stress result is stored
     * @param rIsoStressVector vector where the stress result is stored
     */
    //virtual void CalculatePlasticIsochoricStress( MaterialResponseVariables & rElasticVariables,
    //MPMFlowRule::RadialReturnVariables & rReturnMappingVariables,
    //StressMeasure rStressMeasure,
    //Matrix& rIsoStressMatrix,
    //Vector& rIsoStressVector);




    Vector& GetStressVectorFromMatrix(const Matrix& rStressMatrix,
                                      Vector& rMainStress,
                                      const Matrix& rEigenVectors);

    virtual void CalculateHenckyMainStrain(const Matrix& rCauchyGreeMatrix,
                                           MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
                                           Vector& rMainStrain);

    virtual void CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables,Parameters & rValues, const MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
            Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix);

    virtual void CalculateLeftStretchTensor(Matrix& rLeftStretchTensor, const Matrix& rCauchyGreenMatrix);

    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    bool CheckParameters(Parameters& rValues) override;


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )

        rSerializer.save("mElasticLeftCauchyGreen",mElasticLeftCauchyGreen);
        rSerializer.save("mpMPMFlowRule",mpMPMFlowRule);
        rSerializer.save("mpYieldCriterion",mpYieldCriterion);
        rSerializer.save("mpHardeningLaw",mpHardeningLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )

        rSerializer.load("mElasticLeftCauchyGreen",mElasticLeftCauchyGreen);
        rSerializer.load("mpMPMFlowRule",mpMPMFlowRule);
        rSerializer.load("mpYieldCriterion",mpYieldCriterion);
        rSerializer.load("mpHardeningLaw",mpHardeningLaw);
    }




}; // Class HenckyElasticPlastic3DLaw

} //namespace Kratos

#endif  //KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED

