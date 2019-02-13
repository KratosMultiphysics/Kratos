// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//                   Vicente Mataix
//                   Fernando Rastellini
//  Collaborator:    Lucia Barbu
//

#if !defined(KRATOS_SERIAL_PARALLEL_RULE_OF_MIXTURES_H_INCLUDED)
#define KRATOS_SERIAL_PARALLEL_RULE_OF_MIXTURES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"


namespace Kratos
{
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
/**
 * @class GenericConstitutiveLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief This CL implements the serial-parallel rule of mixtures developed by F.Rastellini 
 * @details
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SerialParallelRuleOfMixturesLaw
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{


    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();
    /// Counted pointer of SerialParallelRuleOfMixturesLaw
    KRATOS_CLASS_POINTER_DEFINITION(SerialParallelRuleOfMixturesLaw);

    /// Definition of the base class
    typedef typename ElasticIsotropic3D BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    SerialParallelRuleOfMixturesLaw()
    {
    }

    /**
    * Constructor.
    */
    SerialParallelRuleOfMixturesLaw(double FiberVolParticipation, const Vector& rParallelDirections)
        : mFiberVolumetricParticipation(FiberVolParticipation), mParallelDirections(rParallelDirections)
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<SerialParallelRuleOfMixturesLaw>(*this);
    }

    // Copy constructor
    SerialParallelRuleOfMixturesLaw(SerialParallelRuleOfMixturesLaw const& rOther)
        : ConstitutiveLaw(rOther), mpMatrixConstitutiveLaw(rOther.mpMatrixConstitutiveLaw), mpFiberConstitutiveLaw(rOther.mpFiberConstitutiveLaw),
        mFiberVolumetricParticipation(rOther.mFiberVolumetricParticipation), mParallelDirections(rOther.mParallelDirections)
    {
    }

    /**
    * Destructor.
    */
    ~SerialParallelRuleOfMixturesLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    bool Has(const Variable<bool>& rThisVariable) override;
    bool Has(const Variable<double>& rThisVariable) override;
    bool Has(const Variable<Vector>& rThisVariable) override;
    bool Has(const Variable<Matrix>& rThisVariable) override;

    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues);


    void IntegrateStrainSerialParallelBehaviour(
        const Vector& rStrainVector,
        Vector& FiberStressVector,
        Vector& MatrixStressVector,
        const Properties& rMaterialProperties,
        ConstitutiveLaw::Parameters& rValues,
        Vector& rSerialStrainMatrix);

    void CalculateSerialParallelProjectionMatrices(
        Matrix& rFiberStressVector,
        Matrix& rMatrixStressVector);

    void InitializeMaterialResponsePK2(Parameters& rValues);  

    void CalculateStrainsOnEachComponent(
        const Vector& rStrainVector,
        const Properties& rMaterialProperties,
        const Matrix& rParallelProjector,
        const Matrix& rSerialProjector,
        const Vector& rSerialStrainMatrix,
        Vector& rStrainVectorMatrix,
        Vector& rStrainVectorFiber);

    void CalculateInitialApproximationSerialStrainMatrix(
        const Vector& rStrainVector,
        const Vector& rPreviousStrainVector,
        const Properties& rMaterialProperties,
        const Matrix& rParallelProjector,
        const Matrix& rSerialProjector,
        Matrix& rConstitutiveTensorMatrixSS,
        Matrix& rConstitutiveTensorFiberSS,
        Vector& rInitialApproximationSerialStrainMatrix);

    void IntegrateStressesOfFiberAndMatrix(
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rMatrixStrainVector,
        const Vector& rFiberStrainVector,
        Vector& rMatrixStressVector,
        Vector& rFiberStressVector);

    void CheckStressEquilibrium(
        const Vector& rStrainVector,
        const Matrix& rSerialProjector,
        const Vector& rMatrixStressVector,
        const Vector& rFiberStressVector,
        Vector& rStressResidual,
        bool& rIsConverged,
        const Matrix& rConstitutiveTensorMatrixSS,
        const Matrix& rConstitutiveTensorFiberSS);

    void CorrectSerialStrainMatrix(
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rResidualStresses,
        Vector& rSerialStrainMatrix,
        const Matrix& rSerialProjector);

    void CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues);
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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

    ConstitutiveLaw::Pointer GetMatrixConstitutiveLaw()
    {
        return mpMatrixConstitutiveLaw;
    }

    void SetMatrixConstitutiveLaw(ConstitutiveLaw::Pointer pMatrixConstitutiveLaw)
    {
        mpMatrixConstitutiveLaw = pMatrixConstitutiveLaw;
    }

    ConstitutiveLaw::Pointer GetFiberConstitutiveLaw()
    {
        return mpFiberConstitutiveLaw;
    }

    void SetFiberConstitutiveLaw(ConstitutiveLaw::Pointer pFiberConstitutiveLaw)
    {
        mpFiberConstitutiveLaw = pFiberConstitutiveLaw;
    }

    int GetNumberOfSerialComponents()
    {
        const int parallel_components = inner_prod(mParallelDirections, mParallelDirections);
        return this->GetStrainSize() - parallel_components;
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

    ConstitutiveLaw::Pointer mpMatrixConstitutiveLaw;
    ConstitutiveLaw::Pointer mpFiberConstitutiveLaw;
    double mFiberVolumetricParticipation;
    Vector mParallelDirections = ZeroVector(6);
    Vector mPreviousStrainVector = ZeroVector(6);
    Vector mPreviousSerialStrainMatrix = ZeroVector(GetNumberOfSerialComponents());

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the elastic tensor
     * @param rElasticityTensor The elastic tensor
     * @param rMaterialProperties The material properties
     */
    void CalculateElasticMatrix(
        Matrix& rElasticityTensor,
        const Properties& rMaterialProperties);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("MatrixConstitutiveLaw", mpMatrixConstitutiveLaw);
        rSerializer.save("FiberConstitutiveLaw", mpFiberConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("MatrixConstitutiveLaw", mpMatrixConstitutiveLaw);
        rSerializer.load("FiberConstitutiveLaw", mpFiberConstitutiveLaw);
    }

    ///@}

}; // Class SerialParallelRuleOfMixturesLaw

} // namespace Kratos

#endif
