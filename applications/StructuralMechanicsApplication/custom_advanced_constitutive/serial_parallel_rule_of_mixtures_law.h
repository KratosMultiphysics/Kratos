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
 * @class SerialParallelRuleOfMixturesLaw
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

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;


    /**
     * This method computes the stresses in the matrix and fiber according to the Serial-Parallel RoM
     * @param rStrainVector The total strain of the composite
     * @param FiberStressVector the Stress of the Fiber
     * @param MatrixStressVector the Stress of the Matrix
     * @param rMaterialProperties the Properties instance of the current element
     * @param rValues the needed parameters for the CL calculation
     * @param rSerialStrainMatrix the serial component of the matrix strain vector
     */
    void IntegrateStrainSerialParallelBehaviour(
        const Vector& rStrainVector,
        Vector& FiberStressVector,
        Vector& MatrixStressVector,
        const Properties& rMaterialProperties,
        ConstitutiveLaw::Parameters& rValues,
        Vector& rSerialStrainMatrix);

    /**
     * This method computes the projection tensors that divide the serial & paralle behaviours of the Strain/Stress
     * @param rParallelProjector The Parallel behaviour projector
     * @param rSerialProjector The Serial behaviour projector
     */
    void CalculateSerialParallelProjectionMatrices(
        Matrix& rParallelProjector,
        Matrix& rSerialProjector);

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(Parameters& rValues) override;  

    /**
     * This method computes the strain vector in the fiber and matrix according to the total 
     * strain and the serial strain of the matrix
     * @param rStrainVector The total strain of the composite
     * @param rParallelProjector The Parallel behaviour projector
     * @param rSerialProjector The Serial behaviour projector
     * @param rSerialStrainMatrix the serial component of the matrix strain vector
     * @param rStrainVectorMatrix the strain vector of the matrix
     * @param rStrainVectorFiber  the strain vector of the fiber
     */
    void CalculateStrainsOnEachComponent(
        const Vector& rStrainVector,
        const Matrix& rParallelProjector,
        const Matrix& rSerialProjector,
        const Vector& rSerialStrainMatrix,
        Vector& rStrainVectorMatrix,
        Vector& rStrainVectorFiber);

    /**
     * This method computes the initial aproximation of the Newton-Raphson procedure
     * regarding the serial strain of the matrix
     * @param rStrainVector The total strain of the composite
     * @param rPreviousStrainVector The total strain of the composite of the previous step
     * @param rMaterialProperties the Properties instance of the current element
     * @param rParallelProjector The Parallel behaviour projector
     * @param rSerialProjector The Serial behaviour projector
     * @param rConstitutiveTensorMatrixSS the serial-serial components of the constitutive tensor of the matrix
     * @param rConstitutiveTensorFiberSS  the serial-serial components of the constitutive tensor of the fiber
     * @param rInitialApproximationSerialStrainMatrix  initial aproximation of the serial strain of the matrix
     */
    void CalculateInitialApproximationSerialStrainMatrix(
        const Vector& rStrainVector,
        const Vector& rPreviousStrainVector,
        const Properties& rMaterialProperties,
        const Matrix& rParallelProjector,
        const Matrix& rSerialProjector,
        Matrix& rConstitutiveTensorMatrixSS,
        Matrix& rConstitutiveTensorFiberSS,
        Vector& rInitialApproximationSerialStrainMatrix);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param rValues the needed parameters for the CL calculation
     * @param rMatrixStrainVector the strain vector of the matrix
     * @param rMatrixStressVector  the stress vector of the matrix
     * @param rFiberStressVector  the stress vector of the fiber
     */
    void IntegrateStressesOfFiberAndMatrix(
        ConstitutiveLaw::Parameters& rValues,
        Vector rMatrixStrainVector,
        Vector rFiberStrainVector,
        Vector& rMatrixStressVector,
        Vector& rFiberStressVector);

    /**
     * This method checks wether the serial stresses are in equilibrium
     * @param rStrainVector The total strain of the composite
     * @param rSerialProjector The Serial behaviour projector
     * @param rMatrixStressVector  the stress vector of the matrix
     * @param rFiberStressVector  the stress vector of the fiber
     * @param rStressResidual  the stress residual between the serial stresses
     * @param rIsConverged  boolean indicator, true if equilibrium is achieved
     * @param rConstitutiveTensorMatrixSS the serial-serial components of the constitutive tensor of the matrix
     * @param rConstitutiveTensorFiberSS  the serial-serial components of the constitutive tensor of the fiber
     */
    void CheckStressEquilibrium(
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rStrainVector,
        const Matrix& rSerialProjector,
        const Vector& rMatrixStressVector,
        const Vector& rFiberStressVector,
        Vector& rStressResidual,
        bool& rIsConverged,
        const Matrix& rConstitutiveTensorMatrixSS,
        const Matrix& rConstitutiveTensorFiberSS);

    /**
     * This method updates the serial strain of the matrix in order to reach equilibrium
     * @param rValues the needed parameters for the CL calculation
     * @param rResidualStresses  the stress residual between the serial stresses
     * @param rSerialStrainMatrix the serial component of the matrix strain vector
     * @param rSerialProjector The Serial behaviour projector
     */
    void CorrectSerialStrainMatrix(
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rResidualStresses,
        Vector& rSerialStrainMatrix,
        const Matrix& rSerialProjector);

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues);

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }
    
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

    /**
     * @brief This method the constitutive law of the matrix material
     */
    ConstitutiveLaw::Pointer GetMatrixConstitutiveLaw()
    {
        return mpMatrixConstitutiveLaw;
    }

    /**
     * @brief This method sets the constitutive law of the matrix material
     */
    void SetMatrixConstitutiveLaw(ConstitutiveLaw::Pointer pMatrixConstitutiveLaw)
    {
        mpMatrixConstitutiveLaw = pMatrixConstitutiveLaw;
    }

    /**
     * @brief This method the constitutive law of the fiber material
     */
    ConstitutiveLaw::Pointer GetFiberConstitutiveLaw()
    {
        return mpFiberConstitutiveLaw;
    }

    /**
     * @brief This method sets the constitutive law of the fiber material
     */
    void SetFiberConstitutiveLaw(ConstitutiveLaw::Pointer pFiberConstitutiveLaw)
    {
        mpFiberConstitutiveLaw = pFiberConstitutiveLaw;
    }

    /**
     * @brief This method returns the number of directions
     * with serial behaviour (iso-stress behaviour)
     */
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
