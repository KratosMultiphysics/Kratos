// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Vicente Mataix
//
//
//

#pragma once

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
 * @brief This CL implements the serial-parallel rule of mixtures detailed in Cornejo et al. "Methodology for the analysis of post-tensioned structures using a constitutive serial-parallel rule of mixtures"
 * DOI: https://doi.org/10.1016/j.compstruct.2018.05.123
 * The SP-RoM is able to work in 2D and in 3D.
 * @details
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) SerialParallelRuleOfMixturesLaw
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// The node definition
    using NodeType = Node;

    /// The geometry definition
    using GeometryType = Geometry<NodeType>;

    using BaseType = ConstitutiveLaw;

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
        : BaseType(rOther), mpMatrixConstitutiveLaw(rOther.mpMatrixConstitutiveLaw), mpFiberConstitutiveLaw(rOther.mpFiberConstitutiveLaw),
        mFiberVolumetricParticipation(rOther.mFiberVolumetricParticipation), mParallelDirections(rOther.mParallelDirections) , 
        mPreviousStrainVector(rOther.mPreviousStrainVector) , mPreviousSerialStrainMatrix(rOther.mPreviousSerialStrainMatrix) , mIsPrestressed(rOther.mIsPrestressed) 
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
        KRATOS_DEBUG_ERROR_IF(mpMatrixConstitutiveLaw->WorkingSpaceDimension() != mpFiberConstitutiveLaw->WorkingSpaceDimension()) << "The WorkingSpaceDimension of the fiber and matrix mismatch..." << std::endl;
        return mpMatrixConstitutiveLaw->WorkingSpaceDimension();
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        KRATOS_DEBUG_ERROR_IF(mpMatrixConstitutiveLaw->GetStrainSize() != mpFiberConstitutiveLaw->GetStrainSize()) << "The GetStrainSize of the fiber and matrix mismatch..." << std::endl;
        return mpMatrixConstitutiveLaw->GetStrainSize();

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
     * @brief Returns the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

    /**
     * @brief Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(const Variable<int>& rThisVariable, int& rValue) override;

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
     * @brief Sets the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override;

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
     * @brief Calculates the value of a specified variable (bool)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& CalculateValue(
        Parameters& rParameterValues,
        const Variable<bool>& rThisVariable,
        bool& rValue) override;

    /**
     * @brief Calculates the value of a specified variable (integer)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    int& CalculateValue(
        Parameters& rParameterValues,
        const Variable<int>& rThisVariable,
        int& rValue) override;

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
        Vector& rSerialStrainMatrix,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy);

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
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */

    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */

    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * This method computes the Green-Lagrange strain
     * @see Parameters
     */
    void CalculateGreenLagrangeStrain(Parameters &rValues);

    /**
     * This method computes the Almansi strain
     * @see Parameters
     */
    void CalculateAlmansiStrain(Parameters &rValues);

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
        Vector& rStrainVectorFiber,
        ConstitutiveLaw::Parameters& rValues,
        const int Iteration = 1);

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
        Vector& rInitialApproximationSerialStrainMatrix,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure& rStressMeasure);

    /**
     * This method computes the stresses of the matrix/fiber according to its own CL
     * @param rValues the needed parameters for the CL calculation
     * @param rMatrixStrainVector the strain vector of the matrix
     * @param rMatrixStressVector  the stress vector of the matrix
     * @param rFiberStressVector  the stress vector of the fiber
     */
    void IntegrateStressesOfFiberAndMatrix(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rMatrixStrainVector,
        Vector& rFiberStrainVector,
        Vector& rMatrixStressVector,
        Vector& rFiberStressVector,
        const ConstitutiveLaw::StressMeasure& rStressMeasure);

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
        const Matrix& rSerialProjector,
        const ConstitutiveLaw::StressMeasure& rStressMeasure);

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy);

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        if (mpMatrixConstitutiveLaw->RequiresInitializeMaterialResponse() || mpFiberConstitutiveLaw->RequiresInitializeMaterialResponse()) {
            return true;
        }
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        if (mpMatrixConstitutiveLaw->RequiresFinalizeMaterialResponse() || mpFiberConstitutiveLaw->RequiresFinalizeMaterialResponse()) {
            return true;
        }
        return false;
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

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
        return GetStrainSize() - parallel_components;
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
    Vector mParallelDirections   = ZeroVector(GetStrainSize());
    Vector mPreviousStrainVector = ZeroVector(GetStrainSize());
    Vector mPreviousSerialStrainMatrix = ZeroVector(GetNumberOfSerialComponents());
    bool mIsPrestressed = false;

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
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
        rSerializer.save("MatrixConstitutiveLaw", mpMatrixConstitutiveLaw);
        rSerializer.save("FiberConstitutiveLaw", mpFiberConstitutiveLaw);
        rSerializer.save("FiberVolumetricParticipation", mFiberVolumetricParticipation);
        rSerializer.save("ParallelDirections", mParallelDirections);
        rSerializer.save("PreviousStrainVector", mPreviousStrainVector);
        rSerializer.save("PreviousSerialStrainMatrix", mPreviousSerialStrainMatrix);
        rSerializer.save("IsPrestressed", mIsPrestressed);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
        rSerializer.load("MatrixConstitutiveLaw", mpMatrixConstitutiveLaw);
        rSerializer.load("FiberConstitutiveLaw", mpFiberConstitutiveLaw);
        rSerializer.load("FiberVolumetricParticipation", mFiberVolumetricParticipation);
        rSerializer.load("ParallelDirections", mParallelDirections);
        rSerializer.load("PreviousStrainVector", mPreviousStrainVector);
        rSerializer.load("PreviousSerialStrainMatrix", mPreviousSerialStrainMatrix);
        rSerializer.load("IsPrestressed", mIsPrestressed);
    }

    ///@}

}; // Class SerialParallelRuleOfMixturesLaw

} // namespace Kratos
