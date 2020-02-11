// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//


#if !defined(KRATOS_GENERIC_TOTAL_LAGRANGIAN_FEMDEM_H_INCLUDED )
#define  KRATOS_GENERIC_TOTAL_LAGRANGIAN_FEMDEM_H_INCLUDED


// System includes


// External include

// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

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
 * @class GenericTotalLagrangianFemDemElement
 * @ingroup StructuralMechanicsApplication
 * @brief Total Lagrangian element for 2D and 3D geometries.
 * @details Implements a total Lagrangian definition for structural analysis. This works for arbitrary geometries in 2D and 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericTotalLagrangianFemDemElement
    : public BaseSolidElement
{
public:
    ///@name Type Definitions
    ///@{

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef BaseSolidElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// We define the dimension
    static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

    /// We define the number of edges
    static constexpr SizeType NumberOfEdges = (TDim == 3) ? 6 : 3;

    /// Counted pointer of GenericTotalLagrangianFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericTotalLagrangianFemDemElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenericTotalLagrangianFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GenericTotalLagrangianFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    GenericTotalLagrangianFemDemElement(GenericTotalLagrangianFemDemElement const& rOther)
        : BaseSolidElement(rOther)
    {};

    /// Destructor.
    ~GenericTotalLagrangianFemDemElement() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    // int Check(const ProcessInfo& rCurrentProcessInfo) override;

    //std::string Info() const;

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Total Lagrangian FEMDEM Solid Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Total Lagrangian FEMDEM Solid Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    Vector mThresholds;                 // Stress mThreshold on edge
    Vector mDamages;                    // Converged Damage on each edge

    double mThreshold = 0.0;            // Converged Threshold
    double mDamage = 0.0;               // Converged Damage
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    GenericTotalLagrangianFemDemElement() : BaseSolidElement()
    {
    }

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) override;

    // ************** Methods to compute the tangent constitutive tensor via numerical derivation ************** 
    /**
     * this computes the Tangent tensor via numerical derivation (perturbations)
     */
    void CalculateTangentTensor(Matrix& rTangentTensor, const Vector& rStrainVectorGP, const Vector& rStressVectorGP, const Matrix& rDeformationGradientGP, const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);
    void CalculateTangentTensorSecondOrder(Matrix& rTangentTensor, const Vector& rStrainVectorGP, const Vector& rStressVectorGP, const Matrix& rDeformationGradientGP, const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);

    /**
     * this computes the perturbation to the strain
     */
    void CalculatePerturbation(const Vector& rStrainVectorGP, double& rPerturbation, const int Component);

    /**
     * this perturbates the strain vector
     */
    void PerturbateStrainVector(Vector& rPerturbedStrainVector, const Vector& rStrainVectorGP, const double Perturbation, const int Component);

    /**
     * this integrated the perturbed strain
     */
    void IntegratePerturbedStrain(Vector& rPerturbedStressVector, const Vector& rPerturbedStrainVector, const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);

    /**
     * this assings the components to the tangent tensor
     */
    void AssignComponentsToTangentTensor(Matrix& rTangentTensor, const Vector& rDeltaStress, const double Perturbation, const int Component);

    /**
     * this perturbates F
     */
    void PerturbateDeformationGradient(Matrix& rPerturbedDeformationGradient, const Matrix& rDeformationGradientGP, const double Perturbation, const int ComponentI, const int ComponentJ);

    /**
     * this gets the voigt index for a set of components
     */
    int CalculateVoigtIndex(const SizeType VoigtSize, const int ComponentI, const int ComponentJ);

    /**
     * this computes the damage of the FE
     */
    double CalculateElementalDamage(const Vector& rEdgeDamages);
    double CalculateElementalDamage3D(const Vector& rEdgeDamages);
    double CalculateElementalDamage2D(const Vector& rEdgeDamages);

    /**
     * this integrates the constitutive law
     */
    void IntegrateStressDamageMechanics(double& rThreshold,double& rDamage, const Vector& rStrainVector,
        const Vector& rStressVector, const int Edge, const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues, bool& rIsDamaging);

    /**
     * this computes the elements that share an edge -> fills the mEdgeNeighboursContainer
     */
    void ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo);

    /**
     * this computes the elements that share an edge -> fills the mEdgeNeighboursContainer
     */
    void AuxComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo);

    /**
     * this returns the elements that share an edge -> gets the mEdgeNeighboursContainer
     */
    std::vector<Element*> GetEdgeNeighbourElements(const int edge) {return mEdgeNeighboursContainer[edge];}

    /**
     * this storages the mEdgeNeighboursContainer
     */
    void SaveEdgeNeighboursContainer(const std::vector<std::vector<Element*>>& rtoSave) {mEdgeNeighboursContainer = rtoSave;}

    /**
     * this sets the numbering for several purposes
     * at the edges 
     */
    void SetNodeIndexes(Matrix& rMatrix)
    {
        rMatrix.resize(6, 2);
        rMatrix(0, 0) = 0; rMatrix(0, 1) = 1; rMatrix(1, 0) = 0;
        rMatrix(1, 1) = 2; rMatrix(2, 0) = 0; rMatrix(2, 1) = 3;
        rMatrix(3, 0) = 1; rMatrix(3, 1) = 2; rMatrix(4, 0) = 1;
        rMatrix(4, 1) = 3; rMatrix(5, 0) = 2; rMatrix(5, 1) = 3;
    }

    /**
     * this imposes the damage/threshold to be equal
     * at the edges 
     */
    void InitializeInternalVariablesAfterMapping();

    /**
     * this computes the average vector on the edge for a certain variable
     */
    void CalculateAverageVariableOnEdge(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);
    void CalculateAverageVariableOnEdge2D(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);
    void CalculateAverageVariableOnEdge3D(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);

    /**
     * this evaluates the constitutive law
     */
    void CalculateEquivalentStress(const array_1d<double, VoigtSize>& rPredictiveStressVector, const Vector& rStrainVector, double& rEquivalentStress, ConstitutiveLaw::Parameters& rValues);

    /**
     * this gets the initial threshold of the yield surface
     */
    void GetInitialUniaxialThreshold(ConstitutiveLaw::Parameters& rValues, double& rThreshold);

    /**
     * this computes the damage parameter "A"
     */
    void CalculateDamageParameter(ConstitutiveLaw::Parameters& rValues, double& rAParameter, const double CharacteristicLength);

    /**
     * this computes the CharacteristicLength of the element
     */
    double CalculateCharacteristicLength(GenericTotalLagrangianFemDemElement *pCurrentElement);

    /**
     * this computes the damage according to a exp softening
     */
    void CalculateExponentialDamage(double& rDamage, const double DamageParameter, const double UniaxialStress, const double InitialThrehsold);

    /**
     * this computes the Green-Lagrange Strain vector from F
     */
    void CalculateGreenLagrangeStrainVector(Vector& rStrainVector, const Matrix& rF);

    std::size_t GetStrainSize() const;

    void GetValueOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rValues,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rOutput,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;
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

    /**
     * this is called in the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called at the end of each solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    /**
     * @brief This method computes the deformation matrix B
     * @param rB The deformation matrix
     * @param rF The deformation gradient
     * @param rDN_DX The gradient derivative of the shape function
     */
    void CalculateB(Matrix& rB, Matrix const& rF, const Matrix& rDN_DX);

    void Calculate2DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);

    void Calculate3DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX);

    void CalculateStress(Vector& rStrain,
                         std::size_t IntegrationPoint,
                         Vector& rStress,
                         ProcessInfo const& rCurrentProcessInfo);

    void CalculateStress(Matrix const& rF,
                         std::size_t IntegrationPoint,
                         Vector& rStress,
                         ProcessInfo const& rCurrentProcessInfo);

    void CalculateStrain(Matrix const& rF,
                         std::size_t IntegrationPoint,
                         Vector& rStrain,
                         ProcessInfo const& rCurrentProcessInfo);

    void CalculateShapeSensitivity(ShapeParameter Deriv,
                                   Matrix& rDN_DX0,
                                   Matrix& rDN_DX0_Deriv,
                                   Matrix& rF_Deriv,
                                   double& rDetJ0_Deriv,
                                   std::size_t IntegrationPointIndex);

    void CalculateBSensitivity(Matrix const& rDN_DX,
                               Matrix const& rF,
                               Matrix const& rDN_DX_Deriv,
                               Matrix const& rF_Deriv,
                               Matrix& rB_Deriv);

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    // Vector to storage the neigh elements sharing a certain edge
    std::vector<std::vector<Element*>> mEdgeNeighboursContainer;

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class GenericTotalLagrangianFemDemElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_H_INCLUDED  defined
