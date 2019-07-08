//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED


// System includes


// External include

// Project includes

#include "custom_elements/solid_elements/small_displacement_element.hpp"
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
 * @class GenericSmallStrainFemDemElement
 * @ingroup FemToDemApplication
 * @brief Small Displacement element for the 2D and 3D cases
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericSmallStrainFemDemElement 
    : public SmallDisplacementElement // Derived Element from SolidMechanics
{
public:
    ///@name Type Definitions
    ///@{

    ///definition of element type
    typedef Element ElementType;

    ///base type: an GeometricalObject that automatically has a unique number
    typedef GeometricalObject BaseType;

    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef GeometryData GeometryDataType;

    /// We define the dimension
    static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

    /// We define the number of edges
    static constexpr SizeType NumberOfEdges = (TDim == 3) ? 6 : 3;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Counted pointer of GenericSmallStrainFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericSmallStrainFemDemElement);

    ///@}
    ///@name Life Cycle
    ///@{

	/// Default constructors
	GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
	GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	GenericSmallStrainFemDemElement(GenericSmallStrainFemDemElement const &rOther);

	/// Destructor.
	virtual ~GenericSmallStrainFemDemElement();

	/// Assignment operator.
	GenericSmallStrainFemDemElement &operator=(GenericSmallStrainFemDemElement const &rOther);
	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const override;
	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

	GenericSmallStrainFemDemElement()
	{
	}

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
     * this is called for non-linear analysis at the end of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this computes the elements that share an edge -> fills the mEdgeNeighboursContainer
     */
    void ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo);
    void AuxComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo);
    std::vector<Element*> GetEdgeNeighbourElements(const int edge) {return mEdgeNeighboursContainer[edge];}

    /**
     * this storages the mEdgeNeighboursContainer
     */
    void SaveEdgeNeighboursContainer(const std::vector<std::vector<Element*>>& rtoSave) {mEdgeNeighboursContainer = rtoSave;}

	void SetNodeIndexes(Matrix& rMatrix) // Defines the numbering of the edges with the corresponding nodes
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
     * this saves the converged values with the later non-conv values
     */
    void UpdateDataBase();

    /**
     * this computes the damage of the FE
     */
    double CalculateElementalDamage(const Vector& rEdgeDamages);
    double CalculateElementalDamage3D(const Vector& rEdgeDamages);
    double CalculateElementalDamage2D(const Vector& rEdgeDamages);

    /**
     * this computes the average vector on the edge for a certain variable
     */
	void CalculateAverageVariableOnEdge(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);
    void CalculateAverageVariableOnEdge2D(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);
    void CalculateAverageVariableOnEdge3D(const Element* pCurrentElement, const Variable<Vector> ThisVariable, Vector& rAverageStress, const int edge);

    /**
     * this integrates the constitutive law
     */
    void IntegrateStressDamageMechanics(double& rThreshold,double& rDamage, const Vector& rStrainVector,
        const Vector& rStressVector, const int Edge, const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues, bool& rIsDamaging);

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
    double CalculateCharacteristicLength(GenericSmallStrainFemDemElement *pCurrentElement);

    /**
     * this computes VolumeForce of the element
     */
    Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

    /**
     * this computes the damage according to a exp softening
     */
    void CalculateExponentialDamage(double& rDamage, const double DamageParameter, const double UniaxialStress, const double InitialThrehsold);


    // ************** Methods to compute the tangent constitutive tensor via numerical derivation ************** 
    /**
     * this computes the Tangent tensor via numerical derivation (perturbations)
     */
    void CalculateTangentTensor(Matrix& rTangentTensor,const Vector& rStrainVectorGP,const Vector& rStressVectorGP,const Matrix& rElasticMatrix, ConstitutiveLaw::Parameters& rValues);

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
    // *****************************************************************

    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specializations of element must specify the actual interface to the integration points!
     * Note, that these functions expect a std::vector of values for the specified variable type that
     * contains a value for each integration point!
     * SetValueOnIntegrationPoints: set the values for given Variable.
     * GetValueOnIntegrationPoints: get the values for given Variable.
     * these methods are: OPTIONAL
     */
	void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
	void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
	void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate variables on Integration points.
     * This gives access to variables computed in the constitutive law on each integration point.
     * Specialisations of element must specify the actual interface to the integration points!
     * Note, that these functions expect a std::vector of values for the specified variable type that
     * contains a value for each integration point!
     * CalculateValueOnIntegrationPoints: calculates the values of given Variable.
     * these methods are: OPTIONAL
     */
	void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
	void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
	void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

protected:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

	Vector mNonConvergedThresholds;     // Equivalent stress
	Vector mThresholds;                 // Stress mThreshold on edge
	Vector mDamages;                    // Converged Damage on each edge
	Vector mNonConvergedDamages;        // Damages at edges of "i" iteration
	double mThreshold = 0.0;            // Converged Threshold
	double mDamage = 0.0;               // Converged Damage

    // Vector to storage the neigh elements sharing a certain edge
    std::vector<std::vector<Element*>> mEdgeNeighboursContainer;

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
}; // Class GenericSmallStrainFemDemElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif