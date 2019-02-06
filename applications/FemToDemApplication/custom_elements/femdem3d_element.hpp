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

#if !defined(KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED

#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{
class FemDem3DElement : public SmallDisplacementElement // Derived Element from SolidMechanics
{

  public:
  	static constexpr double tolerance = std::numeric_limits<double>::epsilon();
	/// Default constructors
	FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

	FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	FemDem3DElement(FemDem3DElement const &rOther);

	/// Destructor.
	virtual ~FemDem3DElement();

	/// Assignment operator.
	FemDem3DElement &operator=(FemDem3DElement const &rOther);

	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	FemDem3DElement()
	{
	}

	void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo);
	void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo);
	void InitializeNonLinearIteration(ProcessInfo &CurrentProcessInfo);
	void CalculateConstitutiveMatrix(Matrix &rConstitutiveMatrix, const double rYoungModulus,
									 const double rPoissonCoefficient);
	void CalculateDN_DX(Matrix &rDN_DX, int PointNumber);
	void CalculateInfinitesimalStrain(Vector &rStrainVector, const Matrix &rDN_DX);
	void CalculateStressVector(Vector &rStressVector, const Matrix &rConstitutiveMAtrix, const Vector &rInfinitesimalStrainVector);
	void CalculatePrincipalStresses(Vector &PrincipalStressVector, const Vector &StressVector);
	void FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo);
	void CalculateOnIntegrationPoints(const Variable<Vector> &rVariable, std::vector<Vector> &rOutput, const ProcessInfo &rCurrentProcessInfo);
	void CalculateOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rOutput, const ProcessInfo &rCurrentProcessInfo);
	void CalculateOnIntegrationPoints(const Variable<Matrix> &rVariable, std::vector<Matrix> &rOutput, const ProcessInfo &rCurrentProcessInfo);
	void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
							  ProcessInfo &rCurrentProcessInfo) override;
	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
	void AverageVector(Vector &rAverageVector, const Vector &v, const Vector &w);
	void GetValueOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues,
									 const ProcessInfo &rCurrentProcessInfo);
	void GetValueOnIntegrationPoints(const Variable<Vector> &rVariable,
									 std::vector<Vector> &rValues,
									 const ProcessInfo &rCurrentProcessInfo);
	void GetValueOnIntegrationPoints(const Variable<Matrix> &rVariable,
									 std::vector<Matrix> &rValues,
									 const ProcessInfo &rCurrentProcessInfo);
	void Get2MaxValues(Vector &MaxValues, double a, double b, double c);
	void Get2MinValues(Vector &MaxValues, double a, double b, double c);
	void IntegrateStressDamageMechanics(double& rThreshold,
										double &rDamage,
										const Vector &rStrainVector,
										const Vector &rStressVector,
										int Edge,
										double Length,
										bool& rIsDamaging);
	void ModifiedMohrCoulombCriterion(double& rThreshold, double &Damage, const Vector &StressVector, int cont, double L_char, bool& rIsDamaging);
	void RankineCriterion(double& rThreshold, double &Damage, const Vector &StressVector, int cont, double L_char, bool& rIsDamaging);
	void DruckerPragerCriterion(double& rThreshold, double &Damage, const Vector &StressVector, int cont, double L_char, bool& rIsDamaging);
	void SimoJuCriterion(double& rThreshold, double &Damage, const Vector &StrainVector, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void RankineFragileLaw(double& rThreshold, double &Damage, const Vector &StressVector, int cont, double L_char, bool& rIsDamaging);
	void ElasticLaw(double& rThreshold, double &Damage, const Vector &StressVector, int cont, double L_char, bool& rIsDamaging);


	// Stress Invariants in 3D
	double CalculateI1Invariant(const Vector& rStressVector);
	double CalculateI2Invariant(const Vector& rStressVector);
	double CalculateI3Invariant(const Vector& rStressVector);
	void CalculateDeviatorVector(Vector& rDeviator, const Vector& rStressVector, const double I1);
	double CalculateJ2Invariant(const Vector& rDeviator);
	double CalculateJ3Invariant(const Vector& rDeviator);
	void CalculateIntegratedStressVector(Vector &rIntegratedStressVector, const Vector rStressVector, const double Damage)
	{
		noalias(rIntegratedStressVector) = (1.0 - Damage) * rStressVector;
	}
	// Lode's angle
	double CalculateLodeAngle(double J2, double J3);
	// Converged values
	void SetThreshold(double af, int cont) { mThresholds[cont] = af; }
	double GetThreshold(int cont) { return mThresholds[cont]; }
	Vector GetThresholds() { return mThresholds; }
	Vector GetDamages() { return mDamages; }
	void SetThreshold(double af) { mThreshold = af; }
	double GetThreshold() { return mThreshold; }
	void SetConvergedDamage(double af) { mDamage = af; }
	double GetDamage() { return mDamage; }
	void SetConvergedDamages(double af, int cont) { mDamages[cont] = af; }
	double GetConvergedDamages(int cont) { return mDamages[cont]; }
	// Non Converged values
	void SetNonConvergedDamages(double af, int cont) { mNonConvergedDamages[cont] = af; }
	double GetNonConvergedDamages(int cont) { return mNonConvergedDamages[cont]; }
	// Characteristic length Calculations
	Vector CalculateCharacteristicLengths();
	void CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo);
	Vector &CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN);
	// Functions to calculate the Constitutive tangent tensor by numerical derivation
	double GetMaxValue(Vector Strain);
	double GetMaxAbsValue(const Vector& rVector);
	double GetMinAbsValue(const Vector& rVector);
	void CalculateDeformationMatrix(Matrix &rB, const Matrix &rDN_DX);
	double CalculateElementalDamage(const Vector& EdgeDamages);
	// Fills mEdgeNeighboursContainer
	void ComputeEdgeNeighbours(ProcessInfo &rCurrentProcessInfo);
	// Storages mEdgeNeighboursContainer
	void SaveEdgeNeighboursContainer(const std::vector<std::vector<Element *>>& toSave) { mEdgeNeighboursContainer = toSave; }
	std::vector<Element *> GetEdgeNeighbourElements(int edge) { return mEdgeNeighboursContainer[edge]; }
	void CalculateAverageStressOnEdge(Vector &AverageVector, const std::vector<Element *>& VectorOfElems);
	void CalculateAverageStrainOnEdge(Vector &AverageVector, const std::vector<Element *>& VectorOfElems);
	void SetNodeIndexes(Matrix &M) // Defines the numbering of the edges with the corresponding nodes
	{
		M.resize(6, 2);

		M(0, 0) = 0;
		M(0, 1) = 1;
		M(1, 0) = 0;
		M(1, 1) = 2;
		M(2, 0) = 0;
		M(2, 1) = 3;
		M(3, 0) = 1;
		M(3, 1) = 2;
		M(4, 0) = 1;
		M(4, 1) = 3;
		M(5, 0) = 2;
		M(5, 1) = 3;
	}
	double GetNumberOfEdges() {return mNumberOfEdges;}
	void SetValueOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;
	void SetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

	// methods to compute numerical tangent tensor
	void CalculateTangentTensor(
		Matrix& TangentTensor,
		const Vector& rStrainVectorGP,
		const Vector& rStressVectorGP,
		const Matrix& rElasticMatrix);
	void CalculatePerturbation(
		const Vector& rStrainVectorGP,
		double& rPerturbation,
		const int Component);
	void PerturbateStrainVector(
		Vector& rPerturbedStrainVector,
		const Vector& rStrainVectorGP,
		const double Perturbation,
		const int Component);
	void IntegratePerturbedStrain(
		Vector& rPerturbedStressVector,
		const Vector& rPerturbedStrainVector,
		const Matrix& rElasticMatrix);
	void AssignComponentsToTangentTensor(
		Matrix& rTangentTensor,
		const Vector& rDeltaStress,
		const double Perturbation,
		const int Component);
	void InitializeInternalVariablesAfterMapping();
	void UpdateDataBase();
	void CalculateExponentialDamage(
		double& rDamage,
		const double DamageParameter,
		const double UniaxialStress,
		const double InitialThrehsold);

  protected:

	int mNumberOfEdges;
	Vector mNonConvergedThresholds;   // Equivalent stress
	Vector mThresholds; // Stress mThreshold on edge
	Vector mDamages; // Converged mDamage on each edge
	Vector mNonConvergedDamages; // mDamages on edges of "i" iteration
	double mThreshold = 0.0;
	double mDamage = 0.0; // Converged mDamage

  private:
	// Vector to storage the neigh elements sharing a certain edge
	std::vector<std::vector<Element*>> mEdgeNeighboursContainer;

}; // Class FemDem3DElement

} // Namespace Kratos
#endif // KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED  defined