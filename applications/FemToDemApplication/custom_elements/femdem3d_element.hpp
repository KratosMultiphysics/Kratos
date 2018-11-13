//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED

#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{
class FemDem3DElement : public SmallDisplacementElement // Derived Element from SolidMechanics
{

  public:
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

	// *************** Methods Alejandro Cornejo ***************
	//**********************************************************
	static constexpr double tolerance = std::numeric_limits<double>::epsilon();

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
							  ProcessInfo &rCurrentProcessInfo);

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

	void IntegrateStressDamageMechanics(Vector &rIntegratedStress,
										double &rDamage, const Vector &StrainVector, const Vector &StressVector, int cont, double L_char);

	void ModifiedMohrCoulombCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void RankineCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void DruckerPragerCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void SimoJuCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StrainVector, const Vector &StressVector, int cont, double L_char);
	void RankineFragileLaw(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);

	void TangentModifiedMohrCoulombCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);

	// Stress Invariants in 3D
	double CalculateI1Invariant(Vector StressVector);
	double CalculateI2Invariant(const Vector StressVector);
	double CalculateI3Invariant(const Vector StressVector);
	void CalculateDeviatorVector(Vector &rDeviator, const Vector StressVector, const double I1);
	double CalculateJ2Invariant(const Vector Deviator);
	double CalculateJ3Invariant(const Vector Deviator);

	void CalculateIntegratedStressVector(Vector &rIntegratedStressVector, const Vector rStressVector, const double Damage)
	{
		rIntegratedStressVector = (1 - Damage) * rStressVector;
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

	void SetConvergedEquivalentStress(double af) { mF_sigma = af; }
	double GetConvergedEquivalentStress() { return mF_sigma; }

	void SetConvergedEquivalentStress(double af, int cont) { mF_sigmas[cont] = af; }
	double GetConvergedEquivalentStress(int cont) { return mF_sigmas[cont]; }

	void SetConvergedDamages(double af, int cont) { mDamages[cont] = af; }
	double GetConvergedDamages(int cont) { return mDamages[cont]; }

	// Non Converged values
	void SetNonConvergedDamages(double af, int cont) { mNonConvergedDamages[cont] = af; }
	double GetNonConvergedDamages(int cont) { return mNonConvergedDamages[cont]; }

	void SetNonConvergedDamages(double af) { mNonConvergedDamage = af; }
	double GetNonConvergedDamage() { return mNonConvergedDamage; }

	void SetNonConvergedEquivalentStress(double af, int cont) { mNonConvergedFsigmas[cont] = af; }
	double GetNonConvergedEquivalentStress(int cont) { return mNonConvergedFsigmas[cont]; }

	void SetNonConvergedEquivalentStress(double af) { mNonConvergedFsigma = af; }
	double GetNonConvergedEquivalentStress() { return mNonConvergedFsigma; }

	void ResetNonConvergedVars()
	{
		this->SetNonConvergedDamages(0.0);
		this->SetNonConvergedEquivalentStress(0.0);

		for (int cont = 0; cont < 6; cont++)
		{
			this->SetNonConvergedDamages(0, cont);
			this->SetNonConvergedEquivalentStress(0, cont);
		}
	}

	// Characteristic length Calculations
	void Set_l_char(double af, int cont) { mL_char[cont] = af; }
	double Get_l_char(int cont) { return mL_char[cont]; }
	void CalculateLchar();

	// Auxiliar functions...
	void IterationPlus() { iteration++; }
	int GetIteration() { return iteration; }
	void SetToZeroIteration() { iteration = 0; }
	//void AssignSmoothedStress(Element& Elem);

	void CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo);
	Vector &CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN);

	// Functions to calculate the Constitutive tangent tensor by numerical derivation
	double GetMaxValue(Vector Strain);
	double GetMaxAbsValue(Vector Strain);
	double GetMinAbsValue(Vector Strain);

	void SetStressVector(Vector toStressVector)
	{
		toStressVector.resize(6);
		mStressVector = toStressVector;
	}
	Vector GetStressVector() { return mStressVector; }

	void SetStrainVector(Vector toStrainVector)
	{
		toStrainVector.resize(6);
		mStrainVector = toStrainVector;
	}
	Vector GetStrainVector() { return mStrainVector; }

	void SetIntegratedStressVector(Vector toIntegratedStressVector)
	{
		toIntegratedStressVector.resize(6);
		mIntegratedStressVector = toIntegratedStressVector;
	}
	Vector GetIntegratedStressVector() { return mIntegratedStressVector; }

	void SetBMatrix(Matrix toBMatrix) { mB = toBMatrix; }
	Matrix GetBMatrix() { return mB; }

	void CalculateDeformationMatrix(Matrix &rB, const Matrix &rDN_DX);

	// Fills mEdgeNeighboursContainer
	void ComputeEdgeNeighbours(ProcessInfo &rCurrentProcessInfo);

	// Storages mEdgeNeighboursContainer
	void SaveEdgeNeighboursContainer(std::vector<std::vector<Element *>> toSave) { mEdgeNeighboursContainer = toSave; }
	std::vector<Element *> GetEdgeNeighbourElements(int edge) { return mEdgeNeighboursContainer[edge]; }

	void CalculateAverageStressOnEdge(Vector &AverageVector, const std::vector<Element *> VectorOfElems);
	void CalculateAverageStrainOnEdge(Vector &AverageVector, const std::vector<Element *> VectorOfElems);
	void AddDEMContactForces(Vector &NodalRHS);

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

	double CalculateElementalDamage(const Vector &EdgeDamages);
	double GetNumberOfEdges() {return mNumberOfEdges;}

	void SetValueOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

	void SetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

  protected:

	 // Each component == Each edge
	 int mNumberOfEdges;
	 Vector mF_sigmas;   // Equivalent stress
	 Vector mThresholds; // Stress mThreshold on edge
	 Vector mDamages; // Converged mDamage on each edge
	 
	 Vector mNonConvergedDamages; // mDamages on edges of "i" iteration
	 Vector mNonConvergedFsigmas; // Equivalent stress of "i" iteration

  private:

	int iteration = 0;
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

	double mThreshold = 0.0;
	double mF_sigma = 0.0;
	double mDamage = 0.0;			 // Converged mDamage
	double mNonConvergedFsigma = 0.0;
	double mNonConvergedDamage = 0.0; // mDamage of the element of "i" iteration

	Vector mL_char; // Characteristic length on each edge

	Vector mStressVector = ZeroVector(voigt_size);
	Vector mStrainVector = ZeroVector(voigt_size);
	Vector mIntegratedStressVector = ZeroVector(voigt_size);
	Matrix mB = ZeroMatrix(voigt_size, dimension * number_of_nodes);

	// Vector to storage the neigh elements sharing a certain edge
	std::vector<std::vector<Element*>> mEdgeNeighboursContainer;

}; // Class FemDem3DElement

} // Namespace Kratos
#endif // KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED  defined