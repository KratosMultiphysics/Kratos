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

#if !defined(KRATOS_FEMDEM2D_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM2D_ELEMENT_H_INCLUDED

#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{

class FemDem2DElement : public SmallDisplacementElement // Derived Element from SolidMechanics
{

  public:
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();
	/// Default constructors
	FemDem2DElement(IndexType NewId, GeometryType::Pointer pGeometry);

	FemDem2DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	FemDem2DElement(FemDem2DElement const &rOther);

	/// Destructor.
	virtual ~FemDem2DElement();

	/// Assignment operator.
	FemDem2DElement &operator=(FemDem2DElement const &rOther);

	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	FemDem2DElement()
	{
	}

	// *************** Methods Alejandro Cornejo ***************
	//**********************************************************
	void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override;
	void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo) override;
	void InitializeNonLinearIteration(ProcessInfo &CurrentProcessInfo) override;
	void CalculateConstitutiveMatrix(Matrix &rConstitutiveMatrix, const double &rYoungModulus,
									 const double &rPoissonCoefficient);

	void CalculateDN_DX(Matrix &rDN_DX, int PointNumber);

	void CalculateInfinitesimalStrain(Vector &rStrainVector, const Matrix &rDN_DX);

	void CalculateStressVector(Vector &rStressVector, const Matrix &rConstitutiveMAtrix, const Vector &rInfinitesimalStrainVector);

	void CalculatePrincipalStress(Vector &PrincipalStressVector, const Vector StressVector);

	void FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo) override;

	void CalculateOnIntegrationPoints(const Variable<Vector> &rVariable, std::vector<Vector> &rOutput, const ProcessInfo &rCurrentProcessInfo) override;
	void CalculateOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rOutput, const ProcessInfo &rCurrentProcessInfo) override;
	void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
							  ProcessInfo &rCurrentProcessInfo) override;
	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

	void AverageVector(Vector &rAverageVector, const Vector &v, const Vector &w);

	void GetValueOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues,
									 const ProcessInfo &rCurrentProcessInfo) override;

	void GetValueOnIntegrationPoints(const Variable<Vector> &rVariable,
									 std::vector<Vector> &rValues,
									 const ProcessInfo &rCurrentProcessInfo) override;

	void Get2MaxValues(Vector &MaxValues, const double a, const double b, const double c);
	void Get2MinValues(Vector &MaxValues, const double a, const double b, const double c);
	void InitializeInternalVariablesAfterMapping();

	void IntegrateStressDamageMechanics(double &rThreshold,
										double &Damage, 
										const Vector& StrainVector, 
										const Vector& StressVector, 
										const int cont, 
										const double L_char, 
										bool& rIsDamaging);

	void ModifiedMohrCoulombCriterion(double& rThreshold, double &rDamage, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void RankineCriterion(double& rThreshold, double &rDamage, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void DruckerPragerCriterion(double& rThreshold, double &rDamage, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void SimoJuCriterion(double& rThreshold, double &rDamage, const Vector &StrainVector, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void RankineFragileLaw(double& rThreshold, double &rDamage, const Vector &StressVector, const int cont, const double L_char, bool& rIsDamaging);
	void ElasticLaw(double& rThreshold,double &rDamage, const Vector &rStressVector, const int Edge, const double Length, bool& rIsDamaging);

	void CalculateExponentialDamage(
		double& rDamage,
		const double DamageParameter,
		const double UniaxialStress,
		const double InitialThrehsold);
		
	// Stress Invariants in 2D
	double CalculateI1Invariant(double sigma1, double sigma2);
	double CalculateJ2Invariant(double sigma1, double sigma2);
	double CalculateJ3Invariant(double sigma1, double sigma2, double I1);

	void CalculateIntegratedStressVector(Vector &rIntegratedStressVector, const Vector& rStressVector, const double Damage)
	{
		noalias(rIntegratedStressVector) = (1.0 - Damage) * rStressVector;
	}

	// Lode's angle
	double CalculateLodeAngle(double J2, double J3);
	void UpdateDataBase();
	void CalculateAverageStressOnEdge(const Element* Neighbour, const Element* CurrentElement, Vector& rAverageStress);
	void CalculateAverageStrainOnEdge(const Element* Neighbour, const Element* CurrentElement, Vector& rAverageStrain);

	// Converged values
	void SetThreshold(double af, int cont) { mThresholds[cont] = af; }
	double GetThreshold(int cont) { return mThresholds[cont]; }

	Vector GetThresholds() { return mThresholds; }
	Vector GetDamages() { return mDamages; }

	void SetThreshold(double af) { mThreshold = af; }
	double GetThreshold() { return mThreshold; }

	void SetConvergedDamage(double af) { mDamage = af; }
	double GetConvergedDamage() { return mDamage; }

	void SetConvergedDamages(double af, int cont) { mDamages[cont] = af; }
	double GetConvergedDamages(int cont) { return mDamages[cont]; }

	// Non Converged values
	void SetNonConvergedDamages(double af, int cont) { mNonConvergedDamages[cont] = af; }
	double GetNonConvergedDamages(int cont) { return mNonConvergedDamages[cont]; }


	// Characteristic length Calculations
	double CalculateCharacteristicLength(FemDem2DElement *CurrentElement, const Element &NeibElement, int cont);

	void CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo);
	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
	Vector &CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN);
	double CalculateElementalDamage(const Vector& rEdgeDamages);

	// Functions to calculate the Constitutive tangent tensor by numerical derivation
	double GetMaxValue(const Vector& rValues);
	double GetMaxAbsValue(const Vector& rValues);
	double GetMinAbsValue(const Vector& rValues);


	void CalculateDeformationMatrix(Matrix &rB, const Matrix &rDN_DX);

	void SetValueOnIntegrationPoints(
		const Variable<double> &rVariable, std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

	void SetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable, std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

	void CalculateTangentTensor(
		Matrix& TangentTensor,
		const Vector& rStrainVectorGP,
		const Vector& rStressVectorGP,
		const Matrix& rElasticMatrix);
	void CalculateSecondOrderTangentTensor(
		Matrix& TangentTensor,
		const Vector& rStrainVectorGP,
		const Vector& rStressVectorGP,
		const Matrix& rElasticMatrix);
	void CalculateSecondOrderCentralDifferencesTangentTensor(
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
	void AssignComponentsToSecondOrderTangentTensor(
		Matrix& rTangentTensor,
		const Vector& rGaussPointStress,
		const Vector& rPerturbedStress,
		const Vector& rTwicePerturbedStress,
		const double Perturbation,
		const int Component);
	void AssignComponentsToSecondOrderCentralDifferencesTangentTensor(
		Matrix& rTangentTensor,
		const Vector& rGaussPointStress,
		const Vector& rPerturbedStress,
		const Vector& rTwicePerturbedStress,
		const double Perturbation,
		const int Component);

       private:
	int mNumberOfEdges;
	// Each component == Each edge
	Vector mThresholds; // Stress mThreshold on edge
	Vector mNonConvergedThresholds; // Stress mThreshold on edge
	double mThreshold = 0.0;
	Vector mDamages; // Converged mDamage on each edge
	double mDamage = 0.0;			 // Converged mDamage
	Vector mNonConvergedDamages; // mDamages on edges of "i" iteration

}; // Class FemDem2DElement

} // Namespace Kratos
#endif 
