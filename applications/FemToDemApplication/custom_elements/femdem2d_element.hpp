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

#if !defined(KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED)
#define KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{

class FemDem2DElement : public SmallDisplacementElement // Derived Element from SolidMechanics
{

  public:
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

	void AverageVector(Vector &rAverageVector, const Vector &v, const Vector &w);

	void GetValueOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues,
									 const ProcessInfo &rCurrentProcessInfo) override;

	void GetValueOnIntegrationPoints(const Variable<Vector> &rVariable,
									 std::vector<Vector> &rValues,
									 const ProcessInfo &rCurrentProcessInfo) override;

	void Get2MaxValues(Vector &MaxValues, double a, double b, double c);
	void Get2MinValues(Vector &MaxValues, double a, double b, double c);

	void IntegrateStressDamageMechanics(Vector &rIntegratedStress,
										double &Damage, const Vector StrainVector, const Vector StressVector, int cont, double L_char);

	void ModifiedMohrCoulombCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void RankineCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void DruckerPragerCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);
	void SimoJuCriterion(Vector &rIntegratedStress, double &Damage, const Vector &StrainVector, const Vector &StressVector, int cont, double L_char);
	void RankineFragileLaw(Vector &rIntegratedStress, double &Damage, const Vector &StressVector, int cont, double L_char);

	// Stress Invariants in 2D
	double CalculateI1Invariant(double sigma1, double sigma2);
	double CalculateJ2Invariant(double sigma1, double sigma2);
	double CalculateJ3Invariant(double sigma1, double sigma2, double I1);

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
	double GetConvergedDamage() { return mDamage; }

	void SetConverged_f_sigma(double af) { mF_sigma = af; }
	double GetConverged_f_sigma() { return mF_sigma; }

	void SetConvergedEquivalentStresses(double af, int cont) { mF_sigmas[cont] = af; }
	double GetConvergedEquivalentStresses(int cont) { return mF_sigmas[cont]; }

	void SetConvergedDamages(double af, int cont) { mDamages[cont] = af; }
	double GetConvergedDamages(int cont) { return mDamages[cont]; }

	// Non Converged values
	void Set_NonConvergeddamages(double af, int cont) { mNonConvergedDamages[cont] = af; }
	double GetNonConvergedDamages(int cont) { return mNonConvergedDamages[cont]; }

	void Set_NonConvergeddamage(double af) { mNonConvergedDamage = af; }
	double Get_NonConvergeddamage() { return mNonConvergedDamage; }

	void SetNonConvergedEquivalentStress(double af, int cont) { mNonConvergedFsigmas[cont] = af; }
	double GetNonConvergedEquivalentStress(int cont) { return mNonConvergedFsigmas[cont]; }

	void SetNonConvergedEquivalentStress(double af) { mNonConvergedFsigma = af; }
	double GetNonConvergedEquivalentStress() { return mNonConvergedFsigma; }

	void ResetNonConvergedVars()
	{
		this->Set_NonConvergeddamage(0.0);
		this->SetNonConvergedEquivalentStress(0.0);
		for (unsigned int cont = 0; cont < 3; cont++) {
			this->Set_NonConvergeddamages(0, cont);
			this->SetNonConvergedEquivalentStress(0, cont);
		}
	}

	// Characteristic length Calculations
	void SetCharacteristicLength(double af, int cont) { mL_char[cont] = af; }
	double Get_l_char(int cont) { return mL_char[cont]; }
	double CalculateLchar(FemDem2DElement *CurrentElement, const Element &NeibElement, int cont);

	// Auxiliar functions...
	void IterationPlus() { iteration++; }
	int GetIteration() { return iteration; }
	void SetToZeroIteration() { iteration = 0; }

	void CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo);
	Vector &CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN);

	// Functions to calculate the Constitutive tangent tensor by numerical derivation
	double GetMaxValue(Vector Strain);
	double GetMaxAbsValue(Vector Strain);
	double GetMinAbsValue(Vector Strain);

	void SetIntegratedStressVector(Vector toIntegratedStressVector)
	{
		toIntegratedStressVector.resize(3);
		mIntegratedStressVector = toIntegratedStressVector;
	}
	
	Vector GetIntegratedStressVector() { return mIntegratedStressVector; }

	void SetBMatrix(Matrix toBMatrix)
	{
		toBMatrix.resize(3, 6);
		mB = toBMatrix;
	}

	Matrix GetBMatrix() { return mB; }

	void CalculateDeformationMatrix(Matrix &rB, const Matrix &rDN_DX);

	void SetValueOnIntegrationPoints(
		const Variable<double> &rVariable, std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

	void SetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable, std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo) override;

       private:
	int iteration = 0;

	// Each component == Each edge
	Vector mF_sigmas = ZeroVector(3);   // Mohr-Coulomb equivalent stress
	Vector mThresholds = ZeroVector(3); // Stress mThreshold on edge

	double mThreshold = 0.0;
	double mF_sigma = 0.0;

	Vector mDamages = ZeroVector(3); // Converged mDamage on each edge
	double mDamage = 0.0;			 // Converged mDamage

	Vector mNonConvergedDamages = ZeroVector(3); // mDamages on edges of "i" iteration
	Vector mNonConvergedFsigmas = ZeroVector(3); // Equivalent stress of "i" iteration

	double mNonConvergedFsigma = 0.0;
	double mNonConvergedDamage = 0.0; // mDamage of the element of "i" iteration

	Vector mL_char = ZeroVector(3); // Characteristic length on each edge

	Vector mIntegratedStressVector = ZeroVector(3);
	Matrix mB = ZeroMatrix(3, 6);

}; // Class FemDem2DElement

} // Namespace Kratos
#endif // KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED  defined
