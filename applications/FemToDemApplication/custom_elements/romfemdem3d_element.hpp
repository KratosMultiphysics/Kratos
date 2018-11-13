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

#if !defined(KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED)
#define KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/femdem3d_element.hpp"

namespace Kratos
{

class RomFemDem3DElement : public FemDem3DElement
{
  public:
	/// Default constructors
	RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

	RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	RomFemDem3DElement(RomFemDem3DElement const &rOther);

	/// Destructor.
	virtual ~RomFemDem3DElement();

	/// Assignment operator.
	RomFemDem3DElement &operator=(RomFemDem3DElement const &rOther);

	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	RomFemDem3DElement()
	{
	}

	void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo);
	void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo);
	void CalculatePredictiveStresses(const Vector &StrainVector);
	void CalculateAverageStressOnEdge(Vector &rAverageVector, const std::vector<Element *> VectorOfElems);
	Vector &CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN);
	void CalculateOnIntegrationPoints(const Variable<Matrix> &rVariable, std::vector<Matrix> &rOutput, const ProcessInfo &rCurrentProcessInfo);
	void GetValueOnIntegrationPoints(const Variable<Matrix> &rVariable, std::vector<Matrix> &rValues, const ProcessInfo &rCurrentProcessInfo);
	void GetValueOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues,
									 const ProcessInfo &rCurrentProcessInfo);
	void CalculateOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rOutput,
									  const ProcessInfo &rCurrentProcessInfo);

	void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo);

	// Functions for plasticity of the steel
	double GetKp() { return mKp; }
	void SetKp(double toKp) { mKp = toKp; }
	double GetCapap() { return mCapap; }
	void SetCapap(double tomCapap) { mCapap = tomCapap; }
	Vector GetPlasticDeformation() { return mPlasticDeformation; }
	void SetPlasticDeformation(Vector toPlasticDeformation) { mPlasticDeformation = toPlasticDeformation; }

	// non converged values
	double GetNonConvergedKp() { return mNonConvergedKp; }
	void SetNonConvergedKp(double toNonConvergedKp) { mNonConvergedKp = toNonConvergedKp; }
	double GetNonConvergedCapap() { return mNonConvergedCapap; }
	void SetNonConvergedCapap(double tomNonConvergedCapap) { mNonConvergedCapap = tomNonConvergedCapap; }
	Vector GetNonConvergedPlasticDeformation() { return mNonConvergedPlasticDeformation; }
	void SetNonConvergedPlasticDeformation(Vector toNonConvergedPlasticDeformation) { mNonConvergedPlasticDeformation = toNonConvergedPlasticDeformation; }

	void ResetNonConvergedVarsPlast()
	{
		Vector Zero = ZeroVector(6);
		SetNonConvergedKp(0.0);
		SetNonConvergedCapap(0.0);
		SetNonConvergedPlasticDeformation(Zero);
	}

	void UpdateAndSaveInternalVariables()
	{
		this->SetKp(this->GetNonConvergedKp());
		this->SetCapap(this->GetNonConvergedCapap());
		this->SetPlasticDeformation(this->GetNonConvergedPlasticDeformation());
	}

	// Methods
	void IntegrateStressPlasticity(Vector &rIntegratedStress, const Vector &PredictiveStress, const Matrix &C);
	void CalculatePlasticParameters(const Vector &StressVector, double &rYield, double &rKp,
									double &rPlasticDenominator, Vector &rFluxVector, double &Capap, const Vector &PlasticStrainIncr, const Matrix &C);
	void VonMisesYieldCriterion(const Vector &StressVector, Vector &rDeviator, double &ryield, double &rJ2);
	void CalculateFluxVector(const Vector &StressVector, const Vector &rDeviator, const double J2, Vector &rFluxVector);
	void CalculateRFactors(const Vector &StressVector, double &r0, double &r1);
	void CalculatePlasticDissipation(const Vector &PredictiveSress, const double r0, const double r1,
									 const Vector &PlasticStrainInc, double &Capap, Vector &rHCapa);
	void CalculateEquivalentStressThreshold(const double Capap, const double r0, const double r1,
											double &rEquivalentStressThreshold, double &rSlope);
	void LinearCalculateThreshold(const double Capap, const double Gf, double &rEqThreshold, double &rSlope);
	void ExponentialCalculateThreshold(const double Capap, const double Gf, double &rEqThreshold, double &rSlope);
	void HardSoftCalculateThreshold(const double Capap, const double Gf, double &rEqThreshold, double &rSlope);
	void CalculateHardeningParameter(const Vector &FluxVector, const double SlopeThreshold,
									 const Vector &HCapa, double &HardeningParam);
	void CalculatePlasticDenominator(const Vector &FluxVector, const Matrix &ElasticConstMatrix,
									 const double HardeningParam, double &rPlasticDenominator);

  private:
	// Plasticity variables
	double mKp = 0.0;	// Eq threshold
	double mCapap = 0.0; // Normalized Plastic Dissipation
	Vector mPlasticDeformation = ZeroVector(6);

	double mNonConvergedKp = 0.0;
	double mNonConvergedCapap = 0.0;
	Vector mNonConvergedPlasticDeformation = ZeroVector(6);
};

} // Namespace Kratos
#endif // KRATOS_ROMFEMDEM3D_ELEMENT_H_INCLUDED  defined