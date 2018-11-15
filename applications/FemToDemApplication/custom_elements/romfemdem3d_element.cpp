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

#include "includes/define.h"
#include <string>
#include "includes/constitutive_law.h"
#include "custom_constitutive/zarate_law.hpp"
#include "femdem3d_element.hpp"
#include "romfemdem3d_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
//#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: FemDem3DElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: FemDem3DElement(NewId, pGeometry, pProperties)
{
	//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(RomFemDem3DElement const &rOther)
	: FemDem3DElement(rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
	//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

RomFemDem3DElement &RomFemDem3DElement::operator=(RomFemDem3DElement const &rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

	FemDem3DElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer RomFemDem3DElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	//NEEDED TO CREATE AN ELEMENT
	return Element::Pointer(new RomFemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer RomFemDem3DElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
	//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

	RomFemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	return Element::Pointer(new RomFemDem3DElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

RomFemDem3DElement::~RomFemDem3DElement()
{
}

void RomFemDem3DElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
	//*****************************
	KRATOS_TRY

	bool is_active = true;
	if (this->IsDefined(ACTIVE))
	{
		is_active = this->Is(ACTIVE);
	}

	// Inactive elements can have negative determinant of the Jacobian
	if (is_active == true)
	{
		//1.-Initialize sizes for the system components:
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		Vector StrainVector(voigt_size);
		noalias(StrainVector) = ZeroVector(voigt_size);
		Vector StressVector(voigt_size);
		noalias(StressVector) = ZeroVector(voigt_size);
		Matrix ConstitutiveMatrix(voigt_size, voigt_size);
		noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
		Matrix B(voigt_size, dimension * number_of_nodes);
		noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
		Matrix DN_DX(number_of_nodes, dimension);
		noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

		//deffault values for the infinitessimal theory
		double detF = 1;
		Matrix F(dimension, dimension);
		noalias(F) = identity_matrix<double>(dimension);

		//3.-Calculate elemental system:

		//reading integration points
		const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//get the shape functions [N] (for the order of the default integration method)
		const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

		//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

		//calculate delta position (here coincides with the current displacement)
		Matrix DeltaPosition(number_of_nodes, dimension);
		noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

		//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0].resize(dimension, dimension, false);
		noalias(J[0]) = ZeroMatrix(dimension, dimension);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

		// Loop Over Integration Points
		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			Matrix InvJ(dimension, dimension);
			noalias(InvJ) = ZeroMatrix(dimension, dimension);
			double detJ = 0;
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

			if (detJ < 0)
			{
				this->Set(ACTIVE, false); // element alone inside a crack
				detJ = fabs(detJ);
			}

			if (detJ < 0)
				KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

			//compute cartesian derivatives for this integration point  [dN/dx_n]
			noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

			//set shape functions for this integration point
			Vector N = row(Ncontainer, PointNumber);

			//b.-compute infinitessimal strainof the composite
			this->CalculateInfinitesimalStrain(StrainVector, DN_DX);
			this->SetValue(STRAIN_VECTOR, StrainVector);

			// Compute predictive stresses for the concrete and steel
			this->CalculatePredictiveStresses(StrainVector);

			this->CalculateDeformationMatrix(B, DN_DX);
			this->SetBMatrix(B);
		}
	}

	KRATOS_CATCH("")
}

void RomFemDem3DElement::CalculatePredictiveStresses(const Vector &StrainVector)
{
	const double Ec = this->GetProperties()[YOUNG_MODULUS];
	const double Es = this->GetProperties()[YOUNG_MODULUS_STEEL];
	const double nuc = this->GetProperties()[POISSON_RATIO];
	const double nus = this->GetProperties()[POISSON_RATIO_STEEL];

	//const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int voigt_size = dimension * (dimension + 1) * 0.5;
	Matrix ConstitutiveMatrixConcrete = ZeroMatrix(voigt_size, voigt_size);
	Matrix ConstitutiveMatrixSteel = ZeroMatrix(voigt_size, voigt_size);

	// Elastic C
	this->CalculateConstitutiveMatrix(ConstitutiveMatrixConcrete, Ec, nuc);
	this->CalculateConstitutiveMatrix(ConstitutiveMatrixSteel, Es, nus);

	Vector StressVectorConcrete = prod(ConstitutiveMatrixConcrete, StrainVector);

	Vector StressVectorSteel;
	if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
	{
		Vector ElasticStrainVector = StrainVector - this->GetPlasticDeformation(); // E-Ep
		StressVectorSteel = prod(ConstitutiveMatrixSteel, ElasticStrainVector);
	}
	else
		StressVectorSteel = ZeroVector(voigt_size);

	// Predictive Stresses
	this->SetValue(CONCRETE_STRESS_VECTOR, StressVectorConcrete);
	this->SetValue(STEEL_STRESS_VECTOR, StressVectorSteel);
}

void RomFemDem3DElement::CalculateAverageStressOnEdge(Vector &rAverageVector, const std::vector<Element *> VectorOfElems)
{
	// Only averages the stress over the concrete part!!!!
	Vector CurrentElementStress = this->GetValue(CONCRETE_STRESS_VECTOR);
	rAverageVector = CurrentElementStress;
	int counter = 0;

	for (int elem = 0; elem < VectorOfElems.size(); elem++)
	{
		// Only take into account the active elements
		bool is_active = true;
		if (VectorOfElems[elem]->IsDefined(ACTIVE))
		{
			is_active = VectorOfElems[elem]->Is(ACTIVE);
		}

		if (is_active == true)
		{
			rAverageVector += VectorOfElems[elem]->GetValue(CONCRETE_STRESS_VECTOR);
			counter++;
		}
	}
	rAverageVector /= (counter + 1);
}

void RomFemDem3DElement::CalculateLocalSystem(
	MatrixType &rLeftHandSideMatrix,
	VectorType &rRightHandSideVector,
	ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
	unsigned int system_size = number_of_nodes * dimension;
	if (rLeftHandSideMatrix.size1() != system_size)
		rLeftHandSideMatrix.resize(system_size, system_size, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

	if (rRightHandSideVector.size() != system_size)
		rRightHandSideVector.resize(system_size, false);
	noalias(rRightHandSideVector) = ZeroVector(system_size);

	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0].resize(dimension, dimension, false);
	noalias(J[0]) = ZeroMatrix(dimension, dimension);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
	{
		const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
		Vector N = row(Ncontainer, PointNumber);

		double detJ = 0;
		Matrix InvJ(dimension, dimension);
		noalias(InvJ) = ZeroMatrix(dimension, dimension);
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;
		const Matrix &B = this->GetBMatrix();
		Vector IntegratedStressVectorConcrete = ZeroVector(voigt_size);
		Vector DamagesOnEdges = ZeroVector(6);

		// Loop over edges of the element
		for (int edge = 0; edge < 6; edge++)
		{
			std::vector<Element *> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);

			Vector AverageStressVectorConcrete, AverageStrainVectorConcrete, IntegratedStressVectorOnEdge;
			this->CalculateAverageStressOnEdge(AverageStressVectorConcrete, EdgeNeighbours);
			this->CalculateAverageStrainOnEdge(AverageStrainVectorConcrete, EdgeNeighbours);

			double DamageEdge = 0.0;
			double Lchar = this->Get_l_char(edge);

			// Integrate the stress on edge DAMAGE
			this->IntegrateStressDamageMechanics(IntegratedStressVectorOnEdge, DamageEdge,
												 AverageStrainVectorConcrete, AverageStressVectorConcrete, edge, Lchar);

			this->SetNonConvergedDamages(DamageEdge, edge);
			DamagesOnEdges[edge] = DamageEdge;

		} // End loop over edges

		// Compute elemental damage
		double damage_element = this->CalculateElementalDamage(DamagesOnEdges);
		if (damage_element >= 0.999)
		{
			damage_element = 0.999;
		}
		this->SetNonConvergedDamages(damage_element);

		const Vector &StressVectorConcrete = this->GetValue(CONCRETE_STRESS_VECTOR);
		noalias(IntegratedStressVectorConcrete) = (1.0 - damage_element) * StressVectorConcrete;
		this->SetIntegratedStressVector(IntegratedStressVectorConcrete);

		// Linear elastic const matrix concrete
		Matrix ConstitutiveMatrixConcrete = ZeroMatrix(voigt_size, voigt_size);
		const double Ec = this->GetProperties()[YOUNG_MODULUS];
		const double nuc = this->GetProperties()[POISSON_RATIO];
		this->CalculateConstitutiveMatrix(ConstitutiveMatrixConcrete, Ec, nuc);

		// Linear elastic const matrix steel
		Matrix ConstitutiveMatrixSteel = ZeroMatrix(voigt_size, voigt_size);
		const double Es = this->GetProperties()[YOUNG_MODULUS_STEEL];
		const double nus = this->GetProperties()[POISSON_RATIO_STEEL];
		this->CalculateConstitutiveMatrix(ConstitutiveMatrixSteel, Es, nus);

		const double k = this->GetProperties()[STEEL_VOLUMETRIC_PART];

		Matrix CompositeTangentMatrix = k * ConstitutiveMatrixSteel + (1.0 - k) * (1.0 - damage_element) * ConstitutiveMatrixConcrete;

		noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight * Matrix(prod(CompositeTangentMatrix, B))); // LHS

		Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

		// RHS Volumetric load
		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++)
			{
				rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
			}
		}

		//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
		Vector SteelStressVector = this->GetValue(STEEL_STRESS_VECTOR);
		Vector IntegratedSteelVector = ZeroVector(voigt_size);

		// Apply plasticity steel
		this->IntegrateStressPlasticity(IntegratedSteelVector, SteelStressVector, ConstitutiveMatrixSteel);
		this->SetValue(STEEL_STRESS_VECTOR, IntegratedSteelVector);

		Vector CompositeStressVector = k * IntegratedSteelVector + (1.0 - k) * IntegratedStressVectorConcrete;
		noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), CompositeStressVector);

		// Add nodal DEM forces
		Vector NodalRHS = ZeroVector(system_size);
		this->AddDEMContactForces(NodalRHS);

		// Add nodal contact forces from the DEM
		noalias(rRightHandSideVector) += NodalRHS;
	}
	KRATOS_CATCH("")
	//*****************************
}

Vector &RomFemDem3DElement::CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN)
{
	KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	if (rVolumeForce.size() != dimension)
		rVolumeForce.resize(dimension, false);

	noalias(rVolumeForce) = ZeroVector(dimension);

	for (unsigned int j = 0; j < number_of_nodes; j++)
	{
		if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION))
		{ // it must be checked once at the begining only
			array_1d<double, 3> &VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
			for (unsigned int i = 0; i < dimension; i++)
				rVolumeForce[i] += rN[j] * VolumeAcceleration[i];
		}
	}
	double k = this->GetProperties()[STEEL_VOLUMETRIC_PART];
	rVolumeForce *= (GetProperties()[DENSITY] * (1.0 - k) + GetProperties()[DENSITY_STEEL] * k);

	return rVolumeForce;

	KRATOS_CATCH("")
}

// 	TENSOR VARIABLES
void RomFemDem3DElement::CalculateOnIntegrationPoints(
	const Variable<Matrix> &rVariable,
	std::vector<Matrix> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	const unsigned int &integration_points_number = GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
	const unsigned int &dimension = GetGeometry().WorkingSpaceDimension();

	if (rOutput[0].size2() != dimension)
		rOutput[0].resize(dimension, dimension, false);

	if (rVariable == CONCRETE_STRESS_TENSOR)
	{
		rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(CONCRETE_STRESS_VECTOR));
	}
	else if (rVariable == STEEL_STRESS_TENSOR)
	{
		if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
		{
			rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STEEL_STRESS_VECTOR));
		}
		else
		{
			Matrix dummy;
			rOutput[0] = dummy;
		}
	}
	else if (rVariable == STRAIN_TENSOR)
	{
		rOutput[0] = MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
	}
	else if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED)
	{
		rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetIntegratedStressVector());
	}
}

// Tensor variables
void RomFemDem3DElement::GetValueOnIntegrationPoints(
	const Variable<Matrix> &rVariable,
	std::vector<Matrix> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == STRAIN_TENSOR)
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
	else if (rVariable == STEEL_STRESS_TENSOR)
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
	else if (rVariable == CONCRETE_STRESS_TENSOR)
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
	else if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED)
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
}

// **** plasticity methods *****
void RomFemDem3DElement::IntegrateStressPlasticity(
	Vector &rIntegratedStress,
	const Vector &PredictiveStress,
	const Matrix &C)
{ // ecuua
	int iteration = 0;
	const int iter_max = 9000;

	double Kp, Capap;
	Vector PlasticStrain;
	bool is_converged = false;

	// Get Converged values from the prev step
	Kp = this->GetKp();
	Capap = this->GetCapap();
	PlasticStrain = this->GetPlasticDeformation();

	double Yield, PlasticDenominator;
	Vector FluxVector = ZeroVector(6), PlasticStrainIncr = ZeroVector(6);

	// Compute Plastic variables
	this->CalculatePlasticParameters(PredictiveStress, Yield, Kp,
									 PlasticDenominator, FluxVector, Capap,
									 PlasticStrainIncr, C);

	double F = Yield - Kp;

	if (F <= std::abs(1.0e-8 * Kp))
	{ // Elastic
		rIntegratedStress = PredictiveStress;
		this->SetNonConvergedKp(Kp);
		this->SetNonConvergedCapap(Capap);
		this->SetNonConvergedPlasticDeformation(PlasticStrain);

		this->SetValue(EQUIVALENT_STRESS_VM, Yield);
	}
	else
	{ // Plastic case
		noalias(rIntegratedStress) = PredictiveStress;
		double DLambda;
		Vector DS = ZeroVector(6), DESIG = ZeroVector(6);

		while (is_converged == false && iteration <= iter_max)
		{
			DLambda = F * PlasticDenominator;
			// if (DLambda < 0.0)
			// 	DLambda = 0.0;

			noalias(PlasticStrainIncr) = DLambda * FluxVector;
			noalias(PlasticStrain) += PlasticStrainIncr;
			noalias(DS) = prod(C, PlasticStrainIncr);
			noalias(DESIG) -= DS;
			noalias(rIntegratedStress) -= DS;

			this->CalculatePlasticParameters(rIntegratedStress, Yield, Kp,
											 PlasticDenominator, FluxVector,
											 Capap, PlasticStrainIncr, C);

			F = Yield - Kp;

			if (F < std::abs(1.0e-8 * Kp)) // Has converged
			{
				is_converged = true;
				// Update Int Vars
				this->SetNonConvergedKp(Kp);
				this->SetNonConvergedCapap(Capap);
				this->SetNonConvergedPlasticDeformation(PlasticStrain);
				this->SetValue(EQUIVALENT_STRESS_VM, Yield);
			}
			else
				iteration++;
		}
		if (iteration == iter_max)
			KRATOS_ERROR << "Reached Max iterations inside Plasticity Loop" << std::endl;
	}
}

void RomFemDem3DElement::VonMisesYieldCriterion(
	const Vector &StressVector,
	Vector &rDeviator,
	double &ryield,
	double &rJ2)
{
	const double I1 = this->CalculateI1Invariant(StressVector);

	rDeviator = StressVector;
	const double Pmean = I1 / 3.0;

	rDeviator[0] -= Pmean;
	rDeviator[1] -= Pmean;
	rDeviator[2] -= Pmean;

	rJ2 = 0.5 * (rDeviator[0] * rDeviator[0] + rDeviator[1] * rDeviator[1] + rDeviator[2] * rDeviator[2]) +
		  (rDeviator[3] * rDeviator[3] + rDeviator[4] * rDeviator[4] + rDeviator[5] * rDeviator[5]);

	ryield = std::sqrt(3.0 * rJ2);
}

void RomFemDem3DElement::CalculatePlasticParameters(
	const Vector &StressVector,
	double &rYield,
	double &rKp,
	double &rPlasticDenominator,
	Vector &rFluxVector,
	double &rCapap,
	const Vector &PlasticStrainIncr,
	const Matrix &C)
{ // modaux
	Vector Deviator = ZeroVector(6), HCapa = ZeroVector(6);
	double J2 = 0.0, r0 = 0.0, r1 = 0.0, Slope = 0.0, HardeningParam = 0.0;

	this->VonMisesYieldCriterion(StressVector, Deviator, rYield, J2);
	this->CalculateFluxVector(StressVector, Deviator, J2, rFluxVector);
	this->CalculateRFactors(StressVector, r0, r1);
	this->CalculatePlasticDissipation(StressVector, r0, r1, PlasticStrainIncr, rCapap, HCapa);
	this->CalculateEquivalentStressThreshold(rCapap, r0, r1, rKp, Slope);
	this->CalculateHardeningParameter(rFluxVector, Slope, HCapa, HardeningParam);
	this->CalculatePlasticDenominator(rFluxVector, C, HardeningParam, rPlasticDenominator);
}

void RomFemDem3DElement::CalculateFluxVector(
	const Vector &StressVector,
	const Vector &rDeviator,
	const double J2,
	Vector &rFluxVector)
{
	// Only valid for Von Mises Yield Surf
	Vector AuxVec = ZeroVector(6);
	const double denomJ2 = 1.0 / (2.0 * std::sqrt(J2));

	for (int i = 0; i < AuxVec.size(); i++)
	{
		AuxVec[i] = rDeviator[i] * denomJ2;
	}

	AuxVec[3] *= 2.0;
	AuxVec[4] *= 2.0;
	AuxVec[5] *= 2.0;

	rFluxVector = std::sqrt(3.0) * AuxVec;
}

void RomFemDem3DElement::CalculateRFactors(const Vector &StressVector, double &r0, double &r1)
{
	Vector PrincipalStresses = ZeroVector(3);
	this->CalculatePrincipalStresses(PrincipalStresses, StressVector);

	double suma = 0.0, sumb = 0.0, sumc = 0.0;
	Vector SA = ZeroVector(3), SB = ZeroVector(3), SC = ZeroVector(3);

	for (int i = 0; i < 3; i++)
	{
		SA[i] = std::abs(PrincipalStresses[i]);
		SB[i] = 0.5 * (PrincipalStresses[i] + SA[i]);
		SC[i] = 0.5 * (-PrincipalStresses[i] + SA[i]);

		suma += SA[i];
		sumb += SB[i];
		sumc += SC[i];
	}

	if (suma != 0.0)
	{
		r0 = sumb / suma;
		r1 = sumc / suma;
	}
	else
	{
		r0 = sumb;
		r1 = sumc;
	}
}

void RomFemDem3DElement::CalculatePlasticDissipation(
	const Vector &PredictiveSress,
	const double r0,
	const double r1,
	const Vector &PlasticStrainInc,
	double &rCapap,
	Vector &rHCapa)
{
	const double E = this->GetProperties()[YOUNG_MODULUS_STEEL];
	const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];
	const double ft = this->GetProperties()[YIELD_STRESS_T_STEEL];
	const double n = fc / ft;
	const double Gf = this->GetProperties()[FRACTURE_ENERGY_STEEL];
	const double Gfc = Gf * std::pow(n, 2);
	const double Volume = this->GetGeometry().Volume();
	const double l_char = std::pow(Volume, 1.0 / 3.0);

	const double gf = Gf / l_char;
	const double gfc = Gfc / l_char;

	const double hlim = 2.0 * E * gfc / std::pow(fc, 2);
	if (l_char > hlim)
		KRATOS_THROW_ERROR(std::invalid_argument, " Characteristic length lower than minimum ", l_char)

	double Const0 = 0.0, Const1 = 0.0;
	if (gf > 0.000001)
		Const0 = r0 / gf;
	if (gfc > 0.000001)
		Const1 = r1 / gfc;

	const double Const = Const0 + Const1;
	double Dcapa = 0.0;

	for (int i = 0; i < 6; i++)
	{
		rHCapa[i] = Const * PredictiveSress[i];
		Dcapa += rHCapa[i] * PlasticStrainInc[i];
	}

	if (Dcapa<0.0 | Dcapa> 1.0)
		Dcapa = 0.0;
	rCapap += Dcapa;

	if (rCapap >= 1.0)
		rCapap = 0.9999;
}

void RomFemDem3DElement::CalculateEquivalentStressThreshold(
	const double Capap,
	const double r0,
	const double r1,
	double &rEquivalentStressThreshold,
	double &rSlope)
{
	const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];
	const double ft = this->GetProperties()[YIELD_STRESS_T_STEEL];
	const double n = fc / ft;
	Vector G = ZeroVector(2), EqTrhesholds = ZeroVector(2), Slopes = ZeroVector(2);
	G[0] = this->GetProperties()[FRACTURE_ENERGY_STEEL];
	G[1] = n * n * G[0];

	const int HardCurve = this->GetProperties()[HARDENING_LAW];

	for (int i = 0; i < 2; i++) // tension and compression curves
	{
		switch (HardCurve)
		{
		case 1:
			this->LinearCalculateThreshold(Capap, G[i], EqTrhesholds[i], Slopes[i]);
			break;

		case 2:
			this->ExponentialCalculateThreshold(Capap, G[i], EqTrhesholds[i], Slopes[i]);
			break;

		case 3:
			this->HardSoftCalculateThreshold(Capap, G[i], EqTrhesholds[i], Slopes[i]);
			break;
		
		default:
			KRATOS_ERROR << "Hardening law not defined..." << std::endl;
		}
	}

	rEquivalentStressThreshold = r0 * EqTrhesholds[0] + r1 * EqTrhesholds[1];
	rSlope = rEquivalentStressThreshold * ((r0 * Slopes[0] / EqTrhesholds[0]) + (r1 * Slopes[1] / EqTrhesholds[1]));
}

void RomFemDem3DElement::LinearCalculateThreshold(
	const double Capap,
	const double Gf,
	double &rEqThreshold,
	double &rSlope)
{ // Linear softening case!!
	const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];

	rEqThreshold = fc * std::sqrt(1.0 - Capap);
	rSlope = -0.5 * (std::pow(fc, 2.0) / (rEqThreshold));
}

void RomFemDem3DElement::ExponentialCalculateThreshold(
	const double Capap,
	const double Gf,
	double &rEqThreshold,
	double &rSlope)
{ // Exponential softening case!!
	const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];

	rEqThreshold = fc * (1.0 - Capap);
	rSlope = -0.5 * fc;
}

void RomFemDem3DElement::HardSoftCalculateThreshold(
	const double PlasticDissipation,
	const double Gf,
	double &rEqThreshold,
	double &rSlope)
{																						// Linear Hardening followed by exp softening
	const double initial_threshold = this->GetProperties()[YIELD_STRESS_C_STEEL];		// sikma
	const double peak_stress = this->GetProperties()[MAXIMUM_STRESS];					// sikpi
	const double peak_stress_position = this->GetProperties()[MAXIMUM_STRESS_POSITION]; // cappi [0,1]

	if (PlasticDissipation < 1.0)
	{
		const double Ro = std::sqrt(1.0 - initial_threshold / peak_stress);
		double alpha = std::log((1.0 - (1.0 - Ro) * (1.0 - Ro)) / ((3.0 - Ro) * (1.0 + Ro) * peak_stress_position));
		alpha = std::exp(alpha / (1.0 - peak_stress_position));
		const double Phi = std::pow((1.0 - Ro), 2) + ((3.0 - Ro) * (1.0 + Ro) * PlasticDissipation * (std::pow(alpha, (1.0 - PlasticDissipation))));

		rEqThreshold = peak_stress * (2.0 * std::sqrt(Phi) - Phi);
		
		rSlope = peak_stress * ((1.0 / std::sqrt(Phi)) - 1.0) * (3.0 - Ro) * (1.0 + Ro) * (std::pow(alpha, (1.0 - PlasticDissipation))) *
				 (1.0 - std::log(alpha) * PlasticDissipation);
	}
	else 
	{
		KRATOS_ERROR << "The Plastic Dissipation is greater that 1.0 ..." << std::endl;
	}
}

void RomFemDem3DElement::CalculateHardeningParameter(
	const Vector &FluxVector,
	const double SlopeThreshold,
	const Vector &HCapa,
	double &rHardeningParam)
{
	rHardeningParam = -SlopeThreshold;
	double aux = 0.0;

	for (int i = 0; i < 6; i++)
	{
		aux += HCapa[i] * FluxVector[i];
	}
	if (aux != 0.0)
		rHardeningParam *= aux;
}

void RomFemDem3DElement::CalculatePlasticDenominator(
	const Vector &FluxVector,
	const Matrix &ElasticConstMatrix,
	const double HardeningParam,
	double &rPlasticDenominator)
{ // only for isotropic hardening
	double A1 = 0.0, A2 = 0.0, A3 = 0.0;

	const Vector Dvect = prod(FluxVector, ElasticConstMatrix);

	for (int i = 0; i < 6; i++)
	{
		A1 += Dvect[i] * FluxVector[i];
	}

	A3 = HardeningParam;
	rPlasticDenominator = 1.0 / (A1 + A2 + A3);
}

void RomFemDem3DElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
	double current_equivalent_stress = 0.0, damage_element = 0.0;

	//Loop over edges
	for (int cont = 0; cont < 3; cont++)
	{
		this->SetConvergedDamages(this->GetNonConvergedDamages(cont), cont);
		this->SetConvergedEquivalentStress(this->GetNonConvergedEquivalentStress(cont), cont);
		current_equivalent_stress = this->GetConvergedEquivalentStress(cont);
		if (current_equivalent_stress > this->GetThreshold(cont))
		{
			this->SetThreshold(current_equivalent_stress, cont);
		}
	} // End Loop over edges

	damage_element = this->GetNonConvergedDamage();
	this->SetConvergedDamage(damage_element);

	if (damage_element >= 0.98)
	{
		if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
		{
			if (this->GetCapap() > 0.98)
			{
				this->Set(ACTIVE, false);
				double old_threshold = this->GetValue(STRESS_THRESHOLD);
				this->SetValue(INITIAL_THRESHOLD, old_threshold);
			}
		}
		else
		{
			this->Set(ACTIVE, false);
			double old_threshold = this->GetValue(STRESS_THRESHOLD);
			this->SetValue(INITIAL_THRESHOLD, old_threshold);
		}
	}

	this->ResetNonConvergedVars();
	this->SetToZeroIteration();

	// computation of the equivalent damage threshold and damage of the element for AMR mapping
	Vector thresholds = this->GetThresholds();

	Vector TwoMinValues;
	this->Get2MaxValues(TwoMinValues, thresholds[0], thresholds[1], thresholds[2]); // todo ojo con la funcion modificada
	double EqThreshold = 0.5 * (TwoMinValues[0] + TwoMinValues[1]);					// El menor o mayor?? TODO

	this->SetValue(STRESS_THRESHOLD, EqThreshold); // AMR
	this->SetThreshold(EqThreshold);
	this->SetValue(DAMAGE_ELEMENT, damage_element);

	// plasticity
	this->UpdateAndSaveInternalVariables();
	this->ResetNonConvergedVarsPlast();

	// Reset the nodal force flag for the next time step
	Geometry<Node<3>> &NodesElement = this->GetGeometry();
	for (int i = 0; i < 3; i++)
	{
		#pragma omp critical
		{
			NodesElement[i].SetValue(NODAL_FORCE_APPLIED, false);
		}
	}
}

// Double values
void RomFemDem3DElement::GetValueOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD || rVariable == PLASTIC_DISSIPATION_CAPAP)
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
}

// DOUBLE VARIABLES
void RomFemDem3DElement::CalculateOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT)
	{
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
		{
			rOutput[PointNumber] = double(this->GetValue(DAMAGE_ELEMENT));
		}
	}

	if (rVariable == IS_DAMAGED)
	{
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
		{
			rOutput[PointNumber] = double(this->GetValue(IS_DAMAGED));
		}
	}

	if (rVariable == STRESS_THRESHOLD)
	{
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
		{
			rOutput[PointNumber] = double(this->GetValue(STRESS_THRESHOLD));
		}
	}

	if (rVariable == PLASTIC_DISSIPATION_CAPAP)
	{
		rOutput.resize(1);
		if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
		{
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
			{
				rOutput[PointNumber] = double(this->GetCapap());
			}
		}
		else
		{
			double dummy = 0.0;
			rOutput[0] = dummy;
		}
	}
}
} // namespace Kratos