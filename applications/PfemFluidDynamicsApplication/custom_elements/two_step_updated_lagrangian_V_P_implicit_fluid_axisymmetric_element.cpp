//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

/**
 * TODO:
 * external forces to be updated for axisymmetric case... check ZT1
 */

// System includes

// External includes
#include <math.h>

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_axisymmetric_element.h"
#include "custom_utilities/element_utilities.hpp"
//#include "includes/cfd_variables.h"

namespace Kratos {

/*
   * public TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim> functions
   */

template <unsigned int TDim>
Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const {

	TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

	NewElement.SetData(this->GetData());
	NewElement.SetFlags(this->GetFlags());

	return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement(NewElement));
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::Initialize() {
	KRATOS_TRY;

	// If we are restarting, the constitutive law will be already defined
	if (mpConstitutiveLaw == nullptr) {
		const Properties& r_properties = this->GetProperties();
		KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
		    << "In initialization of Element " << this->Info() << ": No CONSTITUTIVE_LAW defined for property "
		    << r_properties.Id() << "." << std::endl;
		mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();
	}

	KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) {}

template <unsigned int TDim>
int TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) {
	KRATOS_TRY;

	// Base class checks for positive Jacobian and Id > 0
	int ierr = Element::Check(rCurrentProcessInfo);
	if (ierr != 0)
		return ierr;

	// Check that all required variables have been registered
	if (VELOCITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "VELOCITY Key is 0. Check that the application was correctly registered.", "");
	if (ACCELERATION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "ACCELERATION Key is 0. Check that the application was correctly registered.", "");
	if (PRESSURE.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "PRESSURE Key is 0. Check that the application was correctly registered.", "");
	if (BODY_FORCE.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "BODY_FORCE Key is 0. Check that the application was correctly registered.", "");
	if (DENSITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "DENSITY Key is 0. Check that the application was correctly registered.", "");
	if (DYNAMIC_VISCOSITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.", "");
	if (DELTA_TIME.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument,
		                   "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");

	// Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
	for (unsigned int i = 0; i < this->GetGeometry().size(); ++i) {
		if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing BODY_FORCE variable on solution step data for node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing DENSITY variable on solution step data for node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
			KRATOS_THROW_ERROR(std::invalid_argument,
			                   "missing DYNAMIC_VISCOSITY variable on solution step data for node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
		    this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
		    this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ",
			                   this->GetGeometry()[i].Id());
		if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ",
			                   this->GetGeometry()[i].Id());
	}

	// If this is a 2D problem, check that nodes are in XY plane
	if (this->GetGeometry().WorkingSpaceDimension() == 2) {
		for (unsigned int i = 0; i < this->GetGeometry().size(); ++i) {
			if (this->GetGeometry()[i].Z() != 0.0)
				KRATOS_THROW_ERROR(std::invalid_argument,
				                   "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
		}
	}

	// Consitutive law checks
	const auto& r_properties = this->GetProperties();
	const auto& r_geometry = this->GetGeometry();
	const SizeType dimension = r_geometry.WorkingSpaceDimension();
	mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);

	// Verify that the constitutive law exists
	KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
	    << "Constitutive law not provided for property " << r_properties.Id() << std::endl;

	// Verify that the constitutive law has the correct dimension
	const SizeType strain_size = r_properties.GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

	if (dimension == 2) {
		KRATOS_ERROR_IF(strain_size < 3 || strain_size > 4)
		    << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 "
		       "(el id = ) "
		    << this->Id() << std::endl;
	} else {
		KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! "
		                                         "expected strain size is 6 (el id = ) "
		                                      << this->Id() << std::endl;
	}

	// Check constitutive law
	return mpConstitutiveLaw->Check(r_properties, r_geometry, rCurrentProcessInfo);

	return ierr;

	KRATOS_CATCH("");
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::ComputeBulkReductionCoefficient(
    MatrixType MassMatrix,
    MatrixType StiffnessMatrix,
    double& meanValueStiff,
    double& bulkCoefficient,
    double timeStep) {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
	IndexType FirstRow = 0;
	IndexType FirstCol = 0;
	double meanValueMass = 0;
	double countStiff = 0;
	double countMass = 0;
	for (SizeType j = 0; j < NumNodes; ++j) {
		for (SizeType i = 0; i < NumNodes; ++i) {
			meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol));
			meanValueStiff += fabs(StiffnessMatrix(FirstRow, FirstCol + 1));
			meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol));
			meanValueStiff += fabs(StiffnessMatrix(FirstRow + 1, FirstCol + 1));
			// Update Counter
			countStiff += 4.0;

			meanValueMass += fabs(MassMatrix(FirstRow, FirstCol));
			meanValueMass += fabs(MassMatrix(FirstRow + 1, FirstCol + 1));
			// Update Counter
			countMass += 2.0;

			FirstRow += 2;
		}
		FirstRow = 0;
		FirstCol += 2;
	}

	meanValueStiff /= countStiff;
	meanValueMass /= countMass;

	if (meanValueMass != 0 && meanValueStiff != 0) {
		bulkCoefficient = meanValueMass * 4 / (timeStep * meanValueStiff);
	} else {
		std::cout << " DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << " coordinates " << this->GetGeometry()[0].X() << " " << this->GetGeometry()[0].Y() << std::endl;
		std::cout << " MeanValueMass=" << meanValueMass;
		std::cout << "\t MeanValueMaterial= " << meanValueStiff;
		bulkCoefficient = timeStep;
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::ComputeBulkMatrixLump(
    Matrix& BulkMatrix,
    const double Weight) {

	for (SizeType i = 0; i < this->GetGeometry().PointsNumber(); ++i) {
		// LHS contribution
		double Mij = Weight / (1.0 + TDim);
		BulkMatrix(i, i) += Mij;
	}
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::ComputeBoundLHSMatrix(
    Matrix& BoundLHSMatrix,
    const ShapeFunctionsType& rN,
    const double Weight) {

	GeometryType& rGeom = this->GetGeometry();
	const double coeff = 1.0 / 3.0;

	if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)) {
		if (rGeom[0].IsNot(INLET))
			BoundLHSMatrix(0, 0) += Weight * coeff;
		if (rGeom[1].IsNot(INLET))
			BoundLHSMatrix(1, 1) += Weight * coeff;
	}
	if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {
		if (rGeom[0].IsNot(INLET))
			BoundLHSMatrix(0, 0) += Weight * coeff;
		if (rGeom[2].IsNot(INLET))
			BoundLHSMatrix(2, 2) += Weight * coeff;
	}
	if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {
		if (rGeom[1].IsNot(INLET))
			BoundLHSMatrix(1, 1) += Weight * coeff;
		if (rGeom[2].IsNot(INLET))
			BoundLHSMatrix(2, 2) += Weight * coeff;
	}
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::ComputeBoundRHSVectorComplete(
    VectorType& BoundRHSVector,
    const double TimeStep,
    const double BoundRHSCoeffAcc,
    const double BoundRHSCoeffDev,
    const VectorType SpatialDefRate) {

    GeometryType& rGeom = this->GetGeometry();
	const double coeff = 1.0 / 3.0;
	const double timeFactor = 0.5 / TimeStep;

	if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)) {
		array_1d<double, 3> AccA(3, 0.0);
		array_1d<double, 3> AccB(3, 0.0);
		array_1d<double, 3> MeanAcc(3, 0.0);
		array_1d<double, 3> NormalVector(3, 0.0);

		this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

		double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

		noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(MeanAcc) = 0.5 * (AccA + AccB);

		const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

		if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
			BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

		if (rGeom[1].IsNot(INLET))
			BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
	}

	if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {

		array_1d<double, 3> AccA(3, 0.0);
		array_1d<double, 3> AccB(3, 0.0);
		array_1d<double, 3> MeanAcc(3, 0.0);
		array_1d<double, 3> NormalVector(3, 0.0);

		this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

		double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

		noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(MeanAcc) = 0.5 * (AccA + AccB);

		const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

		if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
			BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

		if (rGeom[2].IsNot(INLET))
			BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
	}

	if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {

		array_1d<double, 3> AccA(3, 0.0);
		array_1d<double, 3> AccB(3, 0.0);
		array_1d<double, 3> MeanAcc(3, 0.0);
		array_1d<double, 3> NormalVector(3, 0.0);

		this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

		double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

		noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(MeanAcc) = 0.5 * (AccA + AccB);

		const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

		if (rGeom[1].IsNot(INLET))
			BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

		if (rGeom[2].IsNot(INLET))
			BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
	}
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::ComputeBoundRHSVector(
    VectorType& BoundRHSVector,
    const ShapeFunctionsType& rN,
    const double TimeStep,
    const double BoundRHSCoeffAcc,
    const double BoundRHSCoeffDev) {

    GeometryType& rGeom = this->GetGeometry();
	array_1d<double, 3> AccA(3, 0.0);
	array_1d<double, 3> AccB(3, 0.0);

	const double factor = 0.5 / TimeStep;
	const double one_third = 1.0 / 3.0;

	if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)) {
		noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);

        const array_1d<double, 3>& NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
		const array_1d<double, 3>& NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
		if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
			BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
		if (rGeom[1].IsNot(INLET))
			BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
	}
	if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {
		noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
		const array_1d<double, 3>& NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
		const array_1d<double, 3>& NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
		if (rGeom[0].IsNot(INLET))
			BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
		if (rGeom[2].IsNot(INLET))
			BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
	}
	if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)) {
		noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
		noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
		const array_1d<double, 3>& NormalA = rGeom[1].FastGetSolutionStepValue(NORMAL);
		const array_1d<double, 3>& NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
		if (rGeom[1].IsNot(INLET))
			BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
		if (rGeom[2].IsNot(INLET))
			BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::CalculateTauFIC(
    double& Tau,
    double ElemSize,
    const double Density,
    const double Viscosity,
    const ProcessInfo& rCurrentProcessInfo) {

    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
	if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME)) {
		DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
	}

	double MeanVelocity = 0;
	this->CalcMeanVelocity(MeanVelocity, 0);

	Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

	if (MeanVelocity == 0) {
		Tau = 0;
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::AddStabilizationMatrixLHS(
    MatrixType& rLeftHandSideMatrix,
    Matrix& BulkAccMatrix,
    const ShapeFunctionsType& rN,
    const double Weight) {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

	if (BulkAccMatrix.size1() != NumNodes)
		BulkAccMatrix.resize(NumNodes, NumNodes);

	BulkAccMatrix = ZeroMatrix(NumNodes, NumNodes);
	for (SizeType i = 0; i < NumNodes; ++i) {
		// LHS contribution
		for (SizeType j = 0; j < NumNodes; ++j) {
			double Mij = 0.0;
			Mij = Weight * rN[i] * rN[j];
			BulkAccMatrix(i, j) += Mij;
		}
	}
	rLeftHandSideMatrix += BulkAccMatrix;
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::ComputeStabLaplacianMatrix(
    MatrixType& StabLaplacianMatrix,
    const ShapeFunctionDerivativesType& rDN_DX,
    const double Weight) {

	// LHS contribution
	const SizeType NumNodes = this->GetGeometry().PointsNumber();
	for (SizeType i = 0; i < NumNodes; ++i) {
		for (SizeType j = 0; j < NumNodes; ++j) {
			double Lij = 0.0;
			for (SizeType d = 0; d < TDim; ++d) {
				Lij += rDN_DX(i, d) * rDN_DX(j, d);
			}
			StabLaplacianMatrix(i, j) += Weight * Lij;
		}
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::AddStabilizationNodalTermsLHS(
    MatrixType& rLeftHandSideMatrix,
    const double Tau,
    const double Weight,
    const ShapeFunctionDerivativesType& rDN_DX,
    const SizeType i) {

    // LHS contribution
	const SizeType NumNodes = this->GetGeometry().PointsNumber();
	for (SizeType j = 0; j < NumNodes; ++j) {
		double Lij = 0.0;
		for (SizeType d = 0; d < TDim; ++d) {
			Lij += rDN_DX(i, d) * rDN_DX(j, d);
		}
		Lij *= Tau;

		rLeftHandSideMatrix(i, j) += Weight * Lij;
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::AddStabilizationNodalTermsRHS(
    VectorType& rRightHandSideVector,
    const double Tau,
    const double Density,
    const double Weight,
    const ShapeFunctionDerivativesType& rDN_DX,
    const SizeType i) {

	double RHSi = 0;
	if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION)) { // it must be checked once at the begining only
		array_1d<double, 3>& VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

		for (SizeType d = 0; d < TDim; ++d) {
			RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
		}
	}
	rRightHandSideVector[i] += Weight * RHSi;
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::InitializeElementalVariables(ElementalVariables& rElementalVariables) {
	KRATOS_TRY;

	const unsigned int voigtsize = 4; //SWITCH_TO_AXISYM

    rElementalVariables.voigtsize = voigtsize;

	rElementalVariables.ConstitutiveMatrix = ZeroMatrix(voigtsize, voigtsize);

	rElementalVariables.DetFgrad = 1.0;

	rElementalVariables.DetFgradVel = 1.0;

	//rElementalVariables.DeviatoricInvariant = 1.0;

	rElementalVariables.EquivalentStrainRate = 1.0;

	rElementalVariables.VolumetricDefRate = 1.0;

	rElementalVariables.SpatialDefRate = ZeroVector(voigtsize);

	//rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize, false);

	//noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);

	rElementalVariables.Fgrad = ZeroMatrix(TDim, TDim);

	rElementalVariables.InvFgrad = ZeroMatrix(TDim, TDim);

	rElementalVariables.FgradVel = ZeroMatrix(TDim, TDim);

	rElementalVariables.InvFgradVel = ZeroMatrix(TDim, TDim);

	rElementalVariables.SpatialVelocityGrad = ZeroMatrix(TDim, TDim);

	rElementalVariables.MeanPressure = 0;

	rElementalVariables.CurrentTotalCauchyStress = ZeroVector(voigtsize);

	rElementalVariables.UpdatedTotalCauchyStress = ZeroVector(voigtsize);

	rElementalVariables.CurrentDeviatoricCauchyStress = ZeroVector(voigtsize);

	rElementalVariables.UpdatedDeviatoricCauchyStress = ZeroVector(voigtsize);

    // MZ added
    rElementalVariables.B = ZeroMatrix(voigtsize, TDim);

	KRATOS_CATCH("");
}

template <>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::CalcElasticPlasticCauchySplitted(
    ElementalVariables& rElementalVariables,
    double TimeStep,
    unsigned int g,
    const ProcessInfo& rCurrentProcessInfo,
    double& Density,
    double& DeviatoricCoeff,
    double& VolumetricCoeff) {

    // reset previous UpdatedDeviatoricCauchyStress and UpdatedTotalCauchyStress
    // this should be not required since this two vectors are overrided each time
    rElementalVariables.UpdatedDeviatoricCauchyStress = ZeroVector(rElementalVariables.voigtsize);
    rElementalVariables.UpdatedTotalCauchyStress = ZeroVector(rElementalVariables.voigtsize);

	mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
	auto constitutive_law_values =
	    ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

	Flags& constitutive_law_options = constitutive_law_values.GetOptions();
	constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
	constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

	const Vector& r_shape_functions = row((this->GetGeometry()).ShapeFunctionsValues(GeometryData::GI_GAUSS_2), g);
	constitutive_law_values.SetShapeFunctionsValues(r_shape_functions);
	constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
	constitutive_law_values.SetStressVector(rElementalVariables.UpdatedDeviatoricCauchyStress);
	constitutive_law_values.SetConstitutiveMatrix(rElementalVariables.ConstitutiveMatrix);

	// Temporary workaround, to be updated
	auto r_geometry = this->GetGeometry();
	r_geometry[0].SetValue(THETA_MOMENTUM, 0.5);

	mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

	rElementalVariables.UpdatedTotalCauchyStress[0] = rElementalVariables.UpdatedDeviatoricCauchyStress[0] + rElementalVariables.MeanPressure;
	rElementalVariables.UpdatedTotalCauchyStress[1] = rElementalVariables.UpdatedDeviatoricCauchyStress[1] + rElementalVariables.MeanPressure;
	rElementalVariables.UpdatedTotalCauchyStress[2] = rElementalVariables.UpdatedDeviatoricCauchyStress[2] + rElementalVariables.MeanPressure;
    rElementalVariables.UpdatedTotalCauchyStress[3] = rElementalVariables.UpdatedDeviatoricCauchyStress[3];

	const double time_step = rCurrentProcessInfo[DELTA_TIME];
	const double bulk_modulus = this->GetProperties()[BULK_MODULUS];

	DeviatoricCoeff = rElementalVariables.ConstitutiveMatrix(rElementalVariables.voigtsize - 1, rElementalVariables.voigtsize - 1);
    VolumetricCoeff = bulk_modulus * time_step;
	Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

	this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
	this->mMaterialVolumetricCoefficient = VolumetricCoeff;
	this->mMaterialDensity = Density;
    //KRATOS_WATCH(DeviatoricCoeff);
    //KRATOS_WATCH(VolumetricCoeff);
    //KRATOS_WATCH(Density);
}

template <>
bool TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>::CalcCompleteStrainRateAxisym(
    ElementalVariables& rElementalVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const ShapeFunctionDerivativesType& rDN_DX,
    const double theta,
    unsigned int PointNumber) {

    const GeometryType &rGeometry = GetGeometry();
    const auto& N_values = rGeometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);//mThisIntegrationMethod
    const Vector& r_shape_functions = row(N_values, PointNumber);

    double current_radius = PfemElementUtilities::CalculateRadius(r_shape_functions,rGeometry,Current); // try const, but not sure if it works
    double initial_radius = PfemElementUtilities::CalculateRadius(r_shape_functions,rGeometry,Initial);

    PfemElementUtilities::CalculateLinearDeformationMatrix(rElementalVariables.B, rGeometry, rDN_DX, rElementalVariables.voigtsize, PointNumber, r_shape_functions, current_radius);

	double current_radial_velocity = 0.0;
    for (unsigned int i = 0; i < rGeometry.PointsNumber(); i++) {
		current_radial_velocity += r_shape_functions[i] * rGeometry[i].FastGetSolutionStepValue(VELOCITY_X, 0);
	}

	bool computeElement = true;
	const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    const SizeType NumNodes = rGeometry.PointsNumber();
	const SizeType LocalSize = dimension * NumNodes;
	VectorType NodePosition = ZeroVector(LocalSize);
	VectorType VelocityValues = ZeroVector(LocalSize);
	VectorType RHSVelocities = ZeroVector(LocalSize);
	this->GetPositions(NodePosition, rCurrentProcessInfo, theta);
	this->GetVelocityValues(RHSVelocities, 0);
	RHSVelocities *= theta;
	this->GetVelocityValues(VelocityValues, 1);
	RHSVelocities += VelocityValues * (1.0 - theta);

    rElementalVariables.Fgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.FgradVel = ZeroMatrix(dimension, dimension);
    for (SizeType i = 0; i < dimension; i++) {
		for (SizeType j = 0; j < dimension; j++) {
			for (SizeType k = 0; k < NumNodes; k++) {
				rElementalVariables.Fgrad(i, j) += NodePosition[dimension * k + i] * rDN_DX(k, j);
				rElementalVariables.FgradVel(i, j) += RHSVelocities[dimension * k + i] * rDN_DX(k, j);
			}
		}
	}

	//Inverse
    rElementalVariables.InvFgrad = ZeroMatrix(dimension,dimension);
    rElementalVariables.DetFgrad = 0.0; // should be not strictly required
	rElementalVariables.DetFgrad = MathUtils<double>::Det(rElementalVariables.Fgrad);
    MathUtils<double>::InvertMatrix2(rElementalVariables.Fgrad,
	                                 rElementalVariables.InvFgrad,
	                                 rElementalVariables.DetFgrad);
    rElementalVariables.DetFgrad *= current_radius / initial_radius;
    //KRATOS_WATCH(rElementalVariables.DetFgrad);

	//rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension, false); // test_SpatialVelocityGrad
	rElementalVariables.SpatialVelocityGrad = ZeroMatrix(dimension,dimension);
    rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

	rElementalVariables.VolumetricDefRate = 0;
	for (SizeType i = 0; i < dimension; i++) {
		rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
	}
    rElementalVariables.VolumetricDefRate += current_radial_velocity/ current_radius;//rElementalVariables.SpatialVelocityGrad(2, 2);// SWITCH_TO_AXISYM

    rElementalVariables.SpatialDefRate = ZeroVector(rElementalVariables.voigtsize); //should be not strictly required
    rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
    rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
    rElementalVariables.SpatialDefRate[2] = current_radial_velocity/current_radius;
    rElementalVariables.SpatialDefRate[3] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));

    // temporary workaround - start
	//GeometryType& r_Geometry = this->GetGeometry();
	////const SizeType NumNodes = this->GetGeometry().PointsNumber();
	//SizeType FirstRow = 0;
	//for (SizeType i = 0; i < NumNodes; ++i) {
	//	r_Geometry[i].FastGetSolutionStepValue(YIELD_SHEAR, 0) += rElementalVariables.SpatialDefRate[0];//r_Geometry[i].FastGetSolutionStepValue
    //    r_Geometry[i].FastGetSolutionStepValue(YIELDED, 0) += rElementalVariables.SpatialDefRate[1];//rGeometry[i].FastGetSolutionStepValue
	//	r_Geometry[i].FastGetSolutionStepValue(ALPHA_PARAMETER, 0) += rElementalVariables.SpatialDefRate[2];
    //    FirstRow += 2;
	//}
	// temporary workaround - start

    // TODO: to be checked for axisym case
	rElementalVariables.EquivalentStrainRate = sqrt((2.0 * rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
	                                                 2.0 * rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
	                                                 4.0 * rElementalVariables.SpatialDefRate[3] * rElementalVariables.SpatialDefRate[3]));

	// cross checking with AF axisym implementation
	Vector AxiSymMDGreenLagrangeMaterial = ZeroVector(rElementalVariables.voigtsize);
	AxiSymMDGreenLagrangeMaterial[0] = rElementalVariables.FgradVel(0, 0) * rElementalVariables.Fgrad(0, 0) + rElementalVariables.FgradVel(1, 0) * rElementalVariables.Fgrad(1, 0);
	AxiSymMDGreenLagrangeMaterial[1] = rElementalVariables.FgradVel(1, 1) * rElementalVariables.Fgrad(1, 1) + rElementalVariables.FgradVel(0, 1) * rElementalVariables.Fgrad(0, 1);
	AxiSymMDGreenLagrangeMaterial[3] = (rElementalVariables.FgradVel(0, 0) * rElementalVariables.Fgrad(0, 1) + rElementalVariables.FgradVel(1, 0) * rElementalVariables.Fgrad(1, 1) +
	                                    rElementalVariables.FgradVel(0, 1) * rElementalVariables.Fgrad(0, 0) + rElementalVariables.FgradVel(1, 1) * rElementalVariables.Fgrad(1, 0)) *
	                                   0.5;

	Vector AxiSymSpatialDefRate = ZeroVector(rElementalVariables.voigtsize);
	AxiSymSpatialDefRate[0] = (rElementalVariables.InvFgrad(0, 0) * AxiSymMDGreenLagrangeMaterial[0] * rElementalVariables.InvFgrad(0, 0) +
	                           rElementalVariables.InvFgrad(1, 0) * AxiSymMDGreenLagrangeMaterial[3] * rElementalVariables.InvFgrad(0, 0) * 2 +
	                           rElementalVariables.InvFgrad(1, 0) * AxiSymMDGreenLagrangeMaterial[1] * rElementalVariables.InvFgrad(1, 0));

	AxiSymSpatialDefRate[1] = (rElementalVariables.InvFgrad(0, 1) * AxiSymMDGreenLagrangeMaterial[0] * rElementalVariables.InvFgrad(0, 1) +
	                           rElementalVariables.InvFgrad(0, 1) * AxiSymMDGreenLagrangeMaterial[3] * rElementalVariables.InvFgrad(1, 1) * 2 +
	                           rElementalVariables.InvFgrad(1, 1) * AxiSymMDGreenLagrangeMaterial[1] * rElementalVariables.InvFgrad(1, 1));

	AxiSymSpatialDefRate[2] = current_radial_velocity / current_radius;

    AxiSymSpatialDefRate[3] = (rElementalVariables.InvFgrad(0, 0) * AxiSymMDGreenLagrangeMaterial[0] * rElementalVariables.InvFgrad(0, 1) +
	                           rElementalVariables.InvFgrad(0, 0) * AxiSymMDGreenLagrangeMaterial[3] * rElementalVariables.InvFgrad(1, 1) +
	                           rElementalVariables.InvFgrad(1, 0) * AxiSymMDGreenLagrangeMaterial[3] * rElementalVariables.InvFgrad(0, 1) +
	                           rElementalVariables.InvFgrad(1, 0) * AxiSymMDGreenLagrangeMaterial[1] * rElementalVariables.InvFgrad(1, 1));
    //KRATOS_WATCH(AxiSymSpatialDefRate);
    //KRATOS_WATCH(rElementalVariables.SpatialDefRate);
    //KRATOS_WATCH("+++++");

	return computeElement;
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {

	GeometryType& rGeom = this->GetGeometry();
	const unsigned int NumNodes = rGeom.PointsNumber();

	// Check sizes and initialize
	if (rLeftHandSideMatrix.size1() != NumNodes)
		rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

	rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

	if (rRightHandSideVector.size() != NumNodes)
		rRightHandSideVector.resize(NumNodes);

	rRightHandSideVector = ZeroVector(NumNodes);

	// Shape functions and integration points
	ShapeFunctionDerivativesArrayType DN_DX;
	Matrix NContainer;
	VectorType GaussWeights;
	this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
	const unsigned int NumGauss = GaussWeights.size();

	const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
	const double theta = this->GetThetaContinuity();
	const double ElemSize = this->ElementSize();

	ElementalVariables rElementalVariables;
	this->InitializeElementalVariables(rElementalVariables);

	const double maxViscousValueForStabilization = 0.1;
	const double Density = this->mMaterialDensity;
	double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
	double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;

	if (DeviatoricCoeff > maxViscousValueForStabilization) {
		DeviatoricCoeff = maxViscousValueForStabilization;
	}

	double Tau = 0;
	this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

	double totalVolume = 0.0;
	bool computeElement = false;
	double GaussWeight = 0.0;
    double current_radius = 0.0;
	// Loop on integration points
	for (unsigned int g = 0; g < NumGauss; g++) {
        this->InitializeElementalVariables(rElementalVariables);// it seems to be not necessary
		const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
        const ShapeFunctionsType& N = row(NContainer, g);
        current_radius = PfemElementUtilities::CalculateRadius(N,rGeom, Current);

        computeElement = this->CalcCompleteStrainRateAxisym(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta, g);

		GaussWeight = GaussWeights[g] * 2.0 * Globals::Pi * current_radius * rElementalVariables.DetFgrad;
		totalVolume += GaussWeight;

		if (computeElement == true) {

			double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);

			this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

			double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
			double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);

            this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

			double StabLaplacianWeight = Tau * GaussWeight;
			this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);

			for (SizeType i = 0; i < NumNodes; ++i) {
				// RHS contribution
				// Velocity divergence
				rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
				this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
			}
		}
	}
    //KRATOS_WATCH(totalVolume); // ok, checked!
	if (computeElement == true) {

		VectorType PressureValues = ZeroVector(NumNodes);
		VectorType PressureValuesForRHS = ZeroVector(NumNodes);
		this->GetPressureValues(PressureValuesForRHS, 0);
		//the LHS matrix up to now just contains the laplacian term and the bound term
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);

		this->GetPressureValues(PressureValues, 1);
		noalias(PressureValuesForRHS) += -PressureValues;
		MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
		MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes, NumNodes);
		double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
		double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;

		this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
		noalias(rLeftHandSideMatrix) += BulkMatrix;
		noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);

		this->GetPressureVelocityValues(PressureValues, 0);
		noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
		noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
		this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
		noalias(rLeftHandSideMatrix) += BulkMatrix;
		noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);

	} else {
		double lumpedBulkCoeff = totalVolume * Tau * Density / (TimeStep * VolumetricCoeff);// to be const
		MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes, NumNodes);
		this->ComputeBulkMatrixLump(BulkVelMatrixLump, lumpedBulkCoeff);
		noalias(rLeftHandSideMatrix) += BulkVelMatrixLump;
		VectorType PressureValues = ZeroVector(NumNodes);
		VectorType PressureValuesForRHS = ZeroVector(NumNodes);
		this->GetPressureValues(PressureValuesForRHS, 0);
		this->GetPressureValues(PressureValues, 1);
		noalias(PressureValuesForRHS) += -PressureValues;
		noalias(rRightHandSideVector) -= prod(BulkVelMatrixLump, PressureValuesForRHS);
	}
    //KRATOS_WATCH(rRightHandSideVector);
    //KRATOS_WATCH(rLeftHandSideMatrix);
    //KRATOS_WATCH("+++++++++++++++++++++++++++++");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::GetPressureAccelerationValues(Vector& rValues,
                                                                                              const int Step) {
	GeometryType& rGeom = this->GetGeometry();
	const SizeType NumNodes = rGeom.PointsNumber();

	if (rValues.size() != NumNodes)
		rValues.resize(NumNodes);

	for (SizeType i = 0; i < NumNodes; ++i) {
		rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION, Step);
	}
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::CalculateLocalMomentumEquations(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

    KRATOS_TRY;

	GeometryType& rGeom = this->GetGeometry();
	const unsigned int NumNodes = rGeom.PointsNumber();
	const unsigned int LocalSize = 2 * NumNodes;

	MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);
	MatrixType StiffnessMatrix = ZeroMatrix(LocalSize, LocalSize);

	// Check sizes and initialize
	if (rLeftHandSideMatrix.size1() != LocalSize)
		rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

	rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

	if (rRightHandSideVector.size() != LocalSize)
		rRightHandSideVector.resize(LocalSize);

	rRightHandSideVector = ZeroVector(LocalSize);

	// Shape functions and integration points
	ShapeFunctionDerivativesArrayType DN_DX;
	Matrix NContainer;
	VectorType GaussWeights;
	this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
	const unsigned int NumGauss = GaussWeights.size();
	const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

	const double theta = this->GetThetaMomentum();

	ElementalVariables rElementalVariables;
	this->InitializeElementalVariables(rElementalVariables);

	double totalVolume = 0;
	double MeanValueMass = 0;
	double Density = 0.0;
	double DeviatoricCoeff = 0;
	double VolumetricCoeff = 0;
    double current_radius = 0.0;

    //mz test
    rElementalVariables.test_rhs = ZeroVector(4);

	// Loop on integration points
	for (unsigned int g = 0; g < NumGauss; g++) {
        //this->InitializeElementalVariables(rElementalVariables);// required?!?
        const ShapeFunctionsType& N = row(NContainer, g);
        current_radius = PfemElementUtilities::CalculateRadius(N,rGeom, Current);
		double GaussWeight = GaussWeights[g] * 2.0 * Globals::Pi * current_radius;

		const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

		double Pressure = 0;
		double OldPressure = 0;

		this->EvaluateInPoint(Pressure, PRESSURE, N, 0);

		this->EvaluateInPoint(OldPressure, PRESSURE, N, 1);

		rElementalVariables.MeanPressure = OldPressure * (1 - theta) + Pressure * theta;

        bool computeElement = this->CalcCompleteStrainRateAxisym(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta, g);
        totalVolume += GaussWeight;// not sure if the mass matrix requires detFgrad
        GaussWeight *= rElementalVariables.DetFgrad;
        //KRATOS_WATCH(rElementalVariables.DetFgrad);


        //KRATOS_WATCH(rElementalVariables.SpatialDefRate);
		this->CalcElasticPlasticCauchySplitted(rElementalVariables, TimeStep, g, rCurrentProcessInfo, Density,
		                                       DeviatoricCoeff, VolumetricCoeff);

		if (computeElement == true) {

            this->AddInternalForces(rRightHandSideVector, rDN_DX, rElementalVariables, GaussWeight);

            this->AddExternalForces(rRightHandSideVector, Density, N, GaussWeight);

            // temporary workaround - start
			//SizeType FirstRow = 0;
			//array_1d<double, 3> VolumeAcceleration(3, 0.0);
            //Vector rRHSVector = ZeroVector(LocalSize);
			//this->EvaluateInPoint(VolumeAcceleration, VOLUME_ACCELERATION, N);
			//for (SizeType i = 0; i < NumNodes; ++i) {
			//	if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
			//		for (SizeType d = 0; d < 2; ++d) {
			//			rRHSVector[FirstRow + d] += GaussWeight * Density * N[i] * VolumeAcceleration[d];
			//		}
			//	}
			//	this->GetGeometry()[i].SetValue(YIELDED,0) = rRHSVector[FirstRow];
			//	FirstRow += 2;
			//}
            // temporary workaround - end

			this->ComputeCompleteTangentTerm(rElementalVariables, StiffnessMatrix, rDN_DX, DeviatoricCoeff, VolumetricCoeff, theta, GaussWeight);
		}
	}
    //KRATOS_WATCH(rElementalVariables.test_rhs);
    //KRATOS_WATCH(rRightHandSideVector);

	double lumpedDynamicWeight = totalVolume * Density;
	this->ComputeLumpedMassMatrix(MassMatrix, lumpedDynamicWeight, MeanValueMass);

	double BulkReductionCoefficient = 1.0;
	double MeanValueStiffness = 0.0;
	this->ComputeBulkReductionCoefficient(MassMatrix, StiffnessMatrix, MeanValueStiffness, BulkReductionCoefficient, TimeStep);

    //totalVolume = 0.0;
	if (BulkReductionCoefficient != 1.0) {
		VolumetricCoeff *= MeanValueMass * 2.0 / (TimeStep * MeanValueStiffness);
        //double reduction_coefficient = MeanValueMass * 2.0 / (TimeStep * MeanValueStiffness);
		StiffnessMatrix = ZeroMatrix(LocalSize, LocalSize);

        double GaussWeight =0.0;
		for (unsigned int g = 0; g < NumGauss; g++) {
            const ShapeFunctionsType& N = row(NContainer, g);
            const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
            double current_radius = PfemElementUtilities::CalculateRadius(N,rGeom, Current);

            // this line has been added since now the stiffness matrix is calculated as BT*C*B and B must be calculated for each GP in CalcCompleteStrainRateAxisym
            this->CalcCompleteStrainRateAxisym(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta, g);

            GaussWeight = GaussWeights[g] * 2.0 * Globals::Pi * current_radius;
			GaussWeight *= rElementalVariables.DetFgrad;
            //totalVolume += GaussWeight;
            this->ComputeCompleteTangentTerm(rElementalVariables, StiffnessMatrix, rDN_DX, DeviatoricCoeff, VolumetricCoeff, theta, GaussWeight);
		}
        //KRATOS_WATCH(totalVolume);//ok, checked
	}

	// Add residual of previous iteration to RHS
	VectorType VelocityValues = ZeroVector(LocalSize);
	VectorType AccelerationValues = ZeroVector(LocalSize);

	//2nd order
	this->GetAccelerationValues(AccelerationValues, 0);
	this->GetVelocityValues(VelocityValues, 0);
	noalias(AccelerationValues) += -2.0 * VelocityValues / TimeStep;
	this->GetVelocityValues(VelocityValues, 1);
	noalias(AccelerationValues) += 2.0 * VelocityValues / TimeStep; //these are negative accelerations
	noalias(rRightHandSideVector) += prod(MassMatrix, AccelerationValues);
	noalias(rLeftHandSideMatrix) += StiffnessMatrix + MassMatrix * 2 / TimeStep;
    KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<TDim>::CalculateDeformationGradient(
    const Matrix& rDN_DX,
    Matrix& rF,
    Matrix& rDeltaPosition,
    double& rCurrentRadius,
    double& rReferenceRadius) {

	KRATOS_TRY

	const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
	const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

	rF = identity_matrix<double>(3);//SWITCH_TO_AXISYM

	if (dimension == 2) {

		for (unsigned int i = 0; i < number_of_nodes; i++) {
			rF(0, 0) += rDeltaPosition(i, 0) * rDN_DX(i, 0);
			rF(0, 1) += rDeltaPosition(i, 0) * rDN_DX(i, 1);
			rF(1, 0) += rDeltaPosition(i, 1) * rDN_DX(i, 0);
			rF(1, 1) += rDeltaPosition(i, 1) * rDN_DX(i, 1);
		}

		rF(2, 2) = rCurrentRadius / rReferenceRadius;//SWITCH_TO_AXISYM

	} else if (dimension == 3) {

		std::cout << " AXISYMMETRIC case and 3D is not possible " << std::endl;
	} else {

		KRATOS_THROW_ERROR(std::invalid_argument, "something is wrong with the dimension", "");
	}

	KRATOS_CATCH("")
}

template class TwoStepUpdatedLagrangianVPImplicitFluidAxisymmetricElement<2>;

} // namespace Kratos
