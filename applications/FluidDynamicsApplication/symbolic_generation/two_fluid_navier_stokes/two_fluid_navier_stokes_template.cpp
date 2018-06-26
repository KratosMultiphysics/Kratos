//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//

#include "two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"


namespace Kratos {


///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::TwoFluidNavierStokes(IndexType NewId,
    GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
TwoFluidNavierStokes<TElementData>::~TwoFluidNavierStokes() {}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(IndexType NewId, NodesArrayType const& ThisNodes, Properties::Pointer pProperties) const
{
    return Kratos::make_shared<TwoFluidNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer TwoFluidNavierStokes<TElementData>::Create(IndexType NewId, GeometryType::Pointer pGeom, Properties::Pointer pProperties) const
{
    return Kratos::make_shared<TwoFluidNavierStokes>(NewId, pGeom, pProperties);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateLocalSystem(
	MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo) {

	// Resize and intialize output
	if (rLeftHandSideMatrix.size1() != LocalSize)
		rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

	if (rRightHandSideVector.size() != LocalSize)
		rRightHandSideVector.resize(LocalSize, false);

	noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
	noalias(rRightHandSideVector) = ZeroVector(LocalSize);

	if (TElementData::ElementManagesTimeIntegration) {
		TElementData data;
		data.Initialize(*this, rCurrentProcessInfo);

		if (data.IsCut()) {
			GeometryType::Pointer p_geom = this->pGetGeometry();
            Matrix shape_functions_pos;
            Matrix shape_functions_neg;
            Matrix shape_functions_enr_pos;
            Matrix shape_functions_enr_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_pos;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_neg;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_pos;
            GeometryType::ShapeFunctionsGradientsType shape_derivatives_enr_neg;
            ComputeSplitting(
                data,
                shape_functions_pos,
                shape_functions_neg,
                shape_functions_enr_pos,
                shape_functions_enr_neg,
                shape_derivatives_pos,
			    shape_derivatives_neg,
                shape_derivatives_enr_pos,
                shape_derivatives_enr_neg);

			if (data.NumberOfDivisions == 1) {
				//cases exist when the element is not subdivided due to the characteristics of the provided distance
				//in this cases the element is treated as AIR or FLUID depending on the side
				Vector gauss_weights; 
                Matrix shape_functions; 
                ShapeFunctionDerivativesArrayType shape_derivatives; 
                this->CalculateGeometryData( 
                    gauss_weights, shape_functions, shape_derivatives); 
                const unsigned int number_of_gauss_points = gauss_weights.size(); 
				array_1d<double,NumNodes> Ncenter;
				for(unsigned int i=0; i<NumNodes; i++) Ncenter[i]=0.25;
				const double dgauss = inner_prod(data.Distance, Ncenter);
				for (unsigned int g = 0; g < number_of_gauss_points; g++) {
					data.UpdateGeometryValues(g,gauss_weights[g], row(shape_functions, g),
						shape_derivatives[g]);

                    data.CalculateDensityAtGaussPoint();

					if (dgauss > 0.0)
						data.CalculateAirMaterialResponse();
					else
						this->CalculateMaterialResponse(data);				

					this->AddTimeIntegratedSystem(
						data, rLeftHandSideMatrix, rRightHandSideVector);
				}
			}

			else {
				MatrixType Vtot;
				MatrixType Htot;
				MatrixType Kee_tot;
				VectorType rhs_ee_tot;
				Vtot.resize(NumNodes*(Dim + 1), NumNodes, false);
				Htot.resize(NumNodes, NumNodes*(Dim + 1), false);
				Kee_tot.resize(NumNodes, NumNodes, false);
				rhs_ee_tot.resize(NumNodes, false);

                noalias(Vtot) = ZeroMatrix(NumNodes*(Dim + 1), NumNodes);
                noalias(Htot) = ZeroMatrix(NumNodes, NumNodes*(Dim + 1));
                noalias(Kee_tot) = ZeroMatrix(NumNodes, NumNodes);
                noalias(rhs_ee_tot) = ZeroVector(NumNodes);

                for (unsigned int g_pos = 0; g_pos < data.w_gauss_pos_side.size(); g_pos ++) {
                    data.UpdateGeometryValues(
                        g_pos,
                        data.w_gauss_pos_side[g_pos],
                        row(shape_functions_pos, g_pos),
				        shape_derivatives_pos[g_pos],
                        row(shape_functions_enr_pos,g_pos),
                        shape_derivatives_enr_pos[g_pos]);
                    
                    data.CalculateDensityAtGaussPoint();
                    data.CalculateAirMaterialResponse();
                    this->AddTimeIntegratedSystem(
				 		data, rLeftHandSideMatrix, rRightHandSideVector);
                    
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

                for (unsigned int g_neg = 0; g_neg < data.w_gauss_neg_side.size(); g_neg ++) {
                    data.UpdateGeometryValues(
                        g_neg,
                        data.w_gauss_neg_side[g_neg],
                        row(shape_functions_neg, g_neg),
				        shape_derivatives_neg[g_neg],
                        row(shape_functions_enr_neg,g_neg),
                        shape_derivatives_enr_neg[g_neg]);
                    
                    data.CalculateDensityAtGaussPoint();
                    this->CalculateMaterialResponse(data);
                    this->AddTimeIntegratedSystem(
				 		data, rLeftHandSideMatrix, rRightHandSideVector);
                    
                    ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
                }

				CondenseEnrichment(data, rLeftHandSideMatrix,rRightHandSideVector,Htot,Vtot,Kee_tot, rhs_ee_tot);
			}
		}

		else {
			//Get Shape function data 
            Vector gauss_weights; 
            Matrix shape_functions; 
            ShapeFunctionDerivativesArrayType shape_derivatives; 
            this->CalculateGeometryData( 
              gauss_weights, shape_functions, shape_derivatives); 
            const unsigned int number_of_gauss_points = gauss_weights.size(); 
			// Iterate over integration points to evaluate local contribution
			for (unsigned int g = 0; g < number_of_gauss_points; g++) {

				data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),
					shape_derivatives[g]);

                data.CalculateDensityAtGaussPoint();

				if (data.IsAir())
					data.CalculateAirMaterialResponse();
				else
					this->CalculateMaterialResponse(data);

				this->AddTimeIntegratedSystem(
					data, rLeftHandSideMatrix, rRightHandSideVector);
			}

		}

		
	}
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::CalculateRightHandSide(
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo)
{
	MatrixType tmp;
	CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int TwoFluidNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) {

    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
std::string TwoFluidNavierStokes<TElementData>::Info() const {
    std::stringstream buffer;
    buffer << "TwoFluidNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::PrintInfo(
    std::ostream& rOStream) const {
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void TwoFluidNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}


template <>
void TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> >::ComputeGaussPointLHSContribution(
	TwoFluidNavierStokesData<2, 3>& rData, MatrixType& rLHS) {
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_2D

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}


template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointLHSContribution(
	TwoFluidNavierStokesData<3, 4>& rData, MatrixType& rLHS) {

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto& lhs = rData.lhs;

    //substitute_lhs_3D

    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;


}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointRHSContribution(
	TwoFluidNavierStokesData<2, 3>& rData, VectorType& rRHS) {

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& vmesh = rData.MeshVelocity;
    const auto& vconv = v - vmesh;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_2D

    noalias(rRHS) += rData.Weight * rhs;

}


template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointRHSContribution(
	TwoFluidNavierStokesData<3, 4>& rData, VectorType& rRHS) {

    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const auto& v = rData.Velocity;
    const auto& vn = rData.Velocity_OldStep1;
    const auto& vnn = rData.Velocity_OldStep2;
    const auto& vmesh = rData.MeshVelocity;
    const auto& vconv = v - vmesh;
    const auto& f = rData.BodyForce;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    auto& rhs = rData.rhs;

    //substitute_rhs_3D

    noalias(rRHS) += rData.Weight * rhs;
}


template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeGaussPointEnrichmentContributions(
	TwoFluidNavierStokesData<2, 3>& rData,
    MatrixType& rV,
    MatrixType& rH,
    MatrixType& rKee,
	VectorType& rRHS_ee) {

	const double rho = rData.Density;
	const double mu = rData.EffectiveViscosity;

	const double h = rData.ElementSize;

	const double dt = rData.DeltaTime;
	const double bdf0 = rData.bdf0;
	const double bdf1 = rData.bdf1;
	const double bdf2 = rData.bdf2;

	const double dyn_tau = rData.DynamicTau;

	const auto& v = rData.Velocity;
	const auto& vn = rData.Velocity_OldStep1;
	const auto& vnn = rData.Velocity_OldStep2;
	const auto& vmesh = rData.MeshVelocity;
	const auto& vconv = v - vmesh;
	const auto& f = rData.BodyForce;
	const auto& p = rData.Pressure;

	// Get shape function values
	const auto& N = rData.N;
	const auto& DN = rData.DN_DX;
	const auto& Nenr = rData.Nenr;
	const auto& DNenr = rData.DN_DXenr;

	// Stabilization parameters
	constexpr double stab_c1 = 4.0;
	constexpr double stab_c2 = 2.0;

    auto& V = rData.V;
    auto& H = rData.H;
    auto& Kee = rData.Kee;
    auto& rhs_ee = rData.rhs_ee;

	array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    //substitute_enrichment_V_2D

    //substitute_enrichment_H_2D

    //substitute_enrichment_Kee_2D

    //substitute_enrichment_rhs_ee_2D

    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;

}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeGaussPointEnrichmentContributions(
	TwoFluidNavierStokesData<3, 4>& rData,
    MatrixType& rV,
    MatrixType& rH,
    MatrixType& rKee,
	VectorType& rRHS_ee) {

	const double rho = rData.Density;
	const double mu = rData.EffectiveViscosity;

	const double h = rData.ElementSize;

	const double dt = rData.DeltaTime;
	const double bdf0 = rData.bdf0;
	const double bdf1 = rData.bdf1;
	const double bdf2 = rData.bdf2;

	const double dyn_tau = rData.DynamicTau;

	const auto& v = rData.Velocity;
	const auto& vn = rData.Velocity_OldStep1;
	const auto& vnn = rData.Velocity_OldStep2;
	const auto& vmesh = rData.MeshVelocity;
	const auto& vconv = v - vmesh;
	const auto& f = rData.BodyForce;
	const auto& p = rData.Pressure;

	// Get shape function values
	const auto& N = rData.N;
	const auto& DN = rData.DN_DX;
	const auto& Nenr = rData.Nenr;
	const auto& DNenr = rData.DN_DXenr;

	// Stabilization parameters
	constexpr double stab_c1 = 4.0;
	constexpr double stab_c2 = 2.0;

    auto& V = rData.V;
    auto& H = rData.H;
    auto& Kee = rData.Kee;
    auto& rhs_ee = rData.rhs_ee;

	array_1d<double, NumNodes> penr = ZeroVector(NumNodes); //penriched is considered to be zero as we do not want to store it

    //substitute_enrichment_V_3D

    //substitute_enrichment_H_3D

    //substitute_enrichment_Kee_3D

    //substitute_enrichment_rhs_ee_3D

    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeSplitting(
	TwoFluidNavierStokesData<3, 4>& rData,
	MatrixType& rShapeFunctionsPos,
    MatrixType& rShapeFunctionsNeg,
    MatrixType& rEnrichedShapeFunctionsPos,
    MatrixType& rEnrichedShapeFunctionsNeg,
	GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
	GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
	GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg)
{
    Matrix Enr_Neg_Interp = ZeroMatrix(NumNodes,NumNodes);
    Matrix Enr_Pos_Interp = ZeroMatrix(NumNodes,NumNodes);
	for (unsigned int i = 0; i < NumNodes; i++)
	{
        if (rData.Distance[i] > 0.0)
            Enr_Neg_Interp(i,i) = 1.0;
        else {
            Enr_Pos_Interp(i,i) = 1.0;
        }
	}
    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, rData.Distance);
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);
        
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, Enr_Pos_Interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, Enr_Neg_Interp);
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++) {
        rEnrichedShapeDerivativesPos[i] = prod(Enr_Pos_Interp, rShapeDerivativesPos[i]);
    }
    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++) {
        rEnrichedShapeDerivativesNeg[i] = prod(Enr_Neg_Interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = p_modified_sh_func->pGetSplittingUtil()->mDivisionsNumber;
}

template <>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeSplitting(
	TwoFluidNavierStokesData<2, 3>& rData,
	MatrixType& rShapeFunctionsPos,
    MatrixType& rShapeFunctionsNeg,
    MatrixType& rEnrichedShapeFunctionsPos,
    MatrixType& rEnrichedShapeFunctionsNeg,
	GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesPos,
	GeometryType::ShapeFunctionsGradientsType& rShapeDerivativesNeg,
    GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesPos,
	GeometryType::ShapeFunctionsGradientsType& rEnrichedShapeDerivativesNeg)
{
        Matrix Enr_Neg_Interp = ZeroMatrix(NumNodes,NumNodes);
    Matrix Enr_Pos_Interp = ZeroMatrix(NumNodes,NumNodes);
	for (unsigned int i = 0; i < NumNodes; i++)
	{
        if (rData.Distance[i] > 0.0)
            Enr_Neg_Interp(i,i) = 1.0;
        else {
            Enr_Pos_Interp(i,i) = 1.0;
        }
	}
    // Construct the modified shape fucntions utility
    GeometryType::Pointer p_geom = this->pGetGeometry();
    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
    p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, rData.Distance);
    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsPos,
        rShapeDerivativesPos,
        rData.w_gauss_pos_side,
        GeometryData::GI_GAUSS_2);

    // Call the negative side modified shape functions calculator
    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rShapeFunctionsNeg,
        rShapeDerivativesNeg,
        rData.w_gauss_neg_side,
        GeometryData::GI_GAUSS_2);
        
    rEnrichedShapeFunctionsPos = prod(rShapeFunctionsPos, Enr_Pos_Interp);
    rEnrichedShapeFunctionsNeg = prod(rShapeFunctionsNeg, Enr_Neg_Interp);
    rEnrichedShapeDerivativesPos = rShapeDerivativesPos;
    rEnrichedShapeDerivativesNeg = rShapeDerivativesNeg;

    for (unsigned int i = 0; i < rShapeDerivativesPos.size(); i++) {
        rEnrichedShapeDerivativesPos[i] = prod(Enr_Pos_Interp, rShapeDerivativesPos[i]);
    }
    for (unsigned int i = 0; i < rShapeDerivativesNeg.size(); i++) {
        rEnrichedShapeDerivativesNeg[i] = prod(Enr_Neg_Interp, rShapeDerivativesNeg[i]);
    }

    rData.NumberOfDivisions = p_modified_sh_func->pGetSplittingUtil()->mDivisionsNumber;
}

template<>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::CondenseEnrichment(
	TwoFluidNavierStokesData<2, 3>& rData,
	Matrix& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
    MatrixType& Htot,
	MatrixType& Vtot,
	MatrixType& Kee_tot,
	VectorType& rhs_ee_tot)
{
    const double min_area_ratio = -1e-6;

	double positive_volume = 0.0;
	double negative_volume = 0.0;
	for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); igauss_pos++)
	{
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
	}

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); igauss_neg++)
	{
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
	}
	const double Vol = positive_volume + negative_volume;



	double max_diag = 0.0;
	for (unsigned int k = 0; k<Dim + 1; k++)
		if (fabs(Kee_tot(k, k)) > max_diag) max_diag = fabs(Kee_tot(k, k));
	if (max_diag == 0) max_diag = 1.0;



	if (positive_volume / Vol < min_area_ratio)
	{
		for (unsigned int i = 0; i<Dim + 1; i++)
		{
			if (rData.Distance[i] >= 0.0)
			{
				Kee_tot(i, i) += 1000.0*max_diag;
			}
		}
	}
	if (negative_volume / Vol < min_area_ratio)
	{
		for (unsigned int i = 0; i<Dim + 1; i++)
		{
			if (rData.Distance[i] < 0.0)
			{
				Kee_tot(i, i) += 1000.0*max_diag;
			}
		}
	}

	//"weakly" impose continuity
	for (unsigned int i = 0; i<Dim; i++)
	{
		const double di = fabs(rData.Distance[i]);

		for (unsigned int j = i + 1; j<Dim + 1; j++)
		{
			const double dj = fabs(rData.Distance[j]);

			if (rData.Distance[i] * rData.Distance[j] < 0.0) //cut edge
			{
				double sum_d = di + dj;
				double Ni = dj / sum_d;
				double Nj = di / sum_d;

				double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
				Kee_tot(i, i) += penalty_coeff * Ni*Ni;
				Kee_tot(i, j) -= penalty_coeff * Ni*Nj;
				Kee_tot(j, i) -= penalty_coeff * Nj*Ni;
				Kee_tot(j, j) += penalty_coeff * Nj*Nj;

			}
		}
	}
	//add to LHS enrichment contributions
	MatrixType inverse_diag;
	inverse_diag.resize(NumNodes, NumNodes, false);
    double det;
    MathUtils<double>::InvertMatrix(Kee_tot, inverse_diag, det);

	BoundedMatrix<double, 4, 16> tmp = prod(inverse_diag, Htot);
	noalias(rLeftHandSideMatrix) -= prod(Vtot, tmp);

	const array_1d<double, 4> tmp2 = prod(inverse_diag, rhs_ee_tot);
	noalias(rRightHandSideVector) -= prod(Vtot, tmp2);

}

template<>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::CondenseEnrichment(
	TwoFluidNavierStokesData<3, 4>& rData,
	Matrix& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	MatrixType& Htot,
	MatrixType& Vtot,
	MatrixType& Kee_tot,
	VectorType& rhs_ee_tot)
{
	const double min_area_ratio = -1e-6;

	double positive_volume = 0.0;
	double negative_volume = 0.0;
	for (unsigned int igauss_pos = 0; igauss_pos < rData.w_gauss_pos_side.size(); igauss_pos++)
	{
        positive_volume += rData.w_gauss_pos_side[igauss_pos];
	}

    for (unsigned int igauss_neg = 0; igauss_neg < rData.w_gauss_neg_side.size(); igauss_neg++)
	{
        negative_volume += rData.w_gauss_neg_side[igauss_neg];
	}
	const double Vol = positive_volume + negative_volume;



	double max_diag = 0.0;
	for (unsigned int k = 0; k<Dim + 1; k++)
		if (fabs(Kee_tot(k, k)) > max_diag) max_diag = fabs(Kee_tot(k, k));
	if (max_diag == 0) max_diag = 1.0;



	if (positive_volume / Vol < min_area_ratio)
	{
		for (unsigned int i = 0; i<Dim + 1; i++)
		{
			if (rData.Distance[i] >= 0.0)
			{
				Kee_tot(i, i) += 1000.0*max_diag;
			}
		}
	}
	if (negative_volume / Vol < min_area_ratio)
	{
		for (unsigned int i = 0; i<Dim + 1; i++)
		{
			if (rData.Distance[i] < 0.0)
			{
				Kee_tot(i, i) += 1000.0*max_diag;
			}
		}
	}

	//"weakly" impose continuity
	for (unsigned int i = 0; i<Dim; i++)
	{
		const double di = fabs(rData.Distance[i]);

		for (unsigned int j = i + 1; j<Dim + 1; j++)
		{
			const double dj = fabs(rData.Distance[j]);

			if (rData.Distance[i] * rData.Distance[j] < 0.0) //cut edge
			{
				double sum_d = di + dj;
				double Ni = dj / sum_d;
				double Nj = di / sum_d;

				double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
				Kee_tot(i, i) += penalty_coeff * Ni*Ni;
				Kee_tot(i, j) -= penalty_coeff * Ni*Nj;
				Kee_tot(j, i) -= penalty_coeff * Nj*Ni;
				Kee_tot(j, j) += penalty_coeff * Nj*Nj;

			}
		}
	}
	//add to LHS enrichment contributions
	MatrixType inverse_diag;
	inverse_diag.resize(NumNodes, NumNodes, false);
    double det;
    MathUtils<double>::InvertMatrix(Kee_tot, inverse_diag, det);

	BoundedMatrix<double, 4, 16> tmp = prod(inverse_diag, Htot);
	noalias(rLeftHandSideMatrix) -= prod(Vtot, tmp);

	const array_1d<double, 4> tmp2 = prod(inverse_diag, rhs_ee_tot);
	noalias(rRightHandSideVector) -= prod(Vtot, tmp2);

}



template< class TElementData >
void TwoFluidNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
	using BaseType = FluidElement<TElementData>;
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template< class TElementData >
void TwoFluidNavierStokes<TElementData>::load(Serializer& rSerializer)
{
	using BaseType = FluidElement<TElementData>;
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> >;
template class TwoFluidNavierStokes< TwoFluidNavierStokesData<3, 4> >;






}
