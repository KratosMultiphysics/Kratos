//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Ruben Zorrilla
//

#include "two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "utilities/enrichment_utilities_duplicate_dofs.h"


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
		// Get Shape function data
		Vector gauss_weights;
		Matrix shape_functions;
		ShapeFunctionDerivativesArrayType shape_derivatives;
		this->CalculateGeometryData(
			gauss_weights, shape_functions, shape_derivatives);
		const unsigned int number_of_gauss_points = gauss_weights.size();

		TElementData data;
		data.Initialize(*this, rCurrentProcessInfo);

		for (unsigned int i = 0; i < NumNodes; i++)
		{
		    if(data.Distance[i] > 0)
				data.NumPositiveNodes++;
		    else
				data.NumNegativeNodes++;
		}

		if (data.IsCut()) {

			std::vector< MatrixType > enriched_shape_derivatives;
			MatrixType enriched_shape_functions;
			data.NumberOfDivisions = ComputeSplitting(data, shape_functions, shape_derivatives, enriched_shape_derivatives, enriched_shape_functions);

			if (data.NumberOfDivisions == 1) {
				//cases exist when the element is like not subdivided due to the characteristics of the provided distance
				//in this cases the element is treated as AIR or FLUID depending on the side
				array_1d<double,NumNodes> Ncenter;
				for(unsigned int i=0; i<NumNodes; i++) Ncenter[i]=0.25;
				const double dgauss = inner_prod(data.Distance, Ncenter);
				for (unsigned int g = 0; g < number_of_gauss_points; g++) {
					data.UpdateGeometryValues(gauss_weights[g],
						row(shape_functions, g),
						shape_derivatives[g]);

					this->CalculateMaterialPropertiesAtGaussPoint(data);

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

				for (unsigned int g = 0; g < data.PartitionsSigns.size(); g++) {
					data.UpdateGeometryValues(data.PartitionsVolumes[g],
						row(shape_functions, g),
						shape_derivatives[0], //is constant
						row(enriched_shape_functions, g),
						enriched_shape_derivatives[g]);
					const double dgauss = inner_prod(data.Distance, data.N);

					this->CalculateMaterialPropertiesAtGaussPoint(data);

					if (dgauss > 0.0)
						data.CalculateAirMaterialResponse();
					else
						this->CalculateMaterialResponse(data);

					this->AddTimeIntegratedSystem(
						data, rLeftHandSideMatrix, rRightHandSideVector);

					ComputeGaussPointEnrichmentContributions(data, Vtot, Htot, Kee_tot, rhs_ee_tot);
				}

				CondenseEnrichment(data, rLeftHandSideMatrix,rRightHandSideVector,Htot,Vtot,Kee_tot, rhs_ee_tot);

			}
		}

		else {
			// Iterate over integration points to evaluate local contribution
			for (unsigned int g = 0; g < number_of_gauss_points; g++) {

				data.UpdateGeometryValues(gauss_weights[g], row(shape_functions, g),
					shape_derivatives[g]);

				this->CalculateMaterialPropertiesAtGaussPoint(data);

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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             C(0,2)*DN(0,0);
const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
const double clhs3 =             pow(DN(0,0), 2);
const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 =             stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             bdf0*rho;
const double clhs9 =             N[0]*rho;
const double clhs10 =             DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs11 =             N[0]*bdf0 + clhs10;
const double clhs12 =             pow(rho, 2);
const double clhs13 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs14 =             1.0/(clhs6*rho/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 =             1.0*N[0]*clhs12*clhs13*clhs14;
const double clhs16 =             1.0*clhs10*clhs12*clhs14;
const double clhs17 =             pow(N[0], 2)*clhs8 + clhs10*clhs9 + clhs11*clhs15 + clhs11*clhs16;
const double clhs18 =             C(0,1)*DN(0,1) + clhs1;
const double clhs19 =             C(1,2)*DN(0,1);
const double clhs20 =             C(2,2)*DN(0,0) + clhs19;
const double clhs21 =             DN(0,0)*clhs7;
const double clhs22 =             DN(0,1)*clhs21;
const double clhs23 =             1.0*clhs13*clhs14;
const double clhs24 =             1.0*clhs14*rho;
const double clhs25 =             clhs10*clhs24;
const double clhs26 =             -N[0] + clhs23*clhs9 + clhs25;
const double clhs27 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs28 =             C(0,2)*DN(1,0);
const double clhs29 =             C(2,2)*DN(1,1) + clhs28;
const double clhs30 =             DN(0,0)*DN(1,0);
const double clhs31 =             clhs30*clhs7;
const double clhs32 =             N[0]*bdf0*rho;
const double clhs33 =             N[1]*clhs32;
const double clhs34 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs35 =             N[1]*bdf0;
const double clhs36 =             clhs34 + clhs35;
const double clhs37 =             clhs15*clhs36 + clhs16*clhs36 + clhs33 + clhs34*clhs9;
const double clhs38 =             C(0,1)*DN(1,1) + clhs28;
const double clhs39 =             C(1,2)*DN(1,1);
const double clhs40 =             C(2,2)*DN(1,0) + clhs39;
const double clhs41 =             DN(1,1)*clhs21;
const double clhs42 =             DN(0,0)*N[1];
const double clhs43 =             DN(1,0)*N[0];
const double clhs44 =             1.0*clhs13*clhs14*rho;
const double clhs45 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs46 =             C(0,2)*DN(2,0);
const double clhs47 =             C(2,2)*DN(2,1) + clhs46;
const double clhs48 =             DN(0,0)*DN(2,0);
const double clhs49 =             clhs48*clhs7;
const double clhs50 =             N[2]*clhs32;
const double clhs51 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs52 =             N[2]*bdf0 + clhs51;
const double clhs53 =             clhs15*clhs52 + clhs16*clhs52 + clhs50 + clhs51*clhs9;
const double clhs54 =             C(0,1)*DN(2,1) + clhs46;
const double clhs55 =             C(1,2)*DN(2,1);
const double clhs56 =             C(2,2)*DN(2,0) + clhs55;
const double clhs57 =             DN(2,1)*clhs21;
const double clhs58 =             DN(0,0)*N[2];
const double clhs59 =             DN(2,0)*N[0];
const double clhs60 =             C(0,1)*DN(0,0) + clhs19;
const double clhs61 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs62 =             pow(DN(0,1), 2);
const double clhs63 =             C(0,1)*DN(1,0) + clhs39;
const double clhs64 =             DN(0,1)*clhs7;
const double clhs65 =             DN(1,0)*clhs64;
const double clhs66 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs67 =             DN(0,1)*DN(1,1);
const double clhs68 =             clhs67*clhs7;
const double clhs69 =             DN(0,1)*N[1];
const double clhs70 =             DN(1,1)*N[0];
const double clhs71 =             C(0,1)*DN(2,0) + clhs55;
const double clhs72 =             DN(2,0)*clhs64;
const double clhs73 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs74 =             DN(0,1)*DN(2,1);
const double clhs75 =             clhs7*clhs74;
const double clhs76 =             DN(0,1)*N[2];
const double clhs77 =             DN(2,1)*N[0];
const double clhs78 =             clhs11*clhs24;
const double clhs79 =             N[0] + clhs78;
const double clhs80 =             1.0*clhs14;
const double clhs81 =             clhs24*clhs36;
const double clhs82 =             clhs80*(clhs30 + clhs67);
const double clhs83 =             clhs24*clhs52;
const double clhs84 =             clhs80*(clhs48 + clhs74);
const double clhs85 =             N[1]*rho;
const double clhs86 =             1.0*N[1]*clhs12*clhs13*clhs14;
const double clhs87 =             1.0*clhs12*clhs14*clhs34;
const double clhs88 =             clhs10*clhs85 + clhs11*clhs86 + clhs11*clhs87 + clhs33;
const double clhs89 =             clhs24*clhs34;
const double clhs90 =             pow(DN(1,0), 2);
const double clhs91 =             pow(N[1], 2)*clhs8 + clhs34*clhs85 + clhs36*clhs86 + clhs36*clhs87;
const double clhs92 =             DN(1,0)*clhs7;
const double clhs93 =             DN(1,1)*clhs92;
const double clhs94 =             -N[1] + clhs23*clhs85 + clhs89;
const double clhs95 =             DN(1,0)*DN(2,0);
const double clhs96 =             clhs7*clhs95;
const double clhs97 =             N[2]*rho;
const double clhs98 =             clhs35*clhs97;
const double clhs99 =             clhs51*clhs85 + clhs52*clhs86 + clhs52*clhs87 + clhs98;
const double clhs100 =             DN(2,1)*clhs92;
const double clhs101 =             DN(1,0)*N[2];
const double clhs102 =             DN(2,0)*N[1];
const double clhs103 =             pow(DN(1,1), 2);
const double clhs104 =             DN(2,0)*clhs7;
const double clhs105 =             DN(1,1)*clhs104;
const double clhs106 =             DN(1,1)*DN(2,1);
const double clhs107 =             clhs106*clhs7;
const double clhs108 =             DN(1,1)*N[2];
const double clhs109 =             DN(2,1)*N[1];
const double clhs110 =             N[1] + clhs81;
const double clhs111 =             clhs80*(clhs106 + clhs95);
const double clhs112 =             1.0*N[2]*clhs12*clhs13*clhs14;
const double clhs113 =             1.0*clhs12*clhs14*clhs51;
const double clhs114 =             clhs10*clhs97 + clhs11*clhs112 + clhs11*clhs113 + clhs50;
const double clhs115 =             clhs24*clhs51;
const double clhs116 =             clhs112*clhs36 + clhs113*clhs36 + clhs34*clhs97 + clhs98;
const double clhs117 =             pow(DN(2,0), 2);
const double clhs118 =             pow(N[2], 2)*clhs8 + clhs112*clhs52 + clhs113*clhs52 + clhs51*clhs97;
const double clhs119 =             DN(2,1)*clhs104;
const double clhs120 =             -N[2] + clhs115 + clhs23*clhs97;
const double clhs121 =             pow(DN(2,1), 2);
const double clhs122 =             N[2] + clhs83;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs17 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs18 + DN(0,1)*clhs20 + clhs22;
            lhs(0,2)=DN(0,0)*clhs26;
            lhs(0,3)=DN(0,0)*clhs27 + DN(0,1)*clhs29 + clhs31 + clhs37;
            lhs(0,4)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + clhs41;
            lhs(0,5)=DN(1,0)*clhs25 - clhs42 + clhs43*clhs44;
            lhs(0,6)=DN(0,0)*clhs45 + DN(0,1)*clhs47 + clhs49 + clhs53;
            lhs(0,7)=DN(0,0)*clhs54 + DN(0,1)*clhs56 + clhs57;
            lhs(0,8)=DN(2,0)*clhs25 + clhs44*clhs59 - clhs58;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs60 + clhs22;
            lhs(1,1)=DN(0,0)*clhs20 + DN(0,1)*clhs61 + clhs17 + clhs62*clhs7;
            lhs(1,2)=DN(0,1)*clhs26;
            lhs(1,3)=DN(0,0)*clhs29 + DN(0,1)*clhs63 + clhs65;
            lhs(1,4)=DN(0,0)*clhs40 + DN(0,1)*clhs66 + clhs37 + clhs68;
            lhs(1,5)=DN(1,1)*clhs25 + clhs44*clhs70 - clhs69;
            lhs(1,6)=DN(0,0)*clhs47 + DN(0,1)*clhs71 + clhs72;
            lhs(1,7)=DN(0,0)*clhs56 + DN(0,1)*clhs73 + clhs53 + clhs75;
            lhs(1,8)=DN(2,1)*clhs25 + clhs44*clhs77 - clhs76;
            lhs(2,0)=DN(0,0)*clhs79;
            lhs(2,1)=DN(0,1)*clhs79;
            lhs(2,2)=clhs80*(clhs3 + clhs62);
            lhs(2,3)=DN(0,0)*clhs81 + clhs43;
            lhs(2,4)=DN(0,1)*clhs81 + clhs70;
            lhs(2,5)=clhs82;
            lhs(2,6)=DN(0,0)*clhs83 + clhs59;
            lhs(2,7)=DN(0,1)*clhs83 + clhs77;
            lhs(2,8)=clhs84;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs31 + clhs88;
            lhs(3,1)=DN(1,0)*clhs18 + DN(1,1)*clhs20 + clhs65;
            lhs(3,2)=DN(0,0)*clhs89 + clhs42*clhs44 - clhs43;
            lhs(3,3)=DN(1,0)*clhs27 + DN(1,1)*clhs29 + clhs7*clhs90 + clhs91;
            lhs(3,4)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + clhs93;
            lhs(3,5)=DN(1,0)*clhs94;
            lhs(3,6)=DN(1,0)*clhs45 + DN(1,1)*clhs47 + clhs96 + clhs99;
            lhs(3,7)=DN(1,0)*clhs54 + DN(1,1)*clhs56 + clhs100;
            lhs(3,8)=DN(2,0)*clhs89 - clhs101 + clhs102*clhs44;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs60 + clhs41;
            lhs(4,1)=DN(1,0)*clhs20 + DN(1,1)*clhs61 + clhs68 + clhs88;
            lhs(4,2)=DN(0,1)*clhs89 + clhs44*clhs69 - clhs70;
            lhs(4,3)=DN(1,0)*clhs29 + DN(1,1)*clhs63 + clhs93;
            lhs(4,4)=DN(1,0)*clhs40 + DN(1,1)*clhs66 + clhs103*clhs7 + clhs91;
            lhs(4,5)=DN(1,1)*clhs94;
            lhs(4,6)=DN(1,0)*clhs47 + DN(1,1)*clhs71 + clhs105;
            lhs(4,7)=DN(1,0)*clhs56 + DN(1,1)*clhs73 + clhs107 + clhs99;
            lhs(4,8)=DN(2,1)*clhs89 - clhs108 + clhs109*clhs44;
            lhs(5,0)=DN(1,0)*clhs78 + clhs42;
            lhs(5,1)=DN(1,1)*clhs78 + clhs69;
            lhs(5,2)=clhs82;
            lhs(5,3)=DN(1,0)*clhs110;
            lhs(5,4)=DN(1,1)*clhs110;
            lhs(5,5)=clhs80*(clhs103 + clhs90);
            lhs(5,6)=DN(1,0)*clhs83 + clhs102;
            lhs(5,7)=DN(1,1)*clhs83 + clhs109;
            lhs(5,8)=clhs111;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs114 + clhs49;
            lhs(6,1)=DN(2,0)*clhs18 + DN(2,1)*clhs20 + clhs72;
            lhs(6,2)=DN(0,0)*clhs115 + clhs44*clhs58 - clhs59;
            lhs(6,3)=DN(2,0)*clhs27 + DN(2,1)*clhs29 + clhs116 + clhs96;
            lhs(6,4)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + clhs105;
            lhs(6,5)=DN(1,0)*clhs115 + clhs101*clhs44 - clhs102;
            lhs(6,6)=DN(2,0)*clhs45 + DN(2,1)*clhs47 + clhs117*clhs7 + clhs118;
            lhs(6,7)=DN(2,0)*clhs54 + DN(2,1)*clhs56 + clhs119;
            lhs(6,8)=DN(2,0)*clhs120;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs60 + clhs57;
            lhs(7,1)=DN(2,0)*clhs20 + DN(2,1)*clhs61 + clhs114 + clhs75;
            lhs(7,2)=DN(0,1)*clhs115 + clhs44*clhs76 - clhs77;
            lhs(7,3)=DN(2,0)*clhs29 + DN(2,1)*clhs63 + clhs100;
            lhs(7,4)=DN(2,0)*clhs40 + DN(2,1)*clhs66 + clhs107 + clhs116;
            lhs(7,5)=DN(1,1)*clhs115 + clhs108*clhs44 - clhs109;
            lhs(7,6)=DN(2,0)*clhs47 + DN(2,1)*clhs71 + clhs119;
            lhs(7,7)=DN(2,0)*clhs56 + DN(2,1)*clhs73 + clhs118 + clhs121*clhs7;
            lhs(7,8)=DN(2,1)*clhs120;
            lhs(8,0)=DN(2,0)*clhs78 + clhs58;
            lhs(8,1)=DN(2,1)*clhs78 + clhs76;
            lhs(8,2)=clhs84;
            lhs(8,3)=DN(2,0)*clhs81 + clhs101;
            lhs(8,4)=DN(2,1)*clhs81 + clhs108;
            lhs(8,5)=clhs111;
            lhs(8,6)=DN(2,0)*clhs122;
            lhs(8,7)=DN(2,1)*clhs122;
            lhs(8,8)=clhs80*(clhs117 + clhs121);


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

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             C(0,3)*DN(0,0);
const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 =             C(0,5)*DN(0,0);
const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 =             stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs9*h/stab_c1 + mu;
const double clhs11 =             bdf0*rho;
const double clhs12 =             N[0]*rho;
const double clhs13 =             DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs14 =             N[0]*bdf0 + clhs13;
const double clhs15 =             pow(rho, 2);
const double clhs16 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs17 =             1.0/(clhs9*rho/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 =             1.0*N[0]*clhs15*clhs16*clhs17;
const double clhs19 =             1.0*clhs13*clhs15*clhs17;
const double clhs20 =             pow(N[0], 2)*clhs11 + clhs12*clhs13 + clhs14*clhs18 + clhs14*clhs19;
const double clhs21 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs22 =             C(1,3)*DN(0,1);
const double clhs23 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs22;
const double clhs24 =             C(3,5)*DN(0,0);
const double clhs25 =             C(4,5)*DN(0,2);
const double clhs26 =             C(1,5)*DN(0,1) + clhs24 + clhs25;
const double clhs27 =             DN(0,0)*clhs10;
const double clhs28 =             DN(0,1)*clhs27;
const double clhs29 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs30 =             C(3,4)*DN(0,1);
const double clhs31 =             C(2,3)*DN(0,2) + clhs24 + clhs30;
const double clhs32 =             C(2,5)*DN(0,2);
const double clhs33 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs32;
const double clhs34 =             DN(0,2)*clhs27;
const double clhs35 =             1.0*clhs16*clhs17;
const double clhs36 =             1.0*clhs17*rho;
const double clhs37 =             clhs13*clhs36;
const double clhs38 =             -N[0] + clhs12*clhs35 + clhs37;
const double clhs39 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs40 =             C(0,3)*DN(1,0);
const double clhs41 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs40;
const double clhs42 =             C(0,5)*DN(1,0);
const double clhs43 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs42;
const double clhs44 =             DN(0,0)*DN(1,0);
const double clhs45 =             clhs10*clhs44;
const double clhs46 =             N[0]*bdf0*rho;
const double clhs47 =             N[1]*clhs46;
const double clhs48 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs49 =             N[1]*bdf0 + clhs48;
const double clhs50 =             clhs12*clhs48 + clhs18*clhs49 + clhs19*clhs49 + clhs47;
const double clhs51 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs40;
const double clhs52 =             C(1,3)*DN(1,1);
const double clhs53 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs52;
const double clhs54 =             C(3,5)*DN(1,0);
const double clhs55 =             C(4,5)*DN(1,2);
const double clhs56 =             C(1,5)*DN(1,1) + clhs54 + clhs55;
const double clhs57 =             DN(1,1)*clhs27;
const double clhs58 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs42;
const double clhs59 =             C(3,4)*DN(1,1);
const double clhs60 =             C(2,3)*DN(1,2) + clhs54 + clhs59;
const double clhs61 =             C(2,5)*DN(1,2);
const double clhs62 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs61;
const double clhs63 =             DN(1,2)*clhs27;
const double clhs64 =             DN(0,0)*N[1];
const double clhs65 =             DN(1,0)*N[0];
const double clhs66 =             1.0*clhs16*clhs17*rho;
const double clhs67 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs68 =             C(0,3)*DN(2,0);
const double clhs69 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs68;
const double clhs70 =             C(0,5)*DN(2,0);
const double clhs71 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs70;
const double clhs72 =             DN(0,0)*DN(2,0);
const double clhs73 =             clhs10*clhs72;
const double clhs74 =             N[2]*clhs46;
const double clhs75 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs76 =             N[2]*bdf0;
const double clhs77 =             clhs75 + clhs76;
const double clhs78 =             clhs12*clhs75 + clhs18*clhs77 + clhs19*clhs77 + clhs74;
const double clhs79 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs68;
const double clhs80 =             C(1,3)*DN(2,1);
const double clhs81 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs80;
const double clhs82 =             C(3,5)*DN(2,0);
const double clhs83 =             C(4,5)*DN(2,2);
const double clhs84 =             C(1,5)*DN(2,1) + clhs82 + clhs83;
const double clhs85 =             DN(2,1)*clhs27;
const double clhs86 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs70;
const double clhs87 =             C(3,4)*DN(2,1);
const double clhs88 =             C(2,3)*DN(2,2) + clhs82 + clhs87;
const double clhs89 =             C(2,5)*DN(2,2);
const double clhs90 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs89;
const double clhs91 =             DN(2,2)*clhs27;
const double clhs92 =             DN(0,0)*N[2];
const double clhs93 =             DN(2,0)*N[0];
const double clhs94 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs95 =             C(0,3)*DN(3,0);
const double clhs96 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs95;
const double clhs97 =             C(0,5)*DN(3,0);
const double clhs98 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs97;
const double clhs99 =             DN(0,0)*DN(3,0);
const double clhs100 =             clhs10*clhs99;
const double clhs101 =             N[3]*clhs46;
const double clhs102 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs103 =             N[3]*bdf0 + clhs102;
const double clhs104 =             clhs101 + clhs102*clhs12 + clhs103*clhs18 + clhs103*clhs19;
const double clhs105 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs95;
const double clhs106 =             C(1,3)*DN(3,1);
const double clhs107 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs106;
const double clhs108 =             C(3,5)*DN(3,0);
const double clhs109 =             C(4,5)*DN(3,2);
const double clhs110 =             C(1,5)*DN(3,1) + clhs108 + clhs109;
const double clhs111 =             DN(3,1)*clhs27;
const double clhs112 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs97;
const double clhs113 =             C(3,4)*DN(3,1);
const double clhs114 =             C(2,3)*DN(3,2) + clhs108 + clhs113;
const double clhs115 =             C(2,5)*DN(3,2);
const double clhs116 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs115;
const double clhs117 =             DN(3,2)*clhs27;
const double clhs118 =             DN(0,0)*N[3];
const double clhs119 =             DN(3,0)*N[0];
const double clhs120 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs22;
const double clhs121 =             C(0,4)*DN(0,0) + clhs25 + clhs30;
const double clhs122 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs123 =             C(1,4)*DN(0,1);
const double clhs124 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs123;
const double clhs125 =             pow(DN(0,1), 2);
const double clhs126 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs123;
const double clhs127 =             C(2,4)*DN(0,2);
const double clhs128 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs127;
const double clhs129 =             DN(0,1)*clhs10;
const double clhs130 =             DN(0,2)*clhs129;
const double clhs131 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs52;
const double clhs132 =             C(0,4)*DN(1,0) + clhs55 + clhs59;
const double clhs133 =             DN(1,0)*clhs129;
const double clhs134 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs135 =             C(1,4)*DN(1,1);
const double clhs136 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs135;
const double clhs137 =             DN(0,1)*DN(1,1);
const double clhs138 =             clhs10*clhs137;
const double clhs139 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs135;
const double clhs140 =             C(2,4)*DN(1,2);
const double clhs141 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs140;
const double clhs142 =             DN(1,2)*clhs129;
const double clhs143 =             DN(0,1)*N[1];
const double clhs144 =             DN(1,1)*N[0];
const double clhs145 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs80;
const double clhs146 =             C(0,4)*DN(2,0) + clhs83 + clhs87;
const double clhs147 =             DN(2,0)*clhs129;
const double clhs148 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs149 =             C(1,4)*DN(2,1);
const double clhs150 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs149;
const double clhs151 =             DN(0,1)*DN(2,1);
const double clhs152 =             clhs10*clhs151;
const double clhs153 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs149;
const double clhs154 =             C(2,4)*DN(2,2);
const double clhs155 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs154;
const double clhs156 =             DN(2,2)*clhs129;
const double clhs157 =             DN(0,1)*N[2];
const double clhs158 =             DN(2,1)*N[0];
const double clhs159 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs106;
const double clhs160 =             C(0,4)*DN(3,0) + clhs109 + clhs113;
const double clhs161 =             DN(3,0)*clhs129;
const double clhs162 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs163 =             C(1,4)*DN(3,1);
const double clhs164 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs163;
const double clhs165 =             DN(0,1)*DN(3,1);
const double clhs166 =             clhs10*clhs165;
const double clhs167 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs163;
const double clhs168 =             C(2,4)*DN(3,2);
const double clhs169 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs168;
const double clhs170 =             DN(3,2)*clhs129;
const double clhs171 =             DN(0,1)*N[3];
const double clhs172 =             DN(3,1)*N[0];
const double clhs173 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs32;
const double clhs174 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs127;
const double clhs175 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs176 =             pow(DN(0,2), 2);
const double clhs177 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs61;
const double clhs178 =             DN(0,2)*clhs10;
const double clhs179 =             DN(1,0)*clhs178;
const double clhs180 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs140;
const double clhs181 =             DN(1,1)*clhs178;
const double clhs182 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs183 =             DN(0,2)*DN(1,2);
const double clhs184 =             clhs10*clhs183;
const double clhs185 =             DN(0,2)*N[1];
const double clhs186 =             DN(1,2)*N[0];
const double clhs187 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs89;
const double clhs188 =             DN(2,0)*clhs178;
const double clhs189 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs154;
const double clhs190 =             DN(2,1)*clhs178;
const double clhs191 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs192 =             DN(0,2)*DN(2,2);
const double clhs193 =             clhs10*clhs192;
const double clhs194 =             DN(0,2)*N[2];
const double clhs195 =             DN(2,2)*N[0];
const double clhs196 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs115;
const double clhs197 =             DN(3,0)*clhs178;
const double clhs198 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs168;
const double clhs199 =             DN(3,1)*clhs178;
const double clhs200 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs201 =             DN(0,2)*DN(3,2);
const double clhs202 =             clhs10*clhs201;
const double clhs203 =             DN(0,2)*N[3];
const double clhs204 =             DN(3,2)*N[0];
const double clhs205 =             clhs14*clhs36;
const double clhs206 =             N[0] + clhs205;
const double clhs207 =             1.0*clhs17;
const double clhs208 =             clhs36*clhs49;
const double clhs209 =             clhs207*(clhs137 + clhs183 + clhs44);
const double clhs210 =             clhs36*clhs77;
const double clhs211 =             clhs207*(clhs151 + clhs192 + clhs72);
const double clhs212 =             clhs103*clhs36;
const double clhs213 =             clhs207*(clhs165 + clhs201 + clhs99);
const double clhs214 =             N[1]*rho;
const double clhs215 =             1.0*N[1]*clhs15*clhs16*clhs17;
const double clhs216 =             1.0*clhs15*clhs17*clhs48;
const double clhs217 =             clhs13*clhs214 + clhs14*clhs215 + clhs14*clhs216 + clhs47;
const double clhs218 =             clhs36*clhs48;
const double clhs219 =             pow(DN(1,0), 2);
const double clhs220 =             pow(N[1], 2)*clhs11 + clhs214*clhs48 + clhs215*clhs49 + clhs216*clhs49;
const double clhs221 =             DN(1,0)*clhs10;
const double clhs222 =             DN(1,1)*clhs221;
const double clhs223 =             DN(1,2)*clhs221;
const double clhs224 =             -N[1] + clhs214*clhs35 + clhs218;
const double clhs225 =             DN(1,0)*DN(2,0);
const double clhs226 =             clhs10*clhs225;
const double clhs227 =             N[1]*bdf0*rho;
const double clhs228 =             N[2]*clhs227;
const double clhs229 =             clhs214*clhs75 + clhs215*clhs77 + clhs216*clhs77 + clhs228;
const double clhs230 =             DN(2,1)*clhs221;
const double clhs231 =             DN(2,2)*clhs221;
const double clhs232 =             DN(1,0)*N[2];
const double clhs233 =             DN(2,0)*N[1];
const double clhs234 =             DN(1,0)*DN(3,0);
const double clhs235 =             clhs10*clhs234;
const double clhs236 =             N[3]*clhs227;
const double clhs237 =             clhs102*clhs214 + clhs103*clhs215 + clhs103*clhs216 + clhs236;
const double clhs238 =             DN(3,1)*clhs221;
const double clhs239 =             DN(3,2)*clhs221;
const double clhs240 =             DN(1,0)*N[3];
const double clhs241 =             DN(3,0)*N[1];
const double clhs242 =             pow(DN(1,1), 2);
const double clhs243 =             DN(1,1)*clhs10;
const double clhs244 =             DN(1,2)*clhs243;
const double clhs245 =             DN(2,0)*clhs243;
const double clhs246 =             DN(1,1)*DN(2,1);
const double clhs247 =             clhs10*clhs246;
const double clhs248 =             DN(2,2)*clhs243;
const double clhs249 =             DN(1,1)*N[2];
const double clhs250 =             DN(2,1)*N[1];
const double clhs251 =             DN(3,0)*clhs243;
const double clhs252 =             DN(1,1)*DN(3,1);
const double clhs253 =             clhs10*clhs252;
const double clhs254 =             DN(3,2)*clhs243;
const double clhs255 =             DN(1,1)*N[3];
const double clhs256 =             DN(3,1)*N[1];
const double clhs257 =             pow(DN(1,2), 2);
const double clhs258 =             DN(1,2)*clhs10;
const double clhs259 =             DN(2,0)*clhs258;
const double clhs260 =             DN(2,1)*clhs258;
const double clhs261 =             DN(1,2)*DN(2,2);
const double clhs262 =             clhs10*clhs261;
const double clhs263 =             DN(1,2)*N[2];
const double clhs264 =             DN(2,2)*N[1];
const double clhs265 =             DN(3,0)*clhs258;
const double clhs266 =             DN(3,1)*clhs258;
const double clhs267 =             DN(1,2)*DN(3,2);
const double clhs268 =             clhs10*clhs267;
const double clhs269 =             DN(1,2)*N[3];
const double clhs270 =             DN(3,2)*N[1];
const double clhs271 =             N[1] + clhs208;
const double clhs272 =             clhs207*(clhs225 + clhs246 + clhs261);
const double clhs273 =             clhs207*(clhs234 + clhs252 + clhs267);
const double clhs274 =             N[2]*rho;
const double clhs275 =             1.0*N[2]*clhs15*clhs16*clhs17;
const double clhs276 =             1.0*clhs15*clhs17*clhs75;
const double clhs277 =             clhs13*clhs274 + clhs14*clhs275 + clhs14*clhs276 + clhs74;
const double clhs278 =             clhs36*clhs75;
const double clhs279 =             clhs228 + clhs274*clhs48 + clhs275*clhs49 + clhs276*clhs49;
const double clhs280 =             pow(DN(2,0), 2);
const double clhs281 =             pow(N[2], 2)*clhs11 + clhs274*clhs75 + clhs275*clhs77 + clhs276*clhs77;
const double clhs282 =             DN(2,0)*clhs10;
const double clhs283 =             DN(2,1)*clhs282;
const double clhs284 =             DN(2,2)*clhs282;
const double clhs285 =             -N[2] + clhs274*clhs35 + clhs278;
const double clhs286 =             DN(2,0)*DN(3,0);
const double clhs287 =             clhs10*clhs286;
const double clhs288 =             N[3]*rho;
const double clhs289 =             clhs288*clhs76;
const double clhs290 =             clhs102*clhs274 + clhs103*clhs275 + clhs103*clhs276 + clhs289;
const double clhs291 =             DN(3,1)*clhs282;
const double clhs292 =             DN(3,2)*clhs282;
const double clhs293 =             DN(2,0)*N[3];
const double clhs294 =             DN(3,0)*N[2];
const double clhs295 =             pow(DN(2,1), 2);
const double clhs296 =             DN(2,1)*clhs10;
const double clhs297 =             DN(2,2)*clhs296;
const double clhs298 =             DN(3,0)*clhs296;
const double clhs299 =             DN(2,1)*DN(3,1);
const double clhs300 =             clhs10*clhs299;
const double clhs301 =             DN(3,2)*clhs296;
const double clhs302 =             DN(2,1)*N[3];
const double clhs303 =             DN(3,1)*N[2];
const double clhs304 =             pow(DN(2,2), 2);
const double clhs305 =             DN(2,2)*clhs10;
const double clhs306 =             DN(3,0)*clhs305;
const double clhs307 =             DN(3,1)*clhs305;
const double clhs308 =             DN(2,2)*DN(3,2);
const double clhs309 =             clhs10*clhs308;
const double clhs310 =             DN(2,2)*N[3];
const double clhs311 =             DN(3,2)*N[2];
const double clhs312 =             N[2] + clhs210;
const double clhs313 =             clhs207*(clhs286 + clhs299 + clhs308);
const double clhs314 =             1.0*N[3]*clhs15*clhs16*clhs17;
const double clhs315 =             1.0*clhs102*clhs15*clhs17;
const double clhs316 =             clhs101 + clhs13*clhs288 + clhs14*clhs314 + clhs14*clhs315;
const double clhs317 =             clhs102*clhs36;
const double clhs318 =             clhs236 + clhs288*clhs48 + clhs314*clhs49 + clhs315*clhs49;
const double clhs319 =             clhs288*clhs75 + clhs289 + clhs314*clhs77 + clhs315*clhs77;
const double clhs320 =             pow(DN(3,0), 2);
const double clhs321 =             pow(N[3], 2)*clhs11 + clhs102*clhs288 + clhs103*clhs314 + clhs103*clhs315;
const double clhs322 =             DN(3,0)*clhs10;
const double clhs323 =             DN(3,1)*clhs322;
const double clhs324 =             DN(3,2)*clhs322;
const double clhs325 =             -N[3] + clhs288*clhs35 + clhs317;
const double clhs326 =             pow(DN(3,1), 2);
const double clhs327 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs328 =             pow(DN(3,2), 2);
const double clhs329 =             N[3] + clhs212;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs20;
            lhs(0,1)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + DN(0,2)*clhs26 + clhs28;
            lhs(0,2)=DN(0,0)*clhs29 + DN(0,1)*clhs31 + DN(0,2)*clhs33 + clhs34;
            lhs(0,3)=DN(0,0)*clhs38;
            lhs(0,4)=DN(0,0)*clhs39 + DN(0,1)*clhs41 + DN(0,2)*clhs43 + clhs45 + clhs50;
            lhs(0,5)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + DN(0,2)*clhs56 + clhs57;
            lhs(0,6)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs62 + clhs63;
            lhs(0,7)=DN(1,0)*clhs37 - clhs64 + clhs65*clhs66;
            lhs(0,8)=DN(0,0)*clhs67 + DN(0,1)*clhs69 + DN(0,2)*clhs71 + clhs73 + clhs78;
            lhs(0,9)=DN(0,0)*clhs79 + DN(0,1)*clhs81 + DN(0,2)*clhs84 + clhs85;
            lhs(0,10)=DN(0,0)*clhs86 + DN(0,1)*clhs88 + DN(0,2)*clhs90 + clhs91;
            lhs(0,11)=DN(2,0)*clhs37 + clhs66*clhs93 - clhs92;
            lhs(0,12)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs100 + clhs104;
            lhs(0,13)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs110 + clhs111;
            lhs(0,14)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs117;
            lhs(0,15)=DN(3,0)*clhs37 - clhs118 + clhs119*clhs66;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs120 + DN(0,2)*clhs121 + clhs28;
            lhs(1,1)=DN(0,0)*clhs23 + DN(0,1)*clhs122 + DN(0,2)*clhs124 + clhs10*clhs125 + clhs20;
            lhs(1,2)=DN(0,0)*clhs31 + DN(0,1)*clhs126 + DN(0,2)*clhs128 + clhs130;
            lhs(1,3)=DN(0,1)*clhs38;
            lhs(1,4)=DN(0,0)*clhs41 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs133;
            lhs(1,5)=DN(0,0)*clhs53 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs138 + clhs50;
            lhs(1,6)=DN(0,0)*clhs60 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142;
            lhs(1,7)=DN(1,1)*clhs37 - clhs143 + clhs144*clhs66;
            lhs(1,8)=DN(0,0)*clhs69 + DN(0,1)*clhs145 + DN(0,2)*clhs146 + clhs147;
            lhs(1,9)=DN(0,0)*clhs81 + DN(0,1)*clhs148 + DN(0,2)*clhs150 + clhs152 + clhs78;
            lhs(1,10)=DN(0,0)*clhs88 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs156;
            lhs(1,11)=DN(2,1)*clhs37 - clhs157 + clhs158*clhs66;
            lhs(1,12)=DN(0,0)*clhs96 + DN(0,1)*clhs159 + DN(0,2)*clhs160 + clhs161;
            lhs(1,13)=DN(0,0)*clhs107 + DN(0,1)*clhs162 + DN(0,2)*clhs164 + clhs104 + clhs166;
            lhs(1,14)=DN(0,0)*clhs114 + DN(0,1)*clhs167 + DN(0,2)*clhs169 + clhs170;
            lhs(1,15)=DN(3,1)*clhs37 - clhs171 + clhs172*clhs66;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs121 + DN(0,2)*clhs173 + clhs34;
            lhs(2,1)=DN(0,0)*clhs26 + DN(0,1)*clhs124 + DN(0,2)*clhs174 + clhs130;
            lhs(2,2)=DN(0,0)*clhs33 + DN(0,1)*clhs128 + DN(0,2)*clhs175 + clhs10*clhs176 + clhs20;
            lhs(2,3)=DN(0,2)*clhs38;
            lhs(2,4)=DN(0,0)*clhs43 + DN(0,1)*clhs132 + DN(0,2)*clhs177 + clhs179;
            lhs(2,5)=DN(0,0)*clhs56 + DN(0,1)*clhs136 + DN(0,2)*clhs180 + clhs181;
            lhs(2,6)=DN(0,0)*clhs62 + DN(0,1)*clhs141 + DN(0,2)*clhs182 + clhs184 + clhs50;
            lhs(2,7)=DN(1,2)*clhs37 - clhs185 + clhs186*clhs66;
            lhs(2,8)=DN(0,0)*clhs71 + DN(0,1)*clhs146 + DN(0,2)*clhs187 + clhs188;
            lhs(2,9)=DN(0,0)*clhs84 + DN(0,1)*clhs150 + DN(0,2)*clhs189 + clhs190;
            lhs(2,10)=DN(0,0)*clhs90 + DN(0,1)*clhs155 + DN(0,2)*clhs191 + clhs193 + clhs78;
            lhs(2,11)=DN(2,2)*clhs37 - clhs194 + clhs195*clhs66;
            lhs(2,12)=DN(0,0)*clhs98 + DN(0,1)*clhs160 + DN(0,2)*clhs196 + clhs197;
            lhs(2,13)=DN(0,0)*clhs110 + DN(0,1)*clhs164 + DN(0,2)*clhs198 + clhs199;
            lhs(2,14)=DN(0,0)*clhs116 + DN(0,1)*clhs169 + DN(0,2)*clhs200 + clhs104 + clhs202;
            lhs(2,15)=DN(3,2)*clhs37 - clhs203 + clhs204*clhs66;
            lhs(3,0)=DN(0,0)*clhs206;
            lhs(3,1)=DN(0,1)*clhs206;
            lhs(3,2)=DN(0,2)*clhs206;
            lhs(3,3)=clhs207*(clhs125 + clhs176 + clhs5);
            lhs(3,4)=DN(0,0)*clhs208 + clhs65;
            lhs(3,5)=DN(0,1)*clhs208 + clhs144;
            lhs(3,6)=DN(0,2)*clhs208 + clhs186;
            lhs(3,7)=clhs209;
            lhs(3,8)=DN(0,0)*clhs210 + clhs93;
            lhs(3,9)=DN(0,1)*clhs210 + clhs158;
            lhs(3,10)=DN(0,2)*clhs210 + clhs195;
            lhs(3,11)=clhs211;
            lhs(3,12)=DN(0,0)*clhs212 + clhs119;
            lhs(3,13)=DN(0,1)*clhs212 + clhs172;
            lhs(3,14)=DN(0,2)*clhs212 + clhs204;
            lhs(3,15)=clhs213;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs217 + clhs45;
            lhs(4,1)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + DN(1,2)*clhs26 + clhs133;
            lhs(4,2)=DN(1,0)*clhs29 + DN(1,1)*clhs31 + DN(1,2)*clhs33 + clhs179;
            lhs(4,3)=DN(0,0)*clhs218 + clhs64*clhs66 - clhs65;
            lhs(4,4)=DN(1,0)*clhs39 + DN(1,1)*clhs41 + DN(1,2)*clhs43 + clhs10*clhs219 + clhs220;
            lhs(4,5)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + DN(1,2)*clhs56 + clhs222;
            lhs(4,6)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs62 + clhs223;
            lhs(4,7)=DN(1,0)*clhs224;
            lhs(4,8)=DN(1,0)*clhs67 + DN(1,1)*clhs69 + DN(1,2)*clhs71 + clhs226 + clhs229;
            lhs(4,9)=DN(1,0)*clhs79 + DN(1,1)*clhs81 + DN(1,2)*clhs84 + clhs230;
            lhs(4,10)=DN(1,0)*clhs86 + DN(1,1)*clhs88 + DN(1,2)*clhs90 + clhs231;
            lhs(4,11)=DN(2,0)*clhs218 - clhs232 + clhs233*clhs66;
            lhs(4,12)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs235 + clhs237;
            lhs(4,13)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs110 + clhs238;
            lhs(4,14)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs239;
            lhs(4,15)=DN(3,0)*clhs218 - clhs240 + clhs241*clhs66;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs120 + DN(1,2)*clhs121 + clhs57;
            lhs(5,1)=DN(1,0)*clhs23 + DN(1,1)*clhs122 + DN(1,2)*clhs124 + clhs138 + clhs217;
            lhs(5,2)=DN(1,0)*clhs31 + DN(1,1)*clhs126 + DN(1,2)*clhs128 + clhs181;
            lhs(5,3)=DN(0,1)*clhs218 + clhs143*clhs66 - clhs144;
            lhs(5,4)=DN(1,0)*clhs41 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs222;
            lhs(5,5)=DN(1,0)*clhs53 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs10*clhs242 + clhs220;
            lhs(5,6)=DN(1,0)*clhs60 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs244;
            lhs(5,7)=DN(1,1)*clhs224;
            lhs(5,8)=DN(1,0)*clhs69 + DN(1,1)*clhs145 + DN(1,2)*clhs146 + clhs245;
            lhs(5,9)=DN(1,0)*clhs81 + DN(1,1)*clhs148 + DN(1,2)*clhs150 + clhs229 + clhs247;
            lhs(5,10)=DN(1,0)*clhs88 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs248;
            lhs(5,11)=DN(2,1)*clhs218 - clhs249 + clhs250*clhs66;
            lhs(5,12)=DN(1,0)*clhs96 + DN(1,1)*clhs159 + DN(1,2)*clhs160 + clhs251;
            lhs(5,13)=DN(1,0)*clhs107 + DN(1,1)*clhs162 + DN(1,2)*clhs164 + clhs237 + clhs253;
            lhs(5,14)=DN(1,0)*clhs114 + DN(1,1)*clhs167 + DN(1,2)*clhs169 + clhs254;
            lhs(5,15)=DN(3,1)*clhs218 - clhs255 + clhs256*clhs66;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs121 + DN(1,2)*clhs173 + clhs63;
            lhs(6,1)=DN(1,0)*clhs26 + DN(1,1)*clhs124 + DN(1,2)*clhs174 + clhs142;
            lhs(6,2)=DN(1,0)*clhs33 + DN(1,1)*clhs128 + DN(1,2)*clhs175 + clhs184 + clhs217;
            lhs(6,3)=DN(0,2)*clhs218 + clhs185*clhs66 - clhs186;
            lhs(6,4)=DN(1,0)*clhs43 + DN(1,1)*clhs132 + DN(1,2)*clhs177 + clhs223;
            lhs(6,5)=DN(1,0)*clhs56 + DN(1,1)*clhs136 + DN(1,2)*clhs180 + clhs244;
            lhs(6,6)=DN(1,0)*clhs62 + DN(1,1)*clhs141 + DN(1,2)*clhs182 + clhs10*clhs257 + clhs220;
            lhs(6,7)=DN(1,2)*clhs224;
            lhs(6,8)=DN(1,0)*clhs71 + DN(1,1)*clhs146 + DN(1,2)*clhs187 + clhs259;
            lhs(6,9)=DN(1,0)*clhs84 + DN(1,1)*clhs150 + DN(1,2)*clhs189 + clhs260;
            lhs(6,10)=DN(1,0)*clhs90 + DN(1,1)*clhs155 + DN(1,2)*clhs191 + clhs229 + clhs262;
            lhs(6,11)=DN(2,2)*clhs218 - clhs263 + clhs264*clhs66;
            lhs(6,12)=DN(1,0)*clhs98 + DN(1,1)*clhs160 + DN(1,2)*clhs196 + clhs265;
            lhs(6,13)=DN(1,0)*clhs110 + DN(1,1)*clhs164 + DN(1,2)*clhs198 + clhs266;
            lhs(6,14)=DN(1,0)*clhs116 + DN(1,1)*clhs169 + DN(1,2)*clhs200 + clhs237 + clhs268;
            lhs(6,15)=DN(3,2)*clhs218 - clhs269 + clhs270*clhs66;
            lhs(7,0)=DN(1,0)*clhs205 + clhs64;
            lhs(7,1)=DN(1,1)*clhs205 + clhs143;
            lhs(7,2)=DN(1,2)*clhs205 + clhs185;
            lhs(7,3)=clhs209;
            lhs(7,4)=DN(1,0)*clhs271;
            lhs(7,5)=DN(1,1)*clhs271;
            lhs(7,6)=DN(1,2)*clhs271;
            lhs(7,7)=clhs207*(clhs219 + clhs242 + clhs257);
            lhs(7,8)=DN(1,0)*clhs210 + clhs233;
            lhs(7,9)=DN(1,1)*clhs210 + clhs250;
            lhs(7,10)=DN(1,2)*clhs210 + clhs264;
            lhs(7,11)=clhs272;
            lhs(7,12)=DN(1,0)*clhs212 + clhs241;
            lhs(7,13)=DN(1,1)*clhs212 + clhs256;
            lhs(7,14)=DN(1,2)*clhs212 + clhs270;
            lhs(7,15)=clhs273;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs277 + clhs73;
            lhs(8,1)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + DN(2,2)*clhs26 + clhs147;
            lhs(8,2)=DN(2,0)*clhs29 + DN(2,1)*clhs31 + DN(2,2)*clhs33 + clhs188;
            lhs(8,3)=DN(0,0)*clhs278 + clhs66*clhs92 - clhs93;
            lhs(8,4)=DN(2,0)*clhs39 + DN(2,1)*clhs41 + DN(2,2)*clhs43 + clhs226 + clhs279;
            lhs(8,5)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + DN(2,2)*clhs56 + clhs245;
            lhs(8,6)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs62 + clhs259;
            lhs(8,7)=DN(1,0)*clhs278 + clhs232*clhs66 - clhs233;
            lhs(8,8)=DN(2,0)*clhs67 + DN(2,1)*clhs69 + DN(2,2)*clhs71 + clhs10*clhs280 + clhs281;
            lhs(8,9)=DN(2,0)*clhs79 + DN(2,1)*clhs81 + DN(2,2)*clhs84 + clhs283;
            lhs(8,10)=DN(2,0)*clhs86 + DN(2,1)*clhs88 + DN(2,2)*clhs90 + clhs284;
            lhs(8,11)=DN(2,0)*clhs285;
            lhs(8,12)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs287 + clhs290;
            lhs(8,13)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs110 + clhs291;
            lhs(8,14)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs292;
            lhs(8,15)=DN(3,0)*clhs278 - clhs293 + clhs294*clhs66;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs120 + DN(2,2)*clhs121 + clhs85;
            lhs(9,1)=DN(2,0)*clhs23 + DN(2,1)*clhs122 + DN(2,2)*clhs124 + clhs152 + clhs277;
            lhs(9,2)=DN(2,0)*clhs31 + DN(2,1)*clhs126 + DN(2,2)*clhs128 + clhs190;
            lhs(9,3)=DN(0,1)*clhs278 + clhs157*clhs66 - clhs158;
            lhs(9,4)=DN(2,0)*clhs41 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs230;
            lhs(9,5)=DN(2,0)*clhs53 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs247 + clhs279;
            lhs(9,6)=DN(2,0)*clhs60 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs260;
            lhs(9,7)=DN(1,1)*clhs278 + clhs249*clhs66 - clhs250;
            lhs(9,8)=DN(2,0)*clhs69 + DN(2,1)*clhs145 + DN(2,2)*clhs146 + clhs283;
            lhs(9,9)=DN(2,0)*clhs81 + DN(2,1)*clhs148 + DN(2,2)*clhs150 + clhs10*clhs295 + clhs281;
            lhs(9,10)=DN(2,0)*clhs88 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs297;
            lhs(9,11)=DN(2,1)*clhs285;
            lhs(9,12)=DN(2,0)*clhs96 + DN(2,1)*clhs159 + DN(2,2)*clhs160 + clhs298;
            lhs(9,13)=DN(2,0)*clhs107 + DN(2,1)*clhs162 + DN(2,2)*clhs164 + clhs290 + clhs300;
            lhs(9,14)=DN(2,0)*clhs114 + DN(2,1)*clhs167 + DN(2,2)*clhs169 + clhs301;
            lhs(9,15)=DN(3,1)*clhs278 - clhs302 + clhs303*clhs66;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs121 + DN(2,2)*clhs173 + clhs91;
            lhs(10,1)=DN(2,0)*clhs26 + DN(2,1)*clhs124 + DN(2,2)*clhs174 + clhs156;
            lhs(10,2)=DN(2,0)*clhs33 + DN(2,1)*clhs128 + DN(2,2)*clhs175 + clhs193 + clhs277;
            lhs(10,3)=DN(0,2)*clhs278 + clhs194*clhs66 - clhs195;
            lhs(10,4)=DN(2,0)*clhs43 + DN(2,1)*clhs132 + DN(2,2)*clhs177 + clhs231;
            lhs(10,5)=DN(2,0)*clhs56 + DN(2,1)*clhs136 + DN(2,2)*clhs180 + clhs248;
            lhs(10,6)=DN(2,0)*clhs62 + DN(2,1)*clhs141 + DN(2,2)*clhs182 + clhs262 + clhs279;
            lhs(10,7)=DN(1,2)*clhs278 + clhs263*clhs66 - clhs264;
            lhs(10,8)=DN(2,0)*clhs71 + DN(2,1)*clhs146 + DN(2,2)*clhs187 + clhs284;
            lhs(10,9)=DN(2,0)*clhs84 + DN(2,1)*clhs150 + DN(2,2)*clhs189 + clhs297;
            lhs(10,10)=DN(2,0)*clhs90 + DN(2,1)*clhs155 + DN(2,2)*clhs191 + clhs10*clhs304 + clhs281;
            lhs(10,11)=DN(2,2)*clhs285;
            lhs(10,12)=DN(2,0)*clhs98 + DN(2,1)*clhs160 + DN(2,2)*clhs196 + clhs306;
            lhs(10,13)=DN(2,0)*clhs110 + DN(2,1)*clhs164 + DN(2,2)*clhs198 + clhs307;
            lhs(10,14)=DN(2,0)*clhs116 + DN(2,1)*clhs169 + DN(2,2)*clhs200 + clhs290 + clhs309;
            lhs(10,15)=DN(3,2)*clhs278 - clhs310 + clhs311*clhs66;
            lhs(11,0)=DN(2,0)*clhs205 + clhs92;
            lhs(11,1)=DN(2,1)*clhs205 + clhs157;
            lhs(11,2)=DN(2,2)*clhs205 + clhs194;
            lhs(11,3)=clhs211;
            lhs(11,4)=DN(2,0)*clhs208 + clhs232;
            lhs(11,5)=DN(2,1)*clhs208 + clhs249;
            lhs(11,6)=DN(2,2)*clhs208 + clhs263;
            lhs(11,7)=clhs272;
            lhs(11,8)=DN(2,0)*clhs312;
            lhs(11,9)=DN(2,1)*clhs312;
            lhs(11,10)=DN(2,2)*clhs312;
            lhs(11,11)=clhs207*(clhs280 + clhs295 + clhs304);
            lhs(11,12)=DN(2,0)*clhs212 + clhs294;
            lhs(11,13)=DN(2,1)*clhs212 + clhs303;
            lhs(11,14)=DN(2,2)*clhs212 + clhs311;
            lhs(11,15)=clhs313;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs100 + clhs316;
            lhs(12,1)=DN(3,0)*clhs21 + DN(3,1)*clhs23 + DN(3,2)*clhs26 + clhs161;
            lhs(12,2)=DN(3,0)*clhs29 + DN(3,1)*clhs31 + DN(3,2)*clhs33 + clhs197;
            lhs(12,3)=DN(0,0)*clhs317 + clhs118*clhs66 - clhs119;
            lhs(12,4)=DN(3,0)*clhs39 + DN(3,1)*clhs41 + DN(3,2)*clhs43 + clhs235 + clhs318;
            lhs(12,5)=DN(3,0)*clhs51 + DN(3,1)*clhs53 + DN(3,2)*clhs56 + clhs251;
            lhs(12,6)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs62 + clhs265;
            lhs(12,7)=DN(1,0)*clhs317 + clhs240*clhs66 - clhs241;
            lhs(12,8)=DN(3,0)*clhs67 + DN(3,1)*clhs69 + DN(3,2)*clhs71 + clhs287 + clhs319;
            lhs(12,9)=DN(3,0)*clhs79 + DN(3,1)*clhs81 + DN(3,2)*clhs84 + clhs298;
            lhs(12,10)=DN(3,0)*clhs86 + DN(3,1)*clhs88 + DN(3,2)*clhs90 + clhs306;
            lhs(12,11)=DN(2,0)*clhs317 + clhs293*clhs66 - clhs294;
            lhs(12,12)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs10*clhs320 + clhs321;
            lhs(12,13)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs110 + clhs323;
            lhs(12,14)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs324;
            lhs(12,15)=DN(3,0)*clhs325;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs120 + DN(3,2)*clhs121 + clhs111;
            lhs(13,1)=DN(3,0)*clhs23 + DN(3,1)*clhs122 + DN(3,2)*clhs124 + clhs166 + clhs316;
            lhs(13,2)=DN(3,0)*clhs31 + DN(3,1)*clhs126 + DN(3,2)*clhs128 + clhs199;
            lhs(13,3)=DN(0,1)*clhs317 + clhs171*clhs66 - clhs172;
            lhs(13,4)=DN(3,0)*clhs41 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs238;
            lhs(13,5)=DN(3,0)*clhs53 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs253 + clhs318;
            lhs(13,6)=DN(3,0)*clhs60 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs266;
            lhs(13,7)=DN(1,1)*clhs317 + clhs255*clhs66 - clhs256;
            lhs(13,8)=DN(3,0)*clhs69 + DN(3,1)*clhs145 + DN(3,2)*clhs146 + clhs291;
            lhs(13,9)=DN(3,0)*clhs81 + DN(3,1)*clhs148 + DN(3,2)*clhs150 + clhs300 + clhs319;
            lhs(13,10)=DN(3,0)*clhs88 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs307;
            lhs(13,11)=DN(2,1)*clhs317 + clhs302*clhs66 - clhs303;
            lhs(13,12)=DN(3,0)*clhs96 + DN(3,1)*clhs159 + DN(3,2)*clhs160 + clhs323;
            lhs(13,13)=DN(3,0)*clhs107 + DN(3,1)*clhs162 + DN(3,2)*clhs164 + clhs10*clhs326 + clhs321;
            lhs(13,14)=DN(3,0)*clhs114 + DN(3,1)*clhs167 + DN(3,2)*clhs169 + clhs327;
            lhs(13,15)=DN(3,1)*clhs325;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs121 + DN(3,2)*clhs173 + clhs117;
            lhs(14,1)=DN(3,0)*clhs26 + DN(3,1)*clhs124 + DN(3,2)*clhs174 + clhs170;
            lhs(14,2)=DN(3,0)*clhs33 + DN(3,1)*clhs128 + DN(3,2)*clhs175 + clhs202 + clhs316;
            lhs(14,3)=DN(0,2)*clhs317 + clhs203*clhs66 - clhs204;
            lhs(14,4)=DN(3,0)*clhs43 + DN(3,1)*clhs132 + DN(3,2)*clhs177 + clhs239;
            lhs(14,5)=DN(3,0)*clhs56 + DN(3,1)*clhs136 + DN(3,2)*clhs180 + clhs254;
            lhs(14,6)=DN(3,0)*clhs62 + DN(3,1)*clhs141 + DN(3,2)*clhs182 + clhs268 + clhs318;
            lhs(14,7)=DN(1,2)*clhs317 + clhs269*clhs66 - clhs270;
            lhs(14,8)=DN(3,0)*clhs71 + DN(3,1)*clhs146 + DN(3,2)*clhs187 + clhs292;
            lhs(14,9)=DN(3,0)*clhs84 + DN(3,1)*clhs150 + DN(3,2)*clhs189 + clhs301;
            lhs(14,10)=DN(3,0)*clhs90 + DN(3,1)*clhs155 + DN(3,2)*clhs191 + clhs309 + clhs319;
            lhs(14,11)=DN(2,2)*clhs317 + clhs310*clhs66 - clhs311;
            lhs(14,12)=DN(3,0)*clhs98 + DN(3,1)*clhs160 + DN(3,2)*clhs196 + clhs324;
            lhs(14,13)=DN(3,0)*clhs110 + DN(3,1)*clhs164 + DN(3,2)*clhs198 + clhs327;
            lhs(14,14)=DN(3,0)*clhs116 + DN(3,1)*clhs169 + DN(3,2)*clhs200 + clhs10*clhs328 + clhs321;
            lhs(14,15)=DN(3,2)*clhs325;
            lhs(15,0)=DN(3,0)*clhs205 + clhs118;
            lhs(15,1)=DN(3,1)*clhs205 + clhs171;
            lhs(15,2)=DN(3,2)*clhs205 + clhs203;
            lhs(15,3)=clhs213;
            lhs(15,4)=DN(3,0)*clhs208 + clhs240;
            lhs(15,5)=DN(3,1)*clhs208 + clhs255;
            lhs(15,6)=DN(3,2)*clhs208 + clhs269;
            lhs(15,7)=clhs273;
            lhs(15,8)=DN(3,0)*clhs210 + clhs293;
            lhs(15,9)=DN(3,1)*clhs210 + clhs302;
            lhs(15,10)=DN(3,2)*clhs210 + clhs310;
            lhs(15,11)=clhs313;
            lhs(15,12)=DN(3,0)*clhs329;
            lhs(15,13)=DN(3,1)*clhs329;
            lhs(15,14)=DN(3,2)*clhs329;
            lhs(15,15)=clhs207*(clhs320 + clhs326 + clhs328);


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

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 =             rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs8 =             crhs3 + crhs7;
const double crhs9 =             stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs10 =             crhs8*(crhs9*h/stab_c1 + mu);
const double crhs11 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs12 =             N[0]*crhs11*rho;
const double crhs13 =             1.0/(crhs9*rho/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs14 =             1.0*crhs13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs6);
const double crhs15 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs16 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs17 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs18 =             rho*(crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs5*crhs7);
const double crhs19 =             1.0*crhs13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs16 + crhs17 + crhs18);
const double crhs20 =             N[1]*crhs11*rho;
const double crhs21 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs22 =             N[2]*crhs11*rho;
const double crhs23 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs10 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs6 - crhs12*crhs14 - crhs14*crhs15;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs10 - DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 - N[0]*crhs18 - crhs12*crhs19 - crhs15*crhs19;
            rhs[2]=-DN(0,0)*crhs14 - DN(0,1)*crhs19 - N[0]*crhs8;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs10 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs6 - crhs14*crhs20 - crhs14*crhs21;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs10 - DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 - N[1]*crhs18 - crhs19*crhs20 - crhs19*crhs21;
            rhs[5]=-DN(1,0)*crhs14 - DN(1,1)*crhs19 - N[1]*crhs8;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs10 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs6 - crhs14*crhs22 - crhs14*crhs23;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs10 - DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 - N[2]*crhs18 - crhs19*crhs22 - crhs19*crhs23;
            rhs[8]=-DN(2,0)*crhs14 - DN(2,1)*crhs19 - N[2]*crhs8;


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

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 =             DN(0,0)*v(0,0);
const double crhs4 =             DN(1,0)*v(1,0);
const double crhs5 =             DN(2,0)*v(2,0);
const double crhs6 =             DN(3,0)*v(3,0);
const double crhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs10 =             rho*(crhs7*(crhs3 + crhs4 + crhs5 + crhs6) + crhs8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs11 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs12 =             DN(0,1)*v(0,1);
const double crhs13 =             DN(1,1)*v(1,1);
const double crhs14 =             DN(2,1)*v(2,1);
const double crhs15 =             DN(3,1)*v(3,1);
const double crhs16 =             crhs11 + crhs12 + crhs13 + crhs14 + crhs15 + crhs3 + crhs4 + crhs5 + crhs6;
const double crhs17 =             stab_c2*sqrt(pow(crhs7, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs18 =             crhs16*(crhs17*h/stab_c1 + mu);
const double crhs19 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs20 =             N[0]*crhs19*rho;
const double crhs21 =             1.0/(crhs17*rho/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 =             1.0*crhs21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs10 + crhs2);
const double crhs23 =             rho*(DN(0,0)*crhs7 + DN(0,1)*crhs8 + DN(0,2)*crhs9);
const double crhs24 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs25 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs26 =             rho*(crhs7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs8*(crhs12 + crhs13 + crhs14 + crhs15) + crhs9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs27 =             1.0*crhs21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs24 + crhs25 + crhs26);
const double crhs28 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs29 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs30 =             rho*(crhs11*crhs9 + crhs7*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs31 =             1.0*crhs21*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs28 + crhs29 + crhs30);
const double crhs32 =             N[1]*crhs19*rho;
const double crhs33 =             rho*(DN(1,0)*crhs7 + DN(1,1)*crhs8 + DN(1,2)*crhs9);
const double crhs34 =             N[2]*crhs19*rho;
const double crhs35 =             rho*(DN(2,0)*crhs7 + DN(2,1)*crhs8 + DN(2,2)*crhs9);
const double crhs36 =             N[3]*crhs19*rho;
const double crhs37 =             rho*(DN(3,0)*crhs7 + DN(3,1)*crhs8 + DN(3,2)*crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - crhs20*crhs22 - crhs22*crhs23;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs18 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 - crhs20*crhs27 - crhs23*crhs27;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs18 - DN(0,2)*stress[2] + N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs20*crhs31 - crhs23*crhs31;
            rhs[3]=-DN(0,0)*crhs22 - DN(0,1)*crhs27 - DN(0,2)*crhs31 - N[0]*crhs16;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - crhs22*crhs32 - crhs22*crhs33;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs18 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 - crhs27*crhs32 - crhs27*crhs33;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs18 - DN(1,2)*stress[2] + N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs31*crhs32 - crhs31*crhs33;
            rhs[7]=-DN(1,0)*crhs22 - DN(1,1)*crhs27 - DN(1,2)*crhs31 - N[1]*crhs16;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - crhs22*crhs34 - crhs22*crhs35;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs18 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 - crhs27*crhs34 - crhs27*crhs35;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs18 - DN(2,2)*stress[2] + N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs31*crhs34 - crhs31*crhs35;
            rhs[11]=-DN(2,0)*crhs22 - DN(2,1)*crhs27 - DN(2,2)*crhs31 - N[2]*crhs16;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs18 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs10 - N[3]*crhs2 - crhs22*crhs36 - crhs22*crhs37;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs18 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 - crhs27*crhs36 - crhs27*crhs37;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs18 - DN(3,2)*stress[2] + N[3]*crhs28 - N[3]*crhs29 - N[3]*crhs30 - crhs31*crhs36 - crhs31*crhs37;
            rhs[15]=-DN(3,0)*crhs22 - DN(3,1)*crhs27 - DN(3,2)*crhs31 - N[3]*crhs16;


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
	const auto& stress = rData.ShearStress;

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

    const double cV0 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double cV1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cV2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cV3 =             1.0/(rho*stab_c2*sqrt(pow(cV1, 2) + pow(cV2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV4 =             1.0*DNenr(0,0)*cV0*cV3*rho;
const double cV5 =             DN(0,0)*cV1 + DN(0,1)*cV2;
const double cV6 =             1.0*DNenr(0,0)*cV3*rho;
const double cV7 =             1.0*DNenr(1,0)*cV0*cV3*rho;
const double cV8 =             1.0*DNenr(1,0)*cV3*rho;
const double cV9 =             1.0*DNenr(2,0)*cV0*cV3*rho;
const double cV10 =             1.0*DNenr(2,0)*cV3*rho;
const double cV11 =             1.0*DNenr(0,1)*cV0*cV3*rho;
const double cV12 =             1.0*DNenr(0,1)*cV3*rho;
const double cV13 =             1.0*DNenr(1,1)*cV0*cV3*rho;
const double cV14 =             1.0*DNenr(1,1)*cV3*rho;
const double cV15 =             1.0*DNenr(2,1)*cV0*cV3*rho;
const double cV16 =             1.0*DNenr(2,1)*cV3*rho;
const double cV17 =             1.0*cV3;
const double cV18 =             DN(1,0)*cV1 + DN(1,1)*cV2;
const double cV19 =             DN(2,0)*cV1 + DN(2,1)*cV2;
            V(0,0)=-DN(0,0)*Nenr[0] + N[0]*cV4 + cV5*cV6;
            V(0,1)=-DN(0,0)*Nenr[1] + N[0]*cV7 + cV5*cV8;
            V(0,2)=-DN(0,0)*Nenr[2] + N[0]*cV9 + cV10*cV5;
            V(1,0)=-DN(0,1)*Nenr[0] + N[0]*cV11 + cV12*cV5;
            V(1,1)=-DN(0,1)*Nenr[1] + N[0]*cV13 + cV14*cV5;
            V(1,2)=-DN(0,1)*Nenr[2] + N[0]*cV15 + cV16*cV5;
            V(2,0)=cV17*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            V(2,1)=cV17*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            V(2,2)=cV17*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            V(3,0)=-DN(1,0)*Nenr[0] + N[1]*cV4 + cV18*cV6;
            V(3,1)=-DN(1,0)*Nenr[1] + N[1]*cV7 + cV18*cV8;
            V(3,2)=-DN(1,0)*Nenr[2] + N[1]*cV9 + cV10*cV18;
            V(4,0)=-DN(1,1)*Nenr[0] + N[1]*cV11 + cV12*cV18;
            V(4,1)=-DN(1,1)*Nenr[1] + N[1]*cV13 + cV14*cV18;
            V(4,2)=-DN(1,1)*Nenr[2] + N[1]*cV15 + cV16*cV18;
            V(5,0)=cV17*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            V(5,1)=cV17*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            V(5,2)=cV17*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            V(6,0)=-DN(2,0)*Nenr[0] + N[2]*cV4 + cV19*cV6;
            V(6,1)=-DN(2,0)*Nenr[1] + N[2]*cV7 + cV19*cV8;
            V(6,2)=-DN(2,0)*Nenr[2] + N[2]*cV9 + cV10*cV19;
            V(7,0)=-DN(2,1)*Nenr[0] + N[2]*cV11 + cV12*cV19;
            V(7,1)=-DN(2,1)*Nenr[1] + N[2]*cV13 + cV14*cV19;
            V(7,2)=-DN(2,1)*Nenr[2] + N[2]*cV15 + cV16*cV19;
            V(8,0)=cV17*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            V(8,1)=cV17*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            V(8,2)=cV17*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cH2 =             DN(0,0)*cH0 + DN(0,1)*cH1 + N[0]*bdf0;
const double cH3 =             1.0/(rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH4 =             1.0*DNenr(0,0)*cH3*rho;
const double cH5 =             1.0*DNenr(0,1)*cH3*rho;
const double cH6 =             1.0*cH3;
const double cH7 =             DN(1,0)*cH0 + DN(1,1)*cH1 + N[1]*bdf0;
const double cH8 =             DN(2,0)*cH0 + DN(2,1)*cH1 + N[2]*bdf0;
const double cH9 =             1.0*DNenr(1,0)*cH3*rho;
const double cH10 =             1.0*DNenr(1,1)*cH3*rho;
const double cH11 =             1.0*DNenr(2,0)*cH3*rho;
const double cH12 =             1.0*DNenr(2,1)*cH3*rho;
            H(0,0)=cH2*cH4;
            H(0,1)=cH2*cH5;
            H(0,2)=cH6*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1));
            H(0,3)=cH4*cH7;
            H(0,4)=cH5*cH7;
            H(0,5)=cH6*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1));
            H(0,6)=cH4*cH8;
            H(0,7)=cH5*cH8;
            H(0,8)=cH6*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1));
            H(1,0)=cH2*cH9;
            H(1,1)=cH10*cH2;
            H(1,2)=cH6*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1));
            H(1,3)=cH7*cH9;
            H(1,4)=cH10*cH7;
            H(1,5)=cH6*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1));
            H(1,6)=cH8*cH9;
            H(1,7)=cH10*cH8;
            H(1,8)=cH6*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1));
            H(2,0)=cH11*cH2;
            H(2,1)=cH12*cH2;
            H(2,2)=cH6*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1));
            H(2,3)=cH11*cH7;
            H(2,4)=cH12*cH7;
            H(2,5)=cH6*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1));
            H(2,6)=cH11*cH8;
            H(2,7)=cH12*cH8;
            H(2,8)=cH6*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1));
const double cKee3 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2));
            Kee(1,2)=cKee3;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee3;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2));


    const double crhs_ee0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs_ee1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs_ee2 =             1.0/(rho*stab_c2*sqrt(pow(crhs_ee0, 2) + pow(crhs_ee1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee3 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + crhs_ee0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + crhs_ee1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs_ee4 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + crhs_ee0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs_ee1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)));
            rhs_ee[0]=-crhs_ee2*(DNenr(0,0)*crhs_ee3 + DNenr(0,1)*crhs_ee4);
            rhs_ee[1]=-crhs_ee2*(DNenr(1,0)*crhs_ee3 + DNenr(1,1)*crhs_ee4);
            rhs_ee[2]=-crhs_ee2*(DNenr(2,0)*crhs_ee3 + DNenr(2,1)*crhs_ee4);


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
	const auto& stress = rData.ShearStress;

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

    const double cV0 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double cV1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cV2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cV3 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cV4 =             1.0/(rho*stab_c2*sqrt(pow(cV1, 2) + pow(cV2, 2) + pow(cV3, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cV5 =             1.0*DNenr(0,0)*cV0*cV4*rho;
const double cV6 =             DN(0,0)*cV1 + DN(0,1)*cV2 + DN(0,2)*cV3;
const double cV7 =             1.0*DNenr(0,0)*cV4*rho;
const double cV8 =             1.0*DNenr(1,0)*cV0*cV4*rho;
const double cV9 =             1.0*DNenr(1,0)*cV4*rho;
const double cV10 =             1.0*DNenr(2,0)*cV0*cV4*rho;
const double cV11 =             1.0*DNenr(2,0)*cV4*rho;
const double cV12 =             1.0*DNenr(3,0)*cV0*cV4*rho;
const double cV13 =             1.0*DNenr(3,0)*cV4*rho;
const double cV14 =             1.0*DNenr(0,1)*cV0*cV4*rho;
const double cV15 =             1.0*DNenr(0,1)*cV4*rho;
const double cV16 =             1.0*DNenr(1,1)*cV0*cV4*rho;
const double cV17 =             1.0*DNenr(1,1)*cV4*rho;
const double cV18 =             1.0*DNenr(2,1)*cV0*cV4*rho;
const double cV19 =             1.0*DNenr(2,1)*cV4*rho;
const double cV20 =             1.0*DNenr(3,1)*cV0*cV4*rho;
const double cV21 =             1.0*DNenr(3,1)*cV4*rho;
const double cV22 =             1.0*DNenr(0,2)*cV0*cV4*rho;
const double cV23 =             1.0*DNenr(0,2)*cV4*rho;
const double cV24 =             1.0*DNenr(1,2)*cV0*cV4*rho;
const double cV25 =             1.0*DNenr(1,2)*cV4*rho;
const double cV26 =             1.0*DNenr(2,2)*cV0*cV4*rho;
const double cV27 =             1.0*DNenr(2,2)*cV4*rho;
const double cV28 =             1.0*DNenr(3,2)*cV0*cV4*rho;
const double cV29 =             1.0*DNenr(3,2)*cV4*rho;
const double cV30 =             1.0*cV4;
const double cV31 =             DN(1,0)*cV1 + DN(1,1)*cV2 + DN(1,2)*cV3;
const double cV32 =             DN(2,0)*cV1 + DN(2,1)*cV2 + DN(2,2)*cV3;
const double cV33 =             DN(3,0)*cV1 + DN(3,1)*cV2 + DN(3,2)*cV3;
            V(0,0)=-DN(0,0)*Nenr[0] + N[0]*cV5 + cV6*cV7;
            V(0,1)=-DN(0,0)*Nenr[1] + N[0]*cV8 + cV6*cV9;
            V(0,2)=-DN(0,0)*Nenr[2] + N[0]*cV10 + cV11*cV6;
            V(0,3)=-DN(0,0)*Nenr[3] + N[0]*cV12 + cV13*cV6;
            V(1,0)=-DN(0,1)*Nenr[0] + N[0]*cV14 + cV15*cV6;
            V(1,1)=-DN(0,1)*Nenr[1] + N[0]*cV16 + cV17*cV6;
            V(1,2)=-DN(0,1)*Nenr[2] + N[0]*cV18 + cV19*cV6;
            V(1,3)=-DN(0,1)*Nenr[3] + N[0]*cV20 + cV21*cV6;
            V(2,0)=-DN(0,2)*Nenr[0] + N[0]*cV22 + cV23*cV6;
            V(2,1)=-DN(0,2)*Nenr[1] + N[0]*cV24 + cV25*cV6;
            V(2,2)=-DN(0,2)*Nenr[2] + N[0]*cV26 + cV27*cV6;
            V(2,3)=-DN(0,2)*Nenr[3] + N[0]*cV28 + cV29*cV6;
            V(3,0)=cV30*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            V(3,1)=cV30*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            V(3,2)=cV30*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            V(3,3)=cV30*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            V(4,0)=-DN(1,0)*Nenr[0] + N[1]*cV5 + cV31*cV7;
            V(4,1)=-DN(1,0)*Nenr[1] + N[1]*cV8 + cV31*cV9;
            V(4,2)=-DN(1,0)*Nenr[2] + N[1]*cV10 + cV11*cV31;
            V(4,3)=-DN(1,0)*Nenr[3] + N[1]*cV12 + cV13*cV31;
            V(5,0)=-DN(1,1)*Nenr[0] + N[1]*cV14 + cV15*cV31;
            V(5,1)=-DN(1,1)*Nenr[1] + N[1]*cV16 + cV17*cV31;
            V(5,2)=-DN(1,1)*Nenr[2] + N[1]*cV18 + cV19*cV31;
            V(5,3)=-DN(1,1)*Nenr[3] + N[1]*cV20 + cV21*cV31;
            V(6,0)=-DN(1,2)*Nenr[0] + N[1]*cV22 + cV23*cV31;
            V(6,1)=-DN(1,2)*Nenr[1] + N[1]*cV24 + cV25*cV31;
            V(6,2)=-DN(1,2)*Nenr[2] + N[1]*cV26 + cV27*cV31;
            V(6,3)=-DN(1,2)*Nenr[3] + N[1]*cV28 + cV29*cV31;
            V(7,0)=cV30*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            V(7,1)=cV30*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            V(7,2)=cV30*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            V(7,3)=cV30*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            V(8,0)=-DN(2,0)*Nenr[0] + N[2]*cV5 + cV32*cV7;
            V(8,1)=-DN(2,0)*Nenr[1] + N[2]*cV8 + cV32*cV9;
            V(8,2)=-DN(2,0)*Nenr[2] + N[2]*cV10 + cV11*cV32;
            V(8,3)=-DN(2,0)*Nenr[3] + N[2]*cV12 + cV13*cV32;
            V(9,0)=-DN(2,1)*Nenr[0] + N[2]*cV14 + cV15*cV32;
            V(9,1)=-DN(2,1)*Nenr[1] + N[2]*cV16 + cV17*cV32;
            V(9,2)=-DN(2,1)*Nenr[2] + N[2]*cV18 + cV19*cV32;
            V(9,3)=-DN(2,1)*Nenr[3] + N[2]*cV20 + cV21*cV32;
            V(10,0)=-DN(2,2)*Nenr[0] + N[2]*cV22 + cV23*cV32;
            V(10,1)=-DN(2,2)*Nenr[1] + N[2]*cV24 + cV25*cV32;
            V(10,2)=-DN(2,2)*Nenr[2] + N[2]*cV26 + cV27*cV32;
            V(10,3)=-DN(2,2)*Nenr[3] + N[2]*cV28 + cV29*cV32;
            V(11,0)=cV30*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            V(11,1)=cV30*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            V(11,2)=cV30*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            V(11,3)=cV30*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            V(12,0)=-DN(3,0)*Nenr[0] + N[3]*cV5 + cV33*cV7;
            V(12,1)=-DN(3,0)*Nenr[1] + N[3]*cV8 + cV33*cV9;
            V(12,2)=-DN(3,0)*Nenr[2] + N[3]*cV10 + cV11*cV33;
            V(12,3)=-DN(3,0)*Nenr[3] + N[3]*cV12 + cV13*cV33;
            V(13,0)=-DN(3,1)*Nenr[0] + N[3]*cV14 + cV15*cV33;
            V(13,1)=-DN(3,1)*Nenr[1] + N[3]*cV16 + cV17*cV33;
            V(13,2)=-DN(3,1)*Nenr[2] + N[3]*cV18 + cV19*cV33;
            V(13,3)=-DN(3,1)*Nenr[3] + N[3]*cV20 + cV21*cV33;
            V(14,0)=-DN(3,2)*Nenr[0] + N[3]*cV22 + cV23*cV33;
            V(14,1)=-DN(3,2)*Nenr[1] + N[3]*cV24 + cV25*cV33;
            V(14,2)=-DN(3,2)*Nenr[2] + N[3]*cV26 + cV27*cV33;
            V(14,3)=-DN(3,2)*Nenr[3] + N[3]*cV28 + cV29*cV33;
            V(15,0)=cV30*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            V(15,1)=cV30*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            V(15,2)=cV30*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            V(15,3)=cV30*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cH0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cH1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cH2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cH3 =             DN(0,0)*cH0 + DN(0,1)*cH1 + DN(0,2)*cH2 + N[0]*bdf0;
const double cH4 =             1.0/(rho*stab_c2*sqrt(pow(cH0, 2) + pow(cH1, 2) + pow(cH2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cH5 =             1.0*DNenr(0,0)*cH4*rho;
const double cH6 =             1.0*DNenr(0,1)*cH4*rho;
const double cH7 =             1.0*DNenr(0,2)*cH4*rho;
const double cH8 =             1.0*cH4;
const double cH9 =             DN(1,0)*cH0 + DN(1,1)*cH1 + DN(1,2)*cH2 + N[1]*bdf0;
const double cH10 =             DN(2,0)*cH0 + DN(2,1)*cH1 + DN(2,2)*cH2 + N[2]*bdf0;
const double cH11 =             DN(3,0)*cH0 + DN(3,1)*cH1 + DN(3,2)*cH2 + N[3]*bdf0;
const double cH12 =             1.0*DNenr(1,0)*cH4*rho;
const double cH13 =             1.0*DNenr(1,1)*cH4*rho;
const double cH14 =             1.0*DNenr(1,2)*cH4*rho;
const double cH15 =             1.0*DNenr(2,0)*cH4*rho;
const double cH16 =             1.0*DNenr(2,1)*cH4*rho;
const double cH17 =             1.0*DNenr(2,2)*cH4*rho;
const double cH18 =             1.0*DNenr(3,0)*cH4*rho;
const double cH19 =             1.0*DNenr(3,1)*cH4*rho;
const double cH20 =             1.0*DNenr(3,2)*cH4*rho;
            H(0,0)=cH3*cH5;
            H(0,1)=cH3*cH6;
            H(0,2)=cH3*cH7;
            H(0,3)=cH8*(DN(0,0)*DNenr(0,0) + DN(0,1)*DNenr(0,1) + DN(0,2)*DNenr(0,2));
            H(0,4)=cH5*cH9;
            H(0,5)=cH6*cH9;
            H(0,6)=cH7*cH9;
            H(0,7)=cH8*(DN(1,0)*DNenr(0,0) + DN(1,1)*DNenr(0,1) + DN(1,2)*DNenr(0,2));
            H(0,8)=cH10*cH5;
            H(0,9)=cH10*cH6;
            H(0,10)=cH10*cH7;
            H(0,11)=cH8*(DN(2,0)*DNenr(0,0) + DN(2,1)*DNenr(0,1) + DN(2,2)*DNenr(0,2));
            H(0,12)=cH11*cH5;
            H(0,13)=cH11*cH6;
            H(0,14)=cH11*cH7;
            H(0,15)=cH8*(DN(3,0)*DNenr(0,0) + DN(3,1)*DNenr(0,1) + DN(3,2)*DNenr(0,2));
            H(1,0)=cH12*cH3;
            H(1,1)=cH13*cH3;
            H(1,2)=cH14*cH3;
            H(1,3)=cH8*(DN(0,0)*DNenr(1,0) + DN(0,1)*DNenr(1,1) + DN(0,2)*DNenr(1,2));
            H(1,4)=cH12*cH9;
            H(1,5)=cH13*cH9;
            H(1,6)=cH14*cH9;
            H(1,7)=cH8*(DN(1,0)*DNenr(1,0) + DN(1,1)*DNenr(1,1) + DN(1,2)*DNenr(1,2));
            H(1,8)=cH10*cH12;
            H(1,9)=cH10*cH13;
            H(1,10)=cH10*cH14;
            H(1,11)=cH8*(DN(2,0)*DNenr(1,0) + DN(2,1)*DNenr(1,1) + DN(2,2)*DNenr(1,2));
            H(1,12)=cH11*cH12;
            H(1,13)=cH11*cH13;
            H(1,14)=cH11*cH14;
            H(1,15)=cH8*(DN(3,0)*DNenr(1,0) + DN(3,1)*DNenr(1,1) + DN(3,2)*DNenr(1,2));
            H(2,0)=cH15*cH3;
            H(2,1)=cH16*cH3;
            H(2,2)=cH17*cH3;
            H(2,3)=cH8*(DN(0,0)*DNenr(2,0) + DN(0,1)*DNenr(2,1) + DN(0,2)*DNenr(2,2));
            H(2,4)=cH15*cH9;
            H(2,5)=cH16*cH9;
            H(2,6)=cH17*cH9;
            H(2,7)=cH8*(DN(1,0)*DNenr(2,0) + DN(1,1)*DNenr(2,1) + DN(1,2)*DNenr(2,2));
            H(2,8)=cH10*cH15;
            H(2,9)=cH10*cH16;
            H(2,10)=cH10*cH17;
            H(2,11)=cH8*(DN(2,0)*DNenr(2,0) + DN(2,1)*DNenr(2,1) + DN(2,2)*DNenr(2,2));
            H(2,12)=cH11*cH15;
            H(2,13)=cH11*cH16;
            H(2,14)=cH11*cH17;
            H(2,15)=cH8*(DN(3,0)*DNenr(2,0) + DN(3,1)*DNenr(2,1) + DN(3,2)*DNenr(2,2));
            H(3,0)=cH18*cH3;
            H(3,1)=cH19*cH3;
            H(3,2)=cH20*cH3;
            H(3,3)=cH8*(DN(0,0)*DNenr(3,0) + DN(0,1)*DNenr(3,1) + DN(0,2)*DNenr(3,2));
            H(3,4)=cH18*cH9;
            H(3,5)=cH19*cH9;
            H(3,6)=cH20*cH9;
            H(3,7)=cH8*(DN(1,0)*DNenr(3,0) + DN(1,1)*DNenr(3,1) + DN(1,2)*DNenr(3,2));
            H(3,8)=cH10*cH18;
            H(3,9)=cH10*cH19;
            H(3,10)=cH10*cH20;
            H(3,11)=cH8*(DN(2,0)*DNenr(3,0) + DN(2,1)*DNenr(3,1) + DN(2,2)*DNenr(3,2));
            H(3,12)=cH11*cH18;
            H(3,13)=cH11*cH19;
            H(3,14)=cH11*cH20;
            H(3,15)=cH8*(DN(3,0)*DNenr(3,0) + DN(3,1)*DNenr(3,1) + DN(3,2)*DNenr(3,2));


    const double cKee0 =             1.0/(rho*stab_c2*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cKee1 =             cKee0*(DNenr(0,0)*DNenr(1,0) + DNenr(0,1)*DNenr(1,1) + DNenr(0,2)*DNenr(1,2));
const double cKee2 =             cKee0*(DNenr(0,0)*DNenr(2,0) + DNenr(0,1)*DNenr(2,1) + DNenr(0,2)*DNenr(2,2));
const double cKee3 =             cKee0*(DNenr(0,0)*DNenr(3,0) + DNenr(0,1)*DNenr(3,1) + DNenr(0,2)*DNenr(3,2));
const double cKee4 =             cKee0*(DNenr(1,0)*DNenr(2,0) + DNenr(1,1)*DNenr(2,1) + DNenr(1,2)*DNenr(2,2));
const double cKee5 =             cKee0*(DNenr(1,0)*DNenr(3,0) + DNenr(1,1)*DNenr(3,1) + DNenr(1,2)*DNenr(3,2));
const double cKee6 =             cKee0*(DNenr(2,0)*DNenr(3,0) + DNenr(2,1)*DNenr(3,1) + DNenr(2,2)*DNenr(3,2));
            Kee(0,0)=cKee0*(pow(DNenr(0,0), 2) + pow(DNenr(0,1), 2) + pow(DNenr(0,2), 2));
            Kee(0,1)=cKee1;
            Kee(0,2)=cKee2;
            Kee(0,3)=cKee3;
            Kee(1,0)=cKee1;
            Kee(1,1)=cKee0*(pow(DNenr(1,0), 2) + pow(DNenr(1,1), 2) + pow(DNenr(1,2), 2));
            Kee(1,2)=cKee4;
            Kee(1,3)=cKee5;
            Kee(2,0)=cKee2;
            Kee(2,1)=cKee4;
            Kee(2,2)=cKee0*(pow(DNenr(2,0), 2) + pow(DNenr(2,1), 2) + pow(DNenr(2,2), 2));
            Kee(2,3)=cKee6;
            Kee(3,0)=cKee3;
            Kee(3,1)=cKee5;
            Kee(3,2)=cKee6;
            Kee(3,3)=cKee0*(pow(DNenr(3,0), 2) + pow(DNenr(3,1), 2) + pow(DNenr(3,2), 2));


    const double crhs_ee0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs_ee1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs_ee2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs_ee3 =             1.0/(rho*stab_c2*sqrt(pow(crhs_ee0, 2) + pow(crhs_ee1, 2) + pow(crhs_ee2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs_ee4 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + DNenr(0,0)*penr[0] + DNenr(1,0)*penr[1] + DNenr(2,0)*penr[2] + DNenr(3,0)*penr[3] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + crhs_ee0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + crhs_ee1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs_ee2*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs_ee5 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + DNenr(0,1)*penr[0] + DNenr(1,1)*penr[1] + DNenr(2,1)*penr[2] + DNenr(3,1)*penr[3] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + crhs_ee0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs_ee1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + crhs_ee2*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs_ee6 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + DNenr(0,2)*penr[0] + DNenr(1,2)*penr[1] + DNenr(2,2)*penr[2] + DNenr(3,2)*penr[3] - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + crhs_ee0*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs_ee1*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs_ee2*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)));
            rhs_ee[0]=-crhs_ee3*(DNenr(0,0)*crhs_ee4 + DNenr(0,1)*crhs_ee5 + DNenr(0,2)*crhs_ee6);
            rhs_ee[1]=-crhs_ee3*(DNenr(1,0)*crhs_ee4 + DNenr(1,1)*crhs_ee5 + DNenr(1,2)*crhs_ee6);
            rhs_ee[2]=-crhs_ee3*(DNenr(2,0)*crhs_ee4 + DNenr(2,1)*crhs_ee5 + DNenr(2,2)*crhs_ee6);
            rhs_ee[3]=-crhs_ee3*(DNenr(3,0)*crhs_ee4 + DNenr(3,1)*crhs_ee5 + DNenr(3,2)*crhs_ee6);


    noalias(rV) += rData.Weight * V;
    noalias(rH) += rData.Weight * H;
    noalias(rKee) += rData.Weight * Kee;
    noalias(rRHS_ee) += rData.Weight * rhs_ee;
}

template <>
unsigned int TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeSplitting(
	TwoFluidNavierStokesData<3, 4>& rData,
	MatrixType& rShapeFunctions,
	ShapeFunctionDerivativesArrayType& rShapeDerivatives,
	std::vector< MatrixType>& rEnrichedShapeDerivatives,
	MatrixType& rEnrichedShapeFunctions)
{
	MatrixType coords(NumNodes, Dim);
	//fill coordinates
	for (unsigned int i = 0; i < NumNodes; i++)
	{
	    const array_1d<double, Dim > & xyz = this->GetGeometry()[i].Coordinates();
	    for (unsigned int j = 0; j < Dim; j++)
	        coords(i, j) = xyz[j];
	}
	Vector distances(NumNodes);
	for (int i = 0; i < NumNodes; i++)
		distances[i] = rData.Distance[i];
	unsigned int ndivisions = EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, rShapeDerivatives[0], distances, rData.PartitionsVolumes, rShapeFunctions, rData.PartitionsSigns, rEnrichedShapeDerivatives, rEnrichedShapeFunctions);
	return ndivisions;

}

template <>
unsigned int TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeSplitting(
	TwoFluidNavierStokesData<2, 3>& rData,
	MatrixType& rShapeFunctions,
	ShapeFunctionDerivativesArrayType& rShapeDerivatives,
	std::vector< MatrixType>& rEnrichedShapeDerivatives,
	MatrixType& rEnrichedShapeFunctions)
{

	MatrixType coords(NumNodes, Dim);

	//fill coordinates
	for (unsigned int i = 0; i < NumNodes; i++)
	{
		const array_1d<double, Dim > & xyz = this->GetGeometry()[i].Coordinates();
		for (unsigned int j = 0; j < Dim; j++)
			coords(i, j) = xyz[j];
	}

	Vector distances(NumNodes);
	for (int i = 0; i < NumNodes; i++)
		distances[i] = rData.Distance[i];
	unsigned int ndivisions = EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, rShapeDerivatives[0], distances, rData.PartitionsVolumes, rShapeFunctions, rData.PartitionsSigns, rEnrichedShapeDerivatives, rEnrichedShapeFunctions);
	return ndivisions;

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
    for (unsigned int igauss = 0; igauss < rData.PartitionsVolumes.size(); igauss++)
    {
        double wGauss = rData.PartitionsVolumes[igauss];

        if(rData.PartitionsSigns[igauss] >= 0) //check positive and negative volume
            positive_volume += wGauss;
        else
            negative_volume += wGauss;
    }
    const double Vol = positive_volume + negative_volume;



    double max_diag = 0.0;
    for(unsigned int k=0; k<Dim+1; k++)
        if(fabs(Kee_tot(k,k) ) > max_diag) max_diag = fabs(Kee_tot(k,k) );
    if(max_diag == 0) max_diag = 1.0;
    


    if(positive_volume/Vol < min_area_ratio)
    {
        for(unsigned int i=0; i<Dim+1; i++)
        {
            if(rData.Distance[i] >= 0.0)
            {
                Kee_tot(i,i) += 1000.0*max_diag;
            }
        }
    }
    if(negative_volume/Vol < min_area_ratio)
    {
        for(unsigned int i=0; i<Dim+1; i++)
        {
            if(rData.Distance[i] < 0.0)
            {
                Kee_tot(i,i) += 1000.0*max_diag;
            }
        }
    }

    //"weakly" impose continuity
    for(unsigned int i=0; i<Dim; i++)
    {
        const double di = fabs(rData.Distance[i]);

        for(unsigned int j=i+1; j<Dim+1; j++)
        {
            const double dj =  fabs(rData.Distance[j]);

            if(rData.Distance[i]* rData.Distance[j] < 0.0) //cut edge
            {
                double sum_d = di+dj;
                double Ni = dj/sum_d;
                double Nj = di/sum_d;

                double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
                Kee_tot(i,i) += penalty_coeff * Ni*Ni;
                Kee_tot(i,j) -= penalty_coeff * Ni*Nj;
                Kee_tot(j,i) -= penalty_coeff * Nj*Ni;
                Kee_tot(j,j) += penalty_coeff * Nj*Nj;

            }
        }
    }
	//add to LHS enrichment contributions
	MatrixType inverse_diag;
	inverse_diag.resize(NumNodes, NumNodes, false);
	bool inversion_successful = InvertMatrix<>(Kee_tot, inverse_diag);

    if(!inversion_successful )
    {
        KRATOS_WATCH(rData.Distance)
        KRATOS_WATCH(positive_volume/Vol)
        KRATOS_WATCH(negative_volume/Vol)
        KRATOS_WATCH(Kee_tot)
        KRATOS_THROW_ERROR(std::logic_error,"error in the inversion of the enrichment matrix for element ",this->Id());
    }

    const boost::numeric::ublas::bounded_matrix<double,4,16> tmp = prod(inverse_diag,Htot);
    noalias(rLeftHandSideMatrix) -= prod(Vtot,tmp);

    const array_1d<double,4> tmp2 = prod(inverse_diag, rhs_ee_tot);
    noalias(rRightHandSideVector) -= prod(Vtot,tmp2);

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
	for (unsigned int igauss = 0; igauss < rData.PartitionsVolumes.size(); igauss++)
	{
		double wGauss = rData.PartitionsVolumes[igauss];

		if (rData.PartitionsSigns[igauss] >= 0) //check positive and negative volume
			positive_volume += wGauss;
		else
			negative_volume += wGauss;
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
	bool inversion_successful = InvertMatrix<>(Kee_tot, inverse_diag);

	if (!inversion_successful)
	{
		KRATOS_WATCH(rData.Distance)
			KRATOS_WATCH(positive_volume / Vol)
			KRATOS_WATCH(negative_volume / Vol)
			KRATOS_WATCH(Kee_tot)
			KRATOS_THROW_ERROR(std::logic_error, "error in the inversion of the enrichment matrix for element ", this->Id());
	}

	const boost::numeric::ublas::bounded_matrix<double, 4, 16> tmp = prod(inverse_diag, Htot);
	noalias(rLeftHandSideMatrix) -= prod(Vtot, tmp);

	const array_1d<double, 4> tmp2 = prod(inverse_diag, rhs_ee_tot);
	noalias(rRightHandSideVector) -= prod(Vtot, tmp2);

}


template<>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::CalculateMaterialPropertiesAtGaussPoint(TwoFluidNavierStokesData<2, 3>& rData)
{
	double dist = 0.0;
	for (unsigned int i = 0; i < NumNodes; i++)
		dist += rData.N[i] * rData.Distance[i];

	double navg = 0.0;
	double density = 0.0;
	double viscosity = 0.0;
	for (unsigned int i = 0; i < NumNodes; i++)
	{
		if (dist * rData.Distance[i] > 0.0)
		{
			navg += 1.0;
			density += this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
			viscosity += this->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		}
	}

	rData.Density = density / navg;
	rData.DynamicViscosity = viscosity / navg;

	const double c_smag = this->GetValue(C_SMAGORINSKY);
	if (c_smag > 0.0)
	{
		unsigned int strain_size = 3;

		rData.ComputeStrain(strain_size);
		Vector& S = rData.StrainRate;
		double strain_rate_norm = std::sqrt(2.*S[0] * S[0] + 2.*S[1] * S[1] + S[2] * S[2]);

		double length_scale = c_smag*rData.ElementSize;
		length_scale *= length_scale; // square
		rData.EffectiveViscosity = rData.DynamicViscosity + 2.0*length_scale*strain_rate_norm;
	}
	else rData.EffectiveViscosity = rData.DynamicViscosity;

}
template<>
void TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::CalculateMaterialPropertiesAtGaussPoint(TwoFluidNavierStokesData<3, 4>& rData)
{

	double dist = 0.0;
	for (unsigned int i = 0; i < NumNodes; i++)
		dist += rData.N[i] * rData.Distance[i];

	double navg = 0.0;
	double density = 0.0;
	double viscosity = 0.0;
	for (unsigned int i = 0; i < NumNodes; i++)
	{
		if (dist * rData.Distance[i] > 0.0)
		{
			navg += 1.0;
			density += this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
			viscosity += this->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		}
	}

	rData.Density = density / navg;
	rData.DynamicViscosity = viscosity / navg;

	const double c_smag = this->GetValue(C_SMAGORINSKY);
	if (c_smag > 0.0)
	{
		unsigned int strain_size = 6;

		rData.ComputeStrain(strain_size);
		Vector& S = rData.StrainRate;
		double strain_rate_norm = std::sqrt(2.*S[0] * S[0] + 2.*S[1] * S[1] + 2.*S[2] * S[2] +
				S[3] * S[3] + S[4] * S[4] + S[5] * S[5]);

		double length_scale = c_smag*rData.ElementSize;
		length_scale *= length_scale; // square
		rData.EffectiveViscosity = rData.DynamicViscosity + 2.0*length_scale*strain_rate_norm;
	}
	else rData.EffectiveViscosity = rData.DynamicViscosity;
}


template< class TElementData>
template< class T>
bool TwoFluidNavierStokes<TElementData>::InvertMatrix(const T& input, T& inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    T A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if (res != 0)
        return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
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
