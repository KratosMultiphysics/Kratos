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

	//this->InitializeGeometryData(data);

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

		//ESTO SE PUEDE METER EN EL DATA
		unsigned int npos=0, nneg=0;
		for (unsigned int i = 0; i < NumNodes; i++)
		{
		    if(data.Distance[i] > 0)
				data.NumPositiveNodes++;
		    else
				data.NumNegativeNodes++;
		}

		if (data.IsCut()) {
			VectorType volumes;
			VectorType signs(6); //ATTENTION: this shall be initialized of size 6
			std::vector< MatrixType > DNenr;
			MatrixType Nenr;
			data.NumberOfDivisions = ComputeSplitting(data, shape_functions, volumes, signs,  DNenr, Nenr);
			if (data.NumberOfDivisions == 1) {
				//cases exist when the element is like not subdivided due to the characteristics of the provided distance
				//in this cases the element is treated as AIR or FLUID depending on the side
				array_1d<double,NumNodes> Ncenter;
				for(unsigned int i=0; i<NumNodes; i++) Ncenter[i]=0.25;
				const double dgauss = inner_prod(data.Distance, Ncenter);
				if (dgauss > 0)
					data.CalculateAirMaterialResponse();
				else
					this->CalculateMaterialResponse(data);

				//ojito con above, no se que me estoy perdiendo

				this->AddTimeIntegratedSystem(
					data, rLeftHandSideMatrix, rRightHandSideVector);

			}
		}

		else {
			// Iterate over integration points to evaluate local contribution
			for (unsigned int g = 0; g < number_of_gauss_points; g++) {

				data.UpdateGeometryValues(gauss_weights[g], row(shape_functions, g),
					shape_derivatives[g]);
				if (data.IsAir())
					data.CalculateAirMaterialResponse();
				else
					this->CalculateMaterialResponse(data);
				//ojito con above, no se que me estoy perdiendo

				this->AddTimeIntegratedSystem(
					data, rLeftHandSideMatrix, rRightHandSideVector);
			}

		}

		
	}
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
    const auto& pn = rData.Pressure_OldStep1;
    const auto& pnn = rData.Pressure_OldStep2;
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
    const auto& pn = rData.Pressure_OldStep1;
    const auto& pnn = rData.Pressure_OldStep2;
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
	const auto& pn = rData.Pressure_OldStep1;
	const auto& pnn = rData.Pressure_OldStep2;
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
	const auto& pn = rData.Pressure_OldStep1;
	const auto& pnn = rData.Pressure_OldStep2;
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
unsigned int TwoFluidNavierStokes<TwoFluidNavierStokesData<3, 4>>::ComputeSplitting(
	TwoFluidNavierStokesData<3, 4>& rData,
	MatrixType& rShapeFunctions,
	VectorType& rVolumes,
	VectorType& rSigns,
	std::vector< MatrixType>& rDNenr,
	MatrixType& rNenr)
{

	MatrixType coords(NumNodes, Dim);

	//fill coordinates
	for (unsigned int i = 0; i < NumNodes; i++)
	{
	    const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
	    for (unsigned int j = 0; j < Dim; j++)
	        coords(i, j) = xyz[j];
	}

	unsigned int ndivisions = 1; // EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, rData.DN_DX, rData, rVolumes, rShapeFunctions, rSigns, rDNenr, rNenr);
	return ndivisions;

}

template <>
unsigned int TwoFluidNavierStokes<TwoFluidNavierStokesData<2, 3>>::ComputeSplitting(
	TwoFluidNavierStokesData<2, 3>& rData,
	MatrixType& rShapeFunctions,
	VectorType& rVolumes,
	VectorType& rSigns,
	std::vector< MatrixType>& rDNenr,
	MatrixType& rNenr)
{

	MatrixType coords(NumNodes, Dim);

	//fill coordinates
	for (unsigned int i = 0; i < NumNodes; i++)
	{
		const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
		for (unsigned int j = 0; j < Dim; j++)
			coords(i, j) = xyz[j];
	}

	unsigned int ndivisions = 1; // EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, rData.DN_DX, rData, rVolumes, rShapeFunctions, rSigns, rDNenr, rNenr);
	return ndivisions;

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


//void NavierStokesEnr2D::ComputeGaussPointEnrichmentContributions(
//    boost::numeric::ublas::bounded_matrix<double,3,9>& H,
//    boost::numeric::ublas::bounded_matrix<double,9,3>& V,
//    boost::numeric::ublas::bounded_matrix<double,3,3>&  Kee,
//    array_1d<double,3>& rhs_ee,
//    const element_data<3,2>& data,
//    const array_1d<double,3>& distances,
//    const array_1d<double,3>& Nenr,
//    const boost::numeric::ublas::bounded_matrix<double,3,3>& DNenr 
//    )
//    {
//        const int nnodes = 3;
//        const int dim = 2;
//        const int strain_size = 3;
//        
//        const double rho = inner_prod(data.N, data.rho);        // Density
//        const double nu = inner_prod(data.N, data.nu);          // Kinematic viscosity
//        const double h = data.h;                                // Characteristic element size
//
//        const double& bdf0 = data.bdf0;
//        const double& bdf1 = data.bdf1;
//        const double& bdf2 = data.bdf2;
//        const double& delta_t = data.delta_t;
//        const double& dyn_tau_coeff = data.dyn_tau_coeff;
//        const double& tau1_coeff = data.tau1_coeff;
//
//        const bounded_matrix<double,nnodes,dim>& v = data.v;
//        const bounded_matrix<double,nnodes,dim>& vn = data.vn;
//        const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
//        const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
//        const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
//        const bounded_matrix<double,nnodes,dim>& f = data.f;
//        const array_1d<double,nnodes>& p = data.p;
//        const array_1d<double,strain_size>& stress = data.stress;
//        
//        // Get constitutive matrix 
//        // const Matrix& C = data.C;
//        
//        // Get shape function values
//        const array_1d<double,nnodes>& N = data.N;
//        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
//        
//        const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
//        
//        const double vconv_norm = norm_2(vconv_gauss);
//                
//        // Stabilization parameters
//        const double tau1 = tau1_coeff/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*rho*nu)/(h*h));
//        // const double tau2 = (rho*nu) + 0.5*h*vconv_norm;
//        
//        // Auxiliary variables used in the calculation of the RHS
//        const array_1d<double,dim> f_gauss = prod(trans(f), N);
//        const array_1d<double,dim> grad_p = prod(trans(DN), p);
//        //~ const double p_gauss = inner_prod(N,p);
//        
//        //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
//        //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
//        //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);;
//        
//        array_1d<double,3> penr = ZeroVector(3); //penriched is considered to be zero as we do not want to store it
//
//        //substitute_enrichment_V
//        
//        //substitute_enrichment_H
//        
//        //substitute_enrichment_Kee
//        
//        //substitute_enrichment_rhs_ee
//        
//        
//    }





}
