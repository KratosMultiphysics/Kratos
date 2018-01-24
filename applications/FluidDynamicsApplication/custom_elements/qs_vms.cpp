//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "qs_vms.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId):
    FluidElement<TElementData>(NewId)
{}

template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, const NodesArrayType& ThisNodes):
    FluidElement<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, GeometryType::Pointer pGeometry):
    FluidElement<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    FluidElement<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
QSVMS<TElementData>::~QSVMS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer QSVMS<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<QSVMS>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer QSVMS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<QSVMS>(NewId, pGeom, pProperties);
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(const Variable<double>& rVariable,
    double& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void QSVMS<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(const Variable<Vector>& rVariable,
    Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void QSVMS<TElementData>::Calculate(const Variable<Matrix>& rVariable,
    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int QSVMS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Extra variables
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    // Output variables (for Calculate() functions)
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_PRESSURE);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::GetValueOnIntegrationPoints(Variable<array_1d<double, 3 > > const& rVariable,
                                            std::vector<array_1d<double, 3 > >& rValues,
                                            ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            data.UpdateGeometryValues(GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

            this->SubscaleVelocity(data, rCurrentProcessInfo, rValues[g]);
        }
    }
    else if (rVariable == VORTICITY)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            this->IntegrationPointVorticity(ShapeDerivatives[g],rValues[g]);
        }
    }
}


template< class TElementData >
void QSVMS<TElementData>::GetValueOnIntegrationPoints(Variable<double> const& rVariable,
                                            std::vector<double>& rValues,
                                            ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_PRESSURE)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            data.UpdateGeometryValues(GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

            this->SubscalePressure(data,rCurrentProcessInfo,rValues[g]);
        }

    }
    else if (rVariable == Q_VALUE)
    {
		Vector GaussWeights;
		Matrix ShapeFunctions;
		ShapeFunctionDerivativesArrayType ShapeDerivatives;
		this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
		const unsigned int NumGauss = GaussWeights.size();

		rValues.resize(NumGauss);
		Matrix GradVel;

		// Loop on integration points
		for (unsigned int g = 0; g < NumGauss; g++)
		{
			GradVel = ZeroMatrix(Dim,Dim);
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

			// Compute velocity gradient
			for (unsigned int i=0; i < Dim; ++i)
				for (unsigned int j=0; j < Dim; ++j)
					for (unsigned int iNode=0; iNode < NumNodes; ++iNode)
					{
						array_1d<double,3>& Vel =
							this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
						GradVel(i,j) += Vel[i] * rDN_DX(iNode,j);
					}

			// Compute Q-value
			double qval = 0.0;
			for (unsigned int i=0; i < Dim; ++i)
				for (unsigned int j=0; j < Dim; ++j)
					qval += GradVel(i,j) * GradVel(j,i);

			qval *= -0.5;
			rValues[g] = qval;
		}
	}
	else if (rVariable == VORTICITY_MAGNITUDE)
	{
		Vector GaussWeights;
		Matrix ShapeFunctions;
		ShapeFunctionDerivativesArrayType ShapeDerivatives;
		this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
		const unsigned int NumGauss = GaussWeights.size();

		rValues.resize(NumGauss);
		
  		// Loop on integration points
		for (unsigned int g = 0; g < NumGauss; g++)
		{
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
			array_1d<double,3> Vorticity(3,0.0);

            this->IntegrationPointVorticity(rDN_DX,Vorticity);

			rValues[g] = sqrt(Vorticity[0] * Vorticity[0] + Vorticity[1] * Vorticity[1]
					+ Vorticity[2] * Vorticity[2]);
		}
	}
}

template <class TElementData>
void QSVMS<TElementData>::GetValueOnIntegrationPoints(Variable<array_1d<double, 6>> const& rVariable,
                                                    std::vector<array_1d<double, 6>>& rValues,
                                                    ProcessInfo const& rCurrentProcessInfo)
{
}

template <class TElementData>
void QSVMS<TElementData>::GetValueOnIntegrationPoints(Variable<Vector> const& rVariable,
                                                    std::vector<Vector>& rValues,
                                                    ProcessInfo const& rCurrentProcessInfo)
{
}

template <class TElementData>
void QSVMS<TElementData>::GetValueOnIntegrationPoints(Variable<Matrix> const& rVariable,
                                                    std::vector<Matrix>& rValues,
                                                    ProcessInfo const& rCurrentProcessInfo)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template< class TElementData >
std::string QSVMS<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "QSVMS #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void QSVMS<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "QSVMS" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

template< class TElementData >
void QSVMS<TElementData>::ASGSMomentumResidual(
    TElementData& rData,
    array_1d<double,3> &rMomentumRes)
{
    const GeometryType rGeom = this->GetGeometry();

    Vector AGradN;
    array_1d<double, 3> convective_velocity =
        this->Interpolate(rData.Velocity, rData.N) -
        this->Interpolate(rData.MeshVelocity, rData.N);
    
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    double density = this->Interpolate(rData.Density,rData.N);
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);

        for (unsigned int d = 0; d < Dim; d++)
        {
            rMomentumRes[d] += density * ( rData.N[i]*(r_body_forces(i,d) - rAcc[d]) - AGradN[i]*r_velocities(i,d)) - rData.DN_DX(i,d)*r_pressures[i];
        }
    }
}


template< class TElementData >
void QSVMS<TElementData>::ASGSMassResidual(
    TElementData& rData,
    double &rMomentumRes)
{
    this->MassProjTerm(rData,rMomentumRes);
}


template< class TElementData >
void QSVMS<TElementData>::OSSMomentumResidual(
    TElementData& rData,
    array_1d<double,3> &rMomentumRes)
{
    this->MomentumProjTerm(rData,rMomentumRes);

    array_1d<double,3> momentum_projection = this->Interpolate(rData.MomentumProjection,rData.N);
    for (unsigned int d = 0; d < Dim; d++)
        rMomentumRes[d] -= momentum_projection[d];
}


template< class TElementData >
void QSVMS<TElementData>::OSSMassResidual(
    TElementData& rData,
    double &rMassRes)
{
    this->MassProjTerm(rData,rMassRes);
    double mass_projection = this->Interpolate(rData.MassProjection,rData.N);
    rMassRes -= mass_projection;
}


template< class TElementData >
void QSVMS<TElementData>::MomentumProjTerm(
    TElementData& rData,
    array_1d<double,3> &rMomentumRHS)
{
        array_1d<double, 3> convective_velocity =
        this->Interpolate(rData.Velocity, rData.N) -
        this->Interpolate(rData.MeshVelocity, rData.N);
    
    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    double density = this->Interpolate(rData.Density,rData.N);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        for (unsigned int d = 0; d < Dim; d++)
        {
            rMomentumRHS[d] += density * ( rData.N[i]*(rData.BodyForce(i,d) /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d)) - rData.DN_DX(i,d)*rData.Pressure[i];
        }
    }
}


template< class TElementData >
void QSVMS<TElementData>::MassProjTerm(
    TElementData& rData,
    double &rMassRHS)
{
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        for (unsigned int d = 0; d < Dim; d++)
            rMassRHS -= rData.DN_DX(i,d)*rData.Velocity(i,d);
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {

    // Call specialized implementation (it is on a helper class to avoid partial template specialization problems)
    Internals::SpecializedAddTimeIntegratedSystem<TElementData,
        TElementData::ElementManagesTimeIntegration>::AddSystem(this, rData,
        rLHS, rRHS);
}

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {
        KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;
    }

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {
        KRATOS_ERROR << "AddTimeIntegratedRHS is not implemented." << std::endl;
    }

template< class TElementData >
void QSVMS<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    // Interpolate nodal data on the integration point
    double ElemSize = this->ElementSize();

    double density = this->Interpolate(rData.Density,rData.N);
    double dynamic_viscosity = this->EffectiveViscosity(rData,ElemSize);
    array_1d<double,3> body_force = this->Interpolate(rData.BodyForce,rData.N);
    array_1d<double,3> momentum_projection = this->Interpolate(rData.MomentumProjection,rData.N);
    double mass_projection = this->Interpolate(rData.MassProjection,rData.N);

    double TauOne;
    double TauTwo;
    array_1d<double, 3> convective_velocity =
        this->Interpolate(rData.Velocity, rData.N) -
        this->Interpolate(rData.MeshVelocity, rData.N);
        
    this->CalculateStaticTau(rData,density,dynamic_viscosity,convective_velocity,ElemSize,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);
    
    // Multiplying some quantities by density to have correct units
    body_force *= density; // Force per unit of volume
    AGradN *= density; // Convective term is always multiplied by density

    // Temporary containers
    double K,G,PDivV;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            //K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = rData.N[i]*AGradN[j];
            K += AGradN[i]*TauOne*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
            K *= rData.Weight;

            // q-p stabilization block (initialize result)
            double laplacian = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < Dim; d++)
            //    K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < Dim; d++)
            {
                //K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);
                rLHS(row+d,col+d) += K;

                // v * Grad(p) block
                G = TauOne * AGradN[i] * rData.DN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(row+d,col+Dim) += rData.Weight * (G - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(col+Dim,row+d) += rData.Weight * (G + PDivV);

                // q-p stabilization block
                laplacian += rData.DN_DX(i,d) * rData.DN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)

                for (unsigned int e = 0; e < Dim; e++) // Stabilization: Div(v) * TauTwo * Div(u)
                    rLHS(row+d,col+e) += rData.Weight*TauTwo*rData.DN_DX(i,d)*rData.DN_DX(j,e);
            }

            // Write q-p term
            rLHS(row+Dim,col+Dim) += rData.Weight*TauOne*laplacian;
        }

        // RHS terms
        double forcing = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*BodyForce
            rRHS[row+d] += rData.Weight * TauOne * AGradN[i] * ( body_force[d] - momentum_projection[d]); // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            rRHS[row+d] -= rData.Weight * TauTwo * rData.DN_DX(i,d) * mass_projection;
            forcing += rData.DN_DX(i, d) * (body_force[d] - momentum_projection[d]);
        }
        rRHS[row + Dim] += rData.Weight * TauOne * forcing; // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    Internals::AddViscousTerm<Dim>(dynamic_viscosity,rData.Weight,rData.DN_DX,rLHS);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    double density = this->Interpolate(rData.Density,rData.N);

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;
            const double Mij = rData.Weight * density * rData.N[i] * rData.N[j];
            for (unsigned int d = 0; d < Dim; d++)
                rMassMatrix(row+d,col+d) += Mij;
        }
    }

    /* Note on OSS and full projection: Riccardo says that adding the terms provided by
     * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble whith the Bossak scheme:
     * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
     * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
     */
    if ( rData.UseOSS != 1.0 )
        this->AddMassStabilization(rData,rMassMatrix);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    double ElemSize = this->ElementSize();

    double density = this->Interpolate(rData.Density,rData.N);
    double dynamic_viscosity = this->EffectiveViscosity(rData,ElemSize);

    double TauOne;
    double TauTwo;
    array_1d<double,3> convective_velocity = this->Interpolate(rData.Velocity,rData.N) - this->Interpolate(rData.MeshVelocity,rData.N);
    this->CalculateStaticTau(rData,density,dynamic_viscosity,convective_velocity,ElemSize,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= density; // Convective term is always multiplied by density

    // Temporary container
    double K;
    double weight = rData.Weight * TauOne * density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            K = weight * AGradN[i] * rData.N[j];

            for (unsigned int d = 0; d < Dim; d++)
            {
                rMassMatrix(row+d,col+d) += K;
                rMassMatrix(row+Dim,col+d) += weight*rData.DN_DX(i,d)*rData.N[j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
double QSVMS<TElementData>::EffectiveViscosity(
    TElementData& rData, double ElementSize) {
    
    double c_s = rData.CSmagorinsky;
    double viscosity = this->Interpolate(rData.DynamicViscosity, rData.N);

    if (c_s != 0.0) {
        const double density = this->Interpolate(rData.Density, rData.N);
        const auto& r_velocities = rData.Velocity;
        const auto& r_dndx = rData.DN_DX;

        // Calculate Symetric gradient
        MatrixType strain_rate = ZeroMatrix(Dim, Dim);
        for (unsigned int n = 0; n < NumNodes; ++n) {
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    strain_rate(i, j) +=
                        0.5 * (r_dndx(n, j) * r_velocities(n, i) +
                                  r_dndx(n, i) * r_velocities(n, j));
        }

        // Norm of symetric gradient
        double strain_rate_norm = 0.0;
        for (unsigned int i = 0; i < Dim; ++i)
            for (unsigned int j = 0; j < Dim; ++j)
                strain_rate_norm += strain_rate(i, j) * strain_rate(i, j);
        strain_rate_norm = sqrt(2.0 * strain_rate_norm);

        // Nu_sgs = (c_s * Delta)^2 * (2*Sij*Sij)^(1/2)
        viscosity +=
            density * c_s * c_s * ElementSize * ElementSize * strain_rate_norm;
    }

    return viscosity;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::CalculateStaticTau(
    const TElementData& rData,
    double Density,
    double DynamicViscosity,
    const array_1d<double,3> &Velocity,
    double ElemSize,
    double &TauOne,
    double &TauTwo)
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double InvTau = c1 * DynamicViscosity / (ElemSize*ElemSize) + Density * ( rData.DynamicTau/rData.DeltaTime + c2 * velocity_norm / ElemSize );
    TauOne = 1.0/InvTau;
    TauTwo = DynamicViscosity + c2 * Density * velocity_norm * ElemSize / c1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    GeometryType& rGeom = this->GetGeometry();
    VectorType MomentumRHS = ZeroVector(NumNodes * Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        data.UpdateGeometryValues(GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

        array_1d<double, 3> MomentumRes(3, 0.0);
        double MassRes = 0.0;

        this->MomentumProjTerm(data, MomentumRes);
        this->MassProjTerm(data,MassRes);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double W = data.Weight*data.N[i];
            unsigned int Row = i*Dim;
            for (unsigned int d = 0; d < Dim; d++)
                MomentumRHS[Row+d] += W*MomentumRes[d];
            NodalArea[i] += W;
            MassRHS[i] += W*MassRes;
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    unsigned int Row = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = rGeom[i].FastGetSolutionStepValue(ADVPROJ);
        for (unsigned int d = 0; d < Dim; ++d)
            rMomValue[d] += MomentumRHS[Row++];
        rGeom[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        rGeom[i].UnSetLock(); // Free the node for other threads
    }
}

template< class TElementData >
void QSVMS<TElementData>::SubscaleVelocity(
    TElementData& rData,
    const ProcessInfo &rProcessInfo,
    array_1d<double,3> &rVelocitySubscale)
{
    double ElemSize = this->ElementSize();
    double dynamic_viscosity = this->EffectiveViscosity(rData,ElemSize);

    double TauOne;
    double TauTwo;
    double density = this->Interpolate(rData.Density,rData.N);
    array_1d<double,3> convective_velocity = this->Interpolate(rData.Velocity,rData.N) - this->Interpolate(rData.MeshVelocity,rData.N);
    this->CalculateStaticTau(rData,density,dynamic_viscosity,convective_velocity,ElemSize,TauOne,TauTwo);

    array_1d<double,3> Residual(3,0.0);

    if (rData.UseOSS != 1.0)
        this->ASGSMomentumResidual(rData,Residual);
    else
        this->OSSMomentumResidual(rData,Residual);

    rVelocitySubscale = TauOne*Residual;
}

template< class TElementData >
void QSVMS<TElementData>::SubscalePressure(
        TElementData& rData,
        const ProcessInfo& rProcessInfo,
        double &rPressureSubscale)
{
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double ElemSize = this->ElementSize();
    double dynamic_viscosity = this->EffectiveViscosity(rData,ElemSize);

    double TauOne;
    double TauTwo;
    double density = this->Interpolate(rData.Density,rData.N);
    array_1d<double, 3> convective_velocity =
        this->Interpolate(rData.Velocity, rData.N) -
        this->Interpolate(rData.MeshVelocity, rData.N);
    this->CalculateStaticTau(rData, density, dynamic_viscosity, convective_velocity, ElemSize,
                             TauOne, TauTwo);

    double Residual = 0.0;

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMassResidual(rData,Residual);
    else
        this->OSSMassResidual(rData,Residual);

    rPressureSubscale = TauTwo*Residual;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void QSVMS<TElementData>::save(Serializer& rSerializer) const
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void QSVMS<TElementData>::load(Serializer& rSerializer)
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Internals
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace Internals {

///////////////////////////////////////////////////////////////////////////////////////////////////
// For Dim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template <>
void AddViscousTerm<2>(double DynamicViscosity,
                       double GaussWeight,
                       const Kratos::Matrix& rDN_DX,
                       Kratos::Matrix& rLHS)
{
    double weight = GaussWeight * DynamicViscosity;

    constexpr double four_thirds = 4.0 / 3.0;
    constexpr double minus_two_thirds = -2.0 / 3.0;

    const unsigned int num_nodes = rDN_DX.size1();
    const unsigned int block_size = 3;

    for (unsigned int a = 0; a < num_nodes; ++a)
    {
        unsigned int row = a*block_size;
        for (unsigned int b = 0; b < num_nodes; ++b)
        {
            unsigned int col = b*block_size;

            // First row
            rLHS(row,col) += weight * ( four_thirds * rDN_DX(a,0) * rDN_DX(b,0) + rDN_DX(a,1) * rDN_DX(b,1) );
            rLHS(row,col+1) += weight * ( minus_two_thirds * rDN_DX(a,0) * rDN_DX(b,1) + rDN_DX(a,1) * rDN_DX(b,0) );

            // Second row
            rLHS(row+1,col) += weight * ( minus_two_thirds * rDN_DX(a,1) * rDN_DX(b,0) + rDN_DX(a,0) * rDN_DX(b,1) );
            rLHS(row+1,col+1) += weight * ( four_thirds * rDN_DX(a,1) * rDN_DX(b,1) + rDN_DX(a,0) * rDN_DX(b,0) );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For Dim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template <>
void AddViscousTerm<3>(double DynamicViscosity,
                       double GaussWeight,
                       const Kratos::Matrix& rDN_DX,
                       Kratos::Matrix& rLHS)
{
    double weight = GaussWeight * DynamicViscosity;

    constexpr double one_third = 1.0 / 3.0;
    constexpr double minus_two_thirds = -2.0 / 3.0;

    const unsigned int num_nodes = rDN_DX.size1();
    const unsigned int block_size = 4;

    unsigned int row(0),col(0);

    for (unsigned int i = 0; i < num_nodes; ++i)
    {
        row = i*block_size;
        for (unsigned int j = 0; j < num_nodes; ++j)
        {
            col = j*block_size;
            // (dN_i/dx_k dN_j/dx_k)
            const double diag =  rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

            // First row
            rLHS(row,col) += weight * ( one_third * rDN_DX(i,0) * rDN_DX(j,0) + diag );
            rLHS(row,col+1) += weight * ( minus_two_thirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );
            rLHS(row,col+2) += weight * ( minus_two_thirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0) );

            // Second row
            rLHS(row+1,col) += weight * ( minus_two_thirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rLHS(row+1,col+1) += weight * ( one_third * rDN_DX(i,1) * rDN_DX(j,1) + diag );
            rLHS(row+1,col+2) += weight * ( minus_two_thirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1) );

            // Third row
            rLHS(row+2,col) += weight * ( minus_two_thirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2) );
            rLHS(row+2,col+1) += weight * ( minus_two_thirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2) );
            rLHS(row+2,col+2) += weight * ( one_third * rDN_DX(i,2) * rDN_DX(j,2) + diag );
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For Standard data: Time integration is not available
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, false>::AddSystem(
    QSVMS<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
    KRATOS_TRY;
    KRATOS_ERROR << "Trying to use time-integrated element functions with a "
                    "data type that does not know previous time step data"
                 << std::endl;
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Specialized time integration
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, true>::AddSystem(
    QSVMS<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
        Matrix mass_matrix = ZeroMatrix(rLHS.size1(),rLHS.size2());
        Matrix velocity_lhs = ZeroMatrix(rLHS.size1(),rLHS.size2());

        pElement->AddVelocitySystem(rData,velocity_lhs,rRHS);
        pElement->AddMassLHS(rData,mass_matrix);

        noalias(rLHS) += rData.bdf0*mass_matrix + velocity_lhs;
        
        Vector values = ZeroVector(rRHS.size());
        Vector acceleration = ZeroVector(rRHS.size());

        int LocalIndex = 0;
        const auto& r_velocities = rData.Velocity;
        const auto& r_velocities_step1 = rData.Velocity_OldStep1;
        const auto& r_velocities_step2 = rData.Velocity_OldStep2;
        const auto& r_pressures = rData.Pressure;

        for (unsigned int i = 0; i < TElementData::NumNodes; ++i) {
            for (unsigned int d = 0; d < TElementData::Dim; ++d)  {
                values[LocalIndex] = r_velocities(i,d);
                // Velocity Dofs
                acceleration[LocalIndex] = rData.bdf0*r_velocities(i,d);
                acceleration[LocalIndex] += rData.bdf1*r_velocities_step1(i,d);
                acceleration[LocalIndex] += rData.bdf2*r_velocities_step2(i,d);
                ++LocalIndex;
            }
            values[LocalIndex] = r_pressures[i];
            ++LocalIndex;
        }

        noalias(rRHS) -= prod(velocity_lhs,values);
        noalias(rRHS) -= prod(mass_matrix,acceleration);
}

} // namespace Internals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class QSVMS< QSVMSData<2,3> >;
template class QSVMS< QSVMSData<3,4> >;

template class QSVMS< TimeIntegratedQSVMSData<2,3> >;
template class QSVMS< TimeIntegratedQSVMSData<3,4> >;


} // namespace Kratos
