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

#include "dss.h"
#include "fluid_element_data.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
DSS<TElementData>::DSS(IndexType NewId):
    FluidElement<TElementData>(NewId)
{}

template< class TElementData >
DSS<TElementData>::DSS(IndexType NewId, const NodesArrayType& ThisNodes):
    FluidElement<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
DSS<TElementData>::DSS(IndexType NewId, GeometryType::Pointer pGeometry):
    FluidElement<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
DSS<TElementData>::DSS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    FluidElement<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
DSS<TElementData>::~DSS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer DSS<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Element::Pointer(new DSS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}


template< class TElementData >
Element::Pointer DSS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Element::Pointer(new DSS(NewId, pGeom, pProperties));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int DSS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    // Generic geometry check
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);

    // Check that required variables are registered

    // Nodal data
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that required variables have been correctly registered.","");
    if(PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that required variables have been correctly registered.","");
    if(BODY_FORCE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(VISCOSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(MESH_VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check that the application was correctly registered.","");

    // Variables used in computing projections
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check that the application was correctly registered.","");
    if(NODAL_AREA.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check that the application was correctly registered.","");
    if(ADVPROJ.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ADVPROJ Key is 0. Check that the application was correctly registered.","");
    if(DIVPROJ.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DIVPROJ Key is 0. Check that the application was correctly registered.","");

    // Output variables (for Calculate() functions)
    if(SUBSCALE_VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"SUBSCALE_VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(SUBSCALE_PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"SUBSCALE_PRESSURE Key is 0. Check that the application was correctly registered.","");

    // Elemental data
    if(C_SMAGORINSKY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"C_SMAGORINSKY Key is 0. Check that the application was correctly registered.","");

    // Process Info data
    // none at the moment, consider OSS_SWITCH


    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        // Check that required nodal variables are included in SolutionStepData
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());

        if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(NODAL_AREA) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(ADVPROJ) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing ADVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DIVPROJ) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing DIVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());

        // Check that required dofs exist
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
    {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if (this->GetGeometry()[i].Z() != 0.0)
                KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
        }
    }

    return out;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DSS<TElementData>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                            std::vector<array_1d<double, 3 > >& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
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

        TElementData Data(this->GetGeometry());

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
            this->SubscaleVelocity(Data,g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
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
void DSS<TElementData>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                            std::vector<double>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
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

        TElementData Data(this->GetGeometry());

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            this->SubscalePressure(Data,g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
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
			GradVel = ZeroMatrix(TElementData::Dim,TElementData::Dim);
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

			// Compute velocity gradient
			for (unsigned int i=0; i < TElementData::Dim; ++i)
				for (unsigned int j=0; j < TElementData::Dim; ++j)
					for (unsigned int iNode=0; iNode < TElementData::NumNodes; ++iNode)
					{
						array_1d<double,3>& Vel =
							this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
						GradVel(i,j) += Vel[i] * rDN_DX(iNode,j);
					}

			// Compute Q-value
			double qval = 0.0;
			for (unsigned int i=0; i < TElementData::Dim; ++i)
				for (unsigned int j=0; j < TElementData::Dim; ++j)
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

			if(TElementData::Dim == 2)
			{
				for (unsigned int iNode = 0; iNode < TElementData::NumNodes; iNode++)
				{
					array_1d<double,3>& Vel =
						this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
					Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
				}
			}
			else
			{
				for (unsigned int iNode = 0; iNode < this->GetGeometry().size(); iNode++)
				{
					array_1d<double,3>& Vel =
						this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
					Vorticity[0] += Vel[2] * rDN_DX(iNode,1) - Vel[1] * rDN_DX(iNode,2);
					Vorticity[1] += Vel[0] * rDN_DX(iNode,2) - Vel[2] * rDN_DX(iNode,0);
					Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
				}
			}

			rValues[g] = sqrt(Vorticity[0] * Vorticity[0] + Vorticity[1] * Vorticity[1]
					+ Vorticity[2] * Vorticity[2]);
		}
	}
}


template< class TElementData >
void DSS<TElementData>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                            std::vector<array_1d<double, 6 > >& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< class TElementData >
void DSS<TElementData>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                            std::vector<Vector>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< class TElementData >
void DSS<TElementData>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                            std::vector<Matrix>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< class TElementData >
std::string DSS<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "DSS #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void DSS<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DSS" << TElementData::Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DSS<TElementData>::ASGSMomentumResidual(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    array_1d<double,3> &rMomentumRes)
{
    const GeometryType rGeom = this->GetGeometry();
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
        {
            rMomentumRes[d] += Density * ( rN[i]*(rData.BodyForce(i,d) - rAcc[d]) - AGradN[i]*rData.Velocity(i,d)) - rDN_DX(i,d)*rData.Pressure[i];
        }
    }
}


template< class TElementData >
void DSS<TElementData>::ASGSMassResidual(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    double &rMomentumRes)
{
    this->MassProjTerm(rData,GaussIndex,rN,rDN_DX,rMomentumRes);
}


template< class TElementData >
void DSS<TElementData>::OSSMomentumResidual(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    array_1d<double,3> &rMomentumRes)
{
    this->MomentumProjTerm(rData,GaussIndex,rN,rDN_DX,rMomentumRes);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMomentumRes[d] -= rN[i]*rData.MomentumProjection(i,d);
    }
}


template< class TElementData >
void DSS<TElementData>::OSSMassResidual(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    double &rMassRes)
{
    this->MassProjTerm(rData,GaussIndex,rN,rDN_DX,rMassRes);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        rMassRes -= rN[i]*rData.MassProjection[i];
    }
}


template< class TElementData >
void DSS<TElementData>::MomentumProjTerm(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    array_1d<double,3> &rMomentumRHS)
{
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);


    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        for (unsigned int d = 0; d < TElementData::Dim; d++)
        {
            rMomentumRHS[d] += Density * ( rN[i]*(rData.BodyForce(i,d) /*- BossakAcc[d]*/ /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d)) - rDN_DX(i,d)*rData.Pressure[i];
        }
    }
}


template< class TElementData >
void DSS<TElementData>::MassProjTerm(
    const TElementData& rData,
    double GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    double &rMassRHS)
{
    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMassRHS -= rDN_DX(i,d)*rData.Velocity(i,d);
    }

}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DSS<TElementData>::AddSystemTerms(
    const TElementData& rData,
    unsigned int GaussIndex,
    double GaussWeight,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    const ProcessInfo &rProcessInfo,
    MatrixType &rLHS,
    VectorType &rRHS)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> BodyForce(3,0.0);
    this->EvaluateInPoint(BodyForce,rData.BodyForce,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);

    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rData,rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj(3,0.0);
    double MassProj = 0;
    this->EvaluateInPoint(MomentumProj,rData.MomentumProjection,rN);
    this->EvaluateInPoint(MassProj,rData.MassProjection,rN);
    
    // Multiplying some quantities by density to have correct units
    Viscosity *= Density; // Dynamic viscosity
    BodyForce *= Density; // Force per unit of volume
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary containers
    double K,L,G,PDivV,qF;


    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        Row = i*TElementData::BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < TElementData::NumNodes; j++)
        {
            Col = j*TElementData::BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            //K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = rN[i]*AGradN[j];
            K += AGradN[i]*TauOne*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
            K *= GaussWeight;

            // q-p stabilization block (reset result)
            L = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < TElementData::Dim; d++)
            //    K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < TElementData::Dim; d++)
            {
                //K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);
                rLHS(Row+d,Col+d) += K;

                // v * Grad(p) block
                G = TauOne * AGradN[i] * rDN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rDN_DX(i,d) * rN[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TElementData::Dim) += GaussWeight * (G - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(Col+TElementData::Dim,Row+d) += GaussWeight * (G + PDivV);

                // q-p stabilization block
                L += rDN_DX(i,d) * rDN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)

                for (unsigned int e = 0; e < TElementData::Dim; e++) // Stabilization: Div(v) * TauTwo * Div(u)
                    rLHS(Row+d,Col+e) += GaussWeight*TauTwo*rDN_DX(i,d)*rDN_DX(j,e);
            }

            // Write q-p term
            rLHS(Row+TElementData::Dim,Col+TElementData::Dim) += GaussWeight*TauOne*L;
        }

        // RHS terms
        qF = 0.0;
        for (unsigned int d = 0; d < TElementData::Dim; ++d)
        {
            rRHS[Row+d] += GaussWeight * rN[i] * BodyForce[d]; // v*BodyForce
            rRHS[Row+d] += GaussWeight * TauOne * AGradN[i] * ( BodyForce[d] - MomentumProj[d]); // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            rRHS[Row+d] -= GaussWeight * TauTwo * rDN_DX(i,d) * MassProj;
            qF += rDN_DX(i, d) * (BodyForce[d] - MomentumProj[d]);
        }
        rRHS[Row + TElementData::Dim] += GaussWeight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    this->AddViscousTerm(Viscosity,GaussWeight,rDN_DX,rLHS);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DSS<TElementData>::AddMassTerms(
    const TElementData& rData,
    double GaussWeight,
    const ShapeFunctionsType &rN,
    MatrixType &rMassMatrix)
{
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    unsigned int Row = 0;
    unsigned int Col = 0;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        Row = i*TElementData::BlockSize;
        for (unsigned int j = 0; j < TElementData::NumNodes; j++)
        {
            Col = j*TElementData::BlockSize;
            const double Mij = GaussWeight * Density * rN[i] * rN[j];
            for (unsigned int d = 0; d < TElementData::Dim; d++)
                rMassMatrix(Row+d,Col+d) += Mij;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DSS<TElementData>::AddMassStabilization(
    const TElementData& rData,
    unsigned int GaussIndex,
    double GaussWeight,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    const ProcessInfo &rProcessInfo,
    MatrixType &rMassMatrix)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);

    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rData,rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary container
    double K;
    double W = GaussWeight * TauOne * Density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        Row = i*TElementData::BlockSize;

        for (unsigned int j = 0; j < TElementData::NumNodes; j++)
        {
            Col = j*TElementData::BlockSize;

            K = W * AGradN[i] * rN[j];

            for (unsigned int d = 0; d < TElementData::Dim; d++)
            {
                rMassMatrix(Row+d,Col+d) += K;
                rMassMatrix(Row+TElementData::Dim,Col+d) += W*rDN_DX(i,d)*rN[j];
            }
        }
    }
}

template< class TElementData >
void DSS<TElementData>::CalculateProjections()
{
    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    GeometryType& rGeom = this->GetGeometry();
    VectorType MomentumRHS = ZeroVector(TElementData::NumNodes*TElementData::Dim);
    VectorType MassRHS = ZeroVector(TElementData::NumNodes);
    VectorType NodalArea = ZeroVector(TElementData::NumNodes);

    TElementData Data(rGeom);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        array_1d<double,3> MomentumRes(3,0.0);
        double MassRes = 0.0;

        this->MomentumProjTerm(Data,g,rN,rDN_DX,MomentumRes);
        this->MassProjTerm(Data,g,rN,rDN_DX,MassRes);

        for (unsigned int i = 0; i < TElementData::NumNodes; i++)
        {
            double W = GaussWeight*rN[i];
            unsigned int Row = i*TElementData::Dim;
            for (unsigned int d = 0; d < TElementData::Dim; d++)
                MomentumRHS[Row+d] += W*MomentumRes[d];
            NodalArea[i] += W;
            MassRHS[i] += W*MassRes;
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    unsigned int Row = 0;
    for (SizeType i = 0; i < TElementData::NumNodes; ++i)
    {
        rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = rGeom[i].FastGetSolutionStepValue(ADVPROJ);
        for (unsigned int d = 0; d < TElementData::Dim; ++d)
            rMomValue[d] += MomentumRHS[Row++];
        rGeom[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        rGeom[i].UnSetLock(); // Free the node for other threads
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS< FluidElementData<2,3> >::AddViscousTerm(double DynamicViscosity,
                            double GaussWeight,
                            const ShapeFunctionDerivativesType &rDN_DX,
                            MatrixType &rLHS)
{
    double Weight = GaussWeight * DynamicViscosity;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int Row(0),Col(0);

    for (unsigned int a = 0; a < FluidElementData<2,3>::NumNodes; ++a)
    {
        Row = a*FluidElementData<2,3>::BlockSize;
        for (unsigned int b = 0; b < FluidElementData<2,3>::NumNodes; ++b)
        {
            Col = b*FluidElementData<2,3>::BlockSize;

            // First Row
            rLHS(Row,Col) += Weight * ( FourThirds * rDN_DX(a,0) * rDN_DX(b,0) + rDN_DX(a,1) * rDN_DX(b,1) );
            rLHS(Row,Col+1) += Weight * ( nTwoThirds * rDN_DX(a,0) * rDN_DX(b,1) + rDN_DX(a,1) * rDN_DX(b,0) );

            // Second Row
            rLHS(Row+1,Col) += Weight * ( nTwoThirds * rDN_DX(a,1) * rDN_DX(b,0) + rDN_DX(a,0) * rDN_DX(b,1) );
            rLHS(Row+1,Col+1) += Weight * ( FourThirds * rDN_DX(a,1) * rDN_DX(b,1) + rDN_DX(a,0) * rDN_DX(b,0) );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS< FluidElementData<3,4> >::AddViscousTerm(double DynamicViscosity,
                            double GaussWeight,
                            const ShapeFunctionDerivativesType &rDN_DX,
                            MatrixType &rLHS)
{
    double Weight = GaussWeight * DynamicViscosity;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int Row(0),Col(0);

    for (unsigned int i = 0; i < FluidElementData<3,4>::NumNodes; ++i)
    {
        Row = i*FluidElementData<3,4>::BlockSize;
        for (unsigned int j = 0; j < FluidElementData<3,4>::NumNodes; ++j)
        {
            Col = j*FluidElementData<3,4>::BlockSize;
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag =  rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

            // First Row
            rLHS(Row,Col) += Weight * ( OneThird * rDN_DX(i,0) * rDN_DX(j,0) + Diag );
            rLHS(Row,Col+1) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );
            rLHS(Row,Col+2) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0) );

            // Second Row
            rLHS(Row+1,Col) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rLHS(Row+1,Col+1) += Weight * ( OneThird * rDN_DX(i,1) * rDN_DX(j,1) + Diag );
            rLHS(Row+1,Col+2) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1) );

            // Third Row
            rLHS(Row+2,Col) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2) );
            rLHS(Row+2,Col+1) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2) );
            rLHS(Row+2,Col+2) += Weight * ( OneThird * rDN_DX(i,2) * rDN_DX(j,2) + Diag );
        }
    }
}


template< class TElementData >
void DSS<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    unsigned int GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    const ProcessInfo &rProcessInfo,
    array_1d<double,3> &rVelocitySubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rData,rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    array_1d<double,3> Residual(3,0.0);

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMomentumResidual(rData,GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMomentumResidual(rData,GaussIndex,rN,rDN_DX,Residual);

    rVelocitySubscale = TauOne*Residual;
}

template< class TElementData >
void DSS<TElementData>::SubscalePressure(
    const TElementData& rData,
    unsigned int GaussIndex,
    const ShapeFunctionsType &rN,
    const ShapeFunctionDerivativesType &rDN_DX,
    const ProcessInfo &rProcessInfo,
    double &rPressureSubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,rData.Density,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(rData,rN,ConvVel);

    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rData,rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    double Residual = 0.0;

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMassResidual(rData,GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMassResidual(rData,GaussIndex,rN,rDN_DX,Residual);

    rPressureSubscale = TauTwo*Residual;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void DSS<TElementData>::save(Serializer& rSerializer) const
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void DSS<TElementData>::load(Serializer& rSerializer)
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS< FluidElementData<2,3> >;
template class DSS< FluidElementData<3,4> >;

}
