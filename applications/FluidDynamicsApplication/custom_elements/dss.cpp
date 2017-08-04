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
void DSS<TElementData>::ASGSMomentumResidual(double GaussIndex,
                                     const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     array_1d<double,3> &rMomentumRes)
{
    const GeometryType rGeom = this->GetGeometry();
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rBodyForce = rGeom[i].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double Press = rGeom[i].FastGetSolutionStepValue(PRESSURE);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
        {
            rMomentumRes[d] += Density * ( rN[i]*(rBodyForce[d] - rAcc[d]) - AGradN[i]*rVel[d]) - rDN_DX(i,d)*Press;
        }
    }
}


template< class TElementData >
void DSS<TElementData>::ASGSMassResidual(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 double &rMomentumRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);
}


template< class TElementData >
void DSS<TElementData>::OSSMomentumResidual(double GaussIndex,
                                    const ShapeFunctionsType &rN,
                                    const ShapeFunctionDerivativesType &rDN_DX,
                                    array_1d<double,3> &rMomentumRes)
{
    this->MomentumProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);

    const GeometryType& rGeom = this->GetGeometry();

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rProj = rGeom[i].FastGetSolutionStepValue(ADVPROJ);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMomentumRes[d] -= rN[i]*rProj[d];
    }
}


template< class TElementData >
void DSS<TElementData>::OSSMassResidual(double GaussIndex,
                                const ShapeFunctionsType &rN,
                                const ShapeFunctionDerivativesType &rDN_DX,
                                double &rMassRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMassRes);

    const GeometryType& rGeom = this->GetGeometry();

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const double Proj = rGeom[i].FastGetSolutionStepValue(DIVPROJ);
        rMassRes -= rN[i]*Proj;
    }
}


template< class TElementData >
void DSS<TElementData>::MomentumProjTerm(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 array_1d<double,3> &rMomentumRHS)
{
    const GeometryType& rGeom = this->GetGeometry();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rBodyForce = rGeom[i].FastGetSolutionStepValue(BODY_FORCE);
        //const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        //const array_1d<double,3>& rAccOld = rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
        //array_1d<double,3> BossakAcc = 1.3*rAcc - 0.3*rAccOld;
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double Press = rGeom[i].FastGetSolutionStepValue(PRESSURE);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
        {
            rMomentumRHS[d] += Density * ( rN[i]*(rBodyForce[d] /*- BossakAcc[d]*/ /*- rAcc[d]*/) - AGradN[i]*rVel[d]) - rDN_DX(i,d)*Press;
        }
    }
}


template< class TElementData >
void DSS<TElementData>::MassProjTerm(double GaussIndex,
                             const ShapeFunctionsType &rN,
                             const ShapeFunctionDerivativesType &rDN_DX,
                             double &rMassRHS)
{
    const GeometryType& rGeom = this->GetGeometry();

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMassRHS -= rDN_DX(i,d)*rVel[d];
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< class TElementData >
void DSS<TElementData>::ResolvedConvectiveVelocity(array_1d<double,3> &rConvVel, const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();

    array_1d<double,3> NodeVel = rGeom[0].FastGetSolutionStepValue(VELOCITY);
    NodeVel -= rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY);
    rConvVel = rN[0] * NodeVel;

    for (unsigned int i = 1; i < TElementData::NumNodes; i++)
    {
        NodeVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        NodeVel -= rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        rConvVel += rN[i] * NodeVel;
    }
}


template< class TElementData >
void DSS<TElementData>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
                                       const ShapeFunctionsType &rN,
                                       const array_1d<double,3> &rSubscaleVel)
{
    GeometryType& rGeom = this->GetGeometry();

    rConvVel = rSubscaleVel;

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
        rConvVel += rN[i]*(rGeom[i].FastGetSolutionStepValue(VELOCITY)-rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY));
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
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> BodyForce(3,0.0);
    this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);
    //this->BodyForceTest(rProcessInfo,rN,BodyForce);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj(3,0.0);
    double MassProj = 0;
    this->EvaluateInPoint(MomentumProj,ADVPROJ,rN);
    this->EvaluateInPoint(MassProj,DIVPROJ,rN);

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

    // Modulated gradient diffusion
    //this->ModulatedGradientDiffusion(rLHS,rDN_DX,GaussWeight);
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
    KRATOS_WATCH("NOTE TO SELF: ERROR HERE")
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
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // Multiplying some quantities by density to have correct units
    //Viscosity *= Density; // Dynamic viscosity
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
void DSS<TElementData>::SubscaleVelocity(unsigned int GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 const ProcessInfo &rProcessInfo,
                                 array_1d<double,3> &rVelocitySubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    array_1d<double,3> Residual(3,0.0);

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMomentumResidual(GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMomentumResidual(GaussIndex,rN,rDN_DX,Residual);

    rVelocitySubscale = TauOne*Residual;
}

template< class TElementData >
void DSS<TElementData>::SubscalePressure(unsigned int GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 const ProcessInfo &rProcessInfo,
                                 double &rPressureSubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    double Residual = 0.0;

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMassResidual(GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMassResidual(GaussIndex,rN,rDN_DX,Residual);

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
