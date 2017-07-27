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

#include "fluid_element.h"
#include "fluid_element_data.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId):
    Element(NewId)

{}

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
FluidElement<TElementData>::~FluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "Attempting to Create base FluidElement<" << TElementData::Dim << "> instances, but this is an abstract element." << std::endl;
}


template< class TElementData >
void FluidElement<TElementData>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                 VectorType &rRightHandSideVector,
                                                 ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TElementData::Dim+1);

    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
    rRightHandSideVector = ZeroVector(LocalSize);
}


template< class TElementData >
void FluidElement<TElementData>::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TElementData::Dim+1);

    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
}


template< class TElementData >
void FluidElement<TElementData>::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TElementData::Dim+1);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rRightHandSideVector = ZeroVector(LocalSize);
}


template< class TElementData >
void FluidElement<TElementData>::CalculateLocalVelocityContribution(MatrixType &rDampMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes*(TElementData::Dim+1);

    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rDampMatrix = ZeroMatrix(LocalSize,LocalSize);
    rRightHandSideVector = ZeroVector(LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        this->AddSystemTerms(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rDampMatrix,rRightHandSideVector);
    }

    // Rewrite local contribution into residual form (A*dx = b - A*x)
    VectorType U = ZeroVector(LocalSize);
    int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double,3> &rVel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < TElementData::Dim; ++d) // Velocity Dofs
            U[LocalIndex++] = rVel[d];
        U[LocalIndex++] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
    }

    noalias(rRightHandSideVector) -= prod(rDampMatrix, U);
}


template< class TElementData >
void FluidElement<TElementData>::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes*(TElementData::Dim+1);

    // Resize and intialize output
    if( rMassMatrix.size1() != LocalSize )
        rMassMatrix.resize(LocalSize,LocalSize,false);

    rMassMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        this->AddMassTerms(GaussWeight,rN,rMassMatrix);

        /* Note on OSS and full projection: Riccardo says that adding the terms provided by
         * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
         * could help reduce the non-linearity of the coupling between projection and u,p
         * However, leaving them on gives a lot of trouble whith the Bossak scheme:
         * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
         * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
         */
        if ( rCurrentProcessInfo[OSS_SWITCH] != 1.0 )
            this->AddMassStabilization(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rMassMatrix);
    }
}


template< class TElementData >
void FluidElement<TElementData>::Calculate(const Variable<double> &rVariable,
                          double &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{

}


template< class TElementData >
void FluidElement<TElementData>::Calculate(const Variable<array_1d<double,3> > &rVariable,
                          array_1d<double,3> &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == ADVPROJ)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        VectorType MomentumRHS = ZeroVector(NumNodes*TElementData::Dim);
        VectorType MassRHS = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const double GaussWeight = GaussWeights[g];
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            array_1d<double,3> MomentumRes(3,0.0);
            double MassRes = 0.0;

            this->MomentumProjTerm(g,rN,rDN_DX,MomentumRes);
            this->MassProjTerm(g,rN,rDN_DX,MassRes);

            for (unsigned int i = 0; i < NumNodes; i++)
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
        for (SizeType i = 0; i < NumNodes; ++i)
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
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void FluidElement< FluidElementData<2,3> >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = 3*NumNodes;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void FluidElement< FluidElementData<2,3> >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
     const unsigned int NumNodes = rGeom.PointsNumber();
     const unsigned int LocalSize = 3*NumNodes;

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void FluidElement< FluidElementData<3,4> >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = 4*NumNodes;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void FluidElement< FluidElementData<3,4> >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
     const unsigned int NumNodes = rGeom.PointsNumber();
     const unsigned int LocalSize = 4*NumNodes;

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TElementData::Dim+1);

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< class TElementData >
void FluidElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TElementData::Dim+1);

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rValues[Index++] = rAcc[d];
        rValues[Index++] = 0.0; // skip pressure Dof
    }
}


template< class TElementData >
GeometryData::IntegrationMethod FluidElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}


template< class TElementData >
void FluidElement<TElementData>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
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

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
            this->SubscaleVelocity(g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
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
void FluidElement<TElementData>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            this->SubscalePressure(g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
        }

    }
    else if (rVariable == Q_VALUE)
    {
		Vector GaussWeights;
		Matrix ShapeFunctions;
		ShapeFunctionDerivativesArrayType ShapeDerivatives;
		this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
		const unsigned int NumNodes = this->GetGeometry().PointsNumber();
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
					for (unsigned int iNode=0; iNode < NumNodes; ++iNode)
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
		const unsigned int NumNodes = this->GetGeometry().PointsNumber();
		const unsigned int NumGauss = GaussWeights.size();

		rValues.resize(NumGauss);
		
  		// Loop on integration points
		for (unsigned int g = 0; g < NumGauss; g++)
		{
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
			array_1d<double,3> Vorticity(3,0.0);

			if(TElementData::Dim == 2)
			{
				for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
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
void FluidElement<TElementData>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                            std::vector<array_1d<double, 6 > >& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< class TElementData >
void FluidElement<TElementData>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                            std::vector<Vector>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< class TElementData >
void FluidElement<TElementData>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                            std::vector<Matrix>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int FluidElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);

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
std::string FluidElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidElement #" << Id();
    return buffer.str();
}


template< class TElementData >
void FluidElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluidElement" << TElementData::Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                      Matrix &rNContainer,
                                      ShapeFunctionDerivativesArrayType &rDN_DX)
{
    const GeometryData::IntegrationMethod IntMethod = this->GetIntegrationMethod();
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(IntMethod);

    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,IntMethod);

    rNContainer.resize(NumGauss,NumNodes,false);
    rNContainer = rGeom.ShapeFunctionsValues(IntMethod);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(IntMethod);

    rGaussWeights.resize(NumGauss,false);

    for (unsigned int g = 0; g < NumGauss; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

/** Calculate characteristic element length.
 * @return Minimum element height
 */
template< class TElementData >
double FluidElement<TElementData>::ElementSize()
{
    KRATOS_TRY;
    GeometryType& rGeom = this->GetGeometry();
    GeometryData::KratosGeometryFamily GeoFamily = rGeom.GetGeometryFamily();

    switch (GeoFamily)
    {
    case GeometryData::Kratos_Triangle:
    {
        /* Calculate node-edge distances */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();

        // node 0, edge 12
        double nx = -(y20-y10);
        double ny = x20-x10;
        double Hsq = x10*nx + y10*ny;
        Hsq *= Hsq / (nx*nx + ny*ny);

        // node 1, edge 20
        nx = -y20;
        ny = x20;
        double hsq = x10*nx + y10*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;

        // node 2, edge 10
        nx = -y10;
        ny = x10;
        hsq = x20*nx + y20*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Quadrilateral:
    {
        /* Calculate node-edge distances, assuming parallel faces */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();

        // node 0, edge 12
        double nx = -(y20-y10);
        double ny = x20-y10;
        double Hsq = x10*nx + y10*ny;
        Hsq *= Hsq / std::sqrt(nx*nx + ny*ny);

        // node 0, edge 23
        nx = -(y30-y20);
        ny = x30-y20;
        double hsq = x20*nx + y20*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Tetrahedra:
    {
        /* Calculate distances between each node and the opposite face */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();
        double z20 = rGeom[2].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        // face 123
        double nx = (y30-y10)*(z20-z10) - (z30-z10)*(y20-y10);
        double ny = (z30-z10)*(x20-x10) - (x30-x10)*(z20-z10);
        double nz = (x30-x10)*(y20-y10) - (y30-y10)*(x20-x10);
        double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
        Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2

        // face 230
        nx = y30*z20 - z30*y20;
        ny = z30*x20 - x30*z20;
        nz = x30*y20 - y30*x20;
        double hsq = x10*nx + y10*ny + z10*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 301
        nx = y10*z30 - z10*y30;
        ny = z10*x30 - x10*z30;
        nz = x10*y30 - y10*x30;
        hsq = x20*nx + y20*ny + z20*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 012
        nx = y10*z20 - z10*y20;
        ny = z10*x20 - x10*z20;
        nz = x10*y20 - y10*x20;
        hsq = x30*nx + y30*ny + z30*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Hexahedra:
    {
        /* Numbering assumes bottom face nodes 0123, top nodes 4567
         * considering the distance between a few face--node pairs:
         * face node
         *  034  1
         *  014  3
         *  013  4
         * This assumes parallel faces!
         * (otherwise all nodes should be checked against opposite faces)
         */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        double x40 = rGeom[4].X() - rGeom[0].X();
        double y40 = rGeom[4].Y() - rGeom[0].Y();
        double z40 = rGeom[4].Z() - rGeom[0].Z();


        // Face 034
        double nx = y30*z40 - z30*y40;
        double ny = z30*x40 - x30*z40;
        double nz = x30*y40 - y30*x40;
        double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
        Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2

        // face 014
        nx = y10*z40 - z10*y40;
        ny = z10*x40 - x10*z40;
        nz = x10*y40 - y10*x40;
        double hsq = x30*nx + y30*ny + z30*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 013
        nx = y10*z30 - z10*y30;
        ny = z10*x30 - x10*z30;
        nz = x10*y30 - y10*x30;
        hsq = x40*nx + y40*ny + z40*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    default:
    {
        KRATOS_THROW_ERROR(std::invalid_argument,"FluidElement::ElementSize not implemented for this geometry type","");
        return 0;
    }
    }
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::CalculateStaticTau(double Density,
                                   double KinematicVisc,
                                   const array_1d<double,3> &Velocity,
                                   double ElemSize,
                                   const ProcessInfo& rProcessInfo,
                                   double &TauOne,
                                   double &TauTwo)
{
    const double c1 = 8.0;
    const double c2 = 2.0;

    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TElementData::Dim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( c1 * KinematicVisc / (ElemSize*ElemSize) + c2 * VelNorm / ElemSize );
    TauOne = 1.0/InvTau;
    TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::ASGSMomentumResidual(double GaussIndex,
                                     const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     array_1d<double,3> &rMomentumRes)
{
    const GeometryType rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < NumNodes; i++)
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
void FluidElement<TElementData>::ASGSMassResidual(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 double &rMomentumRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);
}


template< class TElementData >
void FluidElement<TElementData>::OSSMomentumResidual(double GaussIndex,
                                    const ShapeFunctionsType &rN,
                                    const ShapeFunctionDerivativesType &rDN_DX,
                                    array_1d<double,3> &rMomentumRes)
{
    this->MomentumProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rProj = rGeom[i].FastGetSolutionStepValue(ADVPROJ);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMomentumRes[d] -= rN[i]*rProj[d];
    }
}


template< class TElementData >
void FluidElement<TElementData>::OSSMassResidual(double GaussIndex,
                                const ShapeFunctionsType &rN,
                                const ShapeFunctionDerivativesType &rDN_DX,
                                double &rMassRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMassRes);

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const double Proj = rGeom[i].FastGetSolutionStepValue(DIVPROJ);
        rMassRes -= rN[i]*Proj;
    }
}


template< class TElementData >
void FluidElement<TElementData>::MomentumProjTerm(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 array_1d<double,3> &rMomentumRHS)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < NumNodes; i++)
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
void FluidElement<TElementData>::MassProjTerm(double GaussIndex,
                             const ShapeFunctionsType &rN,
                             const ShapeFunctionDerivativesType &rDN_DX,
                             double &rMassRHS)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rMassRHS -= rDN_DX(i,d)*rVel[d];
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< class TElementData >
double FluidElement<TElementData>::EffectiveViscosity(const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     double ElemSize,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    const FluidElement* const_this = static_cast<const FluidElement*>(this);
    double Csmag = const_this->GetValue(C_SMAGORINSKY);

    double KinViscosity = 0.0;
    this->EvaluateInPoint(KinViscosity,VISCOSITY,rN);

    if (Csmag != 0.0 )
    {
        const unsigned int NumNodes = this->GetGeometry().PointsNumber();

        // Calculate Symetric gradient
        MatrixType S = ZeroMatrix(TElementData::Dim,TElementData::Dim);
        for (unsigned int n = 0; n < NumNodes; ++n)
        {
            const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < TElementData::Dim; ++i)
                for (unsigned int j = 0; j < TElementData::Dim; ++j)
                    S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
        }

        // Norm of symetric gradient
        double NormS = 0.0;
        for (unsigned int i = 0; i < TElementData::Dim; ++i)
            for (unsigned int j = 0; j < TElementData::Dim; ++j)
                NormS += S(i,j) * S(i,j);
        NormS = sqrt(2.0*NormS);

        // Nu_sgs = (Csmag * Delta)^2 * (2*Sij*Sij)^(1/2)
        KinViscosity += Csmag * Csmag * ElemSize * ElemSize * NormS;
    }

    return KinViscosity;
}


template< class TElementData >
void FluidElement<TElementData>::ResolvedConvectiveVelocity(array_1d<double,3> &rConvVel, const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    array_1d<double,3> NodeVel = rGeom[0].FastGetSolutionStepValue(VELOCITY);
    NodeVel -= rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY);
    rConvVel = rN[0] * NodeVel;

    for (unsigned int i = 1; i < NumNodes; i++)
    {
        NodeVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        NodeVel -= rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        rConvVel += rN[i] * NodeVel;
    }
}


template< class TElementData >
void FluidElement<TElementData>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
                                       const ShapeFunctionsType &rN,
                                       const array_1d<double,3> &rSubscaleVel)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    rConvVel = rSubscaleVel;

    for (unsigned int i = 0; i < NumNodes; i++)
        rConvVel += rN[i]*(rGeom[i].FastGetSolutionStepValue(VELOCITY)-rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    if(rResult.size() != NumNodes) rResult.resize(NumNodes,false);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < TElementData::Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

template< class TElementData >
void FluidElement<TElementData>::SubscaleVelocity(unsigned int GaussIndex,
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
void FluidElement<TElementData>::SubscalePressure(unsigned int GaussIndex,
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


template<>
void FluidElement< FluidElementData<2,3> >::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[2] += rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0];
    }
}


template<>
void FluidElement< FluidElementData<3,4> >::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[0] += rDN_DX(i,1)*rVelocity[2] - rDN_DX(i,2)*rVelocity[1];
        rVorticity[1] += rDN_DX(i,2)*rVelocity[0] - rDN_DX(i,0)*rVelocity[2];
        rVorticity[2] += rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void FluidElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
}


template< class TElementData >
void FluidElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class FluidElement< FluidElementData<2,3> >;
template class FluidElement< FluidElementData<3,4> >;

}
