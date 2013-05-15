#include "fractional_step_discontinuous.h"

namespace Kratos {







template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
        rLeftHandSideMatrix.resize(NumNodes,NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
        rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    Matrix NContainer;
    ShapeFunctionDerivativesArrayType DN_DX;
	VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
	const SizeType NumGauss = NContainer.size1();

    // Stabilization parameters
    double ElemSize = this->ElementSize();
    double TauOne;
    double TauTwo;

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& N = row(NContainer,g);
        const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

        // Evaluate required variables at the integration point
        double Density;
//        double Viscosity;
//        array_1d<double,3> Velocity(3,0.0);
//        array_1d<double,3> MeshVelocity(3,0.0);
        array_1d<double,3> BodyForce(3,0.0);
        array_1d<double,3> MomentumProjection(3,0.0);

        this->EvaluateInPoint(Density,DENSITY,N);
//        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
//        this->EvaluateInPoint(Velocity,VELOCITY,N);
//        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
//        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection,PRESS_PROJ,N);

//        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
//        double OldPressure;
//        this->EvaluateInPoint(OldPressure,PRESSURE,N,0);

        array_1d<double,TDim> OldPressureGradient(TDim,0.0);
        this->EvaluateGradientInPoint(OldPressureGradient,PRESSURE,rDN_DX);

//        // For ALE: convective velocity
//        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Stabilization parameters
        array_1d<double,3> ConvVel(3,0.0);
        this->EvaluateConvVelocity(ConvVel,N);
        double Viscosity = this->EffectiveViscosity(N,rDN_DX,ElemSize,rCurrentProcessInfo);
        this->CalculateTau(TauOne,TauTwo,ElemSize,ConvVel,Density,Viscosity,rCurrentProcessInfo);

//        // Evaluate convection operator Velocity * Grad(N)
//        Vector UGradN(NumNodes);
//        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

/*        double DivU;
        this->EvaluateDivergenceInPoint(DivU,VELOCITY,rDN_DX);*/
		array_1d<double,3> ugauss;
		this->EvaluateInPoint(ugauss,VELOCITY,N);
		for (SizeType i = 0; i < NumNodes; ++i)
		{
			double aux = 0.0;
			for (SizeType d = 0; d < TDim; ++d)
				aux += ugauss[d]*rDN_DX(i,d);
		
			rRightHandSideVector[i] += GaussWeight * aux;
		}

        // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
        const double LaplacianCoeff = 1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]) ;

        // Add convection, stabilization and RHS contributions to the local system equation
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // LHS contribution
            for (SizeType j = 0; j < NumNodes; ++j)
            {
                double Lij = 0.0;
                for (SizeType d = 0; d < TDim; ++d)
                    Lij += rDN_DX(i,d) * rDN_DX(j,d);
                Lij *= (LaplacianCoeff + TauOne);

                rLeftHandSideMatrix(i,j) += GaussWeight * Lij;
            }

            // RHS contribution

            // Velocity divergence
            //double RHSi = - N[i] * DivU;
			double RHSi = 0.0;
            for (SizeType d = 0; d < TDim; ++d)
            {
//                double Conv = UGradN[0] * rGeom[0].FastGetSolutionStepValue(VELOCITY)[d];
//                for (SizeType j = 1; j < NumNodes; ++j) Conv += UGradN[i] * rGeom[i].FastGetSolutionStepValue(VELOCITY)[d];
                // Momentum stabilization
                RHSi += rDN_DX(i,d) * TauOne * ( Density  * ( BodyForce[d]/* - Conv*/ ) - OldPressureGradient[d] - MomentumProjection[d] );
            }

            rRightHandSideVector[i] += GaussWeight * RHSi;
        }
    }
}

/*
 * protected FractionalStepDiscontinuous<TDim> functions
 */

template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
                                           Matrix& rNContainer,
                                           Vector& rGaussWeights)
{
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_2);
    rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2),false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    /*
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);

    // Initialize arrays to proper size
    rDN_DX.resize(NumGauss);
    rDetJ.resize(NumGauss);

    const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( GeometryData::GI_GAUSS_2 );

    // Temporary container for inverse of J
    Matrix InvJ;

    GeometryType::JacobiansType J;
    rGeom.Jacobian( J, GeometryData::GI_GAUSS_2 );

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        // calculate inverse of the jacobian and its determinant
        MathUtils<double>::InvertMatrix( J[g], InvJ, rDetJ[g] );

        // calculate the shape function derivatives in global coordinates
        rDN_DX[g].resize(NumNodes,TDim);
        noalias( rDN_DX[g] ) = prod( DN_De[g], InvJ );
    }*/
}







///*
// * FractionalStepDiscontinuous<TDim> static members
// */

//template< unsigned int TDim >
//static const FractionalStepDiscontinuous<TDim>::ShapeFunctionsType FractionalStepDiscontinuous<TDim>::msNg = FractionalStepDiscontinuous<TDim>::InitializeShapeFunctions();

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class FractionalStepDiscontinuous<2>;
template class FractionalStepDiscontinuous<3>;

}
