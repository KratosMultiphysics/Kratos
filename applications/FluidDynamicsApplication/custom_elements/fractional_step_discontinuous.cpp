#include "fractional_step_discontinuous.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {



template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::AddMomentumSystemTerms(Matrix& rLHSMatrix,
                                                  Vector& rRHSVector,
                                                  const double Density,
                                                  const Vector& rConvOperator,
                                                  const array_1d<double,3>& rBodyForce,
                                                  const double OldPressure,
                                                  const double TauOne,
                                                  const double TauTwo,
                                                  const array_1d<double,3>& rMomentumProjection,
                                                  const double MassProjection,
                                                  const ShapeFunctionsType& rN,
                                                  const ShapeFunctionDerivativesType& rDN_DX,
                                                  const double Weight)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;
    
    array_1d<double,TDim+1> OldPressures;
    for (SizeType i = 0; i < NumNodes; ++i) 
      OldPressures[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
    
    array_1d<double,TDim> grad_p = prod(trans(rDN_DX),OldPressures);

    for (SizeType i = 0; i < NumNodes; ++i)
    {

        // Build RHS
        for (SizeType d = 0; d < TDim; ++d)
        {
            // Body force
            double RHSi = Density * rN[i] * rBodyForce[d];
            // Pressure gradient (integrated by parts)
            //RHSi += rDN_DX(i,d) * OldPressure;
	    RHSi -= rN[i]*grad_p[d];
	    
            // Momentum Stabilization
            RHSi += Density * rConvOperator[i] * TauOne * ( /*rBodyForce[d] + rOldPressureGradient[d]*/ - rMomentumProjection[d] );
            // Mass Stabilization
            RHSi -= rDN_DX(i,d) * TauTwo * MassProjection;

            rRHSVector[FirstRow+d] += Weight * RHSi;
        }

        // Build LHS
        for (SizeType j = 0; j < NumNodes; ++j)
        {
            // Convective term
            double Kij = Density * rN[i] * rConvOperator[j];

            // Streamline stabilization
            Kij += Density * rConvOperator[i] * TauOne * Density * rConvOperator[j];

            Kij *= Weight;

            for (SizeType d = 0; d < TDim; ++d)
                rLHSMatrix(FirstRow + d,FirstCol + d) += Kij;

            // Mass-GLS (TauTwo) stabiliziation term
            for (SizeType m = 0; m < TDim; ++m)
                for (SizeType n = 0; n < TDim; ++n)
                    rLHSMatrix(FirstRow+m,FirstCol+n) += Weight * rDN_DX(i,m) * TauTwo * rDN_DX(j,n);

            FirstCol += TDim;
        }
        FirstRow += TDim;
        FirstCol = 0;
    }
}



template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
        rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
        rRightHandSideVector.resize(NumNodes,false);

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
		//TauOne = 0.0;

//        // Evaluate convection operator Velocity * Grad(N)
//        Vector UGradN(NumNodes);
//        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

        //double DivU;
        //this->EvaluateDivergenceInPoint(DivU,VELOCITY,rDN_DX);
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
	
	//
	bool split_element = this->GetValue(SPLIT_ELEMENT);
	if (split_element == true)
	{
		const array_1d<double,3>& vel = this->GetValue(EMBEDDED_VELOCITY);
		const Vector& distances = this->GetValue(ELEMENTAL_DISTANCES);
			
		double Volume_tot;
		boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DXcontinuous;
		array_1d<double, 4 > Ncontinuous;
		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DXcontinuous, Ncontinuous, Volume_tot);
		
		array_1d<double,TDim> grad_d = prod(trans(DN_DXcontinuous),distances);
		grad_d /= norm_2(grad_d);

		double vn = grad_d[0]*vel[0];
		for(unsigned int i=1; i<TDim; i++) vn+=grad_d[i]*vel[i];
		
		//KRATOS_WATCH(medge_areas);

		//loop on edges
		unsigned int edge_counter = 0.0;
		for(unsigned int i=0; i<TDim; i++)
		{
			for(unsigned int j=i+1; j<TDim+1; j++)
			{
				const double aux = medge_areas[edge_counter]*vn;
				if( distances[i]*distances[j] < 0.0) //cut edge
				{           
					if(distances[i] > 0)
					{
						rRightHandSideVector[i] += aux;
						rRightHandSideVector[j] -= aux;
					}
					else
					{
						rRightHandSideVector[i] -= aux;
						rRightHandSideVector[j] += aux;
					}
				}
				edge_counter++;
			}
		}
		
		//FSI STRUCTURAL MASS TERM
		if( rCurrentProcessInfo.Has(DENSITY) )
                {
                  const double structural_density = rCurrentProcessInfo[DENSITY];
                  edge_counter = 0.0;
                  for(unsigned int i=0; i<TDim; i++)
                  {
                      for(unsigned int j=i+1; j<TDim+1; j++)
                      {
                              if( distances[i]*distances[j] < 0.0) //cut edge
                              {          
                                      //double li = fabs(distances[i]) + fabs(distances[j]);
                                      double Nj = 1.0; //distances[i]/li;
                                      double Ni = 1.0; // - Nj;
                                      
                                      const double aux = 1.0 / (structural_density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]) * medge_areas[edge_counter];
                                      rLeftHandSideMatrix(i,i)  += Ni * aux;
                                      rLeftHandSideMatrix(j,j)  += Nj * aux;                             
                              }
                              edge_counter++;
                      }
                  }                  
                  
                }

	}
}















template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::Calculate(const Variable<array_1d<double,3> > &rVariable,
                                     array_1d<double,3> &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == ADVPROJ)
    {
        double Tmp = 0.0;
        this->Calculate(DIVPROJ,Tmp,rCurrentProcessInfo);
    }
    else if (rVariable == CONV_PROJ)
    {
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const unsigned int LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType ConvTerm = ZeroVector(LocalSize);
        VectorType PresTerm = ZeroVector(LocalSize);
        VectorType DivTerm = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = GaussWeights[g];

            for (unsigned int i = 0; i < NumNodes; i++)
                NodalArea[i] += N[i] * GaussWeight;

            this->CalculateProjectionRHS(ConvTerm,PresTerm,DivTerm,N,DN_DX[g],GaussWeight);
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int RowIndex = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rConvVal = this->GetGeometry()[i].FastGetSolutionStepValue(CONV_PROJ);
            array_1d<double,3>& rPresVal = this->GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rConvVal[d] += ConvTerm[RowIndex];
                rPresVal[d] += PresTerm[RowIndex];
                ++RowIndex;
            }
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += DivTerm[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
    else if (rVariable == VELOCITY)
    {
        GeometryType& rGeom = this->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const SizeType LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType NodalVelCorrection = ZeroVector(LocalSize);
	
	

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; ++g)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

            double Density;

            this->EvaluateInPoint(Density,DENSITY,N);

            const double Coeff = GaussWeights[g] / ( Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0] );
	    
	    array_1d<double,TDim> DeltaPressureGradient(TDim,0.0);
	    this->EvaluateGradientInPoint(DeltaPressureGradient,PRESSURE_OLD_IT,rDN_DX);

            // Calculate contribution to the gradient term (RHS)
//             double DeltaPressure;
//             this->EvaluateInPoint(DeltaPressure,PRESSURE_OLD_IT,N);

            SizeType RowIndex = 0;

            for (SizeType i = 0; i < NumNodes; ++i)
            {
                for (SizeType d = 0; d < TDim; ++d)
                {
//                     NodalVelCorrection[RowIndex++] += Coeff * rDN_DX(i,d) * DeltaPressure;
		    NodalVelCorrection[RowIndex++] -= Coeff * N(i) * DeltaPressureGradient[d];
                }
            }
        }

        SizeType Index = 0;

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < TDim; ++d)
            {
                rTemp[d] += NodalVelCorrection[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }
}


template< unsigned int TDim >
void FractionalStepDiscontinuous<TDim>::Calculate(const Variable<double> &rVariable,
                                     double &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DIVPROJ)
    {
        const GeometryType& rGeom = this->GetGeometry();
        const SizeType NumNodes = rGeom.PointsNumber();
        const unsigned int LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        ShapeFunctionDerivativesArrayType DN_DX;
        Matrix NContainer;
        VectorType GaussWeights;
        this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
        const unsigned int NumGauss = GaussWeights.size();

        VectorType MomentumRHS = ZeroVector(LocalSize);
        VectorType MassRHS = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = GaussWeights[g];

            for (unsigned int i = 0; i < NumNodes; i++)
                NodalArea[i] += N[i] * GaussWeight;

            this->CalculateProjectionRHS(MomentumRHS,MassRHS,N,DN_DX[g],GaussWeight);
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int RowIndex = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rMomValue = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
            for (unsigned int d = 0; d < TDim; ++d)
                rMomValue[d] += MomentumRHS[RowIndex++];
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
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
	bool split_element = this->GetValue(SPLIT_ELEMENT);
/*	KRATOS_WATCH(this->Id())
	KRATOS_WATCH(this->GetValue(SPLIT_ELEMENT))*/
	if (split_element == true)
	{
		
		array_1d<double,(TDim+1)> distances = this->GetValue(ELEMENTAL_DISTANCES);
		unsigned int npos = 0, nneg=0;
		for(unsigned int i=0; i<TDim+1; i++)
			if(distances[i] >= 0) npos++;
			else nneg++;
		
		if(nneg !=0 && npos != 0)
		{
			boost::numeric::ublas::bounded_matrix<double,3*(TDim-1), (TDim+1)> Nenriched;
			array_1d<double,(3*(TDim-1))> volumes;
			boost::numeric::ublas::bounded_matrix<double,(TDim+1), TDim > coords;
			boost::numeric::ublas::bounded_matrix<double, 3*(TDim-1), (TDim+1) > Ngauss;
			array_1d<double,(3*(TDim-1))> signs;
			std::vector< Matrix > gauss_gradients(3*(TDim-1));

			boost::numeric::ublas::bounded_matrix<double, (TDim+1), TDim > DN_DX;
			array_1d<double, (TDim+1) > N;
				
			

			for (unsigned int i = 0; i < TDim+1; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				for (unsigned int j = 0; j < TDim; j++) coords(i, j) = xyz[j];
			}

			for (unsigned int i = 0; i < 3*(TDim-1); i++) gauss_gradients[i].resize(TDim+1, TDim, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions= DiscontinuousShapeFunctionsUtilities::CalculateDiscontinuousShapeFunctions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched,medge_areas);

			if(rGaussWeights.size() != ndivisions)
				rGaussWeights.resize(ndivisions,false);
			if(rNContainer.size1() != ndivisions || rNContainer.size2() != TDim+1)	
				rNContainer.resize(ndivisions,TDim+1,false);
			rDN_DX.resize(ndivisions);
			
//KRATOS_WATCH(this->Id())
			for (unsigned int g = 0; g < ndivisions; g++)
			{
				for(unsigned int d=0; d<TDim+1; d++)
				  rNContainer(g,d) = Nenriched(g,d);
				rGaussWeights[g] = volumes[g];
				rDN_DX[g] = gauss_gradients[g];			
/*			KRATOS_WATCH(g)
			KRATOS_WATCH(rDN_DX[g])
			KRATOS_WATCH(rGaussWeights[g])*/
			}
//			KRATOS_WATCH(rNContainer)
			

		}
		else
		{
			const GeometryType& rGeom = this->GetGeometry();
			Vector DetJ;
			rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_2);
			rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
			const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

			rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2),false);

			for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
				rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
		}
	}
	else
	{
		const GeometryType& rGeom = this->GetGeometry();
		Vector DetJ;
		rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_2);
		rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
		const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

		rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2),false);

		for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
			rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
	}
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
