//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#include "behr_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void BehrWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void BehrWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void BehrWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void BehrWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/// Computes the Gauss pt. LHS contribution
/**
* @param lhs_gauss reference to the local LHS matrix
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs_gauss,
const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    noalias(lhs_gauss) = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);
}

/// Computes the Gauss pt. RHS contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
const ConditionDataStruct& data)
{
    // Initialize the local RHS
    const unsigned int LocalSize = TDim+1;
    noalias(rhs_gauss) = ZeroVector(TNumNodes*LocalSize);

    // Gauss pt. Neumann BC contribution
    this->ComputeRHSNeumannContribution(rhs_gauss, data);

    // Gauss pt. outlet inflow prevention contribution
    if (this->Is(OUTLET))
    {
        this->ComputeRHSOutletInflowContribution(rhs_gauss, data);
    }
}

/// Computes the condition RHS Neumann BC contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
                                                                              const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Add Neumann BC contribution
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const double pext = rGeom[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

        for (unsigned int j=0; j<TNumNodes; ++j)
        {
            unsigned int row = j*LocalSize;
            for (unsigned int d=0; d<TDim; ++d)
            {
                rhs_gauss[row+d] -= data.wGauss*data.N[j]*data.N[i]*pext*data.Normal[d];
            }
        }
    }
}

/// Computes the condition RHS outlet inflow prevention contribution
/**
* @param rhs_gauss reference to the local RHS vector
* @param data Gauss pt. data structure
*/
template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
                                                                                   const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Compute Gauss pt. density, velocity norm and velocity projection
    double rhoGauss = 0.0;
    array_1d<double, 3> vGauss = ZeroVector(3);
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const double& rRho = rGeom[i].FastGetSolutionStepValue(DENSITY);
        const array_1d<double, 3>& rVelNode = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        rhoGauss += data.N[i]*rRho;
        vGauss += data.N[i]*rVelNode;
    }

    const double vGaussProj = inner_prod(vGauss, data.Normal);
    const double vGaussSquaredNorm = std::pow(vGauss[0],2) + std::pow(vGauss[1],2) + std::pow(vGauss[2],2);

    // Add outlet inflow prevention contribution
    const double delta = data.delta;
    const double U_0 = data.charVel;
    const double S_0 = 0.5*(1-tanh(vGaussProj/(U_0*delta)));

    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        unsigned int row = i*LocalSize;
        for (unsigned int d=0; d<TDim; ++d)
        {
            rhs_gauss[row+d] += data.wGauss*data.N[i]*0.5*rhoGauss*vGaussSquaredNorm*S_0*data.Normal[d];
        }
    }
}


/// Computes the 2D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void BehrWallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}


/// Computes the 3D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void BehrWallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::Initialize()
{
	KRATOS_TRY;
	const array_1d<double,3>& rNormal = this->GetValue(NORMAL);
     // std::cout << rNormal << std::endl;
	if (norm_2(rNormal) == 0.0){
	    std::cout << "error on condition -> " << this->Id() << std::endl;
	    KRATOS_THROW_ERROR(std::logic_error, "NORMAL must be calculated before using this condition","");
	}
	if (mInitializeWasPerformed){
        KRATOS_ERROR_IF(mpParentElement == NULL) << "Parent as NULL pointer (Check 1)";
	    return;
	}

	mInitializeWasPerformed = true;
	double EdgeLength;
	array_1d<double,3> Edge;
	GeometryType& rGeom = this->GetGeometry();
	WeakPointerVector<Element> ElementCandidates;
    
    // std::cout << "before finding candidates" << std::flush;
	for (SizeType i = 0; i < TNumNodes; i++){
         // std::cout << "Looking for condidates for node " << i << " ----- " << std::flush;
		WeakPointerVector<Element>& rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);

         // Check that the condition has candidate parent elements
         KRATOS_ERROR_IF(rNodeElementCandidates.size() == 0) <<
             "Condition " << this->Id() << " has no candidate parent elements.\n" <<
             "Check that the FindNodalNeighboursProcess has been executed ( e.g in Solver.Initialize() )";
         // std::cout << "The size of candidates in " << rNodeElementCandidates.size() << std::flush;

		for (SizeType j = 0; j < rNodeElementCandidates.size(); j++){

            ElementCandidates.push_back(rNodeElementCandidates(j));
            // std::cout << "1 candidate found" << std::flush;
		}
	}

	std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
	for (SizeType i=0; i < TNumNodes; i++){
		NodeIds[i] = rGeom[i].Id();
         // std::cout << rGeom[i].Id() << std::flush;
	}

	std::sort(NodeIds.begin(), NodeIds.end());

	for (SizeType i=0; i < ElementCandidates.size(); i++)
	{
		GeometryType& rElemGeom = ElementCandidates[i].GetGeometry();
		ElementNodeIds.resize(rElemGeom.PointsNumber());

		for (SizeType j=0; j < rElemGeom.PointsNumber(); j++){
			ElementNodeIds[j] = rElemGeom[j].Id();
             // std::cout << rElemGeom[j].Id() << std::flush;
		}

		std::sort(ElementNodeIds.begin(), ElementNodeIds.end());
		if ( std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()) )
		{
            // std::cout << "Parent was found " << i << std::endl;
			mpParentElement = &ElementCandidates[i];

            std::cout << mpParentElement->Info() << std::endl;   /// works
            
            // // find the additional point
            // IndexType myAdditionalNodeId = 10;
            // for (auto iter = ElementNodeIds.begin(); iter != ElementNodeIds.end(); iter++){
            //     if ( std::find(NodeIds.begin(), NodeIds.end(), *iter) == NodeIds.end() ){
            //         myAdditionalNodeId = *iter;
            //         break;
            //     } 
            // }

            // if (myAdditionalNodeId != 10){
            //     for (SizeType i = 0; i < rElemGeom.PointsNumber(); i++){
            //         if ( rElemGeom[i].Id() == myAdditionalNodeId ){
            //             mAdditionalNodeFromParent = i;
            //             std::cout << "Additional node from parent element " << mAdditionalNodeFromParent << " -- " << std::endl;
            //             break;
            //         }
            //     }
            // } else {
            //     std::cout << "The node of the parent element was not found" << std::endl;
            //     KRATOS_THROW_ERROR(std::logic_error, "Condition cannot find complete parent element","");
            // }

			// Edge = rElemGeom[1].Coordinates() - rElemGeom[0].Coordinates();
			// mMinEdgeLength = Edge[0]*Edge[0];
			// for (SizeType d=1; d < TDim; d++){
			// 	mMinEdgeLength += Edge[d]*Edge[d];
			// }
			// for (SizeType j=2; j < rElemGeom.PointsNumber(); j++){
			// 	for(SizeType k=0; k < j; k++){
			// 		Edge = rElemGeom[j].Coordinates() - rElemGeom[k].Coordinates();
			// 		EdgeLength = Edge[0]*Edge[0];
			// 		for (SizeType d = 1; d < TDim; d++){
			// 			EdgeLength += Edge[d]*Edge[d];
			// 		}
			// 		mMinEdgeLength = (EdgeLength < mMinEdgeLength) ? EdgeLength : mMinEdgeLength;
			// 	}
			// }
			// mMinEdgeLength = sqrt(mMinEdgeLength);
            //  // Leave immediately one the parent was found

            KRATOS_ERROR_IF(mpParentElement == NULL) << "Parent as NULL pointer (Check2)";
			return;
		}
	}
	std::cout << "PARENT NOT FOUND : error in condition -> " << this->Id() << std::endl;
	KRATOS_THROW_ERROR(std::logic_error, "Condition cannot find parent element","");
	KRATOS_CATCH("");
}




/// Calculates the LHS condition contributions                                           ----- TO DO -----------
/**
 * Clones the selected element variables, creating a new one
 * @param rLeftHandSideMatrix reference to the LHS matrix
 * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
 */

template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::CalculateLeftHandSide(  MatrixType& rLeftHandSideMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY
    
    // GeometryType& rGeom = this->GetGeometry();     // one edge or one face

    // std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
    // std::vector< Matrix > NodalNormals(TNumNodes);

    // // Retrieving the nodal normal vectors and storing them in matrix format
	// for (SizeType i=0; i < TNumNodes; i++){
	// 	NodeIds[i] = rGeom[i].Id();
    //     NodalNormals[i].resize(3,1);

    //     double sum = 0.0;
    //     for (int j = 0; j < 3; j++)
    //     {
    //         NodalNormals[i](j,0) = rGeom[i].FastGetSolutionStepValue(NORMAL)(j);
    //         sum += NodalNormals[i](j,0)*NodalNormals[i](j,0);
    //     }
    //     // normalization
    //     for (int j = 0; j < 3; j++)
    //     {
    //         NodalNormals[i](j,0) /= sqrt(sum);
    //     }

	// }

    // // KRATOS_WATCH( NodalNormals )

    // // Computation of NodalMultMatrix = ( [I] - (na)(na) )
    // std::vector<MatrixType> NodalMultMatrix(TNumNodes);
    // const MatrixType auxIdentMatrix = identity_matrix<double>(3);
    // MatrixType auxMatrix = zero_matrix<double>(3,3);
    
    // for (SizeType i=0; i < TNumNodes; i++){
    //     auxMatrix = prod( NodalNormals[i], trans(NodalNormals[i]) );
    //     NodalMultMatrix[i] = auxIdentMatrix - auxMatrix;
    // }

    // MatrixType MultTimesNormal = zero_matrix<double>(3,1);
    // MatrixType BaseLHSMatrix = zero_matrix<double>(rLeftHandSideMatrix.size1(),rLeftHandSideMatrix.size2());


    // for (unsigned int nnodes = 0; nnodes < TNumNodes; nnodes++){
    //     // assignment of the normal vector for the current node
    //     for (unsigned int i = 0; i < 3; i++){
    //         MultTimesNormal(i,0) = data.Normal(i);
    //     }

    //     // Multiplication with the MultMatrix of the current node
    //     MultTimesNormal = prod( NodalMultMatrix[nnodes], MultTimesNormal );

    //     // KRATOS_WATCH( MultTimesNormal )

    //     if (nnodes == 0){
    //         for (unsigned int row = 0; row < 3; row++){
    //             BaseLHSMatrix(row,2) = MultTimesNormal(row,0);
    //         } 
    //     }
    //     else if (nnodes == 1){
    //         for (unsigned int row = 3; row < 6; row++){
    //             BaseLHSMatrix(row,5) = MultTimesNormal(row,0);
    //         } 
    //     }
    // }

    // // KRATOS_WATCH( BaseLHSMatrix )

    // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    // const unsigned int NumGauss = IntegrationPoints.size();
    // Vector GaussPtsJDet = ZeroVector(NumGauss);
    // rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    
    // const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    // MatrixType Nvert = zero_matrix<double> (6,1);
    // MatrixType Nhori = zero_matrix<double> (1,6);

    // const unsigned int LocalSize = TDim+1;
    // MatrixType lhs_gauss = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);

    // noalias(rLeftHandSideMatrix) = zero_matrix<double>(rLeftHandSideMatrix.size1(),rLeftHandSideMatrix.size2());

    // // Loop on gauss points
    // for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    // {
    //     lhs_gauss = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);

    //     data.N = row(Ncontainer, igauss);
    //     const double J = GaussPtsJDet[igauss];
    //     const double wGauss = J * IntegrationPoints[igauss].Weight();

    //     Nvert(2,0) = data.N[0];
    //     Nvert(5,0) = data.N[1];

    //     Nhori(0,2) = data.N[0];
    //     Nhori(0,5) = data.N[1];

    //     lhs_gauss = prod( Nvert, Nhori );

    //     lhs_gauss = prod( BaseLHSMatrix, lhs_gauss );

    //     noalias(rLeftHandSideMatrix) += wGauss * lhs_gauss * 0.0;
    // }

    // // const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    // // std::cout << "N = " << N << std::endl;
    // // this->GetGeometry().ShapeFunctionsLocalGradients
    // // const ShapeFunctionDerivativesType& DN;

    // // std::cout << "use count of weak pointer " << mpParentElement.use_count() << std::flush;
    // // ElementPointerType mParElement = mpParentElement.lock();
    // // auto rParGeom = mpParentElement->pGetGeometry();

    // // ShapeFunctionsGradientsType DNelem;

    // // rParGeom->ShapeFunctionsIntegrationPointsGradients(DNelem, GeometryData::GI_GAUSS_2);
    // // std::cout << "DNelem = " << DNelem << std::endl;

    // // KRATOS_WATCH( rLeftHandSideMatrix )

    // KRATOS_CATCH("")
}





/// Calculates the RHS condition contributions
/**
 * Clones the selected element variables, creating a new one
 * @param rRightHandSideVector reference to the RHS matrix
 * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
 */

template<unsigned int TDim, unsigned int TNumNodes>
void BehrWallCondition<TDim,TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();

    std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
    std::vector< Matrix > NodalNormals(TNumNodes);
    std::vector< Matrix > NodalEntriesRHS(TNumNodes);
    Matrix sigmaVoigtMatrix;
    Matrix stressTimesNormal;

    // Retrieve the nodal normal vectors, normalize, and store them in matrix
	for (SizeType nnode=0; nnode < TNumNodes; nnode++){
		NodeIds[nnode] = rGeom[nnode].Id();
        NodalNormals[nnode].resize(3,1);

        double sumOfSquares = 0.0;
        for (int j = 0; j < 3; j++)
        {
            NodalNormals[nnode](j,0) = rGeom[nnode].FastGetSolutionStepValue(NORMAL)(j);
            sumOfSquares += NodalNormals[nnode](j,0) * NodalNormals[nnode](j,0);
        }
        // normalization
        for (int j = 0; j < 3; j++)
        {
            NodalNormals[nnode](j,0) /= sqrt(sumOfSquares);
        }
	}

    // KRATOS_WATCH( NodalNormals )


    // Computation of NodalMultMatrix = ( [I] - (na)(na) ) (for all nodes)
    std::vector<MatrixType> NodalProjectionMatrix(TNumNodes);
    const MatrixType auxIdentMatrix = identity_matrix<double>(3);
    MatrixType auxMatrix = zero_matrix<double>(3,3);
    
    for (SizeType nnode=0; nnode < TNumNodes; nnode++){
        auxMatrix = prod( NodalNormals[nnode], trans(NodalNormals[nnode]) );
        NodalProjectionMatrix[nnode] = auxIdentMatrix - auxMatrix;
    }

    // KRATOS_WATCH( NodalProjectionMatrix )


    // FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal, voigt_normal_projection_matrix)
    MatrixType normalVoigtMatrix = zero_matrix<double>(3,3);

    if ( TNumNodes == 2 ){
        if (normalVoigtMatrix.size2() != 3){
            normalVoigtMatrix.resize( 3, 3, false );
        }
        normalVoigtMatrix = zero_matrix<double>(3,3);
        sigmaVoigtMatrix = zero_matrix<double>(3,1);
        normalVoigtMatrix(0,0) = data.Normal(0);
        normalVoigtMatrix(1,1) = data.Normal(1);
        normalVoigtMatrix(0,2) = data.Normal(1);
        normalVoigtMatrix(1,2) = data.Normal(0);
    } else if ( TNumNodes == 3 ){
        if (normalVoigtMatrix.size2() != 6){
            normalVoigtMatrix.resize( 3, 6, false );
        }
        normalVoigtMatrix = zero_matrix<double>(3,6);
        sigmaVoigtMatrix = zero_matrix<double>(6,1);
        normalVoigtMatrix(0,0) = data.Normal(0);
        normalVoigtMatrix(1,1) = data.Normal(1);
        normalVoigtMatrix(2,2) = data.Normal(2);
        normalVoigtMatrix(0,3) = data.Normal(1);
        normalVoigtMatrix(0,5) = data.Normal(2);
        normalVoigtMatrix(1,3) = data.Normal(0);
        normalVoigtMatrix(1,4) = data.Normal(2);
        normalVoigtMatrix(2,4) = data.Normal(1);
        normalVoigtMatrix(2,5) = data.Normal(0);
    } 

    // KRATOS_WATCH( normalVoigtMatrix )

    // checking for existance of the parent
    if (mpParentElement != NULL){
        std::cout << mpParentElement->Info() << std::endl;    ////   PROBLEM
    } else {
        std::cout << "NULL was found..." << std::endl;
    }
    
    // linearization in terms of the ShearStress ( not considered as f(u^{i+1} )
    VectorType ShearStressOfElement(3, 0.0);
    mpParentElement->Calculate( FLUID_STRESS, ShearStressOfElement, rCurrentProcessInfo );

    // adding PRESSURE to the viscous shear stresses
    std::vector< Matrix > CompleteNodalSigma(TNumNodes);

    if ( TNumNodes == 2 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            if ( CompleteNodalSigma[nnode].size1() != 3 ){
                CompleteNodalSigma[nnode].resize(3,1);
            }
            CompleteNodalSigma[nnode](0,0) = ShearStressOfElement[0] + rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](1,0) = ShearStressOfElement[1] + rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](2,0) = ShearStressOfElement[2]; // no pressure in shear component
        }
    } else if ( TNumNodes == 3 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            if ( CompleteNodalSigma[nnode].size1() != 6 ){
                CompleteNodalSigma[nnode].resize(6,1);
            }
            CompleteNodalSigma[nnode](0,0) = ShearStressOfElement[0] + rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](1,0) = ShearStressOfElement[1] + rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](2,0) = ShearStressOfElement[2] + rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](3,0) = ShearStressOfElement[3]; // no pressure in shear component
            CompleteNodalSigma[nnode](4,0) = ShearStressOfElement[4]; // no pressure in shear component
            CompleteNodalSigma[nnode](5,0) = ShearStressOfElement[5]; // no pressure in shear component
        }
    }

    // KRATOS_WATCH( CompleteNodalSigma )

    // retrieving Gauss integration point data
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    Matrix CompleteSigmaInterpolated;
    CompleteSigmaInterpolated.resize(3,1);

    // Loop on nodes
    for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){

        if (NodalEntriesRHS[nnode].size1() != TNumNodes+1){
            NodalEntriesRHS[nnode].resize(TNumNodes+1, 1);
        }

        NodalEntriesRHS[nnode] = zero_matrix<double>(TNumNodes+1, 1);
        // Loop on gauss points
        for(unsigned int igauss = 0; igauss<NumGauss; igauss++){

            data.N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            data.wGauss = J * IntegrationPoints[igauss].Weight();

            // interpolation
            if (TNumNodes == 2){
                CompleteSigmaInterpolated = data.N[0] * prod( normalVoigtMatrix, CompleteNodalSigma[0] ) + 
                                            data.N[1] * prod( normalVoigtMatrix, CompleteNodalSigma[1] );
            } else if (TNumNodes == 3){
                CompleteSigmaInterpolated = data.N[0] * prod( normalVoigtMatrix, CompleteNodalSigma[0] ) +
                                            data.N[1] * prod( normalVoigtMatrix, CompleteNodalSigma[1] ) +
                                            data.N[2] * prod( normalVoigtMatrix, CompleteNodalSigma[2] );
            }
            
            NodalEntriesRHS[nnode] += ( data.wGauss * data.N(nnode) * CompleteSigmaInterpolated );
        }
        NodalEntriesRHS[nnode] = prod( NodalProjectionMatrix[nnode], NodalEntriesRHS[nnode] );
    }

    for (unsigned int entry = 0; entry < TNumNodes+1; entry++){

        for (unsigned int node = 0; node < TNumNodes; node++){
            rRightHandSideVector(entry + node*(TNumNodes+1) ) = -NodalEntriesRHS[node](entry,0);
        }
    }

    // KRATOS_WATCH( rRightHandSideVector )

    KRATOS_CATCH("")
}


template class BehrWallCondition<2,2>;
template class BehrWallCondition<3,3>;

} // namespace Kratos
