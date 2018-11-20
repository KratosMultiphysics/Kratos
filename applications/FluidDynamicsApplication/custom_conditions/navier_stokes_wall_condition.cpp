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

#include "navier_stokes_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{


/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void NavierStokesWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
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
void NavierStokesWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
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



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);
    if (rLeftHandSideMatrix.size1() != MatrixSize)
        rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!
    if (rRightHandSideVector.size() != MatrixSize)
        rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

    // Struct to pass around the data
    ConditionDataStruct data;
    // Allocate memory needed
    array_1d<double,MatrixSize> rhs_gauss;
    BoundedMatrix<double,MatrixSize, MatrixSize> lhs_gauss;

    // LHS and RHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
    
    // Compute condition unit normal vector
    this->CalculateNormal(data.Normal); //this already contains the area
    const double A = norm_2(data.Normal);
    data.Normal /= A;
    
    // Store the outlet inflow prevention constants in the data structure
    data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
    const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
    data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];
    
    // Gauss point information
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    
    // Loop on gauss points
    for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    {
        data.N = row(Ncontainer, igauss);
        const double J = GaussPtsJDet[igauss];
        data.wGauss = J * IntegrationPoints[igauss].Weight();
        ComputeGaussPointRHSContribution(rhs_gauss, data);
        ComputeGaussPointLHSContribution(lhs_gauss, data);
        noalias(rLeftHandSideMatrix) += lhs_gauss;
        noalias(rRightHandSideVector) += rhs_gauss;
    }
    
    // Adding the BEHR2004 contribution if a slip BC is detected
    // Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663
    if (this->Is(SLIP)){

        Element parent;
        FindParent( parent );

        MatrixType rBehrSlipLeftHandSideMatrix = ZeroMatrix(MatrixSize,MatrixSize);
        VectorType rBehrSlipRightHandSideVector = ZeroVector(MatrixSize); 
        CalculateBehrSlipLeftHandSideContribution(rBehrSlipLeftHandSideMatrix, rCurrentProcessInfo, data, parent); 
        CalculateBehrSlipRightHandSideContribution(rBehrSlipRightHandSideVector, rCurrentProcessInfo, data, parent);
        noalias(rLeftHandSideMatrix) += rBehrSlipLeftHandSideMatrix;
        noalias(rRightHandSideVector) += rBehrSlipRightHandSideVector; 
    }

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

    if (rLeftHandSideMatrix.size1() != MatrixSize)
        rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!
    
    // LHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

    if (rRightHandSideVector.size() != MatrixSize)
        rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!
    
    // Struct to pass around the data
    ConditionDataStruct data;
    // Allocate memory needed
    array_1d<double,MatrixSize> rhs_gauss;
    // Loop on gauss points
    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
    
    // Compute condition normal
    this->CalculateNormal(data.Normal); //this already contains the area
    const double A = norm_2(data.Normal);
    data.Normal /= A;

    // Store the outlet inflow prevention constants in the data structure
    data.delta = 1e-2; // TODO: Decide if this constant should be fixed or not
    const ProcessInfo& rProcessInfo = rCurrentProcessInfo; // const to avoid race conditions on data_value_container access/initialization
    data.charVel = rProcessInfo[CHARACTERISTIC_VELOCITY];
    
    // Gauss point information
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    
    for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    {
        data.N = row(Ncontainer, igauss);
        const double J = GaussPtsJDet[igauss];
        data.wGauss = J * IntegrationPoints[igauss].Weight();
        ComputeGaussPointRHSContribution(rhs_gauss, data);
        noalias(rRightHandSideVector) += rhs_gauss;
    }

    KRATOS_CATCH("")
}


/// Condition check
/**
 * @param rCurrentProcessInfo reference to the ProcessInfo
 */
template<unsigned int TDim, unsigned int TNumNodes>
int NavierStokesWallCondition<TDim,TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0
    if (Check != 0) {
        return Check;
    }
    else {
        // Check that all required variables have been registered
        if(VELOCITY.Key() == 0)
            KRATOS_ERROR << "VELOCITY Key is 0. Check if the application was correctly registered.";
        if(MESH_VELOCITY.Key() == 0)
            KRATOS_ERROR << "MESH_VELOCITY Key is 0. Check if the application was correctly registered.";
        if(ACCELERATION.Key() == 0)
            KRATOS_ERROR << "ACCELERATION Key is 0. Check if the application was correctly registered.";
        if(PRESSURE.Key() == 0)
            KRATOS_ERROR << "PRESSURE Key is 0. Check if the application was correctly registered.";
        if(DENSITY.Key() == 0)
            KRATOS_ERROR << "DENSITY Key is 0. Check if the application was correctly registered.";
        if(DYNAMIC_VISCOSITY.Key() == 0)
            KRATOS_ERROR << "DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.";
        if(EXTERNAL_PRESSURE.Key() == 0)
            KRATOS_ERROR << "EXTERNAL_PRESSURE Key is 0. Check if the application was correctly registered.";
        
        // Checks on nodes
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR << "missing PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_ERROR << "missing MESH_VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_ERROR << "missing ACCELERATION variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                KRATOS_ERROR << "missing EXTERNAL_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_ERROR << "missing VELOCITY component degree of freedom on node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_ERROR << "missing PRESSURE component degree of freedom on node " << this->GetGeometry()[i].Id();
        }

        return Check;
    }

    KRATOS_CATCH("");
}   


/**
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void NavierStokesWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
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



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::FindParent(Element& rParent){

    KRATOS_TRY;

    if ( this->Is(SLIP) ){
    
        const array_1d<double,3>& rNormal = this->GetValue(NORMAL);

        if (norm_2(rNormal) == 0.0){
            KRATOS_ERROR << "NORMAL must be calculated before using this condition. Error on condition " << this->Id() << std::endl;
        }

        array_1d<double,3> Edge;
        GeometryType& rGeom = this->GetGeometry();
        WeakPointerVector<Element> ElementCandidates;
        
        for (SizeType i = 0; i < TNumNodes; i++){
            WeakPointerVector<Element>& rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);

            // Check that the condition has candidate parent elements
            KRATOS_ERROR_IF(rNodeElementCandidates.size() == 0) <<
                "Condition " << this->Id() << " has no candidate parent elements.\n" <<
                "Check that the FindNodalNeighboursProcess has been executed ( e.g in Solver.Initialize() )" << std::endl;

            for (SizeType j = 0; j < rNodeElementCandidates.size(); j++){

                ElementCandidates.push_back(rNodeElementCandidates(j));
            }
        }

        std::vector<IndexType> NodeIds(TNumNodes), ElementNodeIds;
        for (SizeType i=0; i < TNumNodes; i++){
            NodeIds[i] = rGeom[i].Id();
        }

        std::sort(NodeIds.begin(), NodeIds.end());

        for (SizeType i=0; i < ElementCandidates.size(); i++)
        {
            GeometryType& rElemGeom = ElementCandidates[i].GetGeometry();
            ElementNodeIds.resize(rElemGeom.PointsNumber());

            for (SizeType j=0; j < rElemGeom.PointsNumber(); j++){
                ElementNodeIds[j] = rElemGeom[j].Id();
            }

            std::sort(ElementNodeIds.begin(), ElementNodeIds.end());
            if ( std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()) )
            {
                rParent = ElementCandidates[i];
                return;
            }
        }
        KRATOS_ERROR << "PARENT NOT FOUND : error in condition -> " << this->Id() << std::endl;
    }

    KRATOS_CATCH("");
}



template <>
void NavierStokesWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
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




template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs_gauss,
const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    noalias(lhs_gauss) = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
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



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
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



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
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
void NavierStokesWallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.0;

}


/// Computes the 3D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void NavierStokesWallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
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
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateBehrSlipLeftHandSideContribution(  MatrixType& rLeftHandSideMatrix,
                                                                                            ProcessInfo& rCurrentProcessInfo,
                                                                                            ConditionDataStruct& rDataStruct,
                                                                                            Element& rParentElement )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();

    std::vector< Matrix > NodalNormals(TNumNodes);

    // Retrieving the nodal normal vectors and storing them in matrix format
	for (SizeType nnode = 0; nnode < TNumNodes; nnode++){

        NodalNormals[nnode].resize(3,1);

        double sum = 0.0;
        for (int j = 0; j < 3; j++){
            NodalNormals[nnode](j,0) = rGeom[nnode].FastGetSolutionStepValue(NORMAL)(j);
            sum += NodalNormals[nnode](j,0)*NodalNormals[nnode](j,0);
        }
        for (int j = 0; j < 3; j++){
            NodalNormals[nnode](j,0) /= sqrt(sum);
        }
    }

    MatrixType BaseLHSMatrix = zero_matrix<double>( TNumNodes*(TNumNodes+1) , TNumNodes*(TNumNodes+1) );
    MatrixType ProjectionLHSMatrix = zero_matrix<double>( TNumNodes*(TNumNodes+1) , TNumNodes*(TNumNodes+1) );
    MatrixType BaseLHSMatrixGPcontribution = zero_matrix<double>( TNumNodes*(TNumNodes+1) , TNumNodes*(TNumNodes+1) );

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);


    // Loop on gauss points
    for(unsigned int igauss = 0; igauss < NumGauss; igauss++) {
        
        BaseLHSMatrixGPcontribution = ZeroMatrix( (TNumNodes + 1)*TNumNodes , (TNumNodes + 1)*TNumNodes );

        array_1d<double, TNumNodes> N = row(Ncontainer, igauss);
        const double J = GaussPtsJDet[igauss];
        const double wGauss = J * IntegrationPoints[igauss].Weight();

        for (unsigned int lineBlock = 0; lineBlock < TNumNodes; lineBlock++){
            for(unsigned int rowBlock = 0; rowBlock < TNumNodes; rowBlock++){
                for(unsigned int i = 0; i < TNumNodes; i++){

                    BaseLHSMatrixGPcontribution( lineBlock * (TNumNodes + 1) + i , rowBlock * (TNumNodes + 1) + TNumNodes ) 
                    = rDataStruct.Normal[i] * N[lineBlock] * N[rowBlock] * wGauss;

                }
            }
        }
 
        BaseLHSMatrix += BaseLHSMatrixGPcontribution;
    }

    // Computation of NodalMultMatrix = ( [I] - (na)(na) ) for each node
    std::vector<MatrixType> NodalMultMatrix(TNumNodes);
    const MatrixType auxIdentMatrix = identity_matrix<double>(3);
    MatrixType auxMatrix = zero_matrix<double>(3,3);
    
    for (SizeType nnode=0; nnode < TNumNodes; nnode++){
        auxMatrix = prod( NodalNormals[nnode], trans(NodalNormals[nnode]) );
        NodalMultMatrix[nnode] = auxIdentMatrix - auxMatrix;
    }

    for(unsigned int nnode = 0; nnode < TNumNodes; nnode++){
        for( unsigned int i = 0; i < 3; i++){
            for( unsigned int j = 0; j < 3; j++){

                const unsigned int istart = nnode * (TNumNodes+1);
                const unsigned int jstart = nnode * (TNumNodes+1);
                ProjectionLHSMatrix(istart + i, jstart + j) = NodalMultMatrix[nnode](i,j);

            }
        }
    }

    rLeftHandSideMatrix = prod( ProjectionLHSMatrix, BaseLHSMatrixGPcontribution );

    KRATOS_CATCH("");
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateBehrSlipRightHandSideContribution( VectorType& rRightHandSideVector,
                                                                                            const ProcessInfo& rCurrentProcessInfo,
                                                                                            const ConditionDataStruct& rDataStruct,
                                                                                            Element& rParentElement )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();

    std::vector< IndexType > NodeIds(TNumNodes), ElementNodeIds;
    std::vector< Matrix > NodalNormals(TNumNodes);
    std::vector< Matrix > NodalEntriesRHS(TNumNodes);
    Matrix sigmaVoigtMatrix;
    Matrix stressTimesNormal;

    // Retrieve the nodal normal vectors, normalize, and store them in matrix
	for (unsigned int nnode=0; nnode < TNumNodes; nnode++){
		NodeIds[nnode] = rGeom[nnode].Id();
        NodalNormals[nnode].resize(3,1);

        double sumOfSquares = 0.0;
        for (unsigned int j = 0; j < 3; j++)
        {
            NodalNormals[nnode](j,0) = rGeom[nnode].FastGetSolutionStepValue(NORMAL)(j); 
            sumOfSquares += NodalNormals[nnode](j,0) * NodalNormals[nnode](j,0);
        }
        // normalization
        for (unsigned int j = 0; j < 3; j++)
        {
            NodalNormals[nnode](j,0) /= sqrt(sumOfSquares);
        }
	}

    // Computation of NodalMultMatrix = ( [I] - (na)(na) ) (for all nodes)
    std::vector<MatrixType> NodalProjectionMatrix(TNumNodes);
    const MatrixType auxIdentMatrix = identity_matrix<double>(3);
    
    for (unsigned int nnode=0; nnode < TNumNodes; nnode++){
        NodalProjectionMatrix[nnode] = auxIdentMatrix - prod( NodalNormals[nnode], trans(NodalNormals[nnode]) );
    }

    MatrixType normalVoigtMatrix = zero_matrix<double>(3,3);

    if ( TNumNodes == 2 ){
        if (normalVoigtMatrix.size2() != 3){
            normalVoigtMatrix.resize( 3, 3, false );
        }
        normalVoigtMatrix = zero_matrix<double>(3,3);
        sigmaVoigtMatrix = zero_matrix<double>(3,1);
        normalVoigtMatrix(0,0) = rDataStruct.Normal(0);
        normalVoigtMatrix(1,1) = rDataStruct.Normal(1);
        normalVoigtMatrix(0,2) = rDataStruct.Normal(1);
        normalVoigtMatrix(1,2) = rDataStruct.Normal(0);
    } else if ( TNumNodes == 3 ){
        if (normalVoigtMatrix.size2() != 6){
            normalVoigtMatrix.resize( 3, 6, false );
        }
        normalVoigtMatrix = zero_matrix<double>(3,6);
        sigmaVoigtMatrix = zero_matrix<double>(6,1);
        normalVoigtMatrix(0,0) = rDataStruct.Normal(0);
        normalVoigtMatrix(1,1) = rDataStruct.Normal(1);
        normalVoigtMatrix(2,2) = rDataStruct.Normal(2);
        normalVoigtMatrix(0,3) = rDataStruct.Normal(1);
        normalVoigtMatrix(0,5) = rDataStruct.Normal(2);
        normalVoigtMatrix(1,3) = rDataStruct.Normal(0);
        normalVoigtMatrix(1,4) = rDataStruct.Normal(2);
        normalVoigtMatrix(2,4) = rDataStruct.Normal(1);
        normalVoigtMatrix(2,5) = rDataStruct.Normal(0);
    } 

    VectorType ShearStressOfElement( (TDim-1)*3 , 0.0);
    ShearStressOfElement = ZeroVector( (TDim-1)*3 );
    
    rParentElement.Calculate( FLUID_STRESS, ShearStressOfElement, rCurrentProcessInfo );

    // adding PRESSURE to the viscous shear stresses
    std::vector< Matrix > CompleteNodalSigma(TNumNodes);

    if ( TNumNodes == 2 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            if ( CompleteNodalSigma[nnode].size1() != 3 ){
                CompleteNodalSigma[nnode].resize(3,1);
            }
            CompleteNodalSigma[nnode](0,0) = ShearStressOfElement[0] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](1,0) = ShearStressOfElement[1] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](2,0) = ShearStressOfElement[2]; // no pressure in shear component

#ifdef KRATOS_DEBUG
            if ( abs( ShearStressOfElement[0] ) > 0.001 || abs( ShearStressOfElement[1] ) > 0.001 ){
                KRATOS_WARNING("Behr Contribution in SLIP condition") << "The normal components of the viscous stress are still present" << std::endl;
            }
#endif
        }
    } else if ( TNumNodes == 3 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            if ( CompleteNodalSigma[nnode].size1() != 6 ){
                CompleteNodalSigma[nnode].resize(6,1);
            }
            CompleteNodalSigma[nnode](0,0) = ShearStressOfElement[0] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](1,0) = ShearStressOfElement[1] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](2,0) = ShearStressOfElement[2] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode](3,0) = ShearStressOfElement[3];  // no pressure in shear component
            CompleteNodalSigma[nnode](4,0) = ShearStressOfElement[4];  // no pressure in shear component
            CompleteNodalSigma[nnode](5,0) = ShearStressOfElement[5];  // no pressure in shear component

#ifdef KRATOS_DEBUG
            if ( abs( ShearStressOfElement[0] ) > 0.001 || abs( ShearStressOfElement[1] ) > 0.001 || abs( ShearStressOfElement[2] ) > 0.001 ){
                KRATOS_WARNING("Behr Contribution in SLIP condition") << "The normal components of the viscous stress are still present" << std::endl;
            }
#endif
        }
    }

    // retrieving Gauss integration point data
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    Matrix CompleteSigmaInterpolated;
    CompleteSigmaInterpolated = ZeroMatrix(3,1);

    // Loop on nodes
    for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){

        if (NodalEntriesRHS[nnode].size1() != TNumNodes+1){
            NodalEntriesRHS[nnode].resize(TNumNodes+1, 1);
        }

        NodalEntriesRHS[nnode] = zero_matrix<double>(3, 1);
        
        // Loop on gauss points
        for(unsigned int igauss = 0; igauss < NumGauss; igauss++){

            array_1d<double, TNumNodes> N = row(Ncontainer, igauss);
            const double J = GaussPtsJDet[igauss];
            const double wGauss = J * IntegrationPoints[igauss].Weight();

            CompleteSigmaInterpolated = ZeroMatrix(3,1);

            for( unsigned int comp = 0; comp < TNumNodes; comp++){

                CompleteSigmaInterpolated += N[comp] * prod( normalVoigtMatrix, CompleteNodalSigma[comp] );
            }
            
            NodalEntriesRHS[nnode] += ( wGauss * N(nnode) * CompleteSigmaInterpolated );
        }

        NodalEntriesRHS[nnode] = prod( NodalProjectionMatrix[nnode], NodalEntriesRHS[nnode] );
    }

    for (unsigned int entry = 0; entry < 3; entry++){

        for (unsigned int node = 0; node < TNumNodes; node++){
            rRightHandSideVector(entry + node*(TNumNodes+1) ) = NodalEntriesRHS[node](entry,0);
        }
    }

    KRATOS_CATCH("")
}


template class NavierStokesWallCondition<2,2>;
template class NavierStokesWallCondition<3,3>;

} // namespace Kratos
