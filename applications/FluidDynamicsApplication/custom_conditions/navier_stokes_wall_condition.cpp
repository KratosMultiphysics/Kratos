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
#include "includes/checks.h"

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

    if ( this->Is(SLIP) ){
        // finding parent element to retrieve viscous stresses which are later stored in "data"
        WeakPointerVector<Element> parentElement = this->GetValue( NEIGHBOUR_ELEMENTS );
        KRATOS_ERROR_IF( parentElement.size() > 1 ) << "A condition was assigned more than one parent element." << std::endl;
        KRATOS_ERROR_IF( parentElement.size() == 0 ) << "A condition was NOT assigned a parent element. "
        << "This leads to errors for the slip condition [BEHR2004] "
        << "Please execute the check_and_prepare_model_process_fluid process." << std::endl;

        Element& parent = parentElement[0];
        data.ViscousStress = ZeroVector( 3*(TDim-1) );
        parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    }

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

    if ( this->Is(SLIP) ){
        // finding parent element to retrieve viscous stresses which are later stored in "data"
        WeakPointerVector<Element> parentElement = this->GetValue( NEIGHBOUR_ELEMENTS );
        KRATOS_ERROR_IF( parentElement.size() > 1 ) << "A condition was assigned more than one parent element." << std::endl;
        KRATOS_ERROR_IF( parentElement.size() == 0 ) << "A condition was NOT assigned a parent element. "
        << "This leads to errors for the slip condition [BEHR2004] "
        << "Please execute the check_and_prepare_model_process_fluid process." << std::endl;

        Element& parent = parentElement[0];
        data.ViscousStress = ZeroVector( 3*(TDim-1) );
        parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    }

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
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(MESH_VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
        KRATOS_CHECK_VARIABLE_KEY(PRESSURE)
        KRATOS_CHECK_VARIABLE_KEY(DENSITY)
        KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY)
        KRATOS_CHECK_VARIABLE_KEY(EXTERNAL_PRESSURE)

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
    lhs_gauss = ZeroMatrix(TNumNodes*LocalSize, TNumNodes*LocalSize);

    // contribution to avoid tangential components in the residual (BEHR2004)
    // Adding the BEHR2004 contribution if a slip BC is detected
    // Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663
    if (this->Is(SLIP)){

        Matrix lhs_gauss_behr = ZeroMatrix(TNumNodes*LocalSize,TNumNodes*LocalSize);
        ComputeGaussPointBehrSlipLHSContribution( lhs_gauss_behr, data );
        lhs_gauss += lhs_gauss_behr;
    }
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
    if (this->Is(OUTLET)){
        this->ComputeRHSOutletInflowContribution(rhs_gauss, data);
    }

    // contribution to avoid tangential components in the residual (BEHR2004)
    // Adding the BEHR2004 contribution if a slip BC is detected
    // Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663
    if (this->Is(SLIP)){

        Vector rhs_gauss_behr = ZeroVector(TNumNodes*LocalSize);
        ComputeGaussPointBehrSlipRHSContribution( rhs_gauss_behr, data );
        rhs_gauss += rhs_gauss_behr;
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
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointBehrSlipLHSContribution(  Matrix& rLeftHandSideMatrix,
                                                                                           const ConditionDataStruct& rDataStruct )
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();

    // Retrieve the nodal consistent normal vectors, normalize, and store them for each node
    std::vector<array_1d<double,3>> NodalNormals(TNumNodes);
	for (unsigned int nnode=0; nnode < TNumNodes; nnode++){

        NodalNormals[nnode] = rGeom[nnode].FastGetSolutionStepValue(NORMAL);
        double sumOfSquares = 0.0;
        for (unsigned int j = 0; j < 3; j++){
            sumOfSquares += NodalNormals[nnode][j] * NodalNormals[nnode][j];
        }
        NodalNormals[nnode] = ( 1.0 / sqrt(sumOfSquares) ) *  NodalNormals[nnode];
	}

    MatrixType BaseLHSMatrix = zero_matrix<double>( TNumNodes*(TNumNodes+1) , TNumNodes*(TNumNodes+1) );
    MatrixType ProjectionLHSMatrix = zero_matrix<double>( TNumNodes*(TNumNodes+1) , TNumNodes*(TNumNodes+1) );

    // Loop on gauss points removed (!!!)
    BaseLHSMatrix = ZeroMatrix( (TNumNodes + 1)*TNumNodes , (TNumNodes + 1)*TNumNodes );

    array_1d<double, TNumNodes> N = rDataStruct.N;
    const double wGauss = rDataStruct.wGauss;

    for (unsigned int lineBlock = 0; lineBlock < TNumNodes; lineBlock++){
        for(unsigned int rowBlock = 0; rowBlock < TNumNodes; rowBlock++){
            for(unsigned int i = 0; i < TNumNodes; i++){

                BaseLHSMatrix( lineBlock * (TNumNodes + 1) + i , rowBlock * (TNumNodes + 1) + TNumNodes )
                = rDataStruct.Normal[i] * N[lineBlock] * N[rowBlock] * wGauss;
            }
        }
    }

    // Computation of NodalProjectionMatrix = ( [I] - (na)(na) ) for all nodes and store it
    std::vector< BoundedMatrix<double, 3, 3> > NodalProjectionMatrix(TNumNodes);
    for(unsigned int node = 0; node < TNumNodes; node++){
        FluidElementUtilities<3>::SetTangentialProjectionMatrix( NodalNormals[node], NodalProjectionMatrix[node] );
    }

    for(unsigned int nnode = 0; nnode < TNumNodes; nnode++){
        for( unsigned int i = 0; i < 3; i++){
            for( unsigned int j = 0; j < 3; j++){

                const unsigned int istart = nnode * (TNumNodes+1);
                const unsigned int jstart = nnode * (TNumNodes+1);
                ProjectionLHSMatrix(istart + i, jstart + j) = NodalProjectionMatrix[nnode](i,j);
            }
        }
    }
    rLeftHandSideMatrix = prod( ProjectionLHSMatrix, BaseLHSMatrix );

    KRATOS_CATCH("");
}




template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointBehrSlipRHSContribution(   VectorType& rRightHandSideVector,
                                                                                            const ConditionDataStruct& rDataStruct )
{
    KRATOS_TRY
    const unsigned int voigtSize = 3 * (TDim-1);

    GeometryType& rGeom = this->GetGeometry();

    // Retrieve the nodal consistent normal vectors, normalize, and store them for each node
    std::vector<array_1d<double,3>> NodalNormals(TNumNodes);
	for (unsigned int nnode=0; nnode < TNumNodes; nnode++){

        NodalNormals[nnode] = rGeom[nnode].FastGetSolutionStepValue(NORMAL);
        double sumOfSquares = 0.0;
        for (unsigned int j = 0; j < 3; j++){
            sumOfSquares += NodalNormals[nnode][j] * NodalNormals[nnode][j];
        }

        NodalNormals[nnode] = ( 1.0 / sqrt(sumOfSquares) ) *  NodalNormals[nnode];
	}

    // Computation of NodalProjectionMatrix = ( [I] - (na)(na) ) for all nodes and store it
    std::vector< BoundedMatrix<double, 3, 3> > NodalProjectionMatrix(TNumNodes);
    for(unsigned int node = 0; node < TNumNodes; node++){
        FluidElementUtilities<3>::SetTangentialProjectionMatrix( NodalNormals[node], NodalProjectionMatrix[node] );
    }

    // Computation of a matrix to replace [S] n by [A(n)] [Svoigt]
    BoundedMatrix<double, TNumNodes, voigtSize> conditionNormalForVoigt = ZeroMatrix( TNumNodes, voigtSize );
    FluidElementUtilities<3>::VoigtTransformForProduct( rDataStruct.Normal, conditionNormalForVoigt );


    // Computing the full stress for the nodes (still in Voigt notation)
    Vector ShearStressOfElement( voigtSize, 0.0);

    // step 1: Retrieving viscous pressure from the element (constant in element assumed)
    ShearStressOfElement = rDataStruct.ViscousStress;


    // step 2: adding pressure (different at nodes) to the viscous shear stresses
    std::vector< array_1d< double, voigtSize > > CompleteNodalSigma(TNumNodes);

    if ( TNumNodes == 2 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            CompleteNodalSigma[nnode] = ZeroVector( voigtSize );
            CompleteNodalSigma[nnode][0] = ShearStressOfElement[0] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode][1] = ShearStressOfElement[1] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode][2] = ShearStressOfElement[2]; // no pressure in shear component
#ifdef KRATOS_DEBUG
            if ( abs( ShearStressOfElement[0] ) > 0.001 || abs( ShearStressOfElement[1] ) > 0.001 ){
                KRATOS_WARNING("Behr Contribution in SLIP condition") << "The normal components of the viscous stress are still present" << std::endl;
            }
#endif
        }
    } else if ( TNumNodes == 3 ){
        for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){
            CompleteNodalSigma[nnode] = ZeroVector( voigtSize );
            CompleteNodalSigma[nnode][0] = ShearStressOfElement[0] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode][1] = ShearStressOfElement[1] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode][2] = ShearStressOfElement[2] - rGeom[nnode].FastGetSolutionStepValue(PRESSURE);
            CompleteNodalSigma[nnode][3] = ShearStressOfElement[3];  // no pressure in shear component
            CompleteNodalSigma[nnode][4] = ShearStressOfElement[4];  // no pressure in shear component
            CompleteNodalSigma[nnode][5] = ShearStressOfElement[5];  // no pressure in shear component
#ifdef KRATOS_DEBUG
            if ( abs( ShearStressOfElement[0] ) > 0.001 || abs( ShearStressOfElement[1] ) > 0.001 || abs( ShearStressOfElement[2] ) > 0.001 ){
                KRATOS_WARNING("Behr Contribution in SLIP condition") << "The normal components of the viscous stress are still present" << std::endl;
            }
#endif
        }
    }

    Vector CompleteSigmaInterpolated;
    CompleteSigmaInterpolated = ZeroVector(3);

    std::vector<array_1d<double,TNumNodes+1>> NodalEntriesRHS(TNumNodes);

    // Loop all nodal contributions
    for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){

        NodalEntriesRHS[nnode] = zero_vector<double>(3);
        const array_1d<double, TNumNodes> N = rDataStruct.N;
        const double wGauss = rDataStruct.wGauss;

        CompleteSigmaInterpolated = ZeroVector(3);
        for( unsigned int comp = 0; comp < TNumNodes; comp++){

            CompleteSigmaInterpolated += N[comp] * prod( conditionNormalForVoigt, CompleteNodalSigma[comp] );
        }

        NodalEntriesRHS[nnode] = ( wGauss * N(nnode) * CompleteSigmaInterpolated );
        NodalEntriesRHS[nnode] = prod( NodalProjectionMatrix[nnode], NodalEntriesRHS[nnode] );
    }

    for (unsigned int node = 0; node < TNumNodes; node++){
        for (unsigned int entry = 0; entry < TNumNodes; entry++){
            rRightHandSideVector( node*(TNumNodes+1) + entry ) = NodalEntriesRHS[node][entry];
        }
    }

    KRATOS_CATCH("")
}

template class NavierStokesWallCondition<2,2>;
template class NavierStokesWallCondition<3,3>;

} // namespace Kratos
