#include "swimming_DEM_application.h"
#include "calculate_fluid_fraction_gradient_element.h"

namespace Kratos
{


template<unsigned int TDim, unsigned int TNumNodes>
Element::Pointer ComputeFluidFractionGradient<TDim, TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ComputeFluidFractionGradient>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template<unsigned int TDim, unsigned int TNumNodes>
Element::Pointer ComputeFluidFractionGradient<TDim, TNumNodes>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ComputeFluidFractionGradient>(NewId, pGeom, pProperties);
}



/// Calculate the element's local contribution to the system for the current step.
template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = TDim * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) = 0.0;
        }
        rRightHandSideVector(i) = 0.0;
    }

    CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const {

    const unsigned int LocalSize = TDim * TNumNodes;
    unsigned int LocalIndex = 0;
    unsigned int local_position = this->GetGeometry()[0].GetDofPosition(FLUID_FRACTION_GRADIENT_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(FLUID_FRACTION_GRADIENT_X, local_position).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(FLUID_FRACTION_GRADIENT_Y, local_position + 1).EquationId();
        if (TDim == 3) rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(FLUID_FRACTION_GRADIENT_Z, local_position + 2).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const {

    const unsigned int LocalSize = TDim * TNumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;
    unsigned int local_position = this->GetGeometry()[0].GetDofPosition(FLUID_FRACTION_GRADIENT_X);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(FLUID_FRACTION_GRADIENT_X,local_position);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(FLUID_FRACTION_GRADIENT_Y,local_position+1);
        if (TDim == 3) rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(FLUID_FRACTION_GRADIENT_Z,local_position+2);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeFluidFractionGradient<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    KRATOS_ERROR_IF(this->GetGeometry().size() != TNumNodes)<< "Wrong number of nodes for element" << this->Id() << std::endl;

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(FLUID_FRACTION_GRADIENT, this->GetGeometry()[i])
    }
    return 0;

    KRATOS_CATCH("");
}


template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                               const double Mass)
{
    unsigned int DofIndex = 0;
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rLHSMatrix(DofIndex, DofIndex) += Mass;
            ++DofIndex;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix,
        const array_1d<double,TNumNodes>& rShapeFunc,
        const double Weight)
{
    const unsigned int BlockSize = TDim;

    double Coef = Weight;

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        unsigned int row = i*BlockSize;
        // Loop over columns
        for (unsigned int j = 0; j < TNumNodes; ++j)
        {
            unsigned int col = j*BlockSize;
            const double K = Coef * rShapeFunc[i] * rShapeFunc[j];
            for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
            {
                rLHSMatrix(row+d, col+d) += K;
            }

        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = TDim * TNumNodes;

    // Resize and set to zero
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Get the element's geometric parameters

    // Add 'classical' mass matrix (lumped)
    if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == 1){
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        DenseVector<Matrix> shape_derivatives;
        Geometry<Node >& geom = this->GetGeometry();
        double Area = geom.Volume();
        GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        auto number_integration_points = geom.IntegrationPointsNumber(integration_method);
        Vector gauss_weights = ZeroVector(number_integration_points);
        Matrix shape_functions = ZeroMatrix(number_integration_points,TNumNodes);
        const std::vector<IntegrationPoint<3>>& IntegrationPoints = geom.IntegrationPoints(integration_method);
        Vector DetJ;
        geom.ShapeFunctionsIntegrationPointsGradients(shape_derivatives,DetJ,integration_method);
        shape_functions = geom.ShapeFunctionsValues(integration_method);

        for (unsigned int g = 0; g < number_integration_points; g++){
                gauss_weights[g] = DetJ[g] * IntegrationPoints[g].Weight();
                for (unsigned int i = 0; i < TNumNodes; ++i){
                    for (unsigned int j = 0; j < TNumNodes; ++j){
                        for (unsigned int d = 0; d < TDim; ++d)
                            DN_DX(j,d) += gauss_weights[g] * shape_functions(g,i) * shape_derivatives[g](j,d);
                    }
                }
        }
        double Coeff = Area / TNumNodes; //Optimize!
        this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);
    }

    else {
        // Add 'consistent' mass matrix
        MatrixType NContainer;
        ShapeFunctionDerivativesArrayType DN_DXContainer;
        VectorType GaussWeights;
        this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
        const SizeType NumGauss = NContainer.size1();

        for (SizeType g = 0; g < NumGauss; g++){
            const double GaussWeight = GaussWeights[g];
            const ShapeFunctionsType& Ng = row(NContainer, g);
            this->AddConsistentMassMatrixContribution(rMassMatrix, Ng, GaussWeight);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::CalculateRHS(VectorType& F, const ProcessInfo& rCurrentProcessInfo)
{
    // Get the element's geometric parameters
    array_1d<double, TNumNodes> N;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX;

    MatrixType NContainer;
    ShapeFunctionDerivativesArrayType DN_DXContainer;
    VectorType GaussWeights;
    this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
    const SizeType NumGauss = NContainer.size1();

    for (SizeType g = 0; g < NumGauss; g++){
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& Ng = row(NContainer, g);
        this->AddIntegrationPointRHSContribution(F, Ng, DN_DXContainer[g], GaussWeight);
    }

}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
        Matrix& rNContainer,
        Vector& rGaussWeights)
{
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::IntegrationMethod::GI_GAUSS_5);
    rNContainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_5);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_5);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_5), false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_5); g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::EvaluateInPoint(double& rResult,
                             const Variable< double >& rVariable,
                             const array_1d< double, TNumNodes >& rShapeFunc)
{
    // Compute the weighted value of the nodal variable in the (Gauss) Point
    rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
    for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
        rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeFluidFractionGradient<TDim, TNumNodes>::AddIntegrationPointRHSContribution(VectorType& F,
                             const array_1d<double, TNumNodes>& rShapeFunc,
                             const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    //double fluid_fraction;
    //this->EvaluateInPoint(fluid_fraction, FLUID_FRACTION, rShapeFunc);
    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < TNumNodes; ++iNodeB){

        for (unsigned int dj = 0; dj < TDim; ++dj){
            double value = 0.0;
            for (unsigned int iNodeA = 0; iNodeA < TNumNodes; ++iNodeA){
                double& r_nodal_fluid_fraction = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(FLUID_FRACTION);
                //for (unsigned int di = 0; di < TDim; ++di){
                    value += rShapeFunc[iNodeB] * rShapeDeriv(iNodeA, dj) * r_nodal_fluid_fraction;
                //}
            }
            F[LocalIndex++] += Coef * value;

        }
    }
}
// Explicit instantiations
template class ComputeFluidFractionGradient<2, 3>;
template class ComputeFluidFractionGradient<3, 4>;
template class ComputeFluidFractionGradient<2, 4>;
template class ComputeFluidFractionGradient<2, 9>;
template class ComputeFluidFractionGradient<3, 8>;
template class ComputeFluidFractionGradient<3, 27>;
} // namespace Kratos
