// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter,
//                   Vahid Galavi
//

// System includes

// External includes

// Project includes
#include "custom_elements/geo_linear_truss_element_base.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::GeoTrussElementLinearBase(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::GeoTrussElementLinearBase(IndexType NewId,
                                                                      GeometryType::Pointer pGeometry,
                                                                      PropertiesType::Pointer pProperties)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementLinearBase<TDim, TNumNodes>::Create(IndexType NewId,
                                                                    NodesArrayType const& rThisNodes,
                                                                    PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_intrusive<GeoTrussElementLinearBase>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementLinearBase<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                    GeometryType::Pointer pGeom,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElementLinearBase>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::~GeoTrussElementLinearBase()
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CreateElementStiffnessMatrix(MatrixType& rLocalStiffnessMatrix,
                                                                              const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    this->CalculateElasticStiffnessMatrix(rLocalStiffnessMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::AddPrestressLinear(VectorType& rRightHandSideVector)
{
    KRATOS_TRY

    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }
    const double A = this->GetProperties()[CROSS_AREA];
    const double N = prestress * A;

    // internal force vectors
    FullDofVectorType f_local = ZeroVector(TDim * TNumNodes);
    f_local[0]                = -1.00 * N;
    f_local[TDim]             = 1.00 * N;
    rRightHandSideVector -= prod(transformation_matrix, f_local);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rRightHandSideVector = ZeroVector(TDim * TNumNodes);

    FullDofVectorType internal_forces = ZeroVector(TDim * TNumNodes);
    this->UpdateInternalForces(internal_forces, rCurrentProcessInfo);

    noalias(rRightHandSideVector) -= internal_forces;

    AddPrestressLinear(rRightHandSideVector);

    // add bodyforces
    FullDofVectorType GlobalBodyForces;
    this->CalculateBodyForces(GlobalBodyForces);
    noalias(rRightHandSideVector) += GlobalBodyForces;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);
    // creating LHS
    this->CreateElementStiffnessMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType& integration_points =
        this->GetGeometry().IntegrationPoints();

    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == FORCE) {
        BoundedVector<double, TDim> truss_forces = ZeroVector(TDim);
        const double                A            = this->GetProperties()[CROSS_AREA];

        double prestress = 0.00;
        if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
        Vector temp_strain = ZeroVector(1);
        Vector temp_stress = ZeroVector(1);
        temp_strain[0]     = CalculateLinearStrain();
        Values.SetStrainVector(temp_strain);
        Values.SetStressVector(temp_stress);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        truss_forces[0] = (temp_stress[0] + prestress) * A;

        rOutput[0] = truss_forces;
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                              std::vector<Vector>& rOutput,
                                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        this->GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }
    if (rVariable == STRAIN) {
        Vector Strain = ZeroVector(TDim);
        Strain[0]     = CalculateLinearStrain();
        rOutput[0]    = Strain;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementLinearBase<3, 2>::WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates)
{
    KRATOS_TRY

    static constexpr unsigned int DIM       = 3;
    static constexpr unsigned int NUM_NODES = 2;

    rReferenceCoordinates    = ZeroVector(DIM * NUM_NODES);
    rReferenceCoordinates[0] = this->GetGeometry()[0].X0();
    rReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
    rReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
    rReferenceCoordinates[3] = this->GetGeometry()[1].X0();
    rReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
    rReferenceCoordinates[5] = this->GetGeometry()[1].Z0();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <>
void GeoTrussElementLinearBase<2, 2>::WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates)
{
    KRATOS_TRY

    static constexpr unsigned int DIM       = 2;
    static constexpr unsigned int NUM_NODES = 2;

    rReferenceCoordinates = ZeroVector(DIM * NUM_NODES);

    rReferenceCoordinates[0] = this->GetGeometry()[0].X0();
    rReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
    rReferenceCoordinates[2] = this->GetGeometry()[1].X0();
    rReferenceCoordinates[3] = this->GetGeometry()[1].Y0();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateLinearStrain()
{
    Vector current_disp = ZeroVector(TDim * TNumNodes);
    this->GetValuesVector(current_disp);
    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    current_disp = prod(Matrix(trans(transformation_matrix)), current_disp);

    double length_0;
    if constexpr (TDim == 2) {
        length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    } else if constexpr (TDim == 3) {
        length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    } else {
        KRATOS_ERROR << "Dimension of truss element should be either 2D or 3D" << std::endl;
    }

    const double e = (current_disp[TDim] - current_disp[0]) / length_0;

    return e;
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::UpdateInternalForces(FullDofVectorType& rInternalForces,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0]     = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    Vector temp_internal_stresses = ZeroVector(TDim * TNumNodes);
    temp_internal_stresses[0]     = -1.0 * temp_stress[0];
    temp_internal_stresses[TDim]  = temp_stress[0];

    rInternalForces = temp_internal_stresses * this->GetProperties()[CROSS_AREA];

    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0]     = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template class GeoTrussElementLinearBase<2, 2>;
template class GeoTrussElementLinearBase<3, 2>;

} // namespace Kratos.
