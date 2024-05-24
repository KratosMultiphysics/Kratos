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
#include "custom_elements/geo_cable_element.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

#include "includes/checks.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
GeoCableElement<TDim, TNumNodes>::GeoCableElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTrussElement<TDim, TNumNodes>(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoCableElement<TDim, TNumNodes>::GeoCableElement(IndexType               NewId,
                                                  GeometryType::Pointer   pGeometry,
                                                  PropertiesType::Pointer pProperties)
    : GeoTrussElement<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoCableElement<TDim, TNumNodes>::Create(IndexType               NewId,
                                                          NodesArrayType const&   rThisNodes,
                                                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_intrusive<GeoCableElement>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoCableElement<TDim, TNumNodes>::Create(IndexType               NewId,
                                                          GeometryType::Pointer   pGeom,
                                                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCableElement>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoCableElement<TDim, TNumNodes>::~GeoCableElement()
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCableElement<TDim, TNumNodes>::CreateElementStiffnessMatrix(MatrixType& rLocalStiffnessMatrix,
                                                                    const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    if (mIsCompressed) {
        rLocalStiffnessMatrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);
    } else {
        this->CalculateElasticStiffnessMatrix(rLocalStiffnessMatrix, rCurrentProcessInfo);

        FullDofMatrixType K_geo;
        this->CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

        rLocalStiffnessMatrix += K_geo;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCableElement<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    rRightHandSideVector = ZeroVector(TDim * TNumNodes);

    BoundedVector<double, TDim * TNumNodes> internal_forces = ZeroVector(TDim * TNumNodes);
    UpdateInternalForces(internal_forces, rCurrentProcessInfo);

    if (!mIsCompressed) {
        noalias(rRightHandSideVector) -= internal_forces;
    }

    // add bodyforces
    BoundedVector<double, TDim * TNumNodes> GlobalBodyForces;
    this->CalculateBodyForces(GlobalBodyForces);
    noalias(rRightHandSideVector) += GlobalBodyForces;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCableElement<TDim, TNumNodes>::UpdateInternalForces(BoundedVector<double, TDim * TNumNodes>& rInternalForces,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double A  = this->GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0]     = this->CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    mInternalStresses = temp_stress;
    temp_stress += mInternalStressesFinalizedPrevious;

    const double normal_force = ((temp_stress[0] + prestress) * l * A) / L0;

    mIsCompressed = false;
    if ((normal_force < 0.00) && (std::abs(l - L0) > numerical_limit)) {
        mIsCompressed = true;
    }

    // internal force vectors
    BoundedVector<double, TDim * TNumNodes> f_local = ZeroVector(TDim * TNumNodes);
    f_local[0]                                      = -1.00 * normal_force;
    f_local[TDim]                                   = 1.00 * normal_force;
    rInternalForces                                 = ZeroVector(TDim * TNumNodes);
    noalias(rInternalForces)                        = prod(transformation_matrix, f_local);
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCableElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                    std::vector<array_1d<double, 3>>& rOutput,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == FORCE) {
        GeoTrussElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        if (rOutput[0][0] < 0.0) {
            rOutput[0] = ZeroVector(3);
        }
    }
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoCableElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                    std::vector<Vector>& rOutput,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR || rVariable == CAUCHY_STRESS_VECTOR ||
        rVariable == PK2_STRESS_VECTOR) {
        GeoTrussElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        if (rOutput[0][0] < 0.0) {
            rOutput[0] = ZeroVector(TDim);
        }
    }
}

//--------------------------------------------------------------------------------------------
template class GeoCableElement<2, 2>;
template class GeoCableElement<3, 2>;

} // namespace Kratos.
