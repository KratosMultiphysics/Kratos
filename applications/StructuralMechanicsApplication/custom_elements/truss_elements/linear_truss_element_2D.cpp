// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "linear_truss_element_2D.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                if constexpr (NNodes == 2) {
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                } else {
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
                }
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::InitializeMaterial()
{
    KRATOS_TRY
    const auto &r_props = GetProperties();

    if (r_props[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = r_props[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
Element::Pointer LinearTrussElement2D<TNNodes>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTrussElement2D<TNNodes>::Pointer p_new_elem = Kratos::make_intrusive<LinearTrussElement2D<TNNodes>>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    IndexType local_index = 0;

    if (rResult.size() != SystemSize)
        rResult.resize(SystemSize, false);

    const IndexType xpos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    rElementalDofList.resize(SystemSize);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * DofsPerNode;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);
    if constexpr (NNodes == 2) {
        rN[0] = 0.5 * (1.0 - xi);
        rN[2] = 0.5 * (1.0 + xi);
    } else { // 3N
        rN[0] = 0.5 * xi * (xi - 1.0);
        rN[2] = (1.0 - std::pow(xi, 2));
        rN[4] = 0.5 * xi * (xi + 1.0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rdN_dX,
    const double Length,
    const double xi
    ) const
{
    if (rdN_dX.size() != SystemSize)
        rdN_dX.resize(SystemSize, false);
    if constexpr (NNodes == 2) {
        const double inverse_l = 1.0 / Length;
        rdN_dX[0] = -inverse_l;
        rdN_dX[2] = inverse_l;
    } else { // 3N
        rdN_dX[0] = xi - 0.5;
        rdN_dX[2] = -2.0 * xi;
        rdN_dX[4] = xi + 0.5;
        rdN_dX *= 2.0 / Length;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetNodalValuesVector(SystemSizeBoundedArrayType& rNodalValues) const
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
array_1d<double, 3> LinearTrussElement2D<TNNodes>::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
)
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
int LinearTrussElement2D<TNNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template class LinearTrussElement2D<2>;
template class LinearTrussElement2D<3>;

} // Namespace Kratos
