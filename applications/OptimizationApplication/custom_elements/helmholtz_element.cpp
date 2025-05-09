//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/helmholtz_solid_shape_data_container.h"
#include "custom_elements/data_containers/helmholtz_surface_data_container.h"
#include "custom_elements/data_containers/helmholtz_solid_data_container.h"
#include "custom_utilities/entity_calculation_utils.h"
#include "optimization_application_variables.h"

// Include base h
#include "custom_elements/helmholtz_element.h"

namespace Kratos {

//************************************************************************************
//************************************************************************************

template<class TDataContainer>
HelmholtzElement<TDataContainer>::HelmholtzElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry),
      mDataContainer(this->GetGeometry())
{
    // DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

template<class TDataContainer>
HelmholtzElement<TDataContainer>::HelmholtzElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties),
      mDataContainer(this->GetGeometry())
{
    // DO NOT ADD DOFS HERE!!!
}

template<class TDataContainer>
Element::Pointer HelmholtzElement<TDataContainer>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataContainer>
Element::Pointer HelmholtzElement<TDataContainer>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TDataContainer>
Element::Pointer HelmholtzElement<TDataContainer>::Clone(
    IndexType NewId,
    const NodesArrayType& rThisNodes) const
{
    KRATOS_TRY

    HelmholtzElement::Pointer p_new_cond =
        Kratos::make_intrusive<HelmholtzElement>(
            NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));

    return p_new_cond;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();

    if (rLeftHandSideMatrix.size1() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize); // resetting LHS

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(LocalSize); // resetting RHS

    MatrixType M;
    CalculateMassMatrix(M, rCurrentProcessInfo);

    MatrixType K;
    CalculateStiffnessMatrix(K, rCurrentProcessInfo);

    const bool is_inversed = rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE];

    noalias(rLeftHandSideMatrix) += M;
    if (!is_inversed) {
        noalias(rLeftHandSideMatrix) += K;
    }

    const bool is_integrated_field = rCurrentProcessInfo[HELMHOLTZ_INTEGRATED_FIELD];

    constexpr auto& source_variables_list = TDataContainer::SourceVariablesList;

    BoundedVector<double, NumberOfNodes * DataDimension> nodal_vals;
    IndexType local_index = 0;
    for (IndexType node_element = 0; node_element < NumberOfNodes; node_element++) {
        const auto& r_node = r_geometry[node_element];
        for (IndexType d = 0; d < DataDimension; ++d) {
            nodal_vals[local_index++] = r_node.GetValue(*source_variables_list[d]);
        }
    }

    if (is_integrated_field) {
        IndexType local_index = 0;
        for (IndexType node_element = 0; node_element < NumberOfNodes; node_element++) {
            const auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            for (IndexType d = 0; d < DataDimension; ++d) {
                nodal_vals[local_index++] /= node_weight;
            }
        }
        noalias(rRightHandSideVector) += nodal_vals;
    } else if (is_inversed) {
        noalias(rRightHandSideVector) += prod(K + M, nodal_vals);
    } else {
        noalias(rRightHandSideVector) += prod(M, nodal_vals);
    }

    // apply drichlet BC
    Vector temp;
    GetValuesVector(temp, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp);

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != LocalSize) {
         rValues.resize(LocalSize);
    }

    constexpr auto& r_variables_list = TDataContainer::TargetVariablesList;
    const auto& r_geometry = this->GetGeometry();

    IndexType local_index = 0;
    if constexpr(DataDimension == 1) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            rValues[local_index++] = r_geometry[i].FastGetSolutionStepValue(*r_variables_list[0]);
        }
    } else if constexpr(DataDimension == 3) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const auto& r_node = r_geometry[i];
            for (IndexType d = 0; d < DataDimension; ++d) {
                rValues[local_index++] = r_node.FastGetSolutionStepValue(*r_variables_list[d]);
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported data dimension. [ DataDimension = " << DataDimension
                     << ".\n";
    }
}
//******************************************************************************
//******************************************************************************

template<class TDataContainer>
void HelmholtzElement<TDataContainer>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == ELEMENT_STRAIN_ENERGY) {
        MatrixType K;
        CalculateStiffnessMatrix(K, rCurrentProcessInfo);

        auto& r_geometry = this->GetGeometry();

        const unsigned int NumberOfNodes = r_geometry.size();
        Vector nodal_vals(NumberOfNodes * 3);
        for (unsigned int node_element = 0; node_element < NumberOfNodes; node_element++) {
            nodal_vals[3 * node_element + 0] = r_geometry[node_element].X0();
            nodal_vals[3 * node_element + 1] = r_geometry[node_element].Y0();
            nodal_vals[3 * node_element + 2] = r_geometry[node_element].Z0();
        }

        rOutput = inner_prod(nodal_vals, prod(K, nodal_vals));
    }
    else {
        auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
        parentElement[0].Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0, 0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    if (rResult.size() != LocalSize) {
         rResult.resize(LocalSize);
    }

    constexpr auto& r_variables_list = TDataContainer::TargetVariablesList;
    const auto& r_geometry = this->GetGeometry();

    IndexType local_index = 0;
    if constexpr(DataDimension == 1) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            rResult[local_index++] = r_geometry[i].GetDof(*r_variables_list[0]).EquationId();
        }
    } else if constexpr(DataDimension == 3) {
        const auto pos = this->GetGeometry()[0].GetDofPosition(*r_variables_list[0]);
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const auto& r_node = r_geometry[i];
            for (IndexType d = 0; d < DataDimension; ++d) {
                rResult[local_index++] = r_node.GetDof(*r_variables_list[d], pos + d).EquationId();
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported data dimension. [ DataDimension = " << DataDimension
                     << ".\n";
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
void HelmholtzElement<TDataContainer>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    if (rElementalDofList.size() != LocalSize) {
         rElementalDofList.resize(LocalSize);
    }

    constexpr auto& r_variables_list = TDataContainer::TargetVariablesList;
    const auto& r_geometry = this->GetGeometry();

    IndexType local_index = 0;
    if constexpr(DataDimension == 1) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(*r_variables_list[0]);
        }
    } else if constexpr(DataDimension == 3) {
        const auto pos = this->GetGeometry()[0].GetDofPosition(*r_variables_list[0]);
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const auto& r_node = r_geometry[i];
            for (IndexType d = 0; d < DataDimension; ++d) {
                rElementalDofList[local_index++] = r_node.pGetDof(*r_variables_list[d], pos + d);
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported data dimension. [ DataDimension = " << DataDimension
                     << ".\n";
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<class TDataContainer>
int HelmholtzElement<TDataContainer>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_HELMHOLTZ_INVERSE))
        << "COMPUTE_HELMHOLTZ_INVERSE not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_INTEGRATED_FIELD))
        << "HELMHOLTZ_INTEGRATED_FIELD not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_RADIUS))
        << "HELMHOLTZ_RADIUS not defined in the ProcessInfo!" << std::endl;

    constexpr auto& r_variables_list = TDataContainer::TargetVariablesList;

    const auto& r_geometry = GetGeometry();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const auto& r_node = r_geometry[i];
        for (const auto p_variable : r_variables_list) {
            const auto& r_variable = *p_variable;
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_variable, r_node)
            KRATOS_CHECK_DOF_IN_NODE(r_variable, r_node)
        }
    }

    return Element::Check(rCurrentProcessInfo);

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

template<class TDataContainer>
void HelmholtzElement<TDataContainer>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rMassMatrix.size1() != LocalSize || rMassMatrix.size2() != LocalSize) {
        rMassMatrix.resize(LocalSize, LocalSize, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);

    const auto& r_cond_geom = GetGeometry();
    const auto& integration_method = r_cond_geom.GetDefaultIntegrationMethod();
    const auto& integration_points = r_cond_geom.IntegrationPoints(integration_method);

    VectorType Ws;
    MatrixType Ns;
    EntityCalculationUtils::CalculateElementGaussPointData(Ws, Ns, r_cond_geom, integration_method);

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        const double integration_weight = Ws[point_number];
        const Vector& rN = row(Ns, point_number);

        const auto& NiNj_weight = outer_prod(rN, rN) * integration_weight;

        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const SizeType index_i = i * DataDimension;
            for (IndexType j = 0; j < NumberOfNodes; ++j) {
                const SizeType index_j = j * DataDimension;
                for (IndexType k = 0; k < DataDimension; ++k) {
                    rMassMatrix(index_i + k, index_j + k) += NiNj_weight(i, j);
                }
            }
        }
    }

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

template<class TDataContainer>
void HelmholtzElement<TDataContainer>::CalculateStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();

    if (rStiffnessMatrix.size1() != LocalSize || rStiffnessMatrix.size2() != LocalSize) {
        rStiffnessMatrix.resize(LocalSize, LocalSize, false);
    }

    noalias(rStiffnessMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // reading integration points and local gradients
    const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
    const auto& integration_points = r_geom.IntegrationPoints(integration_method);
    const IndexType number_of_gauss_points = integration_points.size();

    Vector gauss_pts_J_det = ZeroVector(number_of_gauss_points);
    r_geom.DeterminantOfJacobian(gauss_pts_J_det, integration_method);

    typename TDataContainer::ConstantDataContainer constant_data(*this, integration_method, rCurrentProcessInfo);

    for (IndexType i_point = 0; i_point < number_of_gauss_points; ++i_point) {
        const double W = integration_points[i_point].Weight() * gauss_pts_J_det[i_point];
        mDataContainer.AddStiffnessGaussPointContributions(rStiffnessMatrix, W, i_point, constant_data);
    }

    KRATOS_CATCH("");
}

// template instantiations
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSurfaceDataContainer<3, 3, 1>>;
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSurfaceDataContainer<3, 4, 1>>;

template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSurfaceDataContainer<3, 3, 3>>;
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSurfaceDataContainer<3, 4, 3>>;

template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidDataContainer<3, 4, 1>>;
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidDataContainer<3, 8, 1>>;

template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidDataContainer<3, 4, 3>>;
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidDataContainer<3, 8, 3>>;

template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidShapeDataContainer<3, 4>>;
template class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzElement<HelmholtzSolidShapeDataContainer<3, 8>>;

} // Namespace Kratos
