/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Anna Bauer
//                  Thomas Oberbichler
//                  Tobias Teschemacher
*/

// System includes
#include "includes/define.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "iga_truss_element.h"
#include "iga_application_variables.h"

namespace Kratos {

Element::Pointer IgaTrussElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaTrussElement>(NewId, geometry,
        pProperties);
}

void IgaTrussElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaTrussElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void IgaTrussElement::Initialize()
{
    mReferenceBaseVector = GetActualBaseVector();
}

IgaTrussElement::Vector3 IgaTrussElement::GetActualBaseVector() const
{
    const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    Vector3 actual_base_vector = ZeroVector(3);

    for (std::size_t i = 0; i < NumberOfNodes(); i++)
    {
        actual_base_vector[0] += DN_De(0, i) * GetGeometry()[i].X();
        actual_base_vector[1] += DN_De(0, i) * GetGeometry()[i].Y();
        actual_base_vector[2] += DN_De(0, i) * GetGeometry()[i].Z();
    }

    return actual_base_vector;
}

void IgaTrussElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    // get integration data
    
    const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    const Vector3 actual_base_vector = GetActualBaseVector();

    const double reference_a = norm_2(mReferenceBaseVector);
    const double actual_a = norm_2(actual_base_vector);

    const double actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;

    // green-lagrange strain

    const double e11_membrane = 0.5 * (actual_aa - reference_aa);

    // normal force

    const double s11_membrane = prestress * A + e11_membrane * A * E /
        reference_aa;

    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        const std::size_t dof_type_r = GetDofTypeIndex(r);
        const std::size_t shape_index_r = GetShapeIndex(r);

        const double epsilon_var_r = actual_base_vector[dof_type_r] *
            shape_derivatives(shape_index_r, 0) / reference_aa;

        if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                const std::size_t dof_type_s = GetDofTypeIndex(s);
                const std::size_t shape_index_s = GetShapeIndex(s);

                const double epsilon_var_s =
                    actual_base_vector[dof_type_s] *
                    shape_derivatives(shape_index_s, 0) / reference_aa;

                rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r *
                    epsilon_var_s;

                if (dof_type_r == dof_type_s) {
                    const double epsilon_var_rs =
                        shape_derivatives(shape_index_r, 0) *
                        shape_derivatives(shape_index_s, 0) / reference_aa;

                    rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;
                }
            }
        }

        if (ComputeRightHandSide) {
            rRightHandSideVector[r] = -s11_membrane * epsilon_var_r;
        }
    }

    if (ComputeLeftHandSide) {
        rLeftHandSideMatrix *= reference_a * integration_weight;
    }

    if (ComputeRightHandSide) {
        rRightHandSideVector *= reference_a * integration_weight;
    }

    KRATOS_CATCH("")
}

void IgaTrussElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaTrussElement\" #" << Id();
}

} // namespace Kratos