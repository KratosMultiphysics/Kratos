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
#include "iga_edge_cable_element.h"
#include "iga_application_variables.h"

namespace Kratos {              // Geometrie von Menmbrane einarbeiten 

Element::Pointer IgaEdgeCableElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const      //kann ich hier bei Properties feslegen zu welchem Element die Punkte gehören? bzw. neue Eigenschaft hinzufüen?
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaEdgeCableElement>(NewId, geometry,
        pProperties);
}
//Hier sollen nachher zwei Geometrie Blöcke stehen mit Konten mit einem Index für Cable und Knoten mit zwei Indizes für Membrane 

void IgaEdgeCableElement::GetDofList(       // DOFs nur für CPs Cable oder werden auch die von der Membrane benötigt?
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

void IgaEdgeCableElement::EquationIdVector(
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

void IgaEdgeCableElement::Initialize()
{
    mReferenceBaseVector = GetActualBaseVector();
}

 //Neue Definition von Base Vector g1

IgaEdgeCableElement::Vector3 IgaEdgeCableElement::GetActualBaseVector() const
{
    const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // Base Vector Teile (2) und (4) (aufsummiert über CPs Cable)
 
    Vector3 actual_base_vector = ZeroVector(3);

    for (std::size_t k = 0; k < NumberOfNodes(); k++) // k = Number of Nodes Cable
    {
        actual_base_vector[0] += DN_De(k, 0) * GetGeometry()[k].X();        //entspricht (2)
        actual_base_vector[1] += DN_De(k, 0) * GetGeometry()[k].Y();        //entspricht (4)
        actual_base_vector[2] += DN_De(k, 0) * GetGeometry()[k].Z();
    }

    return actual_base_vector;

    // Gesamter Base Vector MEMBRANE

    //ShapeFunction + Derivatives Membrane importieren
    const Vector& ShapeFunctionNsrf = GetValue(SHAPE_FUNCTUIN_VALUES_MEMBRANE);        
    const Matrix& DNsrf_Du = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MEMBRANE_U);
    const Matrix& DNsrf_Dv = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MEMBRANE_V);

    for (std::size_t i=0; i < NumberOfNodes(); i++)
    {
        final_actual_base_vector[0] = ( DNsrf_Du(i, 0) * actual_base_vector[0] + DNsrf_Dv(i, 1) * actual_base_vector[1] ) * GetGeometry()[i].X();
        final_actual_base_vector[0] = ( DNsrf_Du(i, 0) * actual_base_vector[0] + DNsrf_Dv(i, 1) * actual_base_vector[1] ) * GetGeometry()[i].Y();
        final_actual_base_vector[0] = ( DNsrf_Du(i, 0) * actual_base_vector[0] + DNsrf_Dv(i, 1) * actual_base_vector[1] ) * GetGeometry()[i].Z();
    }

        return final_actual_base_vector;
}

void IgaEdgeCableElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    // get integration data
    
    const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
    Matrix& shape_derivatives_u = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MEMBRANE_U);
    Matrix& shape_derivatives_v = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MEMBRANE_V);

    // get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    const Vector3 actual_base_vector = GetActualBaseVector();

    const double reference_a = norm_2(mReferenceBaseVector);
    const double actual_a = actual_base_vector;
    const double final_actual_a = norm_2(final_actual_base_vector);

    const double final_actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;

    // green-lagrange strain

    const double e11_membrane = 0.5 * (final_actual_aa - reference_aa);

    // normal force

    const double s11_membrane = prestress * A + e11_membrane * A * E /
        reference_aa;


    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        const std::size_t dof_type_r = GetDofTypeIndex(r);
        const std::size_t shape_index_r = GetShapeIndex(r);

        const double epsilon_var_r = final_actual_a[dof_type_r] *
            (shape_derivatives_u(shape_index_r, 0) * acutal_a[0] + shape_derivatives_v(shape_index_r, 0) * acutal_a[1]) / reference_aa;

       if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                const std::size_t dof_type_s = GetDofTypeIndex(s);
                const std::size_t shape_index_s = GetShapeIndex(s);

                const double epsilon_var_s =
                    final_actual_a[dof_type_s] *
                    (shape_derivatives_u(shape_index_r, 0) * acutal_a[0] + shape_derivatives_v(shape_index_r, 0) * acutal_a[1]) / reference_aa;

                rLeftHandSideMatrix(r, s) = E * A * epsilon_var_r *
                    epsilon_var_s;
//
//                if (dof_type_r == dof_type_s) {
//                    const double epsilon_var_rs =
//                        shape_derivatives(shape_index_r, 0) *
//                        shape_derivatives(shape_index_s, 0) / reference_aa;
//
//                    rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;
//                }
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

void IgaEdgeCableElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaEdgeCableElement\" #" << Id();
}

} // namespace Kratos