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
//#include "includes/define.h"
//#include "includes/variables.h"

// External includes

// Project includes
#include "iga_edge_cable_element.h"

namespace Kratos {              

Element::Pointer IgaEdgeCableElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const     
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaEdgeCableElement>(NewId, geometry,
        pProperties);
} 

void IgaEdgeCableElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    const int number_of_control_points = NumberOfNodes();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(NumberOfDofs());

    for (unsigned int i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    KRATOS_CATCH("")
}

void IgaEdgeCableElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    const int number_of_control_points = NumberOfNodes();

    if (rResult.size() != 3 * number_of_control_points)
        rResult.resize(3 * number_of_control_points, false);

    for (unsigned int i = 0; i < number_of_control_points; ++i) {
        const unsigned int index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH("")
}

void IgaEdgeCableElement::Initialize()
{
    mReferenceBaseVector = GetActualBaseVector();   
}

array_1d<double, 3> IgaEdgeCableElement::GetActualBaseVector() const
{
    const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    const Vector& t = GetValue(TANGENTS);
       
    array_1d<double, 3> actual_base_vector = ZeroVector(3);

    IgaCurveOnSurfaceUtilities::CalculateTangent(
        GetGeometry(),
        DN_De,
        t,
        actual_base_vector
    );

    return actual_base_vector;
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
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    const Vector& t = GetValue(TANGENTS);

// get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

// compute base vectors

    const array_1d<double, 3> actual_base_vector = GetActualBaseVector();
    
    const double reference_a = norm_2(mReferenceBaseVector);
    const double actual_a = norm_2(actual_base_vector);

    const double actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;

// green-lagrange strain

    const double e11_membrane = 0.5 * (inner_prod(actual_base_vector, actual_base_vector) - inner_prod(mReferenceBaseVector, mReferenceBaseVector));

// normal forcereference_aa

    const double s11_membrane = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector,mReferenceBaseVector); 
    
    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        const std::size_t dof_type_r = r % 3;
        const std::size_t shape_index_r = r / 3;
       
        const double epsilon_var_r = actual_base_vector[dof_type_r] *
            (shape_derivatives(shape_index_r, 0) * t[0] 
            + shape_derivatives(shape_index_r, 1) * t[1]) / inner_prod(mReferenceBaseVector,mReferenceBaseVector);
 
       if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                const std::size_t dof_type_s = s % 3;
                const std::size_t shape_index_s = s / 3;
                const Vector& t = GetValue(TANGENTS);
          
                const double epsilon_var_s =
                    actual_base_vector[dof_type_s] *
                    (shape_derivatives(shape_index_s, 0) * t[0]
                    + shape_derivatives(shape_index_s, 1) * t[1])
                    / inner_prod(mReferenceBaseVector,mReferenceBaseVector);

                rLeftHandSideMatrix(r, s) = 
                E * A * epsilon_var_r *
                    epsilon_var_s;

            
 
                if (dof_type_r == dof_type_s) {
                    const double epsilon_var_rs =
                        (shape_derivatives(shape_index_r, 0) * t[0] + shape_derivatives(shape_index_r, 1) * t[1]) *
                        (shape_derivatives(shape_index_s, 0) * t[0] + shape_derivatives(shape_index_s, 1) * t[1]) /inner_prod(mReferenceBaseVector,mReferenceBaseVector);
                     
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

void IgaEdgeCableElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaEdgeCableElement\" #" << Id();
}

//********************************************************************************************************************
//********************************************************************************************************************
//********************************************************************************************************************

void IgaEdgeCableElement::Calculate(
        const Variable<double>& rVariable,
        double& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        if (rVariable == CABLE_STRESS)
        {
        const double& integration_weight = GetValue(INTEGRATION_WEIGHT);
        Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Vector& t = GetValue(TANGENTS);
        
        const auto& properties = GetProperties();

        const double E = properties[YOUNG_MODULUS];
        const double A = properties[CROSS_AREA];
        const double prestress = properties[PRESTRESS_CAUCHY];

        const array_1d<double, 3> actual_base_vector = GetActualBaseVector();
        const double reference_a = norm_2(mReferenceBaseVector);
        const double actual_a = norm_2(actual_base_vector);

        const double actual_aa = actual_a * actual_a;
        const double reference_aa = reference_a * reference_a;

        // green-lagrange strain
        const double e11_membrane = 0.5 * (actual_aa - reference_aa);

        // normal forcereference_aa
        double principal_stress = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector,mReferenceBaseVector); 
        
        rOutput = principal_stress;
        } 
    }  

} // namespace Kratos