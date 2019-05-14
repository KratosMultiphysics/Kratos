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
{/*
    if(this->Id() == 0)
        {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Vector& t = GetValue(TANGENTS);
        Vector actual_base_vector_update;
        IgaCurveOnSurfaceUtilities::CalculateTangent(
        GetGeometry(),
        DN_De,
        t,
        actual_base_vector_update);

        return actual_base_vector_update;

        mReferenceBaseVector = actual_base_vector_update;
        }
    //rTangentVector = g1 * rTangents[0] + g2 * rTangents[1];
    else
        {
        mReferenceBaseVector = GetActualBaseVector();
        }*/
    mReferenceBaseVector = GetActualBaseVector();
    
}
/*
IgaEdgeCableElement::Vector3 IgaEdgeCableElement::GetReferenceBaseVector() const
{
    Vector reference_base_vector;
//if (this->Id() == 0) {
if (mode == "UPDATED") {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Vector& t = GetValue(TANGENTS);
        array_1d<double, 3> actual_base_vector_update = ZeroVector(3);
        IgaCurveOnSurfaceUtilities::CalculateTangent(
        GetGeometry(),
        DN_De,
        t,
        actual_base_vector_update);

        return actual_base_vector_update;

        reference_base_vector = actual_base_vector_update;
        }
    else {
        reference_base_vector = mReferenceBaseVector;
    }
    return reference_base_vector;
}*/
//KRATOS_WATCH(mReferenceBaseVector)
//Definition von Base Vector 

IgaEdgeCableElement::Vector3 IgaEdgeCableElement::GetActualBaseVector() const
{
  // void CalculateTangent(
  //      const Matrix& rDN_De,
  //      const array_1d<double, 2>& Tangents);{
    
  //  return Tangents()
  //  }
    const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    const Vector& t = GetValue(TANGENTS);
    //const Vector& t = Tangent1;

    //const double t_u = t[0];
    //const double t_v = t[1];
    
    array_1d<double, 3> actual_base_vector = ZeroVector(3);

    IgaCurveOnSurfaceUtilities::CalculateTangent(
        GetGeometry(),
        DN_De,
        t,
        actual_base_vector
    );
//if(this->Id() == 0)
//KRATOS_WATCH(actual_base_vector)

    //for (std::size_t k = 0; k < NumberOfNodes(); k++) // k = Number of Nodes Cable
    //{
    //    actual_base_vector[0] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].X();        
    //    actual_base_vector[1] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].Y();        
    //    actual_base_vector[2] += ( DN_De(k, 0) * t[0] + DN_De(k,1) * t[1] ) * GetGeometry()[k].Z();
    //}

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

//KRATOS_WATCH(integration_weight)

    // get properties

    const auto& properties = GetProperties();

    const double E = properties[YOUNG_MODULUS];
    const double A = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];

//KRATOS_WATCH(E)
//KRATOS_WATCH(A)
//KRATOS_WATCH(prestress)
//KRATOS_WATCH(rLeftHandSideMatrix)

    // compute base vectors

    const Vector3 actual_base_vector = GetActualBaseVector();
    //const Vector3 reference_base_vector = GetReferenceBaseVector();

    const double reference_a = norm_2(mReferenceBaseVector);
   // const double reference_a = norm_2(reference_base_vector);
    const double actual_a = norm_2(actual_base_vector);

    const double actual_aa = actual_a * actual_a;
    const double reference_aa = reference_a * reference_a;
/*
KRATOS_WATCH(actual_base_vector)

KRATOS_WATCH(mReferenceBaseVector)
KRATOS_WATCH(actual_aa)
KRATOS_WATCH(reference_aa)
*/
    // green-lagrange strain

    const double e11_membrane = 0.5 * (actual_aa - reference_aa);

    // normal forcereference_aa

    const double s11_membrane = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector,mReferenceBaseVector); //reference_aa;

//KRATOS_WATCH(prestress)
//KRATOS_WATCH(e11_membrane)
//KRATOS_WATCH(s11_membrane)

    //for (std::size_t k = 0; k < NumberOfDofs(); k++){
    //    const double dof_type_m[NumberOfDofs()];
    //    const std::size_t dof_type_m[ ] = { GetDofTypeIndex(k) };
    
    for (std::size_t r = 0; r < NumberOfDofs(); r++) {
        //const std::size_t dof_type_g[r] = {GetDofTypeIndex(r)};
        //const std::size_t dof_type_r = dof_type_m[r];
        const std::size_t dof_type_r = GetDofTypeIndex(r);
        const std::size_t shape_index_r = GetShapeIndex(r);
        //const double t_u = t[0];
        //const double t_v = t[1];

        //KRATOS_WATCH(r)
        //KRATOS_WATCH(dof_type_r)
        //KRATOS_WATCH(shape_index_r)
        //if(reference_a == actual_a) //NEUNEU 
        //{//NEUNEU
        //const double epsilon_var_r = 0;//NEUNEU
        //}//NEUNEU
        //else//NEUNEU
        //{ //NEUNEU
        const double epsilon_var_r = actual_base_vector[dof_type_r] *
            (shape_derivatives(shape_index_r, 0) * t[0] 
            + shape_derivatives(shape_index_r, 1) * t[1]) / inner_prod(mReferenceBaseVector,mReferenceBaseVector);//reference_aa;
        //} //NEUNEU
        //KRATOS_WATCH(epsilon_var_r)
 
       if (ComputeLeftHandSide) {
            for (std::size_t s = 0; s < NumberOfDofs(); s++) {
                const std::size_t dof_type_s = GetDofTypeIndex(s);
                const std::size_t shape_index_s = GetShapeIndex(s);
                const Vector& t = GetValue(TANGENTS);
                //const double t_u = t[0];
                //const double t_v = t[1];

            //KRATOS_WATCH(s)
            //KRATOS_WATCH(dof_type_s)
            //KRATOS_WATCH(shape_index_s)
 
                const double epsilon_var_s =
                    actual_base_vector[dof_type_s] *
                    (shape_derivatives(shape_index_s, 0) * t[0]
                    + shape_derivatives(shape_index_s, 1) * t[1])
                    / inner_prod(mReferenceBaseVector,mReferenceBaseVector);//reference_aa;

                //KRATOS_WATCH(epsilon_var_s)

                rLeftHandSideMatrix(r, s) = 
                E * A * epsilon_var_r *
                    epsilon_var_s;

            
 
                if (dof_type_r == dof_type_s) {
                    const double epsilon_var_rs =
                        (shape_derivatives(shape_index_r, 0) * t[0] + shape_derivatives(shape_index_r, 1) * t[1]) *
                        (shape_derivatives(shape_index_s, 0) * t[0] + shape_derivatives(shape_index_s, 1) * t[1]) /inner_prod(mReferenceBaseVector,mReferenceBaseVector);// reference_aa;
                    
                //KRATOS_WATCH(epsilon_var_rs)
 
                    rLeftHandSideMatrix(r, s) += s11_membrane * epsilon_var_rs;

                //KRATOS_WATCH(rLeftHandSideMatrix)
                }
            }
        }
        if (ComputeRightHandSide) {
            rRightHandSideVector[r] = -s11_membrane * epsilon_var_r;
            
         //KRATOS_WATCH(rLeftHandSideMatrix)
         //KRATOS_WATCH(rRightHandSideVector)
        }
    }
    //}

//KRATOS_WATCH(reference_a)

    if (ComputeLeftHandSide) {
        rLeftHandSideMatrix *= reference_a * integration_weight;
    }

    if (ComputeRightHandSide) {
        rRightHandSideVector *= reference_a * integration_weight;
    }

//KRATOS_WATCH(rLeftHandSideMatrix)
//KRATOS_WATCH(rRightHandSideVector)
    //if(this->Id() == 0){KRATOS_WATCH(s11_membrane)}

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
        )// override;
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

        const Vector3 actual_base_vector = GetActualBaseVector();
        const double reference_a = norm_2(mReferenceBaseVector);
        const double actual_a = norm_2(actual_base_vector);

        const double actual_aa = actual_a * actual_a;
        const double reference_aa = reference_a * reference_a;

        // green-lagrange strain
        const double e11_membrane = 0.5 * (actual_aa - reference_aa);

        // normal forcereference_aa
        double principal_stress = prestress * A + e11_membrane * A * E / inner_prod(mReferenceBaseVector,mReferenceBaseVector); 
        //double principal_stress = e11_membrane * A * E / inner_prod(mReferenceBaseVector,mReferenceBaseVector); 

        rOutput = principal_stress;
        } 
    }  

} // namespace Kratos