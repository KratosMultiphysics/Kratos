// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/total_lagrangian.h"

namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Constants for the computation of the stress
        const double E = p_elem_prop->GetValue(YOUNG_MODULUS);
        const double NU = p_elem_prop->GetValue(POISSON_RATIO);
        const double c1 = E / (1.00 - NU * NU);
        const double c2 = c1 * NU;
        const double c3 = 0.5* E / (1 + NU);


        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }

        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("TotalLagrangianElement2D3N", 1, element_nodes, p_elem_prop);

        p_element->Initialize(r_model_part.GetProcessInfo());

        // Define a matrix A and impose that the displacement on each node is u = A*x0
        Matrix A0 = ZeroMatrix(3,3);
        A0(0,0) = 2e-2;
        A0(0,1) = 5e-2;
        A0(1,0) = 5e-2;
        A0(1,1) = -1e-2;
        A0(2,2) = 1.0;

        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        double factor = 1.0;
        for(unsigned int step=0; step<5; step++)
        {
            factor+=1.0;
            Matrix A = factor*A0;
            for (auto& r_node : r_model_part.Nodes()){
                noalias(r_node.GetSolutionStepValue(DISPLACEMENT)) = prod(A, r_node.GetInitialPosition());
                r_node.Coordinates() = r_node.GetInitialPosition() + r_node.GetSolutionStepValue(DISPLACEMENT); //here i update the node coordinates
            }

            p_element->InitializeSolutionStep(const_procinfo_ref);
            p_element->InitializeNonLinearIteration(const_procinfo_ref);
            p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
            p_element->FinalizeNonLinearIteration(const_procinfo_ref);
            p_element->FinalizeSolutionStep(const_procinfo_ref);

            // Known A we can compute the deformation gradient analytically as F = A + I and from that the analytical value of the strain
            Matrix F = IdentityMatrix(3,3) + A;
            Matrix strain_mat = 0.5*(prod(trans(F),F)-IdentityMatrix(3,3));
            std::vector<double> reference_strain = {strain_mat(0,0), strain_mat(1,1), 2.0*strain_mat(0,1)};

            Matrix F22 = ZeroMatrix(2,2);
            for(unsigned int i=0; i<2; i++)
                for(unsigned int j=0; j<2; j++)
                    F22(i,j) = F(i,j);

            // Also the stress can be computed analytically as (PK2 stress)
            Vector reference_stress(3);
            reference_stress[0] = c1 * reference_strain[0] + c2 * reference_strain[1];
            reference_stress[1] = c2 * reference_strain[0] + c1 * reference_strain[1];
            reference_stress[2] = c3 * reference_strain[2];

            Matrix reference_stress_tensor =  MathUtils<double>::StressVectorToTensor( reference_stress );

            const double reference_von_mises_pk2 = std::sqrt( 
                                             reference_stress[0]*reference_stress[0] 
                                             + reference_stress[1]*reference_stress[1]
                                             - reference_stress[0]*reference_stress[1]
                                             + 3.0*reference_stress[2]*reference_stress[2]
                                            );

            //compute cauchy stress
            const auto& S = reference_stress_tensor;
            double J = MathUtils<double>::Det(F22);
            Matrix reference_cauchy_stress_tensor = 1/J*prod(F22, Matrix(prod(S,trans(F22))));



            std::vector<Vector> output_strains(1);
            p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,output_strains, r_model_part.GetProcessInfo());

            std::vector<Vector> output_stress(1);
            p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,output_stress, r_model_part.GetProcessInfo());

            std::vector<double> output_von_mises(1);
            p_element->CalculateOnIntegrationPoints(VON_MISES_STRESS,output_von_mises, r_model_part.GetProcessInfo());

            std::vector<Matrix> output_pk2_tensor(1);
            p_element->CalculateOnIntegrationPoints(PK2_STRESS_TENSOR,output_pk2_tensor, r_model_part.GetProcessInfo());

            std::vector<Vector> output_cauchy_vector(1);
            p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR,output_cauchy_vector, r_model_part.GetProcessInfo());

            std::vector<Matrix> output_cauchy_tensor(1);
            p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR,output_cauchy_tensor, r_model_part.GetProcessInfo());

            KRATOS_CHECK_VECTOR_EQUAL(output_strains[0], reference_strain);
            KRATOS_CHECK_VECTOR_EQUAL(output_stress[0], reference_stress);
            KRATOS_CHECK_NEAR((output_von_mises[0]-reference_von_mises_pk2)/reference_von_mises_pk2, 0.0,1e-14);
            KRATOS_CHECK_MATRIX_EQUAL(output_pk2_tensor[0], reference_stress_tensor);

            for(unsigned int i=0; i<reference_cauchy_stress_tensor.size1(); ++i){
                for(unsigned int j=0; j<reference_cauchy_stress_tensor.size2(); ++j){
                    KRATOS_CHECK_NEAR( (output_cauchy_tensor[0](i,j) - reference_cauchy_stress_tensor(i,j))/reference_cauchy_stress_tensor(i,j),0.0,1e-14);
                }
            }





        }
    }
}
}
