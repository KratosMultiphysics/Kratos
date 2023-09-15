// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"
#include "custom_elements/small_displacement.h"
#include "includes/table_accessor.h"

namespace Kratos::Testing
{

    namespace
    {
        void Fill2D3NModelpart(
            ModelPart &rModelPart,
            const std::string& rConstitutiveLawName)
        {
            rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
            rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
            rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);

            // Set the element properties
            auto p_elem_prop = rModelPart.CreateNewProperties(0);
            p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
            p_elem_prop->SetValue(POISSON_RATIO, 0.3);
            p_elem_prop->SetValue(THERMAL_EXPANSION_COEFFICIENT, 7.2e-6);
            p_elem_prop->SetValue(REFERENCE_TEMPERATURE, 0.0);
            const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get(rConstitutiveLawName);
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            // Create the test element
            auto p_node_1 = rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            auto p_node_3 = rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            auto p_node_2 = rModelPart.CreateNewNode(3, 0.0 , 1.0 , 0.0);

            for (auto& r_node : rModelPart.Nodes()){
                r_node.AddDof(DISPLACEMENT_X);
                r_node.AddDof(DISPLACEMENT_Y);
                r_node.AddDof(DISPLACEMENT_Z);
            }

            std::vector<ModelPart::IndexType> element_nodes {1,2,3};
            auto p_element = rModelPart.CreateNewElement("SmallDisplacementElement2D3N", 1, element_nodes, p_elem_prop);
            p_element->Initialize(rModelPart.GetProcessInfo());

            // Set a displacement and temperature
            array_1d<double, 3> aux_disp = ZeroVector(3);
            noalias(p_node_1->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
            noalias(p_node_2->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
            aux_disp(0) = 0.001;
            noalias(p_node_3->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        }
    }

    /**
    * Checks the ThermalLinearPlaneStrain
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementElement2D3NThermalLinearPlaneStrain, KratosConstitutiveLawsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        Fill2D3NModelpart(r_model_part, "ThermalLinearPlaneStrain");

        r_model_part.pGetNode(1)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;
        r_model_part.pGetNode(2)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;
        r_model_part.pGetNode(3)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;

        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        auto p_element = r_model_part.pGetElement(1);
        p_element->InitializeSolutionStep(const_procinfo_ref);
        p_element->InitializeNonLinearIteration(const_procinfo_ref);
        p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
        p_element->FinalizeNonLinearIteration(const_procinfo_ref);
        p_element->FinalizeSolutionStep(const_procinfo_ref);

        std::vector<Vector> output_strains(1);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,output_strains, const_procinfo_ref);

        std::vector<Vector> output_stress(1);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,output_stress, const_procinfo_ref);

        std::vector<double> output_von_mises(1);
        p_element->CalculateOnIntegrationPoints(VON_MISES_STRESS,output_von_mises, const_procinfo_ref);

        Vector reference_stress(3);
        reference_stress(0) = -1627.69;
        reference_stress(1) = -3166.15;
        reference_stress(2) = 0.0;
        Vector reference_strain(3);
        reference_strain(0) = -0.0001232;
        reference_strain(1) = -0.0011232;
        reference_strain(2) = 0.0;
        const double reference_von_mises_pk2 = 2742.33;

        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_strains[0], reference_strain, 1e-4);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_stress[0], reference_stress, 1e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(output_von_mises[0],reference_von_mises_pk2, 1.0e-4);
    }

    /**
    * Checks the ThermalLinearPlaneStress
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementElement2D3NThermalLinearPlaneStress, KratosConstitutiveLawsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        Fill2D3NModelpart(r_model_part, "ThermalLinearPlaneStress");

        r_model_part.pGetNode(1)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;
        r_model_part.pGetNode(2)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;
        r_model_part.pGetNode(3)->FastGetSolutionStepValue(TEMPERATURE) = 120.0;
        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        auto p_element = r_model_part.pGetElement(1);
        p_element->InitializeSolutionStep(const_procinfo_ref);
        p_element->InitializeNonLinearIteration(const_procinfo_ref);
        p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
        p_element->FinalizeNonLinearIteration(const_procinfo_ref);
        p_element->FinalizeSolutionStep(const_procinfo_ref);

        std::vector<Vector> output_strains(1);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,output_strains, const_procinfo_ref);

        std::vector<Vector> output_stress(1);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,output_stress, const_procinfo_ref);

        std::vector<double> output_von_mises(1);
        p_element->CalculateOnIntegrationPoints(VON_MISES_STRESS,output_von_mises, const_procinfo_ref);

        Vector reference_stress(3);
        reference_stress(0) = -270.769;
        reference_stress(1) = -1809.23;
        reference_stress(2) = 0.0;
        Vector reference_strain(3);
        reference_strain(0) = 0.000136;
        reference_strain(1) = -0.000864;
        reference_strain(2) = 0.0;
        const double reference_von_mises_pk2 = 1690.19;

        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_strains[0], reference_strain, 1e-4);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_stress[0], reference_stress, 1e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(output_von_mises[0],reference_von_mises_pk2, 1.0e-4);
    }

    /**
    * Checks the ThermalLinearPlaneStrainWithTable
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementElement2D3NThermalLinearPlaneStrainWithTable, KratosConstitutiveLawsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        Fill2D3NModelpart(r_model_part, "ThermalLinearPlaneStrain");

        Table<double> temp_E_table;
        temp_E_table.insert(-1.0, 2.0e6);
        temp_E_table.insert(3500.0, 2.1e4);
        temp_E_table.insert(1.0e6, 2.1e4);

        auto p_elem_prop = r_model_part.pGetProperties(0);
        p_elem_prop->SetTable(TEMPERATURE, YOUNG_MODULUS, temp_E_table);

        TableAccessor table_accessor = TableAccessor(TEMPERATURE, "node_historical");
        p_elem_prop->SetAccessor(YOUNG_MODULUS, table_accessor.Clone());

        r_model_part.pGetNode(1)->FastGetSolutionStepValue(TEMPERATURE) = 2400.2;
        r_model_part.pGetNode(2)->FastGetSolutionStepValue(TEMPERATURE) = 2000.2;
        r_model_part.pGetNode(3)->FastGetSolutionStepValue(TEMPERATURE) = 1500.2;
        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        auto p_element = r_model_part.pGetElement(1);
        p_element->InitializeSolutionStep(const_procinfo_ref);
        p_element->InitializeNonLinearIteration(const_procinfo_ref);
        p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
        p_element->FinalizeNonLinearIteration(const_procinfo_ref);
        p_element->FinalizeSolutionStep(const_procinfo_ref);

        std::vector<Vector> output_strains(1);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,output_strains, const_procinfo_ref);

        std::vector<Vector> output_stress(1);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,output_stress, const_procinfo_ref);

        std::vector<double> output_von_mises(1);
        p_element->CalculateOnIntegrationPoints(VON_MISES_STRESS,output_von_mises, const_procinfo_ref);

        Vector reference_stress(3);
        reference_stress(0) = -30230.4;
        reference_stress(1) = -30913.2;
        reference_stress(2) = 0.0;
        Vector reference_strain(3);
        reference_strain(0) = -0.0174099;
        reference_strain(1) = -0.0184099;
        reference_strain(2) = 0.0;
        const double reference_von_mises_pk2 = 30577.5;

        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_strains[0], reference_strain, 1e-4);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_stress[0], reference_stress, 1e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(output_von_mises[0],reference_von_mises_pk2, 1.0e-4);
    }

    /**
    * Checks the ThermalLinearIsotropic3DWithTable
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementElement3D4NThermalLinearIsotropic3DWithTable, KratosConstitutiveLawsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THERMAL_EXPANSION_COEFFICIENT, 7.2e-6);
        p_elem_prop->SetValue(REFERENCE_TEMPERATURE, 0.0);

        Table<double> temp_E_table;
        temp_E_table.insert(-1.0, 2.0e6);
        temp_E_table.insert(3500.0, 2.1e4);
        temp_E_table.insert(1.0e6, 2.1e4);
        p_elem_prop->SetTable(TEMPERATURE, YOUNG_MODULUS, temp_E_table);

        TableAccessor table_accessor = TableAccessor(TEMPERATURE, "node_historical");
        p_elem_prop->SetAccessor(YOUNG_MODULUS, table_accessor.Clone());

        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("ThermalElasticIsotropic3D");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.1 , 0.1 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.5 , 0.2 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.9 , 0.4 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.3 , 0.3 , 0.5);

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }

        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementElement3D4N", 1, element_nodes, p_elem_prop);
        p_element->Initialize(r_model_part.GetProcessInfo());

        // Set a displacement and temperature
        array_1d<double, 3> aux_disp = ZeroVector(3);
        noalias(p_node_1->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        noalias(p_node_2->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        noalias(p_node_3->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        aux_disp(0) = 0.001;
        noalias(p_node_4->FastGetSolutionStepValue(DISPLACEMENT)) = aux_disp;
        p_node_1->FastGetSolutionStepValue(TEMPERATURE) = 2400.2;
        p_node_2->FastGetSolutionStepValue(TEMPERATURE) = 2400.2;
        p_node_3->FastGetSolutionStepValue(TEMPERATURE) = 2400.2;
        p_node_4->FastGetSolutionStepValue(TEMPERATURE) = 2400.2;

        Matrix lhs;
        Vector rhs;
        const auto& const_procinfo_ref = r_model_part.GetProcessInfo();
        p_element->InitializeSolutionStep(const_procinfo_ref);
        p_element->InitializeNonLinearIteration(const_procinfo_ref);
        p_element->CalculateLocalSystem(lhs,rhs,const_procinfo_ref);
        p_element->FinalizeNonLinearIteration(const_procinfo_ref);
        p_element->FinalizeSolutionStep(const_procinfo_ref);

        std::vector<Vector> output_strains(1);
        p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,output_strains, r_process_info);

        std::vector<Vector> output_stress(1);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR,output_stress, r_process_info);

        std::vector<double> output_von_mises(1);
        p_element->CalculateOnIntegrationPoints(VON_MISES_STRESS,output_von_mises, r_process_info);

        Vector reference_stress(6);
        reference_stress(0) = -27766.1;
        reference_stress(1) = -27766.1;
        reference_stress(2) = -27766.1;
        reference_stress(3) = 0.0;
        reference_stress(4) = 0.0;
        reference_stress(5) = 494.37;
        Vector reference_strain(6);
        reference_strain(0) = -0.0172814;
        reference_strain(1) = -0.0172814;
        reference_strain(2) = -0.0172814;
        reference_strain(3) = 0.0;
        reference_strain(4) = 0.0;
        reference_strain(5) = 0.002;
        const double reference_von_mises_pk2 =  856.274;

        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_strains[0], reference_strain, 1e-4);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(output_stress[0], reference_stress, 1e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(output_von_mises[0],reference_von_mises_pk2, 1.0e-4);
    }

} // namespace Kratos::Testing
