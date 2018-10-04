//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/euler_2d_law.h"
#include "custom_constitutive/euler_3d_law.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

namespace Kratos {
	namespace Testing {

        /** 
	     * Auxiliar function to generate a triangular element within 
         * a given model part using the constitutive law to be tested.
	     */
        void GenerateTriangle(
            ModelPart& rModelPart,
            const ConstitutiveLaw::Pointer pConstitutiveLaw) {

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);

            // Set the element properties
            Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 3.0e-01);
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConstitutiveLaw);

            // Element creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
            rModelPart.CreateNewElement("Element2D3N", 1, elem_nodes, p_elem_prop);
        }

        /** 
	     * Auxiliar function to generate a tetrahedral element within 
         * a given model part using the constitutive law to be tested.
	     */
        void GenerateTetrahedron(
            ModelPart& rModelPart,
            const ConstitutiveLaw::Pointer pConstitutiveLaw) {

            // Variables addition
            rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);

            // Set the element properties
            Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 3.0e-01);
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConstitutiveLaw);

            // Element creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
            rModelPart.CreateNewElement("Element3D4N", 1, elem_nodes, p_elem_prop);
        }

	    /** 
	     * Checks the Newtonian fluid 2D constitutive law.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(Newtonian2DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		{
            // Declare the constitutive law pointer as well as its required arrays
            const unsigned int strain_size = 3;
            Newtonian2DLaw::Pointer p_cons_law(new Newtonian2DLaw());
            Vector stress_vector = ZeroVector(strain_size);
            Vector strain_vector = ZeroVector(strain_size);
            Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

            // Get the trial element
			ModelPart model_part("Main", 3);
            GenerateTriangle(model_part, p_cons_law);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Set the constitutive law values
            ConstitutiveLaw::Parameters cons_law_values(
                p_element->GetGeometry(), 
                p_element->GetProperties(), 
                model_part.GetProcessInfo());

            // Set constitutive law flags:
            Flags& constitutive_law_options = cons_law_values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            cons_law_values.SetStrainVector(strain_vector);  // Input strain values
            cons_law_values.SetStressVector(stress_vector);  // Output stress values
            cons_law_values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

            // Check computed values
            const double tolerance = 1e-10;

            KRATOS_CHECK_NEAR(c_matrix(0,0),  0.4, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,1), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,0), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,1),  0.4, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,2),  0.0, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,0),  0.0, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,1),  0.0, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,2),  0.3, tolerance);

            KRATOS_CHECK_NEAR(stress_vector(0), 0.0, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(1), 1.8, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(2), 0.3, tolerance);
	    }

	    /** 
	     * Checks the Newtonian fluid 3D constitutive law.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(Newtonian3DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		{
            // Declare the constitutive law pointer as well as its required arrays
            const unsigned int strain_size = 6;
            Newtonian3DLaw::Pointer p_cons_law(new Newtonian3DLaw());
            Vector stress_vector = ZeroVector(strain_size);
            Vector strain_vector = ZeroVector(strain_size);
            Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

            // Get the trial element
            ModelPart model_part("Main", 3);
            GenerateTetrahedron(model_part, p_cons_law);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Set the constitutive law values
            ConstitutiveLaw::Parameters cons_law_values(
                p_element->GetGeometry(), 
                p_element->GetProperties(), 
                model_part.GetProcessInfo());

            // Set constitutive law flags:
            Flags& constitutive_law_options = cons_law_values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            strain_vector(3) = 2.0;
            strain_vector(4) = 3.0;
            strain_vector(5) = 4.0;
            cons_law_values.SetStrainVector(strain_vector);  // Input strain values
            cons_law_values.SetStressVector(stress_vector);  // Output stress values
            cons_law_values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

            // Check computed values
            const double tolerance = 1e-10;

            KRATOS_CHECK_NEAR(c_matrix(0,0),  0.4, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,1), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,2), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,0), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,1),  0.4, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,2), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,0), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,1), -0.2, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,2),  0.4, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(3,3),  0.3, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(4,4),  0.3, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(5,5),  0.3, tolerance);

            KRATOS_CHECK_NEAR(stress_vector(0), -0.2, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(1),  1.6, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(2), -1.4, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(3),  0.6, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(4),  0.9, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(5),  1.2, tolerance);
	    }

        /** 
	     * Checks the Newtonian Two Fluid 3D constitutive law.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(NewtonianTwoFluid3DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		{
            // Declare the constitutive law pointer as well as its required arrays
            const unsigned int nnodes = 4;
            const unsigned int dim = 3;
            const unsigned int strain_size = (dim - 1) * 3;

            NewtonianTwoFluid3DLaw::Pointer p_cons_law(new NewtonianTwoFluid3DLaw());
            Vector stress_vector = ZeroVector(strain_size);
            Vector strain_vector = ZeroVector(strain_size);
            Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

            // Get the trial element
            ModelPart model_part("Main", 3);
            model_part.AddNodalSolutionStepVariable(DISTANCE);
            model_part.AddNodalSolutionStepVariable(VELOCITY);
            model_part.AddNodalSolutionStepVariable(DENSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            GenerateTetrahedron(model_part, p_cons_law);
            Properties::Pointer p_elem_prop = model_part.pGetProperties(0);
            p_elem_prop->SetValue(C_SMAGORINSKY, 0.15);
            Element::Pointer p_element = model_part.pGetElement(1);

            // Set Nodal Values
            const array_1d<double,3> velocity (3, 0.2);  
            Geometry<Node<3>>& geom = p_element->GetGeometry();
            geom[0].GetSolutionStepValue(DISTANCE) = -1.0;
            geom[1].GetSolutionStepValue(DISTANCE) = -1.0;
            geom[2].GetSolutionStepValue(DISTANCE) = -1.0;
            geom[3].GetSolutionStepValue(DISTANCE) = 1.0;

            geom[0].GetSolutionStepValue(VELOCITY) = velocity;
            geom[1].GetSolutionStepValue(VELOCITY) = velocity;
            geom[2].GetSolutionStepValue(VELOCITY) = velocity;
            geom[3].GetSolutionStepValue(VELOCITY) = 2.0*velocity;

            geom[0].GetSolutionStepValue(DENSITY) = 2.0;
            geom[1].GetSolutionStepValue(DENSITY) = 2.0;
            geom[2].GetSolutionStepValue(DENSITY) = 2.0;
            geom[3].GetSolutionStepValue(DENSITY) = 0.5;

            geom[0].GetSolutionStepValue(DYNAMIC_VISCOSITY) = 0.3;
            geom[1].GetSolutionStepValue(DYNAMIC_VISCOSITY) = 0.3;
            geom[2].GetSolutionStepValue(DYNAMIC_VISCOSITY) = 0.3;
            geom[3].GetSolutionStepValue(DYNAMIC_VISCOSITY) = 0.2;

            // Set the constitutive law values
            ConstitutiveLaw::Parameters cons_law_values(
                p_element->GetGeometry(), 
                p_element->GetProperties(), 
                model_part.GetProcessInfo());

            // Set constitutive law flags:
            Flags& constitutive_law_options = cons_law_values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            strain_vector(3) = 2.0;
            strain_vector(4) = 3.0;
            strain_vector(5) = 4.0;
            cons_law_values.SetStrainVector(strain_vector);  // Input strain values
            cons_law_values.SetStressVector(stress_vector);  // Output stress values
            cons_law_values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            // Set Shape Functions Values
            array_1d<double, nnodes> N;            
            BoundedMatrix<double,nnodes, dim> DN_DX;
            double volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
            Vector N_copy = N;
            Matrix DN_DX_copy = DN_DX;
            cons_law_values.SetShapeFunctionsValues(N_copy);
            cons_law_values.SetShapeFunctionsDerivatives(DN_DX_copy);

            p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

            // Check computed values
            const double tolerance = 1e-7;

            KRATOS_CHECK_NEAR(c_matrix(0,0),  0.4100, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,1), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(0,2), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,0), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,1),  0.4100, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(1,2), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,0), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,1), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(2,2),  0.4100, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(3,3),  0.3075, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(4,4),  0.3075, tolerance);
            KRATOS_CHECK_NEAR(c_matrix(5,5),  0.3075, tolerance);

            KRATOS_CHECK_NEAR(stress_vector(0), -0.2050, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(1),  1.6400, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(2), -1.4350, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(3),  0.6150, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(4),  0.9225, tolerance);
            KRATOS_CHECK_NEAR(stress_vector(5),  1.2300, tolerance);
	    }

	    /** 
	     * Checks the Euler fluid 2D constitutive law.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(Euler2DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		{
            // Declare the constitutive law pointer as well as its required arrays
            const unsigned int strain_size = 3;
            Euler2DLaw::Pointer p_cons_law(new Euler2DLaw());
            Vector stress_vector = ZeroVector(strain_size);
            Vector strain_vector = ZeroVector(strain_size);
            Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

            // Create a raw model part
			ModelPart model_part("Main", 3);
            GenerateTriangle(model_part, p_cons_law);
			Element::Pointer p_element = model_part.pGetElement(1);

            // Set the constitutive law values
            ConstitutiveLaw::Parameters cons_law_values(
                p_element->GetGeometry(), 
                p_element->GetProperties(), 
                model_part.GetProcessInfo());

            // Set constitutive law flags:
            Flags& constitutive_law_options = cons_law_values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            cons_law_values.SetStrainVector(strain_vector);  // Input strain values
            cons_law_values.SetStressVector(stress_vector);  // Output stress values
            cons_law_values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

            // Check computed values
            const double tolerance = 1e-10;

            for (unsigned int i = 0; i < strain_size; ++i) {
                for (unsigned int j = 0; j < strain_size; ++j) {
                    KRATOS_CHECK_NEAR(c_matrix(i,j), 0.0, tolerance);
                }
            }

            for (unsigned int i = 0; i < strain_size; ++i) {
                KRATOS_CHECK_NEAR(stress_vector(i), 0.0, tolerance);
            }
	    }

	    /** 
	     * Checks the Euler fluid 3D constitutive law.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(Euler3DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		{
            // Declare the constitutive law pointer as well as its required arrays
            const unsigned int strain_size = 6;
            Euler3DLaw::Pointer p_cons_law(new Euler3DLaw());
            Vector stress_vector = ZeroVector(strain_size);
            Vector strain_vector = ZeroVector(strain_size);
            Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

            // Create a raw model part
			ModelPart model_part("Main", 3);
            GenerateTetrahedron(model_part, p_cons_law);
			Element::Pointer p_element = model_part.pGetElement(1);

            // Set the constitutive law values
            ConstitutiveLaw::Parameters cons_law_values(
                p_element->GetGeometry(), 
                p_element->GetProperties(), 
                model_part.GetProcessInfo());

            // Set constitutive law flags:
            Flags& constitutive_law_options = cons_law_values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            strain_vector(3) = 2.0;
            strain_vector(4) = 3.0;
            strain_vector(5) = 4.0;
            cons_law_values.SetStrainVector(strain_vector);  // Input strain values
            cons_law_values.SetStressVector(stress_vector);  // Output stress values
            cons_law_values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            p_cons_law->CalculateMaterialResponseCauchy(cons_law_values);

            // Check computed values
            const double tolerance = 1e-10;

            for (unsigned int i = 0; i < strain_size; ++i) {
                for (unsigned int j = 0; j < strain_size; ++j) {
                    KRATOS_CHECK_NEAR(c_matrix(i,j), 0.0, tolerance);
                }
            }

            for (unsigned int i = 0; i < strain_size; ++i) {
                KRATOS_CHECK_NEAR(stress_vector(i), 0.0, tolerance);
            }
	    }
	} // namespace Testing
}  // namespace Kratos.
