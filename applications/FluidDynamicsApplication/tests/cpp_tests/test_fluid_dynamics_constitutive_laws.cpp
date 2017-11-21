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
#include "includes/properties.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

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

            // Create a raw model part
			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);

			// Set the element properties
			Properties::Pointer p_elem_prop = modelPart.pGetProperties(0);
			p_elem_prop->SetValue(DENSITY, 1000.0);
			p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 3.0e-01);
			p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

			// Element creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("EmbeddedNavierStokes2D3N", 1, elemNodes, p_elem_prop);
			Element::Pointer pElement = modelPart.pGetElement(1);

            // Set the constitutive law values
            ConstitutiveLaw::Parameters Values(pElement->GetGeometry(), pElement->GetProperties(), modelPart.GetProcessInfo());
        
            // Set the shape function values on an hypotetic centered Gauss pt.
            array_1d<double, 3> shape_func_values;
            shape_func_values(0) = 1.0/3.0;
            shape_func_values(1) = 1.0/3.0;
            shape_func_values(2) = 1.0/3.0;
            Values.SetShapeFunctionsValues(shape_func_values);

            // Set constitutive law flags:
            Flags& constitutive_law_options = Values.GetOptions();
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            // Set the constitutive arrays
            strain_vector(0) = 3.0;
            strain_vector(1) = 6.0;
            strain_vector(2) = 1.0;
            Values.SetStrainVector(strain_vector);  // Input strain values
            Values.SetStressVector(stress_vector);  // Output stress values
            Values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

            p_cons_law->CalculateMaterialResponseCauchy(Values);

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

	    // /** 
	    //  * Checks the Euler fluid 2D constitutive law.
	    //  */
	    // KRATOS_TEST_CASE_IN_SUITE(Euler2DConstitutiveLaw, FluidDynamicsApplicationFastSuite)
		// {
        //     // Declare the constitutive law pointer as well as its required arrays
        //     const unsigned int strain_size = 3;
        //     Newtonian2DLaw::Pointer p_cons_law(new Newtonian2DLaw());
        //     Vector stress_vector = ZeroVector(strain_size);
        //     Vector strain_vector = ZeroVector(strain_size);
        //     Matrix c_matrix = ZeroMatrix(strain_size, strain_size);

        //     // Create a raw model part
		// 	ModelPart modelPart("Main");
		// 	modelPart.SetBufferSize(3);

		// 	// Variables addition
		// 	modelPart.AddNodalSolutionStepVariable(DENSITY);
		// 	modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);

		// 	// Set the element properties
		// 	Properties::Pointer p_elem_prop = modelPart.pGetProperties(0);
		// 	p_elem_prop->SetValue(DENSITY, 1000.0);
		// 	p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 3.0e-01);
		// 	p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

		// 	// Element creation
		// 	modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
		// 	modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
		// 	modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
		// 	std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
		// 	modelPart.CreateNewElement("EmbeddedNavierStokes2D3N", 1, elemNodes, p_elem_prop);
		// 	Element::Pointer pElement = modelPart.pGetElement(1);

        //     // Set the constitutive law values
        //     ConstitutiveLaw::Parameters Values(pElement->GetGeometry(), pElement->GetProperties(), modelPart.GetProcessInfo());
        
        //     // Set the shape function values on an hypotetic centered Gauss pt.
        //     array_1d<double, 3> shape_func_values;
        //     shape_func_values(0) = 1.0/3.0;
        //     shape_func_values(1) = 1.0/3.0;
        //     shape_func_values(2) = 1.0/3.0;
        //     Values.SetShapeFunctionsValues(shape_func_values);

        //     // Set constitutive law flags:
        //     Flags& constitutive_law_options = Values.GetOptions();
        //     constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
        //     constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //     // Set the constitutive arrays
        //     strain_vector(0) = 3.0;
        //     strain_vector(1) = 6.0;
        //     strain_vector(2) = 1.0;
        //     Values.SetStrainVector(strain_vector);  // Input strain values
        //     Values.SetStressVector(stress_vector);  // Output stress values
        //     Values.SetConstitutiveMatrix(c_matrix); // Output constitutive tensor

        //     p_cons_law->CalculateMaterialResponseCauchy(Values);

        //     // Check computed values
        //     const double tolerance = 1e-10;

        //     KRATOS_CHECK_NEAR(c_matrix(0,0), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(0,1), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(0,2), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(1,0), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(1,1), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(1,2), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(2,0), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(2,1), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(c_matrix(2,2), 0.0, tolerance);

        //     KRATOS_CHECK_NEAR(stress_vector(0), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(stress_vector(1), 0.0, tolerance);
        //     KRATOS_CHECK_NEAR(stress_vector(2), 0.0, tolerance);
	    // }
	} // namespace Testing
}  // namespace Kratos.
