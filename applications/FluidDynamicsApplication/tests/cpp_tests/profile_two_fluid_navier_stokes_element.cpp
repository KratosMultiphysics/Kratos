//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Ruben Zorrilla
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "utilities/builtin_timer.h"


#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"
#include "processes/structured_mesh_generator_process.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

        namespace{
            void InitializeElements(ModelPart& TheModelPart){
                // Define the nodal values
                Matrix vel_original(4,3);
                vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
                vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
                vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
                vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

                Properties::Pointer pElemProp = TheModelPart.pGetProperties(0);

                // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
                for (NodeIteratorType it_node=TheModelPart.NodesBegin(); it_node<TheModelPart.NodesEnd(); ++it_node){
                    it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                    it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                    it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
                }
                
                for(auto& r_element : TheModelPart.Elements()){

                    for(unsigned int i=0; i<4; i++){
                        r_element.GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
                        for(unsigned int k=0; k<3; k++){
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                            r_element.GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                        }
                    }
                    r_element.GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) =  1.0;
                    r_element.GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
                    r_element.GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  1.0;
                    r_element.GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

               }
                
            }

        }

	    // /** Checks the TwoFluidNavierStokes3D4N element
	    //  * Checks the LHS and RHS for a cut element
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ProfileElementTwoFluidNavierStokesPositive3D4N, FluidDynamicsApplicationFastSuite)
		{
			Model current_model;
			ModelPart& modelPart = current_model.CreateModelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
            modelPart.AddNodalSolutionStepVariable(ACCELERATION);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
            Node<3>::Pointer p_point1 = Kratos::make_intrusive<Node<3>>(1, 0.00, 0.00, 0.00);
            Node<3>::Pointer p_point2 = Kratos::make_intrusive<Node<3>>(2, 10.00, 0.00, 0.00);
            Node<3>::Pointer p_point3 = Kratos::make_intrusive<Node<3>>(3, 10.00, 10.00, 0.00);
            Node<3>::Pointer p_point4 = Kratos::make_intrusive<Node<3>>(4, 0.00, 10.00, 0.00);
            Node<3>::Pointer p_point5 = Kratos::make_intrusive<Node<3>>(5, 0.00, 0.00, 10.00);
            Node<3>::Pointer p_point6 = Kratos::make_intrusive<Node<3>>(6, 10.00, 0.00, 10.00);
            Node<3>::Pointer p_point7 = Kratos::make_intrusive<Node<3>>(7, 10.00, 10.00, 10.00);
            Node<3>::Pointer p_point8 = Kratos::make_intrusive<Node<3>>(8, 0.00, 10.00, 10.00);

            Hexahedra3D8<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

            /////SET number_of_divisions to 100 for a realistic benchmark
            Parameters mesher_parameters(R"(
            {
                "number_of_divisions":   10,
                "element_name":     "TwoFluidNavierStokes3D4N"
            })");

		
            StructuredMeshGeneratorProcess(geometry, modelPart, mesher_parameters).Execute();
            
            for(auto& node : modelPart.Nodes()){
                node.AddDof(VELOCITY_X);
                node.AddDof(VELOCITY_Y);
                node.AddDof(VELOCITY_Z);
                node.AddDof(PRESSURE);

            }
            InitializeElements(modelPart);

            // Compute RHS and LHS
            const auto& const_procss_info_ref = modelPart.GetProcessInfo();

            block_for_each(modelPart.Elements(), [&](ModelPart::ElementType& rElement) {
                rElement.Initialize(const_procss_info_ref); // Initialize the element to initialize the constitutive law
            });
            
            block_for_each(modelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
                    rElement.Check(const_procss_info_ref);
            });

            KRATOS_WATCH(modelPart.NumberOfElements())
            for (int i = 0; i < 10; i++)
            {
                BuiltinTimer system_construction_time;      
#pragma omp parallel 
                {
                    Vector RHS = ZeroVector(16);
                    Matrix LHS = ZeroMatrix(16, 16);
#pragma omp for 
                    for (int i = 0; i < static_cast<int>(modelPart.NumberOfElements()); i++)
                    {
                        auto i_element = modelPart.ElementsBegin() + i;
                        i_element->CalculateLocalSystem(LHS, RHS, const_procss_info_ref);
                    }
                }
                std::cout << " time = " << system_construction_time.ElapsedSeconds() << std::endl;
            }
        }
    }
}
