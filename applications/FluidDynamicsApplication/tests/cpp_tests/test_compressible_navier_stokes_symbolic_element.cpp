//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "custom_elements/compressible_navier_stokes.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

        /** Checks the CompressibleNavierStokes2D3N element.
		 * Checks the RHS in case of constant variables (no acceleration).
         * Test for Air at 20 Degrees and Mach 3
		 */
	    KRATOS_TEST_CASE_IN_SUITE(ElementCompressibleNavierStokes2D3NConstant, FluidDynamicsApplicationFastSuite)
		{
			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(MOMENTUM);
			modelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
			modelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);

			modelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

			// Process info creation
            double delta_time = 0.1;
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5/delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
            pElemProp->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
            pElemProp->SetValue(CONDUCTIVITY,  0.0257);
            pElemProp->SetValue(SPECIFIC_HEAT,  718);
            pElemProp->SetValue(HEAT_CAPACITY_RATIO, 1.4);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("CompressibleNavierStokes2D3N", 1, elemNodes, pElemProp);


			Element::Pointer pElement = modelPart.pGetElement(1);
            //std::cout<<*pElement<<std::endl; //Info on element

			// Define the nodal values
			array_1d<double, 3> momentum;
			double R = 287.05;						// J/(kg*K)
			double density = 1.772;                 // density
			double T = 0.00247;                         // Temperature in K
			double velocity = 2.9;
            momentum[0] = velocity*density;        //momentum
            momentum[1] = velocity*density;        //momentum
			momentum[2] = 0.0;        //momentum
			//double Ma = velocity/sqrt(1.4*R*T);
            double total_energy = density*R*T/(1.4-1)+(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2])/(2*density);

            // std::cout<<"\n\nSupersonic Test For Constant Variables"<<std::endl;
            // std::cout<<"\nDensity = "<<density<<std::endl;
			// std::cout<<"\nTemperature = "<<T<<std::endl;
			// std::cout<<"\nMomentum_0 = "<<momentum[0]<<std::endl;
            // std::cout<<"\nMomentum_1 = "<<momentum[1]<<std::endl;
			// std::cout<<"\nMomentum_2 = "<<momentum[2]<<std::endl;
            // std::cout<<"\nMach = "<<Ma<<std::endl;
            // std::cout<<"\nTotal energy:"<<total_energy<<"\n"<<std::endl;

            array_1d<double, 3> f_ext;
            f_ext[0] = 0.0;
            f_ext[1] = 0.0;
            f_ext[2] = 0.0;
            double r = 0.0;         // EXTERNAL_PRESSURE

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = density;
				it_node->FastGetSolutionStepValue(MOMENTUM) = momentum;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
				it_node->FastGetSolutionStepValue(BODY_FORCE) = f_ext;
				it_node->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
			}

			modelPart.CloneTimeStep(delta_time); //generate other time step equal to this one
            modelPart.CloneTimeStep(2*delta_time); //to simulate the zero acceleration

			// Compute RHS and LHS
			Vector RHS = ZeroVector(12);

			pElement->Initialize(); // Initialize the element
			pElement->CalculateRightHandSide(RHS, modelPart.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i){
				sum_RHS += RHS[i]*RHS[i];
                // std::cout<<"RHS["<<i<<"]="<<RHS[i]<<std::endl;
            }
            sum_RHS = sqrt(sum_RHS);
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-9);
        }

		/** Checks the CompressibleNavierStokes2D3N element.
		 * Checks the LHS and RHS stationary solid rigid movements.
         * Test for air at 20 degrees and Mach 3
		 */
        KRATOS_TEST_CASE_IN_SUITE(ElementCompressibleNavierStokes2D3NStationarySupersonic, FluidDynamicsApplicationFastSuite)
		{
			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(MOMENTUM);
			modelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
			modelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);

			modelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

			// Process info creation
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
            pElemProp->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
            pElemProp->SetValue(CONDUCTIVITY,  0.0257);
            pElemProp->SetValue(SPECIFIC_HEAT,  718);
            pElemProp->SetValue(HEAT_CAPACITY_RATIO, 1.4);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("CompressibleNavierStokes2D3N", 1, elemNodes, pElemProp);


			Element::Pointer pElement = modelPart.pGetElement(1);
//             std::cout<<*pElement<<std::endl; //Info on element

			// Define the nodal values
			array_1d<double, 3> momentum;
			array_1d<double, 3> momentum_n;
			double density = 1.772;              //density
			double T = 0.00247;                      // Temperature in K
			double R = 287.05;						// J/(kg*K)
			double velocity = 2.9;
            momentum[0] = velocity*density;        //momentum (velocity = 300 m/s)
            momentum[1] = velocity*density;        //momentum
			momentum[2] = 0.0;        //momentum
			//double Ma = velocity/sqrt(1.4*R*T);
            double total_energy = density*R*T/(1.4-1)+(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2])/(2*density);
            momentum_n[0] = 1.0; momentum_n[1] = 1.0; momentum_n[2] = 1.0; //Casual value for the previous steps

            // std::cout<<"\n\nSupersonic Test For Stationary Rigid Movements"<<std::endl;
            // std::cout<<"\n\nDensity = "<<density<<std::endl;
			// std::cout<<"\nTemperature = "<<T<<std::endl;
			// std::cout<<"\nMomentum_0 = "<<momentum[0]<<std::endl;
            // std::cout<<"\nMomentum_1 = "<<momentum[1]<<std::endl;
			// std::cout<<"\nMomentum_2 = "<<momentum[2]<<std::endl;
            // std::cout<<"\nMach = "<<Ma<<std::endl;
            // std::cout<<"\nTotal energy:"<<total_energy<<std::endl;

            array_1d<double, 3> f_ext;
            f_ext[0] = 0.0;
            f_ext[1] = 0.0;
            f_ext[2] = 0.0;
            double r = 0.0;         // EXTERNAL_PRESSURE

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = density;
				it_node->FastGetSolutionStepValue(DENSITY,1) = 1;
				it_node->FastGetSolutionStepValue(DENSITY,2) = 1;
				it_node->FastGetSolutionStepValue(MOMENTUM) = momentum;
				it_node->FastGetSolutionStepValue(MOMENTUM,1) = momentum_n;
				it_node->FastGetSolutionStepValue(MOMENTUM,2) = momentum_n;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY,1) = 1;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY,2) = 1;
				it_node->FastGetSolutionStepValue(BODY_FORCE) = f_ext;
				it_node->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(12);
			Matrix LHS = ZeroMatrix(12,12);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i)
				sum_RHS += RHS[i]*RHS[i];
            sum_RHS = sqrt(sum_RHS);
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-12);

			// Check rigid movement modes
			Vector a(12);
			Vector rhs(12);
			double sum_rhs;

            // Mode 1 check
            // std::cout<<"\nRigid Body Movement: Mode 1\n"<<std::endl;

			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 0.0;
            a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 0.0;
            a[8] = 1.0; a[9] = 0.0; a[10] = 0.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);

			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);

            // Mode 2 check
            // std::cout<<"\nRigid Body Movement: Mode 2\n"<<std::endl;

			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;
            a[4] = 0.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;
            a[8] = 0.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);

			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);

            // Mode 3 check
            // std::cout<<"\nRigid Body Movement: Mode 3\n"<<std::endl;

			a[0] = 0.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;
            a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;
            a[8] = 0.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);

            // Mode 4 check
            // std::cout<<"\nRigid Body Movement: Mode 4\n"<<std::endl;

			a[0] = 0.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 1.0;
            a[4] = 0.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 1.0;
            a[8] = 0.0; a[9] = 0.0; a[10] = 0.0; a[11] = 1.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);


	    }

	    /** Checks the CompressibleNavierStokes2D3N element.
		 * Checks the LHS and RHS stationary solid rigid movements.
         * Test for Air at 20 degrees and Mach 0.5
		 */
		KRATOS_TEST_CASE_IN_SUITE(ElementCompressibleNavierStokes2D3NStationarySubsonic, FluidDynamicsApplicationFastSuite)
		{
			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(MOMENTUM);
			modelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
            modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
            modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
			modelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);

			modelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

			// Process info creation
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
            pElemProp->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
            pElemProp->SetValue(CONDUCTIVITY,  0.0257);
            pElemProp->SetValue(SPECIFIC_HEAT,  718);
            pElemProp->SetValue(HEAT_CAPACITY_RATIO, 1.4);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("CompressibleNavierStokes2D3N", 1, elemNodes, pElemProp);


			Element::Pointer pElement = modelPart.pGetElement(1);
//             std::cout<<*pElement<<std::endl; //Info on element

			// Define the nodal values
			array_1d<double, 3> momentum;
			array_1d<double, 3> momentum_n;
			double density = 1.772;              //density
			double R = 287.05;						// J/(kg*K)
			double T = 293;                      // Temperature in K
			double velocity = 1.0;
            momentum[0] = velocity*density;        //momentum (velocity = 300 m/s)
            momentum[1] = velocity*density;        //momentum
			momentum[2] = 0.0;        //momentum
			//double Ma = velocity/sqrt(1.4*R*T);
            double total_energy = density*R*T/(1.4-1)+(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2])/(2*density);
            momentum_n[0] = 1.0; momentum_n[1] = 1.0; momentum_n[2] = 1.0; //Casual value for the previous steps

            // std::cout<<"\n\nSubsonic Test For Stationary Rigid Movements"<<std::endl;
            // std::cout<<"\n\nDensity = "<<density<<std::endl;
			// std::cout<<"\nTemperature = "<<T<<std::endl;
			// std::cout<<"\nMomentum_0 = "<<momentum[0]<<std::endl;
            // std::cout<<"\nMomentum_1 = "<<momentum[1]<<std::endl;
			// std::cout<<"\nMomentum_2 = "<<momentum[2]<<std::endl;
            // std::cout<<"\nMach = "<<Ma<<std::endl;
            // std::cout<<"\nTotal energy:"<<total_energy<<std::endl;

            array_1d<double, 3> f_ext;
            f_ext[0] = 0.0;
            f_ext[1] = 0.0;
            f_ext[2] = 0.0;
            double r = 0.0;         // EXTERNAL_PRESSURE

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = density;
				it_node->FastGetSolutionStepValue(DENSITY,1) = 1;
				it_node->FastGetSolutionStepValue(DENSITY,2) = 1;
				it_node->FastGetSolutionStepValue(MOMENTUM) = momentum;
				it_node->FastGetSolutionStepValue(MOMENTUM,1) = momentum_n;
				it_node->FastGetSolutionStepValue(MOMENTUM,2) = momentum_n;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY,1) = 1;
				it_node->FastGetSolutionStepValue(TOTAL_ENERGY,2) = 1;
				it_node->FastGetSolutionStepValue(BODY_FORCE) = f_ext;
				it_node->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(12);
			Matrix LHS = ZeroMatrix(12,12);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i)
				sum_RHS += RHS[i]*RHS[i];
            sum_RHS = sqrt(sum_RHS);
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-10);

			// Check rigid movement modes
			Vector a(12);
			Vector rhs(12);
			double sum_rhs;

            // Mode 1 check
            // std::cout<<"\nRigid Body Movement: Mode 1\n"<<std::endl;

			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 0.0;
            a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 0.0;
            a[8] = 1.0; a[9] = 0.0; a[10] = 0.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);

			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-10);

            // Mode 2 check
            // std::cout<<"\nRigid Body Movement: Mode 2\n"<<std::endl;

			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;
            a[4] = 0.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;
            a[8] = 0.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);

			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-10);

            // Mode 3 check
            // std::cout<<"\nRigid Body Movement: Mode 3\n"<<std::endl;

			a[0] = 0.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;
            a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;
            a[8] = 0.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-10);

            // Mode 4 check
            // std::cout<<"\nRigid Body Movement: Mode 4\n"<<std::endl;

			a[0] = 0.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 1.0;
            a[4] = 0.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 1.0;
            a[8] = 0.0; a[9] = 0.0; a[10] = 0.0; a[11] = 1.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i){
				sum_rhs += rhs[i];
                // std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-10);


	    }
		
		        /** Checks the CompressibleNavierStokes2D3N element.
		 * Checks the RHS in case of constant variables (no acceleration).
         * Test for Air at 20 Degrees and Mach 3
		 */
		 KRATOS_TEST_CASE_IN_SUITE(ElementCompressibleNavierStokes3D4NConstant, FluidDynamicsApplicationFastSuite)
		 {
			 ModelPart modelPart("Main");
			 modelPart.SetBufferSize(3);
 
			 // Variables addition
			 modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			 modelPart.AddNodalSolutionStepVariable(DENSITY);
			 modelPart.AddNodalSolutionStepVariable(MOMENTUM);
			 modelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
			 modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
			 modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			 modelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
			 modelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
			 modelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
 
			 modelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
 
			 // Process info creation
			 double delta_time = 0.1;
			 Vector bdf_coefs(3);
			 bdf_coefs[0] = 3.0/(2.0*delta_time);
			 bdf_coefs[1] = -2.0/delta_time;
			 bdf_coefs[2] = 0.5/delta_time;
			 modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
 
			 // Set the element properties
			 Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			 pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
			 pElemProp->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
			 pElemProp->SetValue(CONDUCTIVITY,  0.0257);
			 pElemProp->SetValue(SPECIFIC_HEAT,  718);
			 pElemProp->SetValue(HEAT_CAPACITY_RATIO, 1.4);
 
			 // Geometry creation
			 modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			 modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			 modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			 modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
			 std::vector<ModelPart::IndexType> elemNodes {1, 2, 3,4};
			 modelPart.CreateNewElement("CompressibleNavierStokes3D4N", 1, elemNodes, pElemProp);
			 
 
			 Element::Pointer pElement = modelPart.pGetElement(1);
			 //std::cout<<*pElement<<std::endl; //Info on element
			 
			 // Define the nodal values
			 array_1d<double, 3> momentum;
			 double R = 287.05;						// J/(kg*K)
			 double density = 1.772;                 // density
			 double T = 0.00247;                         // Temperature in K
			 double velocity = 2.9;   
			 momentum[0] = velocity*density;        //momentum
			 momentum[1] = velocity*density;        //momentum
			 momentum[2] = velocity*density;         //momentum
			 //double Ma = velocity/sqrt(1.4*R*T); 
			 double total_energy = density*R*T/(1.4-1)+(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2])/(2*density);   	
 
			//  std::cout<<"\n\nSupersonic Test For 3D Constant Variables"<<std::endl;
			//  std::cout<<"\nDensity = "<<density<<std::endl;
			//  std::cout<<"\nTemperature = "<<T<<std::endl;
			//  std::cout<<"\nMomentum_0 = "<<momentum[0]<<std::endl;
			//  std::cout<<"\nMomentum_1 = "<<momentum[1]<<std::endl;
			//  std::cout<<"\nMomentum_2 = "<<momentum[2]<<std::endl;
			//  std::cout<<"\nMach = "<<Ma<<std::endl;
			//  std::cout<<"\nTotal energy:"<<total_energy<<"\n"<<std::endl;
						 
			 array_1d<double, 3> f_ext;
			 f_ext[0] = 0.0;
			 f_ext[1] = 0.0;
			 f_ext[2] = 0.0;
			 double r = 0.0;         // EXTERNAL_PRESSURE
			 
			 // Set the nodal values
			 for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				 it_node->FastGetSolutionStepValue(DENSITY) = density;
				 it_node->FastGetSolutionStepValue(MOMENTUM) = momentum;
				 it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
				 it_node->FastGetSolutionStepValue(BODY_FORCE) = f_ext;
				 it_node->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
			 }
			 
			 modelPart.CloneTimeStep(delta_time); //generate other time step equal to this one
			 modelPart.CloneTimeStep(2*delta_time); //to simulate the zero acceleration 
 
			 // Compute RHS and LHS
			 Vector RHS = ZeroVector(20);
 
			 pElement->Initialize(); // Initialize the element 
			 pElement->CalculateRightHandSide(RHS, modelPart.GetProcessInfo());
 
			 // Check obtained RHS
			 double sum_RHS = 0.0;
			 for (unsigned int i=0; i<RHS.size(); ++i){
				 sum_RHS += RHS[i]*RHS[i];
				//  std::cout<<"RHS["<<i<<"]="<<RHS[i]<<std::endl;
			 }
			 sum_RHS = sqrt(sum_RHS);
			 KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-9);
		 }
		
		 /** Checks the CompressibleNavierStokes3D4N element.
		 * Checks the LHS and RHS stationary solid rigid movements.
         * Test for air at 20 degrees and Mach 3
		 */
		 KRATOS_TEST_CASE_IN_SUITE(ElementCompressibleNavierStokes3D4NStationarySupersonic, FluidDynamicsApplicationFastSuite)
		 {
			 ModelPart modelPart("Main");
			 modelPart.SetBufferSize(3);
 
			 // Variables addition
			 modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			 modelPart.AddNodalSolutionStepVariable(DENSITY);
			 modelPart.AddNodalSolutionStepVariable(MOMENTUM);
			 modelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
			 modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
			 modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			 modelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
			 modelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
			 modelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
 
			 modelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
 
			 // Process info creation
			 Vector bdf_coefs(3);
			 bdf_coefs[0] = 0.0;
			 bdf_coefs[1] = 0.0;
			 bdf_coefs[2] = 0.0;
			 modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
 
			 // Set the element properties
			 Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			 pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
			 pElemProp->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
			 pElemProp->SetValue(CONDUCTIVITY,  0.0257);
			 pElemProp->SetValue(SPECIFIC_HEAT,  718);
			 pElemProp->SetValue(HEAT_CAPACITY_RATIO, 1.4);
 
			 // Geometry creation
			 modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			 modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			 modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			 modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
			 std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			 modelPart.CreateNewElement("CompressibleNavierStokes3D4N", 1, elemNodes, pElemProp);
			 
 
			 Element::Pointer pElement = modelPart.pGetElement(1);
 //             std::cout<<*pElement<<std::endl; //Info on element
			 
			 // Define the nodal values
			 array_1d<double, 3> momentum;
			 array_1d<double, 3> momentum_n;
			 double density = 1.772;              //density
			 double T = 0.00247;                      // Temperature in K
			 double R = 287.05;						// J/(kg*K)
			 double velocity = 2.9;     
			 momentum[0] = velocity*density;        //momentum (velocity = 300 m/s)
			 momentum[1] = velocity*density;        //momentum
			 momentum[2] = 0.0;        //momentum
			 //double Ma = velocity/sqrt(1.4*R*T);
			 double total_energy = density*R*T/(1.4-1)+(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2])/(2*density);   	
			 momentum_n[0] = 1.0; momentum_n[1] = 1.0; momentum_n[2] = 1.0; //Casual value for the previous steps 
 
			//  std::cout<<"\n\nSupersonic Test For 3D Stationary Rigid Movements"<<std::endl;
			//  std::cout<<"\n\nDensity = "<<density<<std::endl;
			//  std::cout<<"\nTemperature = "<<T<<std::endl;
			//  std::cout<<"\nMomentum_0 = "<<momentum[0]<<std::endl;
			//  std::cout<<"\nMomentum_1 = "<<momentum[1]<<std::endl;
			//  std::cout<<"\nMomentum_2 = "<<momentum[2]<<std::endl;
			//  std::cout<<"\nMach = "<<Ma<<std::endl;
			//  std::cout<<"\nTotal energy:"<<total_energy<<std::endl;
						 
			 array_1d<double, 3> f_ext;
			 f_ext[0] = 0.0;
			 f_ext[1] = 0.0;
			 f_ext[2] = 0.0;
			 double r = 0.0;         // EXTERNAL_PRESSURE
			 
			 // Set the nodal values
			 for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				 it_node->FastGetSolutionStepValue(DENSITY) = density;
				 it_node->FastGetSolutionStepValue(DENSITY,1) = 1;
				 it_node->FastGetSolutionStepValue(DENSITY,2) = 1;
				 it_node->FastGetSolutionStepValue(MOMENTUM) = momentum;
				 it_node->FastGetSolutionStepValue(MOMENTUM,1) = momentum_n;
				 it_node->FastGetSolutionStepValue(MOMENTUM,2) = momentum_n;
				 it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
				 it_node->FastGetSolutionStepValue(TOTAL_ENERGY,1) = 1;
				 it_node->FastGetSolutionStepValue(TOTAL_ENERGY,2) = 1;
				 it_node->FastGetSolutionStepValue(BODY_FORCE) = f_ext;
				 it_node->FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
			 }
 
			 // Compute RHS and LHS
			 Vector RHS = ZeroVector(20);
			 Matrix LHS = ZeroMatrix(20,20);
 
			 pElement->Initialize(); // Initialize the element to initialize the constitutive law
			 pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());
 
			 // Check obtained RHS
			 double sum_RHS = 0.0;
			 for (unsigned int i=0; i<RHS.size(); ++i)
				 sum_RHS += RHS[i]*RHS[i];
			 sum_RHS = sqrt(sum_RHS);
			 KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-12);
			 
			 // Check rigid movement modes
			 Vector a(20);
			 Vector rhs(20);
			 double sum_rhs;
			 
			 // Mode 1 check
			//  std::cout<<"\nRigid Body Movement: Mode 1\n"<<std::endl;
			 
			 a[0] = 1.0; 	a[1] = 0.0; 	a[2] = 0.0;		a[3] = 0.0;		a[4] = 0.0;
			 a[5] = 1.0;	a[6] = 0.0;		a[7] = 0.0;		a[8] = 0.0; 	a[9] = 0.0;
			 a[10] = 1.0; 	a[11] = 0.0; 	a[12] = 0.0; 	a[13] = 0.0; 	a[14] = 0.0;
			 a[15] = 1.0;	a[16] = 0.0; 	a[17] = 0.0; 	a[18] = 0.0;    a[19] = 0.0;
			 
			 sum_rhs = 0.0;
			 rhs = prod(LHS,a);
			 
			 for (unsigned int i=0; i<rhs.size(); ++i){
				 sum_rhs += rhs[i];
				//  std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
			 }
			 KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);
			 
			 // Mode 2 check
			//  std::cout<<"\nRigid Body Movement: Mode 2\n"<<std::endl;
			 
			 a[0] = 0.0; 	a[1] = 1.0; 	a[2] = 0.0;		a[3] = 0.0;		a[4] = 0.0;
			 a[5] = 0.0;	a[6] = 1.0;		a[7] = 0.0;		a[8] = 0.0; 	a[9] = 0.0;
			 a[10] = 0.0; 	a[11] = 1.0; 	a[12] = 0.0; 	a[13] = 0.0; 	a[14] = 0.0;
			 a[15] = 0.0;	a[16] = 1.0; 	a[17] = 0.0; 	a[18] = 0.0;    a[19] = 0.0;
			 sum_rhs = 0.0;
			 rhs = prod(LHS,a);
			 
			 for (unsigned int i=0; i<rhs.size(); ++i){
				 sum_rhs += rhs[i];
				//  std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
			 }
			 KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);
			 
			 // Mode 3 check
			//  std::cout<<"\nRigid Body Movement: Mode 3\n"<<std::endl;
			 
			 a[0] = 0.0; 	a[1] = 0.0; 	a[2] = 1.0;		a[3] = 0.0;		a[4] = 0.0;
			 a[5] = 0.0;	a[6] = 0.0;		a[7] = 1.0;		a[8] = 0.0; 	a[9] = 0.0;
			 a[10] =0.0; 	a[11] = 0.0; 	a[12] = 1.0; 	a[13] = 0.0; 	a[14] = 0.0;
			 a[15] = 0.0;	a[16] = 0.0; 	a[17] = 1.0; 	a[18] = 0.0;    a[19] = 0.0;
			 sum_rhs = 0.0;
			 rhs = prod(LHS,a);
			 for (unsigned int i=0; i<rhs.size(); ++i){
				 sum_rhs += rhs[i];
				//  std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
			 }
			 KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);
			 
			 // Mode 4 check
			//  std::cout<<"\nRigid Body Movement: Mode 4\n"<<std::endl;
			 
			 a[0] = 0.0; 	a[1] = 0.0; 	a[2] = 0.0;		a[3] = 1.0;		a[4] = 0.0;
			 a[5] = 0.0;	a[6] = 0.0;		a[7] = 0.0;		a[8] = 1.0; 	a[9] = 0.0;
			 a[10] = 0.0; 	a[11] = 0.0; 	a[12] = 0.0; 	a[13] = 1.0; 	a[14] = 0.0;
			 a[15] = 0.0;	a[16] = 0.0; 	a[17] = 0.0; 	a[18] = 1.0;    a[19] = 0.0;
			 sum_rhs = 0.0;
			 rhs = prod(LHS,a);
			 for (unsigned int i=0; i<rhs.size(); ++i){
				 sum_rhs += rhs[i];
				//  std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
			 }
			 KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);

			 // Mode 5 check
			//  std::cout<<"\nRigid Body Movement: Mode 5\n"<<std::endl;
			 
			 a[0] = 0.0; 	a[1] = 0.0; 	a[2] = 0.0;		a[3] = 0.0;		a[4] = 1.0;
			 a[5] = 0.0;	a[6] = 0.0;		a[7] = 0.0;		a[8] = 0.0; 	a[9] = 1.0;
			 a[10] = 0.0; 	a[11] = 0.0; 	a[12] = 0.0; 	a[13] = 0.0; 	a[14] =1.0;
			 a[15] = 0.0;	a[16] = 0.0; 	a[17] = 0.0; 	a[18] = 0.0;    a[19] = 1.0;
			 sum_rhs = 0.0;
			 rhs = prod(LHS,a);
			 for (unsigned int i=0; i<rhs.size(); ++i){
				 sum_rhs += rhs[i];
				//  std::cout<<"\nRHS["<<i<<"]="<<rhs[i]<<std::endl;
			 }
			 KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-9);
			 	 
		}
    }
}
