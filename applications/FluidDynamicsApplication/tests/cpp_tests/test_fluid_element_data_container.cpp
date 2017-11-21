//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

#include "custom_utilities/nodal_data_handler.h"
#include "custom_utilities/fluid_element_data_container.h"

#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos {
    namespace Testing {

        void InitializeTestModelPart(ModelPart& rModelPart)
        {
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);

            Properties::Pointer p_properties = rModelPart.pGetProperties(0);

            // Geometry creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
            rModelPart.CreateNewElement("QSVMS2D3N", 1, element_nodes, p_properties);
        }

        void InitializeCompleteElement(ModelPart& rModelPart, const std::string& rElementName, unsigned int BufferSize)
        {
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            rModelPart.AddNodalSolutionStepVariable(ADVPROJ); // OSS only!
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);
            rModelPart.AddNodalSolutionStepVariable(DENSITY);
            rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            rModelPart.AddNodalSolutionStepVariable(DIVPROJ); // OSS only!

            rModelPart.SetBufferSize(BufferSize);

            Properties::Pointer p_properties = rModelPart.pGetProperties(0);

            // Geometry creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
            rModelPart.CreateNewElement(rElementName, 1, element_nodes, p_properties);

            const double dt = 0.1;
            // Loop starts at 1 because you need one less clone than time steps (JC)
            for (unsigned int i = 1; i < BufferSize; i++) {
                rModelPart.CloneTimeStep(i*dt);
            }

            Element& r_element = *(rModelPart.ElementsBegin());
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

            for (unsigned int i = 0; i < 3; i++) {
                Node<3>& r_node = r_geometry[i];
                r_node.FastGetSolutionStepValue(PRESSURE) = 10.0 * r_node.X();
                r_node.FastGetSolutionStepValue(VELOCITY_X) = r_node.Y();
                r_node.FastGetSolutionStepValue(DENSITY) = 100.0;
                r_node.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = 1.0;
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(FluidElementNodalDataHandler, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeTestModelPart(model_part);

            Element& r_element = *(model_part.ElementsBegin());
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

            for (unsigned int i = 0; i < 3; i++) {
                r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 1.0 + i;
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y) = 3.0 * i;
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_Z) = i - 5.0;
            }

            Matrix NContainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

            NodalDataHandler<double, 3, array_1d<double, 3>> pressure_handler(PRESSURE);
            NodalDataHandler<array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,3,2>> velocity_handler(VELOCITY);

			pressure_handler.Initialize(r_element, model_part.GetProcessInfo());
            velocity_handler.Initialize(r_element, model_part.GetProcessInfo());

            boost::numeric::ublas::matrix_row< Matrix > shape_functions = row(NContainer,0);
            KRATOS_CHECK_NEAR(2.0, pressure_handler.Interpolate(shape_functions, &r_element), 1e-6);

            array_1d<double,3> velocity = velocity_handler.Interpolate(shape_functions, &r_element);
            KRATOS_CHECK_NEAR(0.0, velocity[0], 1e-6);
            KRATOS_CHECK_NEAR(3.0, velocity[1], 1e-6);
            KRATOS_CHECK_NEAR(0.0, velocity[2], 1e-6); // Note: velocity Z is not stored in the 2D handler, so it should return 0

            // Test Check method: success case
            KRATOS_CHECK_EQUAL(pressure_handler.Check(r_element), 0);
            
            // Test Check method: failure case
            NodalDataHandler<array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,3,2>> displacement_handler(DISPLACEMENT);
            KRATOS_CHECK_EXCEPTION_IS_THROWN(displacement_handler.Check(r_element), "Missing DISPLACEMENT variable");
        }

        // In-situ definition of FluidElementDataContainer list for tests
        #define FLUID_ELEMENT_VARIABLES(MACRO_TO_APPLY) \
        MACRO_TO_APPLY(VELOCITY,NodalVector) \
        MACRO_TO_APPLY(PRESSURE,NodalScalar)

        MAKE_FLUID_ELEMENT_DATA_CONTAINER(TestFluidDataContainer, 2, 3, FLUID_ELEMENT_VARIABLES)
        #undef FLUID_ELEMENT_VARIABLES

        KRATOS_TEST_CASE_IN_SUITE(FluidElementGaussPointData, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeTestModelPart(model_part);

            Element& r_element = *(model_part.ElementsBegin());
			Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

            for (unsigned int i = 0; i < 3; i++) {
                r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 1.0 + i;
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y) = 3.0 * i;
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_Z) = i - 5.0;
            }

            // Data container tests go here
            TestFluidDataContainer DataPoint;

            DataPoint.Initialize(r_element,model_part.GetProcessInfo());

            Matrix NContainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
            boost::numeric::ublas::matrix_row< Matrix > shape_functions = row(NContainer,0);
            double interpolated_pressure = DataPoint.GetPRESSURE().Interpolate(shape_functions, &r_element);
            KRATOS_CHECK_NEAR(2.0, interpolated_pressure, 1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(LocalMatrixQSVMS2D3N, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeCompleteElement(model_part,"QSVMS2D3N",2);

            Matrix LHS;
            Vector RHS;

            model_part.ElementsBegin()->CalculateLocalVelocityContribution(LHS,RHS,model_part.GetProcessInfo());

            KRATOS_WATCH(LHS);
            KRATOS_WATCH(RHS);

            KRATOS_CHECK_NEAR(LHS(0,0), 5.533840, 1e-5);
            KRATOS_CHECK_NEAR(LHS(3,5), -0.0243628, 1e-5);
            KRATOS_CHECK_NEAR(LHS(0,8), 0.166667, 1e-5);
            KRATOS_CHECK_NEAR(LHS(8,8), 0.00609399, 1e-5);

            KRATOS_CHECK_NEAR(RHS[0], 0.256372, 1e-5);
            KRATOS_CHECK_NEAR(RHS[2], 0.0609399, 1e-5);
        }
        
        KRATOS_TEST_CASE_IN_SUITE(LocalMatrixSymbolicNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeCompleteElement(model_part,"SymbolicNavierStokes2D3N",3);

            // Set the BDF coefficients (0.1 is the time step in InitializeCompleteElement)
            Vector BdfVector = ZeroVector(3);
            BdfVector[0] = 1.5/0.1;
            BdfVector[1] = -2./0.1;
            BdfVector[2] = 0.5/0.1;
            model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS,BdfVector);
            // Set sound velocity
			model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);

			// Set the element properties
            Node<3>& r_node = *(model_part.Nodes().begin());
			Properties::Pointer p_properties = model_part.pGetProperties(0);
			p_properties->SetValue(DENSITY, r_node.FastGetSolutionStepValue(DENSITY));
            p_properties->SetValue(DYNAMIC_VISCOSITY, r_node.FastGetSolutionStepValue(DYNAMIC_VISCOSITY));
			Newtonian2DLaw::Pointer p_constitutive(new Newtonian2DLaw());
			p_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive);

            Matrix LHS;
            Vector RHS;

            model_part.ElementsBegin()->Initialize();
            model_part.ElementsBegin()->CalculateLocalSystem(LHS,RHS,model_part.GetProcessInfo());

            KRATOS_WATCH(LHS);
            KRATOS_WATCH(RHS);

            KRATOS_CHECK_NEAR(LHS(0,0), 50.9043, 1e-5);
            KRATOS_CHECK_NEAR(LHS(3,5), -0.00938651, 1e-5);
            KRATOS_CHECK_NEAR(LHS(0,8), 0.166667, 1e-5);
            KRATOS_CHECK_NEAR(LHS(8,8), 0.00689309, 1e-5);

            KRATOS_CHECK_NEAR(RHS[0], 0.406135, 1e-5);
            KRATOS_CHECK_NEAR(RHS[2], 0.0689307, 1e-5);
        }
    }
}
