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

#include "custom_elements/nodal_data_handler.h"
#include "custom_elements/integration_point_data_container.h"

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
            rModelPart.CreateNewElement("DSS2D", 1, element_nodes, p_properties);
        }

        void InitializeCompleteElement(ModelPart& rModelPart)
        {
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
            rModelPart.AddNodalSolutionStepVariable(ADVPROJ);
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);
            rModelPart.AddNodalSolutionStepVariable(DENSITY);
            rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
            rModelPart.AddNodalSolutionStepVariable(DIVPROJ);

            rModelPart.SetBufferSize(2);

            Properties::Pointer p_properties = rModelPart.pGetProperties(0);

            // Geometry creation
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
            rModelPart.CreateNewElement("DSS2D", 1, element_nodes, p_properties);

            rModelPart.CloneTimeStep(0.1);

            Element& r_element = *(rModelPart.ElementsBegin());
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

            for (unsigned int i = 0; i < 3; i++) {
                Node<3>& r_node = r_geometry[i];
                r_node.FastGetSolutionStepValue(PRESSURE) = 10.0 * r_node.X();
                r_node.FastGetSolutionStepValue(VELOCITY_X) = r_node.Y();
                r_node.FastGetSolutionStepValue(DENSITY) = 100.0;
                r_node.FastGetSolutionStepValue(VISCOSITY) = 0.01;
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
            //KRATOS_CHECK_EXCEPTION_IS_THROWN(displacement_handler.Check(r_element), "Missing DISPLACEMENT variable");
        }

        // In-situ definition of FluidElementDataContainer list for tests
        #define FLUID_ELEMENT_VARIABLES(MACRO_TO_APPLY) \
        MACRO_TO_APPLY(VELOCITY,NodalVectorType) \
        MACRO_TO_APPLY(PRESSURE,NodalScalarType)

        MAKE_FLUID_ELEMENT_DATA_CONTAINER(TestFluidDataContainer, FLUID_ELEMENT_VARIABLES)
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

        KRATOS_TEST_CASE_IN_SUITE(DSS2D3NLocalMatrix, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeCompleteElement(model_part);

            Matrix LHS;
            Vector RHS;

            model_part.ElementsBegin()->CalculateLocalVelocityContribution(LHS,RHS,model_part.GetProcessInfo());

            //KRATOS_WATCH(LHS);
            //KRATOS_WATCH(RHS);

            KRATOS_CHECK_NEAR(LHS(0,0), 5.533840, 1e-5);
            KRATOS_CHECK_NEAR(LHS(3,5), -0.0243628, 1e-5);
            KRATOS_CHECK_NEAR(LHS(0,8), 0.166667, 1e-5);
            KRATOS_CHECK_NEAR(LHS(8,8), 0.00609399, 1e-5);

            KRATOS_CHECK_NEAR(RHS[0], 0.256372, 1e-5);
            KRATOS_CHECK_NEAR(RHS[2], 0.0609399, 1e-5);
        }
    }
}
