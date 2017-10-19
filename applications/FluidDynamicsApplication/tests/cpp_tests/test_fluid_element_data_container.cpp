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
			rModelPart.CreateNewElement("VMS2D", 1, element_nodes, p_properties);
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

			NodalDataHandler<double, 3, array_1d<double, 3>> PressureHandler(PRESSURE);
			NodalDataHandler<array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,3,2>> VelocityHandler(VELOCITY);

			PressureHandler.Initialize(r_element, model_part.GetProcessInfo());
            VelocityHandler.Initialize(r_element, model_part.GetProcessInfo());

            boost::numeric::ublas::matrix_row< Matrix > shape_functions = row(NContainer,0);
			KRATOS_CHECK_NEAR(2.0, PressureHandler.Interpolate(shape_functions, &r_element), 1e-6);

            array_1d<double,3> velocity = VelocityHandler.Interpolate(shape_functions, &r_element);
            KRATOS_CHECK_NEAR(0.0, velocity[0], 1e-6);
			KRATOS_CHECK_NEAR(3.0, velocity[1], 1e-6);
            KRATOS_CHECK_NEAR(0.0, velocity[2], 1e-6); // Note: velocity Z is not stored in the 2D handler, so it should return 0
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
			KRATOS_CHECK_NEAR(2.0, DataPoint.GetPRESSURE().Interpolate(shape_functions, &r_element), 1e-6);
        }
    }
}
