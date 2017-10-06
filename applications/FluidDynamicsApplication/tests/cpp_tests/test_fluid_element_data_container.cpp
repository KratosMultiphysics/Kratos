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

namespace Kratos {
    namespace Testing {

        void InitializeTestModelpart(ModelPart& rModelPart)
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
            InitializeTestModelpart(model_part);

            Element& r_element = *(model_part.ElementsBegin());
			Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

			for (unsigned int i = 0; i < 3; i++) {
				r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 1.0 + i;
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y) = 3.0 * i;
            }

            Matrix NContainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

			NodalDataHandler<double, 3, array_1d<double, 3>> PressureHandler(PRESSURE);
			NodalDataHandler<array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,2,3>> VelocityHandler(VELOCITY);

			PressureHandler.Initialize(r_element, model_part.GetProcessInfo());

            boost::numeric::ublas::matrix_row< Matrix > shape_functions = row(NContainer,0);
			KRATOS_CHECK_NEAR(2.0, PressureHandler.Interpolate(shape_functions, &r_element), 1e-6);

            array_1d<double,3> velocity = VelocityHandler.Interpolate(shape_functions, &r_element);
            KRATOS_CHECK_NEAR(0.0, velocity[0], 1e-6);
			KRATOS_CHECK_NEAR(3.0, velocity[1], 1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(FluidElementGaussPointData, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeTestModelpart(model_part);

            Element& r_element = *(model_part.ElementsBegin());
			Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

            // Data container tests go here
        }
    }
}
