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

#include "custom_elements/nodal_data_list.h"

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

        KRATOS_TEST_CASE_IN_SUITE(FluidElementGaussPointData, FluidDynamicsApplicationFastSuite)
        {
            ModelPart model_part("Test");
            InitializeTestModelpart(model_part);

            Element& r_element = *(model_part.ElementsBegin());

            //std::vector< DataHandler

            //NodalDataList element_data(r_element.GetGeometry(),DataHandlers);

            //NodalDataList< DataHandler< array_1d<double,3> >, DataHandler< array_1d<double,3> >, DataHandler<double> > elemental_data(r_element.GetGeometry());

        }
    }
}