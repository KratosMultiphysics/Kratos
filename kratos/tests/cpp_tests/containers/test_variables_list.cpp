//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "containers/variables_list_visitor.h"


namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariablesListHas, KratosCoreFastSuite) {
    VariablesList variables_list;
    variables_list.Add(NODAL_AREA);
    variables_list.Add(VELOCITY);


    KRATOS_CHECK(variables_list.Has(NODAL_AREA));
    KRATOS_CHECK(variables_list.Has(VELOCITY));
    KRATOS_CHECK_IS_FALSE(variables_list.Has(DISPLACEMENT));
}

KRATOS_TEST_CASE_IN_SUITE(VariablesListHasComponent, KratosCoreFastSuite) {
    VariablesList variables_list;
    variables_list.Add(DISPLACEMENT);
 
    KRATOS_CHECK_IS_FALSE(variables_list.Has(VELOCITY_X));
    KRATOS_CHECK_IS_FALSE(variables_list.Has(VELOCITY_Y));
    KRATOS_CHECK_IS_FALSE(variables_list.Has(VELOCITY_Z));

    KRATOS_CHECK(variables_list.Has(DISPLACEMENT_X));
    KRATOS_CHECK(variables_list.Has(DISPLACEMENT_Y));
    KRATOS_CHECK(variables_list.Has(DISPLACEMENT_Z));
}


KRATOS_TEST_CASE_IN_SUITE(VariablesListGetDofInfo, KratosCoreFastSuite) {
    VariablesList variables_list;
    auto dof_index = variables_list.AddDof(&DISPLACEMENT_Y, &REACTION_Y);
    KRATOS_CHECK_EQUAL(variables_list.GetDofVariable(dof_index), DISPLACEMENT_Y);
    KRATOS_CHECK_EQUAL(variables_list.pGetDofReaction(dof_index), &REACTION_Y);

}


    class SetModelPartVariableToOne{
            ModelPart& mrModelPart;
        public:
            SetModelPartVariableToOne(ModelPart& TheModelPart) : mrModelPart(TheModelPart){}

            template<typename TVariableType> 
            void operator()(TVariableType const& TheVariable) const {
                for (auto &node : mrModelPart.Nodes())
                {
                    auto& r_value = node.FastGetSolutionStepValue(TheVariable);
                    r_value = 1.0;    
                }                     
            }

            template<typename TDataType> 
            void operator()(Variable<array_1d<TDataType,3>> const& TheVariable) const {
                for (auto &node : mrModelPart.Nodes())
                {
                    auto& r_value = node.FastGetSolutionStepValue(TheVariable);
                    r_value = ScalarVector(3,-1.0);    
                }                     
            }

            void operator()(Variable<Vector> const& TheVaraible) const {} // Skip it for Vector

            void operator()(Variable<Matrix> const& TheVaraible) const {} // Skip it for Matrix
    };  

KRATOS_TEST_CASE_IN_SUITE(VariablesListApplyVisitor, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart &model_part = current_model.CreateModelPart("test");
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);  //not to have an empty var list
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT); //not to have an empty var list
    model_part.AddNodalSolutionStepVariable(VELOCITY);     //not to have an empty var list

    model_part.CreateNewNode(1, 1.0, 2.0, 3.0);
    model_part.CreateNewNode(2, 1.0, 2.0, 3.0);
    model_part.CreateNewNode(3, 1.0, 2.0, 3.0);

    //now create a container of pointers to variables of type VariablesData*
    std::vector<VariableData*> variables{&TEMPERATURE, &VELOCITY, &DISPLACEMENT_X, &DISPLACEMENT_Z};

   VariablesListVisitor<SetModelPartVariableToOne> the_visitor(model_part);

   the_visitor.VisitVariables(variables);

    for (auto& rNode : model_part.Nodes())
    {
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(TEMPERATURE), 1.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(DISPLACEMENT_X), 1.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(DISPLACEMENT_Y), 0.0); //NOT SET!!
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(DISPLACEMENT_Z), 1.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(VELOCITY_X), -1.0);
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(VELOCITY_Y), -1.0); 
        KRATOS_CHECK_EQUAL(rNode.FastGetSolutionStepValue(VELOCITY_Z), -1.0);

    }
}

}
} // namespace Kratos.
