//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

#include <type_traits>

#include "testing/testing.h"
#include "containers/generic_variables_list.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/mpi_serializer.h"

namespace Kratos
{
namespace Testing
{

    class MyTestVisitor : public GenericVariablesList::visitor_base_type
    {
    public:
        //here i pass in the constructor a ModelPart reference...but it is NOt needed
        MyTestVisitor(ModelPart &rModelPart) : mrModelPart(rModelPart)
        {
        }


        //this method is good for both Variable<double> and compoents
        //this template ACCEPTS ANY SCALAR TYPE
        template <class T>
        typename std::enable_if<std::is_scalar<typename T::Type>::value, void>::type 
        operator()(const T& rVar)
        {
            KRATOS_WATCH(rVar.Name())
            for (auto &node : mrModelPart.Nodes())
                {
                    auto& r_value = node.FastGetSolutionStepValue(rVar);
                    r_value = 1.0;
                }
        }

        //THIS OVERLOAD IS NEEDED TO ACCEPT NON-SCALAR TYPES
        template <class T>
        typename std::enable_if<!(std::is_scalar<typename T::Type>::value), void>::type 
        operator()(const T& rVar)
        {
            KRATOS_ERROR << " dummy operator overload for variable : " << rVar.Name() << std::endl;
        }

        //treat specially the case of array_1d variable
        void operator()(const Variable < array_1d<double, 3> >& rVar)
        {
            for (auto &node : mrModelPart.Nodes())
            {
                auto& r_value = node.FastGetSolutionStepValue(rVar);
                for (unsigned int i = 0; i < 3; ++i)
                    r_value[i] = -1.0;
            }
        }

    private:
        ModelPart &mrModelPart;
    };


KRATOS_TEST_CASE_IN_SUITE(GenericVariablesList, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart &mp = current_model.CreateModelPart("test");
    mp.AddNodalSolutionStepVariable(TEMPERATURE);  //not to have an empty var list
    mp.AddNodalSolutionStepVariable(DISPLACEMENT); //not to have an empty var list
    mp.AddNodalSolutionStepVariable(VELOCITY);     //not to have an empty var list

    mp.CreateNewNode(1, 1.0, 2.0, 3.0);
    mp.CreateNewNode(2, 1.0, 2.0, 3.0);
    mp.CreateNewNode(3, 1.0, 2.0, 3.0);




    //now create a generic list of variables
    GenericVariablesList vars;
    vars.AddVariable(TEMPERATURE);
    vars.AddVariable(VELOCITY);
    vars.AddVariable(DISPLACEMENT_X);
    vars.AddVariable(std::string("DISPLACEMENT_Z")); //note that we can add by the name

    //define the functor to be applied onto all the vars types
    //NOTE: this is the part to be user-defined!!
    MyTestVisitor aux(mp);
    vars.ApplyVisitor( aux );
   

    for (auto& rNode : mp.Nodes())
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

} // namespace Testing
} // namespace Kratos
