//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    RAUL BRAVO
//

#if !defined( HROM_MODEL_PART_UTILITY_H_INCLUDED )
#define  HROM_MODEL_PART_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{

    // This utility creates a Hyper-Reduced Model part containing a selected set of elements and conditions
    class HromModelPartUtility
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(HromModelPartUtility);

        HromModelPartUtility(
        ModelPart& rModelPart,
        ModelPart& rModelPart2,
        Vector VectorOfElements,
        Vector VectorOfConditions): mpModelPart(rModelPart), mpModelPart2(rModelPart2)
        {
            Elements = VectorOfElements;
            Conditions = VectorOfConditions;
            //initializing commands
        }

        ~HromModelPartUtility()= default; //destructor

        void DoSomethig(){
            // for(int i=0; i<Elements.size();i++){
            //     KRATOS_WATCH(Elements(i))
            //     KRATOS_WATCH("\n\n\n")
            // }
            // for(int i=0; i<Conditions.size();i++){
            //     KRATOS_WATCH(Conditions(i))
            //     KRATOS_WATCH("\n\n\n")
            // }
            //ModelPart::NodesContainerType::iterator i_boundary_node
            // std::vector<ModelPart::IndexType> aux;
            // // //PointerVector aux;
            // for(auto& cond : mpModelPart.Conditions())
            // {
            //     auto& geom = cond.GetGeometry();
            //     for(unsigned int i=0; i<geom.size(); ++i)
            //     {
            //         aux.push_back(geom(i)->GetId());
            //     }
            // }
            //mpModelPart2.AddNodes(iterator.begin(), iterator.end())
            mpModelPart2.AddNodes(mpModelPart.NodesBegin(), mpModelPart.NodesEnd());

        }


        protected:
            //int some_stuff;  //variables
            ModelPart& mpModelPart;
            ModelPart& mpModelPart2;
            Vector Elements;
            Vector Conditions;
        };



} // namespace Kratos



#endif // HROM_MODEL_PART_UTILITY_H_INCLUDED  defined