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
//
//

// System includes
#include <vector>
#include <algorithm>
// External includes

// Project includes
#include "processes/reorder_and_optimize_modelpart_process.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{
	ReorderAndOptimizeModelPartProcess::ReorderAndOptimizeModelPartProcess(ModelPart& rModelPart, Parameters settings)
		: Process()
		, mrModelPart(rModelPart.GetRootModelPart())
	{

		Parameters default_parameters(R"(
                {
                }
                )");


	}


	void ReorderAndOptimizeModelPartProcess::Execute()
    {
        KRATOS_TRY
        //reorder nodes,#include "spaces/ublas_space.h"elements and conditions so that their Id start in 1 and is consecutive
        IndexPartition<std::size_t>(mrModelPart.Nodes().size()).for_each([&](std::size_t Index){
            (mrModelPart.NodesBegin() + Index)->SetId(Index+1);
        });
        mrModelPart.Nodes().Sort();

        IndexPartition<std::size_t>(mrModelPart.Elements().size()).for_each([&](std::size_t Index){
            (mrModelPart.ElementsBegin() + Index)->SetId(Index+1);
        });
        mrModelPart.Elements().Sort();

        IndexPartition<std::size_t>(mrModelPart.Conditions().size()).for_each([&](std::size_t Index){
            (mrModelPart.ConditionsBegin() + Index)->SetId(Index+1);
        });
        mrModelPart.Conditions().Sort();

        OptimizeOrdering();

        //make a parallel clone of all the nodes
        IndexPartition<std::size_t>(mrModelPart.Nodes().size()).for_each([&](std::size_t Index){
            Node<3>::Pointer& pnode = *(mrModelPart.Nodes().ptr_begin() + Index);
            pnode = pnode->Clone();
        });

        IndexPartition<std::size_t>(mrModelPart.Elements().size()).for_each([&](std::size_t Index){
            Element::Pointer& pelem = *(mrModelPart.Elements().ptr_begin() + Index);
            Geometry<Node<3>>::PointsArrayType tmp;
            const auto& geom = pelem->GetGeometry();
            tmp.reserve(geom.size());
            for(unsigned int k=0; k<geom.size(); ++k)
                tmp.push_back( mrModelPart.Nodes()(geom[k].Id()) );

            auto paux = pelem->Create(pelem->Id(), tmp, pelem->pGetProperties());

            paux->Data() = pelem->Data();

            pelem = paux;
        });

        IndexPartition<std::size_t>(mrModelPart.Conditions().size()).for_each([&](std::size_t Index){
            Condition::Pointer& pcond = *(mrModelPart.Conditions().ptr_begin() + Index);
            Geometry<Node<3>>::PointsArrayType tmp;
            const auto& geom = pcond->GetGeometry();
            tmp.reserve(geom.size());
            for(unsigned int k=0; k<geom.size(); ++k)
                tmp.push_back( mrModelPart.Nodes()(geom[k].Id()));

            auto paux = pcond->Create(pcond->Id(), tmp, pcond->pGetProperties());

            paux->Data() = pcond->Data();

            pcond = paux;
        });

        //actualize pointers within submodelparts
        for(auto& subpart : mrModelPart.SubModelParts())
        {
            ActualizeSubModelPart(subpart);
        }

            //here i do a check
//             for(auto& subpart : mrModelPart.SubModelParts())
//             {
//                 for(const auto& node: subpart.Nodes() )
//                     if(node != mrModelPart.Nodes()[node.Id()])
//                         KRATOS_ERROR << "node in submodelpart with Id " << node.Id() << " is not in the main model part" << std::endl;
//
//                 for(auto it = subpart.Elements().ptr_begin(); it!=subpart.Elements().ptr_end(); it++)
//                 {
//                     for(unsigned int k=0; k<(*it)->GetGeometry().size(); k++)
//                         if(&(*it)->GetGeometry()[k] != &mrModelPart.Nodes()[(*it)->GetGeometry()[k].Id()])
//                             KRATOS_ERROR << *it << " node in element " << (*it)->GetGeometry()[k].Id() << " is not in the main model part" << std::endl;
//                 }
//
//                 for(auto it = subpart.Conditions().ptr_begin(); it!=subpart.Conditions().ptr_end(); it++)
//                 {
//                     for(unsigned int k=0; k<(*it)->GetGeometry().size(); k++)
//                     {
//                         if(&(*it)->GetGeometry()[k] != &mrModelPart.Nodes()[(*it)->GetGeometry()[k].Id()])
//                         {
//
//                             std::cout << "problematic condition Id = " << (*it)->Id() <<std::endl;
//                             std::cout << "problematic condition pointer " << *it << std::endl;
//                             std::cout << "pointer in main model part    " << mrModelPart.Conditions()((*it)->Id());
//                             std::cout << "problematic condition geometry = " << (*it)->GetGeometry() << std::endl;
//                             KRATOS_ERROR << (*it) << "node in condition " << (*it)->GetGeometry()[k].Id() << " is not in the main model part" << std::endl;
//                         }
//                     }
//                 }
//             }





            KRATOS_CATCH("")

        }


        void ReorderAndOptimizeModelPartProcess::ActualizeSubModelPart(ModelPart& subpart)
        {
            //make a parallel clone of all the nodes
            IndexPartition<std::size_t>(subpart.Nodes().size()).for_each([&](std::size_t Index){
                Node<3>::Pointer& pnode = *(subpart.NodesBegin() + Index).base();
                pnode = mrModelPart.Nodes()(pnode->Id());
            });
            subpart.Nodes().Sort();

            IndexPartition<std::size_t>(subpart.Elements().size()).for_each([&](std::size_t Index){
                Element::Pointer& pelem = *(subpart.ElementsBegin() + Index).base();
                pelem = mrModelPart.Elements()(pelem->Id());
            });
            subpart.Elements().Sort();

            IndexPartition<std::size_t>(subpart.Conditions().size()).for_each([&](std::size_t Index){
                Condition::Pointer& pcond = *(subpart.ConditionsBegin() + Index).base();
                pcond = mrModelPart.Conditions()(pcond->Id());
            });

             //actualize pointers within submodelparts
            for(auto& part : subpart.SubModelParts())
            {
                ActualizeSubModelPart(part);
            }
        }

        void ReorderAndOptimizeModelPartProcess::OptimizeOrdering()
        {
            const unsigned int n = mrModelPart.Nodes().size();
            std::vector< std::set< std::size_t > > graph(n);

            for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i)
            {
                const auto& elem = mrModelPart.ElementsBegin() + i;
                const auto& geom = elem->GetGeometry();
                for(unsigned int k=0; k<geom.size(); ++k)
                    for(unsigned int l=0; l<geom.size(); ++l)
                        graph[geom[k].Id()-1].insert(geom[l].Id()-1); //the -1 is because the ids of our nodes start in 1
            }

            for(int i=0; i<static_cast<int>(mrModelPart.Conditions().size()); ++i)
            {
                const auto& cond = mrModelPart.ConditionsBegin() + i;
                const auto& geom = cond->GetGeometry();
                for(unsigned int k=0; k<geom.size(); ++k)
                    for(unsigned int l=0; l<geom.size(); ++l)
                        graph[geom[k].Id()-1].insert(geom[l].Id()-1); //the -1 is because the ids of our nodes start in 1
            }


            unsigned int nnz = 0; //number of nonzeros in the graph
            for(unsigned int i=0; i<n; ++i)
            {
                //if(graph[i].size() == 0) KRATOS_ERROR << "node with Id " << i+1 << "has zero neighbours " << std::endl;
                nnz += graph[i].size();
            }


            CompressedMatrix graph_csr(n,n);
            for(unsigned int i=0; i<n; ++i)
                for(const auto& k : graph[i])
                    graph_csr.push_back(i,k,1);


            //here must do the ordering!!
            std::vector<int> invperm(graph_csr.size1());
            CuthillMcKee<false>().get<CompressedMatrix>(graph_csr,invperm);

            IndexPartition<std::size_t>(mrModelPart.Nodes().size()).for_each([&](std::size_t Index){
                (mrModelPart.NodesBegin() + Index)->SetId(invperm[Index]+1);
            });

            //reorder
            mrModelPart.Nodes().Sort();

            ReorderElements();
        }

        void ReorderAndOptimizeModelPartProcess::ReorderElements()
        {
            // Expects element ids are ordered 1 ... Elements().size().
            std::vector<std::size_t> element_ids(mrModelPart.NumberOfElements());
            std::vector<std::size_t> node_ids(mrModelPart.NumberOfElements());

            block_for_each(mrModelPart.Elements(), [&](Element& rElem){
                element_ids.at(rElem.Id() - 1) = rElem.Id();
                std::size_t node_id = rElem.GetGeometry()[0].Id();
                for (const auto& r_node : rElem.GetGeometry().Points())
                    node_id = std::min(node_id, r_node.Id());
                node_ids.at(rElem.Id() - 1) = node_id;
            });

            std::stable_sort(element_ids.begin(), element_ids.end(),
                      [&node_ids](const std::size_t& i, const std::size_t& j) {
                          return node_ids[i - 1] < node_ids[j - 1];
                      });

            IndexPartition<std::size_t>(element_ids.size()).for_each([&](std::size_t Index){
                (mrModelPart.ElementsBegin() + element_ids[Index] - 1)->SetId(Index + 1);
            });
            mrModelPart.Elements().Sort();
        }

        std::string ReorderAndOptimizeModelPartProcess::Info() const
        {
            return "ReorderAndOptimizeModelPartProcess";
        }

        void ReorderAndOptimizeModelPartProcess::PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
        }

        void ReorderAndOptimizeModelPartProcess::PrintData(std::ostream& rOStream) const
        {
        }

}  // namespace Kratos.
