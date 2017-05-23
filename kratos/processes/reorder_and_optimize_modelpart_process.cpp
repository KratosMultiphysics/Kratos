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
// External includes

// Project includes
#include "processes/reorder_and_optimize_modelpart_process.h"


namespace Kratos
{
	ReorderAndOptimizeModelPartProcess::ReorderAndOptimizeModelPartProcess(ModelPart& rModelPart, Parameters settings)
		: Process()
		, mrModelPart(rModelPart.GetRootModelPart()) {

		Parameters default_parameters(R"(
                {
                }  
                )");
            

	}

	ReorderAndOptimizeModelPartProcess::~ReorderAndOptimizeModelPartProcess() {
            
	}

	void ReorderAndOptimizeModelPartProcess::Execute() 
        {
            KRATOS_TRY

            //reorder nodes Ids so that they start in 1
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); ++i)
            {
                auto pnode = *(mrModelPart.NodesBegin() + i).base();
                pnode->SetId(i+1);
            }
            
            //optimize ordering
            
            //make a parallel clone of all the nodes
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); ++i)
            {
                Node<3>::Pointer& pnode = *(mrModelPart.Nodes().ptr_begin() + i);
                auto paux = (*pnode).Clone();
                pnode.swap(paux);
            }
            
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i)
            {
                Element::Pointer& pelem = *(mrModelPart.Elements().ptr_begin() + i);
                
                PointerVector< Node<3> > tmp;
                const auto& geom = pelem->GetGeometry();
                tmp.reserve(geom.size());
                for(unsigned int k=0; k<geom.size(); ++k)
                {
                    Node<3>::Pointer pn = mrModelPart.Nodes()(geom[k].Id());
                    tmp.push_back( pn );
                }
                auto paux = pelem->Create(pelem->Id(), tmp, pelem->pGetProperties());
                
                //TODO: copy elemental database
                
                pelem.swap(paux);
            }

            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(mrModelPart.Conditions().size()); ++i)
            {
                Condition::Pointer pcond = *(mrModelPart.Conditions().ptr_begin() + i);
                
                PointerVector< Node<3> > tmp;
                const auto& geom = pcond->GetGeometry();
                tmp.reserve(geom.size());
                for(unsigned int k=0; k<geom.size(); ++k)
                    tmp.push_back( mrModelPart.Nodes()(geom[k].Id()));
                
                auto paux = pcond->Create(pcond->Id(), tmp, pcond->pGetProperties());
                
                //TODO: copy elemental database
                
                pcond.swap(paux);
            }
            
            //actualize pointers within submodelparts
            for(auto& subpart : mrModelPart.SubModelParts())
            {
                ActualizeSubModelPart(subpart);
            }
            
            KRATOS_CATCH("")
            
        }
        
        
        void ReorderAndOptimizeModelPartProcess::ActualizeSubModelPart(ModelPart& subpart)
        {
            //make a parallel clone of all the nodes
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(subpart.Nodes().size()); ++i)
            {
                Node<3>::Pointer& pnode = *(subpart.NodesBegin() + i).base();
                auto paux = mrModelPart.Nodes()(pnode->Id());
                pnode.swap(paux);
            }
            
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(subpart.Elements().size()); ++i)
            {
                Element::Pointer& pelem = *(subpart.ElementsBegin() + i).base();              
                auto paux = mrModelPart.Elements()(pelem->Id());
                pelem.swap(paux);
            }
            
            #pragma omp parallel for
            for(int i=0; i<static_cast<int>(subpart.Conditions().size()); ++i)
            {
                Condition::Pointer& pcond = *(subpart.ConditionsBegin() + i).base();            
                auto paux = mrModelPart.Conditions()(pcond->Id());
                pcond.swap(paux);
            }
            
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
                        graph[geom[k].Id()].insert(geom[l].Id());
            }
            
            for(int i=0; i<static_cast<int>(mrModelPart.Conditions().size()); ++i)
            {
                const auto& cond = mrModelPart.ConditionsBegin() + i;              
                const auto& geom = cond->GetGeometry();
                for(unsigned int k=0; k<geom.size(); ++k)
                    for(unsigned int l=0; l<geom.size(); ++l)
                        graph[geom[k].Id()].insert(geom[l].Id());
            }
            
             
            unsigned int nnz = 0;
            for(unsigned int i=0; i<n; ++i)
                nnz += graph[i].size();
            
//             compress_matrix<std::size_t> graph_csr(n,n);
//             for(unsigned int i=0; i<n; ++i)
//                 for(const auto& k : graph[i])
//                     graph_csr.push_back(i,k,0);
                
            //here must do the ordering!!
        
            
            
        }
        
	std::string ReorderAndOptimizeModelPartProcess::Info() const {
		return "ReorderAndOptimizeModelPartProcess";
	}

	void ReorderAndOptimizeModelPartProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void ReorderAndOptimizeModelPartProcess::PrintData(std::ostream& rOStream) const {

	}





}  // namespace Kratos.
