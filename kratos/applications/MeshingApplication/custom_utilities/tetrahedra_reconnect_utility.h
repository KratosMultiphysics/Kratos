//
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <omp.h>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/node_erase_process.h"

#include "u_qualityMetrics.h"
#include "Math3D.h"
#include "u_Types.h"
#include "u_TetraFunctions.h"
#include "u_ShowMetrics.h"
#include "u_ParallelFunctions.h"
#include "u_MeshLoaders.h"
#include "u_elementCluster.h"
#include "u_ProcessTime.h"



namespace Kratos
{
	///@addtogroup ApplicationNameApplication
	///@{

	///@name Kratos Globals
	///@{

	///@}
	///@name Type Definitions
	///@{

	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class TetrahedraReconnectUtility
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TetrahedraReconnectUtility
		KRATOS_CLASS_POINTER_DEFINITION(TetrahedraReconnectUtility);

		///@}
		///@name Life Cycle
		///@{
		// inner Mesh
		TVolumeMesh *m ;
		
		ModelPart refMP ;

		/// Default constructor.
		TetrahedraReconnectUtility(ModelPart& r_model_part) 
		{
			std::cout << "Creating mesh" << "\n";
			m = new TVolumeMesh();
			// Convert to inner format
			innerConvertFromKratos(r_model_part , m );
			refMP = r_model_part;
		
		}

		/// Destructor.
		virtual ~TetrahedraReconnectUtility() {}


		///@}
		///@name Operators
		///@{


		///@}
		///@{
		
		///@brief function innerConvertFromKratos
		/// This function converts from Kratos format, to the inner structure
		///@param r_model_part the input mesh
		///@param m the output mesh
		void innerConvertFromKratos(ModelPart& r_model_part, TVolumeMesh *m)
		{
			std::cout << "Reading nodes"<< "\n";
			//loop on nodes
			for (ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
			{
				float4 fPos ;
				fPos.x = it->X();
				fPos.y = it->Y();
				fPos.z = it->Z();
				TVertex *v = new TVertex(fPos);
				v->id = it->Id();

				m->vertexes->Add(v);
			}
			std::cout << "Reading elements"<< "\n";
			
			for (ModelPart::ElementsContainerType::iterator el_it=r_model_part.ElementsBegin(); el_it!=r_model_part.ElementsEnd(); el_it++)
			{
				Geometry< Node<3> >& geom = el_it->GetGeometry();
				if (geom.size() != 4)
				{
					std::cout << "Invalid element size" <<  el_it->Id();
					continue;
				}
				TVertex *v0 = m->findVertexById(geom[0].Id());
				if (  v0 == NULL )
				{
					std::cout << "Invalid element reference" <<  el_it->Id();
					continue;
				}
				TVertex *v1 = m->findVertexById(geom[1].Id());
				if (  v1 == NULL )
				{
					std::cout << "Invalid element reference" <<  el_it->Id();
					continue;
				}
				TVertex *v2 = m->findVertexById(geom[2].Id());
				if (  v2 == NULL )
				{
					std::cout << "Invalid element reference" <<  el_it->Id();
					continue;
				}
				TVertex *v3 = m->findVertexById(geom[3].Id());
				if (  v3 == NULL )
				{
					std::cout << "Invalid element reference" <<  el_it->Id();
					continue;
				}

				TTetra *t = new TTetra(NULL, v0,v1,v2,v3);
				m->elements->Add(t);
			}

			std::cout << " Number of vertexes read :"<< m->vertexes->Count() <<"\n";
			std::cout << " Number of elements read :"<< m->elements->Count() << "\n";
			m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
			std::cout << " Number of faces read :"<< m->fFaces->Count() << "\n";

		}

		///@brief function innerConvertToKratos
		/// This function converts back to Kratos format
		///@param m the input mesh
		///@param mrModelPart the output mesh
		///@param removeFreeVertexes Free vertexes are not added to the new structure
		void innerConvertToKratos(ModelPart& mrModelPart , TVolumeMesh *m, bool removeFreeVertexes)
		{
			std::cout << "-------------Generating for Kratos----------------" << "\n";

			Element::Pointer pReferenceElement = *(mrModelPart.Elements().begin()).base();

			// Mark elements to delete
			// 0 Not remove
			// 1 Remove
			for (int i=0; i< m->vertexes->Count(); i++)
			{
				m->vertexes->structure[i]->id = i+1;
				if (m->vertexes->structure[i]->elementsList->Count() == 0)
					m->vertexes->elementAt(i)->flag = 1;
				else
					m->vertexes->elementAt(i)->flag = 0;
			}
			std::cout << " Nodes : Mark elements to remove : " << "\n";
			// Create new Nodes
			for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.Nodes().begin() ; i_node != mrModelPart.Nodes().end() ; i_node++)
			{

				TVertex* v = m->findVertexById(i_node->Id());
				if (v == NULL)
				{
					std::cout << "Error at id "<<i_node->Id() << "\n";
					continue;
				}
				i_node->SetValue(ERASE_FLAG ,v->flag == 1);
			}
			std::cout << " Generating Elements for Kratos " << "\n";
			//generate new nodes
			mrModelPart.Elements().clear();
			//add preserved elements to the kratos
			Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(1);
			(mrModelPart.Elements()).reserve(m->elements->Count());

			for (int i=0; i< m->elements->Count() ; i++)
			{

				TTetra* t = (TTetra*)(m->elements->elementAt(i));
				Node<3>::Pointer v0 = mrModelPart.pGetNode(t->vertexes[0]->id);
				Node<3>::Pointer v1 = mrModelPart.pGetNode(t->vertexes[1]->id);
				Node<3>::Pointer v2 = mrModelPart.pGetNode(t->vertexes[2]->id);
				Node<3>::Pointer v3 = mrModelPart.pGetNode(t->vertexes[3]->id);
				if ((v0 == NULL) || (v1==NULL) || (v2 == NULL) || (v3 == NULL))
				{
					std::cout << "Invalid vertex access " << t->id <<"\n";
					continue ;
				}
				Tetrahedra3D4<Node<3> > geom(
					v0,v1,v2,v3
					);

				Element::Pointer p_element = pReferenceElement->Create(i+1, geom, properties);
				(mrModelPart.Elements()).push_back(p_element);
				//KRATOS_WATCH(*p_element);
			}
			std::cout << "Generation OK " << "\n";
			mrModelPart.Elements().Sort();

			if (removeFreeVertexes )
			{
				std::cout << "Removing free vertexes " << "\n";
				(NodeEraseProcess(mrModelPart)).Execute();
			}
			std::cout << "-----------------Generation Finished OK-------------------" << "\n";

			//  TVolumeMesh *testM = new TVolumeMesh();
			//  innerConvertFromKratos( mrModelPart , testM);

		}

		void EvaluateQuality()
		{
		    TetQuality *qt = new TetQuality(m) ;
			qt->refresh();
			qt->print();
			delete qt;
		}

		/**
		* This function performs the meshing optimization by Cluster reconnection.
		*/
		void TestRemovingElements()
		{
			for (int i=0 ;i<100 ; i++)
			{			   
				m->elementsToRemove->Add( m->elements->elementAt(i));
			}
			m->updateRefs();
		}
		
		///@brief function updateNodesPositions
		/// Update only nodes positions without regenerating the structure
		///@param r_model_part the input mesh
		void updateNodesPositions(ModelPart& r_model_part)
		{
			std::cout << "Updating nodes"<< "\n";
			//loop on nodes
			for (ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
			{
				float4 fPos ;
				fPos.x = it->X();
				fPos.y = it->Y();
				fPos.z = it->Z();
				fPos.w = 1.0f;
				
				int id = it->Id();
				TVertex *v = m->findVertexById(id);
				if (v == NULL)
					 v->fPos = fPos;
			}
		}
		///@brief function OptimizeQuality
		///@param r_model_part the input mesh
		///@param iterations amount of iterations to optimize
		///@param processByNode boolean that indicates if processing is done by Node
		///@param processByFace boolean that indicates if processing is done by Face
		///@param processByEdge boolean that indicates if processing is done by Edge
		///@param saveToFile boolean to use as Debug Mode and see intermediate generated meshes
		///@param removeFreeVertexes boolean Removes vertexes that do not have elements referencing them
		///@param evaluateInParallel boolean Activate/Deactivate parallel processing mode
		///@param reinsertNodes boolean Activate/Deactivate to try re inserting removed nodes
		

		void OptimizeQuality(ModelPart& r_model_part, int iterations ,
			bool processByNode, bool processByFace, bool processByEdge,  
			bool saveToFile, bool removeFreeVertexes ,
			bool evaluateInParallel , bool reinsertNodes , bool debugMode)
		{
		    m->vertexes->Sort(sortByID);
			// Save the mesh as generated from Kratos
			if (saveToFile)
			{
				TMeshLoader* ml2 = new TVMWLoader();
				std::string s("out_MeshFromKratos.vwm");				
				ml2->save( s, m);
				delete ml2;
			}
			if (debugMode)
				EvaluateQuality();
			
			std::cout <<"...Start Optimization..." <<"\n";	
			if (evaluateInParallel)
			{
				std::cout <<"Number of active threads"<< omp_get_num_threads() <<"\n";
			}
			if (debugMode)
			{
				std::cout <<"Debug mode is Active" <<"\n";
				startTimers();
			}
			else
			{
				stopTimers();
			}

			for (int iter = 0 ; iter< iterations ; iter ++)
			{	
				//ParallelEvaluateClusterByNode(m,vrelaxQuality);
				if (processByNode)
				{
				    
					if (evaluateInParallel )
					{
					    std::cout <<"...Parallel optimizing by Node. Iteration : "<< iter <<"\n";
						ParallelEvaluateClusterByNode((TVolumeMesh*)(m),diedralAngle);   
					}
					else
					{
					   std::cout <<"...Optimizing by Node. Iteration : "<< iter <<"\n";
					   evaluateClusterByNode( (TVolumeMesh*)(m),5000000,diedralAngle);
					}
					if (debugMode)
						m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
				}

				if (processByFace)
				{
					std::cout <<"...Optimizing by Face. Iteration : "<< iter <<"\n";
					evaluateClusterByFace(m,500000,vrelaxQuality);
					if (debugMode)
						m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
				}

				if (processByEdge)
				{
					
					if (evaluateInParallel )
					{
						std::cout <<"...Parallel optimizing by Edge. Iteration : "<< iter <<"\n";
						ParallelEvaluateClusterByEdge((TVolumeMesh*)(m),diedralAngle);  
						std::cout <<"...End. Iteration : "<< iter <<"\n";
					}
					else
					{
					    std::cout <<"...Optimizing by Edge. Iteration : "<< iter <<"\n";
						evaluateClusterByEdge( (TVolumeMesh*)(m),50000,diedralAngle);
					}
					if (debugMode)
						m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
				}
				
				if (debugMode)
				{
					std :: cout<< "Number of faces:" << m->fFaces->Count() << "\n";
					EvaluateQuality();
					m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
					m->validate(true);				
					std :: cout<< "Number of faces:" << m->fFaces->Count() << "\n";
					showProcessTime();
				}
				// Save the mesh as generated from Kratos
				if (saveToFile)
				{
					TMeshLoader* ml2 = new TVMWLoader();
					std::string s("");					
					s = "out_MeshFromKratos" + intToString(iter)+".vwm";					
					
					// BUG Linux/Windows
					ml2->save(s , m);
					delete ml2;
				}
			}
						
			if (saveToFile)
			{
				TMeshLoader* ml2 = new TVMWLoader();
				std::string s("out_MeshC.vwm");
				ml2->save(s, m);
				delete ml2;
			}

			if (reinsertNodes)
			{
				std::cout <<"...Trying to reinsert nodes..." <<"\n";
				tryToReinsertNodes();
			}
			
		}
        ///@brief function tryToReinsertNodes
		/// Reinsert removed nodes into the structure
		void tryToReinsertNodes()
		{
			int vToR =m->vertexesToRemove->Count();
			int ri = vertexTetraReInsertion(m , m->vertexesToRemove);
			m->updateIndexes(GENERATE_SURFACE | KEEP_ORIG_IDS);
			std :: cout<< " Reinsert vertexes " << ri << " of " <<  vToR <<"\n";
			EvaluateQuality();
			m->validate(true);
			std :: cout<< "........................................"<<"\n";
		}
		///@brief function FinalizeOptimization
		/// Destroy the structure
		void FinalizeOptimization(bool removeFreeVertexes )
		{
			if (m == NULL ) return ;
			std::cout <<"...Output to Kratos Format" <<"\n";
			// Get back in Kratos
			innerConvertToKratos(refMP , m , removeFreeVertexes);
			delete m;
			m = NULL;
		}


		///@}
		///@name Access
		///@{


		///@}
		///@name Inquiry
		///@{


		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			std::stringstream buffer;
			buffer << "TetrahedraReconnectUtility" ;
			return buffer.str();
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "TetrahedraReconnectUtility";
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const {}


		///@}
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables
		///@{


		///@}
		///@name Protected member Variables
		///@{


		///@}
		///@name Protected Operators
		///@{


		///@}
		///@name Protected Operations
		///@{


		///@}
		///@name Protected  Access
		///@{


		///@}
		///@name Protected Inquiry
		///@{


		///@}
		///@name Protected LifeCycle
		///@{


		///@}

	private:
		///@name Static Member Variables
		///@{


		///@}
		///@name Member Variables
		///@{


		///@}
		///@name Private Operators
		///@{


		///@}
		///@name Private Operations
		///@{


		///@}
		///@name Private  Access
		///@{


		///@}
		///@name Private Inquiry
		///@{


		///@}
		///@name Un accessible methods
		///@{

		/// Assignment operator.
		TetrahedraReconnectUtility& operator=(TetrahedraReconnectUtility const& rOther) 
		{
		   return *this;
		}

		/// Copy constructor.
		TetrahedraReconnectUtility(TetrahedraReconnectUtility const& rOther) {}


		///@}

	}; // Class TetrahedraReconnectUtility

	///@}

	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,
		TetrahedraReconnectUtility& rThis) {}

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream,
		const TetrahedraReconnectUtility& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}

	///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED  defined 


