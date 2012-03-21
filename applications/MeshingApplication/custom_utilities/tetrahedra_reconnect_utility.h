<<<<<<< .mine
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


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

#include "u_qualityMetrics.h"
#include "Math3D.h"
#include "u_types.h"
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

    /// Default constructor.
    TetrahedraReconnectUtility() {}

    /// Destructor.
    virtual ~TetrahedraReconnectUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    /// this function evaluates the quality of all the tetrahedras within a given model_part
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
		int numVertexes = m->vertexes->Count();
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

	}

    void innerConvertToKratos(ModelPart& mrModelPart , TVolumeMesh *m)
	{
		std::cout << "Generating nodes for Kratos " << "\n";
		Element rReferenceElement =  mrModelPart.GetElement(1);		
		Element::Pointer p_ele=  mrModelPart.pGetElement(1);		
		if (p_ele ==NULL )
		{
			std::cout << "Invalid element access " << "\n";
			return ;
		}

		//generate new nodes
		    mrModelPart.Elements().clear();
			mrModelPart.Nodes().clear();

            Node < 3 > ::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin())->GetDofs();
            for(int i=0; i< m->vertexes->Count(); i++)
            {
                // Get Vertex
				TVertex* v = m->vertexes->elementAt(i);
				v->id = i+1;
                double x = v->fPos.x;
                double y = v->fPos.y;
                double z = v->fPos.z;

                
				Node<3>::Pointer p_new_node = Node<3>::Pointer(new Node<3>(v->id, x, y, z));
				// Store Ref to Kratos
				v->userData =  (object*)( &p_new_node);

                // Giving model part's variables list to the node
                p_new_node->SetSolutionStepVariablesList(&(mrModelPart.GetNodalSolutionStepVariablesList()));
                p_new_node->SetBufferSize(mrModelPart.NodesBegin()->GetBufferSize());
                for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                {
                    Node < 3 > ::DofType& rDof = *iii;
                    Node < 3 > ::DofType::Pointer p_new_dof = p_new_node->pAddDof(rDof);
                    (p_new_dof)->FreeDof(); //the variables are left as free for the internal node
                }  
				
				mrModelPart.Nodes().push_back(p_new_node);
            }
		std::cout << "Generating Elements for Kratos " << "\n";
            //generate new Elements           
            ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
            Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(0);
			
            for(int i=0; i< m->elements->Count() ; i++)
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

                 Element::Pointer p_element = rReferenceElement.Create(i+1, geom, properties);
	            (mrModelPart.Elements()).push_back(p_element);
            }
          std::cout << "Generation OK " << "\n";

		  TVolumeMesh *testM = new TVolumeMesh();
		  innerConvertFromKratos( mrModelPart , testM);

	}

	/**
         * This function performs the meshing optimization by Cluster reconnection. 
         */
	
	void EvaluateQuality(ModelPart& r_model_part, bool flag , bool saveToFile, bool removeFreeVertexes )
    {
	    std::cout << "Creating mesh" << "\n";
        TVolumeMesh *m = new TVolumeMesh();
		// Convert to inner format
		innerConvertFromKratos(r_model_part , m );
		
		// Save the mesh as generated from Kratos
		if (saveToFile)
		   {
		       TMeshLoader* ml2 = new TVMWLoader();
	           ml2->save("D:/out_MeshFromKratos.vwm" , m);
			   delete ml2;
			}

		std::cout <<"...Optimizing by Face" <<"\n"; 
		TetQuality *qt = new TetQuality(m);
   
		qt->refresh();   qt->print();
		std::cout <<"...Parallel Optimizing by Node" <<"\n"; 
		startProcess("Parallel evaluation");
		    //ParallelEvaluateClusterByNode(m,vrelaxQuality);   
            evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
			m->updateIndexes(0);
		endProcess("Parallel evaluation");
		qt->refresh();   qt->print();
		m->validate(true);

		std::cout <<"...Parallel Optimizing by Node" <<"\n"; 
		startProcess("Parallel evaluation");
		    //ParallelEvaluateClusterByNode(m,vrelaxQuality);   
			evaluateClusterByNode( (TVolumeMesh*)(m),5000000,vrelaxQuality);
			m->updateIndexes(0);
		endProcess("Parallel evaluation");
		qt->refresh();   qt->print();
		m->validate(true);

		if (removeFreeVertexes ) 
			m->removeFreeVertexes();

		showProcessTime();

		if (saveToFile)
		{
			TMeshLoader* ml2 = new TVMWLoader();
			ml2->save("D:/out_MeshC.vwm" , m);
			delete ml2;
		}
			
		std::cout <<"...Output to Kratos Format" <<"\n";         
		// Get back in Kratos
		innerConvertToKratos(r_model_part , m);
		
		delete m;
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
    TetrahedraReconnectUtility& operator=(TetrahedraReconnectUtility const& rOther) {}

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


