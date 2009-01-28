/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-01-21 17:06:57 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_MESH_SUITE_MODELER_H_INCLUDED )
#define  KRATOS_MESH_SUITE_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <iomanip>
#include <fstream>


// External includes 
#include "malla.h"
#include "voronoi.h" 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
//#include "geometries/triangle_2d.h"
#include "punto.h"
#include "nodo.h"

namespace Kratos
{

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
	class MeshSuiteModeler  : public malla
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MeshSuiteModeler
		KRATOS_CLASS_POINTER_DEFINITION(MeshSuiteModeler);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MeshSuiteModeler(){}

		/// Destructor.
		virtual ~MeshSuiteModeler(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		//**************************************************************************************
		//**************************************************************************************

		void GenerateMesh(ModelPart& ThisModelPart, Element const& rReferenceElement)
		{
			KRATOS_TRY

			mk_puntos(true);
			//delaunay(true,false);

			//erasing the elements
			ThisModelPart.Nodes().clear();
			ThisModelPart.Elements().clear();

			//generating in Kratos the new nodes

			ThisModelPart.Nodes().reserve(n.len);
			for(int i = 0 ; i < n.len ; i++)
			{
				Node<3>::Pointer p_node = Kratos::Node<3>::Pointer(new Kratos::Node<3>(0,n[i][0], n[i][1], n[i][2]));
				p_node->SetId(i+1);
				ThisModelPart.Nodes().push_back(p_node);
			}

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);

			
			ThisModelPart.Elements().reserve(e.len);
			for(int i_element = 0 ; i_element < e.len ; i_element++)
			{
				Element::NodesArrayType temp;

				for(int i_node = 0 ; i_node < e[i_element].nv() ; i_node++)
				{
					temp.push_back(ThisModelPart.Nodes()(e[i_element][i_node] + 1));
				}
				Element::Pointer p_element = rReferenceElement.Create(i_element + 1, temp, properties);
				ThisModelPart.Elements().push_back(p_element);
			}
			KRATOS_CATCH("")
		}

		//**************************************************************************************
		//**************************************************************************************
		void FixNodes(ModelPart& ThisModelPart)
		{

			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
			{
				n[i_node->Id()-1].f.set(n_permanente);
			}
		}

		//**************************************************************************************
		//**************************************************************************************
		int UpdateMesh(ModelPart& ThisModelPart)
		{
			// lee los nodos
			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
			{
				std::size_t id = i_node->Id() - 1;
				n[id][0] = i_node->X();
				n[id][1] = i_node->Y();
				n[id][2] = i_node->Z();
			}


			return true;
		}


		//**************************************************************************************
		//**************************************************************************************
		int SetMesh(ModelPart& ThisModelPart)
		{
			KRATOS_TRY
				ini();

			// lee de acuerdo a filext
			array1<nodo> nl; array1<elemento> el;

			nl.clean(); 
			el.clean();

			// lee los nodos
			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
				nl+=nodo(i_node->X(), i_node->Y() , i_node->Z());



			// lee los elementos
			for(ModelPart::ElementIterator i_element = ThisModelPart.ElementsBegin() ; i_element != ThisModelPart.ElementsEnd() ; ++i_element)
			{
				if(i_element->GetGeometry().PointsNumber() == 2)
				{
					elemento temp(e_segmento);
					temp[0] = i_element->GetGeometry()[0].Id() - 1;
					temp[1] = i_element->GetGeometry()[1].Id() - 1;
					// 	    KRATOS_WATCH(temp[0])
					// 	    KRATOS_WATCH(temp[1])
					el+=temp;
				}
				else if(i_element->GetGeometry().PointsNumber() == 3)
				{
					elemento temp(e_triangulo);
					temp[0] = i_element->GetGeometry()[0].Id() - 1;
					temp[1] = i_element->GetGeometry()[1].Id() - 1;
					temp[2] = i_element->GetGeometry()[2].Id() - 1;
					el+=temp;
				}
			}
			//   // nombre y extensión
			//   if (!inombre[0]) {add_error(No_Mesh); return false;}
			//   renombra(inombre);


			if (!nl) {add_error(No_Mesh); return false;}

			int i,j,in,nc=0;
			int &nlen=nl.len;
			n.resize(nlen);

			// bounding box
			pmin=pmax=nl[0];
			for (i=1;i<nlen;i++) {
				if (nl[i].f.es(n_borrado)) continue;
				pmin.set_min(nl[i]); pmax.set_max(nl[i]);
			}

			// hay z?
			bool hayz=pmax[2]-pmin[2]>ERRADM;
			if (!hayz) tipo.set(m_planaxy);
			else tipo.reset(m_planaxy);
			if(!o)
				o=new octree(n,pmin,pmax,hayz ? 3 : 2);
			bool puesto,hayrep=false;
			double d,dmax=0;
			pline map(nlen); map.len=nlen;
			// los nodos frontera van primero y los
			// de h (si hay) despues
			bool haynh=false;
			for (i=0; i<nlen; i++){
				if (nl[i].f.es(n_borrado)) continue;
				if (nl[i].f.es(n_h)) {haynh=true; continue;}
				in=map[i]=n+=nl[i]; // presupongo que es n.len-1
				o->add_no_rep(in,nc,puesto,epsilon);
				if (puesto) continue;
				if (nc<0||nc>=nlen)
				{add_error(Bad_Format); return false;}
				hayrep=true;
				d=n[nc].distancia(nl[i]);
				if (d>dmax) dmax=d;
				n.len--; map[i]=nc;// si no es n.len-1 usar remove
				n[nc].f.set(nl[i].f); // oreo flags
				if (n[nc].h>ERRADM&&nl[i].h>ERRADM)
					set_min(n[nc].h,nl[i].h); // h es el menor
				else
					set_max(n[nc].h,nl[i].h); // h > 0
				n[nc].v=(n[nc].v+nl[i].v)/2;
			}
			// los nodos de h
			if (haynh) for (i=0; i<nlen; i++){
				if (nl[i].f.es(n_borrado)) continue;
				if (nl[i].f.noes(n_h)) continue;
				in=map[i]=n+=nl[i]; // presupongo que es n.len-1
				//    o->add_no_rep(in,nc,puesto,epsilon);
				if (puesto) {nodosh++; if (n[in].h<ERRADM) n[in].h=MAXREAL; continue;}
				n.len--; map[i]=nc;// si no es n.len-1 usar remove
				// no oreo flags
				if (n[nc].h>ERRADM&&nl[i].h>ERRADM)
					set_min(n[nc].h,nl[i].h); // h es el menor
				else
					set_max(n[nc].h,nl[i].h); // h > 0
				if (n[nc].h<ERRADM) n[nc].h=MAXREAL;
				n[nc].v=(n[nc].v+nl[i].v)/2;
			}

			// modifica los elementos y hace eldenod
			if (el) {
				tipo.reset(m_nodos);
				int nv,dim,ie; flagtype fdim;
				e.resize(el);
				for (i=0;i<el.len;i++){
					elemento &ei=el[i]; nv=ei.nv(); dim=ei.dim();
					// verifica nodos h
					for (haynh=false,j=0;j<nv;j++) {
						if (n[ei[j]=map[ei[j]]].f.es(n_h)) haynh=true;
					}
					// verifica si el elm esta repetido
					cpline &en=n[ei[0]].e;
					for (j=0;j<en.len;j++) if (ei==e[en[j]]) break;
					if (j<en.len) { // repetido
						add_warning(Repeated_Elements);
						continue;
					}
					// verifica si el elm tiene nodos repetidos
					for (j=0;j<nv;j++)
						if (ei[j]==ei.npos(j))
							break;
					if (j<nv) continue; // no lo agrega
					// nodos h
					if (haynh){ // no agrega el elemento, solo saca el h de cada nodo
						for (j=0;j<nv;j++) {
							nodo &ni=n[ei[j]];
							if (ni.f.noes(n_h)) continue;
							if (ni.h!=MAXREAL) continue; // ya tiene h
							d=ni.distancia(n[ei.nant(j)]);
							if (d>ERRADM) ni.h=d;
							d=ni.distancia(n[ei.npos(j)]);
							if (d>ERRADM&&d<ni.h) ni.h=d;
						}
						continue;
					}
					// agrega el elemento
					ie=e+=ei;
					for (j=0;j<nv;j++) n[ei[j]].e+=ie;
					eltipos.set(ei.ftipo());
					if (dim==0) fdim.set(m_nodos);
					else if (dim==1) fdim.set(m_lin);
					else if (dim==2) fdim.set(m_sup);
					else if (dim==3) fdim.set(m_vol);
				}
				if (fdim==m_nodos) tipo.set(m_nodos);
				if (fdim==m_lin) tipo.set(m_lin);
				if (fdim==m_sup) tipo.set(m_sup);
				if (fdim==m_vol) tipo.set(m_vol);
			}
			else tipo.set(m_nodos);

			return true;
			KRATOS_CATCH("")
		}

		//*******************************************************************************************
		//*******************************************************************************************
		//void UpdateNodePosition(ModelPart& ThisModelPart)
		//{
		//    KRATOS_TRY
		//
		//        for(int i = 0 ; i < n.len ; i++)
		//        {
		//            Node<3>& r_node = ThisModelPart.Nodes()[i+1];
		//            n[i][0] = r_node.X();
		//            n[i][1] = r_node.Y();
		//            n[i][2] = r_node.Z();
		//        }
		//        KRATOS_CATCH("")
		//}

		//*******************************************************************************************
		//*******************************************************************************************
		//void SetNodalH(ModelPart& ThisModelPart)
		//{
		//    KRATOS_TRY
		//
		//		ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); 
		//        for(int i = 0 ; i < n.len ; i++)
		//        {
		//            Node<3>& r_node = *(in);
		//			double KratosH = r_node.FastGetSolutionStepValue(NODAL_H);

		//			if(r_node.FastGetSolutionStepValue(IS_BOUNDARY) == 0) //only on internal nodes
		//			{
		//				n[i].h = KratosH;
		//			}

		//			//KRATOS_WATCH(n[i].h);
		//			in++;
		//        }
		//        KRATOS_CATCH("")
		//}

		//*******************************************************************************************
		//*******************************************************************************************
		void ReGenerateMesh(
			ModelPart& ThisModelPart , 
			Element const& rReferenceElement, 
			Condition const& rReferenceBoundaryCondition,
			double my_alpha = 1.4)
		{
			KRATOS_TRY

			delaunay(true,false);
			alpha_shape(my_alpha);

			//mk_vecinos_y_frontera();
			mk_nn();
			//mk_frontera();

			ThisModelPart.Elements().clear();

			for(int i = 0 ; i < n.len ; i++)
			{
				Node<3>& r_node = ThisModelPart.Nodes()[i+1];
				r_node.X() = n[i][0];
				r_node.Y() = n[i][1];
				r_node.Z() = n[i][2];
			}

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);

			ModelPart::NodesContainerType& r_model_nodes = ThisModelPart.Nodes();

			//adding domain elements
			for(int i_element = 0 ; i_element < e.len ; i_element++)
			{
				Element::GeometryType temp;
				
				Element::GeometryType::Pointer p_geometry;

				int nv = e[i_element].nv();
				//temp.reserve(nv);
				for(int i_node = 0 ; i_node < nv ; i_node++)
				{
					int index = e[i_element][i_node] + 1;
					temp.push_back(r_model_nodes(index));
				}

				p_geometry = Element::GeometryType::Pointer(new Element::GeometryType(temp));

				Element::Pointer p_element = rReferenceElement.Create(i_element + 1, *p_geometry, properties);

				ThisModelPart.Elements().push_back(p_element);
			}


			//adding boundary elements
			//ThisModelPart.Conditions().clear();
			//for(int i_bound = 0; i_bound < frontera.len ; i_bound++)
			//{
			//	for(int i_element = 0 ; i_element < frontera[i_bound].e.len ; i_element++)
			//	{
			//		Condition::NodesArrayType temp;

			//		for(int i_node = 0 ; i_node < frontera[i_bound].e[i_element].nv() ; i_node++)
			//		{
			//			int index = frontera[i_bound].e[i_element][i_node] + 1;
			//			temp.push_back(r_model_nodes(index));
			//		}

			//		int pos = ThisModelPart.Conditions().size();
			//		Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(pos, temp, properties);

			//		ThisModelPart.Conditions().push_back(p_cond);

			//	}
			//}

			//resetting the boundary flags needed on the frontier
			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
			{
				i_node->GetSolutionStepValue(IS_BOUNDARY) = 0.00;
				i_node->GetSolutionStepValue(IS_FREE_SURFACE) = 0.00;
			}

			//setting the flags needed on the frontier
			for(int i_bound = 0; i_bound < frontera.len ; i_bound++)
			{
				for(int ii = 0 ; ii < frontera[i_bound].n.len ; ii++)
				{
					//attention: nodes identified with 1 are nodes in a larger domain
					//nodes identified with 2 identify nodes which are boundaries of a single isolated element

					if(frontera[i_bound].e.len > 3) //in this case the frontier is "big" as it has more than three edges
					{
						ThisModelPart.Nodes()[frontera[i_bound].n[ii]+1].GetSolutionStepValue(IS_BOUNDARY) = 1;
					}
					else
					{
						ThisModelPart.Nodes()[frontera[i_bound].n[ii]+1].GetSolutionStepValue(IS_BOUNDARY) = 2;
					}

					//attention: nodes which
				}
			}

			//free surface nodes are nodes of the boundary which are not nodes of the structure
			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
			{
				if(i_node->GetSolutionStepValue(IS_BOUNDARY) == 1.00 && i_node->GetSolutionStepValue(IS_STRUCTURE) != 1)
					i_node->GetSolutionStepValue(IS_FREE_SURFACE) = 1.00;
			}

			////fill in the nodal neighbour list
			for(int i = 0 ; i < nn.len ; i++)
			{
				int numb_of_neighb = nn[i].len;
				WeakPointerVector< Node<3> >& neighbours = (r_model_nodes[i+1].GetValues(NEIGHBOUR_NODES));
				//PointerVector< Node<3> >& neighbours = (r_model_nodes[i+1].GetValues(NEIGHBOUR_NODES));

				neighbours.clear();
				neighbours.reserve(numb_of_neighb);
				for(int j = 0; j<numb_of_neighb; j++)
				{
					int ii = nn[i][j]+1;
					//neighbours.push_back(boost::weak_ptr< Node<3> >( r_model_nodes(ii) ) );
					neighbours.push_back( r_model_nodes(ii)  );
				}
			}

			//perform an update to ensure the database is allocated correctly
			for(ModelPart::NodeIterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; ++i_node)
			{
				(i_node->SolutionStepData()).Update();
			}


			KRATOS_CATCH("")
		}

		//*******************************************************************************************
		//*******************************************************************************************
		//void RefineMesh(ModelPart& ThisModelPart)
		//{
		//	KRATOS_TRY
		//		//      voronoi v(this);
		//		//      v.parcial=false;

		//		//      if (!v.delaunay()) // triangulacion (del para que no meta ni saque ptos)
		//		//        return ; //false

		//		//      double alpha=1.2;
		//		//      v.alpha_shape(alpha);// alpha shape

		//		//misma cantidad de nodos porque hice delaunay sin refinar
		//		for(int i=0;i<n.len;i++)  n[i].h=fabs(n[i][0]/10.00 + n[i][1]/10.00) + 0.2; 

		//	//       bool calpha = false; // Don't do alpha_shape
		//	// //       bool calpha = true; // Do alpha_shape
		//	//      if (!v.refina_esferas(calpha)) // triangulacion
		//	//        return ; //false
		//	// //     if (slivers) v.rm_slivers();
		//	//      v.s2e();
		//	//      mk_nn();

		//	delaunay_refinado(true);
		//	ThisModelPart.Nodes().clear();
		//	ThisModelPart.Elements().clear();

		//	for(int i = 0 ; i < n.len ; i++)
		//	{
		//		Node<3>::Pointer p_node = Kratos::Node<3>::Pointer(new Kratos::Node<3>(0,n[i][0], n[i][1], n[i][2]));
		//		//Node<3>::Pointer p_node = n[i].GetKratosNode();
		//		p_node->SetId(i+1);
		//		ThisModelPart.Nodes().push_back(p_node);
		//	}

		//	//Node<3> zero(0,0,0,0);
		//	Properties::Pointer properties(new Properties);

		//	for(int i_element = 0 ; i_element < e.len ; i_element++)
		//	{
		//		Element::GeometryType temp;
		//		Element::GeometryType::Pointer p_geometry;

		//		for(int i_node = 0 ; i_node < e[i_element].nv() ; i_node++)
		//			temp.push_back(ThisModelPart.Nodes()[e[i_element][i_node] + 1]);


		//		//if(e[i_element].tipo() == e_triangulo)
		//		//	p_geometry =  Element::GeometryType::Pointer(new Triangle2D<Node<3> >(temp.Points()));
		//		//else
		//		p_geometry = Element::GeometryType::Pointer(new Element::GeometryType(temp));

		//		Element::Pointer p_element(new Element(i_element + 1, p_geometry, properties));
		//		ThisModelPart.Elements().push_back(p_element);
		//	}
		//	KRATOS_CATCH("")
		//}



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
		virtual std::string Info() const{return "";}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const{}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const{}


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
		MeshSuiteModeler& operator=(MeshSuiteModeler const& rOther);

		/// Copy constructor.
		MeshSuiteModeler(MeshSuiteModeler const& rOther);


		///@}    

	}; // Class MeshSuiteModeler 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MeshSuiteModeler& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MeshSuiteModeler& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_MESH_SUITE_MODELER_H_INCLUDED  defined 


