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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_COPY_MODELER_H_INCLUDED )
#define  KRATOS_COPY_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"

namespace Kratos {

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
    class DuplicateMeshModeler : public Modeler {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of DuplicateMeshModeler
        KRATOS_CLASS_POINTER_DEFINITION(DuplicateMeshModeler);

        typedef Modeler BaseType;

        typedef Point < 3 > PointType;

        typedef Node < 3 > NodeType;

        typedef Geometry<NodeType> GeometryType;

        typedef PointerVector<NodeType> NodesVectorType;

        typedef std::size_t SizeType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// constructor.

        DuplicateMeshModeler(ModelPart& rSourceModelPart) :
        mrModelPart(rSourceModelPart) {
        }

        /// Destructor.

        virtual ~DuplicateMeshModeler() {
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement, Condition const& rReferenceCondition) {
			KRATOS_TRY;

            Timer::Start("Generating Mesh");

            Timer::Start("Generating Nodes");

            rThisModelPart.Nodes().clear();

            // Generate Copy of Nodes
            for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++) {
                rThisModelPart.CreateNewNode(i_node->Id(), *i_node);
            }


            Timer::Stop("Generating Nodes");

            rThisModelPart.PropertiesArray() = mrModelPart.PropertiesArray();

            Timer::Start("Generating Elements");


	//generating the elements
        const std::size_t element_size = rReferenceElement.GetGeometry().size();
        Element::NodesArrayType element_nodes_array(element_size);

	for(ModelPart::ElementsContainerType::iterator i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++)
	{
            Element::GeometryType& geometry = i_element->GetGeometry();
            if(geometry.size() != element_size)
                KRATOS_ERROR(std::invalid_argument, "The given element is not compatible with the reference element", "");

            for(std::size_t i = 0 ; i < element_size ; i++)
                element_nodes_array(i) = rThisModelPart.pGetNode(geometry[i].Id());
            
            Element::Pointer p_element = rReferenceElement.Create(i_element->Id(), element_nodes_array, i_element->pGetProperties());
	    rThisModelPart.Elements().push_back(p_element);
	}

            Timer::Stop("Generating Elements");

             Timer::Start("Generating Conditions");


	//generating the conditions
        const std::size_t condition_size = rReferenceCondition.GetGeometry().size();
        Condition::NodesArrayType condition_nodes_array(condition_size);

	for(ModelPart::ConditionsContainerType::iterator i_condition = mrModelPart.ConditionsBegin(); i_condition != mrModelPart.ConditionsEnd(); i_condition++)
	{
            Condition::GeometryType& geometry = i_condition->GetGeometry();
            if(geometry.size() != condition_size)
                KRATOS_ERROR(std::invalid_argument, "The given condition is not compatible with the reference condition", "");

            for(std::size_t i = 0 ; i < condition_size ; i++)
                condition_nodes_array(i) = rThisModelPart.pGetNode(geometry[i].Id());

            Condition::Pointer p_condition = rReferenceCondition.Create(i_condition->Id(), condition_nodes_array, i_condition->pGetProperties());
	    rThisModelPart.Conditions().push_back(p_condition);
	}

            Timer::Stop("Generating Conditions");

           Timer::Stop("Generating Mesh");

//			//generating the conditions
//			id = 1;
//			for(ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
//			{
//				Geometry< Node<3> >& geom = iii->GetGeometry();
//				double nfree_surf = 0;
//				for(unsigned int k = 0; k<geom.size(); k++)
//					nfree_surf += geom[k].FastGetSolutionStepValue(IS_FREE_SURFACE);
//
//				if(nfree_surf > 1)
//				{
//					Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(id, geom,properties);
//						DestinationModelPart.Conditions().push_back(p_condition);
//						id = id + 1;
//				}
//			}
//			std::cout << "Conditions are generated" << std::endl;

			KRATOS_CATCH("");
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

        virtual std::string Info() const {
            return "DuplicateMeshModeler";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
            rOStream << Info();
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
        }


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

        ModelPart& mrModelPart;

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
        DuplicateMeshModeler & operator=(DuplicateMeshModeler const& rOther);

        /// Copy constructor.
        DuplicateMeshModeler(DuplicateMeshModeler const& rOther);


        ///@}

    }; // Class DuplicateMeshModeler

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream & operator >>(std::istream& rIStream,
            DuplicateMeshModeler& rThis);

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const DuplicateMeshModeler& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


} // namespace Kratos.

#endif // KRATOS_COPY_MODELER_H_INCLUDED  defined


