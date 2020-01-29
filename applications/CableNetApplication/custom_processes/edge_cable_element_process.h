//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//
//



#ifndef EDGE_CABLE_ELEMENT_PROCESS_H
#define EDGE_CABLE_ELEMENT_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "geometries/line_3d_3.h"
#include "custom_geometries/line_3d_n.h"
#include "includes/model_part.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) EdgeCableElementProcess
    : public Process
{
  public:

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(EdgeCableElementProcess);

    /// Constructor.
    EdgeCableElementProcess(ModelPart &rModelPart,
     Parameters InputParameters):mrModelPart(rModelPart),mParameters(InputParameters)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "edge_sub_model_part_name"  : "Structure.example_part",
            "element_type"              : "cable",
            "node_id_order"             : [1,2,3],
            "element_id"                : 1,
            "property_id"               : 1
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);
        KRATOS_CATCH("")
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        this->CreateEdgeCableElement();
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }



    void CreateEdgeCableElement() const
    {
        KRATOS_TRY;
        // !!
        // probably better to create the line equation and calculate the ordering here
        // !!

        const SizeType number_nodes     = mParameters["node_id_order"].size();
        KRATOS_ERROR_IF_NOT(mrModelPart.Nodes().size()==number_nodes)
         << "numbers of nodes in submodel part not consistent with numbers of nodes in process properties"
         << std::endl;


        // get new element id
        const std::size_t new_element_id = mParameters["element_id"].GetInt();

        // create geometric entitity
        std::vector<NodeType::Pointer> element_nodes (number_nodes);
        for (SizeType i=0; i<number_nodes; ++i)
        {
            element_nodes[i] = mrModelPart.pGetNode(mParameters["node_id_order"][i].GetInt());
        }
        Line3DN <NodeType> line_t ( PointerVector<NodeType>{element_nodes} );

        // get properties
        Properties::Pointer p_elem_prop = mrModelPart.pGetProperties(mParameters["property_id"].GetInt());

        // create element
        if (mParameters["element_type"].GetString() == "cable")
        {
            const Element& rElem = KratosComponents<Element>::Get("SlidingCableElement3D3N");
            Element::Pointer pElem = rElem.Create(new_element_id, line_t, p_elem_prop);
            mrModelPart.AddElement(pElem);
        }
        else if (mParameters["element_type"].GetString() == "ring")
        {
            const Element& rElem = KratosComponents<Element>::Get("RingElement3D4N");
            Element::Pointer pElem = rElem.Create(new_element_id, line_t, p_elem_prop);
            mrModelPart.AddElement(pElem);
        }
        else KRATOS_ERROR << "element type :" << mParameters["element_type"].GetString() << " not available for sliding process" << std::endl;




        KRATOS_CATCH("");
    }

  protected:


  private:

    ModelPart& mrModelPart;
    Parameters mParameters;

}; // Class

}; // namespace

#endif
