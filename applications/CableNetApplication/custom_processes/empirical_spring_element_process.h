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


/**
 * @class EmpiricalSpringElementProcess
 *
 * @brief This process creates a spring element w.r.t. to given displacement/load data points
 *
 * @author Klaus B Sautter
 */


#ifndef EMPIRICAL_SPRING_ELEMENT_PROCESS_H
#define EMPIRICAL_SPRING_ELEMENT_PROCESS_H

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
#include "geometries/line_3d_2.h"
#include "includes/model_part.h"
#include "cable_net_application_variables.h"

namespace Kratos
{


class EmpiricalSpringElementProcess
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
    KRATOS_CLASS_POINTER_DEFINITION(EmpiricalSpringElementProcess);

    /// Constructor.
    EmpiricalSpringElementProcess(ModelPart &rModelPart,
     Parameters InputParameters, const DoubleVector FittedPolynomial):mrModelPart(rModelPart),mParameters(InputParameters),mrFittedPoly(FittedPolynomial)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"           : "example_part",
            "computing_model_part_name" : "computing_domain",
            "node_ids"                  : [1,2],
            "element_id"                : 1,
            "property_id"               : 1,
            "displacement_data"         : [0.0,1.0,2.0,3.0],
            "force_data"                : [0.0,1.0,2.0,3.0],
            "polynomial_order"          : 3
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);

        KRATOS_ERROR_IF(mParameters["node_ids"].size()!=2) << "exactly two nodes for each spring needed !" << std::endl;
        KRATOS_ERROR_IF(mParameters["displacement_data"].size()!=mParameters["force_data"].size()) << "only two nodes for each spring allowed !" << std::endl;


        KRATOS_CATCH("")
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        this->CreateEmpiricalSpringElement();
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase

    void CreateEmpiricalSpringElement() const
    {
        KRATOS_TRY;
        Vector polynomial_order = ZeroVector(mrFittedPoly.size());
        for (SizeType i=0;i<mrFittedPoly.size();++i) polynomial_order[i] = mrFittedPoly[i];
        const int number_nodes = 2;

        // get new element id
        const std::size_t new_element_id = mParameters["element_id"].GetInt();

        // create geometric entitity
        std::vector<NodeType::Pointer> element_nodes (number_nodes);
        for (SizeType i=0; i<number_nodes; ++i)
        {
            element_nodes[i] = mrModelPart.pGetNode(mParameters["node_ids"][i].GetInt());
        }
        Line3D2 <NodeType> line_t ( PointerVector<NodeType>{element_nodes} );

        // get properties
        Properties::Pointer p_elem_prop = mrModelPart.pGetProperties(mParameters["property_id"].GetInt());


        p_elem_prop->SetValue(SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL, polynomial_order);


        const Element& rElem = KratosComponents<Element>::Get("EmpiricalSpringElement3D2N");
        Element::Pointer pElem = rElem.Create(new_element_id, line_t, p_elem_prop);
        mrModelPart.AddElement(pElem);


        KRATOS_CATCH("");
    }

  protected:


  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    const DoubleVector mrFittedPoly;

}; // Class

}; // namespace

#endif
