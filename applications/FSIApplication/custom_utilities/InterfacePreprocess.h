/*
 * File:   InterfacePreprocess.h
 * Author: jcotela
 * Co-Author: VMataix
 * Created on 19 January 2010, 11:44
 * Last update on 29 April 2016, 16:54
 */

#if !defined(KRATOS_INTERFACE_PREPROCESS_MAPPER_H_INCLUDED )
#define  KRATOS_INTERFACE_PREPROCESS_MAPPER_H_INCLUDED

#include <iostream>

#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
//#include "geometries/quadrilateral_3d_4.h" // TODO: Add in the future -> Can be replaced with the triangles
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "includes/define.h"
//~ #include "includes/deprecated_variables.h"

namespace Kratos
{

/** \brief InterfacePreprocess Methods
 * Creates Model Parts containing the interface, to be used by AdvancedNMPointsMapper
 */
class InterfacePreprocess
{

public:

    /**
     * It prints the nodes and conditions in the interface, gives an error otherwise there are not
     * @param NodesCounter: Number of nodes in the interface
     * @return CondCounter: Number of conditions in the interface
     */

    void PrintNodesAndConditions(
            const int NodesCounter,
            const int CondCounter
            )
    {
        std::cout << "    " << NodesCounter << " nodes ";
        std::cout << "and " << CondCounter <<  " conditions found." << std::endl;

        // Check that we actually found something
        if( NodesCounter == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument,"No interface nodes found. Please check that nodes on both sides of the interface have been assigned INTERFACE=1.0.","");
        }
        if( CondCounter == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument,"No interface conditions found. Please check that nodes on both sides of the interface have been assigned INTERFACE=1.0 and that the contact surfaces have been assigned conditions.","");
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only linear conditions (just for the 2D cases)
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLineInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart
            )
    {
        KRATOS_TRY;

        const unsigned int dimension = rOriginPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        if (dimension == 2) // By default, but someone can be interested in project values to a BEAM for example
        {
            GenerateLine2DConditions(aux, InterfacePart.Conditions());
        }
        else
        {
            GenerateLine3DConditions(aux, InterfacePart.Conditions());
        }

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular conditions, regardless of what was used while meshing
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangleInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            unsigned int number_nodes = (*cond_it).GetGeometry().PointsNumber();
            if (number_nodes >= 3) 
            {
                bool is_interface = true;  
                for (unsigned int node = 0; node < number_nodes; node++)
                {
                    if ((*cond_it).GetGeometry()[node].Is(INTERFACE) != true)
                    {
                        is_interface = false;
                    }
                }
                if (is_interface == true)
                {
                    aux.push_back( *(cond_it.base()) );
                    CondCounter ++;
                }
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateTriangularConditions(aux, InterfacePart.Conditions());

        KRATOS_CATCH("");

    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 2D lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2DConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get("LineCondition2D2N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 3D lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3DConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get("LineCondition3D2N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of linear triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangularConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get("SurfaceCondition3D3N"); /* Face3D3N */

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType TriId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {

            if (it->GetGeometry().PointsNumber() == 3)
            {
                rTriConds.push_back( rCondition.Create(TriId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 4)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(2));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(0));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 6)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 8)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri5(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(8), it->GetGeometry()(3), it->GetGeometry()(4));
                Triangle3D3< Node<3> > Tri5(it->GetGeometry()(8), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri7(it->GetGeometry()(5), it->GetGeometry()(7), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri8(it->GetGeometry()(0), it->GetGeometry()(8), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri7, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri8, it->pGetProperties()));
            }
        }
    }
};
}

#endif
