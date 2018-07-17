//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti and Pavel Ryzhakov
//


#if !defined(KRATOS_NIST_UTILITIES_INCLUDED )
#define  KRATOS_NIST_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "ULF_application.h"

namespace Kratos
{
class NistUtils
{
public:


    //**********************************************************************************************
    //**********************************************************************************************
    void GenerateModelPart(
        ModelPart& OriginModelPart ,
        ModelPart& DestinationModelPart,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition
    )
    {
        KRATOS_TRY;

        //assigning the nodes to the new model part
        DestinationModelPart.Nodes().clear();
        DestinationModelPart.Nodes() = OriginModelPart.Nodes();

        //generating the elements
        int id = 1;
        Properties::Pointer properties = OriginModelPart.GetMesh().pGetProperties(1);
        for(ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            Element::Pointer p_element = rReferenceElement.Create(id, geom ,properties);
            DestinationModelPart.Elements().push_back(p_element);
            id = id + 1;
        }
        std::cout << "Elements are generated" << std::endl;

        //generating the conditions
        id = 1;
        for(ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            double nfree_surf = 0;
            for(unsigned int k = 0; k<geom.size(); k++)
                nfree_surf += geom[k].FastGetSolutionStepValue(IS_FREE_SURFACE);

            if(nfree_surf > 1)
            {
                Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(id, geom,properties);
                DestinationModelPart.Conditions().push_back(p_condition);
                id = id + 1;
            }
        }
        std::cout << "Conditions are generated" << std::endl;

        KRATOS_CATCH("");
    }


    //**********************************************************************************************
    //**********************************************************************************************
    void ApplyInitialTemperature(
        ModelPart& ThisModelPart ,
        double wall_temperature
    )
    {
        KRATOS_TRY;
        for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                in!=ThisModelPart.NodesEnd(); in++)
        {
            if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 1 &&
                    (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0 )
            {
                in->FastGetSolutionStepValue(TEMPERATURE) = wall_temperature;
            }
        }

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
    double FindFluidLevel(
        ModelPart::NodesContainerType nodes
    )
    {
        KRATOS_TRY;

        double level = 0.00;
        for(ModelPart::NodesContainerType::iterator in = nodes.begin();
                in!=nodes.end(); in++)
        {
            if(  (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0 )
            {
                if( in->Y() > level) level = in->Y();
            }
        }
        return level;

        KRATOS_CATCH("");
    }

private:

};

}  // namespace Kratos.

#endif // KRATOS_NIST_UTILITIES_INCLUDED  defined 


