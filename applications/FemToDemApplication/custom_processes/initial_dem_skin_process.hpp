//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes

#if !defined(KRATOS_INITIAL_DEM_SKIN_PROCESS_HPP_INCLUDED)
#define KRATOS_INITIAL_DEM_SKIN_PROCESS_HPP_INCLUDED

#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{

class InitialDemSkinProcess : public Process
{
  public:
    KRATOS_CLASS_POINTER_DEFINITION(InitialDemSkinProcess);

    // Constructor
    InitialDemSkinProcess(ModelPart &rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /* It will create a submodel part containing nodes to include DEM on them,
    with the DEM radius assigned to each node to be created the DEM afterwards
    The modelpart must include the skinModelpart and damage extrapolated to nodes*/
    void Execute()
    {
        const std::string &name_dem_model_part = "InitialDemSkin";

        if (mrModelPart.HasSubModelPart(name_dem_model_part))
        {
            mrModelPart.RemoveSubModelPart(name_dem_model_part);
        }

        FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(mrModelPart, 4, 4);
        neighbour_finder.Execute();

        mrModelPart.CreateSubModelPart(name_dem_model_part);
        ModelPart* p_auxiliar_model_part = mrModelPart.pGetSubModelPart(name_dem_model_part);

        ModelPart* p_skin_model_part = mrModelPart.pGetSubModelPart("SkinDEMModelPart");

        for (ModelPart::NodeIterator it = (*p_skin_model_part).NodesBegin(); it != (*p_skin_model_part).NodesEnd(); ++it) {
            p_auxiliar_model_part->AddNode(*(it.base()));
        } // InitialDemSkin SubModelPart Filled with nodes

        // Let's assign the DEM radius to those nodes...
        for (ModelPart::NodeIterator it = (*p_auxiliar_model_part).NodesBegin(); it != (*p_auxiliar_model_part).NodesEnd(); ++it) {
            WeakPointerVector<Node<3>> &rneigh = (*it).GetValue(NEIGHBOUR_NODES);
            KRATOS_ERROR_IF(rneigh.size() == 0) << "Nodal neighbours not computed..." << std::endl;
            std::vector<double> radius_is_dems, radius_not_dem;
            double distance, radius_dem, min_radius, min_radius_is_dem, min_radius_no_dem;

            for (int i = 0; i < rneigh.size(); i++)
            {

                distance = this->CalculateDistanceBetweenNodes((*it), rneigh[i]);

                if (rneigh[i].GetValue(DEM_RADIUS) != 0.0)
                {
                    radius_dem = distance - rneigh[i].GetValue(DEM_RADIUS);
                    radius_is_dems.push_back(radius_dem);
                }
                else
                {
                    radius_dem = 0.5 * distance;
                    radius_not_dem.push_back(radius_dem);
                }
            }
            if (radius_is_dems.size() == 0) {
                min_radius = this->GetMinimumValue(radius_not_dem);
                (*it).SetValue(DEM_RADIUS, min_radius);
            } else {
                min_radius_is_dem = this->GetMinimumValue(radius_is_dems);

                if (radius_not_dem.size() != 0)
                    min_radius_no_dem = this->GetMinimumValue(radius_not_dem);
                else
                    min_radius_no_dem = 1000.0;

                if (min_radius_is_dem > 1.5 * min_radius_no_dem)
                    min_radius = min_radius_no_dem;
                else
                    min_radius = min_radius_is_dem;
                (*it).SetValue(DEM_RADIUS, min_radius);
            }
        }
    }

    double CalculateDistanceBetweenNodes(
        const Node<3> &Node1,
        const Node<3> &Node2)
    {
        const double X1 = Node1.X();
        const double X2 = Node2.X();
        const double Y1 = Node1.Y();
        const double Y2 = Node2.Y();
        const double Z1 = Node1.Z();
        const double Z2 = Node2.Z();
        return std::sqrt(std::pow((X1 - X2), 2) + std::pow((Y1 - Y2), 2) + std::pow((Z1 - Z2), 2));
    }

    double GetMinimumValue(std::vector<double> array_values)
    {
        int size = array_values.size();

        double aux = array_values[0];
        for (int i = 1; i < size; i++)
        {
            if (array_values[i] < aux)
                aux = array_values[i];
        }
        return aux;
    }

  protected:
    // Member Variables
    ModelPart &mrModelPart;
};
} // namespace Kratos
#endif
