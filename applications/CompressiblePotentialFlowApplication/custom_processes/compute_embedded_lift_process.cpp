//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez,
//


#include "compute_embedded_lift_process.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
// Constructor for ComputeEmbeddedLiftProcess Process
ComputeEmbeddedLiftProcess::ComputeEmbeddedLiftProcess(ModelPart& rModelPart,
                array_1d<double,3>& rResultForce
                ):
        Process(),
        mrModelPart(rModelPart),
        mrResultForce(rResultForce)
    {
    }

void ComputeEmbeddedLiftProcess::Execute()
{
    KRATOS_TRY;

        double Cl=0;
        double Cd=0;
        double Rz=0;

        for(auto it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); ++it)
        {
            if (it->Is(TO_SPLIT) && it -> Is(ACTIVE)){
                auto geom = it->GetGeometry();
                const unsigned int NumNodes=geom.size();
                array_1d<double,NumNodes> elemental_distance;
                Geometry< Node<3> >::Pointer pgeom = it-> pGetGeometry();

                for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
                    elemental_distance[i_node] = geom[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

                const Vector& r_elemental_distances=elemental_distance;

                Triangle2D3ModifiedShapeFunctions triangle_shape_functions(pgeom, r_elemental_distances);


                array_1d<double,3>  gp_wall = ZeroVector(3);
                Matrix positive_side_interface_sh_func;
                ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_interface_gradients;
                Vector positive_side_interface_weights;
                triangle_shape_functions.ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                    positive_side_interface_sh_func,
                    positive_side_sh_func_interface_gradients,
                    positive_side_interface_weights,
                    GeometryData::GI_GAUSS_1);


                std::vector<Vector> cut_unit_normal;
                std::vector<Vector> gp_cut;
                triangle_shape_functions.ComputePositiveSideInterfaceAreaNormals(cut_unit_normal,GeometryData::GI_GAUSS_1);

                const unsigned int n_gauss = positive_side_interface_weights.size();

                std::vector<double> cp;


                for (unsigned int i_gauss=0;i_gauss<n_gauss;i_gauss++){
                    for (unsigned int i_node=0;i_node<NumNodes;i_node++){
                        array_1d<double,3> coord=geom[i_node].Coordinates();
                        for (unsigned int i_dim=0;i_dim<3;i_dim++){
                            gp_wall[i_dim] += coord[i_dim]*positive_side_interface_sh_func(i_gauss,i_node);//*positive_side_interface_weights(i_gauss);
                        }
                    }
                }

                it->GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());

                double cpressure=cp[0];

                Cl += cpressure*cut_unit_normal[0][1];
                Cd += cpressure*cut_unit_normal[0][0];
                Rz += cpressure*cut_unit_normal[0][2];


                it->SetValue(PRESSURE_COEFFICIENT,cp[0]);
                it->SetValue(NORMAL,cut_unit_normal[0]);
                it->SetValue(BODY_FORCE,gp_wall); //provisional name for cp application point


            }
        }
        mrResultForce[0]=Cd;
        mrResultForce[1]=Cl;
        mrResultForce[2]=Rz;

    KRATOS_CATCH("");
}

ModifiedShapeFunctions::Pointer ComputeEmbeddedLiftProcess::pGetModifiedShapeFunctions<3>(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

}// Namespace Kratos
