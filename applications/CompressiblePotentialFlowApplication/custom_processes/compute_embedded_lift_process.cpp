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
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
// Constructor for ComputeEmbeddedLiftProcess Process
ComputeEmbeddedLiftProcess::ComputeEmbeddedLiftProcess(ModelPart& rModelPart,
                Vector& rResultForce
                ):
        Process(),
        mrModelPart(rModelPart),
        mrResultForce(rResultForce)
    {
    }

void ComputeEmbeddedLiftProcess::Execute()
{
    KRATOS_TRY;

    mrResultForce = ZeroVector(3);

    #pragma omp parallel for
    for(std::size_t i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it=mrModelPart.ElementsBegin()+i;
        if (it->Is(TO_SPLIT) && it -> Is(ACTIVE)){
            auto r_geometry = it->GetGeometry();
            const std::size_t NumNodes = r_geometry.PointsNumber();

            array_1d<double,3> elemental_distances;
            for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
                elemental_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

            const Vector& r_elemental_distances=elemental_distances;
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it->pGetGeometry(), r_elemental_distances);

            // Computing Normal
            std::vector<Vector> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);

            std::vector<double> pressure_coefficient;
            it->GetValueOnIntegrationPoints(PRESSURE_COEFFICIENT,pressure_coefficient,mrModelPart.GetProcessInfo());

            //Storing the local cp and cut normal
            it->SetValue(PRESSURE_COEFFICIENT,pressure_coefficient[0]);
            it->SetValue(NORMAL,cut_normal[0]);

            //Calculating result force as the sum of the pressure contribution of every element
            for (std::size_t i = 0; i<3;i++){
                mrResultForce(i) += pressure_coefficient[0]*cut_normal[0][i];
            }
        }
    }

    KRATOS_CATCH("");
}

ModifiedShapeFunctions::Pointer ComputeEmbeddedLiftProcess::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
    GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();
    switch (geometry_type){
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
        default:
                KRATOS_ERROR << "Only Triangle2D3 geometries are currently implemented. The given geometry was: " << geometry_type;
    }
}

}// Namespace Kratos
