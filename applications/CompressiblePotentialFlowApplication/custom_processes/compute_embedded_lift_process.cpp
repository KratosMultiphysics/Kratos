//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez
//


#include "compute_embedded_lift_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
// Constructor for ComputeEmbeddedLiftProcess Process
template <unsigned int Dim, unsigned int NumNodes>
ComputeEmbeddedLiftProcess<Dim, NumNodes>::ComputeEmbeddedLiftProcess(ModelPart& rModelPart,
                Vector& rResultForce
                ):
        Process(),
        mrModelPart(rModelPart),
        mrResultForce(rResultForce)
    {
    }

template <unsigned int Dim, unsigned int NumNodes>
void ComputeEmbeddedLiftProcess<Dim, NumNodes>::Execute()
{
    KRATOS_TRY;

    mrResultForce = ZeroVector(3);

    //Declaring auxilary variables needed to use with omp
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    #pragma omp parallel for reduction(+:fx,fy,fz)
    for(std::size_t i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem=mrModelPart.ElementsBegin()+i;
        auto r_geometry = it_elem->GetGeometry();

        BoundedVector<double,3> geometry_distances;
        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            geometry_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(geometry_distances);

        if (is_embedded && it_elem->Is(ACTIVE)){
            array_1d<double,3> elemental_distances;
            for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
                elemental_distances[i_node] = r_geometry[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
            }

            const Vector& r_elemental_distances=elemental_distances;
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), r_elemental_distances);

            // Computing Normal
            std::vector<Vector> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);

            std::vector<double> pressure_coefficient;
            it_elem->GetValueOnIntegrationPoints(PRESSURE_COEFFICIENT,pressure_coefficient,mrModelPart.GetProcessInfo());

            //Storing the local cp and cut normal
            it_elem->SetValue(PRESSURE_COEFFICIENT,pressure_coefficient[0]);
            it_elem->SetValue(NORMAL,cut_normal[0]);

            fx += pressure_coefficient[0]*cut_normal[0][0];
            fy += pressure_coefficient[0]*cut_normal[0][1];
            fz += pressure_coefficient[0]*cut_normal[0][2];
        }
    }

    // Storing final result
    mrResultForce[0] = fx;
    mrResultForce[1] = fy;
    mrResultForce[2] = fz;

    KRATOS_CATCH("");
}

template<unsigned int Dim, unsigned int NumNodes>
ModifiedShapeFunctions::Pointer ComputeEmbeddedLiftProcess<Dim, NumNodes>::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
    GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();
    switch (geometry_type){
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
        default:
                KRATOS_ERROR << "Only Triangle2D3 geometries are currently implemented. The given geometry was: " << geometry_type;
    }
}

template class ComputeEmbeddedLiftProcess<2, 3>;
}// Namespace Kratos
