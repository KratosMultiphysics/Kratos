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
                Vector& rResultantForce
                ):
        Process(),
        mrModelPart(rModelPart),
        mrResultantForce(rResultantForce)
    {
    }

template <unsigned int Dim, unsigned int NumNodes>
void ComputeEmbeddedLiftProcess<Dim, NumNodes>::Execute()
{
    KRATOS_TRY;

    mrResultantForce = ZeroVector(3);

    std::ofstream outfile;
    outfile.open("cp_normals.txt");

    //Declaring auxilary variables needed to use with omp
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    // #pragma omp parallel for reduction(+:fx,fy,fz)
    for(int i = 0; i <  static_cast<int>(mrModelPart.NumberOfElements()); ++i) {
        auto it_elem=mrModelPart.ElementsBegin()+i;
        auto r_geometry = it_elem->GetGeometry();

        BoundedVector<double, NumNodes> geometry_distances;
        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            geometry_distances[i_node] = r_geometry[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(geometry_distances);

        if (is_embedded && it_elem->Is(ACTIVE)){
            it_elem->Set(SOLID);
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_distances));

            // Computing Normal
            std::vector<array_1d<double,3>> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);

            Matrix positive_side_sh_func;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
            Vector positive_side_weights;
            pModifiedShFunc -> ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                positive_side_sh_func,
                positive_side_sh_func_gradients,
                positive_side_weights,
            GeometryData::GI_GAUSS_1);

            it_elem->SetValue(VELOCITY_LOWER, cut_normal[0]);

            std::vector<double> pressure_coefficient;
            it_elem->CalculateOnIntegrationPoints(PRESSURE_COEFFICIENT,pressure_coefficient,mrModelPart.GetProcessInfo());

            outfile << it_elem->Id() << " " << it_elem->GetGeometry().Center().X() << " " << it_elem->GetGeometry().Center().Y() << " " << it_elem->GetGeometry().Center().Z() << " " << pressure_coefficient[0] << " " << cut_normal[0][0] << " " << cut_normal[0][1] << " " << cut_normal[0][2] << " " << positive_side_weights(0);
            outfile << "\n";
            //Storing the local cp and cut normal
            it_elem->SetValue(PRESSURE_COEFFICIENT,pressure_coefficient[0]);
            it_elem->SetValue(NORMAL,cut_normal[0]);

            fx += pressure_coefficient[0]*cut_normal[0][0];
            fy += pressure_coefficient[0]*cut_normal[0][1];
            fz += pressure_coefficient[0]*cut_normal[0][2];
        }
    }

    outfile.close();
    // Storing final result
    mrResultantForce[0] = fx;
    mrResultantForce[1] = fy;
    mrResultantForce[2] = fz;

    std::ofstream outfile_solid;
    outfile_solid.open("solid_elements_id.txt");
    for (auto& r_element : mrModelPart.Elements()){
        if(r_element.Is(SOLID)){
            if (r_element.GetGeometry().Center().X() > 0.4) {
                outfile_solid << r_element.Id();
                outfile_solid << " ";
            }
        }
    }
    outfile_solid.close();

    KRATOS_CATCH("");
}

template<>
ModifiedShapeFunctions::Pointer ComputeEmbeddedLiftProcess<2, 3>::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template<>
ModifiedShapeFunctions::Pointer ComputeEmbeddedLiftProcess<3, 4>::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rDistances);
}

template class ComputeEmbeddedLiftProcess<2, 3>;
template class ComputeEmbeddedLiftProcess<3, 4>;
}// Namespace Kratos
