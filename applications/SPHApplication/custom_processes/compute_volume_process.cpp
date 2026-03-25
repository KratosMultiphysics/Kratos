// SPH Application 

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once

#include "custom_processes/compute_volume_process.h"


namespace Kratos
{

void ComputeVolumeProcess::Execute()
{
    KRATOS_TRY

    typedef Geometry<Node>::CoordinatesArrayType CoordinatesArrayType;

    auto& rGeoms = mrThisModelPart.Geometries();
    auto& rNodes = mrThisModelPart.Nodes();
    
    SizeType domain_size = 2;
    
    
    for (auto& node : rNodes){
        node.SetValue(VOLUME, 0.0);
        node.SetValue(BOUNDARY_NORMAL_AREA, ZeroVector(domain_size));
    }

    for (auto& r_geom : rGeoms){

        double volume = r_geom.DomainSize();
        int number_of_points = r_geom.PointsNumber();

        if (!mrThisParameters["structured_mesh"].GetBool()){
            ComputeVolumeUtilities::CalculateBoundaryNetNormal(r_geom);
        }
        
        double volume_aux = volume / number_of_points;
        for (auto& node : r_geom.Points()){
            node.GetValue(VOLUME) += volume_aux;
        }
    }

    if (!mrThisParameters["structured_mesh"].GetBool()){
        ComputeVolumeUtilities::CheckBoundaryNetNormal(rNodes);
    }

    KRATOS_CATCH("")
}

void ComputeVolumeProcess::ExecuteInitialize(){
    this->Execute();
}

}