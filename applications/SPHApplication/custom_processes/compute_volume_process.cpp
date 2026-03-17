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
    auto& rGeoms = mrThisModelPart.Geometries();
    auto& rNodes = mrThisModelPart.Nodes();
    
    
    for (auto& node : rNodes){
        node.SetValue(VOLUME, 0.0);
    }

    for (auto& geom : rGeoms){

        double volume = geom.DomainSize();
        int number_of_points = geom.PointsNumber();
        
        double volume_aux = volume / number_of_points;

        for (auto& node : geom.Points()){
            node.GetValue(VOLUME) += volume_aux;
        }
    }
    
    KRATOS_CATCH("")
}

void ComputeVolumeProcess::ExecuteInitialize(){
    this->Execute();
}

}