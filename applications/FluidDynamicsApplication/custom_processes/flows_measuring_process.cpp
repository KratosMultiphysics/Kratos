//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "flows_measuring_process.h"
#include "custom_utilities/post_process_utilities.h"


namespace Kratos
{


    /// Constructor.
    FlowsMeasuringProcess::FlowsMeasuringProcess(Model& rModel, Kratos::Parameters::Pointer pParameters) : mrModel(rModel)
    {
        Parameters default_parameters( R"(
        {
            "model_part_name" : "default_model_part_name",
            "list_of_outlet_submodelpart_names" : []
        }  )" );

        (*pParameters).ValidateAndAssignDefaults(default_parameters);

        mpModelPartContainingTime = &mrModel.GetModelPart((*pParameters)["model_part_name"].GetString());

        std::ofstream outputfile(mFilename, std::ios_base::out | std::ios_base::app);

        outputfile<<"       ";
        for( unsigned int i=0; i<(*pParameters)["list_of_outlet_submodelpart_names"].size(); ++i) {
            const std::string smp_name = (*pParameters)["list_of_outlet_submodelpart_names"][i].GetString();
            mListOfSubmodelparts.push_back(&mrModel.GetModelPart(smp_name));
            outputfile <<"  "<< smp_name;
        }
        outputfile<<"\n";
    }

    void FlowsMeasuringProcess::ExecuteFinalizeSolutionStep() {

        const double& time = mpModelPartContainingTime->GetProcessInfo()[TIME];
        std::ofstream outputfile(mFilename, std::ios_base::out | std::ios_base::app);
        outputfile << time;

        for(std::vector<ModelPart*>::iterator i_smp = mListOfSubmodelparts.begin(); i_smp != mListOfSubmodelparts.end(); ++i_smp) {
            double flow_of_this_outlet = 0.0;
            PostProcessUtilities().ComputeFlow(**i_smp, flow_of_this_outlet);
            outputfile << "  " << flow_of_this_outlet;
        }// for  (submodelparts)
        outputfile << "\n";
    }
}
