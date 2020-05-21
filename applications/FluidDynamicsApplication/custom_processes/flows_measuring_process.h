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

#ifndef KRATOS_FLOWS_MEASURING_PROCESS_H
#define KRATOS_FLOWS_MEASURING_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class FlowsMeasuringProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FlowsMeasuringProcess
    KRATOS_CLASS_POINTER_DEFINITION(FlowsMeasuringProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FlowsMeasuringProcess(Model& rModel, Kratos::Parameters::Pointer pParameters) : mrModel(rModel)
    {

        Parameters default_parameters( R"(
        {
            "model_part_containing_time_name" : "default_model_part_name",
            "list_of_outlet_submodelpart_names" : []
        }  )" );

        (*pParameters).ValidateAndAssignDefaults(default_parameters);

        mpModelPartContainingTime = &mrModel.GetModelPart((*pParameters)["model_part_containing_time_name"].GetString());

        std::ofstream outputfile(mFilename, std::ios_base::out | std::ios_base::app);

        outputfile<<"       ";
        for( unsigned int i=0; i<(*pParameters)["list_of_outlet_submodelpart_names"].size(); ++i) {
            const std::string smp_name = (*pParameters)["list_of_outlet_submodelpart_names"][i].GetString();
            mListOfSubmodelparts.push_back(&mrModel.GetModelPart(smp_name));
            outputfile <<"  "<< smp_name;
        }
        outputfile<<"\n";

    }

    /// Destructor.
    ~FlowsMeasuringProcess() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteFinalizeSolutionStep() override {

        const double one_third = 1.0/3.0;

        const double& time = mpModelPartContainingTime->GetProcessInfo()[TIME];
        std::ofstream outputfile(mFilename, std::ios_base::out | std::ios_base::app);
        outputfile << time;

        for(std::vector<ModelPart*>::iterator i_smp = mListOfSubmodelparts.begin(); i_smp != mListOfSubmodelparts.end(); ++i_smp) {
            double flow_of_this_outlet = 0.0;
            for(ModelPart::ConditionsContainerType::iterator i_cond = (*i_smp)->ConditionsBegin(); i_cond != (*i_smp)->ConditionsEnd(); ++i_cond) {
                const auto& geom = i_cond->GetGeometry();
                const array_1d<double, 3> tangent_xi  = geom.GetPoint(1) - geom.GetPoint(0);
                const array_1d<double, 3> tangent_eta = geom.GetPoint(2) - geom.GetPoint(0);

                array_1d<double, 3> area_normal;
                MathUtils<double>::CrossProduct(area_normal, tangent_xi, tangent_eta);
                area_normal *= 0.5;

                const double area = MathUtils<double>::Norm3(area_normal);
                array_1d<double,3> unitary_normal;
                noalias(unitary_normal) = area_normal / area;

                double average_velocity = 0.0;

                for (int j=0; j<(int)geom.size(); j++) {
                    average_velocity += one_third * MathUtils<double>::Dot(geom[j].FastGetSolutionStepValue(VELOCITY), unitary_normal);
                }

                const double condition_flow = area * average_velocity;

                flow_of_this_outlet += condition_flow;
            } //for conditions

            outputfile << "  " << flow_of_this_outlet;
        }// for  (submodelparts)
        outputfile << "\n";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "FlowsMeasuringProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "FlowsMeasuringProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    Model& mrModel;
    ModelPart* mpModelPartContainingTime;
    std::vector<ModelPart*> mListOfSubmodelparts;
    std::string mFilename = "flow_data_output.txt";

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    FlowsMeasuringProcess() = delete;

    /// Assignment operator.
    FlowsMeasuringProcess& operator=(FlowsMeasuringProcess const& rOther) = delete;

    /// Copy constructor.
    FlowsMeasuringProcess(FlowsMeasuringProcess const& rOther) = delete;


    ///@}

}; // Class FlowsMeasuringProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_FLOWS_MEASURING_PROCESS_H
