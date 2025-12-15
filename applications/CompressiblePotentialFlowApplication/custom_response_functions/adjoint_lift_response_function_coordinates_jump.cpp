//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on Martin Fusseder work, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_lift_response_function_coordinates_jump.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_processes/compute_nodal_value_process.h"

namespace Kratos
{
    AdjointLiftJumpCoordinatesResponseFunction::AdjointLiftJumpCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointPotentialResponseFunction(rModelPart, ResponseSettings)
    {
        const Parameters default_parameters = Parameters(R"({
            "trailing_edge_model_part_name" : "",
            "reference_chord"               : 0.0,
            "step_size"                     : 0.0,
            "gradient_mode"                 : "semi_analytic",
            "response_type"                 : "adjoint_lift_jump_coordinates"
        })");
        ResponseSettings.ValidateAndAssignDefaults(default_parameters);
        // Reading the reference chord from the parameters
        mReferenceChord = ResponseSettings["reference_chord"].GetDouble();
        const double eps = std::numeric_limits<double>::epsilon();
        KRATOS_ERROR_IF(mReferenceChord < eps)
            << "The reference chord should be larger than 0." << mReferenceChord << std::endl;
        mTrailingEdgeModelPartName = ResponseSettings["trailing_edge_model_part_name"].GetString();
    }

    AdjointLiftJumpCoordinatesResponseFunction::~AdjointLiftJumpCoordinatesResponseFunction(){}

    void AdjointLiftJumpCoordinatesResponseFunction::InitializeSolutionStep() {
        // Get pointer to element that contains the traced node
        this->GetNeighboringElementPointer();
    }

    double AdjointLiftJumpCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        double lift_coefficient = 0.0;
        const array_1d<double, 3> free_stream_velocity = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_VELOCITY);
        double free_stream_velocity_norm = norm_2(free_stream_velocity);
        const int dim = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        if (dim == 2)
        {
            Element tracedElem = rModelPart.GetElement(mTrailingEdgeElementIds[0]);
            unsigned int NumNodes = tracedElem.GetGeometry().size();
            for(IndexType i = 0; i < NumNodes; ++i)
            {
                if(tracedElem.GetGeometry()[i].GetValue(TRAILING_EDGE))
                {
                    double potential = tracedElem.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                    double aux_potential = tracedElem.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                    double potential_jump = std::abs(potential - aux_potential);
                    lift_coefficient = 2.0 * potential_jump / (free_stream_velocity_norm * mReferenceChord);
                }
            }
        } else if (dim == 3)
        {
            KRATOS_ERROR_IF(mTrailingEdgeModelPartName=="") << "AdjointLiftJumpCoordinatesResponseFunction: Define the 'trailing_edge_model_part_name' in 'response_function_settings'." << std::endl;
            const std::vector<std::string> variable_array = {"VELOCITY"};
            ComputeNodalValueProcess ComputeNodalValueProcess(mrModelPart, variable_array);
            ComputeNodalValueProcess.Execute();
            double potential_integral = 0.0;
            for (auto cond_it = mrModelPart.GetModel().GetModelPart(mTrailingEdgeModelPartName).Conditions().ptr_begin(); 
                cond_it != mrModelPart.GetModel().GetModelPart(mTrailingEdgeModelPartName).Conditions().ptr_end(); ++cond_it)
            {
                double length = (*cond_it)->GetGeometry().Area();
                unsigned int NumNodes = (*cond_it)->GetGeometry().size();
                for (IndexType i = 0; i < NumNodes; ++i)
                {
                    double potential      = (*cond_it)->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                    double aux_potential  = (*cond_it)->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                    auto velocity         = (*cond_it)->GetGeometry()[i].GetValue(VELOCITY);
                    double potential_jump = potential - aux_potential;
                    potential_integral   += 0.5 * length * potential_jump;
                }
            }
            lift_coefficient = 2 * potential_integral / (free_stream_velocity_norm * mReferenceChord);
        } else
        {
            KRATOS_ERROR << "Incorrect ModelPart's DOMAIN_SIZE." << std::endl;
        }
        return lift_coefficient;
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
        for (SizeType i = 0; i < mTrailingEdgeElementIds.size(); i++) {
            if (rAdjointElement.Id() == mTrailingEdgeElementIds[i]) {
                unsigned int NumNodes = rAdjointElement.GetGeometry().size();
                const array_1d<double, 3> v_inf = rProcessInfo.GetValue(FREE_STREAM_VELOCITY);
                const int dim = rProcessInfo.GetValue(DOMAIN_SIZE);
                double v_norm = norm_2(v_inf);
                double derivative = 2.0 / (v_norm * mReferenceChord);
                array_1d<double, 2> normal_dir = ZeroVector(2);
                if (dim == 2) {
                    normal_dir[0] = -v_inf[1] / v_norm;
                    normal_dir[1] =  v_inf[0] / v_norm;
                } else {// dim == 3
                    normal_dir[0] = -v_inf[2] / v_norm;
                    normal_dir[1] =  v_inf[0] / v_norm;
                }
                for (IndexType i = 0; i < NumNodes; ++i) {
                    if (rAdjointElement.GetGeometry()[i].GetValue(TRAILING_EDGE)) {
                        rResponseGradient[i + 0*NumNodes] =  derivative * normal_dir[0];
                        rResponseGradient[i + 1*NumNodes] = -derivative * normal_dir[1];
                    }
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }


    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;
        mTrailingEdgeElementIds.clear();
        for (auto elem_it = mrModelPart.GetRootModelPart().GetSubModelPart("trailing_edge_elements_model_part").Elements().ptr_begin(); 
            elem_it != mrModelPart.GetRootModelPart().GetSubModelPart("trailing_edge_elements_model_part").Elements().ptr_end(); ++elem_it)
        {
            if ((*elem_it)->Is(STRUCTURE)){
                mTrailingEdgeElementIds.push_back((*elem_it)->Id());
            }
        }
        if (mTrailingEdgeElementIds.size() == 0) {
            KRATOS_ERROR << "No elements found in submodel part: trailing_edge_elements_model_part" << std::endl;
        }
        KRATOS_CATCH("");
    }
} // namespace Kratos.


