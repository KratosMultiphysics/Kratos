//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nunez, Inigo Lopez
//


// System includes

// External includes

// Project includes
#include "adjoint_far_field_lift_response_function.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
    AdjointLiftFarFieldResponseFunction::AdjointLiftFarFieldResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointPotentialResponseFunction(rModelPart, ResponseSettings)
    {

        const Parameters default_parameters = Parameters(R"(
        {
            "reference_chord"             : 1.0,
            "far_field_model_part_name"   : "",
            "wake_normal"                 : [0,1,0],
            "analyzer"                    : "kratos",
            "response_type"               : "adjoint_lift_far_field",
            "gradient_mode"               : "semi_analytic",
            "step_size"                   : 1e-6
        })" );
        ResponseSettings.ValidateAndAssignDefaults(default_parameters);

        KRATOS_ERROR_IF(ResponseSettings["far_field_model_part_name"].GetString()=="") << "Please define the far "
            "field model part name in the response settings \"far_field_model_part_name\"!" << std::endl;
        mFarFieldModelPartName = ResponseSettings["far_field_model_part_name"].GetString();

        mReferenceChord = ResponseSettings["reference_chord"].GetDouble();
        KRATOS_ERROR_IF(mReferenceChord < std::numeric_limits<double>::epsilon())
            << "The reference chord should be larger than 0." << mReferenceChord << std::endl;

        mWakeNormal = ResponseSettings["wake_normal"].GetVector();
        mStepSize = ResponseSettings["step_size"].GetDouble();

    }

    AdjointLiftFarFieldResponseFunction::~AdjointLiftFarFieldResponseFunction(){}

    void AdjointLiftFarFieldResponseFunction::InitializeSolutionStep() {
        VariableUtils().SetNonHistoricalVariableToZero(NORMAL, mrModelPart.Elements());
        mFreeStreamVelocity = mrModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY];
        double free_stream_velocity_norm = norm_2(mFreeStreamVelocity);
        KRATOS_ERROR_IF(free_stream_velocity_norm<std::numeric_limits<double>::epsilon()) << "Free stream velocity is zero!" << std::endl;
        auto wake_direction = mFreeStreamVelocity/free_stream_velocity_norm;
        double free_stream_velocity_2 = inner_prod(mFreeStreamVelocity, mFreeStreamVelocity);
        auto free_stream_density = mrModelPart.GetProcessInfo()[FREE_STREAM_DENSITY];
        mDynamicPressure = 0.5 * free_stream_velocity_2 * free_stream_density;
        mCurrentLift = this->CalculateValue(mrModelPart);
    }

    double AdjointLiftFarFieldResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        auto& far_field_sub_model_part = rModelPart.GetRootModelPart().GetSubModelPart(mFarFieldModelPartName);
        const auto r_current_process_info = rModelPart.GetProcessInfo();

        array_1d<double, 3> force_coefficient_pres;
        force_coefficient_pres.clear();
        array_1d<double, 3> force_coefficient_vel;
        force_coefficient_vel.clear();

        typedef CombinedReduction<SumReduction<array_1d<double, 3>>,SumReduction<array_1d<double, 3>>> CustomReduction;
        std::tie(force_coefficient_pres, force_coefficient_vel) =
        block_for_each<CustomReduction>(far_field_sub_model_part.Conditions(), [&](Condition& rCondition) {
            auto& r_geometry = rCondition.GetGeometry();
            array_1d<double, 3> local_force_coefficient_pres;
            array_1d<double, 3> local_force_coefficient_vel;

            // Computing normal
            array_1d<double,3> aux_coordinates;
            r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
            const auto normal = r_geometry.Normal(aux_coordinates);
            rCondition.SetValue(NORMAL, normal);

            rCondition.Initialize(r_current_process_info);
            rCondition.FinalizeNonLinearIteration(r_current_process_info);
            double pressure_coefficient = rCondition.GetValue(PRESSURE_COEFFICIENT);
            local_force_coefficient_pres = -normal*pressure_coefficient;

            auto velocity = rCondition.GetValue(VELOCITY);
            double density = rCondition.GetValue(DENSITY);
            double velocity_projection = inner_prod(normal,velocity);
            array_1d<double,3> disturbance = velocity - mFreeStreamVelocity;
            local_force_coefficient_vel = -velocity_projection * disturbance * density;
            return std::make_tuple(local_force_coefficient_pres, local_force_coefficient_vel);
        });

        force_coefficient_pres = force_coefficient_pres / mReferenceChord;
        force_coefficient_vel = force_coefficient_vel/(mDynamicPressure*mReferenceChord);

        auto force_coefficient = force_coefficient_pres + force_coefficient_vel;

        return inner_prod(force_coefficient, mWakeNormal);

        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rResponseGradient = ZeroVector(rResidualGradient.size1());
        auto& r_adjoint_geometry = rAdjointElement.GetGeometry();

        std::size_t far_field_nodes = 0;
        for (std::size_t i_node=0; i_node<r_adjoint_geometry.size(); ++i_node){
            if (r_adjoint_geometry[i_node].GetValue(FAR_FIELD)){
                far_field_nodes++;
            }
        }

        if (far_field_nodes > (r_adjoint_geometry.WorkingSpaceDimension()-1)) {
            auto& r_this_element = mrModelPart.GetElement(rAdjointElement.Id());
            auto& r_geometry = r_this_element.GetGeometry();
            for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){

                double lift = this->ComputeLiftContribution(r_this_element, rProcessInfo);

                if (rAdjointElement.GetValue(WAKE) && (r_geometry[i_node].GetValue(WAKE_DISTANCE) < 0.0)) {
                    r_geometry[i_node].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += mStepSize;
                }
                else{
                    r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += mStepSize;
                }

                double perturbed_lift = this->ComputeLiftContribution(r_this_element, rProcessInfo);

                if (rAdjointElement.GetValue(WAKE) && (r_geometry[i_node].GetValue(WAKE_DISTANCE) < 0.0)) {
                    r_geometry[i_node].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= mStepSize;
                }
                else {
                    r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= mStepSize;
                }

                rResponseGradient[i_node] = (perturbed_lift-lift)/mStepSize;
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }


    void AdjointLiftFarFieldResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointLiftFarFieldResponseFunction::ComputeLiftContribution(Element& rElement, const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        array_1d<double, 3> force_coefficient_pres;
        force_coefficient_pres.clear();
        array_1d<double, 3> force_coefficient_vel;
        force_coefficient_vel.clear();
        const auto normal = rElement.GetValue(NORMAL);

        std::vector<double> pressure_coefficient_vector;
        rElement.CalculateOnIntegrationPoints(PRESSURE_COEFFICIENT,pressure_coefficient_vector, rProcessInfo);
        double pressure_coefficient = pressure_coefficient_vector[0];
        force_coefficient_pres = -normal*pressure_coefficient/ mReferenceChord;;

        std::vector<array_1d<double,3>> velocity_vector;
        rElement.CalculateOnIntegrationPoints(VELOCITY,velocity_vector, rProcessInfo);
        array_1d<double,3> velocity = velocity_vector[0];
        std::vector<double> density_vector;
        rElement.CalculateOnIntegrationPoints(DENSITY,density_vector, rProcessInfo);
        double density = density_vector[0];
        double velocity_projection = inner_prod(normal,velocity);
        array_1d<double,3> disturbance = velocity - mFreeStreamVelocity;
        force_coefficient_vel = -velocity_projection * disturbance * density/(mDynamicPressure*mReferenceChord);

        auto force_coefficient = force_coefficient_pres + force_coefficient_vel;

        return inner_prod(force_coefficient, mWakeNormal);

        KRATOS_CATCH("");
    }
} // namespace Kratos.


