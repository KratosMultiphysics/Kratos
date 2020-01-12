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
#include "adjoint_far_field_lift_response_function.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
    AdjointLiftFarFieldCoordinatesResponseFunction::AdjointLiftFarFieldCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointPotentialResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;

        // Reading the reference chord from the parameters
        mReferenceChord = ResponseSettings["reference_chord"].GetDouble();
        const double eps = std::numeric_limits<double>::epsilon();
        KRATOS_ERROR_IF(mReferenceChord < eps)
            << "The reference chord should be larger than 0." << mReferenceChord << std::endl;

    }

    AdjointLiftFarFieldCoordinatesResponseFunction::~AdjointLiftFarFieldCoordinatesResponseFunction(){}

    void AdjointLiftFarFieldCoordinatesResponseFunction::InitializeSolutionStep() {
        mUnperturbedLift = this->CalculateValue(mrModelPart);
    }

    double AdjointLiftFarFieldCoordinatesResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        ModelPart& root_model_part = rModelPart.GetRootModelPart();
        auto& far_field_model_part = root_model_part.GetSubModelPart("PotentialWallCondition2D_Far_field_Auto1");
        auto free_stream_velocity = root_model_part.GetProcessInfo()[FREE_STREAM_VELOCITY];
        auto free_stream_density = root_model_part.GetProcessInfo()[FREE_STREAM_DENSITY];
        double free_stream_velocity_norm = norm_2(free_stream_velocity);
        KRATOS_ERROR_IF(free_stream_velocity_norm<std::numeric_limits<double>::epsilon()) << "Free stream velocity is zero!" << std::endl;
        auto r_current_process_info = rModelPart.GetProcessInfo();

        array_1d<double, 3> force_coefficient_pres;
        force_coefficient_pres.clear();
        array_1d<double, 3> force_coefficient_vel;
        force_coefficient_vel.clear();
        for(int i = 0; i <  static_cast<int>(far_field_model_part.NumberOfConditions()); ++i) {
            auto it_cond=far_field_model_part.ConditionsBegin()+i;
            auto& r_geometry = it_cond->GetGeometry();
            // Computing normal
            array_1d<double,3> aux_coordinates;
            r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
            const auto normal = r_geometry.Normal(aux_coordinates);
            it_cond->Initialize(r_current_process_info); // TO-DO: Move the initialize outside this loop
            it_cond->FinalizeNonLinearIteration(r_current_process_info);
            double pressure_coefficient = it_cond-> GetValue(PRESSURE_COEFFICIENT);
            force_coefficient_pres -= normal*pressure_coefficient;

            auto velocity = it_cond->GetValue(VELOCITY);
            double density = it_cond->GetValue(DENSITY);
            double velocity_projection = inner_prod(normal,velocity);
            array_1d<double,3> disturbance = velocity - free_stream_velocity;
            force_coefficient_vel -= velocity_projection * disturbance * density;
        }

        force_coefficient_pres = force_coefficient_pres / mReferenceChord;

        auto wake_direction = free_stream_velocity/free_stream_velocity_norm;
        array_1d<double,3> wake_normal;
        wake_normal[0] = -wake_direction[1];
        wake_normal[1] = wake_direction[0];

        double free_stream_velocity_2 = inner_prod(free_stream_velocity, free_stream_velocity);
        double dynamic_pressure = 0.5 * free_stream_velocity_2 * free_stream_density;
        force_coefficient_vel = force_coefficient_vel/(dynamic_pressure*mReferenceChord);

        auto force_coefficient = force_coefficient_pres + force_coefficient_vel;

        double lift =  inner_prod(force_coefficient, wake_normal);

        return lift;

        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;


        rResponseGradient = ZeroVector(rResidualGradient.size1());
        auto this_id = rAdjointElement.Id();
        auto& r_this_element = mrModelPart.GetElement(this_id);
        auto& r_geometry = r_this_element.GetGeometry();
        double epsilon = 1e-9;
        bool is_cond = false;
        std::size_t counter = 0;
        for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){
            if (r_geometry[i_node].Is(INLET)){
                counter++;
            }
        }

        if (counter > r_geometry.WorkingSpaceDimension()-1){
            is_cond = true;
        }

        if (is_cond) {
            for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){

                //Perturbing potential
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += epsilon;
                double lift = this->CalculateValue(mrModelPart);
                rResponseGradient[i_node] = -(lift-mUnperturbedLift)/epsilon;
                //Unperturbing potential
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= epsilon;
            }
        }


        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rResponseGradient = ZeroVector(rResidualGradient.size1());
        auto this_id = rAdjointCondition.Id();
        auto& r_this_condition = mrModelPart.GetCondition(this_id);
        auto& r_geometry = r_this_condition.GetGeometry();
        double epsilon = 1e-9;
        for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){
                //Perturbing potential
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += epsilon;
                double lift = this->CalculateValue(mrModelPart);
                rResponseGradient[i_node] = (lift-mUnperturbedLift)/epsilon;
                //Unperturbing potential
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= epsilon;

        }
        KRATOS_WATCH(rResponseGradient)


        // rResponseGradient = ZeroVector(rResidualGradient.size1());
        // auto this_id = rAdjointCondition.Id();
        // auto& r_this_condition = mrModelPart.GetCondition(this_id);
        // auto& r_geometry = r_this_condition.GetGeometry();
        // double epsilon = 1e-6;
        // for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){
        //     if (r_geometry[i_node].IsNot(MARKER)){
        //         double unperturbed_lift = this->CalculateValue(mrModelPart);
        //         //Perturbing potential
        //         r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += epsilon;
        //         double lift = this->CalculateValue(mrModelPart);
        //         rResponseGradient[i_node] = (lift-unperturbed_lift)/epsilon;
        //         //Unperturbing potential
        //         r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= epsilon;

        //         r_geometry[i_node].Set(MARKER,true);
        //     }

        // }
        // KRATOS_WATCH(rResponseGradient)


        // double epsilon = 1e-9;
        // rResponseGradient = ZeroVector(rResidualGradient.size1());
        // auto this_id = rAdjointCondition.Id();
        // auto& r_this_condition = mrModelPart.GetCondition(this_id);

        // ModelPart& root_model_part = mrModelPart.GetRootModelPart();
        // auto free_stream_velocity = root_model_part.GetProcessInfo()[FREE_STREAM_VELOCITY];
        // auto free_stream_density = root_model_part.GetProcessInfo()[FREE_STREAM_DENSITY];
        // double free_stream_velocity_norm = norm_2(free_stream_velocity);
        // KRATOS_ERROR_IF(free_stream_velocity_norm<std::numeric_limits<double>::epsilon()) << "Free stream velocity is zero!" << std::endl;
        // auto r_current_process_info = mrModelPart.GetProcessInfo();
        // auto wake_direction = free_stream_velocity/free_stream_velocity_norm;
        // double free_stream_velocity_2 = inner_prod(free_stream_velocity, free_stream_velocity);
        // double dynamic_pressure = 0.5 * free_stream_velocity_2 * free_stream_density;

        // array_1d<double,3> wake_normal;
        // wake_normal[0] = -wake_direction[1];
        // wake_normal[1] = wake_direction[0];


        // auto& r_geometry = r_this_condition.GetGeometry();

        // for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node){

        //     if (r_geometry[i_node].Is(INLET)){



        //     array_1d<double, 3> force_coefficient_pres;
        //     force_coefficient_pres.clear();
        //     array_1d<double, 3> force_coefficient_vel;
        //     force_coefficient_vel.clear();
        //     // Computing normal
        //     array_1d<double,3> aux_coordinates;
        //     r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
        //     const auto normal = r_geometry.Normal(aux_coordinates);
        //     r_this_condition.Initialize(r_current_process_info); // TO-DO: Move the initialize outside this loop
        //     r_this_condition.FinalizeNonLinearIteration(r_current_process_info);
        //     double pressure_coefficient = r_this_condition.GetValue(PRESSURE_COEFFICIENT);
        //     force_coefficient_pres = normal*pressure_coefficient;

        //     auto velocity = r_this_condition.GetValue(VELOCITY);
        //     double density = r_this_condition.GetValue(DENSITY);
        //     double velocity_projection = inner_prod(normal,velocity);
        //     array_1d<double,3> disturbance = velocity - free_stream_velocity;
        //     force_coefficient_vel = velocity_projection * disturbance * density;
        //     force_coefficient_pres = force_coefficient_pres / mReferenceChord;
        //     force_coefficient_vel = force_coefficient_vel/(dynamic_pressure*mReferenceChord);

        //     auto force_coefficient = force_coefficient_pres + force_coefficient_vel;

        //     double lift =  inner_prod(force_coefficient, wake_normal);

        //     r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += epsilon;

        //     array_1d<double, 3> force_coefficient_pres_pert;
        //     force_coefficient_pres_pert.clear();
        //     array_1d<double, 3> force_coefficient_vel_pert;
        //     force_coefficient_vel_pert.clear();

        //     // Computing normal
        //     // r_this_condition.Initialize(r_current_process_info); // TO-DO: Move the initialize outside this loop
        //     r_this_condition.FinalizeNonLinearIteration(r_current_process_info);
        //     pressure_coefficient = r_this_condition.GetValue(PRESSURE_COEFFICIENT);
        //     force_coefficient_pres_pert = normal*pressure_coefficient;

        //     velocity = r_this_condition.GetValue(VELOCITY);
        //     density = r_this_condition.GetValue(DENSITY);
        //     velocity_projection = inner_prod(normal,velocity);
        //     disturbance = velocity - free_stream_velocity;
        //     force_coefficient_vel_pert = velocity_projection * disturbance * density;
        //     force_coefficient_pres_pert = force_coefficient_pres_pert / mReferenceChord;
        //     force_coefficient_vel_pert = force_coefficient_vel_pert/(dynamic_pressure*mReferenceChord);

        //     auto force_coefficient_pert = force_coefficient_pres_pert + force_coefficient_vel_pert;

        //     double perturbed_lift =  inner_prod(force_coefficient_pert, wake_normal);

        //     r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= epsilon;

        //     rResponseGradient[i_node] = (perturbed_lift-lift)/epsilon;
        //     }

        // }
        // KRATOS_WATCH(rResponseGradient)

        KRATOS_CATCH("");
    }


    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLiftFarFieldCoordinatesResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        for (auto elem_it = mrModelPart.Elements().ptr_begin(); elem_it != mrModelPart.Elements().ptr_end(); ++elem_it)
        {
            if ((*elem_it)->Is(STRUCTURE)){
                mpNeighboringElement = (*elem_it);
                return;
            }
        }
        KRATOS_ERROR << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");
    }
} // namespace Kratos.


