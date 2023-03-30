//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_CONSTRUCTION_UTILITIES)
#define KRATOS_CONSTRUCTION_UTILITIES

// System includes
#include <cmath>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class ConstructionUtility
{

  public:
    typedef std::size_t IndexType;
    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ConstructionUtility(ModelPart &rMechanicalModelPart, ModelPart &rThermalModelPart,
                        TableType &rTableAmbientTemp, Parameters &rParameters) : mrMechanicalModelPart(rMechanicalModelPart), mrThermalModelPart(rThermalModelPart), mrTableAmbientTemp(rTableAmbientTemp)
    {
        KRATOS_TRY

        // Getting values
        mGravityDirection = rParameters["gravity_direction"].GetString();
        mReferenceCoordinate = mHighestBlockHeight = rParameters["reservoir_bottom_coordinate_in_gravity_direction"].GetDouble();
        mLiftHeight = rParameters["lift_height"].GetDouble();
        mSourceType = rParameters["source_type"].GetString();
        mAging = rParameters["aging"].GetBool();
        mH0 = rParameters["h_0"].GetDouble();
        mActivateSoilPart = rParameters["activate_soil_part"].GetBool();
        mActivateExistingPart = rParameters["activate_existing_part"].GetBool();
        mTimeUnitConverter = mrMechanicalModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];

        if (mActivateSoilPart == true)
        {
            mMechanicalSoilPart = rParameters["mechanical_soil_part"].GetString();
            mThermalSoilPart = rParameters["thermal_soil_part"].GetString();
        }
        if (mActivateExistingPart == true)
        {
            mMechanicalExistingPart = rParameters["mechanical_existing_part"].GetString();
            mThermalExistingPart = rParameters["thermal_existing_part"].GetString();
        }

        if (mSourceType == "NonAdiabatic")
            mAlphaInitial = rParameters["alpha_initial"].GetDouble();

        if (mAging == true)
            mYoungInf = rParameters["young_inf"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ConstructionUtility() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize()
    {
        KRATOS_TRY;

        const int nelements = mrMechanicalModelPart.GetMesh(0).Elements().size();

        mMechanicalLastCondition = mrMechanicalModelPart.GetMesh(0).Conditions().size();
        mThermalLastCondition = mrThermalModelPart.GetMesh(0).Conditions().size();

        const ProcessInfo& r_current_process_info = mrMechanicalModelPart.GetProcessInfo();

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.ElementsBegin();
            mNumNode = el_begin->GetGeometry().PointsNumber();

            // Deactivate all elements of both thermal and mechanical model parts
            VariableUtils().SetFlag(ACTIVE, false, mrThermalModelPart.Elements());
            VariableUtils().SetFlag(ACTIVE, false, mrMechanicalModelPart.Elements());

            // Deactivate all nodes of one of the model parts (thermal) as both model parts share the same nodes
            VariableUtils().SetFlag(ACTIVE, false, mrThermalModelPart.Nodes());
            VariableUtils().SetFlag(SOLID, false, mrThermalModelPart.Nodes());

            block_for_each(mrMechanicalModelPart.Elements(), [&](auto& rMechanicalElement){
                // Initialize non-active elements of the mechanical model part (to avoid segmentation fault)
                rMechanicalElement.Initialize(r_current_process_info);
            });
        }

        // Activation of the existing parts, either the soil or the already built dam ( User must specify each part through the interface)
        if (mActivateSoilPart == true)
        {
            const int soil_nelements = mrMechanicalModelPart.GetSubModelPart(mMechanicalSoilPart).Elements().size();

            if (soil_nelements != 0)
            {
                // Activate elements of the soil part of both model parts
                VariableUtils().SetFlag(ACTIVE, true, mrThermalModelPart.GetSubModelPart(mThermalSoilPart).Elements());
                VariableUtils().SetFlag(ACTIVE, true, mrMechanicalModelPart.GetSubModelPart(mMechanicalSoilPart).Elements());

                // Activate nodes of the soil part of one of the model parts (thermal) as both model parts share the same nodes
                VariableUtils().SetFlag(ACTIVE, true, mrThermalModelPart.GetSubModelPart(mThermalSoilPart).Nodes());
                VariableUtils().SetFlag(SOLID, true, mrThermalModelPart.GetSubModelPart(mThermalSoilPart).Nodes());
            }
        }

        if (mActivateExistingPart == true)
        {
            const int existing_nelements = mrMechanicalModelPart.GetSubModelPart(mMechanicalExistingPart).Elements().size();

            if (existing_nelements != 0)
            {
                // Activate elements of the existing part of both model parts
                VariableUtils().SetFlag(ACTIVE, true, mrThermalModelPart.GetSubModelPart(mThermalExistingPart).Elements());
                VariableUtils().SetFlag(ACTIVE, true, mrMechanicalModelPart.GetSubModelPart(mMechanicalExistingPart).Elements());

                // Activate nodes of the existing part of one of the model parts (thermal) as both model parts share the same nodes
                VariableUtils().SetFlag(ACTIVE, true, mrThermalModelPart.GetSubModelPart(mThermalExistingPart).Nodes());
                VariableUtils().SetFlag(SOLID, true, mrThermalModelPart.GetSubModelPart(mThermalExistingPart).Nodes());
            }
        }

        // Mechanical Conditions
        const int nconditions_mech = mrMechanicalModelPart.GetMesh(0).Conditions().size();

        if (nconditions_mech != 0)
        {
            ModelPart::ConditionsContainerType::iterator cond_begin_mech = mrMechanicalModelPart.ConditionsBegin();

            for (int k = 0; k < nconditions_mech; ++k)
            {
                ModelPart::ConditionsContainerType::iterator it_cond_mech = cond_begin_mech + k;
                const unsigned int number_of_points = (*it_cond_mech).GetGeometry().PointsNumber();
                bool active_condition = true;

                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                {
                    if ((*it_cond_mech).GetGeometry()[i_node].IsNot(ACTIVE))
                    {
                        active_condition = false;
                        break;
                    }
                }
                if (active_condition) it_cond_mech->Set(ACTIVE, true);
                else it_cond_mech->Set(ACTIVE, false);
            }
        }

        // Thermal Conditions
        const int nconditions_thermal = mrThermalModelPart.GetMesh(0).Conditions().size();

        if (nconditions_thermal != 0)
        {
            ModelPart::ConditionsContainerType::iterator cond_begin_thermal = mrThermalModelPart.ConditionsBegin();

            for (int k = 0; k < nconditions_thermal; ++k)
            {
                ModelPart::ConditionsContainerType::iterator it_cond_thermal = cond_begin_thermal + k;
                const unsigned int number_of_points = (*it_cond_thermal).GetGeometry().PointsNumber();
                bool active_condition = true;

                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                {
                    if ((*it_cond_thermal).GetGeometry()[i_node].IsNot(ACTIVE))
                    {
                        active_condition = false;
                        break;
                    }
                }
                if (active_condition) it_cond_thermal->Set(ACTIVE, true);
                else it_cond_thermal->Set(ACTIVE, false);
            }
        }

        // Assign Alpha Initial in case of using Azenha Formulation
        if (mSourceType == "NonAdiabatic")
        {
            if (mAging == false)
            {
                VariableUtils().SetVariable(ALPHA_HEAT_SOURCE, mAlphaInitial, mrThermalModelPart.Nodes());
            }
            else
            {
                VariableUtils().SetVariable(ALPHA_HEAT_SOURCE, mAlphaInitial, mrThermalModelPart.Nodes());
                VariableUtils().SetVariable(NODAL_YOUNG_MODULUS, sqrt(mAlphaInitial) * mYoungInf, mrThermalModelPart.Nodes());
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AssignTimeActivation(std::string ThermalSubModelPartName, int phase, double time_activation, double initial_temperature)
    {
        KRATOS_TRY;

        const int nelements = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements().size();
        int direction;
        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        if (nelements != 0)
        {
            double current_height = mReferenceCoordinate + mLiftHeight * (phase);
            double previous_height = mReferenceCoordinate + mLiftHeight * (phase - 1);

            block_for_each(mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements(), [&](auto& rThermalElement){
                array_1d<double, 3> central_position = rThermalElement.GetGeometry().Center();
                if ((central_position(direction) >= previous_height) && (central_position(direction) <= current_height))
                {
                    const unsigned int number_of_points = rThermalElement.GetGeometry().PointsNumber();
                    for (unsigned int i = 0; i < number_of_points; ++i)
                    {
                        if (rThermalElement.GetGeometry()[i].FastGetSolutionStepValue(TIME_ACTIVATION)==0)
                        {
                            rThermalElement.GetGeometry()[i].FastGetSolutionStepValue(TIME_ACTIVATION) = time_activation * mTimeUnitConverter;
                            rThermalElement.GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE) = rThermalElement.GetGeometry()[i].FastGetSolutionStepValue(PLACEMENT_TEMPERATURE) = initial_temperature;
                        }
                    }
                }
            });
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(std::string ThermalSubModelPartName, std::string MechanicalSubModelPartName, std::string HeatFluxSubModelPartName, std::string HydraulicPressureSubModelPartName, bool thermal_conditions, bool mechanical_conditions, int current_number_of_phase)
    {
        KRATOS_TRY;

        const int nelements_thermal = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements().size();
        const int nelements_mech = mrMechanicalModelPart.GetSubModelPart(MechanicalSubModelPartName).Elements().size();

        int nconditions_thermal = 0;
        int nconditions_mech = 0;
        if (thermal_conditions) nconditions_thermal = mrThermalModelPart.GetSubModelPart(HeatFluxSubModelPartName).Conditions().size();
        if (mechanical_conditions) nconditions_mech = mrMechanicalModelPart.GetSubModelPart(HydraulicPressureSubModelPartName).Conditions().size();

        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        // Getting the value of the table and computing the current height
        double current_height = mReferenceCoordinate + mLiftHeight * current_number_of_phase;

        if (current_height > mHighestBlockHeight) mHighestBlockHeight = current_height;

        if (nelements_thermal != 0)
        {
            // Thermal Elements
            block_for_each(mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements(), [&](auto& rThermalElement){
                array_1d<double, 3> central_position = rThermalElement.GetGeometry().Center();
                if ((central_position(direction) >= (mReferenceCoordinate - mLiftHeight)) && (central_position(direction) <= current_height))
                {
                    rThermalElement.Set(ACTIVE, true);
                }
            });
        }

        if (nelements_mech != 0)
        {
            // Mechanical Elements
            block_for_each(mrMechanicalModelPart.GetSubModelPart(MechanicalSubModelPartName).Elements(), [&](auto& rMechanicalElement){
                array_1d<double, 3> central_position = rMechanicalElement.GetGeometry().Center();
                if ((central_position(direction) >= (mReferenceCoordinate - mLiftHeight)) && (central_position(direction) <= current_height))
                {
                    rMechanicalElement.Set(ACTIVE, true);

                    const unsigned int number_of_points = rMechanicalElement.GetGeometry().PointsNumber();
                    for (unsigned int i = 0; i < number_of_points; i++)
                    {
                        rMechanicalElement.GetGeometry()[i].Set(ACTIVE, true);
                        rMechanicalElement.GetGeometry()[i].Set(SOLID, false);
                    }
                }
            });
        }

        if (thermal_conditions && nconditions_thermal != 0)
        {
            // Thermal Conditions
            block_for_each(mrThermalModelPart.GetSubModelPart(HeatFluxSubModelPartName).Conditions(), [&](auto& rThermalCondition){
                array_1d<double, 3> central_position = rThermalCondition.GetGeometry().Center();
                if ((central_position(direction) >= (mReferenceCoordinate - mLiftHeight)) && (central_position(direction) <= current_height))
                {
                    if (rThermalCondition.IsNot(ACTIVE)) rThermalCondition.Set(ACTIVE, true);
                }
            });
        }

        if (mechanical_conditions && nconditions_mech != 0)
        {
            // Mechanical Conditions
            block_for_each(mrMechanicalModelPart.GetSubModelPart(HydraulicPressureSubModelPartName).Conditions(), [&](auto& rMechanicalCondition){
                array_1d<double, 3> central_position = rMechanicalCondition.GetGeometry().Center();
                if ((central_position(direction) >= (mReferenceCoordinate - mLiftHeight)) && (central_position(direction) <= current_height))
                {
                    if (rMechanicalCondition.IsNot(ACTIVE)) rMechanicalCondition.Set(ACTIVE, true);
                }
            });
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CheckTemperature(Parameters &CheckTemperatureParameters)
    {
        KRATOS_TRY;

        // Getting CheckTemperature Values
        const double maximum_temperature_increment = CheckTemperatureParameters["maximum_temperature_increment"].GetDouble();
        const double maximum_temperature_aux = CheckTemperatureParameters["maximum_temperature"].GetDouble();
        const double minimum_temperature_aux = CheckTemperatureParameters["minimum_temperature"].GetDouble();

        block_for_each(mrThermalModelPart.Nodes(), [&](auto& rNode){
            if (rNode.Is(ACTIVE) && rNode.IsNot(SOLID))
            {
                double maximum_temperature = std::max(rNode.FastGetSolutionStepValue(PLACEMENT_TEMPERATURE) + maximum_temperature_increment, maximum_temperature_aux);
                double minimum_temperature = std::min(rNode.FastGetSolutionStepValue(PLACEMENT_TEMPERATURE), minimum_temperature_aux);
                double current_temperature = rNode.FastGetSolutionStepValue(TEMPERATURE);

                if (current_temperature > maximum_temperature)
                {
                    rNode.FastGetSolutionStepValue(TEMPERATURE) = maximum_temperature;
                }
                else if (current_temperature < minimum_temperature)
                {
                    rNode.FastGetSolutionStepValue(TEMPERATURE) = minimum_temperature;
                }
            }
            else if (rNode.Is(ACTIVE) && rNode.Is(SOLID))
            {
                double maximum_temperature = maximum_temperature_aux;
                double minimum_temperature = minimum_temperature_aux;
                double current_temperature = rNode.FastGetSolutionStepValue(TEMPERATURE);

                if (current_temperature > maximum_temperature)
                {
                    rNode.FastGetSolutionStepValue(TEMPERATURE) = maximum_temperature;
                }
                else if (current_temperature < minimum_temperature)
                {
                    rNode.FastGetSolutionStepValue(TEMPERATURE) = minimum_temperature;
                }
            }
        });

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AfterOutputStep()
    {
        KRATOS_TRY;

        const int nelements = mrThermalModelPart.GetMesh(0).Elements().size();
        const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE];
        std::vector<std::size_t> ConditionNodeIds(Dim);
        if (mNumNode == 8)
            ConditionNodeIds.resize(Dim + 1);

        int last_condition_id = mMechanicalLastCondition + mThermalLastCondition;
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();

            if (Dim == 2)
            {
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    array_1d<double, 3> central_position = it_thermal->GetGeometry().Center();

                    // Elements
                    if ((it_thermal)->IsNot(ACTIVE))
                    {
                        if (central_position(direction) <= mHighestBlockHeight + mLiftHeight)
                        {
                            for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                            {
                                const unsigned int number_of_points = (*it_thermal).GetGeometry().GenerateEdges()[i_edge].PointsNumber();
                                bool active_edge = true;

                                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                                {
                                    if ((*it_thermal).GetGeometry().GenerateEdges()[i_edge][i_node].IsNot(ACTIVE))
                                    {
                                        active_edge = false;
                                        break;
                                    }
                                }
                                if (active_edge)
                                {
                                    for (unsigned int m = 0; m < number_of_points; ++m)
                                    {
                                        ConditionNodeIds[m] = (*it_thermal).GetGeometry().GenerateEdges()[i_edge][m].Id();
                                    }
                                    this->DeactiveFaceHeatFluxStep(ConditionNodeIds);

                                    mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id + 1, 0);
                                    last_condition_id++;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    array_1d<double, 3> central_position = it_thermal->GetGeometry().Center();

                    // Elements
                    if ((it_thermal)->IsNot(ACTIVE))
                    {
                        if (central_position(direction) <= mHighestBlockHeight + mLiftHeight)
                        {
                            for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                            {
                                const unsigned int number_of_points = (*it_thermal).GetGeometry().GenerateFaces()[i_face].PointsNumber();
                                bool active_face = true;

                                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                                {
                                    if ((*it_thermal).GetGeometry().GenerateFaces()[i_face][i_node].IsNot(ACTIVE))
                                    {
                                        active_face = false;
                                        break;
                                    }
                                }
                                if (active_face)
                                {
                                    for (unsigned int m = 0; m < number_of_points; ++m)
                                    {
                                        ConditionNodeIds[m] = (*it_thermal).GetGeometry().GenerateFaces()[i_face][m].Id();
                                    }
                                    this->DeactiveFaceHeatFluxStep(ConditionNodeIds);

                                    mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id + 1, 0);
                                    last_condition_id++;
                                }
                            }
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SearchingFluxes()
    {
        KRATOS_TRY;

        const int nelements = mrMechanicalModelPart.GetMesh(0).Elements().size();
        const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE];
        std::vector<std::size_t> ConditionNodeIds(Dim);
        if (mNumNode == 8)
            ConditionNodeIds.resize(Dim + 1);

        int last_condition_id = mMechanicalLastCondition + mThermalLastCondition;
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();

            if (Dim == 2)
            {
                // Searching for thermal boundary conditions Edges
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    array_1d<double, 3> central_position = it_thermal->GetGeometry().Center();

                    // Elements
                    if ((it_thermal)->IsNot(ACTIVE))
                    {
                        if (central_position(direction) <= mHighestBlockHeight + mLiftHeight)
                        {
                            for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                            {
                                const unsigned int number_of_points = (*it_thermal).GetGeometry().GenerateEdges()[i_edge].PointsNumber();
                                bool active_edge = true;

                                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                                {
                                    if ((*it_thermal).GetGeometry().GenerateEdges()[i_edge][i_node].IsNot(ACTIVE))
                                    {
                                        active_edge = false;
                                        break;
                                    }
                                }
                                if (active_edge)
                                {
                                    for (unsigned int m = 0; m < number_of_points; ++m)
                                    {
                                        ConditionNodeIds[m] = (*it_thermal).GetGeometry().GenerateEdges()[i_edge][m].Id();
                                    }
                                    this->ActiveFaceHeatFluxStep(ConditionNodeIds);

                                    mrThermalModelPart.CreateNewCondition("FluxCondition2D2N", last_condition_id + 1, ConditionNodeIds, 0);
                                    last_condition_id++;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // Searching for thermal boundary conditions
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    array_1d<double, 3> central_position = it_thermal->GetGeometry().Center();

                    // Elements
                    if ((it_thermal)->IsNot(ACTIVE))
                    {
                        if (central_position(direction) <= mHighestBlockHeight + mLiftHeight)
                        {
                            for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                            {
                                const unsigned int number_of_points = (*it_thermal).GetGeometry().GenerateFaces()[i_face].PointsNumber();
                                bool active_face = true;

                                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                                {
                                    if ((*it_thermal).GetGeometry().GenerateFaces()[i_face][i_node].IsNot(ACTIVE))
                                    {
                                        active_face = false;
                                        break;
                                    }
                                }
                                if (active_face)
                                {
                                    for (unsigned int m = 0; m < number_of_points; ++m)
                                    {
                                        ConditionNodeIds[m] = (*it_thermal).GetGeometry().GenerateFaces()[i_face][m].Id();
                                    }
                                    this->ActiveFaceHeatFluxStep(ConditionNodeIds);

                                    if (number_of_points == 3)
                                    {
                                        mrThermalModelPart.CreateNewCondition("FluxCondition3D3N", last_condition_id + 1, ConditionNodeIds, 0);
                                        last_condition_id++;
                                    }
                                    else
                                    {
                                        mrThermalModelPart.CreateNewCondition("FluxCondition3D4N", last_condition_id + 1, ConditionNodeIds, 0);
                                        last_condition_id++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveHeatFluxNoorzai(Parameters &NoorzaiParameters)
    {
        KRATOS_TRY;

        // Getting Noorzai Values
        const double density = NoorzaiParameters["density"].GetDouble();
        const double specific_heat = NoorzaiParameters["specific_heat"].GetDouble();
        const double alpha = NoorzaiParameters["alpha"].GetDouble();
        const double t_max = NoorzaiParameters["t_max"].GetDouble();
        const double time = mrThermalModelPart.GetProcessInfo()[TIME];
        const double delta_time = mrThermalModelPart.GetProcessInfo()[DELTA_TIME];

        block_for_each(mrThermalModelPart.Nodes(), [&](auto& rNode){
            double current_activation_time = time - (rNode.FastGetSolutionStepValue(TIME_ACTIVATION));
            if (current_activation_time >= 0.0 && (rNode.Is(SOLID) == false))
            {
                // Computing the value of heat flux according the time
                double value = density * specific_heat * alpha * t_max * (exp(-alpha * (current_activation_time + 0.5 * delta_time)));
                rNode.FastGetSolutionStepValue(HEAT_FLUX) = value;
            }
        });
        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveHeatFluxAzenha(Parameters &AzenhaParameters)
    {
        KRATOS_TRY;

        if (mAging == false)
        {
            // Getting Azenha Values
            const double activation_energy = AzenhaParameters["activation_energy"].GetDouble();
            const double gas_constant = AzenhaParameters["gas_constant"].GetDouble();
            const double constant_rate = AzenhaParameters["constant_rate"].GetDouble();
            const double q_total = AzenhaParameters["q_total"].GetDouble();
            const double a_coef = AzenhaParameters["A"].GetDouble();
            const double b_coef = AzenhaParameters["B"].GetDouble();
            const double c_coef = AzenhaParameters["C"].GetDouble();
            const double d_coef = AzenhaParameters["D"].GetDouble();

            // Temporal variables
            const double time = mrThermalModelPart.GetProcessInfo()[TIME];
            const double delta_time = mrThermalModelPart.GetProcessInfo()[DELTA_TIME];

            block_for_each(mrThermalModelPart.Nodes(), [&](auto& rNode){
                double current_activation_time = time - (rNode.FastGetSolutionStepValue(TIME_ACTIVATION));
                if (current_activation_time >= 0.0 && (rNode.Is(SOLID) == false))
                {
                    // Computing the current alpha according las step.
                    double current_alpha = ((rNode.FastGetSolutionStepValue(HEAT_FLUX, 1)) / q_total) * delta_time + (rNode.FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                    double f_alpha = a_coef * (pow(current_alpha, 2)) * exp(-b_coef * pow(current_alpha, 3)) + c_coef * current_alpha * exp(-d_coef * current_alpha);

                    // This is neccesary for stopping the addition to the system once the process finish.
                    if (current_alpha >= 1.0)
                    {
                        f_alpha = 0.0;
                        current_alpha = 1.0;
                    }

                    // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                    const double temp_current = rNode.FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                    const double heat_flux = constant_rate * f_alpha * exp((-activation_energy) / (gas_constant * temp_current));

                    // Updating values according the computations
                    rNode.FastGetSolutionStepValue(HEAT_FLUX) = heat_flux;
                    rNode.FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                }
            });
        }
        else
        {
            this->ActiveHeatFluxAzenhaAging(AzenhaParameters);
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    ModelPart &mrMechanicalModelPart;
    ModelPart &mrThermalModelPart;
    int mNumNode;
    std::string mGravityDirection;
    std::string mMechanicalSoilPart;
    std::string mThermalSoilPart;
    std::string mMechanicalExistingPart;
    std::string mThermalExistingPart;
    std::string mSourceType;
    bool mActivateSoilPart;
    bool mActivateExistingPart;
    double mReferenceCoordinate;
    double mHighestBlockHeight;
    double mLiftHeight;
    double mH0;
    double mTimeUnitConverter;
    double mAlphaInitial;
    bool mAging;
    double mYoungInf;
    unsigned int mMechanicalLastCondition;
    unsigned int mThermalLastCondition;
    TableType &mrTableAmbientTemp;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveHeatFluxAzenhaAging(Parameters &AzenhaParameters)
    {
        KRATOS_TRY;

        // Getting Azenha Values
        double activation_energy = AzenhaParameters["activation_energy"].GetDouble();
        double gas_constant = AzenhaParameters["gas_constant"].GetDouble();
        double constant_rate = AzenhaParameters["constant_rate"].GetDouble();
        double q_total = AzenhaParameters["q_total"].GetDouble();
        double a_coef = AzenhaParameters["A"].GetDouble();
        double b_coef = AzenhaParameters["B"].GetDouble();
        double c_coef = AzenhaParameters["C"].GetDouble();
        double d_coef = AzenhaParameters["D"].GetDouble();

        // Tempotal variables
        double time = mrThermalModelPart.GetProcessInfo()[TIME];
        double delta_time = mrThermalModelPart.GetProcessInfo()[DELTA_TIME];

        block_for_each(mrThermalModelPart.Nodes(), [&](auto& rNode){
            double current_activation_time = time - (rNode.FastGetSolutionStepValue(TIME_ACTIVATION));
            if (current_activation_time >= 0.0 && (rNode.Is(SOLID) == false))
            {
                // Computing the current alpha according las step.
                double current_alpha = ((rNode.FastGetSolutionStepValue(HEAT_FLUX, 1)) / q_total) * delta_time + (rNode.FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                double f_alpha = a_coef * (pow(current_alpha, 2)) * exp(-b_coef * pow(current_alpha, 3)) + c_coef * current_alpha * exp(-d_coef * current_alpha);

                // This is neccesary for stopping the addition to the system once the process finish.
                if (current_alpha >= 1.0)
                {
                    f_alpha = 0.0;
                    current_alpha = 1.0;
                }

                // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                const double temp_current = rNode.FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                const double heat_flux = constant_rate * f_alpha * exp((-activation_energy) / (gas_constant * temp_current));

                // Updating values according the computations
                rNode.FastGetSolutionStepValue(HEAT_FLUX) = heat_flux;
                rNode.FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                rNode.FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = sqrt(current_alpha) * mYoungInf;
            }
        });

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveFaceHeatFluxStep(std::vector<IndexType> ConditionNodeIds)
    {
        KRATOS_TRY;

        const unsigned int size = ConditionNodeIds.size();
        const unsigned int nnodes = mrThermalModelPart.GetMesh(0).Nodes().size();
        ModelPart::NodesContainerType::iterator it_begin_thermal = mrThermalModelPart.GetMesh(0).NodesBegin();

        double time = mrThermalModelPart.GetProcessInfo()[TIME];
        time = time / mTimeUnitConverter;
        double ambient_temp = mrTableAmbientTemp.GetValue(time);

        if (size != 0)
        {
            for (unsigned int j = 0; j < size; ++j)
            {
                for (unsigned int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it_thermal = it_begin_thermal + i;

                    if (it_thermal->Id() == ConditionNodeIds[j])
                    {
                        const double temp_current = it_thermal->FastGetSolutionStepValue(TEMPERATURE);
                        const double heat_flux = mH0 * (ambient_temp - temp_current);
                        it_thermal->FastGetSolutionStepValue(FACE_HEAT_FLUX) = heat_flux;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void DeactiveFaceHeatFluxStep(std::vector<IndexType> ConditionNodeIds)
    {
        KRATOS_TRY;

        const unsigned int size = ConditionNodeIds.size();
        const unsigned int nnodes = mrThermalModelPart.GetMesh(0).Nodes().size();
        ModelPart::NodesContainerType::iterator it_begin_thermal = mrThermalModelPart.GetMesh(0).NodesBegin();

        if (size != 0)
        {
            for (unsigned int j = 0; j < size; ++j)
            {
                for (unsigned int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it_thermal = it_begin_thermal + i;

                    if (it_thermal->Id() == ConditionNodeIds[j])
                    {
                        //Setting to 0 the fluxes since in the next step
                        it_thermal->FastGetSolutionStepValue(FACE_HEAT_FLUX) = 0.0;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */
