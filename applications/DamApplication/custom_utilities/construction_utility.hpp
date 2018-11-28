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
        mReferenceCoordinate = rParameters["reservoir_bottom_coordinate_in_gravity_direction"].GetDouble();
        mHeight = rParameters["height_dam"].GetDouble();
        mPhases = rParameters["number_of_phases"].GetInt();
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
        const int nnodes = mrMechanicalModelPart.GetMesh(0).Nodes().size();

        mMechanicalLastCondition = mrMechanicalModelPart.GetMesh(0).Conditions().size();
        mThermalLastCondition = mrThermalModelPart.GetMesh(0).Conditions().size();

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.ElementsBegin();
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();
            mNumNode = el_begin->GetGeometry().PointsNumber();

            #pragma omp parallel for
            for (int k = 0; k < nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                it->Set(ACTIVE, false);
                it_thermal->Set(ACTIVE, false);
            }

            // Same nodes for both computing model part
            ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();
            #pragma omp parallel for
            for (int i = 0; i < nnodes; ++i)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Set(ACTIVE, false);
                it->Set(SOLID, false);
            }
        }

        // Activation of the existing parts, either the soil or the already built dam ( User must specify each part through the interface)
        if (mActivateSoilPart == true)
        {
            const int soil_nelements = mrMechanicalModelPart.GetSubModelPart(mMechanicalSoilPart).Elements().size();
            const int soil_nnodes = mrMechanicalModelPart.GetSubModelPart(mMechanicalSoilPart).Nodes().size();

            if (soil_nelements != 0)
            {
                ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.GetSubModelPart(mMechanicalSoilPart).ElementsBegin();
                ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(mThermalSoilPart).ElementsBegin();
                mNumNode = el_begin->GetGeometry().PointsNumber();

                #pragma omp parallel for
                for (int k = 0; k < soil_nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    it->Set(ACTIVE, true);
                    it_thermal->Set(ACTIVE, true);
                }

                // Same nodes for both computing model part
                ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.GetSubModelPart(mThermalSoilPart).NodesBegin();
                #pragma omp parallel for
                for (int i = 0; i < soil_nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->Set(ACTIVE, true);
                    it->Set(SOLID, true);
                }
            }
        }

        if (mActivateExistingPart == true)
        {
            const int existing_nelements = mrMechanicalModelPart.GetSubModelPart(mMechanicalExistingPart).Elements().size();
            const int existing_nnodes = mrMechanicalModelPart.GetSubModelPart(mMechanicalExistingPart).Nodes().size();

            if (existing_nelements != 0)
            {
                ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.GetSubModelPart(mMechanicalExistingPart).ElementsBegin();
                ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(mThermalExistingPart).ElementsBegin();
                mNumNode = el_begin->GetGeometry().PointsNumber();

                #pragma omp parallel for
                for (int k = 0; k < existing_nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    it->Set(ACTIVE, true);
                    it_thermal->Set(ACTIVE, true);
                }

                // Same nodes for both computing model part
                ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.GetSubModelPart(mThermalExistingPart).NodesBegin();
                #pragma omp parallel for
                for (int i = 0; i < existing_nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->Set(ACTIVE, true);
                    it->Set(SOLID, true);
                }
            }
        }

        // Assign Alpha Initial in case of using Azenha Formulation
        if (mSourceType == "NonAdiabatic")
        {
            if (mAging == false)
            {
                ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();
            #pragma omp parallel for
                for (int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = mAlphaInitial;
                }
            }
            else
            {
                ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();
            #pragma omp parallel for
                for (int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = mAlphaInitial;
                    it->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = sqrt(mAlphaInitial) * mYoungInf;
                }
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
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).ElementsBegin();

            double current_height = mReferenceCoordinate + (mHeight / mPhases) * (phase);
            double previous_height = mReferenceCoordinate + (mHeight / mPhases) * (phase - 1);

            #pragma omp parallel for
            for (int k = 0; k < nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                array_1d<double, 3> central_position = it_thermal->GetGeometry().Center();

                if ((central_position(direction) >= previous_height) && (central_position(direction) <= current_height))
                {

                    const unsigned int number_of_points = it_thermal->GetGeometry().PointsNumber();
                    for (unsigned int i = 0; i < number_of_points; ++i)
                    {
                        if (it_thermal->GetGeometry()[i].FastGetSolutionStepValue(TIME_ACTIVATION)==0)
                        {
                            it_thermal->GetGeometry()[i].FastGetSolutionStepValue(TIME_ACTIVATION) = time_activation * mTimeUnitConverter;
                            it_thermal->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE) = initial_temperature;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(std::string ThermalSubModelPartName, std::string MechanicalSubModelPartName, int current_number_of_phase)
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

        // Getting the value of the table and computing the current height
        double current_height = mReferenceCoordinate + (mHeight / mPhases) * current_number_of_phase;

        if (nelements != 0)
        {
            // ELEMENTS
            ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.GetSubModelPart(MechanicalSubModelPartName).ElementsBegin();
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).ElementsBegin();

            #pragma omp parallel for
            for (int k = 0; k < nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                array_1d<double, 3> central_position = it->GetGeometry().Center();

                if ((central_position(direction) >= mReferenceCoordinate) && (central_position(direction) <= current_height))
                {
                    it->Set(ACTIVE, true);
                    it_thermal->Set(ACTIVE, true);

                    const unsigned int number_of_points = it_thermal->GetGeometry().PointsNumber();
                    for (unsigned int i = 0; i < number_of_points; i++)
                    {
                        it->GetGeometry()[i].Set(ACTIVE, true);
                    }
                }
            }
        }

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

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();

            if (Dim == 2)
            {
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    // Elements
                    if ((it_thermal)->Is(ACTIVE) == false)
                    {
                        for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Edges()[i_edge].PointsNumber();
                            bool active_edge = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if ((*it_thermal).GetGeometry().Edges()[i_edge][i_node].Is(ACTIVE) == true)
                                {
                                    active_edge = false;
                                    break;
                                }
                            }
                            if (active_edge)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Edges()[i_edge][m].Id();
                                }
                                this->DeactiveFaceHeatFluxStep(ConditionNodeIds);

                                mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id + 1, 0);
                                last_condition_id++;
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
                    // Elements
                    if ((it_thermal)->Is(ACTIVE) == false)
                    {
                        for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                            bool active_face = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if ((*it_thermal).GetGeometry().Faces()[i_face][i_node].Is(ACTIVE) == true)
                                {
                                    active_face = false;
                                    break;
                                }
                            }
                            if (active_face)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
                                }
                                this->DeactiveFaceHeatFluxStep(ConditionNodeIds);

                                mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id + 1, 0);
                                last_condition_id++;
                            }
                        }
                    }
                    if ((it_thermal)->Is(ACTIVE) == true)
                    {
                        for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                            bool active_face = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if (!(*it_thermal).GetGeometry().Faces()[i_face][i_node].FastGetSolutionStepValue(DAM_SURFACE_NODE))
                                {
                                    active_face = false;
                                    break;
                                }
                            }
                            if (active_face)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
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

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();

            if (Dim == 2)
            {
                // Searching for thermal boundary conditions Edges
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                    // Elements
                    if ((it_thermal)->Is(ACTIVE) == false)
                    {
                        for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Edges()[i_edge].PointsNumber();
                            bool active_edge = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if ((*it_thermal).GetGeometry().Edges()[i_edge][i_node].Is(ACTIVE) == false)
                                {
                                    active_edge = false;
                                    break;
                                }
                            }
                            if (active_edge)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Edges()[i_edge][m].Id();
                                }
                                this->ActiveFaceHeatFluxStep(ConditionNodeIds);

                                mrThermalModelPart.CreateNewCondition("FluxCondition2D2N", last_condition_id + 1, ConditionNodeIds, 0);
                                last_condition_id++;
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
                    // Elements
                    if ((it_thermal)->Is(ACTIVE) == false)
                    {
                        for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                            bool active_face = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if ((*it_thermal).GetGeometry().Faces()[i_face][i_node].Is(ACTIVE) == false)
                                {
                                    active_face = false;
                                    break;
                                }
                            }
                            if (active_face)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
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
                    if ((it_thermal)->Is(ACTIVE) == true)
                    {
                        for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                        {
                            const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                            bool surface_condition = true;

                            for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                            {
                                if (!(*it_thermal).GetGeometry().Faces()[i_face][i_node].FastGetSolutionStepValue(DAM_SURFACE_NODE))
                                {
                                    surface_condition = false;
                                    break;
                                }
                            }
                            if (surface_condition)
                            {
                                for (unsigned int m = 0; m < number_of_points; ++m)
                                {
                                    ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
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

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveHeatFluxNoorzai(Parameters &NoorzaiParameters)
    {
        KRATOS_TRY;

        const int nnodes = mrThermalModelPart.Nodes().size();

        // Getting Noorzai Values
        double density = NoorzaiParameters["density"].GetDouble();
        double specific_heat = NoorzaiParameters["specific_heat"].GetDouble();
        double alpha = NoorzaiParameters["alpha"].GetDouble();
        double t_max = NoorzaiParameters["t_max"].GetDouble();
        double time = mrThermalModelPart.GetProcessInfo()[TIME];

        ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            double current_activation_time = time - (it->FastGetSolutionStepValue(TIME_ACTIVATION));
            if (current_activation_time >= 0.0 && (it->Is(SOLID) == false))
            {
                // Computing the value of heat flux according the time
                double value = density * specific_heat * alpha * t_max * (exp(-alpha * current_activation_time));
                it->FastGetSolutionStepValue(HEAT_FLUX) = value;
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ActiveHeatFluxAzenha(Parameters &AzenhaParameters)
    {
        KRATOS_TRY;

        if (mAging == false)
        {
            const int nnodes = mrThermalModelPart.Nodes().size();

            // Getting Azenha Values
            double activation_energy = AzenhaParameters["activation_energy"].GetDouble();
            double gas_constant = AzenhaParameters["gas_constant"].GetDouble();
            double constant_rate = AzenhaParameters["constant_rate"].GetDouble();
            double q_total = AzenhaParameters["q_total"].GetDouble();
            double a_coef = AzenhaParameters["A"].GetDouble();
            double b_coef = AzenhaParameters["B"].GetDouble();
            double c_coef = AzenhaParameters["C"].GetDouble();
            double d_coef = AzenhaParameters["D"].GetDouble();

            // Temporal variables
            double time = mrThermalModelPart.GetProcessInfo()[TIME];
            double delta_time = mrThermalModelPart.GetProcessInfo()[DELTA_TIME];

            ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();
            #pragma omp parallel for
            for (int i = 0; i < nnodes; ++i)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                double current_activation_time = time - (it->FastGetSolutionStepValue(TIME_ACTIVATION));
                if (current_activation_time >= 0.0 && (it->Is(SOLID) == false))
                {
                    // Computing the current alpha according las step.
                    double current_alpha = ((it->FastGetSolutionStepValue(HEAT_FLUX, 1)) / q_total) * delta_time + (it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                    double f_alpha = a_coef * (pow(current_alpha, 2)) * exp(-b_coef * pow(current_alpha, 3)) + c_coef * current_alpha * exp(-d_coef * current_alpha);

                    // This is neccesary for stopping the addition to the system once the process finish.
                    if (current_alpha >= 1.0)
                    {
                        f_alpha = 0.0;
                        current_alpha = 1.0;
                    }

                    // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                    const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                    const double heat_flux = constant_rate * f_alpha * exp((-activation_energy) / (gas_constant * temp_current));

                    // Updating values according the computations
                    it->FastGetSolutionStepValue(HEAT_FLUX) = heat_flux;
                    it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                }
            }
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
    double mHeight;
    int mPhases;
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

        const int nnodes = mrThermalModelPart.Nodes().size();

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

        ModelPart::NodesContainerType::iterator it_begin = mrThermalModelPart.NodesBegin();
        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            double current_activation_time = time - (it->FastGetSolutionStepValue(TIME_ACTIVATION));
            if (current_activation_time >= 0.0 && (it->Is(SOLID) == false))
            {
                // Computing the current alpha according las step.
                double current_alpha = ((it->FastGetSolutionStepValue(HEAT_FLUX, 1)) / q_total) * delta_time + (it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                double f_alpha = a_coef * (pow(current_alpha, 2)) * exp(-b_coef * pow(current_alpha, 3)) + c_coef * current_alpha * exp(-d_coef * current_alpha);

                // This is neccesary for stopping the addition to the system once the process finish.
                if (current_alpha >= 1.0)
                {
                    f_alpha = 0.0;
                    current_alpha = 1.0;
                }

                // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                const double heat_flux = constant_rate * f_alpha * exp((-activation_energy) / (gas_constant * temp_current));

                // Updating values according the computations
                it->FastGetSolutionStepValue(HEAT_FLUX) = heat_flux;
                it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                it->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = sqrt(current_alpha) * mYoungInf;
            }
        }

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
