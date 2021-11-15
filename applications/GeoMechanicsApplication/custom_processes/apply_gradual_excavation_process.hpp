// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia,
//                   Aron Noordam,
//                   Vahid Galavi
//
//


// System includes
#include <cmath>
#include <iostream>
#include<string>
#include<limits>

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

#if !defined(KRATOS_GEO_APPLY_GRADUAL_EXCAVATION_PROCESS )
#define  KRATOS_GEO_APPLY_GRADUAL_EXCAVATION_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_processes/apply_excavation_process.hpp"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyGradualExcavationProcess : public ApplyExcavationProcess
{

  public:


    typedef std::size_t IndexType;

    KRATOS_CLASS_POINTER_DEFINITION(ApplyGradualExcavationProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyGradualExcavationProcess(ModelPart&  model_part,
                           Parameters& rParameters) : ApplyExcavationProcess(model_part, rParameters)
    {
        KRATOS_TRY
        //unsigned int TableId = rParameters["table"].GetInt();
        //mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        mGravityDirection = rParameters["gravity_direction"].GetInt();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyGradualExcavationProcess() override{}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        double Time = mr_model_part.GetProcessInfo()[TIME] / mTimeUnitConverter;

        const double start_time = mr_model_part.GetProcessInfo()[START_TIME] / mTimeUnitConverter;
        const double end_time = mr_model_part.GetProcessInfo()[END_TIME] / mTimeUnitConverter;

		Time = Time + start_time;
        std::cout << "Time " << Time << '\n';
        std::cout << "start_time " << start_time << '\n';
        std::cout << "end_time " << end_time << '\n';
        const int nelements = mr_model_part.GetMesh(0).Elements().size();
        const int nnodes = mr_model_part.GetMesh(0).Nodes().size();
        // mModelLastCondition = mr_model_part.GetMesh(0).Conditions().size();

        array_1d<double, 3> Coordinates;

        // calculate the maximum and minimum coordinate and calculate the total height of the model part
        double max_vertical_coordinate = -9999999;
        double min_vertical_coordinate = 9999999;

        ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            noalias(Coordinates) = it->Coordinates();
            if (Coordinates[mGravityDirection] > max_vertical_coordinate)
            {
                max_vertical_coordinate = Coordinates[mGravityDirection];
            }
            if (Coordinates[mGravityDirection] < min_vertical_coordinate)
            {
                min_vertical_coordinate = Coordinates[mGravityDirection];
            }
        }

        const double total_height = max_vertical_coordinate - min_vertical_coordinate;
        double deltaH = (Time - start_time )/(end_time - start_time) * total_height;

		std::cout << "min_vertical_coordinate " << min_vertical_coordinate << '\n';
        std::cout << "total_height " << total_height << '\n';
        std::cout << "deltaH " << deltaH << '\n';

        double activate_coordinate = min_vertical_coordinate + deltaH;

        if (nelements != 0 && mDeactivateSoilPart == true)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();

            #pragma omp parallel for
            for (int k = 0; k < nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                it->Set(ACTIVE, false);
            }

            // Same nodes for both computing model part
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
            #pragma omp parallel for
            for (int i = 0; i < nnodes; ++i)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Set(ACTIVE, false);
                it->Set(SOLID, false);
            }    
        }
        
        // Activation of the existing parts, either the soil or the already built dam ( User must specify each part through the interface)
        if (mDeactivateSoilPart == false)
        {
            if (nelements != 0)
            {
                ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();

                #pragma omp parallel for
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it_elem = el_begin + k;

                    const unsigned int number_of_points = (*it_elem).GetGeometry().PointsNumber();

                    double centre_vertical_coordinate = 0.0;
                    for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                    {
                        centre_vertical_coordinate = ((*it_elem).GetGeometry()[i_node].Coordinates()[mGravityDirection]) + centre_vertical_coordinate;
                    }
                    centre_vertical_coordinate = centre_vertical_coordinate / number_of_points;

                    std::cout << "centre_vertical_coordinate " << centre_vertical_coordinate << '\n';
                    std::cout << "activate_coordinate " << activate_coordinate << '\n';

                    if (centre_vertical_coordinate < activate_coordinate) 
                    { 
                        it_elem->Set(ACTIVE, true);
                    }
                }

                // Same nodes for both computing model part
                ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
                #pragma omp parallel for
                for (int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    if (it->Coordinates()[mGravityDirection] < activate_coordinate)
                    {
                        it->Set(ACTIVE, true);
                        it->Set(SOLID, true);
                    }
                }
            }
        }

        // Conditions
        const int nconditions = mr_model_part.GetMesh(0).Conditions().size();
        if (nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();

            for (int k = 0; k < nconditions; ++k)
            {
                ModelPart::ConditionsContainerType::iterator it_cond = cond_begin + k;
                const unsigned int number_of_points = (*it_cond).GetGeometry().PointsNumber();
                bool active_condition = true;

                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                {
                    if ((*it_cond).GetGeometry()[i_node].IsNot(ACTIVE))
                    {
                        active_condition = false;
                        break;
                    }
                }
                if (active_condition) it_cond->Set(ACTIVE, true);
                else it_cond->Set(ACTIVE, false);
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables
      //TableType::Pointer mpTable;
      double mTimeUnitConverter;
      unsigned int mGravityDirection;

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */