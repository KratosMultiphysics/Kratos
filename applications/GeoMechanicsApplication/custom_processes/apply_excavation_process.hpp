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

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

#if !defined(KRATOS_GEO_APPLY_EXCAVATION_PROCESS )
#define  KRATOS_GEO_APPLY_EXCAVATION_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyExcavationProcess : public Process
{

  public:


    typedef std::size_t IndexType;
    typedef Table<double, double> TableType;

    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyExcavationProcess(ModelPart&  model_part,
                           Parameters& rParameters) : Process(Flags()), mr_model_part(model_part)
    {
        KRATOS_TRY
        mDeactivateSoilPart =  rParameters["deactivate_soil_part"].GetBool();
        mModelPartName      =  rParameters["model_part_name"].GetString();
        

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyExcavationProcess() override{}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const int nelements = mr_model_part.GetMesh(0).Elements().size();
        const int nnodes = mr_model_part.GetMesh(0).Nodes().size();

        if (nelements != 0)
        {
            if (mDeactivateSoilPart == true)
            {
                // Deactivation of the existing parts:
                // ( User must specify each part through the interface)
                ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
                #pragma omp parallel for
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    it->Set(ACTIVE, false);
                    it->ResetConstitutiveLaw();
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
            else
            {
                // Activation of the existing parts:
                // ( User must specify each part through the interface)
                ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
                #pragma omp parallel for
                for (int k = 0; k < nelements; ++k)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    it->Set(ACTIVE, true);
                }

                // Same nodes for both computing model part
                ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();
                #pragma omp parallel for
                for (int i = 0; i < nnodes; ++i)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;
                    it->Set(ACTIVE, true);
                    it->Set(SOLID, true);
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

                // VG: there is a problem in this part
                /*
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
                */

                if (mDeactivateSoilPart == true) it_cond->Set(ACTIVE, false);
                else it_cond->Set(ACTIVE, true);

            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    ModelPart& mr_model_part;
    bool mDeactivateSoilPart;
    std::string mModelPartName;

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */