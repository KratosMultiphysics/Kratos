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

#if !defined(KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS)
#define KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamUpliftCircularConditionLoadProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamUpliftCircularConditionLoadProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamUpliftCircularConditionLoadProcess(ModelPart &rModelPart,
                                          Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "Modify"                                                : true,
                "Gravity_Direction"                                     : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Upstream_Coordinate_first_bracket"                     : [0.0,0.0,0.0],
                "Downstream_Coordinate_first_bracket"                   : [0.0,0.0,0.0],
                "Focus"                                                 : [0.0,0.0,0.0],
                "Spe_weight"                                            : 10000,
                "Water_level"                                           : 0.0,
                "Drains"                                                : false,
                "Height_drain"                                          : 0.0,
                "Distance"                                              : 0.0,
                "Effectiveness"                                         : 0.0,
                "table"                                                 : 0
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Focus"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mSpecific = rParameters["Spe_weight"].GetDouble();

        // Getting the values of the coordinates at first bracket (reference value)
        mUpstream.resize(3, false);
        mUpstream[0] = rParameters["Upstream_Coordinate_first_bracket"][0].GetDouble();
        mUpstream[1] = rParameters["Upstream_Coordinate_first_bracket"][1].GetDouble();
        mUpstream[2] = rParameters["Upstream_Coordinate_first_bracket"][2].GetDouble();

        mDownstream.resize(3, false);
        mDownstream[0] = rParameters["Downstream_Coordinate_first_bracket"][0].GetDouble();
        mDownstream[1] = rParameters["Downstream_Coordinate_first_bracket"][1].GetDouble();
        mDownstream[2] = rParameters["Downstream_Coordinate_first_bracket"][2].GetDouble();

        // Getting the coordinates of the focus (reference value)
        mFocus.resize(3, false);
        mFocus[0] = rParameters["Focus"][0].GetDouble();
        mFocus[1] = rParameters["Focus"][1].GetDouble();
        mFocus[2] = rParameters["Focus"][2].GetDouble();

        // Drains
        mDrain = rParameters["Drains"].GetBool();
        mHeightDrain = rParameters["Height_drain"].GetDouble();
        mDistanceDrain = rParameters["Distance"].GetDouble();
        mEffectivenessDrain = rParameters["Effectiveness"].GetDouble();
        mWaterLevel = rParameters["Water_level"].GetInt();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();

        if (mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamUpliftCircularConditionLoadProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {

        KRATOS_TRY;

        //Defining necessary variables
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        array_1d<double, 3> auxiliar_vector;

        // Gravity direction for computing the hydrostatic pressure
        int direction;
        int radius_comp_1;
        int radius_comp_2;

        if (mGravityDirection == "X")
        {
            direction = 0;
            radius_comp_1 = 1;
            radius_comp_2 = 2;
        }
        else if (mGravityDirection == "Y")
        {
            direction = 1;
            radius_comp_1 = 0;
            radius_comp_2 = 2;
        }
        else
        {
            direction = 2;
            radius_comp_1 = 0;
            radius_comp_2 = 1;
        }

        // Computation of the angle in radians for each bracket
        //double tetha = (mangle*2*pi)/(mnum_brackets*360);

        // Computation of radius for Upstream and Downstream
        double up_radius = norm_2(mFocus - mUpstream);
        double down_radius = norm_2(mFocus - mDownstream);
        double width_dam = up_radius - down_radius;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            double ref_coord = mReferenceCoordinate + mWaterLevel;

            if (mDrain == true)
            {
                double coefficient_effectiveness = 1.0 - mEffectivenessDrain;
                double aux_drain = coefficient_effectiveness * (mWaterLevel - mHeightDrain) * ((width_dam - mDistanceDrain) / width_dam) + mHeightDrain;

#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    auxiliar_vector.resize(3, false);
                    const array_1d<double,3>& r_coordinates = it->Coordinates();
                    noalias(auxiliar_vector) = mFocus - r_coordinates;

                    //// Computing the new coordinates
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1] * auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2] * auxiliar_vector[radius_comp_2]);

                    //// We compute the first part of the uplift law
                    mUpliftPressure = mSpecific * ((ref_coord - aux_drain) - (r_coordinates[direction])) * (1.0 - ((1.0 / mDistanceDrain) * (fabs(current_radius - up_radius)))) + (mSpecific * aux_drain);

                    //// If uplift pressure is greater than the limit we compute the second part and we update the value
                    if (mUpliftPressure <= mSpecific * aux_drain)
                    {
                        mUpliftPressure = (mSpecific * ((mReferenceCoordinate + aux_drain) - (r_coordinates[direction]))) * (1.0 - ((1.0 / (width_dam - mDistanceDrain)) * (fabs(current_radius - (up_radius - mDistanceDrain)))));
                    }

                    if (mUpliftPressure < 0.0)
                    {
                        it->FastGetSolutionStepValue(var) = 0.0;
                    }
                    else
                    {
                        it->FastGetSolutionStepValue(var) = mUpliftPressure;
                    }
                }
            }
            else
            {
#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    auxiliar_vector.resize(3, false);
                    const array_1d<double,3>& r_coordinates = it->Coordinates();
                    noalias(auxiliar_vector) = mFocus - r_coordinates;

                    // Computing the current distance to the focus.
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1] * auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2] * auxiliar_vector[radius_comp_2]);

                    mUpliftPressure = mSpecific * (ref_coord - (r_coordinates[direction])) * (1.0 - (1.0 / width_dam) * (fabs(current_radius - up_radius)));

                    if (mUpliftPressure < 0.0)
                    {
                        it->FastGetSolutionStepValue(var) = 0.0;
                    }
                    else
                    {
                        it->FastGetSolutionStepValue(var) = mUpliftPressure;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        //Defining necessary variables
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        array_1d<double, 3> auxiliar_vector;

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTable->GetValue(time);
        }

        // Gravity direction for computing the hydrostatic pressure
        int direction;
        int radius_comp_1;
        int radius_comp_2;

        if (mGravityDirection == "X")
        {
            direction = 0;
            radius_comp_1 = 1;
            radius_comp_2 = 2;
        }
        else if (mGravityDirection == "Y")
        {
            direction = 1;
            radius_comp_1 = 0;
            radius_comp_2 = 2;
        }
        else
        {
            direction = 2;
            radius_comp_1 = 0;
            radius_comp_2 = 1;
        }

        // Computation of the angle in radians for each bracket
        //double tetha = (mangle*2*pi)/(mnum_brackets*360);

        // Computation of radius for Upstream and Downstream
        double up_radius = norm_2(mFocus - mUpstream);
        double down_radius = norm_2(mFocus - mDownstream);
        double width_dam = up_radius - down_radius;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            double ref_coord = mReferenceCoordinate + mWaterLevel;

            if (mDrain == true)
            {
                double coefficient_effectiveness = 1.0 - mEffectivenessDrain;
                double aux_drain = coefficient_effectiveness * (mWaterLevel - mHeightDrain) * ((width_dam - mDistanceDrain) / width_dam) + mHeightDrain;

#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    auxiliar_vector.resize(3, false);
                    const array_1d<double,3>& r_coordinates = it->Coordinates();
                    noalias(auxiliar_vector) = mFocus - r_coordinates;

                    //// Computing the new coordinates
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1] * auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2] * auxiliar_vector[radius_comp_2]);

                    //// We compute the first part of the uplift law
                    mUpliftPressure = mSpecific * ((ref_coord - aux_drain) - (r_coordinates[direction])) * (1.0 - ((1.0 / mDistanceDrain) * (fabs(current_radius - up_radius)))) + (mSpecific * aux_drain);

                    //// If uplift pressure is greater than the limit we compute the second part and we update the value
                    if (mUpliftPressure <= mSpecific * aux_drain)
                    {
                        mUpliftPressure = (mSpecific * ((mReferenceCoordinate + aux_drain) - (r_coordinates[direction]))) * (1.0 - ((1.0 / (width_dam - mDistanceDrain)) * (fabs(current_radius - (up_radius - mDistanceDrain)))));
                    }

                    if (mUpliftPressure < 0.0)
                    {
                        it->FastGetSolutionStepValue(var) = 0.0;
                    }
                    else
                    {
                        it->FastGetSolutionStepValue(var) = mUpliftPressure;
                    }
                }
            }
            else
            {
#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    auxiliar_vector.resize(3, false);
                    const array_1d<double,3>& r_coordinates = it->Coordinates();
                    noalias(auxiliar_vector) = mFocus - r_coordinates;

                    // Computing the current distance to the focus.
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1] * auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2] * auxiliar_vector[radius_comp_2]);

                    mUpliftPressure = mSpecific * (ref_coord - (r_coordinates[direction])) * (1.0 - (1.0 / width_dam) * (fabs(current_radius - up_radius)));

                    if (mUpliftPressure < 0.0)
                    {
                        it->FastGetSolutionStepValue(var) = 0.0;
                    }
                    else
                    {
                        it->FastGetSolutionStepValue(var) = mUpliftPressure;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamUpliftCircularConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamUpliftCircularConditionLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    ModelPart &mrModelPart;
    std::string mVariableName;
    std::string mGravityDirection;
    double mReferenceCoordinate;
    double mSpecific;
    double mWaterLevel;
    bool mDrain;
    double mHeightDrain;
    double mDistanceDrain;
    double mEffectivenessDrain;
    double mUpliftPressure;
    Vector mUpstream;
    Vector mDownstream;
    Vector mFocus;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamUpliftCircularConditionLoadProcess &operator=(DamUpliftCircularConditionLoadProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamUpliftCircularConditionLoadProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamUpliftCircularConditionLoadProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS defined */
