// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_APPLY_CONSTANT_PHREATIC_LINE_PRESSURE_PROCESS )
#define  KRATOS_GEO_APPLY_CONSTANT_PHREATIC_LINE_PRESSURE_PROCESS

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantPhreaticLinePressureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticLinePressureProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyConstantPhreaticLinePressureProcess(ModelPart& model_part,
                                                Parameters rParameters
                                                ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "is_seepage": false,
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "first_reference_coordinate":           [0.0,1.0,0.0],
                "second_reference_coordinate":          [1.0,0.5,0.0],
                "specific_weight" : 10000.0,
                "pressure_tension_cut_off" : 0.0,
                "table" : [0,1]
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["first_reference_coordinate"];
        rParameters["second_reference_coordinate"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        mIsFixedProvided = rParameters.Has("is_fixed");

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mIsSeepage = rParameters["is_seepage"].GetBool();
        mGravityDirection = rParameters["gravity_direction"].GetInt();
        mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
        if (mGravityDirection == mOutOfPlaneDirection)
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                         << rParameters
                         << std::endl;

        mHorizontalDirection = 0;
        for (unsigned int i=0; i<N_DIM_3D; ++i)
           if (i!=mGravityDirection && i!=mOutOfPlaneDirection) mHorizontalDirection = i;

        mFirstReferenceCoordinate = rParameters["first_reference_coordinate"].GetVector();
        mSecondReferenceCoordinate= rParameters["second_reference_coordinate"].GetVector();

        mMinHorizontalCoordinate = std::min(mFirstReferenceCoordinate[mHorizontalDirection], mSecondReferenceCoordinate[mHorizontalDirection]);
        mMaxHorizontalCoordinate = std::max(mFirstReferenceCoordinate[mHorizontalDirection], mSecondReferenceCoordinate[mHorizontalDirection]);

        if (!(mMaxHorizontalCoordinate > mMinHorizontalCoordinate))
        {
            KRATOS_ERROR << "First and second point on the phreatic line have the same horizontal coordinate"
                         << rParameters
                         << std::endl;
        }

        mSlope = (mSecondReferenceCoordinate[mGravityDirection] - mFirstReferenceCoordinate[mGravityDirection])
                /(mSecondReferenceCoordinate[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]);

        mSpecificWeight = rParameters["specific_weight"].GetDouble();
        if (rParameters.Has("pressure_tension_cut_off"))
          mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
        else
          mPressureTensionCutOff = 0.0;

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantPhreaticLinePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantPhreaticLinePressureProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfNodes() > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
                    const double pressure = CalculatePressure(rNode);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        if (mIsFixedProvided) rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
                    if (mIsFixed) rNode.Fix(var);
                    else if (mIsFixedProvided) rNode.Free(var);
                    const double pressure = CalculatePressure(rNode);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < mPressureTensionCutOff) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                    } else {
                        rNode.FastGetSolutionStepValue(var) = mPressureTensionCutOff;
                    }
                });
            }
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyConstantPhreaticLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantPhreaticLinePressureProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    std::string mVariableName;
    bool mIsFixed;
    bool mIsFixedProvided;
    bool mIsSeepage;
    unsigned int mGravityDirection;
    double mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    Vector3 mFirstReferenceCoordinate;
    Vector3 mSecondReferenceCoordinate;
    double mSlope;
    double mMinHorizontalCoordinate;
    double mMaxHorizontalCoordinate;
    double mPressureTensionCutOff;

    double CalculatePressure(const Node &rNode) const
    {
        double height = 0.0;
        if (rNode.Coordinates()[mHorizontalDirection] >= mMinHorizontalCoordinate && rNode.Coordinates()[mHorizontalDirection] <= mMaxHorizontalCoordinate) {
            height = mSlope * (rNode.Coordinates()[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
        } else if (rNode.Coordinates()[mHorizontalDirection] < mMinHorizontalCoordinate) {
            height = mSlope * (mMinHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
        } else if (rNode.Coordinates()[mHorizontalDirection] > mMaxHorizontalCoordinate) {
            height = mSlope * (mMaxHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
        }

        const double distance = height - rNode.Coordinates()[mGravityDirection];
        const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;
        return pressure;
    }
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyConstantPhreaticLinePressureProcess& operator=(ApplyConstantPhreaticLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantPhreaticLinePressureProcess(ApplyConstantPhreaticLinePressureProcess const& rOther);

}; // Class ApplyConstantPhreaticLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantPhreaticLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantPhreaticLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_CONSTANT_PHREATIC_LINE_PRESSURE_PROCESS defined */
