// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_APPLY_CONSTANT_PHREATIC_SURFACE_PRESSURE_PROCESS )
#define  KRATOS_GEO_APPLY_CONSTANT_PHREATIC_SURFACE_PRESSURE_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantPhreaticSurfacePressureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticSurfacePressureProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyConstantPhreaticSurfacePressureProcess(ModelPart& model_part,
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
                "first_reference_coordinate":          [0.0,1.0,0.0],
                "second_reference_coordinate":         [1.0,0.5,0.0],
                "third_reference_coordinate":          [1.0,0.5,0.0],
                "specific_weight" : 10000.0,
                "pressure_tension_cut_off" : 0.0,
                "table" : [0,1,2]
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["first_reference_coordinate"];
        rParameters["second_reference_coordinate"];
        rParameters["third_reference_coordinate"];

        rParameters["variable_name"];
        rParameters["model_part_name"];

        mIsFixedProvided = rParameters.Has("is_fixed");

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mIsSeepage = rParameters["is_seepage"].GetBool();
        mGravityDirection = rParameters["gravity_direction"].GetInt();
        mFirstReferenceCoordinate = rParameters["first_reference_coordinate"].GetVector();
        mSecondReferenceCoordinate= rParameters["second_reference_coordinate"].GetVector();
        mThirdReferenceCoordinate= rParameters["third_reference_coordinate"].GetVector();

        calculateEquationParameters();

        mSpecificWeight = rParameters["specific_weight"].GetDouble();

        if (rParameters.Has("pressure_tension_cut_off"))
          mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
        else
          mPressureTensionCutOff = 0.0;

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantPhreaticSurfacePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantPhreaticSurfacePressureProcess algorithms.
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
            Vector3 direction=ZeroVector(3);
            direction[mGravityDirection] = 1.0;

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, &direction, this](Node& rNode) {
                    double distance = 0.0;
                    double d = 0.0;
                    for (unsigned int j=0; j < rNode.Coordinates().size(); ++j) {
                        distance += mNormalVector[j] * rNode.Coordinates()[j];
                        d += mNormalVector[j]*direction[j];
                    }
                    distance = -(distance - mEqRHS) / d;

                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        if (mIsFixedProvided) rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, &direction, this](Node& rNode) {
                    if (mIsFixed) rNode.Fix(var);
                    else if (mIsFixedProvided) rNode.Free(var);

                    double distance = 0.0;
                    double d = 0.0;
                    for (unsigned int j=0; j < rNode.Coordinates().size(); ++j) {
                        distance += mNormalVector[j] * rNode.Coordinates()[j];
                        d += mNormalVector[j]*direction[j];
                    }
                    distance = -(distance - mEqRHS) / d;

                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;

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
        return "ApplyConstantPhreaticSurfacePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantPhreaticSurfacePressureProcess";
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
    Vector3 mFirstReferenceCoordinate;
    Vector3 mSecondReferenceCoordinate;
    Vector3 mThirdReferenceCoordinate;
    Vector3 mNormalVector;
    double mEqRHS;
    double mPressureTensionCutOff;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void calculateEquationParameters()
    {
        Vector3 v1;
        Vector3 v2;
        v1 = mSecondReferenceCoordinate - mFirstReferenceCoordinate;
        v2 = mThirdReferenceCoordinate - mFirstReferenceCoordinate;

        mNormalVector[0] = v1[1] * v2[2] - v1[2] * v2[1];
        mNormalVector[1] = v1[2] * v2[0] - v1[0] * v2[2];
        mNormalVector[2] = v1[0] * v2[1] - v1[1] * v2[0];

        double sum=0;
        for (unsigned int i=0; i < mNormalVector.size(); ++i)
        {
            sum += mNormalVector[i]*mNormalVector[i];
        }
        double norm = std::sqrt(sum);
        if (!(norm > 0.0))
            KRATOS_ERROR << "Normal vector to phreatic surface has zero size!"
                         << std::endl;

        mEqRHS = 0.0;
        for (unsigned int i=0; i < mFirstReferenceCoordinate.size(); ++i)
        {
            mEqRHS += mNormalVector[i]*mFirstReferenceCoordinate[i];
        }
    }

    /// Assignment operator.
    ApplyConstantPhreaticSurfacePressureProcess& operator=(ApplyConstantPhreaticSurfacePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantPhreaticSurfacePressureProcess(ApplyConstantPhreaticSurfacePressureProcess const& rOther);

}; // Class ApplyConstantPhreaticSurfacePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantPhreaticSurfacePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantPhreaticSurfacePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_CONSTANT_PHREATIC_SURFACE_PRESSURE_PROCESS defined */
