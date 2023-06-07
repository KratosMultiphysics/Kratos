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

#if !defined(KRATOS_GEO_APPLY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_phreatic_surface_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyPhreaticSurfacePressureTableProcess : public ApplyConstantPhreaticSurfacePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticSurfacePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyPhreaticSurfacePressureTableProcess(ModelPart& model_part,
                                         Parameters rParameters
                                         ) : ApplyConstantPhreaticSurfacePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyPhreaticSurfacePressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyPhreaticSurfacePressureTableProcess algorithms.
    void Execute() override
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfNodes() > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
            const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
            const double deltaH = mpTable->GetValue(Time);

            Vector3 direction=ZeroVector(3);
            direction[mGravityDirection] = 1.0;

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, &direction, &deltaH, this](Node& rNode) {
                    double distance = 0.0;
                    double d = 0.0;
                    for (unsigned int j=0; j < rNode.Coordinates().size(); ++j) {
                        distance += mNormalVector[j]*rNode.Coordinates()[j];
                        d += mNormalVector[j]*direction[j];
                    }
                    distance = -(distance - mEqRHS) / d;
                    distance += deltaH;
                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance ;

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        if (mIsFixedProvided) rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, &direction, &deltaH, this](Node& rNode) {
                    double distance = 0.0;
                    double d = 0.0;
                    for (unsigned int j=0; j < rNode.Coordinates().size(); ++j) {
                        distance += mNormalVector[j]*rNode.Coordinates()[j];
                        d += mNormalVector[j]*direction[j];
                    }
                    distance = -(distance - mEqRHS) / d;
                    distance += deltaH;
                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance ;

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
        return "ApplyPhreaticSurfacePressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyPhreaticSurfacePressureTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    TableType::Pointer mpTable;
    double mTimeUnitConverter;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyPhreaticSurfacePressureTableProcess& operator=(ApplyPhreaticSurfacePressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyPhreaticSurfacePressureTableProcess(ApplyPhreaticSurfacePressureTableProcess const& rOther);
}; // Class ApplyPhreaticSurfacePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyPhreaticSurfacePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyPhreaticSurfacePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS defined */