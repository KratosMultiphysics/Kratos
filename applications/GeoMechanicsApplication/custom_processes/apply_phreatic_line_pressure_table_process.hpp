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

#if !defined(KRATOS_GEO_APPLY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_phreatic_line_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyPhreaticLinePressureTableProcess : public ApplyConstantPhreaticLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyPhreaticLinePressureTableProcess(ModelPart& model_part,
                                         Parameters rParameters
                                         ) : ApplyConstantPhreaticLinePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        for (unsigned int i=0; i < mpTable.size(); ++i) {
            unsigned int TableId = rParameters["table"][i].GetInt();
            if (TableId > 0) {
                mpTable[i] = model_part.pGetTable(TableId);
            } else {
                mpTable[i] = nullptr;
            }
        }

        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyPhreaticLinePressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyPhreaticLinePressureTableProcess algorithms.
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
            array_1d<double, 2> deltaH;
            for (unsigned int i=0; i < mpTable.size(); ++i) {
                if (!mpTable[i]) {
                deltaH[i] = 0.0;
                } else {
                deltaH[i] = mpTable[i]->GetValue(Time);
                }
            }

            array_1d<double, 2> y;
            y[0] = deltaH[0] + mFirstReferenceCoordinate[mGravityDirection];
            y[1] = deltaH[1] + mSecondReferenceCoordinate[mGravityDirection];

            mSlope = (y[1] - y[0])
                    /(mSecondReferenceCoordinate[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]);

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, &y, this](Node<3>& rNode) {
                    const double pressure = CalculatePressure(rNode, y);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, &y, this](Node<3>& rNode) {
                    const double pressure = CalculatePressure(rNode, y);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < mPressureTensionCutOff) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                    } else {
                        rNode.FastGetSolutionStepValue(var) = mPressureTensionCutOff;
                    }
                });
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyPhreaticLinePressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyPhreaticLinePressureTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    array_1d<TableType::Pointer,2> mpTable;
    double mTimeUnitConverter;

    double CalculatePressure(const Node<3> &rNode, const array_1d<double, 2> &y) const
    {
        double height = 0.0;
        if (rNode.Coordinates()[mHorizontalDirection] >= mMinHorizontalCoordinate && rNode.Coordinates()[mHorizontalDirection] <= mMaxHorizontalCoordinate) {
            height = mSlope * (rNode.Coordinates()[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]) + y[0];

        } else if (rNode.Coordinates()[mHorizontalDirection] < mMinHorizontalCoordinate) {
            height = mSlope * (mMinHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + y[0];

        } else if (rNode.Coordinates()[mHorizontalDirection] > mMaxHorizontalCoordinate) {
            height = mSlope * (mMaxHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + y[0];
        }

        const double distance = height - rNode.Coordinates()[mGravityDirection];
        const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;
        return pressure;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyPhreaticLinePressureTableProcess& operator=(ApplyPhreaticLinePressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyPhreaticLinePressureTableProcess(ApplyPhreaticLinePressureTableProcess const& rOther);
}; // Class ApplyPhreaticLinePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyPhreaticLinePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyPhreaticLinePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS defined */