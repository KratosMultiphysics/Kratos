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

#if !defined(KRATOS_GEO_APPLY_BOUNDARY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_BOUNDARY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_boundary_phreatic_line_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryPhreaticLinePressureTableProcess : public ApplyConstantBoundaryPhreaticLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyBoundaryPhreaticLinePressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryPhreaticLinePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyBoundaryPhreaticLinePressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyBoundaryPhreaticLinePressureTableProcess algorithms.
    void Execute() override
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
        double deltaH = mpTable->GetValue(Time);

        const int nNodes = static_cast<int>(mrModelPart.Nodes().size());

        if (nNodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            Vector3 Coordinates;

            #pragma omp parallel for private(Coordinates)
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                noalias(Coordinates) = it->Coordinates();

                double hight = 0.0;
                if (Coordinates[mHorizontalDirection] >= mMinHorizontalCoordinate && Coordinates[mHorizontalDirection] <= mMaxHorizontalCoordinate)
                {
                    hight = mSlope * (Coordinates[mHorizontalDirection] - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
                } else if (Coordinates[mHorizontalDirection] < mMinHorizontalCoordinate)
                {
                    hight = mSlope * (mMinHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
                } else if (Coordinates[mHorizontalDirection] > mMaxHorizontalCoordinate)
                {
                    hight = mSlope * (mMaxHorizontalCoordinate - mFirstReferenceCoordinate[mHorizontalDirection]) + mFirstReferenceCoordinate[mGravityDirection];
                }

                double distance = hight - Coordinates[mGravityDirection];
                distance += deltaH;
                const double pressure = mSpecificWeight * distance;

                if (pressure > 0.0)
                {
                    it->FastGetSolutionStepValue(var) = pressure;
                }
                else
                {
                    it->FastGetSolutionStepValue(var) = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyBoundaryPhreaticLinePressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyBoundaryPhreaticLinePressureTableProcess";
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
    ApplyBoundaryPhreaticLinePressureTableProcess& operator=(ApplyBoundaryPhreaticLinePressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyBoundaryPhreaticLinePressureTableProcess(ApplyBoundaryPhreaticLinePressureTableProcess const& rOther);
}; // Class ApplyBoundaryPhreaticLinePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyBoundaryPhreaticLinePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyBoundaryPhreaticLinePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_BOUNDARY_BOUNDARY_PHREATIC_LINE_PRESSURE_TABLE_PROCESS defined */