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

#if !defined(KRATOS_GEO_APPLY_BOUNDARY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_BOUNDARY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_boundary_phreatic_surface_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryPhreaticSurfacePressureTableProcess : public ApplyConstantBoundaryPhreaticSurfacePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryPhreaticSurfacePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyBoundaryPhreaticSurfacePressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryPhreaticSurfacePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyBoundaryPhreaticSurfacePressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyBoundaryPhreaticSurfacePressureTableProcess algorithms.
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
            Vector3 direction=ZeroVector(3);
            direction[mGravityDirection] = 1.0;

            #pragma omp parallel for private(Coordinates)
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                noalias(Coordinates) = it->Coordinates();

                double distance = 0.0;
                double d = 0.0;
                for (unsigned int j=0; j < Coordinates.size(); ++j)
                {
                    distance += mNormalVector[j]*Coordinates[j];
                    d += mNormalVector[j]*direction[j];
                }
                distance = -(distance - mEqRHS) / d;
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
        return "ApplyBoundaryPhreaticSurfacePressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyBoundaryPhreaticSurfacePressureTableProcess";
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
    ApplyBoundaryPhreaticSurfacePressureTableProcess& operator=(ApplyBoundaryPhreaticSurfacePressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyBoundaryPhreaticSurfacePressureTableProcess(ApplyBoundaryPhreaticSurfacePressureTableProcess const& rOther);
}; // Class ApplyBoundaryPhreaticSurfacePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyBoundaryPhreaticSurfacePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyBoundaryPhreaticSurfacePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_BOUNDARY_BOUNDARY_PHREATIC_SURFACE_PRESSURE_TABLE_PROCESS defined */