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

#if !defined(KRATOS_GEO_APPLY_CONSTANT_BOUNDARY_PHREATIC_LINE_PRESSURE_PROCESS )
#define  KRATOS_GEO_APPLY_CONSTANT_BOUNDARY_PHREATIC_LINE_PRESSURE_PROCESS

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantBoundaryPhreaticLinePressureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantBoundaryPhreaticLinePressureProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyConstantBoundaryPhreaticLinePressureProcess(ModelPart& model_part,
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
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "first_reference_coordinate":           [0.0,1.0,0.0],
                "second_reference_coordinate":          [1.0,0.5,0.0],
                "specific_weight" : 10000.0,
                "table" : 1
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["first_reference_coordinate"];
        rParameters["second_reference_coordinate"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mGravityDirection = rParameters["gravity_direction"].GetInt();
        mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
        if (mGravityDirection == mOutOfPlaneDirection)
        {
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions "
                         << rParameters
                         << std::endl;
        }
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

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantBoundaryPhreaticLinePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantBoundaryPhreaticLinePressureProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const int nNodes = static_cast<int>(mrModelPart.Nodes().size());

        if (nNodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            Vector3 Coordinates;

            #pragma omp parallel for private(Coordinates)
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if (mIsFixed) it->Fix(var);
                else          it->Free(var);

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
                const double pressure = mSpecificWeight * distance ;

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
        return "ApplyConstantBoundaryPhreaticLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantBoundaryPhreaticLinePressureProcess";
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
    unsigned int mGravityDirection;
    unsigned int mHorizontalDirection;
    double mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    Vector3 mFirstReferenceCoordinate;
    Vector3 mSecondReferenceCoordinate;
    double mSlope;
    double mMinHorizontalCoordinate;
    double mMaxHorizontalCoordinate;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyConstantBoundaryPhreaticLinePressureProcess& operator=(ApplyConstantBoundaryPhreaticLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantPhreaticLinePressureProcess(ApplyConstantBoundaryPhreaticLinePressureProcess const& rOther);

}; // Class ApplyConstantBoundaryPhreaticLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantBoundaryPhreaticLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantBoundaryPhreaticLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_CONSTANT_BOUNDARY_PHREATIC_LINE_PRESSURE_PROCESS defined */
