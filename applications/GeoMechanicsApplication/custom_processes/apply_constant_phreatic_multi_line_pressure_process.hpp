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
//                   Jonathan Nuttall

#if !defined(KRATOS_GEO_APPLY_CONSTANT_PHREATIC_MULTI_LINE_PRESSURE_PROCESS )
#define  KRATOS_GEO_APPLY_CONSTANT_PHREATIC_MULTI_LINE_PRESSURE_PROCESS

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantPhreaticMultiLinePressureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticMultiLinePressureProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyConstantPhreaticMultiLinePressureProcess(ModelPart& model_part,
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
                "x_coordinates":           [0.0,1.0],
                "y_coordinates":           [1.0,0.5],
				"z_coordinates":           [0.0,0.0],
                "specific_weight" : 10000.0,
                "pressure_tension_cut_off" : 0.0,
                "table" : [0,1]
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["x_coordinates"];
        rParameters["y_coordinates"];
        rParameters["z_coordinates"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mIsSeepage = rParameters["is_seepage"].GetBool();
        mGravityDirection = rParameters["gravity_direction"].GetInt();
        mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
        if (GravityDirection() == OutOfPlaneDirection())
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                         << rParameters
                         << std::endl;

        mHorizontalDirection = 0;
        for (unsigned int i=0; i<N_DIM_3D; ++i)
           if (i!=GravityDirection() && i != OutOfPlaneDirection()) mHorizontalDirection = i;

        mXCoordinates = rParameters["x_coordinates"].GetVector();
        mYCoordinates = rParameters["y_coordinates"].GetVector();
        mZCoordinates = rParameters["z_coordinates"].GetVector();

        if (HorizontalDirection() == 0) {
            mHorizontalDirectionCoordinates = XCoordinates();
        }
        else if (HorizontalDirection() == 1)
        {
            mHorizontalDirectionCoordinates = YCoordinates();
        }
        else if (HorizontalDirection() == 2)
        {
            mHorizontalDirectionCoordinates = ZCoordinates();
        }

        if (GravityDirection() == 0) {
            mGravityDirectionCoordinates = XCoordinates();
        }
        else if (GravityDirection() == 1)
        {
            mGravityDirectionCoordinates = YCoordinates();
        }
        else if (GravityDirection() == 2)
        {
            mGravityDirectionCoordinates = ZCoordinates();
        }

        if (!std::is_sorted(HorizontalDirectionCoordinates().begin(), HorizontalDirectionCoordinates().end()))
        {
            KRATOS_ERROR << "The Horizontal Elements Coordinates are not ordered."
                         << rParameters
                         << std::endl;
        }

        mSpecificWeight = rParameters["specific_weight"].GetDouble();
        if (rParameters.Has("pressure_tension_cut_off"))
          mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
        else
          mPressureTensionCutOff = 0.0;

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantPhreaticMultiLinePressureProcess() override = default;

    /// Assignment operator.
    ApplyConstantPhreaticMultiLinePressureProcess& operator=(ApplyConstantPhreaticMultiLinePressureProcess const&) = delete;

    /// Copy constructor.
    ApplyConstantPhreaticMultiLinePressureProcess(ApplyConstantPhreaticMultiLinePressureProcess const&) = delete;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfNodes() <= 0) {
            return;
        }

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(VariableName());

        if (IsSeepage()) {
            block_for_each(mrModelPart.Nodes(), [&var, this](Node<3>& rNode) {
                const double pressure = CalculatePressure(rNode);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                    if (IsFixed()) rNode.Fix(var);
                } else {
                    rNode.Free(var);
                }
            });
        } else {
            block_for_each(mrModelPart.Nodes(), [&var, this](Node<3>& rNode) {
                if (IsFixed()) rNode.Fix(var);
                else           rNode.Free(var);
                const double pressure = CalculatePressure(rNode);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < PressureTensionCutOff()) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                } else {
                    rNode.FastGetSolutionStepValue(var) = PressureTensionCutOff();
                }
            });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyConstantPhreaticMultiLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantPhreaticMultiLinePressureProcess";
    }

    const std::string& VariableName() const
    {
        return mVariableName;
    }

    bool IsFixed() const
    {
        return mIsFixed;
    }

    bool IsSeepage() const
    {
        return mIsSeepage;
    }

    unsigned int GravityDirection() const
    {
        return mGravityDirection;
    }

    double SpecificWeight() const
    {
        return mSpecificWeight;
    }

    unsigned int OutOfPlaneDirection() const
    {
        return mOutOfPlaneDirection;
    }

    unsigned int HorizontalDirection() const
    {
        return mHorizontalDirection;
    }

    const Vector& HorizontalDirectionCoordinates() const
    {
        return mHorizontalDirectionCoordinates;
    }

    const Vector& GravityDirectionCoordinates() const
    {
        return mGravityDirectionCoordinates;
    }

    const Vector& XCoordinates() const
    {
        return mXCoordinates;
    }

    const Vector& YCoordinates() const
    {
        return mYCoordinates;
    }

    const Vector& ZCoordinates() const
    {
        return mZCoordinates;
    }

    double PressureTensionCutOff() const
    {
        return mPressureTensionCutOff;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;

    int findIndex(const Node<3>& rNode) const
    {
        int index = 0;
        auto coords = rNode.Coordinates();
        auto noCoordinates = static_cast<int>(HorizontalDirectionCoordinates().size());
    	for ( ; index < noCoordinates; ++index)
        {
	        if (HorizontalDirectionCoordinates()[index] >= coords[HorizontalDirection()])
	        {
                if (index == 0)
                {
                    return 0;
                }
                else
                {
                    return index - 1;
                }
	        }
        }
        return noCoordinates - 1;
        
    }

    double CalculatePressure(const Node<3> &rNode) const
    {
        double height = 0.0;
        int firstPointIndex;
        int secondPointIndex;
        double slope;

        // find nodes in horizontalDirectionCoordinates

        firstPointIndex = findIndex(rNode);
        secondPointIndex = firstPointIndex + 1;

        slope = (GravityDirectionCoordinates()[secondPointIndex] - GravityDirectionCoordinates()[firstPointIndex]) /
            (HorizontalDirectionCoordinates()[secondPointIndex] - HorizontalDirectionCoordinates()[firstPointIndex]);

        height = slope * (rNode.Coordinates()[HorizontalDirection()] - HorizontalDirectionCoordinates()[firstPointIndex]) + GravityDirectionCoordinates()[firstPointIndex];

        const double distance = height - rNode.Coordinates()[GravityDirection()];
        const double pressure = - PORE_PRESSURE_SIGN_FACTOR * SpecificWeight() * distance;
        return pressure;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    std::string mVariableName;
    bool mIsFixed;
    bool mIsSeepage;
    unsigned int mGravityDirection;
    double mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    Vector mHorizontalDirectionCoordinates;
    Vector mGravityDirectionCoordinates;
    Vector mXCoordinates;
    Vector mYCoordinates;
    Vector mZCoordinates;
    double mPressureTensionCutOff;
}; // Class ApplyConstantPhreaticMultiLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantPhreaticMultiLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantPhreaticMultiLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_CONSTANT_PHREATIC_MULTI_LINE_PRESSURE_PROCESS defined */
