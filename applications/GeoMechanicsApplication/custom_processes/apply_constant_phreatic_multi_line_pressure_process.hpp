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
        if (mGravityDirection == mOutOfPlaneDirection)
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                         << rParameters
                         << std::endl;

        mHorizontalDirection = 0;
        for (unsigned int i=0; i<N_DIM_3D; ++i)
           if (i!=mGravityDirection && i!=mOutOfPlaneDirection) mHorizontalDirection = i;

        mXCoordinates = rParameters["x_coordinates"].GetVector();
        mYCoordinates= rParameters["y_coordinates"].GetVector();
        mZCoordinates = rParameters["z_coordinates"].GetVector();

        if (mHorizontalDirection == 0) {
            mHorizontalDirectionCoordinates = mXCoordinates;
        }
        else if (mHorizontalDirection == 1)
        {
            mHorizontalDirectionCoordinates = mYCoordinates;
        }
        else if (mHorizontalDirection == 2)
        {
            mHorizontalDirectionCoordinates = mZCoordinates;
        }

        if (mGravityDirection == 0) {
            mGravityDirectionCoordinates = mXCoordinates;
        }
        else if (mGravityDirection == 1)
        {
            mGravityDirectionCoordinates = mYCoordinates;
        }
        else if (mGravityDirection == 2)
        {
            mGravityDirectionCoordinates = mZCoordinates;
        }

    	bool mHorizontalDirectionCoordindatesSorted = std::is_sorted(mHorizontalDirectionCoordinates.begin(), mHorizontalDirectionCoordinates.end());

    	if (!mHorizontalDirectionCoordindatesSorted)
        {
            KRATOS_ERROR << "The Horizontal Elements Coordinates are not ordered."
                         << rParameters
                         << std::endl;
        }

        mMinHorizontalCoordinate = *std::min_element(mHorizontalDirectionCoordinates.begin(), mHorizontalDirectionCoordinates.end());
        mMaxHorizontalCoordinate = *std::max_element(mHorizontalDirectionCoordinates.begin(), mHorizontalDirectionCoordinates.end());

        mSpecificWeight = rParameters["specific_weight"].GetDouble();
        if (rParameters.Has("pressure_tension_cut_off"))
          mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
        else
          mPressureTensionCutOff = 0.0;

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantPhreaticMultiLinePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantPhreaticMultiLinePressureProcess algorithms.
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
                block_for_each(mrModelPart.Nodes(), [&var, this](Node<3>& rNode) {
                    const double pressure = CalculatePressure(rNode);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, this](Node<3>& rNode) {
                    if (mIsFixed) rNode.Fix(var);
                    else          rNode.Free(var);
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
        return "ApplyConstantPhreaticMultiLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantPhreaticMultiLinePressureProcess";
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
    
    double mMinHorizontalCoordinate;
    double mMaxHorizontalCoordinate;
    double mPressureTensionCutOff;

    int findIndex(const Node<3>& rNode) const
    {
        int index = 0;
    	for(index; index < mHorizontalDirectionCoordinates.size(); index++)
        {
	        if (mHorizontalDirectionCoordinates[index] >= rNode.Coordinates()[mHorizontalDirection])
	        {
	        	break;
	        }
        }
        if (index == 0) 
        {
            return 0;
        }
        else
        {
            return index - 1;
        }
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

        slope = (mGravityDirectionCoordinates[secondPointIndex] - mGravityDirectionCoordinates[firstPointIndex]) /
            (mHorizontalDirectionCoordinates[secondPointIndex] - mHorizontalDirectionCoordinates[firstPointIndex]);

    	height = slope * (rNode.Coordinates()[mHorizontalDirection] - mHorizontalDirectionCoordinates[firstPointIndex]) + mGravityDirectionCoordinates[firstPointIndex];

        const double distance = height - rNode.Coordinates()[mGravityDirection];
        const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * distance;
        return pressure;
    }
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyConstantPhreaticLinePressureProcess& operator=(ApplyConstantPhreaticMultiLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantPhreaticLinePressureProcess(ApplyConstantPhreaticMultiLinePressureProcess const& rOther);

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
