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

#if !defined(KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS)
#define KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamUpliftConditionLoadProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamUpliftConditionLoadProcess);

    typedef Table<double, double> TableType;

    /// Structs for mapping model parts -------------------------------------------------------------------------------------------------------------------------------------------

    struct GaussPoint
    {
        array_1d<double,3> Coordinates;
        double StateVariable;
    };
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamUpliftConditionLoadProcess(ModelPart &rModelPart,
                                  ModelPart &rJointModelPart,
                                  Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart), mrJointModelPart(rJointModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "Modify"                                           : true,
                "joint_group_name"                                 : "PLEASE_CHOOSE_JOINT_GROUP_NAME",
                "Gravity_Direction"                                : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Upstream_Coordinate"                              : [0.0,0.0,0.0],
                "Downstream_Coordinate"                            : [0.0,0.0,0.0],
                "Upstream_Longitudinal_Coordinate"                 : [0.0,0.0,0.0],
                "Spe_weight"                                       : 0.0,
                "Water_level"                                      : 10.0,
                "Water_Table"                                      : 0,
                "Drains"                                           : false,
                "Height_drain"                                     : 0.0,
                "Distance"                                         : 0.0,
                "Effectiveness"                                    : 0.0,
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["Upstream_Coordinate"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mSpecific = rParameters["Spe_weight"].GetDouble();

        // Getting the values of the coordinates (reference value)
        mX0.resize(3, false);
        mX0[0] = rParameters["Upstream_Coordinate"][0].GetDouble();
        mX0[1] = rParameters["Upstream_Coordinate"][1].GetDouble();
        mX0[2] = rParameters["Upstream_Coordinate"][2].GetDouble();

        mX1.resize(3, false);
        mX1[0] = rParameters["Downstream_Coordinate"][0].GetDouble();
        mX1[1] = rParameters["Downstream_Coordinate"][1].GetDouble();
        mX1[2] = rParameters["Downstream_Coordinate"][2].GetDouble();

        mX2.resize(3, false);
        mX2[0] = rParameters["Upstream_Longitudinal_Coordinate"][0].GetDouble();
        mX2[1] = rParameters["Upstream_Longitudinal_Coordinate"][1].GetDouble();
        mX2[2] = rParameters["Upstream_Longitudinal_Coordinate"][2].GetDouble();

        // Drains
        mDrain = rParameters["Drains"].GetBool();
        mHeightDrain = rParameters["Height_drain"].GetDouble();
        mDistanceDrain = rParameters["Distance"].GetDouble();
        mEffectivenessDrain = rParameters["Effectiveness"].GetDouble();
        mWaterLevel = rParameters["Water_level"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["Water_Table"].GetInt();

        if (mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamUpliftConditionLoadProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {

        KRATOS_TRY;

        //Defining necessary variables
        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nelems = mrJointModelPart.GetMesh(0).Elements().size();
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        BoundedMatrix<double, 3, 3> RotationMatrix;

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTable->GetValue(time);
        }

        // && mrModelPart.GetProcessInfo()[TIME] >= 1.0

        // Computing the rotation matrix accoding with the introduced points by the user
        this->CalculateRotationMatrix(RotationMatrix);
        array_1d<double, 3> newCoordinate;
        array_1d<double, 3> auxiliar_vector;
        array_1d<double, 3> reference_vector;

        // Gravity direction for computing the hydrostatic pressure
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        // Computing the reference vector (coordinates)
        reference_vector = prod(RotationMatrix, mX0);

        // Locate Old Gauss Points in cells
        GaussPoint MyGaussPoint;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        array_1d<double,3> AuxLocalCoordinates;

        double JointPosition = 0.0;
        double ref_water_level = mReferenceCoordinate + mWaterLevel;

        if (nelems != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrJointModelPart.ElementsBegin();

            for(int j = 0; j < nelems; j++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = el_begin + j;

                Element::GeometryType& rGeom = it_elem->GetGeometry();
                MyIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> StateVariableVector(NumGPoints);
                it_elem->CalculateOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfo);

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPoint Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyGaussPoint.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // GaussPoint StateVariable (1 broken, !=1 not broken)
                    if (StateVariableVector.empty()) MyGaussPoint.StateVariable = 0.0;
                    else MyGaussPoint.StateVariable = StateVariableVector[GPoint];

                    if (MyGaussPoint.StateVariable == 1.0) // Broken part
                    {
                        if (MyGaussPoint.Coordinates[0] > JointPosition)
                        {
                            JointPosition = MyGaussPoint.Coordinates[0];
                        }
                    }
                }
            }
        }

        KRATOS_WATCH(JointPosition)

        if (nnodes != 0)
        {

            block_for_each(mrModelPart.GetMesh(0).Nodes(), [&](auto& rNode){

                // Node Global Coordinates
                noalias(auxiliar_vector) = rNode.Coordinates();

                newCoordinate = prod(RotationMatrix, auxiliar_vector);

                double uplift_pressure = 0.0;

                if (newCoordinate(0) < (reference_vector(0) + JointPosition)) // Broken part
                {
                    uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction]));
                }
                else // Not broken part
                {
                    if (mDrain == true && JointPosition < (reference_vector(0) + mDistanceDrain)) // If Drain and Joint broken less than the Drain
                    {
                        if (newCoordinate(0) < (reference_vector(0) + mDistanceDrain)) // UpliftPressure before the Drain
                        {
                            uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - 0.8 * ((1.0 / ((reference_vector(0) + mDistanceDrain) - JointPosition)) * (fabs(newCoordinate(0) - JointPosition))));
                        }
                        else // UpliftPressure after the Drain
                        {
                            uplift_pressure = 0.2 * mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - ((1.0 / (mBaseDam - mDistanceDrain)) * (fabs(newCoordinate(0) - (reference_vector(0) + mDistanceDrain)))));
                        }
                    }
                    else // If Not Drain or Joint broken more than the Drain
                    {
                        uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - ((1.0 / (mBaseDam - (JointPosition - reference_vector(0)))) * (fabs(newCoordinate(0) - JointPosition))));
                    }
                }

                if (uplift_pressure < 0.0)
                {
                    rNode.FastGetSolutionStepValue(var) = 0.0;
                }
                else
                {
                    rNode.FastGetSolutionStepValue(var) = uplift_pressure;
                }
            });
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        //Defining necessary variables
        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nelems = mrJointModelPart.GetMesh(0).Elements().size();
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        BoundedMatrix<double, 3, 3> RotationMatrix;

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTable->GetValue(time);
        }

        // && mrModelPart.GetProcessInfo()[TIME] >= 1.0

        // Computing the rotation matrix accoding with the introduced points by the user
        this->CalculateRotationMatrix(RotationMatrix);
        array_1d<double, 3> newCoordinate;
        array_1d<double, 3> auxiliar_vector;
        array_1d<double, 3> reference_vector;

        // Gravity direction for computing the hydrostatic pressure
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        // Computing the reference vector (coordinates)
        reference_vector = prod(RotationMatrix, mX0);

        // Locate Old Gauss Points in cells
        GaussPoint MyGaussPoint;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        array_1d<double,3> AuxLocalCoordinates;

        double JointPosition = 0.0;
        double ref_water_level = mReferenceCoordinate + mWaterLevel;

        if (nelems != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrJointModelPart.ElementsBegin();

            for(int j = 0; j < nelems; j++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = el_begin + j;

                Element::GeometryType& rGeom = it_elem->GetGeometry();
                MyIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> StateVariableVector(NumGPoints);
                it_elem->CalculateOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfo);

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPoint Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyGaussPoint.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // GaussPoint StateVariable (1 broken, !=1 not broken)
                    if (StateVariableVector.empty()) MyGaussPoint.StateVariable = 0.0;
                    else MyGaussPoint.StateVariable = StateVariableVector[GPoint];

                    if (MyGaussPoint.StateVariable == 1.0) // Broken part
                    {
                        if (MyGaussPoint.Coordinates[0] > JointPosition)
                        {
                            JointPosition = MyGaussPoint.Coordinates[0];
                        }
                    }
                }
            }
        }

        KRATOS_WATCH(JointPosition)

        if (nnodes != 0)
        {

            block_for_each(mrModelPart.GetMesh(0).Nodes(), [&](auto& rNode){

                // Node Global Coordinates
                noalias(auxiliar_vector) = rNode.Coordinates();

                newCoordinate = prod(RotationMatrix, auxiliar_vector);

                double uplift_pressure = 0.0;

                if (newCoordinate(0) < (reference_vector(0) + JointPosition)) // Broken part
                {
                    uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction]));
                }
                else // Not broken part
                {
                    if (mDrain == true && JointPosition < (reference_vector(0) + mDistanceDrain)) // If Drain and Joint broken less than the Drain
                    {
                        if (newCoordinate(0) < (reference_vector(0) + mDistanceDrain)) // UpliftPressure before the Drain
                        {
                            uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - 0.8 * ((1.0 / ((reference_vector(0) + mDistanceDrain) - JointPosition)) * (fabs(newCoordinate(0) - JointPosition))));
                        }
                        else // UpliftPressure after the Drain
                        {
                            uplift_pressure = 0.2 * mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - ((1.0 / (mBaseDam - mDistanceDrain)) * (fabs(newCoordinate(0) - (reference_vector(0) + mDistanceDrain)))));
                        }
                    }
                    else // If Not Drain or Joint broken more than the Drain
                    {
                        uplift_pressure = mSpecific * (ref_water_level - (rNode.Coordinates()[direction])) * (1.0 - ((1.0 / (mBaseDam - (JointPosition - reference_vector(0)))) * (fabs(newCoordinate(0) - JointPosition))));
                    }
                }

                if (uplift_pressure < 0.0)
                {
                    rNode.FastGetSolutionStepValue(var) = 0.0;
                }
                else
                {
                    rNode.FastGetSolutionStepValue(var) = uplift_pressure;
                }
            });
        }

        KRATOS_CATCH("");
    }

    void CalculateRotationMatrix(BoundedMatrix<double, 3, 3> &rRotationMatrix)
    {
        KRATOS_TRY;

        //Unitary vector in uplift direction
        array_1d<double, 3> V_uplift;
        V_uplift = (mX1 - mX0);
        mBaseDam = norm_2(V_uplift);
        double inv_norm_uplift = 1.0 / norm_2(V_uplift);
        V_uplift[0] *= inv_norm_uplift;
        V_uplift[1] *= inv_norm_uplift;
        V_uplift[2] *= inv_norm_uplift;

        //Unitary vector in longitudinal direction
        array_1d<double, 3> V_longitudinal;
        V_longitudinal = (mX2 - mX0);
        double inv_norm_longitudinal = 1.0 / norm_2(V_longitudinal);
        V_longitudinal[0] *= inv_norm_longitudinal;
        V_longitudinal[1] *= inv_norm_longitudinal;
        V_longitudinal[2] *= inv_norm_longitudinal;

        //Unitary vector in local z direction
        array_1d<double, 3> V_normal;
        MathUtils<double>::CrossProduct(V_normal, V_uplift, V_longitudinal);

        //Rotation Matrix
        rRotationMatrix(0, 0) = V_uplift[0];
        rRotationMatrix(0, 1) = V_uplift[1];
        rRotationMatrix(0, 2) = V_uplift[2];

        rRotationMatrix(1, 0) = V_longitudinal[0];
        rRotationMatrix(1, 1) = V_longitudinal[1];
        rRotationMatrix(1, 2) = V_longitudinal[2];

        rRotationMatrix(2, 0) = V_normal[0];
        rRotationMatrix(2, 1) = V_normal[1];
        rRotationMatrix(2, 2) = V_normal[2];

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamUpliftConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamUpliftConditionLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    ModelPart &mrModelPart;
    ModelPart &mrJointModelPart;
    std::string mVariableName;
    std::string mGravityDirection;
    double mReferenceCoordinate;
    double mSpecific;
    double mBaseDam;
    double mWaterLevel;
    bool mDrain;
    double mHeightDrain;
    double mDistanceDrain;
    double mEffectivenessDrain;
    Vector mX0;
    Vector mX1;
    Vector mX2;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamUpliftConditionLoadProcess &operator=(DamUpliftConditionLoadProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamUpliftConditionLoadProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamUpliftConditionLoadProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS defined */
