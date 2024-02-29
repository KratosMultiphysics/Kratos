
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_PROPAGATE_FRACTURES_3D_UTILITIES )
#define  KRATOS_PROPAGATE_FRACTURES_3D_UTILITIES

// System includes
#include <fstream>
#include <iostream>
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class FracturePropagation3DUtilities
{

protected:

/// Basic Structs for the utility ---------------------------------------------------------------------------------------------------------------------------------------------

    struct UtilityVariables
    {
        double X_max, X_min, Y_max, Y_min, Z_max, Z_min;
        int NRows, NColumns, NSections;
        double RowSize, ColumnSize, SectionSize;
    };

/// Structs for fracture propagation check ------------------------------------------------------------------------------------------------------------------------------------

    struct FracturePoint
    {
        array_1d<double,3> Coordinates;
        double Damage, Weight, TipDistance;
        Element::Pointer pElement;
    };

    ///------------------------------------------------------------------------------------

    struct Propagation
    {
        int MotherFractureId;
        array_1d<double,3> TopInitCoordinates;
        array_1d<double,3> BotInitCoordinates;
        array_1d<double,3> LeftInitCoordinates;
        array_1d<double,3> RightInitCoordinates;
        array_1d<double,3> TopEndCoordinates;
        array_1d<double,3> BotEndCoordinates;
        array_1d<double,3> LeftEndCoordinates;
        array_1d<double,3> RightEndCoordinates;
        array_1d<double,3> TipCoordinates;
    };

    ///------------------------------------------------------------------------------------

    struct Bifurcation
    {
        int MotherFractureId;
        array_1d<double,3> TopInitCoordinates;
        array_1d<double,3> BotInitCoordinates;
        array_1d<double,3> LeftInitCoordinates;
        array_1d<double,3> RightInitCoordinates;
        array_1d<double,3> TopTopEndCoordinates;
        array_1d<double,3> TopBotEndCoordinates;
        array_1d<double,3> TopLeftEndCoordinates;
        array_1d<double,3> TopRightEndCoordinates;
        array_1d<double,3> TopTipCoordinates;
        array_1d<double,3> BotTopEndCoordinates;
        array_1d<double,3> BotBotEndCoordinates;
        array_1d<double,3> BotLeftEndCoordinates;
        array_1d<double,3> BotRightEndCoordinates;
        array_1d<double,3> BotTipCoordinates;
    };

    ///------------------------------------------------------------------------------------

    struct PropagationGlobalVariables
    {
        std::vector< std::vector< std::vector< std::vector<FracturePoint> > > > FracturePointsCellMatrix;
        ProcessInfo::Pointer pProcessInfo;
        bool PropagateFractures;
        std::vector<Propagation> PropagationVector;
        std::vector<Bifurcation> BifurcationVector;
    };

    ///------------------------------------------------------------------------------------

    struct PropagationLocalVariables
    {
        BoundedMatrix<double,3,3> RotationMatrix;
        array_1d<double,3> TipCoordinates;
        array_1d<double,3> TipLocalCoordinates;
        std::vector<FracturePoint*> TopFrontFracturePoints;
        std::vector<FracturePoint*> BotFrontFracturePoints;
    };

/// Structs for mapping model parts -------------------------------------------------------------------------------------------------------------------------------------------

    struct GaussPointOld
    {
        array_1d<double,3> Coordinates;
        double StateVariable, Weight;
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION( FracturePropagation3DUtilities );

    /// Constructor
    FracturePropagation3DUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~FracturePropagation3DUtilities() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool CheckFracturePropagation (Parameters& rParameters, ModelPart& rModelPart, const bool& move_mesh_flag)
    {
        // Define necessary variables
        UtilityVariables AuxVariables;
        PropagationGlobalVariables PropagationData;

        this->InitializeCheckFracture(PropagationData, AuxVariables, rParameters, rModelPart, move_mesh_flag);

        //Loop for all pre-existing fractures
        for(unsigned int i = 0; i < rParameters["fractures_list"].size(); i++)
        {
            this->CheckFracture(i, PropagationData, AuxVariables, rParameters);
        }

        this->FinalizeCheckFracture(PropagationData, rParameters, rModelPart, move_mesh_flag);

        return PropagationData.PropagateFractures;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void MappingModelParts (Parameters& rParameters, ModelPart& rModelPartOld, ModelPart& rModelPartNew, const bool& move_mesh_flag)
    {
        // Define necessary variables
        UtilityVariables AuxVariables;

        this->InitializeMapping(AuxVariables,rModelPartNew, move_mesh_flag);

        this->NodalVariablesMapping(AuxVariables,rParameters,rModelPartOld,rModelPartNew);

        this->GaussPointStateVariableMapping(AuxVariables,rParameters,rModelPartOld,rModelPartNew);

        if(move_mesh_flag==true)
            this->UpdateMeshPosition(rModelPartNew);
    }

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/// Fracture Propagation Check ------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeCheckFracture(
        PropagationGlobalVariables& rPropagationData,
        UtilityVariables& rAuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPart,
        const bool& move_mesh_flag)
    {
        // Set PropagateFractures bool to false
        rPropagationData.PropagateFractures = false;

        // Move mesh to the original position to work in the reference state
        if(move_mesh_flag==true)
            this->ResetMeshPosition(rModelPart);

        // Locate fracture points inside CellMatrix
        this->SetFracturePoints(rPropagationData, rAuxVariables, rParameters, rModelPart);
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CheckFracture(
        const unsigned int& itFracture,
        PropagationGlobalVariables& rPropagationData,
        const UtilityVariables& AuxVariables,
        Parameters& rParameters)
    {
        PropagationLocalVariables AuxPropagationVariables;

        // Original tip Coordinates
        array_1d<double,3> LeftPointCoordinates;
        array_1d<double,3> RightPointCoordinates;
        for(unsigned int i = 0; i < 3; i++)
        {
            AuxPropagationVariables.TipCoordinates[i] = rParameters["fractures_list"][itFracture]["tip_point"]["coordinates"][i].GetDouble();
            LeftPointCoordinates[i] = rParameters["fractures_list"][itFracture]["left_point"]["coordinates"][i].GetDouble();
            RightPointCoordinates[i] = rParameters["fractures_list"][itFracture]["right_point"]["coordinates"][i].GetDouble();
        }

        // Tip Rotation Matrix
        this->CalculateTipRotationMatrix(AuxPropagationVariables.RotationMatrix,AuxPropagationVariables.TipCoordinates,LeftPointCoordinates,RightPointCoordinates);

        // Tip search area
        const double PropagationLength = rParameters["fracture_data"]["propagation_length"].GetDouble();

        double X_left = AuxPropagationVariables.TipCoordinates[0] - PropagationLength;
        double X_right = AuxPropagationVariables.TipCoordinates[0] + PropagationLength;
        double Y_top = AuxPropagationVariables.TipCoordinates[1] + PropagationLength;
        double Y_bot = AuxPropagationVariables.TipCoordinates[1] - PropagationLength;
        double Z_back = AuxPropagationVariables.TipCoordinates[2] - PropagationLength;
        double Z_front = AuxPropagationVariables.TipCoordinates[2] + PropagationLength;

        int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
        int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
        int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
        int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);
        int Section_back = int((Z_back-AuxVariables.Z_min)/AuxVariables.SectionSize);
        int Section_front = int((Z_front-AuxVariables.Z_min)/AuxVariables.SectionSize);

        if(Column_left < 0) Column_left = 0;
        if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
        if(Row_top < 0) Row_top = 0;
        if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
        if(Section_back < 0) Section_back = 0;
        if(Section_front >= AuxVariables.NSections) Section_front = AuxVariables.NSections-1;

        // Search FracturePoints neighbours around the tip
        std::vector<FracturePoint*> TipNeighbours;
        noalias(AuxPropagationVariables.TipLocalCoordinates) = prod(AuxPropagationVariables.RotationMatrix,AuxPropagationVariables.TipCoordinates);
        array_1d<double,3> OtherLocalCoordinates;
        for(int s = Section_back; s <= Section_front; s++)
        {
            for(int i = Row_top; i <= Row_bot; i++)
            {
                for(int j = Column_left; j<= Column_right; j++)
                {
                    for(unsigned int k = 0; k < rPropagationData.FracturePointsCellMatrix[i][j][s].size(); k++)
                    {
                        FracturePoint& rOtherPoint = rPropagationData.FracturePointsCellMatrix[i][j][s][k];

                        rOtherPoint.TipDistance = sqrt((rOtherPoint.Coordinates[0]-AuxPropagationVariables.TipCoordinates[0])*(rOtherPoint.Coordinates[0]-AuxPropagationVariables.TipCoordinates[0]) +
                                        (rOtherPoint.Coordinates[1]-AuxPropagationVariables.TipCoordinates[1])*(rOtherPoint.Coordinates[1]-AuxPropagationVariables.TipCoordinates[1]) +
                                        (rOtherPoint.Coordinates[2]-AuxPropagationVariables.TipCoordinates[2])*(rOtherPoint.Coordinates[2]-AuxPropagationVariables.TipCoordinates[2]));

                        if(rOtherPoint.TipDistance <= PropagationLength)
                        {
                            TipNeighbours.push_back(&rOtherPoint);

                            noalias(OtherLocalCoordinates) = prod(AuxPropagationVariables.RotationMatrix,rOtherPoint.Coordinates);

                            // FrontFracturePoints
                            if(OtherLocalCoordinates[0] >= AuxPropagationVariables.TipLocalCoordinates[0])
                            {
                                // TopFrontFracturePoints
                                if(OtherLocalCoordinates[2] >= AuxPropagationVariables.TipLocalCoordinates[2])
                                {
                                    AuxPropagationVariables.TopFrontFracturePoints.push_back(&rOtherPoint);
                                }
                                // BotFrontFracturePoints
                                else
                                {
                                    AuxPropagationVariables.BotFrontFracturePoints.push_back(&rOtherPoint);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Calculate Non-local damage around the tip
        double NonlocalDamage = 0.0;
        double WeightingFunctionDenominator = 0.0;
        int NNeighbours = TipNeighbours.size();

        #pragma omp parallel for reduction(+:NonlocalDamage,WeightingFunctionDenominator)
        for(int i = 0; i < NNeighbours; i++)
        {
            const FracturePoint& MyPoint = *(TipNeighbours[i]);
            NonlocalDamage += (MyPoint.Weight)*
                                exp(-4.0*(MyPoint.TipDistance)*(MyPoint.TipDistance)/(PropagationLength*PropagationLength))*
                                (MyPoint.Damage);
            WeightingFunctionDenominator += (MyPoint.Weight)*exp(-4.0*(MyPoint.TipDistance)*(MyPoint.TipDistance)/(PropagationLength*PropagationLength));
        }

        if(WeightingFunctionDenominator > 1.0e-20)
            NonlocalDamage = NonlocalDamage/WeightingFunctionDenominator;
        else
            NonlocalDamage = 0.0;

        // Check fracture propagation
        if(NonlocalDamage >= rParameters["fracture_data"]["propagation_damage"].GetDouble())
            this->PropagateFracture(itFracture,rPropagationData,AuxPropagationVariables,rParameters);
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeCheckFracture(
        const PropagationGlobalVariables& PropagationData,
        Parameters& rParameters,
        ModelPart& rModelPart,
        const bool& move_mesh_flag)
    {
        if(PropagationData.PropagateFractures==false)
        {
            // Move mesh to the current position
            if(move_mesh_flag==true)
                this->UpdateMeshPosition(rModelPart);
        }
        else
        {
            // Write PropagationData.tcl file
            this->WritePropagationData(PropagationData,rParameters);
        }
    }

    ///------------------------------------------------------------------------------------

    void WritePropagationData(
        const PropagationGlobalVariables& PropagationData,
        Parameters& rParameters)
    {

        std::fstream PropDataFile;
        PropDataFile.open ("PropagationData.tcl", std::fstream::out | std::fstream::trunc);
        PropDataFile.precision(12);

        PropDataFile << "set PropagationData [list]" << std::endl;

        // GiDPath
        PropDataFile << "lappend PropagationData \"" << rParameters["fracture_data"]["gid_path"].GetString() << "\"" << std::endl;

        int Id;

        // BodyVolumesDict
        for(unsigned int i = 0; i < rParameters["body_volumes_list"].size(); i++)
        {
            Id = rParameters["body_volumes_list"][i]["id"].GetInt();

            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["body_volumes_list"][i]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["body_volumes_list"][i]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set BodyVolumesDict " << Id << " Groups $Groups" << std::endl;
            PropDataFile << "set Surfaces [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["body_volumes_list"][i]["surfaces"].size(); j++)
            {
                PropDataFile << "lappend Surfaces " << rParameters["body_volumes_list"][i]["surfaces"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set BodyVolumesDict " << Id << " Surfaces $Surfaces" << std::endl;
            PropDataFile << "dict set BodyVolumesDict " << Id << " MeshSize " << rParameters["body_volumes_list"][i]["mesh_size"].GetDouble() << std::endl;
        }
        PropDataFile << "lappend PropagationData $BodyVolumesDict" << std::endl;

        // FracturesDict
        for(unsigned int i = 0; i < rParameters["fractures_list"].size(); i++)
        {
            Id = rParameters["fractures_list"][i]["id"].GetInt();

            // TipPoint
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Id "
                         << rParameters["fractures_list"][i]["tip_point"]["id"].GetInt() << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["tip_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["tip_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["tip_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Coordinates $Coordinates" << std::endl;
            // TopPoint
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Id "
                         << rParameters["fractures_list"][i]["top_point"]["id"].GetInt() << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["top_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["top_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["top_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Coordinates $Coordinates" << std::endl;
            // BotPoint
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Id "
                         << rParameters["fractures_list"][i]["bot_point"]["id"].GetInt() << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["bot_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["bot_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["bot_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Coordinates $Coordinates" << std::endl;
            // LeftPoint
            PropDataFile << "dict set FracturesDict " << Id << " LeftPoint Id "
                         << rParameters["fractures_list"][i]["left_point"]["id"].GetInt() << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["left_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["left_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["left_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " LeftPoint Coordinates $Coordinates" << std::endl;
            // RightPoint
            PropDataFile << "dict set FracturesDict " << Id << " RightPoint Id "
                         << rParameters["fractures_list"][i]["right_point"]["id"].GetInt() << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["right_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["right_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["right_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " RightPoint Coordinates $Coordinates" << std::endl;
            // TopLine
            PropDataFile << "dict set FracturesDict " << Id << " TopLine Id "
                         << rParameters["fractures_list"][i]["top_line"]["id"].GetInt() << std::endl;
            // BotLine
            PropDataFile << "dict set FracturesDict " << Id << " BotLine Id "
                         << rParameters["fractures_list"][i]["bot_line"]["id"].GetInt() << std::endl;
            // LeftLine
            PropDataFile << "dict set FracturesDict " << Id << " LeftLine Id "
                         << rParameters["fractures_list"][i]["left_line"]["id"].GetInt() << std::endl;
            // RightLine
            PropDataFile << "dict set FracturesDict " << Id << " RightLine Id "
                         << rParameters["fractures_list"][i]["right_line"]["id"].GetInt() << std::endl;
            // TopLeftSurface
            PropDataFile << "dict set FracturesDict " << Id << " TopLeftSurface Id "
                         << rParameters["fractures_list"][i]["top_left_surface"]["id"].GetInt() << std::endl;
            PropDataFile << "set Lines [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["top_left_surface"]["lines"].size(); j++)
            {
                PropDataFile << "lappend Lines " << rParameters["fractures_list"][i]["top_left_surface"]["lines"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " TopLeftSurface Lines $Lines" << std::endl;
            // TopRightSurface
            PropDataFile << "dict set FracturesDict " << Id << " TopRightSurface Id "
                         << rParameters["fractures_list"][i]["top_right_surface"]["id"].GetInt() << std::endl;
            PropDataFile << "set Lines [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["top_right_surface"]["lines"].size(); j++)
            {
                PropDataFile << "lappend Lines " << rParameters["fractures_list"][i]["top_right_surface"]["lines"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " TopRightSurface Lines $Lines" << std::endl;
            // BotLeftSurface
            PropDataFile << "dict set FracturesDict " << Id << " BotLeftSurface Id "
                         << rParameters["fractures_list"][i]["bot_left_surface"]["id"].GetInt() << std::endl;
            PropDataFile << "set Lines [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["bot_left_surface"]["lines"].size(); j++)
            {
                PropDataFile << "lappend Lines " << rParameters["fractures_list"][i]["bot_left_surface"]["lines"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BotLeftSurface Lines $Lines" << std::endl;
            // BotRightSurface
            PropDataFile << "dict set FracturesDict " << Id << " BotRightSurface Id "
                         << rParameters["fractures_list"][i]["bot_right_surface"]["id"].GetInt() << std::endl;
            PropDataFile << "set Lines [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["bot_right_surface"]["lines"].size(); j++)
            {
                PropDataFile << "lappend Lines " << rParameters["fractures_list"][i]["bot_right_surface"]["lines"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BotRightSurface Lines $Lines" << std::endl;
            // LeftInterfaceVolume
            PropDataFile << "dict set FracturesDict " << Id << " LeftInterfaceVolume Id "
                         << rParameters["fractures_list"][i]["left_interface_volume"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " LeftInterfaceVolume Layer \""
                         << rParameters["fractures_list"][i]["left_interface_volume"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["left_interface_volume"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["left_interface_volume"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " LeftInterfaceVolume Groups $Groups" << std::endl;
            // RightInterfaceVolume
            PropDataFile << "dict set FracturesDict " << Id << " RightInterfaceVolume Id "
                         << rParameters["fractures_list"][i]["right_interface_volume"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " RightInterfaceVolume Layer \""
                         << rParameters["fractures_list"][i]["right_interface_volume"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["right_interface_volume"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["right_interface_volume"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " RightInterfaceVolume Groups $Groups" << std::endl;
            // BodyVolumes
            PropDataFile << "set BodyVolumes [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["body_volumes"].size(); j++)
            {
                PropDataFile << "lappend BodyVolumes " << rParameters["fractures_list"][i]["body_volumes"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BodyVolumes $BodyVolumes" << std::endl;
        }
        PropDataFile << "lappend PropagationData $FracturesDict" << std::endl;

        // PropagationDict
        PropDataFile << "set PropagationDict [dict create]" << std::endl;
        for(unsigned int i = 0; i < PropagationData.PropagationVector.size(); i++)
        {
            PropDataFile << "dict set PropagationDict " << i << " MotherFractureId "
                         << PropagationData.PropagationVector[i].MotherFractureId << std::endl;

            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].TopInitCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].TopInitCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].TopInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " TopInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].BotInitCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].BotInitCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].BotInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " BotInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].LeftInitCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].LeftInitCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].LeftInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " LeftInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].RightInitCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].RightInitCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].RightInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " RightInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].TopEndCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].TopEndCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].TopEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " TopEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].BotEndCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].BotEndCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].BotEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " BotEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].LeftEndCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].LeftEndCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].LeftEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " LeftEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].RightEndCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].RightEndCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].RightEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " RightEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.PropagationVector[i].TipCoordinates[0]
                         << " " << PropagationData.PropagationVector[i].TipCoordinates[1] << " "
                         << PropagationData.PropagationVector[i].TipCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set PropagationDict " << i << " TipCoordinates $Coordinates" << std::endl;
        }
        PropDataFile << "lappend PropagationData $PropagationDict" << std::endl;

        // BifurcationDict
        PropDataFile << "set BifurcationDict [dict create]" << std::endl;
        for(unsigned int i = 0; i < PropagationData.BifurcationVector.size(); i++)
        {
            PropDataFile << "dict set BifurcationDict " << i << " MotherFractureId "
                         << PropagationData.BifurcationVector[i].MotherFractureId << std::endl;

            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopInitCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopInitCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotInitCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotInitCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].LeftInitCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].LeftInitCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].LeftInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " LeftInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].RightInitCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].RightInitCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].RightInitCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " RightInitCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopTopEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopTopEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopTopEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopTopEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopBotEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopBotEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopBotEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopBotEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopLeftEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopLeftEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopLeftEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopLeftEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopRightEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopRightEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopRightEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopRightEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopTipCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopTipCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopTipCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopTipCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotTopEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotTopEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotTopEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotTopEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotBotEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotBotEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotBotEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotBotEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotLeftEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotLeftEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotLeftEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotLeftEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotRightEndCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotRightEndCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotRightEndCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotRightEndCoordinates $Coordinates" << std::endl;
            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].BotTipCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].BotTipCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].BotTipCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " BotTipCoordinates $Coordinates" << std::endl;
        }
        PropDataFile << "lappend PropagationData $BifurcationDict" << std::endl;

        PropDataFile << "return $PropagationData" << std::endl;

        PropDataFile.close();
    }

/// Mapping -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeMapping(
        UtilityVariables& rAuxVariables,
        ModelPart& rModelPartNew,
        const bool& move_mesh_flag)
    {
        // Move mesh to the original position to work in the reference state
        if(move_mesh_flag==true)
            this->ResetMeshPosition(rModelPartNew);

        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPartNew);
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void NodalVariablesMapping(
        const UtilityVariables& AuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPartOld,
        ModelPart& rModelPartNew)
    {
        // Define ElementOld Cell matrix
        std::vector< std::vector< std::vector< std::vector<Element::Pointer> > > > ElementOldCellMatrix;
        ElementOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++)
            ElementOldCellMatrix[i].resize(AuxVariables.NColumns);
        for(int i = 0; i < AuxVariables.NRows; i++)
        {
            for(int j = 0; j < AuxVariables.NColumns; j++)
            {
                ElementOldCellMatrix[i][j].resize(AuxVariables.NSections);
            }
        }

        // Locate Old Elments in cells
        double X_me;
        double Y_me;
        double Z_me;
        int PointsNumber;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int m = 0; m < NumBodySubModelParts; m++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][m].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(X_me,Y_me,Z_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                double Z_back = itElemOld->GetGeometry().GetPoint(0).Z0();
                double Z_front = Z_back;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();
                    Z_me = itElemOld->GetGeometry().GetPoint(j).Z0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                    if(Z_me > Z_front) Z_front = Z_me;
                    else if(Z_me < Z_back) Z_back = Z_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);
                int Section_back = int((Z_back-AuxVariables.Z_min)/AuxVariables.SectionSize);
                int Section_front = int((Z_front-AuxVariables.Z_min)/AuxVariables.SectionSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;
                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;
                if(Section_back < 0) Section_back = 0;
                else if(Section_back >= AuxVariables.NSections) Section_back = AuxVariables.NSections-1;
                if(Section_front >= AuxVariables.NSections) Section_front = AuxVariables.NSections-1;
                else if(Section_front < 0) Section_front = 0;

                for(int s = Section_back; s <= Section_front; s++)
                {
                    for(int k = Row_top; k <= Row_bot; k++)
                    {
                        for(int l = Column_left; l <= Column_right; l++)
                        {
                            #pragma omp critical
                            {
                                ElementOldCellMatrix[k][l][s].push_back((*(itElemOld.base())));
                            }
                        }
                    }
                }
            }
        }

        unsigned int NumInterfaceSubModelPartsOld = rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"].size();

        // Loop through all InterfaceSubModelParts
        for(unsigned int m = 0; m < NumInterfaceSubModelPartsOld; m++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"][m].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(X_me,Y_me,Z_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                double Z_back = itElemOld->GetGeometry().GetPoint(0).Z0();
                double Z_front = Z_back;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();
                    Z_me = itElemOld->GetGeometry().GetPoint(j).Z0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                    if(Z_me > Z_front) Z_front = Z_me;
                    else if(Z_me < Z_back) Z_back = Z_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);
                int Section_back = int((Z_back-AuxVariables.Z_min)/AuxVariables.SectionSize);
                int Section_front = int((Z_front-AuxVariables.Z_min)/AuxVariables.SectionSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;
                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;
                if(Section_back < 0) Section_back = 0;
                else if(Section_back >= AuxVariables.NSections) Section_back = AuxVariables.NSections-1;
                if(Section_front >= AuxVariables.NSections) Section_front = AuxVariables.NSections-1;
                else if(Section_front < 0) Section_front = 0;

                for(int s = Section_back; s <= Section_front; s++)
                {
                    for(int k = Row_top; k <= Row_bot; k++)
                    {
                        for(int l = Column_left; l <= Column_right; l++)
                        {
                            #pragma omp critical
                            {
                                ElementOldCellMatrix[k][l][s].push_back((*(itElemOld.base())));
                            }
                        }
                    }
                }
            }
        }

        // Locate new nodes inside old elements and interpolate nodal variables
        const int NNodes = static_cast<int>(rModelPartNew.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPartNew.NodesBegin();

        array_1d<double,3> GlobalCoordinates;
        array_1d<double,3> LocalCoordinates;

        #pragma omp parallel for private(X_me,Y_me,Z_me,PointsNumber,GlobalCoordinates,LocalCoordinates)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNodeNew = node_begin + i;

            X_me = itNodeNew->X0();
            Y_me = itNodeNew->Y0();
            Z_me = itNodeNew->Z0();

            int Column = int((X_me-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Row = int((AuxVariables.Y_max-Y_me)/AuxVariables.RowSize);
            int Section = int((Z_me-AuxVariables.Z_min)/AuxVariables.SectionSize);

            if(Column >= AuxVariables.NColumns) Column = AuxVariables.NColumns-1;
            else if(Column < 0) Column = 0;
            if(Row >= AuxVariables.NRows) Row = AuxVariables.NRows-1;
            else if(Row < 0) Row = 0;
            if(Section >= AuxVariables.NSections) Section = AuxVariables.NSections-1;
            else if(Section < 0) Section = 0;

            noalias(GlobalCoordinates) = itNodeNew->Coordinates(); //Coordinates of new nodes are still in the original position
            bool IsInside = false;
            Element::Pointer pElementOld;

            for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column][Section]).size(); m++)
            {
                pElementOld = ElementOldCellMatrix[Row][Column][Section][m];
                IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column][Section]).size(); m++)
                {
                    pElementOld = ElementOldCellMatrix[Row][Column][Section][m];
                    IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
                std::cout << "ERROR!!, NONE OF THE OLD ELEMENTS CONTAINS NODE: " << itNodeNew->Id() << std::endl;

            PointsNumber = pElementOld->GetGeometry().PointsNumber();
            Vector ShapeFunctionsValuesVector(PointsNumber);
            Vector NodalVariableVector(PointsNumber);

            pElementOld->GetGeometry().ShapeFunctionsValues(ShapeFunctionsValuesVector,LocalCoordinates);

            // Interpolation of nodal variables
            if( itNodeNew->IsFixed(DISPLACEMENT_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DISPLACEMENT_X);
                }
                itNodeNew->FastGetSolutionStepValue(DISPLACEMENT_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(VELOCITY_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(VELOCITY_X);
                }
                itNodeNew->FastGetSolutionStepValue(VELOCITY_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(ACCELERATION_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(ACCELERATION_X);
                }
                itNodeNew->FastGetSolutionStepValue(ACCELERATION_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(DISPLACEMENT_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DISPLACEMENT_Y);
                }
                itNodeNew->FastGetSolutionStepValue(DISPLACEMENT_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(VELOCITY_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(VELOCITY_Y);
                }
                itNodeNew->FastGetSolutionStepValue(VELOCITY_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(ACCELERATION_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(ACCELERATION_Y);
                }
                itNodeNew->FastGetSolutionStepValue(ACCELERATION_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(DISPLACEMENT_Z)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DISPLACEMENT_Z);
                }
                itNodeNew->FastGetSolutionStepValue(DISPLACEMENT_Z) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(VELOCITY_Z)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(VELOCITY_Z);
                }
                itNodeNew->FastGetSolutionStepValue(VELOCITY_Z) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(ACCELERATION_Z)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(ACCELERATION_Z);
                }
                itNodeNew->FastGetSolutionStepValue(ACCELERATION_Z) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(LIQUID_PRESSURE)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(LIQUID_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(LIQUID_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(DT_LIQUID_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GaussPointStateVariableMapping(
        const UtilityVariables& AuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPartOld,
        ModelPart& rModelPartNew)
    {
        // Define GaussPointOld Cell matrix
        std::vector< std::vector< std::vector< std::vector<GaussPointOld> > > > BodyGaussPointOldCellMatrix;
        BodyGaussPointOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++)
            BodyGaussPointOldCellMatrix[i].resize(AuxVariables.NColumns);
        for(int i = 0; i < AuxVariables.NRows; i++)
        {
            for(int j = 0; j < AuxVariables.NColumns; j++)
            {
                BodyGaussPointOldCellMatrix[i][j].resize(AuxVariables.NSections);
            }
        }

        std::vector< std::vector< std::vector< std::vector<GaussPointOld> > > > InterfaceGaussPointOldCellMatrix;
        InterfaceGaussPointOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++)
            InterfaceGaussPointOldCellMatrix[i].resize(AuxVariables.NColumns);
        for(int i = 0; i < AuxVariables.NRows; i++)
        {
            for(int j = 0; j < AuxVariables.NColumns; j++)
            {
                InterfaceGaussPointOldCellMatrix[i][j].resize(AuxVariables.NSections);
            }
        }

        // Locate Old Gauss Points in cells
        GaussPointOld MyGaussPointOld;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfoOld = rModelPartOld.GetProcessInfo();
        array_1d<double,3> AuxLocalCoordinates;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();

        // Loop through all OLD BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(MyGaussPointOld,MyIntegrationMethod,AuxLocalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> StateVariableVector(NumGPoints);
                itElem->CalculateOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoOld);
                int Row;
                int Column;
                int Section;

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPointOld Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyGaussPointOld.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // GaussPointOld Weight
                    MyGaussPointOld.Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // GaussPointOld StateVariable and Damage
                    MyGaussPointOld.StateVariable = StateVariableVector[GPoint];

                    // GaussPointOld Row and Column
                    Row = int((AuxVariables.Y_max-MyGaussPointOld.Coordinates[1])/AuxVariables.RowSize);
                    Column = int((MyGaussPointOld.Coordinates[0]-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    Section = int((MyGaussPointOld.Coordinates[2]-AuxVariables.Z_min)/AuxVariables.SectionSize);
                    #pragma omp critical
                    {
                        BodyGaussPointOldCellMatrix[Row][Column][Section].push_back(MyGaussPointOld);
                    }
                }
            }
        }

        unsigned int NumInterfaceSubModelPartsOld = rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"].size();

        // Loop through all OLD InterfaceSubModelParts
        for(unsigned int i = 0; i < NumInterfaceSubModelPartsOld; i++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"][i].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(MyGaussPointOld,MyIntegrationMethod,AuxLocalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> StateVariableVector(NumGPoints);
                itElem->CalculateOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoOld);
                int Row;
                int Column;
                int Section;

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPointOld Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyGaussPointOld.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // GaussPointOld Weight
                    MyGaussPointOld.Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // GaussPointOld StateVariable
                    MyGaussPointOld.StateVariable = StateVariableVector[GPoint];

                    // GaussPointOld Row and Column
                    Row = int((AuxVariables.Y_max-MyGaussPointOld.Coordinates[1])/AuxVariables.RowSize);
                    Column = int((MyGaussPointOld.Coordinates[0]-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    Section = int((MyGaussPointOld.Coordinates[2]-AuxVariables.Z_min)/AuxVariables.SectionSize);
                    #pragma omp critical
                    {
                        InterfaceGaussPointOldCellMatrix[Row][Column][Section].push_back(MyGaussPointOld);
                    }
                }
            }
        }

        // Transfer state variables from old Gauss points to new Gauss Points (nonlocal average)
        const ProcessInfo& CurrentProcessInfoNew = rModelPartNew.GetProcessInfo();
        const double PropagationLength = rParameters["fracture_data"]["propagation_length"].GetDouble();
        array_1d<double,3> AuxGlobalCoordinates;

        // Loop through all NEW BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& SubModelPart = rModelPartNew.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            double DamageThreshold = el_begin->GetProperties()[DAMAGE_THRESHOLD];

            #pragma omp parallel for private(MyIntegrationMethod,AuxLocalCoordinates,AuxGlobalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                std::vector<double> StateVariableVector(NumGPoints);

                // Loop through NEW GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPointNew Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(AuxGlobalCoordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates
                    double X_me = AuxGlobalCoordinates[0];
                    double Y_me = AuxGlobalCoordinates[1];
                    double Z_me = AuxGlobalCoordinates[2];

                    // GaussPointNew search area
                    double X_left = X_me - PropagationLength;
                    double X_right = X_me + PropagationLength;
                    double Y_top = Y_me + PropagationLength;
                    double Y_bot = Y_me - PropagationLength;
                    double Z_back = Z_me - PropagationLength;
                    double Z_front = Z_me + PropagationLength;

                    int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                    int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);
                    int Section_back = int((Z_back-AuxVariables.Z_min)/AuxVariables.SectionSize);
                    int Section_front = int((Z_front-AuxVariables.Z_min)/AuxVariables.SectionSize);

                    if(Column_left < 0) Column_left = 0;
                    if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                    if(Row_top < 0) Row_top = 0;
                    if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                    if(Section_back < 0) Section_back = 0;
                    if(Section_front >= AuxVariables.NSections) Section_front = AuxVariables.NSections-1;

                    // Search GaussPointOld neighbours around the GaussPointNew and compute nonlocal state variable
                    double Numerator = 0.0;
                    double WeightingFunctionDenominator = 0.0;
                    double Distance;
                    for(int s = Section_back; s <= Section_front; s++)
                    {
                        for(int k = Row_top; k <= Row_bot; k++)
                        {
                            for(int l = Column_left; l<= Column_right; l++)
                            {
                                for(unsigned int m = 0; m < BodyGaussPointOldCellMatrix[k][l][s].size(); m++)
                                {
                                    GaussPointOld& rOtherGaussPointOld = BodyGaussPointOldCellMatrix[k][l][s][m];

                                    Distance = sqrt((rOtherGaussPointOld.Coordinates[0]-X_me)*(rOtherGaussPointOld.Coordinates[0]-X_me) +
                                                    (rOtherGaussPointOld.Coordinates[1]-Y_me)*(rOtherGaussPointOld.Coordinates[1]-Y_me) +
                                                    (rOtherGaussPointOld.Coordinates[2]-Z_me)*(rOtherGaussPointOld.Coordinates[2]-Z_me));

                                    if(Distance <= PropagationLength)
                                    {
                                        Numerator += rOtherGaussPointOld.Weight
                                                    *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength))
                                                    *rOtherGaussPointOld.StateVariable;
                                        WeightingFunctionDenominator += rOtherGaussPointOld.Weight
                                                                        *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength));
                                    }
                                }
                            }
                        }
                    }
                    // Save computed stateVariable
                    if(WeightingFunctionDenominator > 0.0)
                        StateVariableVector[GPoint] = Numerator/WeightingFunctionDenominator;
                    else
                        StateVariableVector[GPoint] = DamageThreshold;
                }
                // Set stateVariable of new GaussPoints
                itElem->SetValuesOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoNew);
            }
        }

        unsigned int NumInterfaceSubModelParts = rParameters["fracture_data"]["interface_domain_sub_model_part_list"].size();
        const double PropagationDamage = rParameters["fracture_data"]["propagation_damage"].GetDouble();

        // Loop through all NEW InterfaceSubModelParts
        for(unsigned int i = 0; i < NumInterfaceSubModelParts; i++)
        {
            ModelPart& SubModelPart = rModelPartNew.GetSubModelPart(rParameters["fracture_data"]["interface_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            //double DamageThreshold = el_begin->GetProperties()[DAMAGE_THRESHOLD];

            #pragma omp parallel for private(MyIntegrationMethod,AuxLocalCoordinates,AuxGlobalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                std::vector<double> StateVariableVector(NumGPoints);

                // Loop through NEW GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // GaussPointNew Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(AuxGlobalCoordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates
                    double X_me = AuxGlobalCoordinates[0];
                    double Y_me = AuxGlobalCoordinates[1];
                    double Z_me = AuxGlobalCoordinates[2];

                    // GaussPointNew search area
                    double X_left = X_me - PropagationLength;
                    double X_right = X_me + PropagationLength;
                    double Y_top = Y_me + PropagationLength;
                    double Y_bot = Y_me - PropagationLength;
                    double Z_back = Z_me - PropagationLength;
                    double Z_front = Z_me + PropagationLength;

                    int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                    int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);
                    int Section_back = int((Z_back-AuxVariables.Z_min)/AuxVariables.SectionSize);
                    int Section_front = int((Z_front-AuxVariables.Z_min)/AuxVariables.SectionSize);

                    if(Column_left < 0) Column_left = 0;
                    if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                    if(Row_top < 0) Row_top = 0;
                    if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                    if(Section_back < 0) Section_back = 0;
                    if(Section_front >= AuxVariables.NSections) Section_front = AuxVariables.NSections-1;

                    // Search GaussPointOld neighbours around the GaussPointNew and compute nonlocal state variable
                    double Numerator = 0.0;
                    double WeightingFunctionDenominator = 0.0;
                    double Distance;
                    for(int s = Section_back; s <= Section_front; s++)
                    {
                        for(int k = Row_top; k <= Row_bot; k++)
                        {
                            for(int l = Column_left; l<= Column_right; l++)
                            {
                                for(unsigned int m = 0; m < InterfaceGaussPointOldCellMatrix[k][l][s].size(); m++)
                                {
                                    GaussPointOld& rOtherGaussPointOld = InterfaceGaussPointOldCellMatrix[k][l][s][m];

                                    Distance = sqrt((rOtherGaussPointOld.Coordinates[0]-X_me)*(rOtherGaussPointOld.Coordinates[0]-X_me) +
                                                    (rOtherGaussPointOld.Coordinates[1]-Y_me)*(rOtherGaussPointOld.Coordinates[1]-Y_me) +
                                                    (rOtherGaussPointOld.Coordinates[2]-Z_me)*(rOtherGaussPointOld.Coordinates[2]-Z_me));

                                    if(Distance <= PropagationLength)
                                    {
                                        Numerator += rOtherGaussPointOld.Weight
                                                    *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength))
                                                    *rOtherGaussPointOld.StateVariable;
                                        WeightingFunctionDenominator += rOtherGaussPointOld.Weight
                                                                        *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength));
                                    }
                                }
                            }
                        }
                    }

                    // Save computed stateVariable
                    if(WeightingFunctionDenominator > 0.0)
                    {
                        StateVariableVector[GPoint] = Numerator/WeightingFunctionDenominator;
                        if(StateVariableVector[GPoint] < PropagationDamage)
                        {
                            StateVariableVector[GPoint] = PropagationDamage;
                        }
                    }
                    else
                        StateVariableVector[GPoint] = PropagationDamage;
                }
                // Set stateVariable of new GaussPoints
                itElem->SetValuesOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoNew);
            }
        }
    }

/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ResetMeshPosition(ModelPart& rModelPart)
    {
        // Move mesh to the original position
        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            itNode->X() = itNode->X0();
            itNode->Y() = itNode->Y0();
            itNode->Z() = itNode->Z0();
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void UpdateMeshPosition(ModelPart& rModelPart)
    {
        // Move mesh to the current position
        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            itNode->X() = itNode->X0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_X);
            itNode->Y() = itNode->Y0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_Y);
            itNode->Z() = itNode->Z0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }

private:

/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeCellMatrixDimensions(
        UtilityVariables& rAuxVariables,
        ModelPart& rModelPart)
    {
        // Compute X, Y and Z limits of the current geometry
        unsigned int NumThreads = ParallelUtilities::GetNumThreads();
        std::vector<double> X_max_partition(NumThreads);
        std::vector<double> X_min_partition(NumThreads);
        std::vector<double> Y_max_partition(NumThreads);
        std::vector<double> Y_min_partition(NumThreads);
        std::vector<double> Z_max_partition(NumThreads);
        std::vector<double> Z_min_partition(NumThreads);

        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            X_max_partition[k] = node_begin->X();
            X_min_partition[k] = X_max_partition[k];
            Y_max_partition[k] = node_begin->Y();
            Y_min_partition[k] = Y_max_partition[k];
            Z_max_partition[k] = node_begin->Z();
            Z_min_partition[k] = Z_max_partition[k];

            double X_me, Y_me, Z_me;

            #pragma omp for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                X_me = itNode->X();
                Y_me = itNode->Y();
                Z_me = itNode->Z();

                if( X_me > X_max_partition[k] ) X_max_partition[k] = X_me;
                else if( X_me < X_min_partition[k] ) X_min_partition[k] = X_me;

                if( Y_me > Y_max_partition[k] ) Y_max_partition[k] = Y_me;
                else if( Y_me < Y_min_partition[k] ) Y_min_partition[k] = Y_me;

                if( Z_me > Z_max_partition[k] ) Z_max_partition[k] = Z_me;
                else if( Z_me < Z_min_partition[k] ) Z_min_partition[k] = Z_me;
            }
        }

        rAuxVariables.X_max = X_max_partition[0];
        rAuxVariables.X_min = X_min_partition[0];
        rAuxVariables.Y_max = Y_max_partition[0];
        rAuxVariables.Y_min = Y_min_partition[0];
        rAuxVariables.Z_max = Z_max_partition[0];
        rAuxVariables.Z_min = Z_min_partition[0];

        for(unsigned int i=1; i < NumThreads; i++)
        {
            if(X_max_partition[i] > rAuxVariables.X_max) rAuxVariables.X_max = X_max_partition[i];
            if(X_min_partition[i] < rAuxVariables.X_min) rAuxVariables.X_min = X_min_partition[i];
            if(Y_max_partition[i] > rAuxVariables.Y_max) rAuxVariables.Y_max = Y_max_partition[i];
            if(Y_min_partition[i] < rAuxVariables.Y_min) rAuxVariables.Y_min = Y_min_partition[i];
            if(Z_max_partition[i] > rAuxVariables.Z_max) rAuxVariables.Z_max = Z_max_partition[i];
            if(Z_min_partition[i] < rAuxVariables.Z_min) rAuxVariables.Z_min = Z_min_partition[i];
        }

        // Calculate Average Element Length
        double AverageElementLength = 0.0;

        int NElems = static_cast<int>(rModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for reduction(+:AverageElementLength)
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;

            AverageElementLength += itElem->GetGeometry().Length();
        }

        AverageElementLength = AverageElementLength/NElems;

        // Compute FracturePoints CellMatrix dimensions
        rAuxVariables.NRows = int((rAuxVariables.Y_max-rAuxVariables.Y_min)/AverageElementLength);
        rAuxVariables.NColumns = int((rAuxVariables.X_max-rAuxVariables.X_min)/AverageElementLength);
        rAuxVariables.NSections = int((rAuxVariables.Z_max-rAuxVariables.Z_min)/AverageElementLength);
        if(rAuxVariables.NRows <= 0) rAuxVariables.NRows = 1;
        if(rAuxVariables.NColumns <= 0) rAuxVariables.NColumns = 1;
        if(rAuxVariables.NSections <= 0) rAuxVariables.NSections = 1;
        rAuxVariables.RowSize = (rAuxVariables.Y_max-rAuxVariables.Y_min)/rAuxVariables.NRows;
        rAuxVariables.ColumnSize = (rAuxVariables.X_max-rAuxVariables.X_min)/rAuxVariables.NColumns;
        rAuxVariables.SectionSize = (rAuxVariables.Z_max-rAuxVariables.Z_min)/rAuxVariables.NSections;
    }

/// Fracture Propagation Check ------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetFracturePoints(
        PropagationGlobalVariables& rPropagationData,
        UtilityVariables& rAuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPart)
    {
        // Compute FracturePointsCellMatrix dimensions
        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPart);

        rPropagationData.FracturePointsCellMatrix.resize(rAuxVariables.NRows);
        for(int i = 0; i < rAuxVariables.NRows; i++)
            rPropagationData.FracturePointsCellMatrix[i].resize(rAuxVariables.NColumns);
        for(int i = 0; i < rAuxVariables.NRows; i++)
        {
            for(int j = 0; j < rAuxVariables.NColumns; j++)
            {
                rPropagationData.FracturePointsCellMatrix[i][j].resize(rAuxVariables.NSections);
            }
        }

        // Locate FracturePoints inside CellMatrix
        FracturePoint MyFracturePoint;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        rPropagationData.pProcessInfo = rModelPart.pGetProcessInfo();
        array_1d<double,3> AuxLocalCoordinates;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& BodySubModelPart = rModelPart.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(BodySubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = BodySubModelPart.ElementsBegin();

            // Loop through all body elements
            #pragma omp parallel for private(MyFracturePoint,MyIntegrationMethod,AuxLocalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> DamageVector(NumGPoints);
                itElem->CalculateOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,rCurrentProcessInfo);
                int Row;
                int Column;
                int Section;

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // FracturePoint Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyFracturePoint.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // FracturePoint Weight
                    MyFracturePoint.Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // FracturePoint Damage
                    MyFracturePoint.Damage = DamageVector[GPoint];

                    // FracturePoint Row, Column and Section
                    Row = int((rAuxVariables.Y_max-MyFracturePoint.Coordinates[1])/rAuxVariables.RowSize);
                    Column = int((MyFracturePoint.Coordinates[0]-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    Section = int((MyFracturePoint.Coordinates[2]-rAuxVariables.Z_min)/rAuxVariables.SectionSize);

                    // Element containing the FracturePoint
                    MyFracturePoint.pElement = (*(itElem.base()));

                    #pragma omp critical
                    {
                        rPropagationData.FracturePointsCellMatrix[Row][Column][Section].push_back(MyFracturePoint);
                    }
                }
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void PropagateFracture(
        const unsigned int& itFracture,
        PropagationGlobalVariables& rPropagationData,
        PropagationLocalVariables& rAuxPropagationVariables,
        Parameters& rParameters)
    {
        // Compute Propagation Coordinates
        double TipX = 0.0;
        double TipY = 0.0;
        double TipZ = 0.0;
        double TipDen = 0.0;

        #pragma omp parallel sections reduction(+:TipX,TipY,TipZ,TipDen)
        {
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                {
                    this->AverageTipCoordinates(TipX,TipY,TipZ,TipDen,*(rAuxPropagationVariables.TopFrontFracturePoints[i]));
                }
            }
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                {
                    this->AverageTipCoordinates(TipX,TipY,TipZ,TipDen,*(rAuxPropagationVariables.BotFrontFracturePoints[i]));
                }
            }
        }

        array_1d<double,3> AuxArray1;
        array_1d<double,3> AuxArray2;
        const double PropagationLength = rParameters["fracture_data"]["propagation_length"].GetDouble();
        const double PropagationWidth = rParameters["fracture_data"]["propagation_width"].GetDouble();
        const double PropagationHeight = rParameters["fracture_data"]["propagation_height"].GetDouble();
        int MotherFractureId = rParameters["fractures_list"][itFracture]["id"].GetInt();
        Propagation MyPropagation;

        MyPropagation.MotherFractureId = MotherFractureId;

        MyPropagation.TipCoordinates[0] = TipX/TipDen;
        MyPropagation.TipCoordinates[1] = TipY/TipDen;
        MyPropagation.TipCoordinates[2] = TipZ/TipDen;

        // LeftInitCoordinates
        noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
        AuxArray1[1] += 0.5*PropagationHeight;
        noalias(MyPropagation.LeftInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
        // RightInitCoordinates
        noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
        AuxArray1[1] -= 0.5*PropagationHeight;
        noalias(MyPropagation.RightInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
        // TopInitCoordinates
        noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
        AuxArray1[2] += 0.5*PropagationWidth;
        noalias(MyPropagation.TopInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
        // BotInitCoordinates
        noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
        AuxArray1[2] -= 0.5*PropagationWidth;
        noalias(MyPropagation.BotInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

        // Check straight propagation
        const double CorrectionTol = rParameters["fracture_data"]["correction_tolerance"].GetDouble();
        noalias(AuxArray2) = MyPropagation.TipCoordinates;
        noalias(AuxArray1) = prod(rAuxPropagationVariables.RotationMatrix,AuxArray2);
        AuxArray1[1] = rAuxPropagationVariables.TipLocalCoordinates[1];
        AuxArray1[2] = rAuxPropagationVariables.TipLocalCoordinates[2];
        noalias(AuxArray2) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

        double Distance = sqrt((MyPropagation.TipCoordinates[0]-AuxArray2[0])*(MyPropagation.TipCoordinates[0]-AuxArray2[0])+
                               (MyPropagation.TipCoordinates[1]-AuxArray2[1])*(MyPropagation.TipCoordinates[1]-AuxArray2[1])+
                               (MyPropagation.TipCoordinates[2]-AuxArray2[2])*(MyPropagation.TipCoordinates[2]-AuxArray2[2]));
        if (Distance <= PropagationLength*CorrectionTol)
            noalias(MyPropagation.TipCoordinates) = AuxArray2;

        // Check whether new tip falls inside a valid element
        array_1d<double,3> LocalCoordinates;
        Element::Pointer pElement;
        bool IsInside = false;

        for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
        {
            pElement = rAuxPropagationVariables.TopFrontFracturePoints[i]->pElement;
            IsInside = pElement->GetGeometry().IsInside(MyPropagation.TipCoordinates,LocalCoordinates);
            if(IsInside) break;
        }
        if(IsInside == false)
        {
            for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
            {
                pElement = rAuxPropagationVariables.BotFrontFracturePoints[i]->pElement;
                IsInside = pElement->GetGeometry().IsInside(MyPropagation.TipCoordinates,LocalCoordinates);
                if(IsInside) break;
            }
        }

        double PropagationDistance = 0.1*PropagationLength;
        const double PropagationDamage = rParameters["fracture_data"]["propagation_damage"].GetDouble();
        const ProcessInfo& rCurrentProcessInfo = *(rPropagationData.pProcessInfo);

        if (IsInside == true)
        {
            std::vector<double> DamageVector;
            pElement->CalculateOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,rCurrentProcessInfo);
            unsigned int NumGPoints = DamageVector.size();
            double InvNumGPoints = 1.0/static_cast<double>(NumGPoints);
            double ElementDamage = 0.0;
            for (unsigned int i = 0; i < NumGPoints; i++)
            {
                ElementDamage += DamageVector[i];
            }
            ElementDamage *= InvNumGPoints;
            if (ElementDamage >= 0.5*PropagationDamage)
            {
                // Compute new Tip RotationMatrix
                this->CalculateNewTipRotationMatrix(rAuxPropagationVariables,MyPropagation.TipCoordinates,
                                                    MyPropagation.LeftInitCoordinates,PropagationDistance);
                // LeftEndCoordinates
                noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
                AuxArray1[1] += 0.5*PropagationHeight;
                noalias(MyPropagation.LeftEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
                // RightEndCoordinates
                noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
                AuxArray1[1] -= 0.5*PropagationHeight;
                noalias(MyPropagation.RightEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
                // TopEndCoordinates
                noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
                AuxArray1[2] += 0.5*PropagationWidth;
                noalias(MyPropagation.TopEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
                // BotEndCoordinates
                noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
                AuxArray1[2] -= 0.5*PropagationWidth;
                noalias(MyPropagation.BotEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

                rPropagationData.PropagationVector.push_back(MyPropagation);
                rPropagationData.PropagateFractures = true;
                return;
            }
        }

        // Being here means that the new tip does not fall inside a valid element. We need to check top and bot fracture points
        bool PropagateTop = false;
        bool PropagateBot = false;
        double TopEndX = 0.0;
        double TopEndY = 0.0;
        double TopEndZ = 0.0;
        double TopEndDen = 0.0;
        double BotEndX = 0.0;
        double BotEndY = 0.0;
        double BotEndZ = 0.0;
        double BotEndDen = 0.0;
        array_1d<double,3> GlobalCoordinates;

        #pragma omp parallel sections private(GlobalCoordinates,LocalCoordinates,pElement,IsInside)
        {
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                {
                    this->AverageTipCoordinates(TopEndX,TopEndY,TopEndZ,TopEndDen,*(rAuxPropagationVariables.TopFrontFracturePoints[i]));
                }
                GlobalCoordinates[0] = TopEndX/TopEndDen;
                GlobalCoordinates[1] = TopEndY/TopEndDen;
                GlobalCoordinates[2] = TopEndZ/TopEndDen;

                // Check whether new tip falls inside a valid element
                IsInside = false;
                for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                {
                    pElement = rAuxPropagationVariables.TopFrontFracturePoints[i]->pElement;
                    IsInside = pElement->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates);
                    if(IsInside) break;
                }

                if (IsInside == true)
                {
                    std::vector<double> DamageVector;
                    pElement->CalculateOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,rCurrentProcessInfo);
                    unsigned int NumGPoints = DamageVector.size();
                    double InvNumGPoints = 1.0/static_cast<double>(NumGPoints);
                    double ElementDamage = 0.0;
                    for (unsigned int i = 0; i < NumGPoints; i++)
                    {
                        ElementDamage += DamageVector[i];
                    }
                    ElementDamage *= InvNumGPoints;
                    if (ElementDamage >= PropagationDamage)
                        PropagateTop = true;
                }
            }
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                {
                    this->AverageTipCoordinates(BotEndX,BotEndY,BotEndZ,BotEndDen,*(rAuxPropagationVariables.BotFrontFracturePoints[i]));
                }
                GlobalCoordinates[0] = BotEndX/BotEndDen;
                GlobalCoordinates[1] = BotEndY/BotEndDen;
                GlobalCoordinates[2] = BotEndZ/BotEndDen;

                // Check whether new tip falls inside a valid element
                IsInside = false;
                for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                {
                    pElement = rAuxPropagationVariables.BotFrontFracturePoints[i]->pElement;
                    IsInside = pElement->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates);
                    if(IsInside) break;
                }

                if (IsInside == true)
                {
                    std::vector<double> DamageVector;
                    pElement->CalculateOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,rCurrentProcessInfo);
                    unsigned int NumGPoints = DamageVector.size();
                    double InvNumGPoints = 1.0/static_cast<double>(NumGPoints);
                    double ElementDamage = 0.0;
                    for (unsigned int i = 0; i < NumGPoints; i++)
                    {
                        ElementDamage += DamageVector[i];
                    }
                    ElementDamage *= InvNumGPoints;
                    if (ElementDamage >= PropagationDamage)
                        PropagateBot = true;
                }
            }
        }

        if (PropagateTop == true && PropagateBot == true) // Bifurcation
        {
            Bifurcation MyBifurcation;
            MyBifurcation.MotherFractureId = MotherFractureId;

            MyBifurcation.TopTipCoordinates[0] = TopEndX/TopEndDen;
            MyBifurcation.TopTipCoordinates[1] = TopEndY/TopEndDen;
            MyBifurcation.TopTipCoordinates[2] = TopEndZ/TopEndDen;

            MyBifurcation.BotTipCoordinates[0] = BotEndX/BotEndDen;
            MyBifurcation.BotTipCoordinates[1] = BotEndY/BotEndDen;
            MyBifurcation.BotTipCoordinates[2] = BotEndZ/BotEndDen;

            // LeftInitCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationDistance;
            AuxArray1[1] += 0.5*PropagationHeight;
            noalias(MyBifurcation.LeftInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // RightInitCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationDistance;
            AuxArray1[1] -= 0.5*PropagationHeight;
            noalias(MyBifurcation.RightInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopInitCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationDistance;
            AuxArray1[2] += 0.5*PropagationWidth;
            noalias(MyBifurcation.TopInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotInitCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationDistance;
            AuxArray1[2] -= 0.5*PropagationWidth;
            noalias(MyBifurcation.BotInitCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

            // Compute new TopTip RotationMatrix
            this->CalculateNewTipRotationMatrix(rAuxPropagationVariables,MyBifurcation.TopTipCoordinates,
                                                MyBifurcation.LeftInitCoordinates,PropagationDistance);
            // TopLeftEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] += 0.5*PropagationHeight;
            noalias(MyBifurcation.TopLeftEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopRightEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] -= 0.5*PropagationHeight;
            noalias(MyBifurcation.TopRightEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopTopEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] += 0.5*PropagationWidth;
            noalias(MyBifurcation.TopTopEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopBotEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] -= 0.5*PropagationWidth;
            noalias(MyBifurcation.TopBotEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

            // Compute new BotTip RotationMatrix
            this->CalculateNewTipRotationMatrix(rAuxPropagationVariables,MyBifurcation.BotTipCoordinates,
                                                MyBifurcation.LeftInitCoordinates,PropagationDistance);
            // BotLeftEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] += 0.5*PropagationHeight;
            noalias(MyBifurcation.BotLeftEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotRightEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] -= 0.5*PropagationHeight;
            noalias(MyBifurcation.BotRightEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotTopEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] += 0.5*PropagationWidth;
            noalias(MyBifurcation.BotTopEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotBotEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] -= 0.5*PropagationWidth;
            noalias(MyBifurcation.BotBotEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

            rPropagationData.BifurcationVector.push_back(MyBifurcation);
            rPropagationData.PropagateFractures = true;
            return;
        }
        else if (PropagateTop == true) // Top Propagation
        {
            MyPropagation.TipCoordinates[0] = TopEndX/TopEndDen;
            MyPropagation.TipCoordinates[1] = TopEndY/TopEndDen;
            MyPropagation.TipCoordinates[2] = TopEndZ/TopEndDen;

            // Compute new Tip RotationMatrix
            this->CalculateNewTipRotationMatrix(rAuxPropagationVariables,MyPropagation.TipCoordinates,
                                                MyPropagation.LeftInitCoordinates,PropagationDistance);
            // LeftEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] += 0.5*PropagationHeight;
            noalias(MyPropagation.LeftEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // RightEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] -= 0.5*PropagationHeight;
            noalias(MyPropagation.RightEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] += 0.5*PropagationWidth;
            noalias(MyPropagation.TopEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] -= 0.5*PropagationWidth;
            noalias(MyPropagation.BotEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

            rPropagationData.PropagationVector.push_back(MyPropagation);
            rPropagationData.PropagateFractures = true;
            return;
        }
        else if (PropagateBot == true) // Bot Propagation
        {
            MyPropagation.TipCoordinates[0] = BotEndX/BotEndDen;
            MyPropagation.TipCoordinates[1] = BotEndY/BotEndDen;
            MyPropagation.TipCoordinates[2] = BotEndZ/BotEndDen;

            // Compute new Tip RotationMatrix
            this->CalculateNewTipRotationMatrix(rAuxPropagationVariables,MyPropagation.TipCoordinates,
                                                MyPropagation.LeftInitCoordinates,PropagationDistance);
            // LeftEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] += 0.5*PropagationHeight;
            noalias(MyPropagation.LeftEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // RightEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] -= 0.5*PropagationHeight;
            noalias(MyPropagation.RightEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // TopEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] += 0.5*PropagationWidth;
            noalias(MyPropagation.TopEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            // BotEndCoordinates
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[2] -= 0.5*PropagationWidth;
            noalias(MyPropagation.BotEndCoordinates) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);

            rPropagationData.PropagationVector.push_back(MyPropagation);
            rPropagationData.PropagateFractures = true;
            return;
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateTipRotationMatrix(
        BoundedMatrix<double,3,3>& rRotationMatrix,
        const array_1d<double,3>& TipPointCoordinates,
        const array_1d<double,3>& LeftPointCoordinates,
        const array_1d<double,3>& RightPointCoordinates)
    {
        array_1d<double, 3> CenterPoint;
        noalias(CenterPoint) = 0.5 * (LeftPointCoordinates + RightPointCoordinates);

        // Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = TipPointCoordinates - CenterPoint;
        double inv_norm = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm;
        Vx[1] *= inv_norm;
        Vx[2] *= inv_norm;

        // Vector in local y direction
        array_1d<double, 3> Vy;
        noalias(Vy) = LeftPointCoordinates - CenterPoint;
        // Unitary vector in local z direction
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        inv_norm = 1.0/norm_2(Vz);
        Vz[0] *= inv_norm;
        Vz[1] *= inv_norm;
        Vz[2] *= inv_norm;

        // Unitary vector in local y direction
        MathUtils<double>::CrossProduct( Vy, Vz, Vx);

        // Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];
        rRotationMatrix(0,2) = Vx[2];

        rRotationMatrix(1,0) = Vy[0];
        rRotationMatrix(1,1) = Vy[1];
        rRotationMatrix(1,2) = Vy[2];

        rRotationMatrix(2,0) = Vz[0];
        rRotationMatrix(2,1) = Vz[1];
        rRotationMatrix(2,2) = Vz[2];
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AverageTipCoordinates(
        double& rTipX,
        double& rTipY,
        double& rTipZ,
        double& rTipDenominator,
        const FracturePoint& MyFracturePoint)
    {
        rTipX += MyFracturePoint.Weight * MyFracturePoint.Damage * (MyFracturePoint.Coordinates[0]);
        rTipY += MyFracturePoint.Weight * MyFracturePoint.Damage * (MyFracturePoint.Coordinates[1]);
        rTipZ += MyFracturePoint.Weight * MyFracturePoint.Damage * (MyFracturePoint.Coordinates[2]);
        rTipDenominator += MyFracturePoint.Weight * MyFracturePoint.Damage;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateNewTipRotationMatrix(
        PropagationLocalVariables& rAuxPropagationVariables,
        const array_1d<double,3>& NewTipCoordinates,
        const array_1d<double,3>& LeftInitCoordinates,
        const double& PropagationDistance)
    {
        // Unitary vector in local x direction
        array_1d<double,3> Vx;
        noalias(Vx) = NewTipCoordinates - rAuxPropagationVariables.TipCoordinates;
        double inv_norm = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm;
        Vx[1] *= inv_norm;
        Vx[2] *= inv_norm;

        array_1d<double,3> CenterPoint;
        noalias(CenterPoint) = rAuxPropagationVariables.TipCoordinates + PropagationDistance*(Vx);

        array_1d<double,3> LeftLine;
        noalias(LeftLine) = NewTipCoordinates - LeftInitCoordinates;

        const double Lambda = (Vx[0]*(CenterPoint[0]-LeftInitCoordinates[0])+Vx[1]*(CenterPoint[1]-LeftInitCoordinates[1])+
                            Vx[2]*(CenterPoint[2]-LeftInitCoordinates[2]))/(LeftLine[0]*Vx[0]+LeftLine[1]*Vx[1]+LeftLine[2]*Vx[2]);

        array_1d<double,3> LeftPoint;
        noalias(LeftPoint) = LeftInitCoordinates + Lambda*LeftLine;

        // Vector in local y direction
        array_1d<double,3> Vy;
        noalias(Vy) = LeftPoint - CenterPoint;
        // Unitary vector in local z direction
        array_1d<double,3> Vz;
        MathUtils<double>::CrossProduct(Vz,Vx,Vy);
        inv_norm = 1.0/norm_2(Vz);
        Vz[0] *= inv_norm;
        Vz[1] *= inv_norm;
        Vz[2] *= inv_norm;

        // Unitary vector in local y direction
        MathUtils<double>::CrossProduct(Vy, Vz, Vx);

        // Rotation Matrix
        rAuxPropagationVariables.RotationMatrix(0,0) = Vx[0];
        rAuxPropagationVariables.RotationMatrix(0,1) = Vx[1];
        rAuxPropagationVariables.RotationMatrix(0,2) = Vx[2];

        rAuxPropagationVariables.RotationMatrix(1,0) = Vy[0];
        rAuxPropagationVariables.RotationMatrix(1,1) = Vy[1];
        rAuxPropagationVariables.RotationMatrix(1,2) = Vy[2];

        rAuxPropagationVariables.RotationMatrix(2,0) = Vz[0];
        rAuxPropagationVariables.RotationMatrix(2,1) = Vz[1];
        rAuxPropagationVariables.RotationMatrix(2,2) = Vz[2];

        noalias(rAuxPropagationVariables.TipLocalCoordinates) = prod(rAuxPropagationVariables.RotationMatrix,CenterPoint);
    }

}; // Class FracturePropagation3DUtilities

} // namespace Kratos.

#endif /* KRATOS_PROPAGATE_FRACTURES_3D_UTILITIES defined */
