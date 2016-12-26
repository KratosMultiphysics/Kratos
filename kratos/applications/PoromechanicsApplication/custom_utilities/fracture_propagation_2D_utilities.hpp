//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               July 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_PROPAGATE_FRACTURES_2D_UTILITIES )
#define  KRATOS_PROPAGATE_FRACTURES_2D_UTILITIES

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

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class FracturePropagation2DUtilities
{

protected:

/// Basic Structs for the utility ---------------------------------------------------------------------------------------------------------------------------------------------

    struct UtilityVariables
    {
        double X_max, X_min, Y_max, Y_min;
        int NRows, NColumns;
        double RowSize, ColumnSize;
    };

/// Structs for fracture propagation check ------------------------------------------------------------------------------------------------------------------------------------

    struct FracturePoint
    {
        array_1d<double,2> Coordinates;
        double Damage, Weight, TipDistance, CosAngle;
    };
    
    ///------------------------------------------------------------------------------------

    struct Propagation
    {
        int MotherFractureId;
        array_1d<double,3> TopInitCoordinates;
        array_1d<double,3> BotInitCoordinates;
        array_1d<double,3> TipCoordinates;
    };

    ///------------------------------------------------------------------------------------

    struct Bifurcation
    {
        int MotherFractureId;
        array_1d<double,3> TopInitCoordinates;
        array_1d<double,3> BotInitCoordinates;
        array_1d<double,3> TopTipCoordinates;
        array_1d<double,3> BotTipCoordinates;
    };

    ///------------------------------------------------------------------------------------

    struct GlobalPropagationVariables
    {
        std::vector< std::vector< std::vector<FracturePoint> > > FracturePointsCellMatrix;

        bool PropagateFractures;
        std::vector<Propagation> PropagationVector;
        std::vector<Bifurcation> BifurcationVector;
    };

    ///------------------------------------------------------------------------------------

    struct LocalPropagationVariables
    {
        boost::numeric::ublas::bounded_matrix<double,2,2> RotationMatrix;
        array_1d<double,2> TipCoordinates;
        array_1d<double,2> TipLocalCoordinates;
        std::vector<FracturePoint*> TopFrontFracturePoints;
        std::vector<FracturePoint*> BotFrontFracturePoints;
    };
    
/// Structs for mapping model parts -------------------------------------------------------------------------------------------------------------------------------------------

    struct GaussPointOld
    {
        array_1d<double,3> Coordinates;
        double StateVariable, Damage, Weight;
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION( FracturePropagation2DUtilities );

    /// Constructor
    FracturePropagation2DUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~FracturePropagation2DUtilities() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool CheckFracturePropagation (Parameters& rParameters, ModelPart& rModelPart, const bool& move_mesh_flag)
    {
        // Define necessary variables
        UtilityVariables AuxVariables;
        GlobalPropagationVariables PropagationData;
        
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

        // Note: after the CheckFracturePropagation method the rModelPartOld is in the reference state

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
        GlobalPropagationVariables& rPropagationData,
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
        GlobalPropagationVariables& rPropagationData,
        const UtilityVariables& AuxVariables,
        Parameters& rParameters)
    {
        LocalPropagationVariables AuxPropagationVariables;

        // Tip Coordinates
        for(unsigned int i = 0; i < 2; i++)
            AuxPropagationVariables.TipCoordinates[i] = rParameters["fractures_list"][itFracture]["tip_point"]["coordinates"][i].GetDouble();

        // Tip Rotation Matrix
        this->CalculateTipRotationMatrix(itFracture,AuxPropagationVariables.RotationMatrix,rParameters);

        // Tip search area
        const double PropagationLength = rParameters["fracture_data"]["propagation_length"].GetDouble();

        double X_left = AuxPropagationVariables.TipCoordinates[0] - PropagationLength;
        double X_right = AuxPropagationVariables.TipCoordinates[0] + PropagationLength;
        double Y_top = AuxPropagationVariables.TipCoordinates[1] + PropagationLength;
        double Y_bot = AuxPropagationVariables.TipCoordinates[1] - PropagationLength;

        int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
        int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
        int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
        int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

        if(Column_left < 0) Column_left = 0;
        if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
        if(Row_top < 0) Row_top = 0;
        if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;

        // Search FracturePoints neighbours around the tip
        std::vector<FracturePoint*> TipNeighbours;
        noalias(AuxPropagationVariables.TipLocalCoordinates) = prod(AuxPropagationVariables.RotationMatrix,AuxPropagationVariables.TipCoordinates);
        array_1d<double,2> OtherLocalCoordinates;
        for(int i = Row_top; i <= Row_bot; i++)
        {
            for(int j = Column_left; j<= Column_right; j++)
            {
                for(unsigned int k = 0; k < rPropagationData.FracturePointsCellMatrix[i][j].size(); k++)
                {
                    FracturePoint& rOtherPoint = rPropagationData.FracturePointsCellMatrix[i][j][k];

                    rOtherPoint.TipDistance = sqrt((rOtherPoint.Coordinates[0]-AuxPropagationVariables.TipCoordinates[0])*(rOtherPoint.Coordinates[0]-AuxPropagationVariables.TipCoordinates[0]) +
                                    (rOtherPoint.Coordinates[1]-AuxPropagationVariables.TipCoordinates[1])*(rOtherPoint.Coordinates[1]-AuxPropagationVariables.TipCoordinates[1]));

                    if(rOtherPoint.TipDistance <= PropagationLength)
                    {
                        TipNeighbours.push_back(&rOtherPoint);

                        noalias(OtherLocalCoordinates) = prod(AuxPropagationVariables.RotationMatrix,rOtherPoint.Coordinates);
                        
                        // FrontFracturePoints
                        if(OtherLocalCoordinates[0] >= AuxPropagationVariables.TipLocalCoordinates[0])
                        {
                            // TopFrontFracturePoints
                            if(OtherLocalCoordinates[1] >= AuxPropagationVariables.TipLocalCoordinates[1])
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
        {
            rPropagationData.PropagateFractures = true;
            this->PropagateFracture(itFracture,rPropagationData,AuxPropagationVariables,rParameters);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeCheckFracture(
        const GlobalPropagationVariables& PropagationData,
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
        const GlobalPropagationVariables& PropagationData,
        Parameters& rParameters)
    {

        std::fstream PropDataFile;
        PropDataFile.open ("PropagationData.tcl", std::fstream::out | std::fstream::trunc);
        PropDataFile.precision(12);

        PropDataFile << "set PropagationData [list]" << std::endl;

        // GiDPath
        PropDataFile << "lappend PropagationData \"" << rParameters["fracture_data"]["gid_path"].GetString() << "\"" << std::endl;

        int Id;

        // BodySurfacesDict
        for(unsigned int i = 0; i < rParameters["body_surfaces_list"].size(); i++)
        {
            Id = rParameters["body_surfaces_list"][i]["id"].GetInt();

            PropDataFile << "dict set BodySurfacesDict " << Id << " Layer \""
                         << rParameters["body_surfaces_list"][i]["layer"].GetString() << "\"" << std::endl;

            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["body_surfaces_list"][i]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["body_surfaces_list"][i]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set BodySurfacesDict " << Id << " Groups $Groups" << std::endl;

            PropDataFile << "set Lines [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["body_surfaces_list"][i]["lines"].size(); j++)
            {
                PropDataFile << "lappend Lines " << rParameters["body_surfaces_list"][i]["lines"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set BodySurfacesDict " << Id << " Lines $Lines" << std::endl;
        }
        PropDataFile << "lappend PropagationData $BodySurfacesDict" << std::endl;

        // FracturesDict
        for(unsigned int i = 0; i < rParameters["fractures_list"].size(); i++)
        {
            Id = rParameters["fractures_list"][i]["id"].GetInt();

            // TipPoint
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Id "
                         << rParameters["fractures_list"][i]["tip_point"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Layer \""
                         << rParameters["fractures_list"][i]["tip_point"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["tip_point"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["tip_point"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Groups $Groups" << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["tip_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["tip_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["tip_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TipPoint Coordinates $Coordinates" << std::endl;
            // TopPoint
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Id "
                         << rParameters["fractures_list"][i]["top_point"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Layer \""
                         << rParameters["fractures_list"][i]["top_point"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["top_point"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["top_point"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Groups $Groups" << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["top_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["top_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["top_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TopPoint Coordinates $Coordinates" << std::endl;
            // BotPoint
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Id "
                         << rParameters["fractures_list"][i]["bot_point"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Layer \""
                         << rParameters["fractures_list"][i]["bot_point"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["bot_point"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["bot_point"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Groups $Groups" << std::endl;
            PropDataFile << "set Coordinates \"" << rParameters["fractures_list"][i]["bot_point"]["coordinates"][0].GetDouble()
                         << " " << rParameters["fractures_list"][i]["bot_point"]["coordinates"][1].GetDouble() << " "
                         << rParameters["fractures_list"][i]["bot_point"]["coordinates"][2].GetDouble() << "\"" << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " BotPoint Coordinates $Coordinates" << std::endl;
            // TopLine
            PropDataFile << "dict set FracturesDict " << Id << " TopLine Id "
                         << rParameters["fractures_list"][i]["top_line"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " TopLine Layer \""
                         << rParameters["fractures_list"][i]["top_line"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["top_line"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["top_line"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " TopLine Groups $Groups" << std::endl;
            // BotLine
            PropDataFile << "dict set FracturesDict " << Id << " BotLine Id "
                         << rParameters["fractures_list"][i]["bot_line"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " BotLine Layer \""
                         << rParameters["fractures_list"][i]["bot_line"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["bot_line"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["bot_line"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BotLine Groups $Groups" << std::endl;
            // InterfaceSurface
            PropDataFile << "dict set FracturesDict " << Id << " InterfaceSurface Id "
                         << rParameters["fractures_list"][i]["interface_surface"]["id"].GetInt() << std::endl;
            PropDataFile << "dict set FracturesDict " << Id << " InterfaceSurface Layer \""
                         << rParameters["fractures_list"][i]["interface_surface"]["layer"].GetString() << "\"" << std::endl;
            PropDataFile << "set Groups [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["interface_surface"]["groups"].size(); j++)
            {
                PropDataFile << "lappend Groups \"" << rParameters["fractures_list"][i]["interface_surface"]["groups"][j].GetString() << "\"" << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " InterfaceSurface Groups $Groups" << std::endl;
            // BodySurfaces
            PropDataFile << "set BodySurfaces [list]" << std::endl;
            for(unsigned int j = 0; j < rParameters["fractures_list"][i]["body_surfaces"].size(); j++)
            {
                PropDataFile << "lappend BodySurfaces " << rParameters["fractures_list"][i]["body_surfaces"][j].GetInt() << std::endl;
            }
            PropDataFile << "dict set FracturesDict " << Id << " BodySurfaces $BodySurfaces" << std::endl;
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

            PropDataFile << "set Coordinates \"" << PropagationData.BifurcationVector[i].TopTipCoordinates[0]
                         << " " << PropagationData.BifurcationVector[i].TopTipCoordinates[1] << " "
                         << PropagationData.BifurcationVector[i].TopTipCoordinates[2] << "\"" << std::endl;
            PropDataFile << "dict set BifurcationDict " << i << " TopTipCoordinates $Coordinates" << std::endl;

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
        std::vector< std::vector< std::vector<Element::Pointer> > > ElementOldCellMatrix;
        ElementOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++) ElementOldCellMatrix[i].resize(AuxVariables.NColumns);

        // Locate Old Elments in cells
        double X_me;
        double Y_me;
        int PointsNumber;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();
                
        // Loop through all BodySubModelParts
        for(unsigned int m = 0; m < NumBodySubModelParts; m++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][m].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(X_me,Y_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;

                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;

                for(int k = Row_top; k <= Row_bot; k++)
                {
                    for(int l = Column_left; l <= Column_right; l++)
                    {
                        #pragma omp critical
                        {
                            ElementOldCellMatrix[k][l].push_back((*(itElemOld.base())));
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

            #pragma omp parallel for private(X_me,Y_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;

                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;

                for(int k = Row_top; k <= Row_bot; k++)
                {
                    for(int l = Column_left; l <= Column_right; l++)
                    {
                        #pragma omp critical
                        {
                            ElementOldCellMatrix[k][l].push_back((*(itElemOld.base())));
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
        double Tolerance = 1.0e-5;

        #pragma omp parallel for private(X_me,Y_me,PointsNumber,GlobalCoordinates,LocalCoordinates)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNodeNew = node_begin + i;

            X_me = itNodeNew->X0();
            Y_me = itNodeNew->Y0();

            int Column = int((X_me-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Row = int((AuxVariables.Y_max-Y_me)/AuxVariables.RowSize);

            if(Column >= AuxVariables.NColumns) Column = AuxVariables.NColumns-1;
            else if(Column < 0) Column = 0;
            if(Row >= AuxVariables.NRows) Row = AuxVariables.NRows-1;
            else if(Row < 0) Row = 0;

            noalias(GlobalCoordinates) = itNodeNew->Coordinates(); //Coordinates of new nodes are still in the original position
            bool IsInside = false;
            Element::Pointer pElementOld;

            for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
            {
                pElementOld = ElementOldCellMatrix[Row][Column][m];
                IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
                {
                    pElementOld = ElementOldCellMatrix[Row][Column][m];
                    IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates);
                    PointsNumber = pElementOld->GetGeometry().PointsNumber();
                    if(PointsNumber == 3)
                    {
                        if(((LocalCoordinates[0]+Tolerance)>=0)&&((LocalCoordinates[1]+Tolerance)>=0)&&((LocalCoordinates[1]-Tolerance)<=(1-LocalCoordinates[0])))
                        {
                            IsInside = true;
                            break;
                        }
                    }
                    else if(PointsNumber == 4)
                    {
                        if(((LocalCoordinates[0]+Tolerance)>=-1)&&((LocalCoordinates[0]-Tolerance)<=1)&&((LocalCoordinates[1]+Tolerance)>=-1)&&((LocalCoordinates[1]-Tolerance)<=1))
                        {
                            IsInside = true;
                            break;
                        }
                    }
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
            if( itNodeNew->IsFixed(WATER_PRESSURE)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(WATER_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(WATER_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DT_WATER_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(DT_WATER_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
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
        std::vector< std::vector< std::vector<GaussPointOld> > > BodyGaussPointOldCellMatrix;
        BodyGaussPointOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++) BodyGaussPointOldCellMatrix[i].resize(AuxVariables.NColumns);
        
        std::vector< std::vector< std::vector<GaussPointOld> > > GaussPointOldCellMatrix;
        GaussPointOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++) GaussPointOldCellMatrix[i].resize(AuxVariables.NColumns);

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
                itElem->GetValueOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoOld);
                std::vector<double> DamageVector(NumGPoints);
                itElem->GetValueOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,CurrentProcessInfoOld);
                int Row;
                int Column;

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
                    MyGaussPointOld.Damage = DamageVector[GPoint];

                    // GaussPointOld Row and Column
                    Row = int((AuxVariables.Y_max-MyGaussPointOld.Coordinates[1])/AuxVariables.RowSize);
                    Column = int((MyGaussPointOld.Coordinates[0]-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    #pragma omp critical
                    {
                        BodyGaussPointOldCellMatrix[Row][Column].push_back(MyGaussPointOld);
                        GaussPointOldCellMatrix[Row][Column].push_back(MyGaussPointOld);
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
                MyIntegrationMethod = GeometryData::GI_GAUSS_1;
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<double> StateVariableVector(NumGPoints);
                itElem->GetValueOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoOld);
                int Row;
                int Column;

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
                    MyGaussPointOld.Damage = MyGaussPointOld.StateVariable;

                    // GaussPointOld Row and Column
                    Row = int((AuxVariables.Y_max-MyGaussPointOld.Coordinates[1])/AuxVariables.RowSize);
                    Column = int((MyGaussPointOld.Coordinates[0]-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    #pragma omp critical
                    {
                        GaussPointOldCellMatrix[Row][Column].push_back(MyGaussPointOld);
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

                    // GaussPointNew search area
                    double X_left = X_me - PropagationLength;
                    double X_right = X_me + PropagationLength;
                    double Y_top = Y_me + PropagationLength;
                    double Y_bot = Y_me - PropagationLength;

                    int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                    int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                    if(Column_left < 0) Column_left = 0;
                    if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                    if(Row_top < 0) Row_top = 0;
                    if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;

                    // Search GaussPointOld neighbours around the GaussPointNew and compute nonlocal state variable
                    double Numerator = 0.0;
                    double WeightingFunctionDenominator = 0.0;
                    double Distance;
                    for(int k = Row_top; k <= Row_bot; k++)
                    {
                        for(int l = Column_left; l<= Column_right; l++)
                        {
                            for(unsigned int m = 0; m < BodyGaussPointOldCellMatrix[k][l].size(); m++)
                            {
                                GaussPointOld& rOtherGaussPointOld = BodyGaussPointOldCellMatrix[k][l][m];

                                Distance = sqrt((rOtherGaussPointOld.Coordinates[0]-X_me)*(rOtherGaussPointOld.Coordinates[0]-X_me) +
                                                (rOtherGaussPointOld.Coordinates[1]-Y_me)*(rOtherGaussPointOld.Coordinates[1]-Y_me));

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

                    // Save computed stateVariable
                    if(WeightingFunctionDenominator > 0.0)
                        StateVariableVector[GPoint] = Numerator/WeightingFunctionDenominator;
                    else
                        StateVariableVector[GPoint] = DamageThreshold;
                }
                // Set stateVariable of new GaussPoints
                itElem->SetValueOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoNew);
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
                MyIntegrationMethod = GeometryData::GI_GAUSS_1;
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

                    // GaussPointNew search area
                    double X_left = X_me - PropagationLength;
                    double X_right = X_me + PropagationLength;
                    double Y_top = Y_me + PropagationLength;
                    double Y_bot = Y_me - PropagationLength;

                    int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                    int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                    int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                    if(Column_left < 0) Column_left = 0;
                    if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                    if(Row_top < 0) Row_top = 0;
                    if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;

                    // Search GaussPointOld neighbours around the GaussPointNew and compute nonlocal state variable
                    double Numerator = 0.0;
                    double WeightingFunctionDenominator = 0.0;
                    double Distance;
                    for(int k = Row_top; k <= Row_bot; k++)
                    {
                        for(int l = Column_left; l<= Column_right; l++)
                        {
                            for(unsigned int m = 0; m < GaussPointOldCellMatrix[k][l].size(); m++)
                            {
                                GaussPointOld& rOtherGaussPointOld = GaussPointOldCellMatrix[k][l][m];

                                Distance = sqrt((rOtherGaussPointOld.Coordinates[0]-X_me)*(rOtherGaussPointOld.Coordinates[0]-X_me) +
                                                (rOtherGaussPointOld.Coordinates[1]-Y_me)*(rOtherGaussPointOld.Coordinates[1]-Y_me));

                                if(Distance <= PropagationLength)
                                {
                                    Numerator += rOtherGaussPointOld.Weight
                                                *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength))
                                                *rOtherGaussPointOld.Damage;
                                    WeightingFunctionDenominator += rOtherGaussPointOld.Weight
                                                                    *exp(-4.0*Distance*Distance/(PropagationLength*PropagationLength));
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
                itElem->SetValueOnIntegrationPoints(STATE_VARIABLE,StateVariableVector,CurrentProcessInfoNew);
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
        // Compute X and Y limits of the current geometry
        unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector<double> X_max_partition(NumThreads);
        std::vector<double> X_min_partition(NumThreads);
        std::vector<double> Y_max_partition(NumThreads);
        std::vector<double> Y_min_partition(NumThreads);
        
        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();
        
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            X_max_partition[k] = node_begin->X();
            X_min_partition[k] = X_max_partition[k];
            Y_max_partition[k] = node_begin->Y();
            Y_min_partition[k] = Y_max_partition[k];

            double X_me, Y_me;

            #pragma omp for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                X_me = itNode->X();
                Y_me = itNode->Y();

                if( X_me > X_max_partition[k] ) X_max_partition[k] = X_me;
                else if( X_me < X_min_partition[k] ) X_min_partition[k] = X_me;

                if( Y_me > Y_max_partition[k] ) Y_max_partition[k] = Y_me;
                else if( Y_me < Y_min_partition[k] ) Y_min_partition[k] = Y_me;
            }
        }

        rAuxVariables.X_max = X_max_partition[0];
        rAuxVariables.X_min = X_min_partition[0];
        rAuxVariables.Y_max = Y_max_partition[0];
        rAuxVariables.Y_min = Y_min_partition[0];

        for(unsigned int i=1; i < NumThreads; i++)
        {
            if(X_max_partition[i] > rAuxVariables.X_max) rAuxVariables.X_max = X_max_partition[i];
            if(X_min_partition[i] < rAuxVariables.X_min) rAuxVariables.X_min = X_min_partition[i];
            if(Y_max_partition[i] > rAuxVariables.Y_max) rAuxVariables.Y_max = Y_max_partition[i];
            if(Y_min_partition[i] < rAuxVariables.Y_min) rAuxVariables.Y_min = Y_min_partition[i];
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
        rAuxVariables.NRows = int((rAuxVariables.Y_max-rAuxVariables.Y_min)/(AverageElementLength));
        rAuxVariables.NColumns = int((rAuxVariables.X_max-rAuxVariables.X_min)/(AverageElementLength));
        if(rAuxVariables.NRows <= 0) rAuxVariables.NRows = 1;
        if(rAuxVariables.NColumns <= 0) rAuxVariables.NColumns = 1;
        rAuxVariables.RowSize = (rAuxVariables.Y_max-rAuxVariables.Y_min)/rAuxVariables.NRows;
        rAuxVariables.ColumnSize = (rAuxVariables.X_max-rAuxVariables.X_min)/rAuxVariables.NColumns;
    }

/// Fracture Propagation Check ------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetFracturePoints(
        GlobalPropagationVariables& rPropagationData,
        UtilityVariables& rAuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPart)
    {
        // Compute FracturePointsCellMatrix dimensions
        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPart);

        rPropagationData.FracturePointsCellMatrix.resize(rAuxVariables.NRows);
        for(int i = 0; i < rAuxVariables.NRows; i++) rPropagationData.FracturePointsCellMatrix[i].resize(rAuxVariables.NColumns);

        // Locate FracturePoints inside CellMatrix
        FracturePoint MyFracturePoint;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        array_1d<double,3> AuxGlobalCoordinates;
        array_1d<double,3> AuxLocalCoordinates;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& BodySubModelPart = rModelPart.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(BodySubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = BodySubModelPart.ElementsBegin();
            
            // Loop through all body elements
            #pragma omp parallel for private(MyFracturePoint,MyIntegrationMethod,AuxLocalCoordinates,AuxGlobalCoordinates)
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
                itElem->GetValueOnIntegrationPoints(DAMAGE_VARIABLE,DamageVector,CurrentProcessInfo);
                int Row;
                int Column;

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // FracturePoint Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(AuxGlobalCoordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates
                    MyFracturePoint.Coordinates[0] = AuxGlobalCoordinates[0];
                    MyFracturePoint.Coordinates[1] = AuxGlobalCoordinates[1];

                    // FracturePoint Weight
                    MyFracturePoint.Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // FracturePoint Damage
                    MyFracturePoint.Damage = DamageVector[GPoint];

                    // FracturePoint Row and Column
                    Row = int((rAuxVariables.Y_max-MyFracturePoint.Coordinates[1])/rAuxVariables.RowSize);
                    Column = int((MyFracturePoint.Coordinates[0]-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    #pragma omp critical
                    {
                        rPropagationData.FracturePointsCellMatrix[Row][Column].push_back(MyFracturePoint);
                    }
                }
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void PropagateFracture(
        const unsigned int& itFracture,
        GlobalPropagationVariables& rPropagationData,
        LocalPropagationVariables& rAuxPropagationVariables,
        Parameters& rParameters)
    {
        // Check whether we have bifurcation or simple propagation
        
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                {
                    this->ComputeCosAngle(*(rAuxPropagationVariables.TopFrontFracturePoints[i]),rAuxPropagationVariables);
                }
            }
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                {
                    this->ComputeCosAngle(*(rAuxPropagationVariables.BotFrontFracturePoints[i]),rAuxPropagationVariables);
                }
            }
        }
        
        // Compute bifurcation and propagation factors
        const double PropagationLength = rParameters["fracture_data"]["propagation_length"].GetDouble();
        const double PropagationCosAngle = rParameters["fracture_data"]["propagation_cosangle"].GetDouble();
        double PropagationFactor = 0.0;
        double PropagationFactorDenominator = 0.0;
        double TopBifurcationFactor = 0.0;
        double TopBifurcationFactorDenominator = 0.0;
        double BotBifurcationFactor = 0.0;
        double BotBifurcationFactorDenominator = 0.0;

        #pragma omp parallel sections reduction(+:PropagationFactor,PropagationFactorDenominator)
        {
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                {
                    if(rAuxPropagationVariables.TopFrontFracturePoints[i]->CosAngle > PropagationCosAngle)
                    {
                        this->ComputeBifurcationFactor(PropagationFactor,PropagationFactorDenominator,*(rAuxPropagationVariables.TopFrontFracturePoints[i]),PropagationLength);
                    }
                    else
                    {
                        this->ComputeBifurcationFactor(TopBifurcationFactor,TopBifurcationFactorDenominator,*(rAuxPropagationVariables.TopFrontFracturePoints[i]),PropagationLength);
                    }
                }
            }
            #pragma omp section
            {
                for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                {
                    if(rAuxPropagationVariables.BotFrontFracturePoints[i]->CosAngle > PropagationCosAngle)
                    {
                        this->ComputeBifurcationFactor(PropagationFactor,PropagationFactorDenominator,*(rAuxPropagationVariables.BotFrontFracturePoints[i]),PropagationLength);
                    }
                    else
                    {
                        this->ComputeBifurcationFactor(BotBifurcationFactor,BotBifurcationFactorDenominator,*(rAuxPropagationVariables.BotFrontFracturePoints[i]),PropagationLength);
                    }
                }
            }
        }
        if(PropagationFactorDenominator > 1.0e-20)
            PropagationFactor = PropagationFactor/PropagationFactorDenominator;
        else
            PropagationFactor = 0.0;
        if(TopBifurcationFactorDenominator > 1.0e-20)
            TopBifurcationFactor = TopBifurcationFactor/TopBifurcationFactorDenominator;
        else
            TopBifurcationFactor = 0.0;
        if(BotBifurcationFactorDenominator > 1.0e-20)
            BotBifurcationFactor = BotBifurcationFactor/BotBifurcationFactorDenominator;
        else
            BotBifurcationFactor = 0.0;
        
        // Compute Propagation Coordinates
        const double PropagationDamage = rParameters["fracture_data"]["propagation_damage"].GetDouble();
        const double PropagationWidth = rParameters["fracture_data"]["propagation_width"].GetDouble();
        int MotherFractureId = rParameters["fractures_list"][itFracture]["id"].GetInt();
        array_1d<double,2> AuxArray1;
        array_1d<double,2> AuxArray2;

        // Bifurcation
        if(PropagationFactor < PropagationDamage
            && TopBifurcationFactor > PropagationDamage && BotBifurcationFactor > PropagationDamage)
        {
            Bifurcation MyBifurcation;
            MyBifurcation.MotherFractureId = MotherFractureId;

            double TopEndX = 0.0;
            double TopEndY = 0.0;
            double TopEndDen = 0.0;
            double BotEndX = 0.0;
            double BotEndY = 0.0;
            double BotEndDen = 0.0;

            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                    {
                        this->AverageTipCoordinates(TopEndX,TopEndY,TopEndDen,*(rAuxPropagationVariables.TopFrontFracturePoints[i]));
                    }
                }
                #pragma omp section
                {
                    for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                    {
                        this->AverageTipCoordinates(BotEndX,BotEndY,BotEndDen,*(rAuxPropagationVariables.BotFrontFracturePoints[i]));
                    }
                }
            }
            
            MyBifurcation.TopTipCoordinates[0] = TopEndX/TopEndDen;
            MyBifurcation.TopTipCoordinates[1] = TopEndY/TopEndDen;
            MyBifurcation.TopTipCoordinates[2] = 0.0;
            
            MyBifurcation.BotTipCoordinates[0] = BotEndX/BotEndDen;
            MyBifurcation.BotTipCoordinates[1] = BotEndY/BotEndDen;
            MyBifurcation.BotTipCoordinates[2] = 0.0;
            
            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationWidth;
            AuxArray1[1] += 0.5*PropagationWidth;
            noalias(AuxArray2) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            MyBifurcation.TopInitCoordinates[0] = AuxArray2[0];
            MyBifurcation.TopInitCoordinates[1] = AuxArray2[1];
            MyBifurcation.TopInitCoordinates[2] = 0.0;

            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[0] -= PropagationWidth;
            AuxArray1[1] -= 0.5*PropagationWidth;
            noalias(AuxArray2) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            MyBifurcation.BotInitCoordinates[0] = AuxArray2[0];
            MyBifurcation.BotInitCoordinates[1] = AuxArray2[1];
            MyBifurcation.BotInitCoordinates[2] = 0.0;

            rPropagationData.BifurcationVector.push_back(MyBifurcation);
        }
        else
        // Propagation
        {
            Propagation MyPropagation;
            MyPropagation.MotherFractureId = MotherFractureId;

            double TipX = 0.0;
            double TipY = 0.0;
            double TipDen = 0.0;

            #pragma omp parallel sections reduction(+:TipX,TipY,TipDen)
            {
                #pragma omp section
                {
                    for(unsigned int i = 0; i < rAuxPropagationVariables.TopFrontFracturePoints.size(); i++)
                    {
                        this->AverageTipCoordinates(TipX,TipY,TipDen,*(rAuxPropagationVariables.TopFrontFracturePoints[i]));
                    }
                }
                #pragma omp section
                {
                    for(unsigned int i = 0; i < rAuxPropagationVariables.BotFrontFracturePoints.size(); i++)
                    {
                        this->AverageTipCoordinates(TipX,TipY,TipDen,*(rAuxPropagationVariables.BotFrontFracturePoints[i]));
                    }
                }
            }
            
            MyPropagation.TipCoordinates[0] = TipX/TipDen;
            MyPropagation.TipCoordinates[1] = TipY/TipDen;
            MyPropagation.TipCoordinates[2] = 0.0;

            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] += 0.5*PropagationWidth;
            noalias(AuxArray2) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            MyPropagation.TopInitCoordinates[0] = AuxArray2[0];
            MyPropagation.TopInitCoordinates[1] = AuxArray2[1];
            MyPropagation.TopInitCoordinates[2] = 0.0;

            noalias(AuxArray1) = rAuxPropagationVariables.TipLocalCoordinates;
            AuxArray1[1] -= 0.5*PropagationWidth;
            noalias(AuxArray2) = prod(trans(rAuxPropagationVariables.RotationMatrix),AuxArray1);
            MyPropagation.BotInitCoordinates[0] = AuxArray2[0];
            MyPropagation.BotInitCoordinates[1] = AuxArray2[1];
            MyPropagation.BotInitCoordinates[2] = 0.0;

            rPropagationData.PropagationVector.push_back(MyPropagation);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateTipRotationMatrix(
        const unsigned int& itFracture,
        boost::numeric::ublas::bounded_matrix<double,2,2>& rRotationMatrix,
        Parameters& rParameters)
    {
        array_1d<double,3> TipPointCoordinates;
        array_1d<double,3> TopPointCoordinates;
        array_1d<double,3> BotPointCoordinates;
        for(unsigned int i = 0; i < 3; i++)
        {
            TipPointCoordinates[i] = rParameters["fractures_list"][itFracture]["tip_point"]["coordinates"][i].GetDouble();
            TopPointCoordinates[i] = rParameters["fractures_list"][itFracture]["top_point"]["coordinates"][i].GetDouble();
            BotPointCoordinates[i] = rParameters["fractures_list"][itFracture]["bot_point"]["coordinates"][i].GetDouble();
        }

        array_1d<double, 3> pmid0;
        noalias(pmid0) = 0.5 * (TopPointCoordinates + BotPointCoordinates);

        array_1d<double, 3> Vx;
        noalias(Vx) = TipPointCoordinates - pmid0;
        double inv_norm_x = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm_x;
        Vx[1] *= inv_norm_x;

        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];

        rRotationMatrix(1,0) = -Vx[1];
        rRotationMatrix(1,1) = Vx[0];
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ComputeCosAngle(
    FracturePoint& rMyFracturePoint,
    LocalPropagationVariables& rAuxPropagationVariables)
{
    array_1d<double,2> FractureLineVector;
    array_1d<double,2> LocalFractureLineVector;

    // Compute vector between FracturePoint and Crack Tip
    FractureLineVector[0] = rMyFracturePoint.Coordinates[0] - rAuxPropagationVariables.TipCoordinates[0];
    FractureLineVector[1] = rMyFracturePoint.Coordinates[1] - rAuxPropagationVariables.TipCoordinates[1];

    // Compute angle between FracturePoint and Crack tip
    FractureLineVector[0] = FractureLineVector[0]/rMyFracturePoint.TipDistance;
    FractureLineVector[1] = FractureLineVector[1]/rMyFracturePoint.TipDistance;
    noalias(LocalFractureLineVector) = prod(rAuxPropagationVariables.RotationMatrix,FractureLineVector);

    rMyFracturePoint.CosAngle = LocalFractureLineVector[0];
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeBifurcationFactor(
        double& rBifurcationFactor,
        double& rBifurcationFactorDenominator,
        const FracturePoint& MyFracturePoint,
        const double& PropagationLength)
    {
        rBifurcationFactor += MyFracturePoint.Weight*
                            exp(-4.0*(MyFracturePoint.TipDistance)*(MyFracturePoint.TipDistance)/(PropagationLength*PropagationLength))*MyFracturePoint.Damage;
        rBifurcationFactorDenominator += MyFracturePoint.Weight*
                                        exp(-4.0*(MyFracturePoint.TipDistance)*(MyFracturePoint.TipDistance)/(PropagationLength*PropagationLength));
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AverageTipCoordinates(
        double& rTipX,
        double& rTipY,
        double& rTipDenominator,
        const FracturePoint& MyFracturePoint)
    {
        rTipX += MyFracturePoint.Weight * MyFracturePoint.Damage * (MyFracturePoint.Coordinates[0]);
        rTipY += MyFracturePoint.Weight * MyFracturePoint.Damage * (MyFracturePoint.Coordinates[1]);
        rTipDenominator += MyFracturePoint.Weight * MyFracturePoint.Damage;
    }

}; // Class FracturePropagation2DUtilities

} // namespace Kratos.

#endif /* KRATOS_PROPAGATE_FRACTURES_2D_UTILITIES defined */
