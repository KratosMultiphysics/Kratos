
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


#if !defined(KRATOS_NONLOCAL_DAMAGE_3D_UTILITIES )
#define  KRATOS_NONLOCAL_DAMAGE_3D_UTILITIES

// Application includes
#include "custom_utilities/nonlocal_damage_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class NonlocalDamage3DUtilities : public NonlocalDamageUtilities
{

public:

    typedef NonlocalDamageUtilities::GaussPoint GaussPoint;
    using NonlocalDamageUtilities::mGaussPointList;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    struct Utility3DVariables
    {
        double X_max, X_min, Y_max, Y_min, Z_max, Z_min;
        int NRows, NColumns, NSections;
        double RowSize, ColumnSize, SectionSize;
        
        std::vector< std::vector< std::vector< std::vector<GaussPoint*> > > > GaussPointCellMatrix;
    };
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:

    KRATOS_CLASS_POINTER_DEFINITION( NonlocalDamage3DUtilities );

    /// Default Constructor
    NonlocalDamage3DUtilities() : NonlocalDamageUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~NonlocalDamage3DUtilities() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void SearchGaussPointsNeighbours (Parameters* pParameters, ModelPart& rModelPart) override
    {
        KRATOS_INFO("Nonlocal Damage 3D utility") << "Starting non-local search of neighbours ..." << std::endl;
        
        // Define necessary variables
        Utility3DVariables AuxVariables;

        //Set GaussPoints inside CellMatrix
        this->InitializeNonlocalSearch(AuxVariables,pParameters,rModelPart);

        this->SearchNeighbours(AuxVariables,pParameters,rModelPart);
        
        KRATOS_INFO("Nonlocal Damage 3D utility") << "... search of neighbours completed." << std::endl;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonlocalSearch(
        Utility3DVariables& rAuxVariables,
        Parameters* pParameters,
        ModelPart& rModelPart)
    {
        // Compute GaussPointsCells dimensions
        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPart);
        
        rAuxVariables.GaussPointCellMatrix.resize(rAuxVariables.NRows);
        for(int i = 0; i < rAuxVariables.NRows; i++)
            rAuxVariables.GaussPointCellMatrix[i].resize(rAuxVariables.NColumns);
        for(int i = 0; i < rAuxVariables.NRows; i++)
        {
            for(int j = 0; j < rAuxVariables.NColumns; j++)
                rAuxVariables.GaussPointCellMatrix[i][j].resize(rAuxVariables.NSections);
        }
        
        // Locate GaussPoints inside CellMatrix
        unsigned int NGPoints = 0;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        array_1d<double,3> AuxGlobalCoordinates;
        array_1d<double,3> AuxLocalCoordinates;
        
        Parameters& rParameters = *pParameters;
        unsigned int NumBodySubModelParts = rParameters["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& BodySubModelPart = rModelPart.GetSubModelPart(rParameters["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(BodySubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = BodySubModelPart.ElementsBegin();
            
            // Loop through all body elements
            #pragma omp parallel for private(MyIntegrationMethod,AuxGlobalCoordinates,AuxLocalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                MyIntegrationMethod = itElem->GetIntegrationMethod();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
                unsigned int NumGPoints = IntegrationPoints.size();
                Vector detJContainer(NumGPoints);
                rGeom.DeterminantOfJacobian(detJContainer,MyIntegrationMethod);
                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(NumGPoints);
                itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,CurrentProcessInfo);
                int Row;
                int Column;
                int Section;

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // MyGaussPoint Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(AuxGlobalCoordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // MyGaussPoint Weight
                    double Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // MyGaussPoint Row, Column and Section
                    Row = int((rAuxVariables.Y_max-AuxGlobalCoordinates[1])/rAuxVariables.RowSize);
                    Column = int((AuxGlobalCoordinates[0]-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    Section = int((AuxGlobalCoordinates[2]-rAuxVariables.Z_min)/rAuxVariables.SectionSize);

                    #pragma omp critical
                    {
                        // Push Back GaussPoint (Heap)
                        mGaussPointList.push_back( new GaussPoint(ConstitutiveLawVector[GPoint],AuxGlobalCoordinates,Weight) );
                        rAuxVariables.GaussPointCellMatrix[Row][Column][Section].push_back(mGaussPointList[NGPoints]);
                        NGPoints++;
                    }
                }
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SearchNeighbours(
        Utility3DVariables& rAuxVariables,
        Parameters* pParameters,
        ModelPart& rModelPart)
    {
        int NGPoints = static_cast<int>(mGaussPointList.size());
        double CharacteristicLength = (*pParameters)["characteristic_length"].GetDouble();
        
        // Loop through all Gauss Points
        #pragma omp parallel for
        for(int i = 0; i < NGPoints; i++)
        {
            GaussPoint& rMyGaussPoint = *(mGaussPointList[i]);

            double X_left = rMyGaussPoint.Coordinates[0] - CharacteristicLength;
            double X_right = rMyGaussPoint.Coordinates[0] + CharacteristicLength;
            double Y_top = rMyGaussPoint.Coordinates[1] + CharacteristicLength;
            double Y_bot = rMyGaussPoint.Coordinates[1] - CharacteristicLength;
            double Z_back = rMyGaussPoint.Coordinates[2] - CharacteristicLength;
            double Z_front = rMyGaussPoint.Coordinates[2] + CharacteristicLength;

            int Column_left = int((X_left-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
            int Column_right = int((X_right-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
            int Row_top = int((rAuxVariables.Y_max-Y_top)/rAuxVariables.RowSize);
            int Row_bot = int((rAuxVariables.Y_max-Y_bot)/rAuxVariables.RowSize);
            int Section_back = int((Z_back-rAuxVariables.Z_min)/rAuxVariables.SectionSize);
            int Section_front = int((Z_front-rAuxVariables.Z_min)/rAuxVariables.SectionSize);

            if(Column_left < 0) Column_left = 0;
            if(Column_right >= rAuxVariables.NColumns) Column_right = rAuxVariables.NColumns-1;
            if(Row_top < 0) Row_top = 0;
            if(Row_bot >= rAuxVariables.NRows) Row_bot = rAuxVariables.NRows-1;
            if(Section_back < 0) Section_back = 0;
            if(Section_front >= rAuxVariables.NSections) Section_front = rAuxVariables.NSections-1;

            // Search GaussPoints neighbours
            for(int s = Section_back; s <= Section_front; s++)
            {
                for(int k = Row_top; k <= Row_bot; k++)
                {
                    for(int l = Column_left; l <= Column_right; l++)
                    {
                        for(unsigned int m = 0; m < rAuxVariables.GaussPointCellMatrix[k][l][s].size(); m++)
                        {
                            if ( (&rMyGaussPoint) != rAuxVariables.GaussPointCellMatrix[k][l][s][m] )
                            {
                                GaussPoint& rMyNeighbourPoint = *(rAuxVariables.GaussPointCellMatrix[k][l][s][m]);

                                double Distance = sqrt((rMyNeighbourPoint.Coordinates[0]-rMyGaussPoint.Coordinates[0])*(rMyNeighbourPoint.Coordinates[0]-rMyGaussPoint.Coordinates[0]) +
                                                (rMyNeighbourPoint.Coordinates[1]-rMyGaussPoint.Coordinates[1])*(rMyNeighbourPoint.Coordinates[1]-rMyGaussPoint.Coordinates[1]) +
                                                (rMyNeighbourPoint.Coordinates[2]-rMyGaussPoint.Coordinates[2])*(rMyNeighbourPoint.Coordinates[2]-rMyGaussPoint.Coordinates[2]));

                                if(Distance <= CharacteristicLength)
                                {
                                    #pragma omp critical
                                    {
                                        rMyGaussPoint.NeighbourPoints.push_back(&rMyNeighbourPoint);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeNeighbourDistance(
        double& rDistance,
        const GaussPoint& ReceiverPoint,
        const GaussPoint& SourcePoint) override
    {
        rDistance = sqrt((SourcePoint.Coordinates[0]-ReceiverPoint.Coordinates[0])*(SourcePoint.Coordinates[0]-ReceiverPoint.Coordinates[0]) +
                         (SourcePoint.Coordinates[1]-ReceiverPoint.Coordinates[1])*(SourcePoint.Coordinates[1]-ReceiverPoint.Coordinates[1]) +
                         (SourcePoint.Coordinates[2]-ReceiverPoint.Coordinates[2])*(SourcePoint.Coordinates[2]-ReceiverPoint.Coordinates[2]));
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void ComputeCellMatrixDimensions(
        Utility3DVariables& rAuxVariables,
        ModelPart& rModelPart)
    {
        // Compute X, Y and Z limits of the current geometry
        unsigned int NumThreads = OpenMPUtils::GetNumThreads();
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

        // Compute CellMatrix dimensions
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
    
}; // Class NonlocalDamage3DUtilities

} // namespace Kratos.

#endif /* KRATOS_NONLOCAL_DAMAGE_3D_UTILITIES defined */