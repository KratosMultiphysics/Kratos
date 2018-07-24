
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


#if !defined(KRATOS_NONLOCAL_DAMAGE_2D_UTILITIES )
#define  KRATOS_NONLOCAL_DAMAGE_2D_UTILITIES

// Application includes
#include "custom_utilities/nonlocal_damage_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class NonlocalDamage2DUtilities : public NonlocalDamageUtilities
{

public:

    typedef NonlocalDamageUtilities::GaussPoint GaussPoint;
    using NonlocalDamageUtilities::mGaussPointList;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    struct Utility2DVariables
    {
        double X_max, X_min, Y_max, Y_min;
        int NRows, NColumns;
        double RowSize, ColumnSize;
        
        std::vector< std::vector< std::vector<GaussPoint*> > > GaussPointCellMatrix;
    };
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:

    KRATOS_CLASS_POINTER_DEFINITION( NonlocalDamage2DUtilities );

    /// Default Constructor
    NonlocalDamage2DUtilities() : NonlocalDamageUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~NonlocalDamage2DUtilities() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void SearchGaussPointsNeighbours (Parameters* pParameters, ModelPart& rModelPart) override
    {
        KRATOS_INFO("Nonlocal Damage 2D utility") << "Starting non-local search of neighbours ..." << std::endl;
        
        // Define necessary variables
        Utility2DVariables AuxVariables;

        //Set GaussPoints inside CellMatrix
        this->InitializeNonlocalSearch(AuxVariables,pParameters,rModelPart);

        this->SearchNeighbours(AuxVariables,pParameters,rModelPart);
        
        KRATOS_INFO("Nonlocal Damage 2D utility") << "... search of neighbours completed." << std::endl;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonlocalSearch(
        Utility2DVariables& rAuxVariables,
        Parameters* pParameters,
        ModelPart& rModelPart)
    {
        // Compute GaussPointsCells dimensions
        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPart);
        
        rAuxVariables.GaussPointCellMatrix.resize(rAuxVariables.NRows);
        for(int i = 0; i < rAuxVariables.NRows; i++) rAuxVariables.GaussPointCellMatrix[i].resize(rAuxVariables.NColumns);
        
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

                    // MyGaussPoint Row and Column
                    Row = int((rAuxVariables.Y_max-AuxGlobalCoordinates[1])/rAuxVariables.RowSize);
                    Column = int((AuxGlobalCoordinates[0]-rAuxVariables.X_min)/rAuxVariables.ColumnSize);

                    #pragma omp critical
                    {
                        // Push Back GaussPoint (Heap)
                        mGaussPointList.push_back( new GaussPoint(ConstitutiveLawVector[GPoint],AuxGlobalCoordinates,Weight) );
                        rAuxVariables.GaussPointCellMatrix[Row][Column].push_back(mGaussPointList[NGPoints]);
                        NGPoints++;
                    }
                }
            }
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SearchNeighbours(
        Utility2DVariables& rAuxVariables,
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

            int Column_left = int((X_left-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
            int Column_right = int((X_right-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
            int Row_top = int((rAuxVariables.Y_max-Y_top)/rAuxVariables.RowSize);
            int Row_bot = int((rAuxVariables.Y_max-Y_bot)/rAuxVariables.RowSize);

            if(Column_left < 0) Column_left = 0;
            if(Column_right >= rAuxVariables.NColumns) Column_right = rAuxVariables.NColumns-1;
            if(Row_top < 0) Row_top = 0;
            if(Row_bot >= rAuxVariables.NRows) Row_bot = rAuxVariables.NRows-1;

            // Search GaussPoints neighbours
            for(int j = Row_top; j <= Row_bot; j++)
            {
                for(int k = Column_left; k <= Column_right; k++)
                {
                    for(unsigned int l = 0; l < rAuxVariables.GaussPointCellMatrix[j][k].size(); l++)
                    {
                        if ( (&rMyGaussPoint) != rAuxVariables.GaussPointCellMatrix[j][k][l] )
                        {
                            GaussPoint& rMyNeighbourPoint = *(rAuxVariables.GaussPointCellMatrix[j][k][l]);

                            double Distance = sqrt((rMyNeighbourPoint.Coordinates[0]-rMyGaussPoint.Coordinates[0])*(rMyNeighbourPoint.Coordinates[0]-rMyGaussPoint.Coordinates[0]) +
                                                   (rMyNeighbourPoint.Coordinates[1]-rMyGaussPoint.Coordinates[1])*(rMyNeighbourPoint.Coordinates[1]-rMyGaussPoint.Coordinates[1]));

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

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeNeighbourDistance(
        double& rDistance,
        const GaussPoint& ReceiverPoint,
        const GaussPoint& SourcePoint) override
    {
        rDistance = sqrt((SourcePoint.Coordinates[0]-ReceiverPoint.Coordinates[0])*(SourcePoint.Coordinates[0]-ReceiverPoint.Coordinates[0]) +
                         (SourcePoint.Coordinates[1]-ReceiverPoint.Coordinates[1])*(SourcePoint.Coordinates[1]-ReceiverPoint.Coordinates[1]));
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void ComputeCellMatrixDimensions(
        Utility2DVariables& rAuxVariables,
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

        // Compute CellMatrix dimensions
        rAuxVariables.NRows = int((rAuxVariables.Y_max-rAuxVariables.Y_min)/(AverageElementLength));
        rAuxVariables.NColumns = int((rAuxVariables.X_max-rAuxVariables.X_min)/(AverageElementLength));
        if(rAuxVariables.NRows <= 0) rAuxVariables.NRows = 1;
        if(rAuxVariables.NColumns <= 0) rAuxVariables.NColumns = 1;
        rAuxVariables.RowSize = (rAuxVariables.Y_max-rAuxVariables.Y_min)/rAuxVariables.NRows;
        rAuxVariables.ColumnSize = (rAuxVariables.X_max-rAuxVariables.X_min)/rAuxVariables.NColumns;
    }
    
}; // Class NonlocalDamage2DUtilities

} // namespace Kratos.

#endif /* KRATOS_NONLOCAL_DAMAGE_2D_UTILITIES defined */
