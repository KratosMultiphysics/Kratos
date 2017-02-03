//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               July 2016 $
//   Revision:            $Revision:                 1.0 $
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
    typedef NonlocalDamageUtilities::NeighbourPoint NeighbourPoint;
    typedef NonlocalDamageUtilities::NonlocalPoint NonlocalPoint;
    using NonlocalDamageUtilities::mNonlocalPointList;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    struct Utility2DVariables
    {
        double X_max, X_min, Y_max, Y_min;
        int NRows, NColumns;
        double RowSize, ColumnSize;
        
        std::vector< std::vector< std::vector<GaussPoint> > > GaussPointCellMatrix;
    };
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:

    KRATOS_CLASS_POINTER_DEFINITION( NonlocalDamage2DUtilities );

    /// Default Constructor
    NonlocalDamage2DUtilities() : NonlocalDamageUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~NonlocalDamage2DUtilities() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void SearchGaussPointsNeighbours (Parameters* pParameters, ModelPart& rModelPart)
    {
        std::cout << "Starting non-local search of neighbours ..." << std::endl;
        
        // Define necessary variables
        Utility2DVariables AuxVariables;

        //Set GaussPoints inside CellMatrix
        this->InitializeNonlocalSearch(AuxVariables,pParameters,rModelPart);
                
        this->SearchNeighbours(AuxVariables,pParameters,rModelPart);
        
        std::cout << "... search of neighbours completed." << std::endl;
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
        GaussPoint MyGaussPoint;
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
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
            #pragma omp parallel for private(MyGaussPoint,MyIntegrationMethod,AuxLocalCoordinates)
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
                    rGeom.GlobalCoordinates(MyGaussPoint.Coordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    // MyGaussPoint Weight
                    MyGaussPoint.Weight = detJContainer[GPoint]*IntegrationPoints[GPoint].Weight();

                    // MyGaussPoint ConstitutiveLaw Pointer
                    MyGaussPoint.pConstitutiveLaw = ConstitutiveLawVector[GPoint];

                    // MyGaussPoint Row and Column
                    Row = int((rAuxVariables.Y_max-MyGaussPoint.Coordinates[1])/rAuxVariables.RowSize);
                    Column = int((MyGaussPoint.Coordinates[0]-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    #pragma omp critical
                    {
                        rAuxVariables.GaussPointCellMatrix[Row][Column].push_back(MyGaussPoint);
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
        NeighbourPoint MyNeighbourPoint;
        array_1d<double,3> MyCoordinates;
        array_1d<double,3> AuxLocalCoordinates;
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        Parameters& rParameters = *pParameters;
        double CharacteristicLength = rParameters["characteristic_length"].GetDouble();
        unsigned int NumBodySubModelParts = rParameters["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int i = 0; i < NumBodySubModelParts; i++)
        {
            ModelPart& BodySubModelPart = rModelPart.GetSubModelPart(rParameters["body_domain_sub_model_part_list"][i].GetString());

            int NElems = static_cast<int>(BodySubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = BodySubModelPart.ElementsBegin();
            
            // Loop through all body elements
            #pragma omp parallel for private(MyNeighbourPoint,MyCoordinates,AuxLocalCoordinates)
            for(int j = 0; j < NElems; j++)
            {
                ModelPart::ElementsContainerType::iterator itElem = el_begin + j;

                Element::GeometryType& rGeom = itElem->GetGeometry();
                const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(itElem->GetIntegrationMethod());
                unsigned int NumGPoints = IntegrationPoints.size();
                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(NumGPoints);
                itElem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,ConstitutiveLawVector,CurrentProcessInfo);

                // Loop through GaussPoints
                for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
                {
                    // Coordinates
                    AuxLocalCoordinates[0] = IntegrationPoints[GPoint][0];
                    AuxLocalCoordinates[1] = IntegrationPoints[GPoint][1];
                    AuxLocalCoordinates[2] = IntegrationPoints[GPoint][2];
                    rGeom.GlobalCoordinates(MyCoordinates,AuxLocalCoordinates); //Note: these are the CURRENT global coordinates

                    double X_left = MyCoordinates[0] - CharacteristicLength;
                    double X_right = MyCoordinates[0] + CharacteristicLength;
                    double Y_top = MyCoordinates[1] + CharacteristicLength;
                    double Y_bot = MyCoordinates[1] - CharacteristicLength;

                    int Column_left = int((X_left-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    int Column_right = int((X_right-rAuxVariables.X_min)/rAuxVariables.ColumnSize);
                    int Row_top = int((rAuxVariables.Y_max-Y_top)/rAuxVariables.RowSize);
                    int Row_bot = int((rAuxVariables.Y_max-Y_bot)/rAuxVariables.RowSize);

                    if(Column_left < 0) Column_left = 0;
                    if(Column_right >= rAuxVariables.NColumns) Column_right = rAuxVariables.NColumns-1;
                    if(Row_top < 0) Row_top = 0;
                    if(Row_bot >= rAuxVariables.NRows) Row_bot = rAuxVariables.NRows-1;
                    
                    NonlocalPoint MyPoint;
                    MyPoint.pConstitutiveLaw = ConstitutiveLawVector[GPoint];
                    
                    // Search GaussPoints neighbours
                    for(int k = Row_top; k <= Row_bot; k++)
                    {
                        for(int l = Column_left; l <= Column_right; l++)
                        {
                            for(unsigned int m = 0; m < rAuxVariables.GaussPointCellMatrix[k][l].size(); m++)
                            {
                                const GaussPoint& MyGaussPoint = rAuxVariables.GaussPointCellMatrix[k][l][m];

                                double Distance = sqrt((MyGaussPoint.Coordinates[0]-MyCoordinates[0])*(MyGaussPoint.Coordinates[0]-MyCoordinates[0]) +
                                                (MyGaussPoint.Coordinates[1]-MyCoordinates[1])*(MyGaussPoint.Coordinates[1]-MyCoordinates[1]));

                                if(Distance <= CharacteristicLength)
                                {
                                    MyNeighbourPoint.pConstitutiveLaw = MyGaussPoint.pConstitutiveLaw;
                                    MyNeighbourPoint.Weight = MyGaussPoint.Weight;
                                    MyNeighbourPoint.Distance = Distance;
                                    
                                    MyPoint.NeighbourPoints.push_back(MyNeighbourPoint);
                                }
                            }
                        }
                    }
                    #pragma omp critical
                    {
                        mNonlocalPointList.push_back(MyPoint);
                    }
                }
            }
        }
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
