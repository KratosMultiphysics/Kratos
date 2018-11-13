//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_MAPPING_VARIABLES_PROCESS )
#define  KRATOS_MAPPING_VARIABLES_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "fem_to_dem_application_variables.h"

namespace Kratos
{

//Only for Triangles2D3N or Quadrilateral2D4N

class MappingVariablesProcess : public Process
{

protected:

    struct NodeNew
    {
      Element::GeometryType::PointType& rNode;
      int Row,Column;
       
      NodeNew(Element::GeometryType::PointType& NodeReference, int Row_i, int Column_j) : rNode(NodeReference)
      {
        Row = Row_i;
        Column = Column_j;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct ElementOldCell
    {
      std::vector<Element::Pointer> ElementOldVector;
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointNew
    {
      //ConstitutiveLaw::Pointer pConstitutiveLaw;
      Element::Pointer pElement;
      double X_coord,Y_coord;
      int Row,Column;
      
      GaussPointNew(Element::Pointer ElementPointer, double X, double Y ,int Row_i, int Column_j)
      {
        //pConstitutiveLaw = ConstitutiveLawPointer;
        pElement = ElementPointer;
        X_coord = X;
        Y_coord = Y;
        Row = Row_i;
        Column = Column_j;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointOld
    {
      //ConstitutiveLaw::Pointer pConstitutiveLaw;
      Element::Pointer pElement;
      double X_coord,Y_coord;
      
      GaussPointOld(Element::Pointer ElementPointer, double X, double Y)  // TODO modify CL pointer with Element pointer
      {
        // pConstitutiveLaw = ConstitutiveLawPointer;
        pElement = ElementPointer;
        X_coord = X;
        Y_coord = Y;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointOldCell
    {
      std::vector<GaussPointOld> GaussPointOldVector;
    };
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:
  
    typedef ModelPart::ElementsContainerType ElementsArrayType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Constructor
    MappingVariablesProcess(ModelPart& r_model_part_old, ModelPart& r_model_part_new, std::string imposed_displacement, std::string Mapping_Procedure) :
     mmodel_part_old(r_model_part_old), mmodel_part_new(r_model_part_new)
    {
        mImposedDisplacement = imposed_displacement;
        mMappingProcedure    = Mapping_Procedure;
    }
    
    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~MappingVariablesProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Main Function
    void Execute()
    {
        double X_max = 0.0, X_min = 0.0, Y_max = 0.0, Y_min = 0.0;
    
        double DamageThreshold = 0.0, CharacteristicLength = 0.0;
    
        this->Initialize(X_max, X_min, Y_max, Y_min, CharacteristicLength, DamageThreshold);
    
        this->NodalDisplacementsMapping(X_max, X_min, Y_max, Y_min);
    
        this->GaussPointStateVariableMapping(X_max, X_min, Y_max, Y_min, CharacteristicLength, DamageThreshold);
    
        this->TransferProcessInfoVariables();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mmodel_part_old;
    ModelPart& mmodel_part_new;
    
    std::string mImposedDisplacement;
	std::string mMappingProcedure;
	double mAverageElementLength;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(double& rX_max, double& rX_min, double& rY_max,double& rY_min, double& rCharacteristicLength, double& rDamageThreshold)
    {
        rX_max = mmodel_part_old.NodesBegin()->X0();
        rX_min = rX_max;
        rY_max = mmodel_part_old.NodesBegin()->Y0();
        rY_min = rY_max;

        double X_me, Y_me;

        for(ModelPart::NodeIterator i = mmodel_part_old.NodesBegin(); i != mmodel_part_old.NodesEnd(); ++i)
        {
            X_me = i->X0();
            Y_me = i->Y0();

            if( X_me > rX_max )  rX_max = X_me;
            else if( X_me < rX_min )  rX_min = X_me;

            if( Y_me > rY_max )  rY_max = Y_me;
            else if( Y_me < rY_min )  rY_min = Y_me;

            //Move old mesh to the original position (to work with both meshes in the reference state)
            (i)->X() = (i)->X0();
            (i)->Y() = (i)->Y0();
            (i)->Z() = (i)->Z0();
        }
        
        //rCharacteristicLength = (*(mmodel_part_old.Elements().ptr_begin()))->GetProperties()[CHARACTERISTIC_LENGTH];
        //rDamageThreshold = (*(mmodel_part_old.Elements().ptr_begin()))->GetProperties()[DAMAGE_THRESHOLD];
        rDamageThreshold = (*(mmodel_part_old.Elements().ptr_begin()))->GetValue(STRESS_THRESHOLD);
    }
  
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void NodalDisplacementsMapping(const double& X_max,const double& X_min,const double& Y_max,const double& Y_min)
    {
        double AverageElementLength = 0.0;

        for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
        {
            AverageElementLength += (*it)->GetGeometry().Length();
        }
        AverageElementLength  = AverageElementLength / mmodel_part_old.NumberOfElements();
		mAverageElementLength = AverageElementLength;

    
        int NRows = int((Y_max - Y_min) / AverageElementLength);
        int NColumns = int((X_max - X_min) / AverageElementLength);
    
        double RowSize = (Y_max - Y_min) / NRows;
        double ColumnSize = (X_max - X_min) / NColumns;

        ElementOldCell** ElementOldMatrix;
        ElementOldMatrix = new ElementOldCell*[NRows];
        for(int i = 0; i < NRows; i++)  ElementOldMatrix[i] = new ElementOldCell[NColumns];
                
        double X_left, X_right, Y_top, Y_bot, X_me, Y_me;
        int Column_left, Column_right, Row_top, Row_bot;
        
        // Locate Old Elements in Cells
        for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
        {

			//if ((*it)->Id() == 10)
			//{
			//	KRATOS_WATCH((*it)->GetValue(STRESS_THRESHOLD))
			//	KRATOS_WATCH((*it)->GetValue(STRESS_THRESHOLD))
			//}



            X_left = (*it)->GetGeometry().GetPoint(0).X0();
            X_right = X_left;
            Y_top = (*it)->GetGeometry().GetPoint(0).Y0();
            Y_bot = Y_top;
            for(unsigned int i = 1; i < (*it)->GetGeometry().PointsNumber(); i++)
            {
                X_me = (*it)->GetGeometry().GetPoint(i).X0();
                Y_me = (*it)->GetGeometry().GetPoint(i).Y0();
                
                if( X_me > X_right )  X_right = X_me;
                else if( X_me < X_left )  X_left = X_me;
              
                if( Y_me > Y_top ){Y_top = Y_me;}
                else if( Y_me < Y_bot ){Y_bot = Y_me;}
            }
            
            Column_left = int((X_left - X_min) / ColumnSize);
            Column_right = int((X_right - X_min) / ColumnSize);
            Row_top = int((Y_max - Y_top) / RowSize);
            Row_bot = int((Y_max - Y_bot) / RowSize);

            if(Column_left == NColumns)  Column_left = NColumns - 1;
            if(Column_right == NColumns)  Column_right = NColumns - 1;
            if(Row_top == NRows)  Row_top = NRows - 1;
            if(Row_bot == NRows)  Row_bot = NRows - 1;

            for(int i = Row_top; i <= Row_bot; i++)
            {
                for(int j = Column_left; j<= Column_right; j++)
                {
                    ElementOldMatrix[i][j].ElementOldVector.push_back((*it));
                }
            }
        }
    
        int Row, Column;
        std::vector<NodeNew*> NodeNewVector;
        
        //Locate New Nodes in Cells
        for(ModelPart::NodeIterator i = mmodel_part_new.NodesBegin(); i != mmodel_part_new.NodesEnd(); ++i)
        {
            X_me = i->X0();
            Y_me = i->Y0();

            Row = int((Y_max - Y_me) / RowSize);
            Column = int((X_me - X_min) / ColumnSize);

            if(Column == NColumns) Column = NColumns - 1;
            if(Row == NRows) Row = NRows - 1;
            
            NodeNewVector.push_back(new NodeNew((*i), Row, Column));
        }

        Element::Pointer pElementOld;
        Element::GeometryType::CoordinatesArrayType NodeLocalCoordinates;
        Vector ElementDisplacements,ElementVelocities,ElementAccelerations, ElementShapeFunctions;
        double Tolerance = 1e-4;
        bool IsInside;
        
        //Locate new nodes inside old elements and interpolate displacements.  
        ElementDisplacements = ZeroVector(3);
        ElementVelocities    = ZeroVector(3);
        ElementAccelerations = ZeroVector(3);

        for(unsigned int i = 0; i < NodeNewVector.size(); i++)
        {
            Element::GeometryType::PointType& NodeNew_i = NodeNewVector[i]->rNode;
            Row = NodeNewVector[i] -> Row;
            Column = NodeNewVector[i] -> Column;
            const Element::GeometryType::CoordinatesArrayType NodeGlobalCoordinates = NodeNew_i.Coordinates(); //Coordinates of new nodes are still in the original position
            IsInside = false;
            
            for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
            {
                pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                IsInside = pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates, NodeLocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                if(IsInside)  break;                                                                        
            }
    
            if(IsInside == false) //TODO: cal??
            {
                for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
                {
                    pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                    pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates,NodeLocalCoordinates);
                    if(((NodeLocalCoordinates[0]+Tolerance)>=0)&&((NodeLocalCoordinates[1]+Tolerance)>=0)&&((NodeLocalCoordinates[1]-Tolerance)<=(1-NodeLocalCoordinates[0])))  break;
                }
            }
    
            ElementShapeFunctions = this->TriangleShapeFunctions(NodeLocalCoordinates[0], NodeLocalCoordinates[1]);
    
            if( (NodeNew_i.pGetDof(DISPLACEMENT_X))->IsFixed() == false )
            {
                ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[0];
                ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[0];
                ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[0];
                NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = inner_prod(ElementShapeFunctions,ElementDisplacements);

                if (mmodel_part_old.GetProcessInfo()[IS_DYNAMIC] == 1)  // Mapping of velocities and accelerations
                {
                    ElementVelocities[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(VELOCITY)[0];
                    ElementVelocities[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(VELOCITY)[0];
                    ElementVelocities[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(VELOCITY)[0];
                    NodeNew_i.FastGetSolutionStepValue(VELOCITY)[0] = inner_prod(ElementShapeFunctions,ElementVelocities);

                    ElementAccelerations[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(ACCELERATION)[0];
                    ElementAccelerations[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(ACCELERATION)[0];
                    ElementAccelerations[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(ACCELERATION)[0];
                    NodeNew_i.FastGetSolutionStepValue(ACCELERATION)[0] = inner_prod(ElementShapeFunctions,ElementAccelerations);
                }
            }
            // else if( mImposedDisplacement == "Linearly_Incremented" )
            //     NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = mmodel_part_old.GetProcessInfo()[TIME_STEPS] * NodeNew_i.FastGetSolutionStepValue(IMPOSED_DISPLACEMENT)[0];
            
            if( (NodeNew_i.pGetDof(DISPLACEMENT_Y))->IsFixed() == false )
            {
                ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[1];
                ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[1];
                ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[1];
                NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[1] = inner_prod(ElementShapeFunctions,ElementDisplacements);

                if (mmodel_part_old.GetProcessInfo()[IS_DYNAMIC] == 1)  // Mapping of velocities and accelerations
                {
                    ElementVelocities[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(VELOCITY)[1];
                    ElementVelocities[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(VELOCITY)[1];
                    ElementVelocities[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(VELOCITY)[1];
                    NodeNew_i.FastGetSolutionStepValue(VELOCITY)[1] = inner_prod(ElementShapeFunctions,ElementVelocities);

                    ElementAccelerations[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(ACCELERATION)[1];
                    ElementAccelerations[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(ACCELERATION)[1];
                    ElementAccelerations[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(ACCELERATION)[1];
                    NodeNew_i.FastGetSolutionStepValue(ACCELERATION)[1] = inner_prod(ElementShapeFunctions,ElementAccelerations);
                }
            }
        }
       
    
        //Deallocate memory
        for(int i=0; i < NRows; i++)
            delete[] ElementOldMatrix[i];
        delete[] ElementOldMatrix;
    
        for(unsigned int i = 0; i < NodeNewVector.size(); i++)
            delete NodeNewVector[i];
    
        std::cout << "-- Nodal Values Mapped --" << std::endl;
    
    } //NodalVariablesMapping_End

    //------------------------------------------------------------------------------------

    Vector TriangleShapeFunctions(const double& rx_local, const double& ry_local)
    {
        Vector TSF = ZeroVector(3);
		TSF[0] = 1 - rx_local - ry_local;
        TSF[1] = rx_local;
        TSF[2] = ry_local;
        return TSF;    
    }

    //------------------------------------------------------------------------------------
    
    Vector QuadrilateralShapeFunctions(const double& rx_local, const double& ry_local)
    {
        Vector QSF = ZeroVector(4);
		QSF[0] = (1 - rx_local - ry_local + rx_local*ry_local) / 4;
		QSF[1] = (1 + rx_local - ry_local - rx_local*ry_local) / 4;
		QSF[2] = (1 + rx_local + ry_local + rx_local*ry_local) / 4;
		QSF[3] = (1 - rx_local + ry_local - rx_local*ry_local) / 4;
        return QSF;    
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GaussPointStateVariableMapping(const double& X_max, const double& X_min, const double& Y_max, const double& Y_min, const double& rCharacteristicLength, const double& rDamageThreshold)
    {
        int Row, Column;
        double X_me, Y_me, X_other, Y_other, Distance;
        Element::Pointer Me;
        Element::Pointer Other;
        Vector Trian_GPLocalCoord = ZeroVector(3);
        Element::GeometryType::CoordinatesArrayType GPGlobalCoord;
        std::vector<GaussPointNew*> GaussPointNewVector;

        //int NRows = int((Y_max - Y_min) / (rCharacteristicLength * 2));
        //int NColumns = int((X_max - X_min) / (rCharacteristicLength * 2));
    
        //double RowSize = (Y_max - Y_min) / NRows;
        //double ColumnSize = (X_max - X_min) / NColumns;

        // if it's too big, Nncol could be 0 -> error
		int NRows    = int((Y_max - Y_min) / (22.5 * mAverageElementLength));
		int NColumns = int((X_max - X_min) / (22.5 * mAverageElementLength));

		double RowSize = (Y_max - Y_min) / NRows;
		double ColumnSize = (X_max - X_min) / NColumns;

        GaussPointOldCell** pGaussPointOldMatrix;
		pGaussPointOldMatrix = new GaussPointOldCell*[NRows];
		for (int i = 0;i < NRows;i++)  pGaussPointOldMatrix[i] = new GaussPointOldCell[NColumns];
        
        
        // GP local coordinates
        Trian_GPLocalCoord[0] = 1 / 3;
        Trian_GPLocalCoord[1] = 1 / 3;
        const Element::GeometryType::CoordinatesArrayType GPLocalCoord = Trian_GPLocalCoord;
        
        // Locate old elements in cells
        for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
        {
            //(*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);

			Vector GPGlobalCoord;
			CalculateGlobalBaricenterCoordinates((*it), GPGlobalCoord);
            X_me = GPGlobalCoord[0];
            Y_me = GPGlobalCoord[1];

    //        if ((*it)->Id() == 2)
    //        {
    //  			KRATOS_WATCH((*it)->GetValue(STRESS_THRESHOLD))
    //            KRATOS_WATCH(Y_me)
				//KRATOS_WATCH(GPGlobalCoord[0])
				//KRATOS_WATCH(GPGlobalCoord[1])
    //            KRATOS_WATCH(Y_me)
    //        }

            Row    = int((Y_max - Y_me) / RowSize);
            Column = int((X_me - X_min) / ColumnSize);

            if(Row == NRows) Row = NRows - 1;
            if(Column == NColumns) Column = NColumns - 1;

			//KRATOS_WATCH((*it)->Id())

            pGaussPointOldMatrix[Row][Column].GaussPointOldVector.push_back(GaussPointOld(*it, X_me, Y_me));

			//if ((*it)->Id() == 126)
			//{
			//	KRATOS_WATCH((*it)->GetValue(DAMAGE_ELEMENT))
			//	
			//	KRATOS_WATCH(Row)
			//	KRATOS_WATCH(Column)
			//	KRATOS_WATCH(mAverageElementLength)
			//	KRATOS_WATCH(NRows)
			//	KRATOS_WATCH(NColumns)


			//	//KRATOS_WATCH((*it)->GetValue(DAMAGE_ELEMENT))
			//}


        }


		//KRATOS_WATCH(pGaussPointOldMatrix[0][0].GaussPointOldVector.size())
		//for (int i = 0; i < pGaussPointOldMatrix[0][0].GaussPointOldVector.size(); i++)
		//{
		//	KRATOS_WATCH(pGaussPointOldMatrix[0][0].GaussPointOldVector[0].pElement->Id())
		//}
		

        // Locate New elements in cells
        for(ElementsArrayType::ptr_iterator it = mmodel_part_new.Elements().ptr_begin(); it != mmodel_part_new.Elements().ptr_end(); ++it)
        {
            // (*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);
            
            // X_me = GPGlobalCoord[0];
            // Y_me = GPGlobalCoord[1];

            Vector GPGlobalCoord;
            CalculateGlobalBaricenterCoordinates((*it), GPGlobalCoord);
            X_me = GPGlobalCoord[0];
            Y_me = GPGlobalCoord[1];

            Row = int((Y_max - Y_me) / RowSize);
            Column = int((X_me - X_min) / ColumnSize);

            if(Row == NRows) Row = NRows - 1;
            if(Column == NColumns) Column = NColumns - 1;
            GaussPointNewVector.push_back(new GaussPointNew(*it, X_me, Y_me, Row, Column));

			//KRATOS_WATCH((*it)->Id())
			//if ((*it)->Id() == 7)
			//{
			//	KRATOS_WATCH(Row)
			//	KRATOS_WATCH(Column)
			//	//KRATOS_WATCH(Column)
			//	//KRATOS_WATCH(Column)
			//}


        }  
		
		// testsss
		//KRATOS_WATCH(GaussPointNewVector.size())
		//for (int i = 0; i < GaussPointNewVector.size();i++)
		//{
		//	if (GaussPointNewVector[i]->Row == 0 && GaussPointNewVector[i]->Column == 0)
		//	{
		//		Element::Pointer elem = GaussPointNewVector[i]->pElement;
		//		KRATOS_WATCH(elem->Id())
		//	}
		//}

        if (mMappingProcedure == "Non_Local_Average")
        {
            //Transfer state variables from old Gauss points to new Gauss Points (nonlocal average)
            double IntegrationCoefficient;     //StateVariable, Numerator, WeightingFunctionDenominator;
            double DamageVariable, StressThreshold, DamageNumerator, ThresholdNumerator, WeightingFunctionDenominator;

            // Analize each new element in its own cell 
            for(unsigned int i = 0; i < GaussPointNewVector.size(); i++)
            {
                Me   = GaussPointNewVector[i]->pElement;
                X_me = GaussPointNewVector[i]->X_coord;
                Y_me = GaussPointNewVector[i]->Y_coord;
                Row  = GaussPointNewVector[i]->Row;
                Column = GaussPointNewVector[i]->Column;
                DamageNumerator = 0.0;
                ThresholdNumerator  = 0.0;
                WeightingFunctionDenominator = 0.0;
                
                // Loop over old elements inside the cell
                for(unsigned int j = 0; j < pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size(); j++)
                {
                    Other   = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].pElement;
                    X_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].X_coord;
                    Y_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].Y_coord;
            
                    Distance = sqrt((X_other - X_me) * (X_other - X_me) + (Y_other - Y_me) * (Y_other - Y_me));

                    if(Distance <= rCharacteristicLength) //TODO: es podria calcular amb tots els de la teva cel·la...
                    {
                        IntegrationCoefficient = Other->GetValue(INTEGRATION_COEFFICIENT);
                        DamageVariable = Other->GetValue(DAMAGE_ELEMENT);
                        StressThreshold = Other->GetValue(STRESS_THRESHOLD);
                        
                        // Mapping of the damage and the threshold GP variables
                        DamageNumerator    += IntegrationCoefficient * exp(- 4 * Distance * Distance / (rCharacteristicLength * rCharacteristicLength)) * DamageVariable;
                        ThresholdNumerator += IntegrationCoefficient * exp(- 4 * Distance * Distance / (rCharacteristicLength * rCharacteristicLength)) * StressThreshold;
                        WeightingFunctionDenominator += IntegrationCoefficient * exp(- 4 * Distance * Distance / (rCharacteristicLength * rCharacteristicLength));
                    }
                }
                
                //Search in adjacent cells
                if (sqrt((X_min + ColumnSize*(Column + 1) - X_me)*(X_min + ColumnSize*(Column + 1) - X_me) + (Y_max - RowSize*(Row + 1) - Y_me)*(Y_max - RowSize*(Row + 1) - Y_me)) < rCharacteristicLength)
                {
                    if ((Row + 1) < NRows)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row + 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if ((Column + 1) < NColumns)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column + 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if (((Row + 1) < NRows) && ((Column + 1) < NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row + 1][Column + 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                }
                else if (sqrt((X_min + ColumnSize*(Column)-X_me)*(X_min + ColumnSize*(Column)-X_me) + (Y_max - RowSize*(Row + 1) - Y_me)*(Y_max - RowSize*(Row + 1) - Y_me)) < rCharacteristicLength)
                {
                    if ((Row + 1) < NRows)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row + 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if ((Column - 1) >= 0)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column - 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if (((Row + 1) < NRows) && ((Column - 1) >= 0))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row + 1][Column - 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                }
                else if (sqrt((X_min + ColumnSize*(Column)-X_me)*(X_min + ColumnSize*(Column)-X_me) + (Y_max - RowSize*(Row)-Y_me)*(Y_max - RowSize*(Row)-Y_me)) < rCharacteristicLength)
                {
                    if ((Row - 1) >= 0)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row - 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if ((Column - 1) >= 0)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column - 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if (((Row - 1) >= 0) && ((Column - 1) >= 0))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row - 1][Column - 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                }
                else if (sqrt((X_min + ColumnSize*(Column + 1) - X_me)*(X_min + ColumnSize*(Column + 1) - X_me) + (Y_max - RowSize*(Row)-Y_me)*(Y_max - RowSize*(Row)-Y_me)) < rCharacteristicLength)
                {
                    if ((Row - 1) >= 0)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row - 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if ((Column + 1) < NColumns)  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column + 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    if (((Row - 1) >= 0) && ((Column + 1) < NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row - 1][Column + 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                }
                else
                {
                    if ((int((X_me - X_min + rCharacteristicLength) / ColumnSize) > Column) && ((Column + 1) < NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column + 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    else if ((int((X_me - X_min - rCharacteristicLength) / ColumnSize) < Column) && ((Column - 1) >= 0))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row][Column - 1], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
            
                    if ((int((Y_max - Y_me + rCharacteristicLength) / RowSize) > Row) && ((Row + 1) < NRows))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row + 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                    else if ((int((Y_max - Y_me - rCharacteristicLength) / RowSize) < Row) && ((Row - 1) >= 0))  this->SearchInAdjacentCell(GaussPointNewVector[i], pGaussPointOldMatrix[Row - 1][Column], DamageNumerator,ThresholdNumerator, WeightingFunctionDenominator, rCharacteristicLength);
                }
        
                if (fabs(WeightingFunctionDenominator) < 1e-15) 
                { // TODO-> CHECK WHY
                    StressThreshold = rDamageThreshold;
                    DamageVariable = 0.0;
                }
                else 
                {
                    DamageVariable  = DamageNumerator / WeightingFunctionDenominator;
                    StressThreshold = ThresholdNumerator / WeightingFunctionDenominator;
                }

                Me->SetValue(DAMAGE_ELEMENT, DamageVariable);
                Me->SetValue(STRESS_THRESHOLD, StressThreshold);
            }

        }

        else // Closest_Point_Transfer
        {
            // Loop over new elements
            for(unsigned int i = 0; i < GaussPointNewVector.size(); i++)
            {
                Me     = GaussPointNewVector[i]->pElement;
                X_me   = GaussPointNewVector[i]->X_coord;
                Y_me   = GaussPointNewVector[i]->Y_coord;
                Row    = GaussPointNewVector[i]->Row;
                Column = GaussPointNewVector[i]->Column;

				//KRATOS_WATCH(Me->Id())
				/*if (Me->Id() == 7)
				{
					KRATOS_WATCH(Row)
					KRATOS_WATCH(Column)
					KRATOS_WATCH(Column)
				}*/


                double* Distance  = new double[pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size()];
                double* Damage    = new double[pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size()];
                double* Threshold = new double[pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size()];
                double damage, threshold, MinDistance;
                
                // Loop over old elements inside the cell
                for (unsigned int j = 0; j < pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size(); j++)
                {
                    Other   = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].pElement;
                    X_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].X_coord;
                    Y_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].Y_coord;
            
                    Distance[j]  = sqrt((X_other - X_me) * (X_other - X_me) + (Y_other - Y_me) * (Y_other - Y_me));
                    Damage[j]    = Other->GetValue(DAMAGE_ELEMENT);
                    Threshold[j] = Other->GetValue(STRESS_THRESHOLD);

					//if (Other->Id() == 10)
					//{
					//	KRATOS_WATCH(Threshold[j])
					//	KRATOS_WATCH(Threshold[j])
					//}
					
                    //GPNodeLocalCoordinates[0] = X_me;
                    //GPNodeLocalCoordinates[1] = Y_me;


                    //KRATOS_WATCH(GPNodeGlobalCoordinates)
                    //GPNodeGlobalCoordinates = (X_me, Y_me, 0.0); not working

                    //KRATOS_WATCH(GPIsInside)

                    std::string flag = "any_nodes_is_in";
					//std::string flag = "GP_is_in";

                    if (flag == "GP_is_in")
                    {
                        // check if new GP is inside element
                        bool GPIsInside = false;
                        Element::GeometryType::CoordinatesArrayType GPNodeGlobalCoordinates; 
                        Element::GeometryType::CoordinatesArrayType GPNodeLocalCoordinates;
                        GPNodeGlobalCoordinates = Me->GetGeometry()[0].Coordinates();
                        GPNodeGlobalCoordinates[0] = X_me;
                        GPNodeGlobalCoordinates[1] = Y_me;
                        GPIsInside = Other->GetGeometry().IsInside(GPNodeGlobalCoordinates,GPNodeLocalCoordinates);

                        bool condition_is_active = true;
                        if ((Other)->IsDefined(ACTIVE))
                        {
                            condition_is_active = (Other)->Is(ACTIVE);
                        }
    
                        if (condition_is_active == false) // the inactive old elements have damage 0
                        {
                            Damage[j] = 0.97;
                            Threshold[j] = Other->GetValue(INITIAL_THRESHOLD); // does not work TODO
                            Distance[j] = 1.0e24;  // avoid to map this value 
                        }
                        if (condition_is_active == false && GPIsInside == true)
                        {
                            Me->Set(ACTIVE, false);
                        }
                    }

                    else // "any_nodes_is_in"
                    {
                        Element::GeometryType::CoordinatesArrayType PointNodeGlobalCoordinates;
                        Element::GeometryType::CoordinatesArrayType PointNodeLocalCoordinates;

                        bool any_node_is_inside = false;
                        bool condition_is_active = true;
                        if ((Other)->IsDefined(ACTIVE))
                        {
                            condition_is_active = (Other)->Is(ACTIVE);
                        }
    
                        if (condition_is_active == false) // the inactive old elements have damage 0
                        {
                            Damage[j] = 0.97;
							Threshold[j] = Other->GetValue(INITIAL_THRESHOLD);
							Distance[j] = 1.0e24;  // avoid to map this value 
                        }
						int num_nodes = 0;
                        for (int node = 0; node < 3; node++)
                        {
                            PointNodeGlobalCoordinates = Me->GetGeometry()[node].Coordinates();
                            any_node_is_inside =  Other->GetGeometry().IsInside(PointNodeGlobalCoordinates,PointNodeLocalCoordinates);
                            if (any_node_is_inside == true) num_nodes++;
                        }

                        if (condition_is_active == false && num_nodes >= 1)
                        {
                            Me->Set(ACTIVE, false);
                        }
                    }

                } // end loop over old elements in cell

                // Let's see if Me element is inside any inactive old element in the neighbour cells
                //***********************************************************************************
                //***********************************************************************************
                bool indicator_active = true;
                if ((Me)->IsDefined(ACTIVE))
                {
                    indicator_active = (Me)->Is(ACTIVE); 
                }
                
                // check if active
                if (indicator_active == true) // if it is inactive no need to searh other cells
                {
                    if(Row + 1 < NRows)
                    {
                        for (unsigned int j = 0; j < pGaussPointOldMatrix[Row + 1][Column].GaussPointOldVector.size(); j++)
                        {
                            Other   = pGaussPointOldMatrix[Row + 1][Column].GaussPointOldVector[j].pElement;
                            X_other = pGaussPointOldMatrix[Row + 1][Column].GaussPointOldVector[j].X_coord;
                            Y_other = pGaussPointOldMatrix[Row + 1][Column].GaussPointOldVector[j].Y_coord;
    
                            Element::GeometryType::CoordinatesArrayType PointNodeGlobalCoordinates;
                            Element::GeometryType::CoordinatesArrayType PointNodeLocalCoordinates;
    
                            bool any_node_is_inside = false;
                            bool condition_is_active = true;
                            if ((Other)->IsDefined(ACTIVE))
                            {
                                condition_is_active = (Other)->Is(ACTIVE);
                            }
                            int num_nodes = 0;
                            for (int node = 0; node < 3; node++)
                            {
                                PointNodeGlobalCoordinates = Me->GetGeometry()[node].Coordinates();
                                any_node_is_inside =  Other->GetGeometry().IsInside(PointNodeGlobalCoordinates,PointNodeLocalCoordinates);
                                if (any_node_is_inside == true) num_nodes++;
                            }
                            if (condition_is_active == false && num_nodes >= 1)
                            {
                                Me->Set(ACTIVE, false);
                            }
                        }
                    }
                    // check if active
                    if ((Me)->IsDefined(ACTIVE))
                    {
                        indicator_active = (Me)->Is(ACTIVE); // if it is inactive no need to searh other cells
                    }

                    if (indicator_active == true)
                    {
                        if(Row - 1 >= 0)
                        {
                            for (unsigned int j = 0; j < pGaussPointOldMatrix[Row - 1][Column].GaussPointOldVector.size(); j++)
                            {
                                Other   = pGaussPointOldMatrix[Row - 1][Column].GaussPointOldVector[j].pElement;
                                X_other = pGaussPointOldMatrix[Row - 1][Column].GaussPointOldVector[j].X_coord;
                                Y_other = pGaussPointOldMatrix[Row - 1][Column].GaussPointOldVector[j].Y_coord;
        
                                Element::GeometryType::CoordinatesArrayType PointNodeGlobalCoordinates;
                                Element::GeometryType::CoordinatesArrayType PointNodeLocalCoordinates;
        
                                bool any_node_is_inside = false;
                                bool condition_is_active = true;
                                if ((Other)->IsDefined(ACTIVE))
                                {
                                    condition_is_active = (Other)->Is(ACTIVE);
                                }
                                int num_nodes = 0;
                                for (int node = 0; node < 3; node++)
                                {
                                    PointNodeGlobalCoordinates = Me->GetGeometry()[node].Coordinates();
                                    any_node_is_inside =  Other->GetGeometry().IsInside(PointNodeGlobalCoordinates,PointNodeLocalCoordinates);
                                    if (any_node_is_inside == true) num_nodes++;
                                }
                                if (condition_is_active == false && num_nodes >= 1)
                                {
                                    Me->Set(ACTIVE, false);
                                }
                            }

                            // check if active
                            if ((Me)->IsDefined(ACTIVE))
                            {
                                indicator_active = (Me)->Is(ACTIVE); // if it is inactive no need to searh other cells
                            }

                            if (indicator_active == true)
                            {
                                if(Column - 1 >= 0)
                                {
                                    for (unsigned int j = 0; j < pGaussPointOldMatrix[Row][Column - 1].GaussPointOldVector.size(); j++)
                                    {
                                        Other   = pGaussPointOldMatrix[Row][Column - 1].GaussPointOldVector[j].pElement;
                                        X_other = pGaussPointOldMatrix[Row][Column - 1].GaussPointOldVector[j].X_coord;
                                        Y_other = pGaussPointOldMatrix[Row][Column - 1].GaussPointOldVector[j].Y_coord;
                
                                        Element::GeometryType::CoordinatesArrayType PointNodeGlobalCoordinates;
                                        Element::GeometryType::CoordinatesArrayType PointNodeLocalCoordinates;
                
                                        bool any_node_is_inside = false;
                                        bool condition_is_active = true;
                                        if ((Other)->IsDefined(ACTIVE))
                                        {
                                            condition_is_active = (Other)->Is(ACTIVE);
                                        }
                                        int num_nodes = 0;
                                        for (int node = 0; node < 3; node++)
                                        {
                                            PointNodeGlobalCoordinates = Me->GetGeometry()[node].Coordinates();
                                            any_node_is_inside =  Other->GetGeometry().IsInside(PointNodeGlobalCoordinates,PointNodeLocalCoordinates);
                                            if (any_node_is_inside == true) num_nodes++;
                                        }
                                        if (condition_is_active == false && num_nodes >= 1)
                                        {
                                            Me->Set(ACTIVE, false);
                                        }
                                    }

                                    // check if active
                                    if ((Me)->IsDefined(ACTIVE))
                                    {
                                        indicator_active = (Me)->Is(ACTIVE); // if it is inactive no need to searh other cells
                                    }
                                    
                                    if (indicator_active == true)
                                    {
                                        if(Column + 1 < NColumns)
                                        {
                                            for (unsigned int j = 0; j < pGaussPointOldMatrix[Row][Column + 1].GaussPointOldVector.size(); j++)
                                            {
                                                Other   = pGaussPointOldMatrix[Row][Column + 1].GaussPointOldVector[j].pElement;
                                                X_other = pGaussPointOldMatrix[Row][Column + 1].GaussPointOldVector[j].X_coord;
                                                Y_other = pGaussPointOldMatrix[Row][Column + 1].GaussPointOldVector[j].Y_coord;
                        
                                                Element::GeometryType::CoordinatesArrayType PointNodeGlobalCoordinates;
                                                Element::GeometryType::CoordinatesArrayType PointNodeLocalCoordinates;
                        
                                                bool any_node_is_inside = false;
                                                bool condition_is_active = true;
                                                if ((Other)->IsDefined(ACTIVE))
                                                {
                                                    condition_is_active = (Other)->Is(ACTIVE);
                                                }
                                                int num_nodes = 0;
                                                for (int node = 0; node < 3; node++)
                                                {
                                                    PointNodeGlobalCoordinates = Me->GetGeometry()[node].Coordinates();
                                                    any_node_is_inside =  Other->GetGeometry().IsInside(PointNodeGlobalCoordinates,PointNodeLocalCoordinates);
                                                    if (any_node_is_inside == true) num_nodes++;
                                                }
                                                if (condition_is_active == false && num_nodes >= 1)
                                                {
                                                    Me->Set(ACTIVE, false);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } // end searching neighbour cells
                //***********************************************************************************
                //***********************************************************************************
                

                // Map the values of the closest gauss point on the current cell
                bool condition_is_active = true;
                if ((Me)->IsDefined(ACTIVE))
                {
                    condition_is_active = (Me)->Is(ACTIVE);
                }

                if (condition_is_active)
                {
                    MinDistance = Distance[0];
                    damage      = Damage[0];
                    threshold   = Threshold[0];
                    //KRATOS_WATCH(damage)

                    // Select the closest point old element
                    for (int elem = 1; elem < pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size(); elem++)
                    {
                        if (Distance[elem] < MinDistance)
                        {
                            MinDistance = Distance[elem];
                            damage      = Damage[elem];
                            threshold   = Threshold[elem];
                        }
                    }
                    // Set GP values
                    Me->SetValue(DAMAGE_ELEMENT, damage);
                    Me->SetValue(STRESS_THRESHOLD, threshold);
                }


                delete[] Distance;
                delete[] Damage;
                delete[] Threshold;
            }
        }

        // Deallocate memory
		for (int i = 0;i < NRows;i++)
			delete[] pGaussPointOldMatrix[i];
		delete[] pGaussPointOldMatrix;
    
		for (unsigned int i = 0;i < GaussPointNewVector.size(); i++)
			delete GaussPointNewVector[i];
    
        std::cout << "-- Gauss Point State Variables Mapped --" << std::endl;




    
    } //GaussPointStateVariableMapping_End

    //------------------------------------------------------------------------------------
    
    void SearchInAdjacentCell(const GaussPointNew* pGaussPointNew, const GaussPointOldCell& NeighbourCell, double& rDamageNumerator,double& rThresholdNumerator,
         double& rWeightingFunctionDenominator, const double& rCharacteristicLength)
    {
        double X_me = pGaussPointNew->X_coord;
        double Y_me = pGaussPointNew->Y_coord;
    
        Element::Pointer Other;
		double X_other, Y_other;
		double IntegrationCoefficient, DamageVariable, StressThreshold;
		double Distance;
                    
		for (unsigned int j = 0; j < NeighbourCell.GaussPointOldVector.size(); j++)
        {
			Other = NeighbourCell.GaussPointOldVector[j].pElement;
            X_other = NeighbourCell.GaussPointOldVector[j].X_coord;
            Y_other = NeighbourCell.GaussPointOldVector[j].Y_coord;

			Distance = sqrt((X_other - X_me)*(X_other - X_me) + (Y_other - Y_me)*(Y_other - Y_me));
      
            if(Distance <= rCharacteristicLength) //TODO: es podria calcular amb tots els de les cel·les vehines...
            {
                IntegrationCoefficient = Other->GetValue(INTEGRATION_COEFFICIENT);
                DamageVariable = Other->GetValue(DAMAGE_ELEMENT);
                StressThreshold = Other->GetValue(STRESS_THRESHOLD);

				rDamageNumerator += IntegrationCoefficient*exp(-4 * Distance*Distance / (rCharacteristicLength*rCharacteristicLength))*DamageVariable;
				rThresholdNumerator += IntegrationCoefficient * exp(-4 * Distance * Distance / (rCharacteristicLength * rCharacteristicLength)) * StressThreshold;
				rWeightingFunctionDenominator += IntegrationCoefficient*exp(-4 * Distance*Distance / (rCharacteristicLength*rCharacteristicLength));
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void TransferProcessInfoVariables()
    {
        //Arc-length parameters
        /*mmodel_part_new.GetProcessInfo()[LOAD_FACTOR] = mmodel_part_old.GetProcessInfo()[LOAD_FACTOR];
        mmodel_part_new.GetProcessInfo()[RADIUS_FACTOR] = mmodel_part_old.GetProcessInfo()[RADIUS_FACTOR];*/
        
        //To compute linearly incremented loads
        //mmodel_part_new.GetProcessInfo()[TIME_STEPS] = mmodel_part_old.GetProcessInfo()[TIME_STEPS];
    }


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateGlobalBaricenterCoordinates(Element::Pointer pElement, Vector& rCoordinates)
    {
        rCoordinates = ZeroVector(2);

        Geometry< Node < 3 > >& NodesElem = pElement->GetGeometry();
        Vector X = ZeroVector(3);  // x coord of the nodes
        Vector Y = ZeroVector(3);  // y coord of the nodes

        for (int node = 0; node < 3; node++)
        {
            X[node] = NodesElem[node].X();
            Y[node] = NodesElem[node].Y();
        }

        rCoordinates[0] = (X[0] + X[1] + X[2]) / 3.0;
        rCoordinates[1] = (Y[0] + Y[1] + Y[2]) / 3.0;
    }

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_MAPPING_VARIABLES_PROCESS defined */
