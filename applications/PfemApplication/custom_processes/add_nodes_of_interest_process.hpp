//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:                  $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:       January 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_ADD_NODES_OF_INTEREST_PROCESS_H_INCLUDED )
#define  KRATOS_ADD_NODES_OF_INTEREST_PROCESS_H_INCLUDED


// External includes

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>


// Project includes
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//StepData: NODAL_H, NORMAL, CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) BOUNDARY, 
//          (set)     BOUNDARY(nodes), TO_ERASE(conditions), NEW_ENTITY(conditions,nodes)(set),
//          (modified)  
//          (reset)   TO_SPLIT
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Add Nodes of Interest Process
/** The process ...
    
*/

class AddNodesOfInterestProcess
  : public Process
 {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( AddNodesOfInterestProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef PointerVector<NodeType>        PointsArrayType;
  
    typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
    typedef ConditionsContainerType::iterator                    ConditionIterator;
    typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;
  
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AddNodesOfInterestProcess(ModelPart& rModelPart,
			      ModelerUtilities::MeshingParameters& rRemeshingParameters,
			      int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {     
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~AddNodesOfInterestProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the Process algorithms.
    virtual void Execute() 
    {
		KRATOS_TRY

		if( this->mEchoLevel > 0 ){
			std::cout<<" [ ADD NODES OF INTEREST : "<<std::endl;
			std::cout<<"  "<<mrModelPart.NumberOfNodes()<<" Nodes and "<<mrModelPart.NumberOfConditions()<<" Conditions "<<std::endl;
		}

		if( mrModelPart.Name() != mrRemesh.SubModelPartName )
			std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;
		
		//if (REFINE_INSERT_NODES or REFINE_ADD_NODES) & REFINE_BOUNDARY
		if( mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_INSERT_NODES) || mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES) )
		{
			if(1 == 0)
			{
				std::cout<<"+++"<<std::endl;
				std::cout<<"+++Refine conditions are both met!"<<std::endl;
				std::cout<<"+++"<<std::endl;
			}
			
			// list of new nodes
			std::vector<NodeType::Pointer> list_of_new_nodes;
			
			if(1 == 0)
			{
				std::cout<<"+++"<<std::endl;
				std::cout<<"+++Refine: test"<<mrRemesh.Refine->ConditionalMeshingNodes[0][0]<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables.size())<<std::endl;
				std::cout<<"+++Refine node variables: "<<(mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables.size())<<std::endl;
				std::cout<<"+++Refine: length: "<<mrRemesh.Refine->ConditionalMeshingNodes[0].size()<<std::endl;
				std::cout<<"+++"<<std::endl;
			}
			

			// check the list of nodes and find nodes that are too close

			//release previous nodes and remove from model part
			this->CleanPreviousNodes(mrModelPart);
			
			//generate the new nodes
			this->GenerateNewNodes(mrModelPart, list_of_new_nodes, mrRemesh.Refine->ConditionalMeshingNodes);

			//project variables to new nodes from mesh elements
			this->ProjectVariablesToNewNodes(mrModelPart, list_of_new_nodes);

			//set node variables
			for(std::vector<NodeType::Pointer>::iterator it = list_of_new_nodes.begin(); it != list_of_new_nodes.end(); it++)
				this->SetNewNodeVariables(mrModelPart, (*it));
		
			//add new nodes to model part
			this->SetNodesToModelPart(mrModelPart, list_of_new_nodes);
		
		}

		// echo output
		if( this->mEchoLevel > 0 ){
/*
			std::cout<<"   [ CONDITIONS ( total : "<<mrModelPart.NumberOfConditions()<<" ) ]"<<std::endl;
			std::cout<<"   [ NODES      ( inserted : "<<mrRemesh.Info->InsertedBoundaryNodes<<" total: "<<mrModelPart.NumberOfNodes()<<" ) ]"<<std::endl;

			if( this->mEchoLevel >=1 ){
				mrRemesh.Refine->Info.BoundaryConditionsRefined.EchoStats();
			}
*/
			std::cout<<"   ADD NODES OF INTEREST ]; "<<std::endl;
		}

		KRATOS_CATCH(" ")
    
    }
    

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {	
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }

    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }

    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }

    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "AddNodesOfInterestProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AddNodesOfInterestProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
  
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;
 
    ModelerUtilities::MeshingParameters& mrRemesh;

    ModelerUtilities mModelerUtilities;
  
    int mEchoLevel;
  
    ///@}
    ///@name Protected Operators
    ///@{

	//*******************************************************************************************
    //*******************************************************************************************
	
	void WriteOutputDataToFiles(ModelPart& rModelPart, ModelPart::NodesContainerType::iterator& i_node)
	{
		KRATOS_TRY
		
		//process info
		ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		//std::cout << CurrentProcessInfo.GetValue(IS_RESTARTED) << std::endl;
		
		//create output files for removed points
		std::ofstream outFile;
		const char* outputFileName = "ConditionalMeshingData.csv";
		
/*
		//header
		std::string headerNodes;
		std::string headerGauss;
		headerNodes += "TIME ID X Y Z ";
		
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables.size(); it++ )
		{
			headerNodes += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables.size(); it++ )
		{
			headerNodes += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables.size(); it++ )
		{
			headerNodes += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables.size(); it++ )
		{
			headerNodes += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables[it])) + " ");
		}
		
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables.size(); it++ )
		{
			headerGauss += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables.size(); it++ )
		{
			headerGauss += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables.size(); it++ )
		{
			headerGauss += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables[it])) + " ");
		}
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables.size(); it++ )
		{
			headerGauss += (std::to_string(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables[it])) + " ");
		}

		if(!(CurrentProcessInfo.GetValue(IS_RESTARTED)))
		{
			//write string into file
			outFile.open(outputFileName, std::fstream::out);
			if(outFile.is_open())
			{
				outFile << headerNodes + headerGauss << std::endl;
				outFile.close();
			}
		}
*/
		//initiaite nodal data with time step
		std::string nodalData = std::to_string(rModelPart.GetProcessInfo().GetValue(TIME)) + " ";
		std::string auxString;
		//node ID, X, Y, Z
		nodalData += (std::to_string(i_node->Id()) + " " + std::to_string(i_node->X()) + " " + std::to_string(i_node->Y()) + " " + std::to_string(i_node->Z()) + " ");
		
		//loop over nodal data double Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables.size(); it++ )
		{
			nodalData += (std::to_string(i_node->FastGetSolutionStepValue(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables[it]))) + " ");
		}
		
		//loop over nodal data array1D Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables.size(); it++ )
		{
			ConvertArray1DToString(&(i_node->FastGetSolutionStepValue(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables[it]))), auxString);
			nodalData += auxString;
		}
		
		//loop over nodal data vector Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables.size(); it++ )
		{
			ConvertVectorToString(&(i_node->FastGetSolutionStepValue(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables[it]))), auxString);
			nodalData += auxString;
		}
		
		//loop over nodal data matrix Variables
		for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables.size(); it++ )
		{
			ConvertMatrixToString(&(i_node->FastGetSolutionStepValue(*(mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables[it]))), auxString);
			nodalData += auxString;
		}
		
		
		//find neighbouring elements for Gauss point data
		WeakPointerVector< Element >& neighb_elems = i_node->GetValue(NEIGHBOUR_ELEMENTS);
		
		//process info, integration point number
		//ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin()->GetGeometry();
	    GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
	    unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );
	    
		//averaged (elemental area) array containing double, 1Darray, vector and matrix data
		//each sorted according to input
		std::vector<double> NodesDoubleVariableArray (mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables.size());
	    std::vector<array_1d<double,3> > NodesArray1DVariableArray (mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables.size());
	    std::vector<Vector> NodesVectorVariableArray (mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables.size());
	    std::vector<Matrix> NodesMatrixVariableArray (mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables.size());
		
		std::vector<double> ElementDoubleVariableArray (integration_points_number, 0);
		std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
	    std::vector<Vector> ElementVectorVariableArray (integration_points_number);
	    std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);
/*
	    std::cout << "-----------------Sizes: " << NodesDoubleVariableArray.size() << std::endl;
	    std::cout << "-----------------Sizes: " << NodesArray1DVariableArray.size() << std::endl;
	    std::cout << "-----------------Sizes: " << NodesVectorVariableArray.size() << std::endl;
	    std::cout << "-----------------Sizes: " << NodesMatrixVariableArray.size() << std::endl << std::endl;
	    std::cout << "-----------------Sizes: " << mrRemesh.Refine->ConditionalMeshingNodeVariables.DoubleVariables.size() << std::endl;
	    std::cout << "-----------------Sizes: " << mrRemesh.Refine->ConditionalMeshingNodeVariables.Array1DVariables.size() << std::endl;
	    std::cout << "-----------------Sizes: " << mrRemesh.Refine->ConditionalMeshingNodeVariables.VectorVariables.size() << std::endl;
	    std::cout << "-----------------Sizes: " << mrRemesh.Refine->ConditionalMeshingNodeVariables.MatrixVariables.size() << std::endl;
*/
	    double Area         = 0;
		double ElementArea  = 0;

		//Initialize variables

		//double
		std::fill( NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0);
  
		//Array1D
		for(unsigned int i = 0; i < mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables.size(); i++)  
			NodesArray1DVariableArray[i].clear();

		//Vector
		for(unsigned int i = 0; i < mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables.size(); i++)
			NodesVectorVariableArray[i] = ZeroVector(3);

		//Matrix
		for(unsigned int i = 0; i < mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables.size(); i++)
			NodesMatrixVariableArray[i] = ZeroMatrix(3);
		
		

		//loop over neighbouring elements: obtain weighted results
		for(unsigned int ne=0; ne < neighb_elems.size(); ne++)
		{		    

			ElementArea = neighb_elems[ne].GetGeometry().Area();
			Area += ElementArea;			 

			//double
			for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables.size(); it++)
			{			  
				//elemental value
				neighb_elems[ne].GetValueOnIntegrationPoints(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.DoubleVariables[it]), ElementDoubleVariableArray, CurrentProcessInfo);   

				//add weighted element value
				for(unsigned int j = 0; j < integration_points_number; j++)
				{
					NodesDoubleVariableArray[it] += ElementDoubleVariableArray[j] * ElementArea/integration_points_number;			
					/* //print
					std::cout<<ElementDoubleVariableArray[j]<<std::endl;
					std::cout<<ElementArea/double(integration_points_number)<<std::endl;
					std::cout<<ElementDoubleVariableArray[j] * ElementArea/double(integration_points_number)<<std::endl;
					std::cout<<NodesDoubleVariableArray[it]<<std::endl<<std::endl;
					// */
				}
			}

			//Array1D
			for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables.size(); it++)
			{			  
				//elemental value
				neighb_elems[ne].GetValueOnIntegrationPoints(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.Array1DVariables[it]), ElementArray1DVariableArray, CurrentProcessInfo);
				
				//check size of NodesMatrixVariableArray[it] and resize
				if(NodesArray1DVariableArray[it].size() != ElementArray1DVariableArray[0].size())
					NodesArray1DVariableArray[it].resize( ElementArray1DVariableArray[0].size(), true);
				
				//add weighted element value
				for(unsigned int j = 0; j < integration_points_number; j++)
				{
					MultiSkalarArray1D(&(ElementArray1DVariableArray[j]), ElementArea/integration_points_number);
					NodesArray1DVariableArray[it] += ElementArray1DVariableArray[j];
				}
			}

			//Vector
			for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables.size(); it++)
			{			  
				//elemental value
				neighb_elems[ne].GetValueOnIntegrationPoints(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.VectorVariables[it]), ElementVectorVariableArray, CurrentProcessInfo);
				
				//check size of NodesMatrixVariableArray[it] and resize
				if(NodesVectorVariableArray[it].size() != ElementVectorVariableArray[0].size())
					NodesVectorVariableArray[it].resize( ElementVectorVariableArray[0].size(), true);
				
				//add weighted elemental value
				for(unsigned int j = 0; j < integration_points_number; j++)
				{
					MultiSkalarVector(&(ElementVectorVariableArray[j]), ElementArea/integration_points_number);
					NodesVectorVariableArray[it] += ElementVectorVariableArray[j];
				}
			}
			
			//Matrix
			for(unsigned int it = 0; it < mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables.size(); it++)
			{			  
				//elemental value
				neighb_elems[ne].GetValueOnIntegrationPoints(*(mrRemesh.Refine->ConditionalMeshingGaussVariables.MatrixVariables[it]), ElementMatrixVariableArray, CurrentProcessInfo);
				
				//check size of NodesMatrixVariableArray[it] and resize
				if(NodesMatrixVariableArray[it].size1() != ElementMatrixVariableArray[0].size1() && NodesMatrixVariableArray[it].size2() != ElementMatrixVariableArray[0].size2())
					NodesMatrixVariableArray[it].resize( ElementMatrixVariableArray[0].size1() , ElementMatrixVariableArray[0].size2(), true);
				
				//add weighted elemental value
				for(unsigned int j = 0; j < integration_points_number; j++)
				{
					MultiSkalarMatrix(&(ElementMatrixVariableArray[j]), ElementArea/integration_points_number);
					NodesMatrixVariableArray[it] += ElementMatrixVariableArray[j];
					
					/*/print result
					ConvertMatrixToString(&(NodesMatrixVariableArray[it]), auxString);
					std::cout<<auxString<<std::endl<<std::endl;
					*/
				}
			}
		}

		//divide by element area and write into output array
		if(Area != 0)
		{
			//loop over nodal data double Variables
			for(unsigned int i = 0; i < NodesDoubleVariableArray.size(); i++ )
			{
				NodesDoubleVariableArray[i] /= Area;
				nodalData += (std::to_string(NodesDoubleVariableArray[i]) + " ");
			}
			
			//loop over nodal data array1D Variables
			for(unsigned int i = 0; i < NodesArray1DVariableArray.size(); i++ )
			{
				MultiSkalarArray1D(&(NodesArray1DVariableArray[i]), (1/Area));
				ConvertArray1DToString(&(NodesArray1DVariableArray[i]), auxString);
				nodalData += auxString;
			}
			
			//loop over nodal data vector Variables
			for(unsigned int i = 0; i < NodesVectorVariableArray.size(); i++ )
			{
				MultiSkalarVector(&(NodesVectorVariableArray[i]), (1/Area));
				ConvertVectorToString(&(NodesVectorVariableArray[i]), auxString);
				nodalData += auxString;
			}
			
			//loop over nodal data matrix Variables
			for(unsigned int i = 0; i < NodesMatrixVariableArray.size(); i++ )
			{
				MultiSkalarMatrix(&(NodesMatrixVariableArray[i]), (1/Area));
				ConvertMatrixToString(&(NodesMatrixVariableArray[i]), auxString);
				nodalData += auxString;
			}
		}
		else {
			std::cout << " Error in mapped Gauss results at nodes of interest! Area is equal to 0!" << std::endl;
		}


		//write string into file
		outFile.open(outputFileName, std::fstream::app);
		if(outFile.is_open())
		{
			outFile << nodalData << std::endl;
			outFile.close();
		}
		
		//
		
		KRATOS_CATCH( "" )
	}

	//*******************************************************************************************
    //*******************************************************************************************

	void ConvertArray1DToString(array_1d<double,3>* array1D, std::string& output)
	{
		output = "";
		for (int i = 0; i < (*array1D).size(); i++)
			output += (std::to_string((*array1D)[i]) + " ");
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void ConvertVectorToString(Vector* vec, std::string& output)
	{
		output = "";
		for (int i = 0; i < (*vec).size(); i++)
			output += (std::to_string((*vec)[i]) + " ");
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void ConvertMatrixToString(Matrix* mat, std::string& output)
	{
		output = "";
		for (int i = 0; i < (*mat).size1(); i++)
		{
			for (int j = 0; j < (*mat).size2(); j++)
				output += (std::to_string((*mat)(i,j)) + " ");
		}
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void MultiSkalarArray1D(array_1d<double,3>* array1D, double skalar)
	{
		for (int i = 0; i < (*array1D).size(); i++)
			(*array1D)[i] = (*array1D)[i] * skalar;
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void MultiSkalarVector(Vector* vec, double skalar)
	{
		for (int i = 0; i < (*vec).size(); i++)
			(*vec)[i] = (*vec)[i] * skalar;
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void MultiSkalarMatrix(Matrix* mat, double skalar)
	{
		for (int i = 0; i < (*mat).size1(); i++)
		{
			for (int j = 0; j < (*mat).size2(); j++)
				(*mat)(i,j) = (*mat)(i,j) * skalar;
		}
	}
	
	//*******************************************************************************************
    //*******************************************************************************************

	void CleanPreviousNodes(ModelPart& rModelPart)
    {
		KRATOS_TRY

		//MESH 0 total domain mesh
	    ModelPart::NodesContainerType temporal_nodes;
	    temporal_nodes.reserve(rModelPart.Nodes().size());

	    temporal_nodes.swap(rModelPart.Nodes());

    	for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	    {
			if( (i_node->Is(NEW_ENTITY) && i_node->Is(OLD_ENTITY)))
			{
				if( this->mEchoLevel >=1 )
					std::cout<<" Old 'forced' node releases: ID "<<i_node->Id()<<" ["<<i_node->X()<<", "<<i_node->Y()<<", "<<i_node->Z()<<"]"<<std::endl;

				//write the nodal data to output files before removing the nodes from model part
				this->WriteOutputDataToFiles(rModelPart, i_node);
			}
			else
			{
				(rModelPart.Nodes()).push_back(*(i_node.base()));
			}
	    }

	    rModelPart.Nodes().Sort();

	    KRATOS_CATCH( "" )
    }
    
    //*******************************************************************************************
    //*******************************************************************************************
	
	//from refine_mesh_elements_in_edges_process.hpp

	void GenerateNewNodes(ModelPart& rModelPart,
			std::vector<Node<3>::Pointer>& rListOfNewNodes,
			std::vector<std::array<double, 6> >& rListOfNodesToCreate)
    {
		KRATOS_TRY

		//ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      
		Node<3>::Pointer pNode;

		//center
		double xc = 0;
		double yc = 0;
		double zc = 0;
		
		//time
		double time = mrModelPart.GetProcessInfo().GetValue(TIME);
      
		//radius
		double radius = 0;

		//assign data to dofs
		Node<3>::DofsContainerType& ReferenceDofs = rModelPart.Nodes().front().GetDofs();

		VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();
      
		unsigned int id = ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())) + 1;
		
		//std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())) "<<ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart()))<<std::endl;
		//std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart)) "<<ModelerUtilities::GetMaxNodeId(rModelPart)<<std::endl;

		unsigned int size  = 0;
		unsigned int count = 0;

		//loop over nodes of interest
		for (std::vector<std::array<double, 6> >::iterator i_node = rListOfNodesToCreate.begin() ; i_node != rListOfNodesToCreate.end(); i_node++) //originally ++i_node
		{
			//create a new node
			pNode = boost::make_shared< Node<3> >( id, (*i_node)[0] + ((*i_node)[3])*time, (*i_node)[1] + ((*i_node)[4])*time, (*i_node)[2] + ((*i_node)[5])*time );
			
			//max id has been initialized before loop
			
			//giving model part variables list to the node
			pNode->SetSolutionStepVariablesList(&VariablesList);
			
			//set buffer size
			pNode->SetBufferSize(rModelPart.GetBufferSize());
			
			//generating the dofs; ATM for nodes inside the domain
			for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); i_dof++)
			{
				Node<3>::DofType& rDof = *i_dof;
				Node<3>::DofType::Pointer pNewDof = pNode->pAddDof( rDof );
		    
				(pNewDof)->FreeDof();
			}
		
			//set flags
			pNode->Set(NEW_ENTITY);
			pNode->Set(OLD_ENTITY);
			
			rListOfNewNodes.push_back(pNode);

			if( this->mEchoLevel >=1 )
				std::cout<<" New 'forced' node added: ID "<<pNode->Id()<<" ["<<pNode->X()<<", "<<pNode->Y()<<", "<<pNode->Z()<<"]"<<std::endl;

			id++;  
	    }

		KRATOS_CATCH( "" )
    }
    
    //**************************************************************************
    //**************************************************************************

    void ProjectVariablesToNewNodes(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_new_nodes)			    
    {

		KRATOS_TRY

		//defintions for spatial search
		typedef Node<3>                                  PointType;
		typedef Node<3>::Pointer                  PointPointerType;
		typedef std::vector<PointPointerType>          PointVector;
		typedef PointVector::iterator                PointIterator;
		//typedef std::vector<double>                 DistanceVector;
		typedef std::vector<double>::iterator     DistanceIterator;

		typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
		//defintions for spatial search

		unsigned int  bucket_size = 20;
		KdtreeType    NodesTree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

		//Find out where the new nodes belong to:
		std::vector<double> ShapeFunctionsN;    
		std::vector<VariablesListDataValueContainer> VariablesListVector(list_of_new_nodes.size());

		VariablesList&  variables_list = rModelPart.GetNodalSolutionStepVariablesList();

		//find the center and "radius" of the element
		double  radius = 0;
		Node<3> center(0,0.0,0.0,0.0);

		unsigned int MaximumNumberOfPointsInRadius = list_of_new_nodes.size();
		std::vector<Node<3>::Pointer> PointsInRadius (MaximumNumberOfPointsInRadius);
		std::vector<double>  PointsInRadiusDistances (MaximumNumberOfPointsInRadius);

		//geometry
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();	  
		const unsigned int nds = element_begin->GetGeometry().size();

		std::vector<std::vector<double> > ElementPointCoordinates(nds);
		std::vector<double> PointCoordinates(3);
		std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
		std::fill( ElementPointCoordinates.begin(), ElementPointCoordinates.end(), PointCoordinates );

		//loop over elements
		for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ie++)
		{
			//write element points coordinates into ElementPointCoordinates
			for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
			{	
				PointCoordinates[0] = ie->GetGeometry()[i].X();
				PointCoordinates[1] = ie->GetGeometry()[i].Y();
				PointCoordinates[2] = ie->GetGeometry()[i].Z();

				ElementPointCoordinates[i] = PointCoordinates;
			}

			std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
			MeshDataTransferUtilities DataTransferUtilities;
			//writes the center point of ElementPointCoordintes into PointCoordintes; updated the largest radius
			DataTransferUtilities.CalculateCenterAndSearchRadius( ElementPointCoordinates, PointCoordinates, radius );

			//find all of the new nodes within the radius
			center.X() = PointCoordinates[0];
			center.Y() = PointCoordinates[1];
			center.Z() = PointCoordinates[2];

			double Radius = radius * 1.01;
			int NumberOfPointsInRadius = NodesTree.SearchInRadius (center, Radius, PointsInRadius.begin(), PointsInRadiusDistances.begin(),  MaximumNumberOfPointsInRadius);

			//check if inside and eventually interpolate
			for(std::vector<Node<3>::Pointer>::iterator it_found = PointsInRadius.begin(); it_found != (PointsInRadius.begin() + NumberOfPointsInRadius) ; it_found++)
			{
		
				PointCoordinates[0] = (*it_found)->X();
				PointCoordinates[1] = (*it_found)->Y();
				PointCoordinates[2] = (*it_found)->Z();

				bool is_inside = false;
				//check if new node (PointCoordinates) is inside ElementPointCoordinates and calculates shape functions vector at PointCoordinates
				is_inside = ModelerUtilities::CalculatePosition( ElementPointCoordinates, PointCoordinates, ShapeFunctionsN );

				if(is_inside == true)
				{
					//do the interpolation based on shape functions vector for new node
					double alpha = 1; //1 to interpolate, 0 to leave the original data
					DataTransferUtilities.Interpolate( ie->GetGeometry(), ShapeFunctionsN, variables_list, *(it_found), alpha );		    
				}
			}
		}

		KRATOS_CATCH( "" )
    }
  
    //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewNodeVariables(ModelPart& rModelPart, NodeType::Pointer& pNode)
    {
		KRATOS_TRY
		
		//set variables:
		//Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();

		//set model part
		pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());

		//set nodal_h - whatever this is doing here???
		//pNode->FastGetSolutionStepValue(NODAL_H) = mrRemesh.Refine->CriticalSide; //too small problems     
		//pNode->FastGetSolutionStepValue(NODAL_H) = rGeometry.DomainSize();

		//set normal - for new contact node
		//noalias(pNode->FastGetSolutionStepValue(NORMAL)) = pCondition->GetValue(NORMAL);
	    
		//set original position
		const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
		pNode->X0() = pNode->X() - Displacement[0];
		pNode->Y0() = pNode->Y() - Displacement[1];
		pNode->Z0() = pNode->Z() - Displacement[2];

		/*
		//set contact force
		unsigned int count = 0;
		for( unsigned int i = 0; i<rGeometry.size(); i++ )
	    {
			if( norm_2(rGeometry[i].FastGetSolutionStepValue(CONTACT_FORCE)) == 0 )
				count++;
	    }

		if( count )
			pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();
		*/
		
		//set contact force to zero
		pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();

		KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************
 
    void SetNodesToModelPart(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes)
    {
		KRATOS_TRY

		if(list_of_nodes.size())
		{
			//add new conditions: ( SOLID body model part )
			for(std::vector<NodeType::Pointer>::iterator i_node = list_of_nodes.begin(); i_node!= list_of_nodes.end(); i_node++)
			{
				rModelPart.Nodes().push_back(*(i_node));
			}
			std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())) "<<ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart()))<<std::endl;
			std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart)) "<<ModelerUtilities::GetMaxNodeId(rModelPart)<<std::endl;
		}

		KRATOS_CATCH( "" )
    }
    
    //*******************************************************************************************
    //*******************************************************************************************
 
    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:

    ///@name Private Static Member Variables
    ///@{

    ///@}
    ///@name Private Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{
  
    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{
  
    ///@}
    ///@name Un accessible methods
    ///@{

  
    /// Assignment operator.
    AddNodesOfInterestProcess& operator=(AddNodesOfInterestProcess const& rOther);


    /// this function is a private function


    /// Copy constructor.
    //Process(Process const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AddNodesOfInterestProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AddNodesOfInterestProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_NODES_OF_INTEREST_PROCESS_H_INCLUDED  defined 


