//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"

#include "pfem_base_application_variables.h"

namespace Kratos
{

       KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, NODE_TO_ELEMENT,        0 );
       KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, ELEMENT_TO_NODE,        1 );
       KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, ELEMENT_TO_ELEMENT,     2 );
       KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, MASTER_ELEMENT_TO_NODE, 3 );
       KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, INITIALIZATION,         4 );


  

	//*******************************************************************************************
	//*******************************************************************************************

        void MeshDataTransferUtilities::TransferData(ModelPart& rModelPart,
						     const Element & rReferenceElement,
						     PointPointerVector &list_of_new_centers,
						     std::vector<Geometry<Node<3> > >& list_of_new_vertices,
						     Flags Options,
						     ModelPart::IndexType MeshId)
	  
	{
	  KRATOS_TRY

	  ModelPart::MeshesContainerType Meshes = rModelPart.GetMeshes();
	  
	  
	    if(Options.Is(MeshDataTransferUtilities::NODE_TO_ELEMENT))
	      {
		TransferNodalValuesToElements(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices,MeshId);
	      }
	    else
	      {
		if(Options.Is(MeshDataTransferUtilities::ELEMENT_TO_ELEMENT)){
		  TransferElementalValuesToElements(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices,MeshId);
		}
		else
		  {
		    if(Options.Is(MeshDataTransferUtilities::ELEMENT_TO_NODE)){
		      TransferElementalValuesToNodes(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices,MeshId);
		    }
		  }
	      }
	    
	  KRATOS_CATCH( "" )

	}



	//*******************************************************************************************
	//*******************************************************************************************

        void MeshDataTransferUtilities::InitializeBoundaryData(Condition::Pointer rCondition,
							       const TransferParameters& rTransferVariables)
	{

	  KRATOS_TRY

	  unsigned int StrainSize;
	  const unsigned int dimension = rCondition->GetGeometry().WorkingSpaceDimension();
	  
	  if ( dimension == 2 )
	    StrainSize = 3;
	  else
	    StrainSize = 6;


	  double DoubleVariable = 0;
	  array_1d<double, 3> Array1DVariable;
	  Vector VectorVariable = ZeroVector(StrainSize);
	  Matrix MatrixVariable = identity_matrix<double>( dimension );


	  //double
	  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
	    {
	      rCondition->SetValue(*(rTransferVariables.DoubleVariables[i]),DoubleVariable);
	    }
	  
	  //array_1d
	  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
	    {
	      rCondition->SetValue(*(rTransferVariables.Array1DVariables[i]),Array1DVariable);
	    }
	  
	  //Vector
	  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
	    {
	      rCondition->SetValue(*(rTransferVariables.VectorVariables[i]),VectorVariable);
	    }
	  
	  //Matrix
	  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
	    {
	      rCondition->SetValue(*(rTransferVariables.MatrixVariables[i]),MatrixVariable);
	    }
	  
		    
	  KRATOS_CATCH( "" )

	}



	//*******************************************************************************************
	//*******************************************************************************************

        void MeshDataTransferUtilities::TransferBoundaryData(Condition::Pointer rCurrentCondition, 
							     Condition::Pointer rReferenceCondition, 
							     const TransferParameters& rTransferVariables)
	{
	  KRATOS_TRY

	  //double
	  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
	    {
	      rCurrentCondition->SetValue(*(rTransferVariables.DoubleVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.DoubleVariables[i])) );
	    }
	  
	  //array_1d
	  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
	    {
	      rCurrentCondition->SetValue(*(rTransferVariables.Array1DVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.Array1DVariables[i])) );
	    }
	  
	  //Vector
	  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
	    {
	      rCurrentCondition->SetValue(*(rTransferVariables.VectorVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.VectorVariables[i])) );
	    }
	  
	  //Matrix
	  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
	    {
	      rCurrentCondition->SetValue(*(rTransferVariables.MatrixVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.MatrixVariables[i])) );
	    }
	  
		    
	  KRATOS_CATCH( "" )

	}

	//*******************************************************************************************
	//*******************************************************************************************

        void MeshDataTransferUtilities::TransferBoundaryData(const TransferParameters& rTransferVariables,
							     ModelPart& rModelPart,
							     ModelPart::IndexType MeshId)
	{
	    KRATOS_TRY

	    if(rTransferVariables.Options.Is(MeshDataTransferUtilities::MASTER_ELEMENT_TO_NODE))
	      {
		
		std::cout<<"  TRANSFER MASTER_ELEMENT_TO_NODE "<<std::endl;

		ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		unsigned int StrainSize;
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
		const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();
		
		if ( dimension == 2 )
		  StrainSize = 3;
		else
		  StrainSize = 6;

		double DoubleVariable = 0;
		array_1d<double, 3> Array1DVariable;
		Vector VectorVariable = ZeroVector(StrainSize);
		Matrix MatrixVariable = identity_matrix<double>( dimension );

		
		//initialize to zero all skin master-nodes
		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {
		    //double
		    for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.DoubleVariables[i]),DoubleVariable);
		      }
		    
		    //array_1d
		    for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.Array1DVariables[i]),Array1DVariable);
		      }
		    
		    //Vector
		    for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.VectorVariables[i]),VectorVariable);
		      }

		    //Matrix
		    for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.MatrixVariables[i]),MatrixVariable);
		      }
		    
		  }
		
	      


		//store that information in all body skin if there is a Contact Condition;
		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {

		    if(ic->Is(CONTACT) && ic->Is(ACTIVE)){
		      
		      //std::cout<<" Transfer: Cond: "<<ic->Id()<<" is Active "<<std::endl;

		      Element::ElementType& MasterElement    = ic->GetValue(MASTER_ELEMENTS)[0];
		      Condition::Pointer MasterCondition     = ic->GetValue(MASTER_CONDITION);

		      unsigned int integration_points_number = (MasterElement.pGetGeometry())->IntegrationPointsNumber(MasterElement.GetIntegrationMethod());
		      
		      std::vector<double> DoubleVariableArray (integration_points_number);
		      std::vector<array_1d<double,3> > Array1DVariableArray (integration_points_number);
		      std::vector<Vector> VectorVariableArray (integration_points_number);
		      std::vector<Matrix> MatrixVariableArray (integration_points_number);
		      
		      
		      //double
		      for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
			{
			  MasterElement.GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),DoubleVariableArray,CurrentProcessInfo);

			  //if there is more than one integration point, an average or an interpolation is need		
			  DoubleVariable = DoubleVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      DoubleVariable  += DoubleVariableArray[j];
			    }
			  DoubleVariable *= (1.0/double(integration_points_number));
			  MasterCondition->SetValue(*(rTransferVariables.DoubleVariables[i]),DoubleVariable);
			}
		    
		      //array_1d
		      for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
			{
			  MasterElement.GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),Array1DVariableArray,CurrentProcessInfo);

			  //if there is more than one integration point, an average or an interpolation is need		
			  Array1DVariable = Array1DVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      Array1DVariable += Array1DVariableArray[j];
			    }
			  Array1DVariable *= (1.0/double(integration_points_number));
			  MasterCondition->SetValue(*(rTransferVariables.Array1DVariables[i]),Array1DVariable);
			}
		      
		      //Vector
		      for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
			{
			
			  MasterElement.GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),VectorVariableArray,CurrentProcessInfo);

			  //std::cout<<" TRANSFER EN "<<MasterElement.Id()<<" "<<*rTransferVariables.VectorVariables[i]<<" "<<VectorVariableArray[0]<<std::endl;

			  //if there is more than one integration point, an average or an interpolation is need		
			  VectorVariable = VectorVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      VectorVariable  += VectorVariableArray[j];
			    }
			  VectorVariable *= (1.0/double(integration_points_number));
			  MasterCondition->SetValue(*(rTransferVariables.VectorVariables[i]),VectorVariable);
			}

		      //Matrix
		      for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
			{

			  MasterElement.GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),MatrixVariableArray,CurrentProcessInfo);

			  //std::cout<<" TRANSFER EN "<<MasterElement.Id()<<" "<<*rTransferVariables.MatrixVariables[i]<<" "<<MatrixVariableArray[0]<<std::endl;
			  
			  //if there is more than one integration point, an average or an interpolation is need		
			  MatrixVariable = MatrixVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      MatrixVariable  += MatrixVariableArray[i];
			    }
			  MatrixVariable *= (1.0/double(integration_points_number));
			  MasterCondition->SetValue(*(rTransferVariables.MatrixVariables[i]),MatrixVariable);
			}
		    
     

		      //std::cout<<" MasterCond: "<<MasterCondition->Id()<<" is Active "<<std::endl;		    

		    }
		  }


		    

	      }

	    if(rTransferVariables.Options.Is(MeshDataTransferUtilities::INITIALIZATION))		
	      {
		
		std::cout<<"  TRANSFER INITIALIZE "<<std::endl;

		unsigned int StrainSize;
		
		//I don't use condition because when I create the boundary skin I must assign the 
		//correct geometry to the boundary face, not only some PointsArray
		//It must be done.
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
		const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();
		
		if ( dimension == 2 )
		  StrainSize = 3;
		else
		  StrainSize = 6;


		double DoubleVariable = 0;
		array_1d<double, 3> Array1DVariable;
		Vector VectorVariable = ZeroVector(StrainSize);
		Matrix MatrixVariable = identity_matrix<double>( dimension );

		
		//initialize to zero all skin master-nodes
		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {
		    //double
		    for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.DoubleVariables[i]),DoubleVariable);
		      }
		    
		    //array_1d
		    for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.Array1DVariables[i]),Array1DVariable);
		      }
		    
		    //Vector
		    for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.VectorVariables[i]),VectorVariable);
			//std::cout<<" TRANSFER "<<*rTransferVariables.VectorVariables[i]<<" "<<VectorVariable<<std::endl;
		      }

		    //Matrix
		    for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
		      {
			ic->SetValue(*(rTransferVariables.MatrixVariables[i]),MatrixVariable);
			//std::cout<<" TRANSFER "<<*rTransferVariables.MatrixVariables[i]<<" "<<MatrixVariable<<std::endl;
		      }
		    
		  }

	      }

	  KRATOS_CATCH( "" )

	}

  
      //*******************************************************************************************
      //*******************************************************************************************
	
      void MeshDataTransferUtilities::TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
								    ModelPart& rModelPart,
								    ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
		    
	    //std::cout<<" [ Data Transfer NODE to ELEMENT ] :"<<std::endl;

	    double alpha = 1; //[0,1] //smoothing level of the Jacobian	      

	    Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin(MeshId)->GetGeometry();
	    GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
	    unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );

	    std::vector<double> NodesDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
	    std::vector<Vector> NodesVectorVariableArray (integration_points_number);
	    std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

	    std::vector<double> ElementDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
	    std::vector<Vector> ElementVectorVariableArray (integration_points_number);
	    std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);
		      

	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
		
		Geometry<Node<3> > & rGeometry = (ie)->GetGeometry();
		IntegrationMethod =  rGeometry.GetDefaultIntegrationMethod();
		integration_points_number = rGeometry.IntegrationPointsNumber( IntegrationMethod );

		const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( IntegrationMethod );

		//shape functions
		Vector N;

		//double
		for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
		  {
		    	      
		    //elemental value
		    (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);

		    std::fill(NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0 );

		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			N = row( Ncontainer, j );

			//nodal value
			for( unsigned int k=0 ; k<rGeometry.size(); k++)
			  {
			    NodesDoubleVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i]));
			  }
			
			NodesDoubleVariableArray[j] *= (alpha);
			NodesDoubleVariableArray[j] += (1-alpha) * ElementDoubleVariableArray[j];
						
		      }
		    
	    
		    (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),NodesDoubleVariableArray,CurrentProcessInfo);

		  }

		//array_1d
		for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
		  {
		    	      
		    //elemental value
		    (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);

		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			NodesArray1DVariableArray[j].clear();
		      }

		    
		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			N = row( Ncontainer, j );

			//nodal value
			for( unsigned int k=0 ; k<rGeometry.size(); k++)
			  {
			    NodesArray1DVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i]));
			  }
			
			NodesArray1DVariableArray[j] *= (alpha);
			NodesArray1DVariableArray[j] += (1-alpha) * ElementArray1DVariableArray[j];
						
		      }
		    
	    
		    (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),NodesArray1DVariableArray,CurrentProcessInfo);

		  }
		      
		//Vector
		for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
		  {
		    	      
		    //elemental value
		    (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);


		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			NodesVectorVariableArray[j] = ZeroVector();
		      }

		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			N = row( Ncontainer, j );

			//nodal value
			for( unsigned int k=0 ; k<rGeometry.size(); k++)
			  {
			    NodesVectorVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]));
			  }
			
			NodesVectorVariableArray[j] *= (alpha);
			NodesVectorVariableArray[j] += (1-alpha) * ElementVectorVariableArray[j];
						
		      }
		    
	    
		    (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),NodesVectorVariableArray,CurrentProcessInfo);

		  }


		//Matrix

		for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
		  {
		    	      
		    //elemental value
		    (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);

		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			NodesMatrixVariableArray[j] = ZeroMatrix();
		      }


		    for(unsigned int j=0; j<integration_points_number; j++)
		      {
			N = row( Ncontainer, j );

			//nodal value
			for( unsigned int k=0 ; k<rGeometry.size(); k++)
			  {
			    NodesMatrixVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]));
			  }
			
			NodesMatrixVariableArray[j] *= (alpha);
			NodesMatrixVariableArray[j] += (1-alpha) * ElementMatrixVariableArray[j];
						
		      }
		    
	    
		    (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),NodesMatrixVariableArray,CurrentProcessInfo);

		  }


	      }

	    
	    //std::cout<<" [ Finished NODE to ELEMENT Transfer ]"<<std::endl;
	    
	    KRATOS_CATCH( "" )
	}
      

      //*******************************************************************************************
      //*******************************************************************************************
	
      void MeshDataTransferUtilities::TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
								    const Variable<double>& rCriticalVariable,
								    const double& rCriticalValue,
								    ModelPart& rModelPart,
								    ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
		    

	    //std::cout<<" [ Data Transfer NODE to ELEMENT ] : based on critical values of "<<rCriticalVariable<<std::endl;

	    double alpha = 0.5; //[0,1] //smoothing level of the Jacobian	      

	    Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin(MeshId)->GetGeometry();
	    GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
	    unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );


            std::vector<double> ComputedValues(integration_points_number);
	    double computed_value=0;
	    double critical_value=0;

	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
	     
	    	(ie)->GetValueOnIntegrationPoints(rCriticalVariable,ComputedValues,CurrentProcessInfo);
		computed_value = ComputedValues[0];
		for(unsigned int j=1; j<integration_points_number; j++)
		  {
		    computed_value += ComputedValues[j];
		  }
		
		computed_value *= (ie)->GetGeometry().Area()/double(integration_points_number);

	    	critical_value = rCriticalValue;

	    	if( computed_value > critical_value )
	    	  {
	    	    for(unsigned int i = 0; i<ie->GetGeometry().size(); i++)
	    	      {
	    		(ie)->GetGeometry()[i].Set(TO_REFINE);
	    	      }
	    	  }
	      }

	    std::vector<double> NodesDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
	    std::vector<Vector> NodesVectorVariableArray (integration_points_number);
	    std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

	    std::vector<double> ElementDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
	    std::vector<Vector> ElementVectorVariableArray (integration_points_number);
	    std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);
		
	    int counter = 0;

	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
		
		Geometry<Node<3> > & rGeometry = (ie)->GetGeometry();
		IntegrationMethod =  rGeometry.GetDefaultIntegrationMethod();
		integration_points_number = rGeometry.IntegrationPointsNumber( IntegrationMethod );
		const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( IntegrationMethod );

		//shape functions
		Vector N;

		bool apply_smoothing = false;
		for( unsigned int k=0 ; k<rGeometry.size(); k++)
		  {
		    if(rGeometry[k].Is(TO_REFINE))
		      apply_smoothing = true;
		  }
		
		if( apply_smoothing ){

		  counter++;
		  //double
		  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
		    {
		    	      
		      //elemental value
		      (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);

		      std::fill(NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0 );

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  N = row( Ncontainer, j );

			  //nodal value
			  for( unsigned int k=0 ; k<rGeometry.size(); k++)
			    {
			      NodesDoubleVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i]));
			    }
			
			  NodesDoubleVariableArray[j] *= (alpha);
			  NodesDoubleVariableArray[j] += (1-alpha) * ElementDoubleVariableArray[j];
						
			}
		    
	    
		      (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),NodesDoubleVariableArray,CurrentProcessInfo);

		    }

		  //array_1d
		  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
		    {
		    	      
		      //elemental value
		      (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);
		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  NodesArray1DVariableArray[j].clear();
			}

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  N = row( Ncontainer, j );

			  //nodal value
			  for( unsigned int k=0 ; k<rGeometry.size(); k++)
			    {
			      NodesArray1DVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i]));
			    }
			
			  NodesArray1DVariableArray[j] *= (alpha);
			  NodesArray1DVariableArray[j] += (1-alpha) * ElementArray1DVariableArray[j];
						
			}
		    
	    
		      (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),NodesArray1DVariableArray,CurrentProcessInfo);

		    }
		      
		  //Vector
		  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
		    {
		    	      
		      //elemental value
		      (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  NodesVectorVariableArray[j] = ZeroVector();
			}

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  N = row( Ncontainer, j );

			  //nodal value
			  for( unsigned int k=0 ; k<rGeometry.size(); k++)
			    {
			      NodesVectorVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]));
			    }
			
			  NodesVectorVariableArray[j] *= (alpha);
			  NodesVectorVariableArray[j] += (1-alpha) * ElementVectorVariableArray[j];
						
			}
		    
	    
		      (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),NodesVectorVariableArray,CurrentProcessInfo);

		    }


		  //Matrix

		  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
		    {
		    	      
		      //elemental value
		      (ie)->GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  NodesMatrixVariableArray[j] = ZeroMatrix();
			}

		      for(unsigned int j=0; j<integration_points_number; j++)
			{
			  N = row( Ncontainer, j );

			  //nodal value
			  for( unsigned int k=0 ; k<rGeometry.size(); k++)
			    {
			      NodesMatrixVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]));
			    }
			
			  NodesMatrixVariableArray[j] *= (alpha);
			  NodesMatrixVariableArray[j] += (1-alpha) * ElementMatrixVariableArray[j];
						
			}
		    
	    
		      (ie)->SetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),NodesMatrixVariableArray,CurrentProcessInfo);

		    }

		}
	      }

	    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++)
	      {
	    	in->Reset(TO_REFINE);
	      }

	    //std::cout<<" [ Finished NODE to ELEMENT Transfer ] : ( Performed "<<counter<<" transfers of "<<rModelPart.NumberOfElements(MeshId)<<" possible )"<<std::endl;

	    
	    KRATOS_CATCH( "" )
	}
      

	//*******************************************************************************************
	//*******************************************************************************************
        void MeshDataTransferUtilities::TransferElementalValuesToNodes( const TransferParameters& rTransferVariables,
									ModelPart& rModelPart,
									ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
	     
	    //std::cout<<" [ Data Transfer ELEMENT to NODE ] :"<<std::endl;

	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
	    NodesContainerType& rNodes      = rModelPart.Nodes(MeshId);

	    Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin(MeshId)->GetGeometry();
	    GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
	    unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );

	    std::vector<double> NodesDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
	    std::vector<Vector> NodesVectorVariableArray (integration_points_number);
	    std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

	    std::vector<double> ElementDoubleVariableArray (integration_points_number);
	    std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
	    std::vector<Vector> ElementVectorVariableArray (integration_points_number);
	    std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);
      
	    unsigned int buffer_size = rModelPart.GetBufferSize();
	    VariablesList& variables_list = rModelPart.GetNodalSolutionStepVariablesList();
	    

	    for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	      {

		if(in->IsNot(RIGID)){
		  
		  //fill variables that are non assigned vectors
		  this->FillVectorData(variables_list, *(in.base()));	  

	    	  WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);

	    	  double Area         = 0;
	    	  double ElementArea  = 0;
		  std::fill( NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0);
		  
	    	  for(unsigned int ne=0; ne < neighb_elems.size(); ne++)
	    	    {		    

		      Geometry<Node<3> > & rGeometry   = neighb_elems[ne].GetGeometry();
		      ElementArea = rGeometry.Area();
		      Area += ElementArea;			 
		      
		      //double
		      for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
			{			  
			  //elemental value
			  neighb_elems[ne].GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);   
			  NodesDoubleVariableArray[i] = ElementDoubleVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      NodesDoubleVariableArray[i] += ElementDoubleVariableArray[j];
			    }
			  
			  NodesDoubleVariableArray[i] *= ElementArea/double(integration_points_number);
			}


		      //Array1D
		      for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
			{			  
			  //elemental value
			  neighb_elems[ne].GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);
			  NodesArray1DVariableArray[i] = ElementArray1DVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      NodesArray1DVariableArray[i] += ElementArray1DVariableArray[j];
			    }
			  
			  NodesArray1DVariableArray[i] *= ElementArea/double(integration_points_number);
			}

		      //Vector
		      for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
			{			  
			  //elemental value
			  neighb_elems[ne].GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);
			  NodesVectorVariableArray[i] = ElementVectorVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      NodesVectorVariableArray[i] += ElementVectorVariableArray[j];
			    }
			  
			  NodesVectorVariableArray[i] *= ElementArea/double(integration_points_number);
			}

		      //Matrix
		      for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
			{			  
			  //elemental value
			  neighb_elems[ne].GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);
			  NodesMatrixVariableArray[i] = ElementMatrixVariableArray[0];
			  for(unsigned int j=1; j<integration_points_number; j++)
			    {
			      NodesMatrixVariableArray[i] += ElementMatrixVariableArray[j];
			    }
			  
			  NodesMatrixVariableArray[i] *= ElementArea/double(integration_points_number);
			}

	    	    }
		
	    	  if(Area!=0){
		    //double
		    for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); i++)
		      {			  
			NodesDoubleVariableArray[i] /= Area;

			if( (in)->SolutionStepsDataHas(*(rTransferVariables.DoubleVariables[i]))){
			  (in)->FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i])) = NodesDoubleVariableArray[i];
			}
			else{
			  std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : variable "<<*(rTransferVariables.DoubleVariables[i])<<" was not defined "<<std::endl;
			}
		      }
		    //Array1D
		    for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); i++)
		      {			  
			NodesArray1DVariableArray[i] /= Area;

			if( (in)->SolutionStepsDataHas(*(rTransferVariables.Array1DVariables[i]))){
			  (in)->FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i])) = NodesArray1DVariableArray[i];
			}
			else{
			  std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : variable "<<*(rTransferVariables.Array1DVariables[i])<<" was not defined "<<std::endl;
			}
		      }
		    //Vector
		    for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); i++)
		      {			  
			NodesVectorVariableArray[i] /= Area;

			if( (in)->SolutionStepsDataHas(*(rTransferVariables.VectorVariables[i]))){
			  (in)->FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i])) = NodesVectorVariableArray[i];
			  //fill buffer if empty
			  for(unsigned int step = 1; step<buffer_size; step++)
			    {
			      if((in)->FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step).size() == 0){
				(in)->FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step) = NodesVectorVariableArray[i];
				(in)->FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step).clear();
			      }
			      
			    }
			}
			else{
			  std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : variable "<<*(rTransferVariables.VectorVariables[i])<<" was not defined "<<std::endl;
			}

		      }
		    //Matrix
		    for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); i++)
		      {			  
			NodesMatrixVariableArray[i] /= Area;

			if( (in)->SolutionStepsDataHas(*(rTransferVariables.MatrixVariables[i]))){
			  (in)->FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i])) = NodesMatrixVariableArray[i];
			  //fill buffer if empty
			  for(unsigned int step = 1; step<buffer_size; step++)
			    {
			      if((in)->FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).size1() == 0 &&
				 (in)->FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).size2() == 0 ){
				(in)->FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step) = NodesMatrixVariableArray[i];
				(in)->FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).clear();
			      }
			    }
			}
			else{
			  std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : variable "<<*(rTransferVariables.MatrixVariables[i])<<" was not defined "<<std::endl;
			}

		      }		    
		  }
	    	  else{
	    	    std::cout<<" ERROR TR: Something Wrong in node ["<<(in)->Id()<<"] : Area = 0 (neighbours: "<<neighb_elems.size()<<") "<<std::endl;
		  }

		    
	    	}

	      }
    

	    //std::cout<<" [ Finished ELEMENT to NODE Transfer ]"<<std::endl;
		   
	    KRATOS_CATCH( "" )
        }



	//*******************************************************************************************
	//*******************************************************************************************
	void MeshDataTransferUtilities::TransferNodalValuesToElements(ModelPart& rModelPart,
								      const Element & rReferenceElement,
								      PointPointerVector &list_of_new_centers,
								      std::vector<Geometry<Node<3> > >& list_of_new_vertices,
								      ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
		    
	    std::cout<<" [ Data Transfer NODE to ELEMENT: NOT IMPLEMENTED YET ]"<<std::endl;
	    std::cout<<" [ Finished NODE to ELEMENT Transfer: NOT IMPLEMENTED YET ]"<<std::endl;
	    
	    KRATOS_CATCH( "" )

	}


	//*******************************************************************************************
	//*******************************************************************************************
	void MeshDataTransferUtilities::TransferElementalValuesToNodes(ModelPart& rModelPart,
								       const Element & rReferenceElement,
								       PointPointerVector &list_of_new_centers,
								       std::vector<Geometry<Node<3> > >& list_of_new_vertices,
								       ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY

	    std::cout<<" [ Data Transfer ELEMENT to NODE: NOT IMPLEMENTED YET ]"<<std::endl;
	    std::cout<<" [ Finished ELEMENT to NODE Transfer: NOT IMPLEMENTED YET ]"<<std::endl;

    	    KRATOS_CATCH( "" )

	}


	//*******************************************************************************************
	//*******************************************************************************************
	void MeshDataTransferUtilities::TransferElementalValuesToElements(ModelPart& rModelPart,
									  const Element & rReferenceElement,
									  PointPointerVector &list_of_new_centers,
									  std::vector<Geometry<Node<3> > >& list_of_new_vertices,				       
									  ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY

	    //std::cout<<" [ Data Transfer ELEMENT to ELEMENT ]"<<std::endl;

	    //definitions for spatial search
  	    typedef Node<3>                                  PointType;
	    typedef Node<3>::Pointer                  PointPointerType;
	    typedef std::vector<PointPointerType>   PointPointerVector;
	    typedef std::vector<PointType>             PointTypeVector;
	    typedef PointPointerVector::iterator         PointIterator;
	    typedef std::vector<double>                 DistanceVector;
	    typedef std::vector<double>::iterator     DistanceIterator;
	    typedef Bucket<3, PointType, PointPointerVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
	    typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
	    //definitions for spatial search


	    ElementsContainerType& rPreElements = rModelPart.Elements(MeshId);
		
	    //creating an auxiliary list for the pre integration points
	    PointPointerVector list_of_pre_centers;
	    	    

	    //find the center and "radius" of the element
	    double xc,  yc, zc=0, radius;

	    int moved_transfers = 0;

	    for(ElementsContainerType::iterator ie = rPreElements.begin(); ie != rPreElements.end(); ie++)
	    {
		PointsArrayType& vertices=ie->GetGeometry().Points();

		CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(),
						vertices[1].X(), vertices[1].Y(),
						vertices[2].X(), vertices[2].Y(),
						xc,yc,radius);
		int id= ie->Id();
		PointPointerType p_center = PointPointerType( new PointType(id,xc,yc,zc) );

		// if ((*ie.base())->GetOptions().Is(Element::THERMAL))
		//   std::cout<<" is thermal "<<std::endl;

		// if ((*ie.base())->GetOptions().Is(Element::MECHANICAL))
		//   std::cout<<" is mechanical "<<std::endl;

		list_of_pre_centers.push_back( p_center );

		//std::cout<<" id pre elems "<<list_of_pre_centers.back()->Id()<<std::endl;
	    }

	    //std::cout<<" list of pre centers "<<list_of_pre_centers.size()<<" list of new centers "<<list_of_new_centers.size()<<std::endl;

	    //creating an auxiliary list for the pre integration points
	    unsigned int   bucket_size = 20;
	    KdtreeType     nodes_tree(list_of_pre_centers.begin(),list_of_pre_centers.end(),bucket_size);	

	    PointType work_point(0,0.0,0.0,0.0);
	    double  ResultDistance;

	    //make a loop on temporal elements

	    ElementsContainerType temporal_elements;
	    temporal_elements.reserve(rModelPart.Elements(MeshId).size());
	    
	    temporal_elements.swap(rModelPart.Elements(MeshId));

	    int count=0;
	    for(PointPointerVector::iterator i_center = list_of_new_centers.begin() ; i_center != list_of_new_centers.end() ; i_center++)
	    {
	      count++;
	      //std::cout<<" NUM "<<count<<std::endl;
	      //find the nearest integration point to the new integration point (work_point)
	      PointType& work_point= (**i_center);


	      //std::cout<<" point post "<<work_point.Id()<<" "<<xc<<" "<<yc<<std::endl;
	      PointPointerType result_point = nodes_tree.SearchNearestPoint(work_point,ResultDistance);

	      ElementsContainerType::iterator pe = temporal_elements.find( result_point->Id() );
	      

	      Element::Pointer PreviousElement = Element::Pointer( *pe.base() );

	      PointsArrayType& vertices=PreviousElement->GetGeometry().Points();



	      if(ResultDistance>1e-20 || ResultDistance == 0){
		// std::cout<<"Id: "<<(*i_center)->Id()<<" Connectivities : ["<<list_of_new_vertices[(*i_center)->Id()-1][0].Coordinates()<<", "<<list_of_new_vertices[(*i_center)->Id()-1][1].Coordinates()<<", "<<list_of_new_vertices[(*i_center)->Id()-1][2].Coordinates()<<"] "<<std::endl;

		
		// std::cout<<"Pre Id: "<<PreviousElement->Id()<<" Connectivities : ["<<vertices[0].Coordinates()<<", "<<vertices[1].Coordinates()<<", "<<vertices[2].Coordinates()<<"] "<<std::endl;
		  
		// std::cout<<"New Center :"<<work_point<<", Old Center :"<<*result_point<<"[Distance: "<<ResultDistance<<"]"<<std::endl;

		 moved_transfers ++;
	      }
	      
	      // std::cout<<" Result distance "<<ResultDistance<<" point Id "<<result_point->Id()<<" center Id "<<(*i_center)->Id()<<std::endl;


	      //get transfer variables

	      
	      //Inside of the loop create the element, set the variables, and push_back to the model part.
	      Element::Pointer new_element = PreviousElement->Clone((*i_center)->Id(), list_of_new_vertices[(*i_center)->Id()-1]);	
	      
	      //set transfer variables
	      new_element->SetValue(DOMAIN_LABEL,vertices[0].GetValue(DOMAIN_LABEL)); //DOMAIN_LABEL set as a variable
	      new_element->AssignFlags(*PreviousElement);	

	      //check
	      //new_element->PrintInfo(std::cout);

	      //setting new elements
	      (rModelPart.Elements(MeshId)).push_back(new_element);

		
		//ELEMENT TRANSFER CHECK//
		// std::cout<<" Write Stress [element: "<<new_element->Id()<<"]: ";
		// std::cout.flush();
		// std::vector<Matrix >  StressMatrix;
		// Matrix Stress;
		// Stress.clear();
		// StressMatrix.push_back(Stress);
		// ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		// new_element->GetValueOnIntegrationPoints(CAUCHY_STRESS_TENSOR,StressMatrix,CurrentProcessInfo);
		// std::cout<<StressMatrix[0]<<std::endl;
		//ELEMENT TRANSFER CHECK//
	    }


	    //std::cout<<" [ MOVED TRANSFERS: "<<moved_transfers<<" ]"<<std::endl;
	    //std::cout<<" [ Finished ELEMENT to ELEMENT Transfer ]"<<std::endl;


	    KRATOS_CATCH( "" )
	}
	

	//*******************************************************************************************
	//*******************************************************************************************

	// inline void MeshDataTransferUtilities::CalculateCenterAndSearchRadius(const double x0, const double y0,
	// 							     const double x1, const double y1,
	// 							     const double x2, const double y2,
	// 							     double& xc, double& yc, double& R)

	// {
	//     xc = 0.3333333333333333333*(x0+x1+x2);
	//     yc = 0.3333333333333333333*(y0+y1+y2);

	//     double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
	//     double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
	//     double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

	//     R = R1;
	//     if(R2 > R) R = R2;
	//     if(R3 > R) R = R3;

	//     R = sqrt(R);
	// }


	//*******************************************************************************************
	//*******************************************************************************************

	// inline double MeshDataTransferUtilities::CalculateVol(const double x0, const double y0,
	// 					     const double x1, const double y1,
	// 					     const double x2, const double y2)
	// {
	//     return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
	// }


	//*******************************************************************************************
	//*******************************************************************************************

	// inline bool MeshDataTransferUtilities::CalculatePosition(const double x0, const double y0,
	// 						const double x1, const double y1,
	// 						const double x2, const double y2,
	// 						const double xc, const double yc,
	// 						array_1d<double,3>& N)
	// {
	//     double area = CalculateVol(x0,y0,x1,y1,x2,y2);

	//     if(area < 1e-20)
	//     {
	// 	KRATOS_THROW_ERROR( std::logic_error,"element with zero area found", "" )
	//     }

	//     N[0] = CalculateVol(x1,y1,x2,y2,xc,yc)  / area;
	//     N[1] = CalculateVol(x2,y2,x0,y0,xc,yc)  / area;
	//     N[2] = CalculateVol(x0,y0,x1,y1,xc,yc)  / area;

	//     double tol = 1e-4;
	//     double upper_limit = 1.0+tol;
	//     double lower_limit = -tol;

	//     if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
	// 	return true;

	//     return false;
	// }


      //*******************************************************************************************
      //*******************************************************************************************
  
      VariablesListDataValueContainer MeshDataTransferUtilities::InterpolateVariables( Geometry<Node<3> > &geom,
										       const array_1d<double,3>& N,
										       VariablesList& rVariablesList,
										       Node<3>::Pointer pnode,
										       double& alpha)
      {

	  KRATOS_TRY

	  //Copy Variables List
	  std::cout<<" node["<<pnode->Id()<<"] Data "<<(pnode)->SolutionStepData()<<std::endl;

	  VariablesListDataValueContainer PreviousVariablesListData = (pnode)->SolutionStepData();

	  std::cout<<" CopiedData "<<PreviousVariablesListData<<std::endl;

      
	  Interpolate( geom, N, rVariablesList, pnode, alpha);

	  VariablesListDataValueContainer CurrentVariablesListData = (pnode)->SolutionStepData();

	  (pnode)->SolutionStepData() = PreviousVariablesListData;
      
	  // std::cout<<" PreviousVariables "<<PreviousVariablesListData<<std::endl;
	  // std::cout<<" CurrentVariables "<<CurrentVariablesListData<<std::endl;

	  return CurrentVariablesListData;

	  KRATOS_CATCH( "" )

      }


      //*******************************************************************************************
      //*******************************************************************************************
  
      void MeshDataTransferUtilities::FillVectorData(VariablesList& rVariablesList,
						     Node<3>::Pointer pnode)
      {

	  KRATOS_TRY

	  unsigned int buffer_size = pnode->GetBufferSize();
     
	  for(VariablesList::const_iterator i_variable =  rVariablesList.begin();  i_variable != rVariablesList.end() ; i_variable++)
	    {
	      //std::cout<<" name "<<i_variable->Name()<<std::endl;
	      //std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
	      std::string variable_name = i_variable->Name();
	      if(KratosComponents<Variable<Vector > >::Has(variable_name))
		{
		  //std::cout<<"Vector"<<std::endl;
		  Variable<Vector> variable = KratosComponents<Variable<Vector > >::Get(variable_name);
		  for(unsigned int step = 0; step<buffer_size; step++)
		    {
		      //getting the data of the solution step
		      Vector& node_data = pnode->FastGetSolutionStepValue(variable, step);
		      
		      if( variable == CONSTRAINT_LABELS ){
			std::cout<<" node ["<<pnode->Id()<<"]"<<std::endl;
			std::cout<<" variable "<<variable<<std::endl;
			std::cout<<" node data "<<node_data<<std::endl;
		      }

		      if( node_data.size() == 0 ){
			node_data = ZeroVector(1);	       		      
		      }
		      
		    }
		}
	      
	    }

	    KRATOS_CATCH( "" )

      }


      //*******************************************************************************************
      //*******************************************************************************************
  
      void MeshDataTransferUtilities::Interpolate( Geometry<Node<3> > &geom,
						   const array_1d<double,3>& N,
						   VariablesList& rVariablesList,
						   Node<3>::Pointer pnode,
						   double& alpha)
      {

	  KRATOS_TRY

	  unsigned int buffer_size = pnode->GetBufferSize();
	     
	  if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0){
	    KRATOS_THROW_ERROR( std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR", "" )
	      }
     
	  for(VariablesList::const_iterator i_variable =  rVariablesList.begin();  i_variable != rVariablesList.end() ; i_variable++)
	    {
	      //std::cout<<" name "<<i_variable->Name()<<std::endl;
	      //std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
	      std::string variable_name = i_variable->Name();
	      if(KratosComponents<Variable<double> >::Has(variable_name))
		{
		  //std::cout<<"double"<<std::endl;
		  Variable<double> variable = KratosComponents<Variable<double> >::Get(variable_name);
		  for(unsigned int step = 0; step<buffer_size; step++)
		    {
		      //getting the data of the solution step
		      double& node_data = pnode->FastGetSolutionStepValue(variable, step);
		  
		      double& node0_data = geom[0].FastGetSolutionStepValue(variable, step);
		      double& node1_data = geom[1].FastGetSolutionStepValue(variable, step);
		      double& node2_data = geom[2].FastGetSolutionStepValue(variable, step);
		  
		      if(alpha != 1 )
			node_data = (alpha) * (N[0]*node0_data + N[1]*node1_data + N[2]*node2_data) + (1-alpha) * node_data;
		      else
			node_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
		    }
		}
	      else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
		{
		  //std::cout<<"array1d"<<std::endl;
		  Variable<array_1d<double, 3> > variable = KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name);
		  for(unsigned int step = 0; step<buffer_size; step++)
		    {
		      //getting the data of the solution step
		      array_1d<double, 3>& node_data = pnode->FastGetSolutionStepValue(variable, step);
		  
		      array_1d<double, 3>& node0_data = geom[0].FastGetSolutionStepValue(variable, step);
		      array_1d<double, 3>& node1_data = geom[1].FastGetSolutionStepValue(variable, step);
		      array_1d<double, 3>& node2_data = geom[2].FastGetSolutionStepValue(variable, step);
		  
		      if(alpha != 1 )
			node_data = (alpha) * (N[0]*node0_data + N[1]*node1_data + N[2]*node2_data) + (1-alpha) * node_data;
		      else
			node_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
		    }

		}
	      else if(KratosComponents<Variable<int > >::Has(variable_name))
		{
		  //std::cout<<"int"<<std::endl;
		  //NO INTERPOLATION
		}
	      else if(KratosComponents<Variable<bool > >::Has(variable_name))
		{
		  //std::cout<<"bool"<<std::endl;
		  //NO INTERPOLATION
		}
	      else if(KratosComponents<Variable<Matrix > >::Has(variable_name))
		{
		  //std::cout<<"Matrix"<<std::endl;
		  Variable<Matrix> variable = KratosComponents<Variable<Matrix > >::Get(variable_name);
		  for(unsigned int step = 0; step<buffer_size; step++)
		    {
		      //getting the data of the solution step
		      Matrix& node_data = pnode->FastGetSolutionStepValue(variable, step);
		  
		      Matrix& node0_data = geom[0].FastGetSolutionStepValue(variable, step);
		      Matrix& node1_data = geom[1].FastGetSolutionStepValue(variable, step);
		      Matrix& node2_data = geom[2].FastGetSolutionStepValue(variable, step);
		  
		      if( node_data.size1() > 0 && node_data.size2() ){
			if( node_data.size1() == node0_data.size1() && node_data.size2() == node0_data.size2() &&
			    node_data.size1() == node1_data.size1() && node_data.size2() == node1_data.size2() &&
			    node_data.size1() == node2_data.size1() && node_data.size2() == node2_data.size2()) {
		      
			  if(alpha != 1 )
			    node_data = (alpha) * (N[0]*node0_data + N[1]*node1_data + N[2]*node2_data) + (1-alpha) * node_data;
			  else
			    node_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
			  
			}
		      }
		    }

		}
	      else if(KratosComponents<Variable<Vector > >::Has(variable_name))
		{
		  //std::cout<<"Vector"<<std::endl;
		  Variable<Vector> variable = KratosComponents<Variable<Vector > >::Get(variable_name);
		  for(unsigned int step = 0; step<buffer_size; step++)
		    {
		      //getting the data of the solution step
		      Vector& node_data = pnode->FastGetSolutionStepValue(variable, step);
		  
		      Vector& node0_data = geom[0].FastGetSolutionStepValue(variable, step);
		      Vector& node1_data = geom[1].FastGetSolutionStepValue(variable, step);
		      Vector& node2_data = geom[2].FastGetSolutionStepValue(variable, step);
		  
		      if( variable != CONSTRAINT_LABELS && variable != LOAD_LABELS && variable != MARKER_LABELS )
			{
			  std::cout<<" node ["<<pnode->Id()<<"]"<<std::endl;
			  std::cout<<" variable "<<variable<<std::endl;
			  std::cout<<" node data "<<node_data<<std::endl;
			  if( node_data.size() > 0 ){
			    if( node_data.size() == node0_data.size() &&
				node_data.size() == node1_data.size() &&
				node_data.size() == node2_data.size()) {
			      
			      if(alpha != 1 )
				node_data = (alpha) * (N[0]*node0_data + N[1]*node1_data + N[2]*node2_data) + (1-alpha) * node_data;
			      else
				node_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
			      
			    }
			  }
			}

		    }
		}

	    }

	    KRATOS_CATCH( "" )

	}


	//*******************************************************************************************
	//*******************************************************************************************

	void MeshDataTransferUtilities::InterpolateData( Geometry<Node<3> >& geom,
							 const array_1d<double,3>& N,
							 unsigned int step_data_size,
							 Node<3>::Pointer pnode,
							 double& alpha)
	{

   	    KRATOS_TRY

	    unsigned int buffer_size = pnode->GetBufferSize();

	    if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0)
		KRATOS_THROW_ERROR( std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR", "" )

	    //alpha [0,1] //smoothing level of the interpolation

	    for(unsigned int step = 0; step<buffer_size; step++)
	    {
		//getting the data of the solution step
		double* step_data = (pnode)->SolutionStepData().Data(step);

		double* node0_data = geom[0].SolutionStepData().Data(step);
		double* node1_data = geom[1].SolutionStepData().Data(step);
		double* node2_data = geom[2].SolutionStepData().Data(step);

		//copying this data in the position of the vector we are interested in
		for(unsigned int j= 0; j<step_data_size; j++)
		{
		  step_data[j] = (alpha) * (N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j]) + (1-alpha) * step_data[j];
		}
	    }

	    KRATOS_CATCH( "" )
	}

	//*******************************************************************************************
	//*******************************************************************************************
  
        VariablesListDataValueContainer MeshDataTransferUtilities::InterpolateVariablesData( Geometry<Node<3> >& geom,
											     const array_1d<double,3>& N,
											     unsigned int step_data_size,
											     Node<3>::Pointer pnode,
											     double& alpha)
	  
	{
	    KRATOS_TRY

	    //Copy Variables List
	    VariablesListDataValueContainer PreviousVariablesListData = (pnode)->SolutionStepData();
      
	    InterpolateData( geom, N, step_data_size, pnode, alpha);
	    
	    VariablesListDataValueContainer CurrentVariablesListData = (pnode)->SolutionStepData();
	  
	    (pnode)->SolutionStepData() = PreviousVariablesListData;
      
	    //std::cout<<" PreviousVariables "<<PreviousVariablesListData<<std::endl;
	    //std::cout<<" CurrentVariables "<<CurrentVariablesListData<<std::endl;
	    
	    return CurrentVariablesListData;
	    
	    KRATOS_CATCH( "" )
	    
	}


}  // namespace Kratos.

