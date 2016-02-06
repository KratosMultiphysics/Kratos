//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"

#include "pfem_solid_mechanics_application_variables.h"

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
	    
	}

	//*******************************************************************************************
	//*******************************************************************************************

        void MeshDataTransferUtilities::TransferBoundaryData(ModelPart& rModelPart,
						    Flags Options,
						    ModelPart::IndexType MeshId)
	{

	    if(Options.Is(MeshDataTransferUtilities::MASTER_ELEMENT_TO_NODE))
	      {
		
		ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
		unsigned int StrainSize;
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();
		const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();
		
		if ( dimension == 2 )
		  StrainSize = 3;
		else
		  StrainSize = 6;

		Vector StressVector=ZeroVector(StrainSize);
		Matrix DeformationGradient=identity_matrix<double>( dimension );

		
		//initialize to zero all skin master-nodes
		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {
		    ic->SetValue(CAUCHY_STRESS_VECTOR,StressVector);
		    ic->SetValue(DEFORMATION_GRADIENT,DeformationGradient);
		  }
		
	      


		//store that information in all body skin if there is a Contact Condition;
		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {

		    StressVector=ZeroVector(StrainSize);
		    DeformationGradient=identity_matrix<double>( dimension );

		    if(ic->Is(CONTACT) && ic->Is(ACTIVE)){
		      
		      //std::cout<<" Transfer: Cond: "<<ic->Id()<<" is Active "<<std::endl;

		      Element::ElementType& MasterElement    = ic->GetValue(MASTER_ELEMENTS)[0];
		      Condition::Pointer MasterCondition     = ic->GetValue(MASTER_CONDITION);

		      unsigned int integration_points_number = (MasterElement.pGetGeometry())->IntegrationPointsNumber(MasterElement.GetIntegrationMethod());
		      //step Deformation Gradient F and total CAUCHY_STRESS
		      std::vector<Vector> StressArray (integration_points_number);
		      //StressArray.push_back(StressVector);
	      
		      MasterElement.GetValueOnIntegrationPoints(CAUCHY_STRESS_VECTOR,StressArray,CurrentProcessInfo);
		      std::vector<Matrix> DeformationGradientArray (integration_points_number);
		      //DeformationGradientArray.push_back(DeformationGradient);
	      
		      MasterElement.GetValueOnIntegrationPoints(DEFORMATION_GRADIENT,DeformationGradientArray,CurrentProcessInfo); 
		      
		      //if there is more than one integration point, an average or an interpolation is need			    
		      StressVector        = StressArray[0];
		      DeformationGradient = DeformationGradientArray[0];
		      for(unsigned int i=1; i<integration_points_number; i++)
			{
			  StressVector+=StressArray[i];
			  DeformationGradient+=DeformationGradientArray[i];
			}
		      StressVector        *= (1.0/double(integration_points_number));
		      DeformationGradient *= (1.0/double(integration_points_number));

		      //std::cout<<" MasterCond: "<<MasterCondition->Id()<<" is Active "<<std::endl;
		      MasterCondition->SetValue(CAUCHY_STRESS_VECTOR,StressArray[0]);
		      MasterCondition->SetValue(DEFORMATION_GRADIENT,DeformationGradientArray[0]);
		    

		    }
		  }


		    

	      }

	    if(Options.Is(MeshDataTransferUtilities::INITIALIZATION))		
	      {
		
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

		Vector StressVector=ZeroVector(StrainSize);
		Matrix DeformationGradient=identity_matrix<double>( dimension );

		for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
		  {
		    ic->SetValue(CAUCHY_STRESS_VECTOR,StressVector);
		    ic->SetValue(DEFORMATION_GRADIENT,DeformationGradient);
		  }

	      }

	}

  
      //*******************************************************************************************
      //*******************************************************************************************
	
      void MeshDataTransferUtilities::TransferNodalValuesToElements(const Variable<double>& rVariable,
								    ModelPart& rModelPart,
								    ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
		    
            if( GetEchoLevel() > 0 )
	      std::cout<<" [ Data Transfer NODE to ELEMENT ] :"<<rVariable<<std::endl;

	    double alpha = 1; //[0,1] //smoothing level of the Jacobian	      

            std::vector<double> Jacobians(1);
	    std::vector<double> InitialJacobians(1);

	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
		
		Geometry<Node<3> >& elem_geometry = (ie)->GetGeometry();
		Jacobians[0]    = 0;
		double Jacobian = 0;
		double size = elem_geometry.size();
		for( unsigned int node_i = 0 ; node_i <elem_geometry.size(); node_i++)
		  {
		    Jacobian = elem_geometry[node_i].FastGetSolutionStepValue(rVariable);
		    Jacobians[0] += (Jacobian/size);
		  }


		ie->GetValueOnIntegrationPoints(rVariable,InitialJacobians,CurrentProcessInfo);
		//std::cout<<" Element: "<<ie->Id()<<" Jacobian 0: "<<InitialJacobians[0]<<" Jacobian 1: "<<Jacobians[0]<<std::endl;

		Jacobians[0] = (alpha)*Jacobians[0]+(1-alpha)*InitialJacobians[0];

		ie->SetValueOnIntegrationPoints(rVariable,Jacobians,CurrentProcessInfo);
		
	      }

	    if( GetEchoLevel() > 0 )		    
	      std::cout<<" [ Finished NODE to ELEMENT Transfer ]"<<std::endl;
	    
	    KRATOS_CATCH( "" )
	}
      



      //*******************************************************************************************
      //*******************************************************************************************
	
      void MeshDataTransferUtilities::TransferNodalValuesToElements(const Variable<double>& rVariable,
								    const Variable<double>& rCriticalVariable,
								    const double& rCriticalValue,
								    ModelPart& rModelPart,
								    ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
		    
	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ Data Transfer NODE to ELEMENT ] :"<<rVariable<<" based on critical values of "<<rCriticalVariable<<std::endl;

	    double alpha = 0.5; //[0,1] //smoothing level of the Jacobian	      

            std::vector<double> Jacobians(1);
	    std::vector<double> InitialJacobians(1);
	    double Jacobian = 0;	  
            std::vector<double> ComputedValues(1);
	    double computed_value=0;
	    double critical_value=0;
	    int counter = 0;

	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
		
		ie->GetValueOnIntegrationPoints(rCriticalVariable,ComputedValues,CurrentProcessInfo);
	    
		computed_value = ComputedValues[0] * ie->GetGeometry().Area();
		critical_value = rCriticalValue;

		if( computed_value > critical_value )
		  {
		    for(unsigned int i = 0; i<ie->GetGeometry().size(); i++)
		      {
			ie->GetGeometry()[i].Set(TO_REFINE);
		      }
		  }
	      }

	     
	    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(MeshId); ie != rModelPart.ElementsEnd(MeshId); ie++)
	      {
		  		
		  Geometry<Node<3> >& elem_geometry = (ie)->GetGeometry();
		  double size = elem_geometry.size();

		  Jacobians[0]    = 0;
		  bool apply_smoothing = false;

		  for( unsigned int node_i = 0 ; node_i <elem_geometry.size(); node_i++)
		    {
		      Jacobian = elem_geometry[node_i].FastGetSolutionStepValue(rVariable);
		      Jacobians[0] += (Jacobian/size);
		      
		      if(elem_geometry[node_i].Is(TO_REFINE))
			apply_smoothing = true;
		    }
		
		  if( apply_smoothing ){

		    ie->GetValueOnIntegrationPoints(rVariable,InitialJacobians,CurrentProcessInfo);
		    //std::cout<<" Element: "<<ie->Id()<<" Jacobian 0: "<<InitialJacobians[0]<<" Jacobian 1: "<<Jacobians[0]<<std::endl;
		    Jacobians[0] = (alpha)*Jacobians[0]+(1-alpha)*InitialJacobians[0];

		    ie->SetValueOnIntegrationPoints(rVariable,Jacobians,CurrentProcessInfo);

		    counter ++;
		  }
	      }

		    
	    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(MeshId); in!=rModelPart.NodesEnd(MeshId); in++)
	      {
		in->Reset(TO_REFINE);
	      }

	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ Finished NODE to ELEMENT Transfer ] : ( Performed "<<counter<<" transfers of "<<rModelPart.NumberOfElements(MeshId)<<" possible )"<<std::endl;
	    
	    KRATOS_CATCH( "" )
	}
      


	//*******************************************************************************************
	//*******************************************************************************************
        void MeshDataTransferUtilities::TransferElementalValuesToNodes( const Variable<double>& rVariable,
							       ModelPart& rModelPart,
							       ModelPart::IndexType MeshId)
	{

	    KRATOS_TRY
	     
	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ Data Transfer ELEMENT to NODE ] :"<<rVariable<<std::endl;

	    std::vector<double> Jacobians(1);					
	    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
	    NodesContainerType& rNodes = rModelPart.Nodes(MeshId);

	    for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	      {
		if(in->IsNot(RIGID)){
		  WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);
		  double Jacobian  = 0;
		  double Area      = 0;
		  double elem_area = 0;
		
		  for(unsigned int i = 0; i < neighb_elems.size(); i++)
		    {
		      neighb_elems[i].GetValueOnIntegrationPoints(rVariable,Jacobians,CurrentProcessInfo);
		      Geometry<Node<3> >& elem_geometry = neighb_elems[i].GetGeometry();
		    
		      elem_area = elem_geometry.Area();

		      Area += elem_area;
		      Jacobian += Jacobians[0]*elem_area;

		    }
		
		  if(Area!=0 || Jacobian<=0)
		    Jacobian /= Area;
		  else
		    std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : Area = 0 (neighbours: "<<neighb_elems.size()<<") "<<std::endl;

		  if( in->SolutionStepsDataHas(rVariable))
		    in->FastGetSolutionStepValue(rVariable) = Jacobian;
		  else
		    std::cout<<" ERROR TR: Something Wrong in node ["<<in->Id()<<"] : variable "<<rVariable<<" was not defined "<<std::endl;
		    
		}

	      }
    

	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ Finished ELEMENT to NODE Transfer ]"<<std::endl;
		   
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
		    
	    if( GetEchoLevel() > 0 ){
	      std::cout<<" [ Data Transfer NODE to ELEMENT: NOT IMPLEMENTED YET ]"<<std::endl;
	      std::cout<<" [ Finished NODE to ELEMENT Transfer: NOT IMPLEMENTED YET ]"<<std::endl;
	    }

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

	    if( GetEchoLevel() > 0 ){   
	      std::cout<<" [ Data Transfer ELEMENT to NODE: NOT IMPLEMENTED YET ]"<<std::endl;
	      std::cout<<" [ Finished ELEMENT to NODE Transfer: NOT IMPLEMENTED YET ]"<<std::endl;
	    }

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

	    if( GetEchoLevel() > 0 )
	      std::cout<<" [ Data Transfer ELEMENT to ELEMENT ]"<<std::endl;

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


	    if( GetEchoLevel() > 0 ){
	      std::cout<<" [ MOVED TRANSFERS: "<<moved_transfers<<" ]"<<std::endl;
	      std::cout<<" [ Finished ELEMENT to ELEMENT Transfer ]"<<std::endl;
	    }

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

	void MeshDataTransferUtilities::Interpolate( Geometry<Node<3> >& geom,
					    const array_1d<double,3>& N,
					    unsigned int step_data_size,
					    Node<3>::Pointer pnode)
	{
	    unsigned int buffer_size = pnode->GetBufferSize();

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
		    step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
		}
	    }

	    if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0)
		KRATOS_THROW_ERROR( std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR", "" )

	}

	//*******************************************************************************************
	//*******************************************************************************************
  
        VariablesListDataValueContainer MeshDataTransferUtilities::InterpolateVariables( Triangle2D3<Node<3> >& geom,
											 const array_1d<double,3>& N,
											 unsigned int step_data_size,
											 Node<3>::Pointer pnode)
	  
	{
	  unsigned int buffer_size = pnode->GetBufferSize();
	 
	  //Copy Variables List
	  VariablesListDataValueContainer VariablesListData = (pnode)->SolutionStepData();
	 
	  for(unsigned int step = 0; step<buffer_size; step++)
	    {
	      //getting the data of the solution step
	      double* step_data  = VariablesListData.Data(step);
	     
	      double* node0_data = geom[0].SolutionStepData().Data(step);
	      double* node1_data = geom[1].SolutionStepData().Data(step);
	      double* node2_data = geom[2].SolutionStepData().Data(step);
	     
	      //copying this data in the position of the vector we are interested in
	      for(unsigned int j= 0; j<step_data_size; j++)
		{
		  step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
		}
	    }

	  if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0)
	    KRATOS_THROW_ERROR( std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR", "" )

	  return VariablesListData;
	}


}  // namespace Kratos.

