//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes


// Project includes
#include "custom_utilities/modeler_utilities.hpp"

namespace Kratos
{

  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::SetDomainLabels (ModelPart& rModelPart)
  {

    unsigned int start=0;
    unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
    if(NumberOfMeshes>1) 
      start=1;
      

    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(MeshId) ; i_node != rModelPart.NodesEnd(MeshId) ; i_node++)
	  {

	    i_node->SetValue(DOMAIN_LABEL,MeshId);
	  }
      }

  }

  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::BuildTotalMesh (ModelPart& rModelPart)
  {
    //Mesh Id=0

    std::cout<<" [ OLD TOTAL MESH (Elements: "<<rModelPart.NumberOfElements()<<" Nodes: "<<rModelPart.NumberOfNodes()<<" Conditions: "<<rModelPart.NumberOfConditions()<<" ] "<<std::endl;

    rModelPart.Nodes().clear();
    rModelPart.Elements().clear();

    //contact conditions are located on Mesh_0
    ModelPart::ConditionsContainerType KeepConditions;


    //std::cout<<" [ Number of Meshes "<<rModelPart.GetMeshes().size()-1<<" ]"<<std::endl;

    unsigned int nodeId=1;
    unsigned int elemId=1;
    unsigned int condId=1;

    unsigned int start=0;
    unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
    if(NumberOfMeshes>1) 
      start=1;
      

    for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
      {
	std::cout<<" [ CHILD MESH : ["<<MeshId<<"] (Elements: "<<rModelPart.NumberOfElements(MeshId)<<" Nodes: "<<rModelPart.NumberOfNodes(MeshId)<<" Conditions: "<<rModelPart.NumberOfConditions(MeshId)<<" ] "<<std::endl;


	for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(MeshId) ; i_elem != rModelPart.ElementsEnd(MeshId) ; i_elem++)
	  {
	    PointsArrayType& vertices=i_elem->GetGeometry().Points();
	    for(unsigned int i=0; i<vertices.size(); i++)
	      {
		vertices[i].Set(ENGAGED);
	      }

	    (rModelPart.Elements()).push_back(*(i_elem.base()));	
	    rModelPart.Elements().back().SetId(elemId);
	    elemId+=1;
	  }



	//Clean Nodes when redefining the total mesh:
	const array_1d<double,3> ZeroNormal(3,0.0);
	ModelPart::NodesContainerType temporal_nodes;
	temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	temporal_nodes.swap(rModelPart.Nodes(MeshId));
	  
	for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	  {
	    //i_node->PrintInfo(std::cout);
	    //std::cout<<std::endl;
	    if(i_node->Is(ENGAGED)){
		
	      i_node->Reset(INSERTED); //reset here if the node is labeled as insert 
	      i_node->Reset(REFINE); //reset here if the node is labeled as refine (to not duplicate boundary conditions)
	      i_node->Reset(ENGAGED); 

	      (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
	      (rModelPart.Nodes()).push_back(*(i_node.base()));	
	      rModelPart.Nodes().back().SetId(nodeId);
	      nodeId+=1;

	    }
	    else if (!rModelPart.NumberOfElements(MeshId)){
		
	      (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
	      (rModelPart.Nodes()).push_back(*(i_node.base()));	
	      rModelPart.Nodes().back().SetId(nodeId);
	      nodeId+=1;
	    }
	    else{
	      std::cout<<" NOT ENGAGED NODE "<<i_node->Id()<<std::endl;
	    }
		
	    if(i_node->IsNot(BOUNDARY))
	      {
		noalias(i_node->GetSolutionStepValue(FORCE_CONTACT_NORMAL)) = ZeroNormal;
		noalias(i_node->GetSolutionStepValue(FORCE_CONTACT_TANGENT)) = ZeroNormal;
	      }


	  }
	  
	//rModelPart.Nodes(MeshId).Sort();  

	for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(MeshId) ; i_cond != rModelPart.ConditionsEnd(MeshId) ; i_cond++)
	  {
	    i_cond->Reset(INSERTED); //reset here if the node is inserted
	    KeepConditions.push_back(*(i_cond.base()));
	    KeepConditions.back().SetId(condId);
	    condId+=1;	
	  }
      }


    for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
      {
	if(i_cond->Is(CONTACT)){
	  KeepConditions.push_back(*(i_cond.base()));
	  KeepConditions.back().SetId(condId);
	  condId+=1;
	}
      }
      
    rModelPart.Conditions().swap(KeepConditions);

    //Sort
    rModelPart.Nodes().Sort();
    rModelPart.Elements().Sort();
    rModelPart.Conditions().Sort();
      
    //Unique
    rModelPart.Nodes().Unique();
    rModelPart.Elements().Unique();
    rModelPart.Conditions().Unique();
      
    std::cout<<" [ NEW TOTAL MESH (Elements: "<<rModelPart.NumberOfElements()<<" Nodes: "<<rModelPart.NumberOfNodes()<<" Conditions: "<<rModelPart.NumberOfConditions()<<" ] "<<std::endl;      
 
  }
  
  //*******************************************************************************************
  //*******************************************************************************************


  void ModelerUtilities::CleanRemovedNodes(ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    KRATOS_TRY;

    //MESH 0 total domain mesh
    ModelPart::NodesContainerType temporal_nodes;
    temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	
    temporal_nodes.swap(rModelPart.Nodes(MeshId));
	    
    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
      {
	if( i_node->IsNot(RELEASE) ){
	  (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));	
	}
	else{
	  if( i_node->Is(BOUNDARY) )
	    std::cout<<" BOUNDARY NODE RELEASED "<<i_node->Id()<<std::endl;
	}
      }
	
	
    rModelPart.Nodes(MeshId).Sort();
	

    KRATOS_CATCH("")
      };


  //*******************************************************************************************
  //*******************************************************************************************

 
  bool ModelerUtilities::CheckSubdomain(Geometry<Node<3> >& vertices)
  {

    unsigned int DomainLabel = vertices[0].GetValue(DOMAIN_LABEL); //DOMAIN_LABEL must be set as nodal variable
      
    int samesbd=0;
      
    for(int pn=1; pn<3; pn++)
      {
	if(DomainLabel!=vertices[pn].GetValue(DOMAIN_LABEL))
	  {
	    samesbd++;
	  }
      }
      
      
    if(samesbd>0)
      return false;

    return true;
  }
  

  //*******************************************************************************************
  //*******************************************************************************************


  bool ModelerUtilities::CheckInnerCentre(Geometry<Node<3> >& vertices)
  {
    int nbound=0;

    int Nn=vertices.size();

    for(int pn=0; pn<Nn; pn++)
      {
	if(vertices[pn].Is(BOUNDARY))
	  {
	    nbound++;
	  }
      }


    if(nbound==3)
      {

	//std::cout<<" checked "<<std::endl;
	//Triangle2D3<Node<3> > Geometry (vertices[0], vertices[1], vertices[2]);
	//Center=Geometry.GetBariCentre(Center);

	//Baricenter
	array_1d<double, 3> Center;
	array_1d<double, 3> Normal;

	Center[0] = vertices[0].X() + vertices[1].X() + vertices[2].X();
	Center[1] = vertices[0].Y() + vertices[1].Y() + vertices[2].Y();
	Center[2] = 0;
	Center /=3.0;
	    
	bool outer=true;
	double tolerance = 0.05;
	int numouter=0;

	for(int v=0; v<Nn; v++)
	  {

	    array_1d<double, 3> Corner;
	    Corner[0] = vertices[v].X();
	    Corner[1] = vertices[v].Y();
	    Corner[2] = 0;

	    Normal=vertices[v].FastGetSolutionStepValue(NORMAL); 

	    if(norm_2(Normal))
	      Normal  /= norm_2(Normal);

	    //change position to be the vector from the vertex to the geometry center
	    Corner = Center-Corner;
	    if(norm_2(Corner))
	      Corner/= norm_2(Corner);

	    double projection=inner_prod(Corner,Normal);

	    if(projection>tolerance)
	      {
		numouter++;
	      }
	  }


	if(numouter>0)
	  outer=false;

	return outer; //if is outside to the  body domain returns false
      }


    return true;

  }


  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckOuterCentre(Geometry<Node<3> >& vertices,double &offset_factor)
  {
    int nbound=0;
    bool outer=false;

    int Nn=vertices.size();

    for(int pn=0; pn<Nn; pn++)
      {
	if(vertices[pn].Is(BOUNDARY))
	  {
	    nbound++;
	  }
      }


    if(nbound==3)
      {

	//Triangle2D3<Node<3> > Geometry (vertices[0], vertices[1], vertices[2]);
	//Center=Geometry.GetBariCentre(Center);

	//Baricenter
	array_1d<double, 3> Center;

	array_1d<double, 3> Normal =vertices[0].FastGetSolutionStepValue(NORMAL); 
	double  Shrink             =vertices[0].FastGetSolutionStepValue(SHRINK_FACTOR);   

	double x0 = vertices[0].X()-Normal[0]*Shrink*offset_factor;
	double y0 = vertices[0].Y()-Normal[1]*Shrink*offset_factor;
	    
	Normal = vertices[1].FastGetSolutionStepValue(NORMAL); 
	Shrink = vertices[1].FastGetSolutionStepValue(SHRINK_FACTOR);   

	double x1 = vertices[1].X()-Normal[0]*Shrink*offset_factor;
	double y1 = vertices[1].Y()-Normal[1]*Shrink*offset_factor;

	Normal = vertices[2].FastGetSolutionStepValue(NORMAL); 
	Shrink = vertices[2].FastGetSolutionStepValue(SHRINK_FACTOR);   
	    

	double x2 = vertices[2].X()-Normal[0]*Shrink*offset_factor;
	double y2 = vertices[2].Y()-Normal[1]*Shrink*offset_factor;

	Center[0] = x0 + x1 + x2;
	Center[1] = y0 + y1 + y2;
	Center[2] = 0;
	Center /=3.0;
	    
	double ortho  = 0.15;
	double slope  = 0.25;  //error assumed for some elements in the corners < 45 degrees
	double extra  = 0.95;

	int numouter=0;
	int numextra=0;
	int numcoplanar=0;
	int numsamedirection=0;
	int numorthogonal=0;

	array_1d<double, 3> Coplanar = vertices[0].FastGetSolutionStepValue(NORMAL); 

	if(norm_2(Coplanar))
	  Coplanar/= norm_2(Coplanar);

	//std::cout<<" vertices contact element: [normal(0)] "<<Coplanar<<" ";
	for(int v=0; v<Nn; v++)
	  {
	    Normal.clear();

	    array_1d<double, 3> Corner;
	    Normal = vertices[v].FastGetSolutionStepValue(NORMAL); 
	    Shrink = vertices[v].FastGetSolutionStepValue(SHRINK_FACTOR);   
				
	    Corner[0] = vertices[v].X()-Normal[0]*Shrink*offset_factor;
	    Corner[1] = vertices[v].Y()-Normal[1]*Shrink*offset_factor;
	    Corner[2] = 0;

	    if(norm_2(Normal))
	      Normal  /= norm_2(Normal);

	    //change position to be the vector from the vertex to the geometry center
	    Corner =Corner-Center;

	    if(norm_2(Corner))
	      Corner/= norm_2(Corner);

	    double projection=inner_prod(Corner,Normal);

	    if(projection<slope)
	      {
		numouter++;
	      }
	    else
	      {
		if(projection<extra)
		  numextra++;
	      }
		
	    double coplanar  =inner_prod(Coplanar,Normal);

	    //std::cout<<" V["<<v<<"]: "<<vertices[v]<<" Normal:"<<Normal<<" coplanar "<<fabs(coplanar)<<" < "<<ortho<<std::endl;

	    if(coplanar>0){
	      numsamedirection++;
	    }
		
	    if(coplanar>extra){
	      numcoplanar++;
	    }
		
	    if(fabs(coplanar)<=ortho){
	      numorthogonal++;
	    }

	  }

	// std::cout<<std::endl;
	// std::cout<<"  [ no:"<<numouter<<";ne:"<<numextra<<";nc:"<<numcoplanar<<";ns: "<<numsamedirection<<";no:"<<numorthogonal<<"]"<<std::endl;

	//std::cout<<" numouter "<<numouter<<" numextra "<<numextra<<std::endl;
	if(numouter==3)
	  outer=true;

	if(numouter==2 && numextra==1)
	  outer=true;

	if(numcoplanar==3)
	  outer=false;
	    
	if(numsamedirection==3)
	  outer=false;
	    
	// if(numorthogonal>=1)
	//   outer=false;
	    

      }

    return outer; //if is out to the  body domain returns true

  };

  
  //*******************************************************************************************
  //*******************************************************************************************

  double ModelerUtilities::FindBoundaryH (Node<3>& BoundaryPoint)
  {

    double havg = 0.00;
      
    if((BoundaryPoint.GetValue(NEIGHBOUR_NODES)).size() != 0)
      {
	double xc = BoundaryPoint.X();
	double yc = BoundaryPoint.Y();
	double zc = BoundaryPoint.Z();

	double h_nodes = 0;
	double h = 1000.0;
	for( WeakPointerVector< Node<3> >::iterator i = BoundaryPoint.GetValue(NEIGHBOUR_NODES).begin();
	     i !=  BoundaryPoint.GetValue(NEIGHBOUR_NODES).end(); i++)
	  {
	    if( i->Is(BOUNDARY) ){
	      double x = i->X();
	      double y = i->Y();
	      double z = i->Z();
	      double l = (x-xc)*(x-xc);
	      l += (y-yc)*(y-yc);
	      l += (z-zc)*(z-zc);
				    
	      if(l<h) h = l;
		
	      h = sqrt(h);
	      havg += h;
	      h_nodes += 1;
	    }
	  }
			      
	//calculate average h
	if(h_nodes == 0)
	  KRATOS_ERROR(std::logic_error,"no node has neighbours!!!!","");
			      
	havg /= h_nodes;
      }
			  
    return havg;
  };


  //*******************************************************************************************
  //*******************************************************************************************

  //returns false if it should be removed
  bool ModelerUtilities::ShrankAlphaShape(double alpha_param, Geometry<Node<3> >& vertices,double & offset_factor)
  {
    KRATOS_TRY

      boost::numeric::ublas::bounded_matrix<double,2,2> mJ; //local jacobian
    boost::numeric::ublas::bounded_matrix<double,2,2> mJinv; //inverse jacobian
    array_1d<double,2> mC; //center pos
    array_1d<double,2> mRhs; //center pos
	
    array_1d<double, 3>  Normal=vertices[0].FastGetSolutionStepValue(NORMAL); 
    double  Shrink             =vertices[0].FastGetSolutionStepValue(SHRINK_FACTOR);   

    //if normal not normalized
    //Shrink /=norm_2(Normal);

    double x0 = vertices[0].X()-Normal[0]*Shrink*offset_factor;
    double y0 = vertices[0].Y()-Normal[1]*Shrink*offset_factor;

    Normal = vertices[1].FastGetSolutionStepValue(NORMAL); 
    Shrink = vertices[1].FastGetSolutionStepValue(SHRINK_FACTOR);   

    double x1 = vertices[1].X()-Normal[0]*Shrink*offset_factor;
    double y1 = vertices[1].Y()-Normal[1]*Shrink*offset_factor;

    Normal = vertices[2].FastGetSolutionStepValue(NORMAL); 
    Shrink = vertices[2].FastGetSolutionStepValue(SHRINK_FACTOR);   

    double x2 = vertices[2].X()-Normal[0]*Shrink*offset_factor;
    double y2 = vertices[2].Y()-Normal[1]*Shrink*offset_factor;

    mJ(0,0)=2.0*(x1-x0);
    mJ(0,1)=2.0*(y1-y0);
    mJ(1,0)=2.0*(x2-x0);
    mJ(1,1)=2.0*(y2-y0);

    double detJ = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0);

    mJinv(0,0) =  mJ(1,1);
    mJinv(0,1) = -mJ(0,1);
    mJinv(1,0) = -mJ(1,0);
    mJinv(1,1) =  mJ(0,0);

    bounded_matrix<double,2,2> check;

    //calculate average h
    double h;
    h =  vertices[0].FastGetSolutionStepValue(NODAL_H);
    h += vertices[1].FastGetSolutionStepValue(NODAL_H);
    h += vertices[2].FastGetSolutionStepValue(NODAL_H);
    h *= 0.333333333;

    //calculate average h of boundary faces
    double h_side = FindBoundaryH(vertices[0]) + FindBoundaryH(vertices[1]) + FindBoundaryH(vertices[2]);
    h_side *= 0.333333333;   
    
    if(h_side>h)
      h = h_side;

    // if( h > h_side ){
    //   //std::cout<<" hside "<<h_side<<" h "<<h<<"  h bigger  "<<std::endl;
    //   if( h > 5 * h_side ){
    // 	h = h_side * 5;
    //   }
    //   else{
    // 	h = h_side;
    //   }
    // }
    // else{
    //   //std::cout<<" hside "<<h_side<<" h "<<h<<"  h_side bigger  "<<std::endl;
    //   if( h_side > 5 * h ){
    // 	h *= 5;
    //   }
    //   else{
    // 	h  = h_side;
    //   }
    // }
         
    if(detJ < 1e-6*h*h) 
      {
	// std::cout << "Alpha shape: volume criterion  = " << detJ << std::endl;
	////mark as boundary
	vertices[0].Set(BOUNDARY);
	vertices[1].Set(BOUNDARY);
	vertices[2].Set(BOUNDARY);
	return false;
	    
      }
    else
      {

	double x0_2 = x0*x0 + y0*y0;
	double x1_2 = x1*x1 + y1*y1;
	double x2_2 = x2*x2 + y2*y2;
	    
	//finalizing the calculation of the inverted matrix
	//std::cout<<"MATR INV"<<MatrixInverse(mJ)<<std::endl;
	mJinv /= detJ;
	//calculating the RHS
	mRhs[0] = (x1_2 - x0_2);
	mRhs[1] = (x2_2 - x0_2);
	    
	//calculate position of the center
	noalias(mC) = prod(mJinv,mRhs);
	    
	double radius = sqrt(pow(mC[0]-x0,2)+pow(mC[1]-y0,2));
	    
	// std::cout<<" CircumCenter "<<mC<<std::endl;
	    
	// std::cout<<" Shrank Alpha Shape "<<radius<<" < "<<h*alpha_param<<" ";//<<std::endl;

	double extra_alpha=1.4;

	if (radius < h*alpha_param*extra_alpha)
	  {
	    //std::cout<<" TRUE "<<std::endl;
	    return true;
	  }
	else
	  {
	    //std::cout<<" FALSE "<<std::endl;
	    return false;
	  }
      }


    KRATOS_CATCH("")
      };



  //*******************************************************************************************
  //*******************************************************************************************

  //returns false if it should be removed
  bool ModelerUtilities::AlphaShape(double alpha_param, Geometry<Node<3> >& vertices)
  {
    KRATOS_TRY

      boost::numeric::ublas::bounded_matrix<double,2,2> mJ; //local jacobian
    boost::numeric::ublas::bounded_matrix<double,2,2> mJinv; //inverse jacobian
    array_1d<double,2> mC; //center pos
    array_1d<double,2> mRhs; //center pos

    double x0 = vertices[0].X();
    double x1 = vertices[1].X();
    double x2 = vertices[2].X();

    double y0 = vertices[0].Y();
    double y1 = vertices[1].Y();
    double y2 = vertices[2].Y();

    mJ(0,0)=2.0*(x1-x0);
    mJ(0,1)=2.0*(y1-y0);
    mJ(1,0)=2.0*(x2-x0);
    mJ(1,1)=2.0*(y2-y0);

    double detJ = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0);

    mJinv(0,0) =  mJ(1,1);
    mJinv(0,1) = -mJ(0,1);
    mJinv(1,0) = -mJ(1,0);
    mJinv(1,1) =  mJ(0,0);

    bounded_matrix<double,2,2> check;

    //calculate average h
    double h;
    h =  vertices[0].FastGetSolutionStepValue(NODAL_H);
    h += vertices[1].FastGetSolutionStepValue(NODAL_H);
    h += vertices[2].FastGetSolutionStepValue(NODAL_H);
    h *= 0.333333333;

    if(detJ < 5e-5*h*h)
      {
	////mark as boundary
	vertices[0].Set(BOUNDARY);
	vertices[1].Set(BOUNDARY);
	vertices[2].Set(BOUNDARY);
	return false;
      }
    else
      {

	double x0_2 = x0*x0 + y0*y0;
	double x1_2 = x1*x1 + y1*y1;
	double x2_2 = x2*x2 + y2*y2;

	//finalizing the calculation of the inverted matrix
	//std::cout<<"MATR INV"<<MatrixInverse(mJ)<<std::endl;
	mJinv /= detJ;
	//calculating the RHS
	mRhs[0] = (x1_2 - x0_2);
	mRhs[1] = (x2_2 - x0_2);

	//calculate position of the center
	noalias(mC) = prod(mJinv,mRhs);

	double radius = sqrt(pow(mC[0]-x0,2)+pow(mC[1]-y0,2));

	if (radius < h*alpha_param)
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }
      }


    KRATOS_CATCH("")
      };


  //*******************************************************************************************
  //*******************************************************************************************

  void ModelerUtilities::CheckParticles (ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    int NumberOfNodes = rModelPart.NumberOfNodes();
    std::cout<<" Number of Nodes "<<NumberOfNodes<<std::endl;
    for(int id=1; id<=NumberOfNodes; id++)
      {
	std::cout<<" Check node: "<<id<<std::endl;
	if(rModelPart.Nodes()[id].Is(BOUNDARY)){
	  std::cout<<" Node : "<<id<<" is boundary "<<std::endl;
	}
	else{
	  std::cout<<" Node : "<<id<<" is not boundary "<<std::endl;
	}  

      }
  }


  //*******************************************************************************************
  //*******************************************************************************************
  
  // inline double ModelerUtilities::CalculateSideLength(PointType& P1,PointType& P2)
  // {
  //   return sqrt( (P1.X()-P2.X())*(P1.X()-P2.X()) + (P1.Y()-P2.Y())*(P1.Y()-P2.Y()) );
  // };

  //*******************************************************************************************
  //*******************************************************************************************

  // inline double ModelerUtilities::CalculateCircRadius(Geometry< Node<3> >& rGeometry)
  // {
  //   double L1 = CalculateSideLength (rGeometry[0],rGeometry[1]);
  //   double L2 = CalculateSideLength (rGeometry[1],rGeometry[2]);
  //   double L3 = CalculateSideLength (rGeometry[2],rGeometry[0]);
      
      
  //   double Area = rGeometry.Area();
       

  //   double  Rcrit = Area*2/(L1+L2+L3);

  //   return Rcrit;

  // };

  //*******************************************************************************************
  //*******************************************************************************************

  // inline double ModelerUtilities::CalculateCircRadius(const double x0, const double y0,
  // 						      const double x1, const double y1,
  // 						      const double x2, const double y2,
  // 						      double& Area)
  // {

  //   double L1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  //   double L2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  //   double L3 = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0));     
      
  //   Area = fabs( 0.5 * ( (x0*y1) - (x0*y2) - (x1*y0) + (x1*y2) + (x2*y0) - (x2*y1) ) );
      
  //   // std::cout<< " Area "<<Area<<" L1 "<<L1<<" L2 "<<L2<<" L3 "<<L3<<std::endl;
      
  //   double  Rcrit = Area*2/(L1+L2+L3);

  //   return Rcrit;
  // };



  //*******************************************************************************************
  //*******************************************************************************************

  // inline double ModelerUtilities::CalculateAverageSideLength(const double x0, const double y0,
  // 							   const double x1, const double y1,
  // 							   const double x2, const double y2)
  // {
  //   double length_0 = sqrt( x0*x0 + y0*y0 );
  //   double length_1 = sqrt( x1*x1 + y1*y1 );
  //   double length_2 = sqrt( x2*x2 + y2*y2 );
      
  //   return 0.5*( length_0 + length_1 + length_2 );
  // };





  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckConditionInBox(Condition::Pointer& pCond, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo)
  {
    bool inside = true;
    Vector Point(2);
    Geometry< Node<3> >& rGeom = pCond->GetGeometry();

    for(unsigned int i=0; i<rGeom.size(); i++)
      {
	Point[0] = rGeom[i].X();
	Point[1] = rGeom[i].Y();
	if( !RefiningBox.IsInside(Point,CurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;
  };

  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckElementInBox(Element::Pointer& pElem, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo)
  {
    bool inside = true;
    Vector Point(2);
    Geometry< Node<3> >& rGeom = pElem->GetGeometry();

    for(unsigned int i=0; i<rGeom.size(); i++)
      {
	Point[0] = rGeom[i].X();
	Point[1] = rGeom[i].Y();
	if( !RefiningBox.IsInside(Point,CurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;
  };

  //*******************************************************************************************
  //*******************************************************************************************

  bool ModelerUtilities::CheckVerticesInBox(Geometry<Node<3> >& rGeom, SpatialBoundingBox& RefiningBox, ProcessInfo& CurrentProcessInfo)
  {
    bool inside = true;
    Vector Point(2);

    for(unsigned int i=0; i<rGeom.size(); i++)
      {
	Point[0] = rGeom[i].X();
	Point[1] = rGeom[i].Y();
	if( !RefiningBox.IsInside(Point,CurrentProcessInfo[TIME]) ){
	  inside = false;
	  break;
	}
      }

    return inside;
  };


  //*******************************************************************************************
  //*******************************************************************************************

  Condition::Pointer ModelerUtilities::FindMasterCondition(Condition::Pointer& pCond, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found)
  {
		
    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeom = pCond->GetGeometry();
    boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
    rGeom.NodesInFaces(lpofa);   
		
    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeom "<<rGeom<<std::endl;

    condition_found=false;
    for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      {
	//2D edges:
	if(ic->IsNot(CONTACT)){

	  Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
	  for(unsigned int i=0; i<lpofa.size2();i++)
	    {
	      // std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	      // std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeom[lpofa(1,i)].Id()<<"] ["<<rGeom[lpofa(2,i)].Id()<<"] "<<std::endl;

	      if( (   rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() 
		      && rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id() ) || 
		  (   rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() 
		      && rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id() ) )
		{
		  pMasterCondition= *(ic.base());
		  condition_found=true;
		  break;
		}
			       
	    }

	}
	if(condition_found)
	  {
	    break;
	  }
			
      }

    if(!condition_found) {
      //   //KRATOS_ERROR(std::logic_error, "Boundary Condition NOT FOUND after CONTACT MESHING SEARCH","")
      std::cout<<" Boundary Condition NOT FOUND after CONTACT MESHING SEARCH "<<std::endl;
      //   std::cout<<" rGeom "<<rGeom<<std::endl;

      //   for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      //     {
      //       //2D edges:
			
      //       Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
      for(unsigned int i=0; i<lpofa.size2();i++)
	{
	  // 	  std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	  std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeom[lpofa(1,i)].Id()<<"] ["<<rGeom[lpofa(2,i)].Id()<<"] "<<std::endl;
		
	}

      //     }
    }

    return pMasterCondition;
	
  };
  
  //*******************************************************************************************
  //*******************************************************************************************

  Condition::Pointer ModelerUtilities::FindMasterCondition(Condition::Pointer& pCond, PointType& pSlaveNode, ModelPart::ConditionsContainerType & rModelConditions,bool & condition_found)
  {
      
    Condition::Pointer pMasterCondition;

    Geometry< Node<3> >& rGeom = pCond->GetGeometry();
    boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces
    rGeom.NodesInFaces(lpofa);   
		
    //std::cout<<" lpofa "<<lpofa<<std::endl;
    //std::cout<<" rGeom "<<rGeom<<std::endl;

    condition_found=false;
    for(ModelPart::ConditionsContainerType::iterator ic=rModelConditions.begin(); ic!=rModelConditions.end(); ic++)
      {
	//2D edges:
	if(ic->IsNot(CONTACT)){
				      
	  Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
	  for(unsigned int i=0; i<lpofa.size2();i++)
	    {
	      // std::cout<<" General Conditions IDs ["<<rConditionGeom[0].Id()<<"] ["<<rConditionGeom[1].Id()<<"] "<<std::endl;
	      // std::cout<<" Local Conditions IDs ("<<i<<"):["<<rGeom[lpofa(1,i)].Id()<<"] ["<<rGeom[lpofa(2,i)].Id()<<"] "<<std::endl;

	      if( (   rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() 
		      && rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id() ) || 
		  (   rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() 
		      && rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id() ) )
		{
		  pMasterCondition = *(ic.base());
		  pSlaveNode = rGeom[lpofa(0,i)];
		  //std::cout<<"   Slave_Node: found: "<<rGeom[lpofa(0,i)].Id()<<std::endl;
		  condition_found=true;
		  break;
		}
			       
	    }
	}
	if(condition_found)
	  {
	    break;
	  }
			
      }

    // if(!found) 
    //     KRATOS_ERROR(std::logic_error, "Boundary Condition NOT FOUND after CONTACT MESHING SEARCH","")

    return pMasterCondition;
  };



} // Namespace Kratos

