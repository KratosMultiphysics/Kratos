//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#include "includes/model_part.h"

#if !defined(KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_H_INCLUDED )
#define  KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_H_INCLUDED


/* System includes */


/* External includes */


/* Project includes */
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/// Tool to evaluate the normals on nodes based on the normals of a set of surface conditions
class BoundaryNormalsCalculationUtilities
{
public:
	/**@name Type Definitions */
	/*@{ */
	typedef ModelPart::NodesContainerType               NodesArrayType;
        typedef ModelPart::ElementsContainerType     ElementsContainerType;
	typedef ModelPart::ConditionsContainerType ConditionsContainerType;
  	typedef ModelPart::MeshType                               MeshType;
	/*@} */
	/**@name Life Cycle
	 */
	/*@{ */

	/** Constructor.
	 */


	/** Destructor.
	 */

	/*@} */
	/**@name Operators
	 */
	/*@{ */


	/*@} */
	/**@name Operations */
	/*@{ */

        /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area).
        /** This is done on the base of the Conditions provided which should be
	 * understood as the surface elements of the area of interest.
	 * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
	 * @param dimension Spatial dimension (2 or 3)
	 * @note Use this fuction instead of its overload taking a Conditions array for MPI applications,
	 * as it will take care of communication between partitions.
	 */
 


	/// Calculates the area normal (unitary vector oriented as the normal).

        void CalculateUnitBoundaryNormals(ModelPart& rModelPart, int EchoLevel = 0)
	{

	  mEchoLevel = EchoLevel;
	  unsigned int start = 0;
	  unsigned int NumberOfMeshes = rModelPart.NumberOfMeshes();

	  if(NumberOfMeshes>1) 
	    start=1;

	  //By the way: set meshes options from bools
	  for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	    {
       
	      this->CalculateMeshUnitBoundaryNormals(rModelPart, MeshId, EchoLevel);
		
	      if(mEchoLevel > 1) 
		std::cout<<"   Unit Boundary Normals Set MESH ["<<MeshId<<"] "<<std::endl;			
	    }
	    

	  // For MPI: correct values on partition boundaries
	  rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

	}
	


	/// Calculates the area normal (unitary vector oriented as the normal) and weight the normal to shrink

        void CalculateBoundaryNormals(ModelPart& rModelPart, int EchoLevel = 0)
	{

	  mEchoLevel = EchoLevel;
	  unsigned int start=0;
	  unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
	  if(NumberOfMeshes>1) 
	    start=1;

	  //By the way: set meshes options from bools
	  for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	    {
       
	      this->CalculateMeshBoundaryNormals(rModelPart, MeshId, EchoLevel);
		
	      if(mEchoLevel > 1) 
		std::cout<<"   Boundary Normals Set MESH ["<<MeshId<<"] "<<std::endl;			
	    }
	
	  // For MPI: correct values on partition boundaries
	  rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
	  rModelPart.GetCommunicator().AssembleCurrentData(SHRINK_FACTOR);

	}
		

       /// Calculates the normals of the BOUNDARY for a mesh given by MeshId
       void CalculateMeshUnitBoundaryNormals(ModelPart& rModelPart, int MeshId, int EchoLevel = 0)
	{
	  KRATOS_TRY

	  mEchoLevel = EchoLevel;
		
	  const unsigned int dimension = (rModelPart.pGetMesh())->WorkingSpaceDimension();

	  if( dimension == 2 ) {
	    	
	    ConditionsContainerType& rConditions = rModelPart.Conditions(MeshId);
	    
	    if( rConditions.size() > 0 ){
	      //Reset Body Normal Variables:
	      //this->ResetBodyNormals(rModelPart,MeshId);
	      //Compute New Normals
	      this->CalculateMeshBoundaryNormals(rConditions, EchoLevel);
	    }
	    
	  }
	  else if ( dimension == 3 ){
	    
	    MeshType& rMesh = rModelPart.GetMesh(MeshId);
	    
	    //check which type of element is in the 3D mesh, a surface (false) or volumetric (true) element
	    if( this->CheckVolumetricElements(rMesh, EchoLevel) ){
	      //Reset Body Normal Variables:
	      //this->ResetBodyNormals(rModelPart,MeshId);
	      //Compute New Normals
	      this->CalculateMeshBoundaryNormals(rMesh, EchoLevel);
	    }
	    else{
	      //std::cout<<" Surface Elements Boundary Calculation "<<std::endl;
	      ElementsContainerType& rElements = rModelPart.Elements(MeshId);
	      if( rElements.size() > 0 ){
		//Reset Body Normal Variables:
		//this->ResetBodyNormals(rModelPart,MeshId);
		//Compute New Normals		
		this->CalculateMeshBoundaryNormals(rElements, EchoLevel);
	      }
	    }
	    
	  }
	  
	  //standard assignation // fails in sharp edges angle<90
	  AddNormalsToNodes(rModelPart.Conditions(MeshId));
		

	  KRATOS_CATCH( "" )

	}



       /// Calculates the normals of the BOUNDARY for a mesh given by MeshId
       void CalculateMeshBoundaryNormals(ModelPart& rModelPart, int MeshId, int EchoLevel = 0)
	{
	  KRATOS_TRY

	  mEchoLevel = EchoLevel;
	    		
	  const unsigned int dimension = (rModelPart.pGetMesh())->WorkingSpaceDimension();

	  if( dimension == 2 ) {
	    	
	    ConditionsContainerType& rConditions = rModelPart.Conditions(MeshId);
	    
	    if( rConditions.size() > 0 ){
	      //Reset Body Normal Variables:
	      this->ResetBodyNormals(rModelPart,MeshId); //Otherwise, some spurious normal appear in the interior of the domain due to interpolation
	      //Compute New Normals
	      this->CalculateMeshBoundaryNormals(rConditions, EchoLevel);
	    }

	  }
	  else if ( dimension == 3 ){
	    

	    ConditionsContainerType& rConditions = rModelPart.Conditions(MeshId);	    
	    ElementsContainerType& rElements = rModelPart.Elements(MeshId);

	    if( rConditions.size() > 0 ){
	      //std::cout<<"    Conditions: Surface Elements Boundary Calculation "<<std::endl;
	      //Reset Body Normal Variables:
	      //this->ResetBodyNormals(rModelPart,MeshId);
	      //Compute New Normals
	      this->CalculateMeshBoundaryNormals(rConditions, EchoLevel);
	    }
	    else if ( rElements.size() > 0 ){

	      MeshType& rMesh = rModelPart.GetMesh(MeshId);
	    
	      //check which type of element is in the 3D mesh, a surface (false) or volumetric (true) element
	      if( this->CheckVolumetricElements(rMesh, EchoLevel) ){
		//std::cout<<"    Volumetric Elements Boundary Calculation "<<std::endl;
		//Reset Body Normal Variables:
		//this->ResetBodyNormals(rModelPart,MeshId);
		//Compute New Normals
		this->CalculateMeshBoundaryNormals(rMesh, EchoLevel);
	      }
	      else{
		//std::cout<<"    Surface Elements Boundary Calculation "<<std::endl;
		ElementsContainerType& rElements = rModelPart.Elements(MeshId);
		if( rElements.size() > 0 ){
		  //Reset Body Normal Variables:
		  //this->ResetBodyNormals(rModelPart,MeshId);
		  //Compute New Normals		
		  this->CalculateMeshBoundaryNormals(rElements, EchoLevel);
		}
	      }

	    }
	  }
	  
	  //solid pfem assignation: Unity Normals on nodes and Shrink_Factor on nodes
	  AddNormalsAndShrink(rModelPart,dimension,MeshId);    

	  KRATOS_CATCH( "" )

	}
       


	/*@} */
	/**@name Acces */
	/*@{ */


	/*@} */
	/**@name Inquiry */
	/*@{ */


	/*@} */
	/**@name Friends */
	/*@{ */


	/*@} */

private:
	/**@name Static Member Variables */
	/*@{ */


	/*@} */
	/**@name Member Variables */
	/*@{ */
        int mEchoLevel;
  
  	/*@} */
	/**@name Private Operators*/
	/*@{ */


	/*@} */
	/**@name Private Operations*/
	/*@{ */


	//this function adds the Contribution of one of the geometries
	//to the corresponding nodes
	static void CalculateUnityNormal2D(ConditionsContainerType::iterator it, array_1d<double,3>& An)
	{
		Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

		An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
		An[1] = - (pGeometry[1].X() - pGeometry[0].X());
		An[2] =    0.00;

		array_1d<double,3>& normal = (it)->GetValue(NORMAL);
		noalias(normal) = An/norm_2(An);

		// 				(it)->SetValue(NORMAL,An);
	}

	static void CalculateUnityNormal3D(ConditionsContainerType::iterator it, array_1d<double,3>& An,
					   array_1d<double,3>& v1,array_1d<double,3>& v2 )
	{
		Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

		v1[0] = pGeometry[1].X() - pGeometry[0].X();
		v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
		v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

		v2[0] = pGeometry[2].X() - pGeometry[0].X();
		v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
		v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

		MathUtils<double>::CrossProduct(An,v1,v2);
		An *= 0.5;

		array_1d<double,3>& normal = (it)->GetValue(NORMAL);

		noalias(normal) = An/norm_2(An);
		//noalias((it)->GetValue(NORMAL)) = An;
	}


	//this function adds the Contribution of one of the geometries
	//to the corresponding nodes
	static void CalculateUnityNormal2D(ElementsContainerType::iterator it, array_1d<double,3>& An)
	{
		Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

		An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
		An[1] = - (pGeometry[1].X() - pGeometry[0].X());
		An[2] =    0.00;

		array_1d<double,3>& normal = (it)->GetValue(NORMAL);
		noalias(normal) = An/norm_2(An);

		// 				(it)->SetValue(NORMAL,An);
	}

	static void CalculateUnityNormal3D(ElementsContainerType::iterator it, array_1d<double,3>& An,
					   array_1d<double,3>& v1,array_1d<double,3>& v2 )
	{
		Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

		v1[0] = pGeometry[1].X() - pGeometry[0].X();
		v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
		v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

		v2[0] = pGeometry[2].X() - pGeometry[0].X();
		v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
		v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

		MathUtils<double>::CrossProduct(An,v1,v2);
		An *= 0.5;

		array_1d<double,3>& normal = (it)->GetValue(NORMAL);

		noalias(normal) = An/norm_2(An);
		// 				noalias((it)->GetValue(NORMAL)) = An;
	}


	void ResetBodyNormals(ModelPart& rModelPart, ModelPart::IndexType MeshId=0)
				      
	{
		KRATOS_TRY
		  
		//resetting the normals
		for(NodesArrayType::iterator in =  rModelPart.NodesBegin(MeshId);
		    in !=rModelPart.NodesEnd(MeshId); in++)
		{
		  (in->GetSolutionStepValue(NORMAL)).clear();
		}

		KRATOS_CATCH( "" )
	}

	void CheckBodyNormals(ModelPart& rModelPart, ModelPart::IndexType MeshId=0)
				      
	{
		KRATOS_TRY
		  
		//resetting the normals
		for(NodesArrayType::iterator in =  rModelPart.NodesBegin(MeshId);
		    in !=rModelPart.NodesEnd(MeshId); in++)
		{
		  std::cout<<" ID: "<<in->Id()<<" normal: "<<(in->GetSolutionStepValue(NORMAL))<<std::endl;
		}

		KRATOS_CATCH( "" )
	}

       /// Calculates the "area normal" (vector oriented as the normal with a dimension proportional to the area).
       /** This is done on the base of the Conditions provided which should be
        * understood as the surface elements of the area of interest.
        * @param rConditions A set of conditions defining the "skin" of a model
        * @param dimension Spatial dimension (2 or 3)
        * @note This function is not recommended for distributed (MPI) runs, as
        * the user has to ensure that the calculated normals are assembled between
        * processes. The overload of this function that takes a ModelPart is
        * preferable in this case, as it performs the required communication.
        */
        void CalculateMeshBoundaryNormals(ConditionsContainerType& rConditions, int EchoLevel = 0)
				      
	{
		KRATOS_TRY
		  
	        //resetting the normals
		for(ConditionsContainerType::iterator it =  rConditions.begin();
		    it != rConditions.end(); it++)
		{
			Element::GeometryType& rNodes = it->GetGeometry();
			for(unsigned int in = 0; in<rNodes.size(); in++)
 			    ((rNodes[in]).GetSolutionStepValue(NORMAL)).clear();
		}

		const unsigned int dimension = (rConditions.begin())->WorkingSpaceDimension();

		//calculating the normals and storing on the conditions
		array_1d<double,3> An;
		if(dimension == 2)
		{
			for(ConditionsContainerType::iterator it =  rConditions.begin();
			    it !=rConditions.end(); it++)
			{
			  if(it->IsNot(CONTACT) && it->Is(BOUNDARY) )
			    CalculateUnityNormal2D(it,An);
			}
		
		}
		else if(dimension == 3)
		{
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			for(ConditionsContainerType::iterator it =  rConditions.begin();
			    it !=rConditions.end(); it++)
			{
			  //calculate the normal on the given condition
			  if(it->IsNot(CONTACT) && it->Is(BOUNDARY))
			    CalculateUnityNormal3D(it,An,v1,v2);
			}
		}

		KRATOS_CATCH( "" )
	}


       /// Calculates the "area normal" (vector oriented as the normal with a dimension proportional to the area).
       /** This is done on the base of the Elements provided which should be
        * understood as the surface elements of the area of interest.
        * @param rElements A set of elements defining the "skin" of a model
        * @param dimension Spatial dimension (2 or 3)
        * @note This function is not recommended for distributed (MPI) runs, as
        * the user has to ensure that the calculated normals are assembled between
        * processes. The overload of this function that takes a ModelPart is
        * preferable in this case, as it performs the required communication.
        */
        void CalculateMeshBoundaryNormals(ElementsContainerType& rElements, int EchoLevel = 0)
				      
	{
		KRATOS_TRY

	        //resetting the normals
		for(ElementsContainerType::iterator it =  rElements.begin();
		    it != rElements.end(); it++)
		{
			Element::GeometryType& rNodes = it->GetGeometry();
			for(unsigned int in = 0; in<rNodes.size(); in++)
 			    ((rNodes[in]).GetSolutionStepValue(NORMAL)).clear();
		}

		const unsigned int dimension = (rElements.begin())->WorkingSpaceDimension();

		//calculating the normals and storing on the conditions
		array_1d<double,3> An;
		if(dimension == 2)
		{
			for(ElementsContainerType::iterator it =  rElements.begin();
			    it !=rElements.end(); it++)
			{
			  if(it->IsNot(CONTACT)){
			    it->Set(BOUNDARY);
			    CalculateUnityNormal2D(it,An);
			  }
			}
		}
		else if(dimension == 3)
		{
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			for(ElementsContainerType::iterator it =  rElements.begin();
			    it !=rElements.end(); it++)
			{
			  //calculate the normal on the given condition
			  if(it->IsNot(CONTACT)){
			    it->Set(BOUNDARY);
			    CalculateUnityNormal3D(it,An,v1,v2);
			  }
			}
		}

		KRATOS_CATCH( "" )
	}


       //Check if the mesh has volumetic elements
       bool CheckVolumetricElements(MeshType& rMesh, int EchoLevel = 0)
       {

	  KRATOS_TRY

	  ElementsContainerType& rElements = rMesh.Elements();
	 
	  ElementsContainerType::iterator it =  rElements.begin();
	  
	  if( (it)->GetGeometry().WorkingSpaceDimension() == 3 ){    
	    return true;
	  }
	  else{
	    return false;
	  }

	  KRATOS_CATCH( "" )
	 
       }
  
       /// Calculates the normals of the BOUNDARY nodes of a given MeshId
       //  using a consistent way: SOTO & CODINA
       //  fails in sharp edges angle<90
       void CalculateMeshBoundaryNormals(MeshType& rMesh, int EchoLevel = 0)
       {
	  KRATOS_TRY

	    
	  //Reset normals
	  ModelPart::NodesContainerType&    rNodes = rMesh.Nodes();
	  ModelPart::ElementsContainerType& rElems = rMesh.Elements();
	  
	  bool neighsearch=true;

	  //Check if the neigbours search is already done and set
	  unsigned int number_of_nodes = rElems.begin()->GetGeometry().PointsNumber();
	  for(unsigned int i=0; i<number_of_nodes; i++)
	    if( (rElems.begin()->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS)).size() > 1 )
	      neighsearch=false;
	  
	  std::cout<<" Neighbour Search "<<neighsearch<<std::endl;
	  
	  for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	    {
	      (in->GetSolutionStepValue(NORMAL)).clear();

	      if(neighsearch){
		//*************  Neigbours of nodes search  ************//
		WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
		rE.erase(rE.begin(),rE.end() );
		in->Reset(BOUNDARY);
		//*************  Neigbours of nodes search ************//
	      }
	    }
	  
	  
	  if(neighsearch){
	    //*************  Neigbours of nodes search ************//
	    //add the neighbour elements to all the nodes in the mesh
	    for(ModelPart::ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
	      {
		Element::GeometryType& pGeom = ie->GetGeometry();
		for(unsigned int i = 0; i < pGeom.size(); i++)
		  {
		    (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
		  }
	      }
	    //*************  Neigbours of nodes search ************//
	  }

	  //calculating the normals and storing it on nodes

	  Vector An(3);
	  Element::IntegrationMethod mIntegrationMethod = Element::GeometryDataType::GI_GAUSS_1; //one gauss point
	  int PointNumber = 0; //one gauss point
	  unsigned int dimension;
	  Matrix J;
	  Matrix InvJ;
	  double detJ;
	  Matrix  DN_DX;
	  
	  int not_assigned   = 0;
	  int boundary_nodes = 0;
	  //int boundarycounter=0;
	  for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
	    {
	      An.clear();

	      //if(in->Is(BOUNDARY)){

	      WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);

	      for(WeakPointerVector<Element >::iterator ie= rE.begin(); ie!=rE.end(); ie++)
		{
		  Element::GeometryType& rGeom = ie->GetGeometry();

		  if( rGeom.EdgesNumber() > 1 ){

		    dimension = rGeom.WorkingSpaceDimension();

		    //********** Compute the element integral ******//
		    const Element::GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( mIntegrationMethod );
		    const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( mIntegrationMethod );


		    J.resize( dimension, dimension );
		    J = rGeom.Jacobian( J, PointNumber , mIntegrationMethod );

		    InvJ.clear();
		    detJ=0;
		    //Calculating the inverse of the jacobian and the parameters needed
		    MathUtils<double>::InvertMatrix( J, InvJ, detJ);
		  

		    //Compute cartesian derivatives for one gauss point
		    DN_DX = prod( DN_De[PointNumber] , InvJ );

		    double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;

		  
		    for(unsigned int i = 0; i < rGeom.size(); i++)
		      {
			if(in->Id() == rGeom[i].Id()){

			  for(unsigned int d=0; d<dimension; d++)
			    {		    
			      An[d] += DN_DX(i,d)*IntegrationWeight;
			    }
			}
		      }
		  
		    //********** Compute the element integral ******//

		  }

		  if(norm_2(An)>1e-12){
		    noalias(in->FastGetSolutionStepValue(NORMAL)) = An/norm_2(An);
		    if(neighsearch){
		      in->Set(BOUNDARY);
		    }
		  }
		  else{
		    (in->FastGetSolutionStepValue(NORMAL)).clear();
		    //std::cout<<" ERROR: normal not set "<<std::endl;
		    not_assigned +=1;
		    if(neighsearch){
		      in->Set(BOUNDARY,false);
		    }
		  }

		}
	    
	      if(in->Is(BOUNDARY))
		boundary_nodes +=1;

	      //boundarycounter++;
	      //}
	    }


	  if(mEchoLevel > 0) 
	    std::cout<<"  [ Boundary_Normals  (Mesh Nodes:"<<rNodes.size()<<")[SET:"<<boundary_nodes<<" / NOT_SET:"<<not_assigned<<"] ]"<<std::endl;

	  
	  //std::cout<<" Boundary COUNTER "<<boundarycounter<<std::endl;
	  KRATOS_CATCH( "" )

	}

        //*****************************

	void AddNormalsToNodes(ConditionsContainerType& rConditions)
	{
			
		//adding the normals to the nodes
		for(ConditionsContainerType::iterator it =  rConditions.begin();
		    it !=rConditions.end(); it++)
		{
			Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
			double coeff = 1.00/pGeometry.size();
			const array_1d<double,3>& An = it->GetValue(NORMAL);

			for(unsigned int i = 0; i<pGeometry.size(); i++)
			{
				noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * An;
			}
		}



	}

        //*****************************

        void AddNormalsAndShrink(ModelPart& rModelPart,const unsigned int &dimension, int MeshId=0)
	{

			
		//************ COMPUTE NEIGHBOUR BOUNDARY ELEMENTS

		ModelPart::NodesContainerType&     rNodes = rModelPart.Nodes(MeshId);
		ModelPart::ConditionsContainerType& rCond = rModelPart.Conditions(MeshId);

		std::vector<WeakPointerVector<Condition> > nConditions(rNodes.size()+1);
		std::vector<int> Ids (rModelPart.NumberOfNodes()+1) ;
		std::fill( Ids.begin(), Ids.end(), 0 );

	
		//*************  Neigbours of nodes search ************//
		//add the neighbour conditions to all the nodes in the mesh
		int id=1;
		for(ModelPart::ConditionsContainerType::iterator ic = rCond.begin(); ic!=rCond.end(); ic++)
		  {
		    if(ic->IsNot(CONTACT) && ic->Is(BOUNDARY)){
			
		      Condition::GeometryType& pGeom = ic->GetGeometry();

		      for(unsigned int i = 0; i < pGeom.size(); i++)
			{
			  if( mEchoLevel > 2 ){
			    std::cout<<" Condition ID "<<ic->Id()<<" id "<<id<<std::endl;
			    if(Ids.size()<=pGeom[i].Id())
			      std::cout<<" Shrink node in geom "<<pGeom[i].Id()<<" number of nodes "<<Ids.size()<<std::endl;
			  }
			  
			  if(Ids[pGeom[i].Id()]==0){
			    Ids[pGeom[i].Id()]=id;
 			    nConditions[id].push_back( Condition::WeakPointer( *(ic.base()) ) );
			    id++;
			  }
			  else{
			    
			    nConditions[Ids[pGeom[i].Id()]].push_back( Condition::WeakPointer( *(ic.base()) ) );
			  }
			    
			}
		      
		    }
		  }
			
		
		if( id > 1 ){

		  //**********Set Boundary Nodes Only

		  ModelPart::NodesContainerType::iterator nodes_begin = rNodes.begin();
		  ModelPart::NodesContainerType  BoundaryNodes;

		
		  for(unsigned int i = 0; i<rNodes.size(); i++)
		    {
		      if((nodes_begin + i)->Is(BOUNDARY)){
			BoundaryNodes.push_back( *((nodes_begin+i).base()) ); 
		      }
		    }



		  //********** Construct the normals in order to have a good direction and length for applying a contraction
			

		  std::vector<double> storenorm;
		  std::vector<double> tipnormal;
		
		  ModelPart::NodesContainerType::iterator boundary_nodes_begin = BoundaryNodes.begin();
		  int boundary_nodes = BoundaryNodes.size();
		
		  int not_assigned = 0;

		  int pn=0;
		  // #pragma omp parallel for private(pn,storenorm,tipnormal)
		  for (pn=0; pn<boundary_nodes; pn++)
		    { 

		      double cosmedio=0;
		      double totalnorm=0;
		      double directiontol=0;
		      double numnorm=0;
		      int indepnorm=0,acuteangle=0;
		    
		      array_1d<double,3>&  Normal=(boundary_nodes_begin + pn)->FastGetSolutionStepValue(NORMAL);
		      double&       shrink_factor=(boundary_nodes_begin + pn)->FastGetSolutionStepValue(SHRINK_FACTOR);
				
		      //initialize
		      Normal.clear();
		      shrink_factor=0;

		      unsigned int normals_size = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]].size();
				
		      if( mEchoLevel > 2 )
			std::cout<<" Id "<<Ids[(boundary_nodes_begin + pn)->Id()]<<" normals size "<<normals_size<<std::endl;

		      storenorm.resize(normals_size);
		      std::fill(storenorm.begin(), storenorm.end(), 0 );

		      tipnormal.resize(normals_size);
		      std::fill(tipnormal.begin(), tipnormal.end(), 0 );

	
		      //Check coincident faces-normals
		      for (unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neighbor faces
			{

			  const array_1d<double,3>& AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
   
			  if (storenorm[esnod]!=1){ //if is not marked as coincident
	      
			    if (esnod+1<normals_size){
		
			      double tooclose=0;
			      for (unsigned int esn=esnod+1;esn<normals_size;esn++)//loop over node neighbor faces
				{
		    
				  const array_1d<double,3>& NormalVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esn].GetValue(NORMAL); //conditions

				  directiontol=inner_prod(AuxVector,NormalVector);
									
				  //std::cout<<"  -- > direction tol"<<directiontol<<" AuxVector "<<AuxVector<<" NormalVector "<<NormalVector<<std::endl;

				  if (directiontol>0.995 && storenorm[esn]==0){//to not have coincident normals
				    storenorm[esn]=1;	
				  }
				  else{
				    if (directiontol>0.95 && directiontol<0.995 && storenorm[esn]==0){//to not have close normals
				      tooclose+=1;
				      storenorm[esn]  =2;
				      storenorm[esnod]=2;
				    }
				    else{
				      if(directiontol<-0.005){
					//std::cout<<" acute angle "<<directiontol<<std::endl;
					acuteangle      +=1;
					tipnormal[esn]  +=1;
					tipnormal[esnod]+=1;
				      }
				    }
				  }

				}//end for the esnod for
		    
			      for (unsigned int esn=esnod; esn<normals_size;esn++)//loop over node neighbour faces
				{
				  if(storenorm[esn]==2 && tooclose>0)
				    storenorm[esn]=(1.0/(tooclose+1));   
				}
		
		
			    }	      
		
			    if(storenorm[esnod]!=0){ //!=1
								
			      Normal+=AuxVector*storenorm[esnod]; //outer normal
			      numnorm  +=storenorm[esnod];

			    }
			    else{

			      Normal+=AuxVector; //outer normal
			      numnorm  +=1;

			    }

			    totalnorm+=1;
	      
			      
			  }
	    
			}
	

		      //std::cout<<" Normal "<<Normal<<" ID "<<(boundary_nodes_begin + pn)->Id()<<std::endl;

		      //Modify direction of tip normals in 3D  --------- start ----------

		      if(dimension==3 && numnorm>=3 && acuteangle){
	  
			//std::cout<<" pn "<<pn<<" normal number: "<<numnorm<<" acute angles: "<<acuteangle<<std::endl;
			//std::cout<<"storenormal"<<storenorm<<std::endl;
			//std::cout<<"tipnormal"<<tipnormal<<std::endl;


			std::vector< array_1d<double,3> > N(3);
			indepnorm=0;

			if(numnorm==3){ //Definite solution for 3 planes

			  for (unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neighbor faces
			    {
			      if(tipnormal[esnod]>=1 && storenorm[esnod]==0){ //tip node correction

				const array_1d<double,3>& AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
				N[indepnorm]=AuxVector;
				indepnorm+=1;
			      }
			    }
	        
			}
			else{ //I must select only 3 planes 

			  double maxprojection=0;
			  double projection=0;
			  std::vector<int> pronormal(normals_size); 
			  pronormal.clear();

			  //get the positive biggest projection between planes
			  for (unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neighbor faces
			    {
			      maxprojection=0;
		
			      if(tipnormal[esnod]>=1 && storenorm[esnod]!=1){ 

				const array_1d<double,3>& AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

				for (unsigned int esn=0;esn<normals_size;esn++)//loop over node neighbor faces to check the most coplanar
				  {
				    if(tipnormal[esn]>=1 && storenorm[esn]!=1 && esnod!=esn){ 


				      const array_1d<double,3>& NormalVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esn].GetValue(NORMAL); //conditions
				      projection=inner_prod(AuxVector,NormalVector);
			
				      if(maxprojection<projection && projection>0.5){
					maxprojection=projection;
					pronormal[esnod]=esn+1; //set in pronormal one position added to the real one the most coplanar planes
				      }
				    }
				  }

			      }
		
			    }
	    
			  // 	    std::cout<<" PROJECTIONS "<<pn<<" pronormal"<<pronormal<<std::endl;
		  
			  //get the most obtuse normals 
			  for (unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neighbor faces
			    {
		
			      if(indepnorm<3){
  
				if(tipnormal[esnod]>=1 && storenorm[esnod]!=1){ //tip node correction
		    
				  array_1d<double,3> AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
										
				  if(storenorm[esnod]!=0) //!=1
				    AuxVector*=storenorm[esnod];
		  
				  for (unsigned int esn=esnod+1;esn<normals_size;esn++)//loop over node neighbor faces to check the most coplanar
				    {
				      if(tipnormal[esn]>=1 && storenorm[esn]!=1){ //tip node correction

					if(esnod+1==(unsigned int)pronormal[esn]){

					  const array_1d<double,3>& NormalVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
				    
					  if(storenorm[esnod]!=0){ //!=1
					    AuxVector+=NormalVector*storenorm[esnod];
					  }
					  else{
					    AuxVector+=NormalVector;
					  }
					  AuxVector=AuxVector/norm_2(AuxVector);
					  //std::cout<<" plane "<<esnod<<" esn "<<esn<<std::endl;
					  storenorm[esn]=1; //to not use this plane again
					}
				      }

				    }


				  N[indepnorm]=AuxVector;
				  indepnorm+=1;

				}
		  
			      }
			      // 		else{

			      // 		  std::cout<<" pn "<<pn<<" pronormal"<<pronormal<<std::endl;
			      // 		}

		
			    }
	  
			}

			//std::cout<<" indepnorm "<<indepnorm<<std::endl;

			if(indepnorm<3){

			  for (unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neighbor faces
			    {
			      if(indepnorm==3)
				break;
		
			      if(tipnormal[esnod]>=1 && storenorm[esnod]!=0 && storenorm[esnod]!=1){ //tip node correction
				const array_1d<double,3>& AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
				N[indepnorm]=AuxVector;
				indepnorm+=1;
			      }

			      if(storenorm[esnod]==0 && tipnormal[esnod]==0 && indepnorm>0){ //if only two of the tip normals are acute one is not stored (corrected here 05/09/2011)
				const array_1d<double,3>& AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
				N[indepnorm]=AuxVector;
				indepnorm+=1;
			      }

			    }
			}



			if(indepnorm==3){
			  double aux=inner_prod(N[0],MathUtils<double>::CrossProduct(N[1],N[2]));
			  if(aux!=0){
			    Normal=((MathUtils<double>::CrossProduct(N[1],N[2])+MathUtils<double>::CrossProduct(N[2],N[0])+MathUtils<double>::CrossProduct(N[0],N[1]))/aux); //intersection of three planes
			  }
			  // else{
			  //   std::cout<<" aux "<<aux<<" Normals :";
			  // }
	    
			}
	  
  
			//Check Boundary
			if(norm_2(Normal)>2){ //a bad solution... but too ensure consistency
			  //Normal.write(" TO LARGE NORMAL ");
			  //std::cout<<" modulus 1 "<<Normal.modulus()<<std::endl;
			  Normal=Normal/norm_2(Normal);
			  Normal*=1.2;
			  //Normal.write(" CORRRECTION ");
			  //std::cout<<" modulus 2 "<<Normal.modulus()<<std::endl;
			}

			//Modify direction of tip normals in 3D ----------- end  -----------

		      }
		      else{

			//std::cout<<" Normal "<<Normal<<std::endl;

			if(norm_2(Normal)!=0)
			  Normal=Normal/norm_2(Normal); //normalize normal *this/modulus();
	  
			for(unsigned int esnod=0;esnod<normals_size;esnod++)//loop over node neigbour faces
			  {    
			    if (storenorm[esnod]!=1){
		
			      array_1d<double,3> AuxVector = nConditions[Ids[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
							
			      if(norm_2(AuxVector))
				AuxVector/=norm_2(AuxVector);
		
			      directiontol=inner_prod(AuxVector,Normal);
		
			      cosmedio+=directiontol;
		
			    }
			  }
	  
	  
			if (totalnorm!=0)
			  cosmedio*=(1.0/totalnorm);
	  
			//acute angles are problematic, reduction of the modulus in that cases
			if (cosmedio<=0.3){
			  cosmedio=0.8;
			  //std::cout<<" acute correction "<<std::endl;
			}
	  
			if(cosmedio>3)
			  std::cout<<" cosmedio "<<cosmedio<<std::endl;

			//if (cosmedio!=0) {
			if (cosmedio>0 && cosmedio>1e-3) { //to ensure consistency
			  cosmedio=1.0/cosmedio;
			  Normal*=cosmedio;               //only to put the correct length   
			  //(boundary_nodes_begin + pn)->SetValue(Normal);
			  //std::cout<<pn<<" cosmedio "<<cosmedio<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<std::endl;		     
			}
		      }
		    
		      //std::cout<<(boundary_nodes_begin + pn)->Id()<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<std::endl;		     

		      //Now Normalize Normal and store the Shrink_Factor
		      shrink_factor=norm_2(Normal);

		      if(shrink_factor!=0)
			{
			  if( mEchoLevel > 2 )
			    std::cout<<"Id "<<(boundary_nodes_begin + pn)->Id()<<" shrink_factor "<<shrink_factor<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<" cosmedio "<<cosmedio<<std::endl;
			  Normal/=shrink_factor;

			}
		      else{

			if( mEchoLevel > 1 )
			  std::cout<<"Id "<<(boundary_nodes_begin + pn)->Id()<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<" cosmedio "<<cosmedio<<std::endl;		     
			
			Normal.clear();
			shrink_factor=1;

			//std::cout<<" ERROR: normal shrinkage calculation failed "<<std::endl;
			not_assigned +=1;
		      }

		    
		    }

		
		  if(mEchoLevel > 0) 
		    std::cout<<"  [ Normals_Shrinkage (BoundaryNodes:"<<boundary_nodes<<")[SET:"<<boundary_nodes-not_assigned<<" / NOT_SET:"<<not_assigned<<"] ][Mesh:"<<MeshId<<"]"<<std::endl;


		}
	}


	/*@} */
	/**@name Private  Acces */
	/*@{ */


	/*@} */
	/**@name Private Inquiry */
	/*@{ */


	/*@} */
	/**@name Un accessible methods */
	/*@{ */

	/*@} */

}; /* Class BoundaryNormalsCalculationUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_FACTOR  defined */

