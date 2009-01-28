#include "custom_utilities/AdvancedNMPointsMapper.h"
#include "fsi_application.h"


namespace Kratos
{

	//in the constructor we count the number of conditions n_ic (is equal to the number of 
	//Gausspoints in our case), which do have nodes 
	//with IS_INTERFACE == 1.0
	
	AdvancedNMPointsMapper::AdvancedNMPointsMapper(ModelPart& origin_model_part, ModelPart& dest_model_part)
		:OriginModelPart(origin_model_part), DestinationModelPart(dest_model_part)
	{
		/*std::cout<<std::endl<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
		std::cout<<"AdvancedNMPointsMapper CONSTRUCTOR STARTED!!!!!!!!!!"<<std::endl;*/
		//cout number of Gausspoints
		this->n=0;
		for(ModelPart::ConditionsContainerType::iterator ic = DestinationModelPart.ConditionsBegin();
			ic!=DestinationModelPart.ConditionsEnd(); ic++)
		{
			if (((*ic).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*ic).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*ic).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   n+=3;
		}
		//std::cout<<"NUMBER OF GAUSSPOINTS IS "<<n<<std::endl;

		//create all the necessary arrays of the size n_ic

		//create an array for GaussPoints of the destination 
		DestGaussPoints = new double*[n];
		for (int i=0;i<n;i++)DestGaussPoints[i] = new double[3];	//first three entries - for x,y,z
		for (int i=0;i<n;i++)										//now we initialize the array with zeros
		{	
			for (int j=0; j<3;j++)DestGaussPoints[i][j]=0;
			
		}
		//Create an array for Aerea of Destination Elements
		DestArea = new double[n];
		for (int i=0; i<n; i++)DestArea[i]=0;
		
		//Create an arry for Normal of Destination Conditions
		DestCondsNormals = new double*[n];
		for (int i=0;i<n;i++)DestCondsNormals[i] = new double[3];	//first three entries - for x,y,z of the normal
		
		//here we create an arrays which will store pointers to the origin and destination conditions, 
		//which are neighbours of a i_th entry of the destination array
		OriginConds = new Condition*[n];
		DestinationConds = new Condition*[n];
		
		//this array is invented to check that a normal from the Gauss point of the 
		//destination conditions, crosses only one origin condition
		//when the entry of CHECK is more than one - it means an ambiguous case where 
		//distances have to be compared and only one has to be chosen for projection
		check = new int[n];
		for (int i=0; i<n; i++)check[i]=0;
		
		//this array is intended to store the coordinates of the origin points x, y, z, and 
		//the values of the shape functions at these points which will be needed for 
		//projection (xi, eta). everythings is initialized with 0
		OriginProjectionPoints = new double*[n];
		for(int i=0; i<n; i++)OriginProjectionPoints[i] = new double[5];
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<5; j++)OriginProjectionPoints[i][j] = 0;
		}
		
		//calculating the new normals
		NormalCalculationUtils normal_tools;
		normal_tools.CalculateOnSimplex(DestinationModelPart.Conditions(),3);
		
		//in this loop we calculate and set the Gausspoints, the Area and Normals of the corresponiding element and fill the
		//array of pointers to the corresponding conditions, which will be used to do further loops to save time
		double n_normed;
		int h=0;
		for(ModelPart::ConditionsContainerType::iterator ic = DestinationModelPart.ConditionsBegin();
				ic!=DestinationModelPart.ConditionsEnd(); ic++)	
		{
			if (((*ic).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*ic).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*ic).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )
			{	
				/////later we have to know with GP is near of a node to know
				DestGaussPoints[h][0] = 0.6*(ic->GetGeometry()[0].X())+0.2*(ic->GetGeometry()[1].X())+0.2*(ic->GetGeometry()[2].X() );
				DestGaussPoints[h][1] = 0.6*(ic->GetGeometry()[0].Y())+0.2*(ic->GetGeometry()[1].Y())+0.2*(ic->GetGeometry()[2].Y() );
				DestGaussPoints[h][2] = 0.6*(ic->GetGeometry()[0].Z())+0.2*(ic->GetGeometry()[1].Z())+0.2*(ic->GetGeometry()[2].Z() );
				
				DestGaussPoints[h+1][0] = 0.2*(ic->GetGeometry()[0].X())+0.6*(ic->GetGeometry()[1].X())+0.2*(ic->GetGeometry()[2].X() );
				DestGaussPoints[h+1][1] = 0.2*(ic->GetGeometry()[0].Y())+0.6*(ic->GetGeometry()[1].Y())+0.2*(ic->GetGeometry()[2].Y() );
				DestGaussPoints[h+1][2] = 0.2*(ic->GetGeometry()[0].Z())+0.6*(ic->GetGeometry()[1].Z())+0.2*(ic->GetGeometry()[2].Z() );
				
				DestGaussPoints[h+2][0] = 0.2*(ic->GetGeometry()[0].X())+0.2*(ic->GetGeometry()[1].X())+0.6*(ic->GetGeometry()[2].X() );
				DestGaussPoints[h+2][1] = 0.2*(ic->GetGeometry()[0].Y())+0.2*(ic->GetGeometry()[1].Y())+0.6*(ic->GetGeometry()[2].Y() );
				DestGaussPoints[h+2][2] = 0.2*(ic->GetGeometry()[0].Z())+0.2*(ic->GetGeometry()[1].Z())+0.6*(ic->GetGeometry()[2].Z() );
				
				DestArea[h] = GetArea(*ic); 
				DestArea[h+1] = DestArea[h]; 
				DestArea[h+2] = DestArea[h];

				const double &n_x = ic->GetValue(NORMAL)[0];
				const double &n_y = ic->GetValue(NORMAL)[1];
				const double &n_z = ic->GetValue(NORMAL)[2];
				n_normed = sqrt(pow(n_x,2)+pow(n_y,2)+pow(n_z,2));
				
				DestCondsNormals[h][0] = n_x/n_normed;
				DestCondsNormals[h][1] = n_y/n_normed;
				DestCondsNormals[h][2] = n_z/n_normed;
				
				//DestCondsNormals[h][0] = ic->GetValue(NORMAL)[0];
				//DestCondsNormals[h][1] = ic->GetValue(NORMAL)[1];
				//DestCondsNormals[h][2] = ic->GetValue(NORMAL)[2];

				DestCondsNormals[h+1][0] = DestCondsNormals[h][0];
				DestCondsNormals[h+1][1] = DestCondsNormals[h][1];
				DestCondsNormals[h+1][2] = DestCondsNormals[h][2];

				DestCondsNormals[h+2][0] = DestCondsNormals[h][0];
				DestCondsNormals[h+2][1] = DestCondsNormals[h][1];
				DestCondsNormals[h+2][2] = DestCondsNormals[h][2];

				//Condition& tempcond=(*ic);
				//DestinationConds[i]=&tempcond;
				DestinationConds[h]=&(*ic);
				DestinationConds[h+1]=DestinationConds[h];
				DestinationConds[h+2]=DestinationConds[h];
				
				//std::cout<<"GAUSSPOINT Nr " <<i<<" COORDINATES: "<<DestGaussPoints[i][0]<<","<<DestGaussPoints[i][1]<<","<<DestGaussPoints[i][2]<<","<<std::endl;
				//std::cout<<"GAUSSPOINT Nr " <<i<<" AREA: "<<DestArea[i]<<std::endl;
				//std::cout<<"GAUSSPOINT Nr " <<i<<" NORMAL: "<<DestCondsNormals[i][0]<<","<<DestCondsNormals[i][1]<<","<<DestCondsNormals[i][2]<<","<<std::endl;
				//std::cout<<std::endl;
				
				h+=3;
			}
		}
	//std::cout<<"CONTRUCTOR END"<<std::endl<<"+++++++++++++++++++++++++++++++++++++++++++++"<<std::endl<<std::endl;
	//DestGaussPoints[310][2]=0.3;//fuer plane_nm_str_unstr
	}
	
	AdvancedNMPointsMapper::~AdvancedNMPointsMapper()
	{
		// free memeory
		for (int i=0;i<n;i++) delete DestGaussPoints[i]; delete DestGaussPoints;		
		delete[] DestArea;
		delete[] check;
		for (int i=0;i<n;i++) delete OriginProjectionPoints[i]; delete OriginProjectionPoints;		
		for (int i=0;i<n;i++) delete DestCondsNormals[i]; delete DestCondsNormals;		
		/////for (int i=0;i<n;i++) delete OriginConds[i];	//NO. Dont delete Conditions (references) we would delete original
		delete OriginConds;								//just delete the pointer array
		/////for (int i=0;i<n;i++) delete DestinationConds[i]; //NO. Dont delete Conditions (references) we would delete original
		delete DestinationConds;
	}
	
	//later in template class 
	//void AdvancedNMPointsMapper<TVar>::Map(double SearchRadiusFactor, int max_iter, double tol_iter)
	//void AdvancedNMPointsMapper::Map(double SearchRadiusFactor, int max_iter, double tol_iter)
	//{
	//	setDestGaussPoints();
	//	FindNeighbours(SearchRadiusFactor);
	//	//calc_p_nodal_dest(max_iter, tol_iter);//replace with ScalarMap and VectorMap		
	//}
	

	//store Gausspoints in Datastructure
	//we use kd-tree of library ANN
	void AdvancedNMPointsMapper::FindNeighbours(double SearchRadiusFactor)
	{
		//boost::timer AnnTree_timer;			//timer for construction time of tree
		//first we create a KD-tree, according to the array of the Gauss points of the destination
		////////////////////Build Tree Object/////////////////////////////////////
		int				dim				= 3;					// dimension
		double			eps				= 0;					// error bound
		int				bucketsize		= 4;
		ANNkd_tree*			kdTree;								// search structure
		//we should work with references here!!!!!!!!!!!!!!!!
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ANNpointArray		dataPts;							// data points as ANN-Datatype
		dataPts = annAllocPts(n, dim);							//allocate memory for data points
		for(int i=0;i<n;i++){dataPts[i]=DestGaussPoints[i];}	//Data Points = Gauss Points
		//cout<<"DataPoints"<<endl;
		
		kdTree = new ANNkd_tree(					// build search structure
						dataPts,					// the data points
						n,							// number of points
						dim,						// dimension of space
						bucketsize);

		/////////////////Search Nearest Neigbour in Tree Object ///////////////////////////
		//std::cout << "  ANN Tree constructed in TIME = " << AnnTree_timer.elapsed() << std::endl;	//timer for construction time of tree
		/*std::cout<<"TREE FOR GAUSPOINTS CONSTRUCTED"<<std::endl;
		std::cout<<std::endl;*/
				
		//int j=0;
		
		//std::cout<<"Start Searching Loop Over OriginElements and Check and Set all Neighbours "<<std::endl;
		//double SetPP=0.0;					//time for the function setOriginProjectionPoints
		//double AnnSearch=0.0;				//time for ANN Search
		//boost::timer ss_timer;			//timer for ANN Search And function setOriginProjectionPoints
		//boost::timer AnnSearch_timer;		//timer for pure ANN Search time
		
		for(ModelPart::ConditionsContainerType::iterator icc = OriginModelPart.ConditionsBegin();
				icc!=OriginModelPart.ConditionsEnd(); icc++)	
		{
			//AnnSearch_timer.restart();
			if (((*icc).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*icc).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
				((*icc).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )
			{
				//std::cout<<std::endl<<std::endl;
				//std::cout<<"_______________________________________________"<<std::endl;
				//std::cout<<"ORIGIN ELEMENT  Nr "<<j<<" ORIGIN ELEMENT  Nr "<<j<<" ORIGIN ELEMENT  Nr "<<j<<std::endl;j++;
				//std::cout<<"_______________________________________________"<<std::endl;
				//std::cout<<std::endl;
				
				double midpoint[3]; 
				GetMidPoint(*icc, midpoint);
				
				double sradius = pow(SearchRadiusFactor,2)*GetSquaredRadius(*icc, midpoint);
				
				ANNpoint			queryPt=midpoint;				// query point
				ANNidxArray			nnIdx;					// near neighbor indices
				ANNdistArray		dists;					// near neighbor distances
				
				//cout<<"queryPoint"<<endl;
				/*std::cout<<" MIDPOINT IS: "<<midpoint[0]<<","<<midpoint[1]<<","<<midpoint[2]<<std::endl;
				std::cout<<" SQUARED RADIUS IS: "<<sradius<<std::endl;*/
				
				//first count number of neighbour k in radius
				int k = kdTree->annkFRSearch(
							queryPt,
							sradius,		//SQUARED DISTANCE!!!!!!!!!!!
							0,
							NULL,
							NULL,
							eps);
				//cout<<k<<" nearest neighbours"<<endl;

				nnIdx = new ANNidx[k];						// allocate near neigh indices
				dists = new ANNdist[k];						// allocate near neighbor dists	
				
				//then search these k neighbours
				kdTree->annkFRSearch(
							queryPt,
							sradius,		//SQUARED DISTANCE!!!!!!!!!!!
							k,
							nnIdx,
							dists,
							eps);

				//std::cout<<" THE NEAREST NEIGHBOURS ARE:"<<std::endl;
				//std::cout << " \tNN:\tIndex\tDistance\n";
				//for (int i = 0; i < k; i++) {			// print summary
				//		dists[i] = sqrt(dists[i]);			// unsquare distance
				//		cout << " \t" << i << "\t" << nnIdx[i] << "\t" << dists[i] << "\n";
				//}
				//AnnSearch+=AnnSearch_timer.elapsed();

				//boost::timer setpp_timer; //timer for setOriginProjectionPoints
				//setpp_timer.restart();
				for(int i = 0; i < k; i++)
				{
					//std::cout<<std::endl<<"++++++++ NEIGBOUR Nr "<<i<<" ++++++++ "<<std::endl;
					setOriginProjectionPoints(nnIdx[i], *icc);
				}
				//SetPP+=setpp_timer.elapsed();

				delete [] nnIdx;							// clean things up
				delete [] dists;
			}
		}
	//std::cout << "    pure ANN Search TIME = " << AnnSearch << std::endl;
	//std::cout << "    Set ProjectionPoints TIME = " << SetPP << std::endl;
	//std::cout<<"Searching Loop Over OriginElements Finished "<<std::endl<<std::endl;
	//std::cout << "  ANN Search and SetPP finished in TIME = " << ss_timer.elapsed() << std::endl;	
	delete kdTree;
	annClose();

	checkBadGP();
	}

	
	void AdvancedNMPointsMapper::setOriginProjectionPoints(int i,  Condition& OriginNeighbourConds)
	{

		array_1d<double,3> Proj_local_coords;
		calcLocalPPCoords(Proj_local_coords, i, OriginNeighbourConds);

		if(check[i]==0)
		{
			if ( Proj_local_coords[0]>=0 && Proj_local_coords[1]>=0 &&  (1.0-Proj_local_coords[0]-Proj_local_coords[1])>=0)
			{
				//Global coordinates??????JA UM DISTANCE ZU VGL
				//setPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
				setlocalPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
				check[i]=1;
			}
			else 
			{
				//setlocalPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
				//std::cout<<"ok, check jetzt 2"<<std::endl;
				check[i]=2;
			}
		}
		else if(check[i]==2)
		{
			if ( Proj_local_coords[0]>=0 && Proj_local_coords[1]>=0 &&  (1.0-Proj_local_coords[0]-Proj_local_coords[1])>=0)
			{
				//Global coordinates??????JA UM DISTANCE ZU VGL
				//setPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
				setlocalPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
				check[i]=1;
			}
		}
		else if(check[i]==1)
		{
			if ( Proj_local_coords[0]>=0 && Proj_local_coords[1]>=0 &&  (1.0-Proj_local_coords[0]-Proj_local_coords[1])>=0)
			{
				//std::cout<<"  we have to compare distances, because there is already a Projection Point for this Gauspoint! "<<std::endl;
				//calculate global coordinates, because we have to cslculate distances, in the corresponding condition!!
				setglobalPPCoordinates(OriginProjectionPoints[i][3],OriginProjectionPoints[i][4],i,*OriginConds[i]);
				
				//std::cout<<"global coordinates GP:  "<<DestGaussPoints[i][0]<<"  " <<DestGaussPoints[i][1]<<"  " <<DestGaussPoints[i][2] <<std::endl;
				//std::cout<<"global coordinates PP:  "<<OriginProjectionPoints[i][0]<<"  " <<OriginProjectionPoints[i][1]<<"  " <<OriginProjectionPoints[i][2] <<std::endl;
					
				//calculate distance old pp and gausspoint
				double old_sdistance = GetSquaredDistance(OriginProjectionPoints[i],DestGaussPoints[i]); //Cuidado! the real arrys are longer, but we use the first three entries here
				//std::cout<<"GP index: "<<i<<std::endl;
				//std::cout<<"old dist: "<<sqrt(old_sdistance)<<std::endl;
				//std::cout<<"current dist: "<<Proj_local_coords[2]<<std::endl<<std::endl;
				if(pow(Proj_local_coords[2],2)<old_sdistance)
				{
					//std::cout<<"replace "<<std::endl;
					setlocalPPCoordinates(Proj_local_coords, i, OriginNeighbourConds);
					//set global. we do not need it but we dont want to mix new with old coordinates
					setglobalPPCoordinates(OriginProjectionPoints[i][3],OriginProjectionPoints[i][4],i,OriginNeighbourConds);
					
					//std::cout<<"global coordinates PP new:  "<<OriginProjectionPoints[i][0]<<"  " <<OriginProjectionPoints[i][1]<<"  " <<OriginProjectionPoints[i][2] <<std::endl;
				
					//check[i]=1;
				}
			}
		}
	
		////we construct a normal  n from the baricenter G of the destination element
		////and find the value of the variable of interest (e,g, pressure) at the element of the 
		////origin, at the point, where the normal crosses the origin element
		//
		////global coordinates of the baricneter of the destination element
		//
		//const double &dest_x = DestGaussPoints[i][0];
		//const double &dest_y = DestGaussPoints[i][1];
		//const double &dest_z = DestGaussPoints[i][2];
		////now we want to solve the system of equations of the type	
		////dest_x = xi*x_a_origin + eta*x_b_origin + (1-xi-eta)*x_c_origin - gamma*n_x
		////dest_y = xi*y_a_origin + eta*y_b_origin + (1-xi-eta)*y_c_origin - gamma*n_y
		////and same statement for z coordinate
		////where gamma is the uknown distance between facets and n_x... n_z  are the 
		////components of the normal to the destination element, which is known for every element
		//
		////this system can be simply modified to the one below simply by isolating 
		////xi, eta and gamma
		////
		//// Note!!! that we know the coordinates of the destination point!!
		//// [dest_x-c_x]		| a_x-c_x   b_x-c_x   -n_x |  |xi		|
		//// [dest_y-c_y] =	| a_y-c_y   b_x-c_y   -n_y |* |eta		|
		//// [dest_z-c_z]		| a_z-c_z   b_x-c_z   -n_z |  |1-eta-xi |
		//const double &a_x = OriginNeighbourConds.GetGeometry()[0].X();
		//const double &b_x = OriginNeighbourConds.GetGeometry()[1].X();
		//const double &c_x = OriginNeighbourConds.GetGeometry()[2].X();

		//const double &a_y = OriginNeighbourConds.GetGeometry()[0].Y();
		//const double &b_y = OriginNeighbourConds.GetGeometry()[1].Y();
		//const double &c_y = OriginNeighbourConds.GetGeometry()[2].Y();
		//
		//const double &a_z = OriginNeighbourConds.GetGeometry()[0].Z();
		//const double &b_z = OriginNeighbourConds.GetGeometry()[1].Z();
		//const double &c_z = OriginNeighbourConds.GetGeometry()[2].Z();
		////
		////normal n to ABC
		//
		////const array_1d<double,3>& normal =DestConds[i].GetValue(NORMAL);
		//const double &n_x = DestCondsNormals[i][0];
		//const double &n_y = DestCondsNormals[i][1];
		//const double &n_z = DestCondsNormals[i][2];
		//	//normal[2];
		//
		////matrix of the known values, that needs to be inverted to solve the linear system
		////in order to obtain the local coordinates xi, eta of the projection p
		//Matrix Trafo(3,3);
		//Matrix InvTrafo(3,3);
		//double det;
		//
		//Trafo(0,0) = a_x-c_x;	Trafo(0,1) = b_x-c_x;	Trafo(0,2) = -n_x;
		//Trafo(1,0) = a_y-c_y;	Trafo(1,1) = b_y-c_y;	Trafo(1,2) = -n_y;
		//Trafo(2,0) = a_z-c_z;	Trafo(2,1) = b_z-c_z;	Trafo(2,2) = -n_z;

		////and the righ-hand side
		//
		////create outside or in constructor for timsteps and large displacements, because
		////then we have to execute this function each timestep
		////(outside!!!performance!!!pass array in functions to use it)
		//array_1d<double,3> RHS;
		//RHS[0]=dest_x-c_x;	RHS[1]=dest_y-c_y;	RHS[2]=dest_z-c_z;
		//array_1d<double,3> Proj_local_coords;
		////std::cout<<" EQUATION FOR PROJECTION: "<<std::endl;
		////std::cout<<"   TRAFO MATRIX: ";
		////std::cout<<Trafo<<std::endl;
		//
		//MathUtils<double>::InvertMatrix3(Trafo, InvTrafo, det);	
		////std::cout<<"   Inv TrAFO ";
		////std::cout<<InvTrafo<<std::endl;
		////std::cout<<"   RHS: "<<RHS<<std::endl;
		//noalias(Proj_local_coords) = prod(InvTrafo, RHS);
		////std::cout<<"   LOCAL COORDINATES OF PROJECTION POINT and DISTANCE: "<<Proj_local_coords<<std::endl;
		////now we have the local coordinates of the point on the origin mesh
		////which resulted from the orthogonal projection of point dest, which
		////was the baricenter of the destination element
		////std::cout<<"  Index of Gausspoint: "<<i<<std::endl;
		//if ( Proj_local_coords[0]>=0 && Proj_local_coords[1]>=0 &&  (1.0-Proj_local_coords[0]-Proj_local_coords[1])>=0)
		//{
		//	//std::cout<<" ---> this point lies inside the element !!! so set Projection Point! "<<std::endl;
		//	//we could do this in a better way
		//	if(check[i]==0)
		//	{
		//		//and shape functions values (first two are shape functions, and third is the distance)
		//		OriginProjectionPoints[i][3] = Proj_local_coords[0];
		//		OriginProjectionPoints[i][4] = Proj_local_coords[1];
		//		//OriginProjectionPoints[i][5] = Proj_local_coords[2];//THIS IS THE DISTANCE DISTANCE DISTANCE DISTANCE
		//		//now we substitute the obtained values of shaoefunctions and gamma into original system
		//		//to obtain the global coordinates of the projection point
		//		//we store the x,y,z coordinates of the origin-conds?
		//		OriginProjectionPoints[i][0] = a_x*Proj_local_coords[0]+b_x*Proj_local_coords[1]+c_x*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//									+ Proj_local_coords[2]*n_x;
		//		OriginProjectionPoints[i][1] = a_y*Proj_local_coords[0]+b_y*Proj_local_coords[1]+c_y*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//									+ Proj_local_coords[2]*n_y;
		//		OriginProjectionPoints[i][2] = a_z*Proj_local_coords[0]+b_z*Proj_local_coords[1]+c_z*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//									+ Proj_local_coords[2]*n_z;
		//		//and the respective condition is saved in the array
		//		//std::cout<<" OriginProjectionPoints: 3 global and 2 local: ";
		//		//std::cout<<OriginProjectionPoints[i][0]<<","<<OriginProjectionPoints[i][1]<<","<<OriginProjectionPoints[i][2]<<","<<OriginProjectionPoints[i][3]<<","<<OriginProjectionPoints[i][4]<<std::endl;
		//		OriginConds[i] = &OriginNeighbourConds;
		//		check[i]=1;
		//	}
		//	else
		//	{
		//		//std::cout<<"  we have to compare distances, because there is already a Projection Point for this Gauspoint! "<<std::endl;
		//		double old_sdistance = GetSquaredDistance(OriginProjectionPoints[i],DestGaussPoints[i]); //Cuidado! the real arrys are longer, but we use the first three entries here
		//		if(pow(Proj_local_coords[2],2)<old_sdistance)
		//		{
		//			//std::cout<<"  we replace this Projection Point because it is closer! "<<std::endl;
		//			//do the same like if check would be zero and set check to 1
		//			OriginProjectionPoints[i][3] = Proj_local_coords[0];
		//			OriginProjectionPoints[i][4] = Proj_local_coords[1];
		//			//OriginProjectionPoints[i][5] = Proj_local_coords[2];//THIS IS THE DISTANCE DISTANCE DISTANCE DISTANCE
		//			OriginProjectionPoints[i][0] = a_x*Proj_local_coords[0]+b_x*Proj_local_coords[1]+c_x*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//										+ Proj_local_coords[2]*n_x;
		//			OriginProjectionPoints[i][1] = a_y*Proj_local_coords[0]+b_y*Proj_local_coords[1]+c_y*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//										+ Proj_local_coords[2]*n_y;
		//			OriginProjectionPoints[i][2] = a_z*Proj_local_coords[0]+b_z*Proj_local_coords[1]+c_z*(1.0-Proj_local_coords[0]-Proj_local_coords[1])
		//										+ Proj_local_coords[2]*n_z;
		//			OriginConds[i] = &OriginNeighbourConds;
		//			check[i]=1;
		//		}
		//	}
		//	
		//}
		///*else
		//{
		//std::cout<<" ---> this point lies NOT inside the element !!! "<<std::endl<<std::endl;
		//}*/
	}

	void AdvancedNMPointsMapper::checkBadGP()
	{
		int badPP = 0;
		for(int i=0;i<n;i++) {if(check[i]==2) badPP+=1;}
		std::cout<<"Number of not projected Gauss points "<<badPP<<". (from "<<n<<" Gauss Points) "<<std::endl;
		////////////////////
		double a_node[3];
		double b_node[3];
		double c_node[3];
		double a_dist;
		double b_dist;
		double c_dist;
		double sd;		//(muss grosser als kleinster abstand sein, dh. bei gefunden gp mit check=2)
						//(d. grösser als grösster suchradius (aber midpoint) DH: WIR nehmen grössten)
						//(durchmesser, woher weiss ich den??)
		
		
		//for all GaussPoints with check=2
		for(int i=0;i<n;i++)
		{
			if(check[i]==2)
			{
				//to obtain a start value for sd we assign values before looping
				ModelPart::ConditionsContainerType::iterator condit = OriginModelPart.ConditionsBegin();
				a_node[0] = condit->GetGeometry()[0].X();
				b_node[0] = condit->GetGeometry()[1].X();
				c_node[0] = condit->GetGeometry()[2].X();

				a_node[1] = condit->GetGeometry()[0].Y();
				b_node[1] = condit->GetGeometry()[1].Y();
				c_node[1] = condit->GetGeometry()[2].Y();
				
				a_node[2] = condit->GetGeometry()[0].Z();
				b_node[2] = condit->GetGeometry()[1].Z();
				c_node[2] = condit->GetGeometry()[2].Z();

				a_dist = GetSquaredDistance(a_node, DestGaussPoints[i]);
				b_dist = GetSquaredDistance(b_node, DestGaussPoints[i]);
				c_dist = GetSquaredDistance(c_node, DestGaussPoints[i]);
					
						sd=a_dist;
						setlocalPPCoordinates(1.0, 0.0, i, *condit); //arguments: xi, eta, index, condition
						//std::cout<<"node A "<<a_node<<std::endl;
						//std::cout<<OriginProjectionPoints[i][3]<<OriginProjectionPoints[i][4]<<DestGaussPoints[i] <<std::endl;
						//cout
				if(b_dist<sd)
					{
						sd=b_dist;
						setlocalPPCoordinates(0.0, 1.0, i, *condit);
						//std::cout<<"node B "<<a_node<<std::endl;
						//std::cout<<OriginProjectionPoints[i][3]<<OriginProjectionPoints[i][4] <<std::endl;
						
					}
				if(c_dist<sd)
					{
						sd=c_dist;
						setlocalPPCoordinates(0.0, 0.0, i, *condit);
						//std::cout<<"node C "<<a_node<<std::endl;
						//std::cout<<OriginProjectionPoints[i][3]<<OriginProjectionPoints[i][4] <<std::endl;
						
					}

				for(ModelPart::ConditionsContainerType::iterator icc = OriginModelPart.ConditionsBegin();
				icc!=OriginModelPart.ConditionsEnd(); icc++)	
				{
					if (((*icc).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
						((*icc).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
						((*icc).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )
					{
						//	setOriginPPforShortestDist((*icc), i, a_node, b_node, c_node, a_dist, b_dist, c_dist, sd);
						a_node[0] = icc->GetGeometry()[0].X();
						b_node[0] = icc->GetGeometry()[1].X();
						c_node[0] = icc->GetGeometry()[2].X();

						a_node[1] = icc->GetGeometry()[0].Y();
						b_node[1] = icc->GetGeometry()[1].Y();
						c_node[1] = icc->GetGeometry()[2].Y();
						
						a_node[2] = icc->GetGeometry()[0].Z();
						b_node[2] = icc->GetGeometry()[1].Z();
						c_node[2] = icc->GetGeometry()[2].Z();

						a_dist = GetSquaredDistance(a_node, DestGaussPoints[i]);
						b_dist = GetSquaredDistance(b_node, DestGaussPoints[i]);
						c_dist = GetSquaredDistance(c_node, DestGaussPoints[i]);
							
						if(a_dist<sd)
							{
								sd=a_dist;
								setlocalPPCoordinates(1.0, 0.0, i, *icc);
							}
						if(b_dist<sd)
							{
								sd=b_dist;
								setlocalPPCoordinates(0.0, 1.0, i, *icc);
							}
						if(c_dist<sd)
							{
								sd=c_dist;
								setlocalPPCoordinates(0.0, 0.0, i, *icc);
							}	
					}
				}
			}
		}
	}

	
	void AdvancedNMPointsMapper::calcLocalPPCoords(array_1d<double,3>& Proj_local_coords, int i, Condition& OriginNeighbourConds)
	{
		//we construct a normal  n from the baricenter G of the destination element
		//and find the value of the variable of interest (e,g, pressure) at the element of the 
		//origin, at the point, where the normal crosses the origin element
		
		//global coordinates of the baricneter of the destination element
		
		const double &dest_x = DestGaussPoints[i][0];
		const double &dest_y = DestGaussPoints[i][1];
		const double &dest_z = DestGaussPoints[i][2];
		//now we want to solve the system of equations of the type	
		//dest_x = xi*x_a_origin + eta*x_b_origin + (1-xi-eta)*x_c_origin - gamma*n_x
		//dest_y = xi*y_a_origin + eta*y_b_origin + (1-xi-eta)*y_c_origin - gamma*n_y
		//and same statement for z coordinate
		//where gamma is the uknown distance between facets and n_x... n_z  are the 
		//components of the normal to the destination element, which is known for every element
		
		//this system can be simply modified to the one below simply by isolating 
		//xi, eta and gamma
		//
		// Note!!! that we know the coordinates of the destination point!!
		// [dest_x-c_x]		| a_x-c_x   b_x-c_x   -n_x |  |xi		|
		// [dest_y-c_y] =	| a_y-c_y   b_x-c_y   -n_y |* |eta		|
		// [dest_z-c_z]		| a_z-c_z   b_x-c_z   -n_z |  |1-eta-xi |
		const double &a_x = OriginNeighbourConds.GetGeometry()[0].X();
		const double &b_x = OriginNeighbourConds.GetGeometry()[1].X();
		const double &c_x = OriginNeighbourConds.GetGeometry()[2].X();

		const double &a_y = OriginNeighbourConds.GetGeometry()[0].Y();
		const double &b_y = OriginNeighbourConds.GetGeometry()[1].Y();
		const double &c_y = OriginNeighbourConds.GetGeometry()[2].Y();
		
		const double &a_z = OriginNeighbourConds.GetGeometry()[0].Z();
		const double &b_z = OriginNeighbourConds.GetGeometry()[1].Z();
		const double &c_z = OriginNeighbourConds.GetGeometry()[2].Z();
		//
		//normal n to ABC
		
		//const array_1d<double,3>& normal =DestConds[i].GetValue(NORMAL);
		const double &n_x = DestCondsNormals[i][0];
		const double &n_y = DestCondsNormals[i][1];
		const double &n_z = DestCondsNormals[i][2];

		
		//matrix of the known values, that needs to be inverted to solve the linear system
		//in order to obtain the local coordinates xi, eta of the projection p
		Matrix Trafo(3,3);
		Matrix InvTrafo(3,3);
		double det;
		
		Trafo(0,0) = a_x-c_x;	Trafo(0,1) = b_x-c_x;	Trafo(0,2) = -n_x;
		Trafo(1,0) = a_y-c_y;	Trafo(1,1) = b_y-c_y;	Trafo(1,2) = -n_y;
		Trafo(2,0) = a_z-c_z;	Trafo(2,1) = b_z-c_z;	Trafo(2,2) = -n_z;

		//and the righ-hand side
		
		//create outside or in constructor for timsteps and large displacements, because
		//then we have to execute this function each timestep
		array_1d<double,3> RHS;
		RHS[0]=dest_x-c_x;	RHS[1]=dest_y-c_y;	RHS[2]=dest_z-c_z;
	//	array_1d<double,3> Proj_local_coords;
		//std::cout<<" EQUATION FOR PROJECTION: "<<std::endl;
		//std::cout<<"   TRAFO MATRIX: ";
		//std::cout<<Trafo<<std::endl;
		
		MathUtils<double>::InvertMatrix3(Trafo, InvTrafo, det);	
		//std::cout<<"   Inv TrAFO ";
		//std::cout<<InvTrafo<<std::endl;
		//std::cout<<"   RHS: "<<RHS<<std::endl;
		noalias(Proj_local_coords) = prod(InvTrafo, RHS);
		//shape functions values (first two are shape functions, and third is the distance)
		//std::cout<<"   LOCAL COORDINATES OF PROJECTION POINT and DISTANCE: "<<Proj_local_coords<<std::endl;
		//now we have the local coordinates of the point on the origin mesh
		//which resulted from the orthogonal projection of point dest, which
		//was the baricenter of the destination element
	}

	void AdvancedNMPointsMapper::setlocalPPCoordinates(array_1d<double,3> Proj_local_coords, int i,  Condition& OriginNeighbourConds)
	{
				//local coordinates
				OriginProjectionPoints[i][3] = Proj_local_coords[0];
				OriginProjectionPoints[i][4] = Proj_local_coords[1];
				//OriginProjectionPoints[i][5] = Proj_local_coords[2];//THIS IS THE DISTANCE DISTANCE DISTANCE DISTANCE
				OriginConds[i] = &OriginNeighbourConds;
	}
	
	void AdvancedNMPointsMapper::setlocalPPCoordinates(double xi, double eta, int i,  Condition& OriginNeighbourConds)
	{
				//local coordinates
				OriginProjectionPoints[i][3] = xi;
				OriginProjectionPoints[i][4] = eta;
				//OriginProjectionPoints[i][5] = Proj_local_coords[2];//THIS IS THE DISTANCE DISTANCE DISTANCE DISTANCE
				OriginConds[i] = &OriginNeighbourConds;
	}
	
	void AdvancedNMPointsMapper::setglobalPPCoordinates(double xi, double eta, int i,  Condition& OriginNeighbourConds)
	{
				const double &a_x = OriginNeighbourConds.GetGeometry()[0].X();
				const double &b_x = OriginNeighbourConds.GetGeometry()[1].X();
				const double &c_x = OriginNeighbourConds.GetGeometry()[2].X();

				const double &a_y = OriginNeighbourConds.GetGeometry()[0].Y();
				const double &b_y = OriginNeighbourConds.GetGeometry()[1].Y();
				const double &c_y = OriginNeighbourConds.GetGeometry()[2].Y();
				
				const double &a_z = OriginNeighbourConds.GetGeometry()[0].Z();
				const double &b_z = OriginNeighbourConds.GetGeometry()[1].Z();
				const double &c_z = OriginNeighbourConds.GetGeometry()[2].Z();

//				const double &n_x = DestCondsNormals[i][0];
//				const double &n_y = DestCondsNormals[i][1];
//				const double &n_z = DestCondsNormals[i][2];
				
				OriginProjectionPoints[i][0] = a_x*xi+b_x*eta+c_x*(1.0-xi-eta);
											
				OriginProjectionPoints[i][1] = a_y*xi+b_y*eta+c_y*(1.0-xi-eta);
											
				OriginProjectionPoints[i][2] = a_z*xi+b_z*eta+c_z*(1.0-xi-eta);
	}
	
	

	//area of a condition
	double AdvancedNMPointsMapper::GetArea(Condition& cond)
	{
	const double& x1 = cond.GetGeometry()[0].X(); const double& x2 = cond.GetGeometry()[1].X(); const double& x3 = cond.GetGeometry()[2].X();
	const double& y1 = cond.GetGeometry()[0].Y(); const double& y2 = cond.GetGeometry()[1].Y(); const double& y3 = cond.GetGeometry()[2].Y();
	const double& z1 = cond.GetGeometry()[0].Z(); const double& z2 = cond.GetGeometry()[1].Z(); const double& z3 = cond.GetGeometry()[2].Z();
	//lengthes of the edges of triangle
	//const double& length1 = sqrt( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2) ); const double& length2 = sqrt( pow(x3-x2,2)+pow(y3-y2,2)+pow(z3-z2,2)   ); const double& length3 = sqrt( pow(x1-x3,2)+pow(y1-y3,2)  +pow(z2-z1,2) );
	///const double&	s = 1.0/2.0*(length1+length2+length3);
	//area of an element acc. to Heron
	//return ( sqrt(s*(s-length1)*(s-length2)*(s-length3)) );
	const double l1=x2-x1; const double l2=y2-y1; const double l3=z2-z1; 
	const double r1=x2-x3; const double r2=y2-y3; const double r3=z2-z3; 
	const double k1=l2*r3-r2*l3; const double k2=l3*r1-l1*r3; const double k3=l1*r2-r1*l2; 
	return 0.5*sqrt(pow(k1,2)+pow(k2,2)+pow(k3,2));
	}
	
	//mid-point of the condition //this is not really the midpoint but the gausspoint
	void AdvancedNMPointsMapper::GetMidPoint(Condition& cond, double* midp)
	{
		//actually this is the Gausspoint, NOT the midpoint
		double x = 0.3333333*(cond.GetGeometry()[0].X()+cond.GetGeometry()[1].X()+cond.GetGeometry()[2].X());
		double y = 0.3333333*(cond.GetGeometry()[0].Y()+cond.GetGeometry()[1].Y()+cond.GetGeometry()[2].Y());
		double z = 0.3333333*(cond.GetGeometry()[0].Z()+cond.GetGeometry()[1].Z()+cond.GetGeometry()[2].Z());
		double midpoint[3]={x,y,z};
		//std::cout<<" MIDPPP"<<x<<y<<z<<std::endl;
		midp[0]=midpoint[0];	
		midp[1]=midpoint[1];
		midp[2]=midpoint[2];

	}
	
	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double AdvancedNMPointsMapper::GetSquaredRadius(Condition& cond, double* point)
	{
		double nodeA[3]={cond.GetGeometry()[0].X(),cond.GetGeometry()[0].Y(),cond.GetGeometry()[0].Z()};
		double nodeB[3]={cond.GetGeometry()[1].X(),cond.GetGeometry()[1].Y(),cond.GetGeometry()[1].Z()};
		double nodeC[3]={cond.GetGeometry()[2].X(),cond.GetGeometry()[2].Y(),cond.GetGeometry()[2].Z()};
		double dist1=GetSquaredDistance(nodeA,point);
		double dist2=GetSquaredDistance(nodeB,point);
		double dist3=GetSquaredDistance(nodeC,point);
		return std::max(dist1,std::max(dist2,dist3));
	}

	double AdvancedNMPointsMapper::GetSquaredDistance(double* point1, double* point2)
	{
		return pow(point1[0]-point2[0],2.0)+pow(point1[1]-point2[1],2.0)+pow(point1[2]-point2[2],2.0);
	}


	////this function will work with vector valued functions, such as e.g. velocities, displ.
	///*template <class Vector>
	//void AdvancedNMPointsMapper<Vector>::calc_p_nodal_dest(Condition& my_condition, Vector& vec_dest  )
	//{
	//}*/
	//
	////this function will work with scalar valued functions, such as e.g. pressures
	////typedef double TScalar;
	////template <class Tscalar>
	////void AdvancedNMPointsMapper<Tscalar>::calc_p_nodal_dest(Condition& my_condition, TScalar& p_dest  )
	
	
	
	
//	void AdvancedNMPointsMapper::calc_p_nodal_dest(int max_iter, double tol_iter)
//	{
//		boost::timer Mapini_timer;
//		OriginScalarValue = new double[n];
//		//std::cout<<"START MAPPING"<<std::endl;
//		//std::cout<<"Assign Nodal Area for DestModelPart"<<std::endl;
//		//initialize nodal area with 0
//		for(ModelPart::NodeIterator i = DestinationModelPart.NodesBegin() ; 
//		i != DestinationModelPart.NodesEnd() ; ++i)
//			{
//					//std::cout<<" HELLO1"<<std::endl;
//				if ((i->FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//				(i->FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//				(i->FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   
//				{
//					i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
//				}
//			}
//		////////////////////////////////////////////////////////////////
//		//std::cout<<"Calculate Nodal Area for DestModelPart"<<std::endl;
//		int i=0;
//		for(ModelPart::ConditionsContainerType::iterator ic = DestinationModelPart.ConditionsBegin();
//			ic!=DestinationModelPart.ConditionsEnd(); ic++)
//		{
//				//std::cout<<" HELLO2!"<<std::endl;
//			if (((*ic).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//				((*ic).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//				((*ic).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   
//			{
//				
//					//std::cout<<" HELLO3"<<std::endl;
//					//AUX needs to be defined as KRATOS variable
//					(*ic).GetGeometry()[0].GetValue(AUX)=0.0;
//					(*ic).GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*GetArea(*ic);
//					////we calculate pressure to be added on pressure of destination domain
//					(*ic).GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)=0.0;
//
//					(*ic).GetGeometry()[1].GetValue(AUX)=0.0;
//					(*ic).GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*GetArea(*ic);
//					(*ic).GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)=0.0;
//
//					(*ic).GetGeometry()[2].GetValue(AUX)=0.0;
//					(*ic).GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*GetArea(*ic);
//					(*ic).GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)=0.0;
//
//					OriginScalarValue[i] =	OriginConds[i]->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)*OriginProjectionPoints[i][3]
//									+ OriginConds[i]->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)*OriginProjectionPoints[i][4]
//									+ OriginConds[i]->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)*(1.0-OriginProjectionPoints[i][3]-OriginProjectionPoints[i][4]);
//					i++;
//			}
//		}
//		array_1d<double,3> loc_RHS;
//		array_1d<double,3> x_vect;
//		Matrix M_cons(3,3);
//		std::cout << "  Map Initialization TIME = " << Mapini_timer.elapsed() << std::endl;
//		///////////////////////////////////////////////////////////////////////////////////////
//		boost::timer Mapiter_timer;
//		//std::cout<<"Start Iterations"<<std::endl;
//		for(int k=0;k<max_iter;k++) //Iteration
//		{
//		///////////////////////////////////////////////////////////////////////////////////////
//			//std::cout<<"Iteration Nr. "<<k+1<<std::endl;
//			for(ModelPart::ConditionsContainerType::iterator ic = DestinationModelPart.ConditionsBegin();
//				ic!=DestinationModelPart.ConditionsEnd(); ic++)
//			{
//			
//				if (((*ic).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//					((*ic).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//					((*ic).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   
//				{
//					(*ic).GetGeometry()[0].GetValue(AUX)=0;
//					(*ic).GetGeometry()[1].GetValue(AUX)=0;
//					(*ic).GetGeometry()[2].GetValue(AUX)=0;
//				}
//			}
//			int i=0;
//			for(ModelPart::ConditionsContainerType::iterator ic = DestinationModelPart.ConditionsBegin();
//				ic!=DestinationModelPart.ConditionsEnd(); ic++)
//			{
//				//std::cout<<" HELLO4!"<<std::endl;
//				
//				if (((*ic).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//					((*ic).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//					((*ic).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   
//				{
//					loc_RHS[0]=0.0; loc_RHS[1]=0.0; loc_RHS[2]=0.0;
//				
//					//gather
//					x_vect[0] = (*ic).GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
//					x_vect[1] = (*ic).GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
//					x_vect[2] = (*ic).GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
//
//						
//					array_1d<double,3> p_orig ;
//					p_orig[0]=DestArea[i]*(0.3333333)*OriginScalarValue[i]; p_orig[1]=DestArea[i]*(0.3333333)*OriginScalarValue[i]; p_orig[2]=DestArea[i]*(0.333333)*OriginScalarValue[i];
//					
//					M_cons(0,0)=2;	M_cons(0,1)=1;	M_cons(0,2)=1;
//					M_cons(1,0)=1;	M_cons(1,1)=2;	M_cons(1,2)=1;
//					M_cons(2,0)=1;	M_cons(2,1)=1;	M_cons(2,2)=2;
//
//					M_cons*=DestArea[i]/(12.0);
//					
//					loc_RHS = p_orig - prod(M_cons, x_vect);
//
//					(*ic).GetGeometry()[0].GetValue(AUX)+=loc_RHS[0];
//					(*ic).GetGeometry()[1].GetValue(AUX)+=loc_RHS[1];
//					(*ic).GetGeometry()[2].GetValue(AUX)+=loc_RHS[2];
//				
//					i++;
//				}
//			}
//			//////////////////////////////////////////////////////////////
//			double norm_dx=0.0;
//			int inodes=0;
//			for(ModelPart::NodeIterator id = DestinationModelPart.NodesBegin() ; 
//			id != DestinationModelPart.NodesEnd() ; ++id)
//			{
//				//std::cout<<" HELLO5!"<<std::endl;
//
//				if ((id->FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//						(id->FastGetSolutionStepValue(IS_INTERFACE)==1.0) &&
//						(id->FastGetSolutionStepValue(IS_INTERFACE)==1.0) )   
//				{
//					const double& A_I = id->FastGetSolutionStepValue(NODAL_AREA);
//					id->FastGetSolutionStepValue(PRESSURE)+=id->GetValue(AUX)/A_I;//p+=dp
//					norm_dx+=(id->GetValue(AUX)/A_I)*(id->GetValue(AUX)/A_I);
//					inodes++;
//									
//				}
//			}

//			//std::cout<<norm_dx/inodes<<std::endl;
//			if((norm_dx/inodes)<(tol_iter*tol_iter))break;
//////////////////////////////////////////////////////////////////
//		}
//	std::cout << "  Map Iteration TIME for "<<k+1<<" Iterations = " << Mapiter_timer.elapsed() << std::endl;
//	//std::cout<<"Iterations finished"<<std::endl;
//	}//calc_p_nodal_dest

	
} //namespace Kratos
	
