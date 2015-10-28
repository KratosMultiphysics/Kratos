//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine  $
//   Date:                $Date: 2012-05-24 $
//   Revision:            $Revision: 1.0 $
//


#if !defined(KRATOS_RIGID_FACE_GEOMETRICAL_OBJECT_CONFIGURE_INCLUDED)
#define  KRATOS_RIGID_FACE_GEOMETRICAL_OBJECT_CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>

// Kratos includes
#include "includes/variables.h"
#include "spatial_containers/spatial_search.h"
#include "GeometryFunctions.h"
#include "geometries/geometry.h"

namespace Kratos
{

    
    typedef Geometry<Node < 3 > > GeometryType;
  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

    
template <std::size_t TDimension>
class RigidFaceGeometricalObjectConfigure
{

public:
  
    enum { 
        Dimension = TDimension,
        DIMENSION = TDimension,
        MAX_LEVEL = 16,
        MIN_LEVEL = 2
    };



    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(RigidFaceGeometricalObjectConfigure);
    
	
	 typedef SpatialSearch                                                       SearchType;

    typedef SearchType::PointType                                               PointType;
    typedef PointerVectorSet<GeometricalObject, IndexedObject>::ContainerType   ContainerType;
    typedef PointerVectorSet<GeometricalObject, IndexedObject>                  ElementsContainerType;
    
    typedef SearchType::ElementType                                             ElementType;
    typedef ContainerType::value_type                                           PointerType;
    typedef ContainerType::iterator                                             IteratorType;
    
    typedef PointerVectorSet<GeometricalObject, IndexedObject>::ContainerType   ResultContainerType;

    
    typedef ResultContainerType::iterator                           ResultIteratorType;
    typedef std::vector<double>::iterator                           DistanceIteratorType;
    
    typedef ContactPair<PointerType>                                ContactPairType;
    typedef std::vector<ContactPairType>                            ContainerContactType;
    typedef ContainerContactType::iterator                          IteratorContactType;
    typedef ContainerContactType::value_type                        PointerContactType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    RigidFaceGeometricalObjectConfigure(){};
    virtual ~RigidFaceGeometricalObjectConfigure(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //******************************************************************************************************************

/////Cfeng: For Particle DEM
  static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
        rHighPoint = rLowPoint  = rObject->GetGeometry().GetPoint(0);
        
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -Radius;
            rHighPoint[i] += Radius;
        }
		
		
    }
	
///Cfeng: For FEM conditions
	static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        ///Cfeng:rObject is condition
	
        array_1d<double, 3> Coord;

        double xyz_min[3] = { 1e20,  1e20,  1e20};
        double xyz_max[3] = {-1e20, -1e20, -1e20};

        for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
        {
                Coord = rObject->GetGeometry()[inode].Coordinates();

                xyz_min[0] = (xyz_min[0] > Coord[0]) ? Coord[0] : xyz_min[0];
                xyz_min[1] = (xyz_min[1] > Coord[1]) ? Coord[1] : xyz_min[1];
                xyz_min[2] = (xyz_min[2] > Coord[2]) ? Coord[2] : xyz_min[2];

                xyz_max[0] = (xyz_max[0] < Coord[0]) ? Coord[0] : xyz_max[0];
                xyz_max[1] = (xyz_max[1] < Coord[1]) ? Coord[1] : xyz_max[1];
                xyz_max[2] = (xyz_max[2] < Coord[2]) ? Coord[2] : xyz_max[2];
        }
		
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint [i] = xyz_min[i];
            rHighPoint[i] = xyz_max[i];
        }

        const double domain_size = rObject->GetGeometry().DomainSize();
        
        for(std::size_t i = 0; i < 3; i++)
        {
            if( (rHighPoint[i]-rLowPoint[i]) < 1e-10 * domain_size)  //altura no area
            {
                rHighPoint[i] = rLowPoint[i] + domain_size;
                  //rLowPoint [i];// + domain_size;
            }
        }

    }



    //******************************************************************************************************************

        static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius)
    {
    
      // int NewContactType  = -1; 
      //-1: No contact;
      // 1: Plane;
      // 2: Edge;
      // 3: Point.
      
      bool ContactExists = EasyIntersection(rObj_1, rObj_2, Radius);//, NewContactType);
 
      return ContactExists;
        
    }
       
    static inline bool RigidContact(int method, SphericParticle* rObj_1, DEMWall* rObj_2){ 
    
      bool result = false;
      
      switch (method)
      {
        case 1:
          result = HierarchicalMethod(rObj_1, rObj_2); 
          break;
          
        case 2:
          result = AreaDistribution(rObj_1, rObj_2);
          break;
        
        default:
          KRATOS_THROW_ERROR( std::invalid_argument," error in rigid_face_geometrical_object_configure: Local search resolution method not defined ", "" )
         
      }
      
      return result;
      
    }
    
    static inline bool AreaDistribution(SphericParticle* rObj_1, DEMWall* rObj_2){ //miquel
       
        //ASSUMPTIONS
        //- Fine for small indentations (at least indentation < R)
        //- Triangle elements numbered in anticlockwise fashion. (normals pointing outwards)
        
        //0.INITIALIZATIONS
        int                 casus                   = -1;
        double              Area                    = 0.0;
        double              CoM[3]                  = {0.0}; //centroid of the region
        double              AddCoM[3]               = {0.0}; //cumulative centroid of the region
        double              normal_flag             = 0.0;
        double              Coord[3][3]             = {{0.0},{0.0},{0.0}};
        double              Normal[3]               = {0.0};
        double              LocalCoordSystem[3][3]  = {{0.0},{0.0},{0.0}};
        double              Particle_Coord[3]       = {0.0};
        double              Particle_Radius         = rObj_1->GetRadius();
        double              delta                   = 0.0; //its artificial for the cases where the normal and the Centre Particle are in oposite direction;
        double              distance_cp_to_plane    = 0.0;
        double              CP_V0[3]                = {0.0}; //Vector going from Vertex 0 of the triangle to the Centre of Particle (Particle_Coord)
        double              CC[3]                   = {0.0}; //Centre of the Circle projected on the triangle plane
        double              RC                      = 0.0; //Radius of the projected Circle
        double              RC_squared              = 0.0; //Radius Squared of the projected Circle
        std::vector<int>    vertex_in_vector;
        int                 num_vertex_in           = 0;
        double              tol                     = 1e-8;
        double              tol_RC                  = tol*RC;
        
        //contact local axes
        double              Outwards[3]             = {0.0};
        double              Auxiliar1[3]            = {0.0};
        double              Auxiliar2[3]            = {0.0};
        
        const array_1d<double,3>& acces_to_particle_coordinates = rObj_1->GetGeometry()[0].Coordinates();
        const GeometryType& rElementGeometry                    = rObj_2->GetGeometry();
        
        for (unsigned int index1 = 0; index1<3; index1++){
            
            Particle_Coord[index1] = acces_to_particle_coordinates[index1];   
            
            for (unsigned int index2 = 0; index2<3; index2++){
                Coord[index1][index2] =  rElementGeometry[index1].Coordinates()[index2];
            }
        }
       
        //1.PROJECTION CIRCLE (SPHERE ON TRIANGLE):  Local Axes, Normal, Centre and Radius of the Circle, and Delta indentation
        
        GeometryFunctions::Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem, normal_flag);
        //NOTE: LocalCoordSystem[2] points the sphere while normal is stricly outwards with a counterclockwise sortig of the element nodes
        for (unsigned int index = 0;index<3;index++){
                
            Normal[index] = normal_flag*LocalCoordSystem[2][index];
            CP_V0[index]  = Particle_Coord[index] - Coord[0][index];
            
        }
        
        distance_cp_to_plane  = GeometryFunctions::DotProduct(CP_V0,LocalCoordSystem[2]); //this ignores the fact that a particle can be with an indentation greater than radius (not becouse of great prenetration but becouse of new neighbour with a strong change of normal like in a stair.                  
        delta                 = Particle_Radius - distance_cp_to_plane; 
        RC_squared            = delta*(2*Particle_Radius-delta);
        RC                    = sqrt(RC_squared); 
        
        if(distance_cp_to_plane>Particle_Radius){KRATOS_WATCH("ERROR_NUMERO 22234528383 EN RIGID_FACE_GEO")}//NO HI HA CONTACTE??;
        
        for (unsigned int index = 0;index<3;index++){
                
            CC[index]  = Particle_Coord[index] - LocalCoordSystem[2][index]*distance_cp_to_plane;
        }
                
        //MSICHECK
        if(GeometryFunctions::DotProduct(CP_V0,LocalCoordSystem[2])<0.0)KRATOS_WATCH("ERROR_NUMERO 292928383 EN RIGID_FACE_GEO")
        if(GeometryFunctions::DotProduct(CP_V0,Normal)<0.0)std::cout<<"LA bola amb id "<<rObj_1->Id()<<" esta darrera del triangle amb nodes"<<rObj_2->GetGeometry()[0].Id()<<" "<< rObj_2->GetGeometry()[1].Id() <<" "<<rObj_2->GetGeometry()[2].Id()<<std::endl;
        
        
        //2. VERTICES INSIDE TEST
        
        for (unsigned int index = 0;index<3;index++){        
            double dist_sq = GeometryFunctions::DistanceOfTwoPointSquared(Coord[index],CC);
            if( (dist_sq - RC_squared) < -tol_RC){
                vertex_in_vector.push_back(index);
            }
        }
        
        num_vertex_in = vertex_in_vector.size();
                
        //3. CASE SELECTION
        
        double V0[3]    = {0.0}; double V1[3]  = {0.0}; double V2[3]  = {0.0}; double V0V1[3] = {0.0}; double V0CC[3] = {0.0}; //double proj[3] = {0.0};
        double V0V2[3] = {0.0}; double V2V0[3] = {0.0}; 
        double intersection1[3]     = {0.0};  double intersection2[3]       = {0.0};  
        double normal_outwards1[3]  = {0.0};  double normal_outwards2[3]    = {0.0};
        double direction_edge1[3]   = {0.0};  double direction_edge2[3]     = {0.0};
        double dist_inline1         =  0.0;   double dist_inline2           =  0.0;
        double dist_normal1         =  0.0;   double dist_normal2           =  0.0;
        double pseudo_dist_inline1  =  0.0;   double pseudo_dist_inline2    =  0.0;
        double AreaCS               =  0.0;   double AreaTV                 =  0.0;   double AreaTC = 0.0;  double AreaTriTotal = 0.0;
        double AreaSegC             =  0.0;
        double CoMSC[3]             = {0.0};  double CoMSegC[3]             = {0.0};
        bool   flag                 = false;
        double inv_Area             = 0.0;
        
        if(num_vertex_in == 0){ //All the cases are solved cutting the circular areas out of the edges off.
        
            Area            = KRATOS_M_PI*RC_squared; 
            for (unsigned int index = 0;index<3;index++){        
                AddCoM[index]  = Area*CC[index];
            }
            casus           = 1;
            
            //loop over edges
            for (unsigned int index = 0;index<3;index++){           
                
                for (unsigned int i = 0;i<3;i++){  
                    
                    V0[i] = Coord[index][i];
                    V1[i] = Coord[(index+1)%3][i];
                    V0V1[i] = V1[i]-V0[i];
                    V0CC[i] = CC[i]-V0[i];
                    
                }
                                    
                GeometryFunctions::AreaAndCentroidCircularSegment(CC,RC,tol_RC,V0,V1,Normal,AreaSegC,CoMSegC,flag);
                    
                if(flag){
                    
                    Area = Area - AreaSegC;
                    for (unsigned int index = 0;index<3;index++){  
                        AddCoM[index] = AddCoM[index] - AreaSegC*CoMSegC[index];
                    }
                    casus++;
                }
                    
            } //loop over the edges
            
            inv_Area = 1/Area;
        
            for (unsigned int index = 0; index<3; index++){        

                CoM[index]     = AddCoM[index]*inv_Area;
           
            }   
            
        }//num_vertex_in 0
        
        else if( (num_vertex_in == 1) || (num_vertex_in == 2)){ //this two cases are grouped in one
            
            int index_vertex=-1;
            double CoMTC[3] = {0.0}; double CoMTV[3] = {0.0}; double CoMTT[3] = {0.0};
            
            if(num_vertex_in==2){
                int suma_index = vertex_in_vector[0]+vertex_in_vector[1];
                switch (suma_index){ case 1: index_vertex = 2; break; case 2: index_vertex = 1; break; case 3: index_vertex = 0; break; }
            }
            else{
                index_vertex = vertex_in_vector[0];
            }
            
            for (unsigned int i = 0;i<3;i++){  
                
                V0[i] = Coord[index_vertex][i];
                V1[i] = Coord[(index_vertex+1)%3][i];
                V2[i] = Coord[(index_vertex+2)%3][i];
            
                V0V1[i] = V1[i]-V0[i];
                V0CC[i] = CC[i]-V0[i];
                V0V2[i] = V2[i]-V0[i];
                V2V0[i] = V0[i]-V2[i];

            }
                
            GeometryFunctions::CrossProduct(V0V1,Normal,normal_outwards1); 
            GeometryFunctions::CrossProduct(V2V0,Normal,normal_outwards2);  
                            
            GeometryFunctions::normalize(normal_outwards1);
            GeometryFunctions::normalize(normal_outwards2);
            //MSICHECK borra aixo, pero vov2 es fa servir tambe 
            if(GeometryFunctions::DotProduct(normal_outwards1,V0V2)>0.0) {KRATOS_WATCH("ERROR_NUMERO 292435383 EN RIGID_FACE_GE0 no estic segur que sigui errorO")} //han de ser contraris pk la normal outwards se suposa que es outwards
            
            dist_normal1   = GeometryFunctions::DotProduct(normal_outwards1,V0CC);  
            dist_normal2   = GeometryFunctions::DotProduct(normal_outwards2,V0CC);    

            dist_inline1 = sqrt(RC_squared-dist_normal1*dist_normal1);
            dist_inline2 = sqrt(RC_squared-dist_normal2*dist_normal2);

            if(num_vertex_in==2){ //its diferent if the reference edge is inside or outside
                
                dist_inline1 = -1*dist_inline1;
                dist_inline2 = -1*dist_inline2;
                
            }

            for (unsigned int index = 0;index<3;index++){   
                
                direction_edge1[index] = V0V1[index];
                direction_edge2[index] = V0V2[index];
            }
            
            GeometryFunctions::normalize(direction_edge1);
            GeometryFunctions::normalize(direction_edge2);
            
            pseudo_dist_inline1 = GeometryFunctions::DotProduct(direction_edge1,V0CC); //can be negative if the centre of circle lies out of triangle
            pseudo_dist_inline2 = GeometryFunctions::DotProduct(direction_edge2,V0CC); //can be negative if the centre of circle lies out of triangle
            
            for (unsigned int index = 0;index<3;index++){        
                intersection1[index] = V0[index]+(dist_inline1+pseudo_dist_inline1)*direction_edge1[index];
                intersection2[index] = V0[index]+(dist_inline2+pseudo_dist_inline2)*direction_edge2[index];

            }
            
            GeometryFunctions::AreaAndCentroidCircularSector(CC, RC, intersection1, intersection2, Normal, AreaCS, CoMSC);
            GeometryFunctions::AreaAndCentroidTriangle(CC, intersection1, intersection2, AreaTC, CoMTC);
            GeometryFunctions::AreaAndCentroidTriangle(V0, intersection1, intersection2, AreaTV, CoMTV);
            
            
            if(num_vertex_in==2){
                
                casus = 7;
                
                GeometryFunctions::AreaAndCentroidTriangle(Coord[0],Coord[1],Coord[2],AreaTriTotal,CoMTT);
                
                Area = AreaTriTotal - AreaTV + AreaCS -AreaTC;
                inv_Area = 1/Area;
                
                for (unsigned int index = 0; index<3; index++){        
                    CoM[index]     = inv_Area*(AreaTriTotal*CoMTT[index] - AreaTV*CoMTV[index] + AreaCS*CoMSC[index] - AreaTC*CoMTC[index]);
                }
            
            }//num_vertex_in 2
            else{ //for the case of 1 vertex, we check now if we cross the oposite edge.
                
                GeometryFunctions::AreaAndCentroidCircularSegment(CC,RC,tol_RC,V1,V2,Normal,AreaSegC,CoMSegC,flag);
                    
                if(flag){
                    
                    casus = 6;
                    
                    Area = AreaTV + AreaCS - AreaTC - AreaSegC;
                    inv_Area = 1/Area;
                    for (unsigned int index = 0;index<3;index++){  
                        
                        CoM[index] = inv_Area*(AreaTV*CoMTV[index] + AreaCS*CoMSC[index] - AreaTC*CoMTC[index] - AreaSegC*CoMSegC[index]);   
                    
                    } 
                }//opposite edge crossed
                else{
                
                    casus = 5;
                    Area = AreaTV + AreaCS - AreaTC;
                    inv_Area = 1/Area;
                    for (unsigned int index = 0;index<3;index++){  
                        
                        CoM[index] = inv_Area*(AreaTV*CoMTV[index] + AreaCS*CoMSC[index] - AreaTC*CoMTC[index]);   
                    
                    } 
                }
                
            }//1 vertex
            
        }//1 or 2 vertices in
        else if(num_vertex_in == 3){
            
            GeometryFunctions::AreaAndCentroidTriangle(Coord[0],Coord[1],Coord[2],Area,CoM);
            casus = 8;
            
        }
   
        double Weight[3] = {0.0};
        GeometryFunctions::TriAngleWeight(Coord[0], Coord[1], Coord[2], CoM, Weight);

        if(Area<0.0)
        {
            KRATOS_WATCH("ERROR_NUMERO 88942132432484 EN RIGID_FACE_GE0")
            KRATOS_WATCH(Area)
            return false;
        }
        
        if( fabs(Weight[0] + Weight[1] + Weight[2] - 1.0) >= 1.0e-5 )
        {
            KRATOS_WATCH("ERROR_NUMERO 8894845484 EN RIGID_FACE_GE0")
            KRATOS_WATCH(Weight[0])
            KRATOS_WATCH(Weight[1])
            KRATOS_WATCH(Weight[2])
            KRATOS_WATCH(Area)  
            KRATOS_WATCH(rObj_2->Id()) 
        KRATOS_WATCH(casus)  
        return false;
        }
        
        std::vector<double>& RF_Pram = rObj_1->mNeighbourRigidFacesPram;
        std::size_t ino = RF_Pram.size();
        
        RF_Pram.resize(ino + 16);

        Outwards[0]   = Particle_Coord[0]-CoM[0];
        Outwards[1]   = Particle_Coord[1]-CoM[1];
        Outwards[2]   = Particle_Coord[2]-CoM[2];
        
        Auxiliar2[0]  = Outwards[0]+0.1;
        Auxiliar2[1]  = Outwards[1]+0.2;
        Auxiliar2[2]  = Outwards[2]+0.3;
        
        GeometryFunctions::CrossProduct(Outwards,Auxiliar2,Auxiliar1);
        GeometryFunctions::CrossProduct(Outwards,Auxiliar1,Auxiliar2); 
                
        GeometryFunctions::normalize(Auxiliar1);
        GeometryFunctions::normalize(Auxiliar2);
        GeometryFunctions::normalize(Outwards);
 
        RF_Pram[ino + 0]  = Auxiliar1[0];
        RF_Pram[ino + 1]  = Auxiliar1[1];
        RF_Pram[ino + 2]  = Auxiliar1[2];
        RF_Pram[ino + 3]  = Auxiliar2[0];
        RF_Pram[ino + 4]  = Auxiliar2[1];
        RF_Pram[ino + 5]  = Auxiliar2[2];
        RF_Pram[ino + 6]  = Outwards[0]; //the result points outwards with a counterclockwise sorting of nodes criteria
        RF_Pram[ino + 7]  = Outwards[1];
        RF_Pram[ino + 8]  = Outwards[2];
        
        double dist = GeometryFunctions::DistanceOfTwoPoint(CoM,Particle_Coord);

        RF_Pram[ino + 9]  = dist;
        RF_Pram[ino + 10] = Weight[0];
        RF_Pram[ino + 11] = Weight[1];
        RF_Pram[ino + 12] = Weight[2];
        RF_Pram[ino + 13] = 0.0;  //Weight[3]... no 4 nodes
        RF_Pram[ino + 14] = rObj_2->Id();
        RF_Pram[ino + 15] = 1; //ContactType always plane  

        return true;
        
    } //AreaDistribution (Elastic Contact)
    
    
    
   
    
    static inline bool HierarchicalMethod(SphericParticle* rObj_1, DEMWall* rObj_2) //joaquin
    {
      //Cfeng: rObj_1 is particle,  and rObj_2 is condition

      bool ContactExists = false;

      int ContactType = -1;
      //-1: No contact;
      // 1: Plane;
      // 2: Edge;
      // 3: Point.
      
      double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
      double Weight[4]              = {0.0, 0.0, 0.0, 0.0};
      double DistPToB               = 0.0;
      double Particle_Coord[3]      = {0.0};
      Particle_Coord[0] = rObj_1->GetGeometry()[0].Coordinates()[0];
      Particle_Coord[1] = rObj_1->GetGeometry()[0].Coordinates()[1];
      Particle_Coord[2] = rObj_1->GetGeometry()[0].Coordinates()[2];

      //double rad = Radius;
      double rad = rObj_1->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

      int i_size = rObj_2->GetGeometry().size();
      
      array_1d<double, 3>& coordinates_obj2_0 = rObj_2->GetGeometry()[0].Coordinates();
      array_1d<double, 3>& coordinates_obj2_1 = rObj_2->GetGeometry()[1].Coordinates();
      
      //CONTACT PARTICLE-LINE 2D      
      if(i_size == 2)
      {
        double Coord1[3]     = {0.0};
        double Coord2[3]     = {0.0};
        double tempWeight[2] = {0.0};
        
        Coord1[0] = coordinates_obj2_0[0];
        Coord1[1] = coordinates_obj2_0[1];
        Coord1[2] = coordinates_obj2_0[2];

        Coord2[0] = coordinates_obj2_1[0];
        Coord2[1] = coordinates_obj2_1[1];
        Coord2[2] = coordinates_obj2_1[2];

        ContactExists = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad,
                                                                                  LocalCoordSystem,  tempWeight, DistPToB);
        
        if(ContactExists == true)
        {
            ContactType = 2;
            Weight[0] = tempWeight[0];
            Weight[1] = tempWeight[1];      
        }
        
        else /// particle point contact search
        {
            for(std::size_t inode1 = 0; inode1 < rObj_2->GetGeometry().size(); inode1++)
            {
                Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
                
                ContactExists = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, LocalCoordSystem, DistPToB);
               
                if(ContactExists == true)
                {
                    ContactType = 3;
                    Weight[inode1] = 1.0;           
                    break;
                }
                
            }       
        }
        
         
      }  //if(rObj_2->GetGeometry().size() == 2)
      
      
      //CONTACT PARTICLE-TRIANGLE/QUADRILATERAL 3D
      
      else if(rObj_2->GetGeometry().size() > 2)
      {
          double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };

          ////Cfeng: Triangle
          int FaceNodeTotal = 3;
  
          Coord[0][0] = coordinates_obj2_0[0];
          Coord[0][1] = coordinates_obj2_0[1];
          Coord[0][2] = coordinates_obj2_0[2];
 
          Coord[1][0] = coordinates_obj2_1[0];
          Coord[1][1] = coordinates_obj2_1[1];
          Coord[1][2] = coordinates_obj2_1[2];

          Coord[2][0] = rObj_2->GetGeometry()[2].Coordinates()[0];
          Coord[2][1] = rObj_2->GetGeometry()[2].Coordinates()[1];
          Coord[2][2] = rObj_2->GetGeometry()[2].Coordinates()[2];

          if(rObj_2->GetGeometry().size() == 4)
          {
              Coord[3][0] = rObj_2->GetGeometry()[3].Coordinates()[0];
              Coord[3][1] = rObj_2->GetGeometry()[3].Coordinates()[1];
              Coord[3][2] = rObj_2->GetGeometry()[3].Coordinates()[2];
 
              ////Cfeng: Quadrilateral
              FaceNodeTotal = 4;
           }
           
           /////Particle-Face contact
           ContactExists = GeometryFunctions::JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord, Particle_Coord, rad, LocalCoordSystem, Weight, DistPToB);

           if(ContactExists == true)
           {
               ContactType = 1;
           }

           ///Particle-edge contact and Particle-point
           if(ContactExists == false)
           {
               Weight[0] = Weight[1] = Weight[2] = Weight[3] = 0.0;
            
               bool Exist_Contact[8] = {false};
               double Dist[8] = {0.0};
            
               for(int inode1 = 0; inode1 < FaceNodeTotal; inode1++)
               {
                   int inode2 = (inode1 + 1) % FaceNodeTotal;
                
                   double Coord1[3]     = {0.0};
                   double Coord2[3]     = {0.0};
                
                   Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                   Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                   Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
                
                   Coord2[0] = rObj_2->GetGeometry()[inode2].Coordinates()[0];
                   Coord2[1] = rObj_2->GetGeometry()[inode2].Coordinates()[1];
                   Coord2[2] = rObj_2->GetGeometry()[inode2].Coordinates()[2];
                   
                   Exist_Contact[inode1]   = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad,  Dist[inode1]); //simplified judgement
                   if (Exist_Contact[inode1] == true) {
                       Exist_Contact[inode1+4] = false;
                   }
                   else {
                       Exist_Contact[inode1+4] = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, Dist[inode1+4]); //simplified judgement
                   }
               }
            
               int contact_index = -1;
               double dist = 0.0;
            
               for(int i = 0; i < 8; i++)
               {
                   if(Exist_Contact[i] == true)
                   {
                       if (contact_index == -1 || Dist[i] < dist)
                       {
                           dist = Dist[i];
                           contact_index = i;
                       }
                   }
               }
            
               if (contact_index == 0 || contact_index == 1 || contact_index == 2 || contact_index == 3) //hierarchy on edges
               {
                   ContactType = 2;
                   double Coord1[3]     = {0.0};
                   double Coord2[3]     = {0.0};
                   double tempWeight[2] = {0.0};
                
                   int inode1 = contact_index;
                   int inode2 = (inode1 + 1) % FaceNodeTotal;
                   int inode3 = (inode1 + 2) % FaceNodeTotal;
                
                   Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                   Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                   Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
                
                   Coord2[0] = rObj_2->GetGeometry()[inode2].Coordinates()[0];
                   Coord2[1] = rObj_2->GetGeometry()[inode2].Coordinates()[1];
                   Coord2[2] = rObj_2->GetGeometry()[inode2].Coordinates()[2];
                
                   ContactExists = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad, LocalCoordSystem, tempWeight, DistPToB); //complete judgement
                
                   Weight[inode1] = tempWeight[0];
                   Weight[inode2] = tempWeight[1];
                   Weight[inode3] = 0.0;
                
                   if(FaceNodeTotal == 4)
                   {
                       int inode4 = (inode1 + 3) % FaceNodeTotal;
                       Weight[inode4] = 0.0;
                   }
               }
            
               if (contact_index == 4 || contact_index == 5 || contact_index == 6 || contact_index == 7) //hierarchy on nodes
               {
                   ContactType = 3;
                
                   int inode1 = contact_index-4;
                 
                   double Coord1[3] = {0.0};
                   Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                   Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                   Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
              
                   ContactExists = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, LocalCoordSystem, DistPToB); //complete judgement

                   Weight[inode1] = 1.0;
                
               }
           } //if(ContactExists == false)
       } //if(rObj_2->GetGeometry().size() > 2)
     
       ////Store the neighbour value, real contact rigidFace
       if(ContactExists == true)
       {
            // In geometrical level to search the neighbour, need pointer convert
            std::vector<double>& RF_Pram = rObj_1->mNeighbourRigidFacesPram;     
            
            //NOTE: THIS member has to be cleared every timestep. It is done in: the SphericParticle Function: ComputeNewRigidFaceNeighboursHistoricalData 

            std::size_t ino = RF_Pram.size();
            std::size_t TotalSize = ino / 16;        
            std::size_t isize;
        
              RF_Pram.resize(ino + 16);
              RF_Pram[ino + 0]  = LocalCoordSystem[0][0];//Room for improvement: ¿could we just store LocalCoordSystem[2]? Could we regenerate 0 and 1 randomly when we need the coordinate system?
              RF_Pram[ino + 1]  = LocalCoordSystem[0][1];
              RF_Pram[ino + 2]  = LocalCoordSystem[0][2];
              RF_Pram[ino + 3]  = LocalCoordSystem[1][0];
              RF_Pram[ino + 4]  = LocalCoordSystem[1][1];
              RF_Pram[ino + 5]  = LocalCoordSystem[1][2];
              RF_Pram[ino + 6]  = LocalCoordSystem[2][0];
              RF_Pram[ino + 7]  = LocalCoordSystem[2][1];
              RF_Pram[ino + 8]  = LocalCoordSystem[2][2];
              RF_Pram[ino + 9]  = DistPToB;
              RF_Pram[ino + 10] = Weight[0];
              RF_Pram[ino + 11] = Weight[1];
              RF_Pram[ino + 12] = Weight[2];
              RF_Pram[ino + 13] = Weight[3];  
              RF_Pram[ino + 14] = rObj_2->Id();
              RF_Pram[ino + 15] = ContactType;

              for(isize = 0; isize < TotalSize; isize++)
              {
                  int ino1 = isize * 16;
          
                  if ( RF_Pram[ino1 + 15] == 1 || RF_Pram[ino1 + 15] == 2 || RF_Pram[ino1 + 15] == 3 )
                  {
                      double Vector[3] = {0.0};
                      Vector[0] = RF_Pram[ino1 + 6];
                      Vector[1] = RF_Pram[ino1 + 7];
                      Vector[2] = RF_Pram[ino1 + 8];
              
                      double Normal_dist1 = GeometryFunctions::DotProduct(LocalCoordSystem[2], Vector) * DistPToB;

                      if ( (Normal_dist1-RF_Pram[ino1 + 9])/fabs(RF_Pram[ino1 + 9]) > -1.0e-15 )
                      {
                          RF_Pram[ino + 15] = -1;
                          break;
                      }
              
                      double Normal_dist2 = GeometryFunctions::DotProduct(LocalCoordSystem[2], Vector) * RF_Pram[ino1 + 9];
              
                      if ( (Normal_dist2-DistPToB)/fabs(DistPToB) > -1.0e-15 )
                      {
                          RF_Pram[ino1 + 15] = -1;
                      }
                  }
              } 
          }

      return ContactExists;
      
    }

    //******************************************************************************************************************

     static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
     {
         //Cfeng: rObject is block

    
    array_1d<double, 3> Coord;

    double xyz_min[3] = { 1e20,  1e20,  1e20};
    double xyz_max[3] = {-1e20, -1e20, -1e20};
                
    for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
    {
      Coord = rObject->GetGeometry()[inode].Coordinates();

      xyz_min[0] = (xyz_min[0] > Coord[0]) ? Coord[0] : xyz_min[0];
      xyz_min[1] = (xyz_min[1] > Coord[1]) ? Coord[1] : xyz_min[1];
      xyz_min[2] = (xyz_min[2] > Coord[2]) ? Coord[2] : xyz_min[2];

      xyz_max[0] = (xyz_max[0] < Coord[0]) ? Coord[0] : xyz_max[0];
      xyz_max[1] = (xyz_max[1] < Coord[1]) ? Coord[1] : xyz_max[1];
      xyz_max[2] = (xyz_max[2] < Coord[2]) ? Coord[2] : xyz_max[2];
    }
    
    
    bool intersect = (rLowPoint [0] <= xyz_max[0] && rLowPoint [1] <= xyz_max[1] && rLowPoint [2] <= xyz_max[2] &&
                      rHighPoint[0] >= xyz_min[0] && rHighPoint[1] >= xyz_min[1] && rHighPoint[2] >= xyz_min[2]);
            

        return  intersect;
  
  
    
     }
   
   
  
  static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
        array_1d<double, 3> center_of_particle = rObject->GetGeometry().GetPoint(0);

        double radius = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        bool intersect = (
          floatle(rLowPoint[0]  - radius,center_of_particle[0]) && 
          floatle(rLowPoint[1]  - radius,center_of_particle[1]) && 
          floatle(rLowPoint[2]  - radius,center_of_particle[2]) &&
          floatge(rHighPoint[0] + radius,center_of_particle[0]) && 
          floatge(rHighPoint[1] + radius,center_of_particle[1]) && 
          floatge(rHighPoint[2] + radius,center_of_particle[2]));
      
      
        return  intersect;
    }
  
  
  
  static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        array_1d<double, 3> center_of_particle1 = rObj_1->GetGeometry().GetPoint(0);
        array_1d<double, 3> center_of_particle2 = rObj_2->GetGeometry().GetPoint(0);
    

        distance = sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
                        (center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
                        (center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    }
     
   
   //*******************************************************************************************************************************************************//
   //MIQUEL NEW FUNCTIONS
   //*******************************************************************************************************************************************************
   
    static inline bool EasyIntersection(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius){//, int& ContactType){
        //rObj_1 is particle,  and rObj_2 is condition

        //int ContactType = -1;
        //-1: No contact;
        // 1: Plane;
        // 2: Edge;
        // 3: Point.
        // 4: Edge or Point

        bool ContactExists = false;

        double Particle_Coord[3] = {0.0};
        Particle_Coord[0] = rObj_1->GetGeometry()[0].Coordinates()[0];
        Particle_Coord[1] = rObj_1->GetGeometry()[0].Coordinates()[1];
        Particle_Coord[2] = rObj_1->GetGeometry()[0].Coordinates()[2];

        double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };
        
        array_1d<double, 3>& coordinates_obj2_0 = rObj_2->GetGeometry()[0].Coordinates();
        array_1d<double, 3>& coordinates_obj2_1 = rObj_2->GetGeometry()[1].Coordinates();
        array_1d<double, 3>& coordinates_obj2_2 = rObj_2->GetGeometry()[2].Coordinates();
        
        Coord[0][0] = coordinates_obj2_0[0];
        Coord[0][1] = coordinates_obj2_0[1];
        Coord[0][2] = coordinates_obj2_0[2];

        Coord[1][0] = coordinates_obj2_1[0];
        Coord[1][1] = coordinates_obj2_1[1];
        Coord[1][2] = coordinates_obj2_1[2];

        Coord[2][0] = coordinates_obj2_2[0];
        Coord[2][1] = coordinates_obj2_2[1];
        Coord[2][2] = coordinates_obj2_2[2];

        //TriangleEasyIntersection
        int FaceNodeTotal = 3;
        int facet_size = rObj_2->GetGeometry().size();
        
        if( facet_size==4 ) //for other convex polygons, see: Meyer, Lee, Barr, Desbrun
        {
            Coord[3][0] = rObj_2->GetGeometry()[3].Coordinates()[0];
            Coord[3][1] = rObj_2->GetGeometry()[3].Coordinates()[1];
            Coord[3][2] = rObj_2->GetGeometry()[3].Coordinates()[2];

            FaceNodeTotal = 4;
        }
  
        //Particle-Face contact
        
        double distance_point_to_plane        = 0.0;
        double dummy_local_coord_system[3][3] = { {0.0},{0.0},{0.0} }; //MSIMSI comproba que estigui ben inicialitzat. he copiat duna altre lloc
        
        ContactExists = GeometryFunctions::Fast_JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord,  Particle_Coord, Radius, dummy_local_coord_system, distance_point_to_plane);
            
        //The key here is to see that we only need to check for further contact if, when not having contact with plane, the distance_point_to_plane is lower than the radius. 
        //In this case it might have contact with edges or vertices, otherwise no contact is possible. 
                
        ///Particle-edge contact and Particle-point
        if( (ContactExists == false) && (distance_point_to_plane < Radius ) )
        {         
            ContactExists = QuickClosestEdgeVertexDetermination(/*ContactType ,*/Coord, Particle_Coord, facet_size, Radius);

        }//no plane contact found

        return ContactExists;
        
    }//EasyIntersection
    
    
     static inline bool QuickClosestEdgeVertexDetermination(/*int& ContactType,*/ double Coord[4][3], double Particle_Coord[3], int size, double particle_radius) 
    {                  
        
        std::vector<double> dist_sq(size*2, 0.0); //vertices and edges of the face       
        std::vector<double> dist(size*2, 0.0);
        double base[4][3];
        double base_sq    = 0.0; 
        double tau_sq     = 0.0;
        double base_leng  = 0.0;
        double diag_back  = 0.0;
        double diag_front = 0.0;
        double particle_radius_sq = particle_radius*particle_radius;
        
        double NodeToP[4][3] = { {0.0},{0.0},{0.0},{0.0} };
        
        //VERTICES
        for (int n=0;n<size;n++){ //vertices
        
            
            for (int k= 0; k<3; k++){ //components
            
                NodeToP[n][k] =  Particle_Coord[k] - Coord[n][k];
                dist_sq[n]    += NodeToP[n][k]*NodeToP[n][k];
                
            }
            
            if ( dist_sq[n] <= particle_radius_sq )
            {
                //ContactType = 3;
                return true; 
            }
            
            dist[n] = sqrt(dist_sq[n]);
            
        }
       
        //EDGES
        for (int i=0;i<size;i++)
        {
           
            int j=(i+1)%size;
            base_sq = 0.0;
            tau_sq  = 0.0;
            for (int k=0;k<3;k++)
            {
                base[i][k] = (Coord[j][k] - Coord[i][k]); 
                base_sq += base[i][k]*base[i][k];
            
            }
                    
            base_leng = sqrt(base_sq);
            diag_back  = dist[i];
            diag_front = dist[j];
            
            double a2  =  dist_sq[i];
            double b2  =  base_sq;
            double c2  =  dist_sq[j];
           
            
            if ( (a2 <= b2 + c2) && (c2 <= a2 + b2) ){  //neglecting obtuse triangles.
                
                GeometryFunctions::TauSqCalculation(tau_sq,base_leng,diag_back,diag_front);
                
                dist_sq[i+size] = tau_sq/base_sq; 
                
                if ( dist_sq[i+size] <= particle_radius_sq )
                {         
                    //ContactType=2;
                    return true;
                    
                }
            
            }
            
        }//for each edge starting at node i going to j

        return false;
     
    } //QuickClosestEdgeVertexDetermination
    
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
    virtual std::string Info() const {return " Spatial Containers Configure for RigidFace"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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

    ///@}
    ///@name Protected Operators
    ///@{

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
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
      
    static inline bool floateq(double a, double b) {
        return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
    }
    
    static inline bool floatle(double a, double b) {
        return std::fabs(a - b) < std::numeric_limits<double>::epsilon() || a < b;
    }
    
    static inline bool floatge(double a, double b) {
        return std::fabs(a - b) < std::numeric_limits<double>::epsilon() || a > b;
    }

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
    RigidFaceGeometricalObjectConfigure& operator=(RigidFaceGeometricalObjectConfigure const& rOther);

    /// Copy constructor.
    RigidFaceGeometricalObjectConfigure(RigidFaceGeometricalObjectConfigure const& rOther);

    ///@}

    }; // Class ParticleConfigure

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template <std::size_t TDimension>
    inline std::istream& operator >> (std::istream& rIStream, RigidFaceGeometricalObjectConfigure<TDimension> & rThis){
        return rIStream;
        }

    /// output stream function
    template <std::size_t TDimension>
    inline std::ostream& operator << (std::ostream& rOStream, const RigidFaceGeometricalObjectConfigure<TDimension>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }
        
    ///@}

}   // namespace Kratos.
#endif  /* rigid_face_geometrical_object_configure.h */
