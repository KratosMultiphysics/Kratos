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

namespace Kratos
{

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

        for(std::size_t i = 0; i < 3; i++)
        {
            if( (rHighPoint[i]-rLowPoint[i]) < 1e-10*rObject->GetGeometry().DomainSize())  //altura no area
            {
                rHighPoint[i] = rLowPoint[i] + rObject->GetGeometry().DomainSize();
                  //rLowPoint [i];// + rObject->GetGeometry().DomainSize();
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
       
    
    
    static inline bool RigidContact(SphericParticle* rObj_1, DEMWall* rObj_2)
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
           ContactExists = JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord, Particle_Coord, rad, LocalCoordSystem, Weight, DistPToB);

           if(ContactExists == true)
           {
               ContactType = 1;
           }

           ///Particle-edge contact and Particle-point
           if(ContactExists == false)
           {
               Weight[0] = Weight[1] = Weight[2] = Weight[3] = 0.0;
            
               bool Contact[8] = {false};
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
                
                   Contact[inode1]   = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad,  Dist[inode1]);
                   Contact[inode1+4] = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, Dist[inode1+4]);
               }
            
               int contact = -1;
               double dist = 0.0;
            
               for(int i = 0; i < 8; i++)
               {
                   if(Contact[i] == true)
                   {
                       if (contact == -1 || Dist[i] < dist)
                       {
                           dist = Dist[i];
                           contact = i;
                       }
                   }
               }
            
               if (contact == 0 || contact == 1 || contact == 2 || contact == 3)
               {
                   ContactType = 2;
                   double Coord1[3]     = {0.0};
                   double Coord2[3]     = {0.0};
                   double tempWeight[2] = {0.0};
                
                   int inode1 = contact;
                   int inode2 = (inode1 + 1) % FaceNodeTotal;
                   int inode3 = (inode1 + 2) % FaceNodeTotal;
                
                   Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                   Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                   Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
                
                   Coord2[0] = rObj_2->GetGeometry()[inode2].Coordinates()[0];
                   Coord2[1] = rObj_2->GetGeometry()[inode2].Coordinates()[1];
                   Coord2[2] = rObj_2->GetGeometry()[inode2].Coordinates()[2];
                
                   ContactExists = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad, LocalCoordSystem, tempWeight, DistPToB);
                
                   Weight[inode1] = tempWeight[0];
                   Weight[inode2] = tempWeight[1];
                   Weight[inode3] = 0.0;
                
                   if(FaceNodeTotal == 4)
                   {
                       int inode4 = (inode1 + 3) % FaceNodeTotal;
                       Weight[inode4] = 0.0;
                   }
               }
            
               if (contact == 4 || contact == 5 || contact == 6 || contact == 7)
               {
                   ContactType = 3;
                
                   int inode1 = contact-4;
                 
                   double Coord1[3] = {0.0};
                   Coord1[0] = rObj_2->GetGeometry()[inode1].Coordinates()[0];
                   Coord1[1] = rObj_2->GetGeometry()[inode1].Coordinates()[1];
                   Coord1[2] = rObj_2->GetGeometry()[inode1].Coordinates()[2];
              
                   ContactExists = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, LocalCoordSystem, DistPToB);

                   Weight[inode1] = 1.0;
                
               }
           } //if(ContactExists == false)
       } //if(rObj_2->GetGeometry().size() > 2)
     
       ////Store the neighbour value, real contact rigidFace
       if(ContactExists == true)
       {
            // In geometrical level to search the neighbour, need pointer convert
            std::vector<double>& RF_Pram = rObj_1->mNeighbourRigidFacesPram;

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
        
        if( facet_size==4 )
        {
            Coord[3][0] = rObj_2->GetGeometry()[3].Coordinates()[0];
            Coord[3][1] = rObj_2->GetGeometry()[3].Coordinates()[1];
            Coord[3][2] = rObj_2->GetGeometry()[3].Coordinates()[2];

            //Cfeng: Quadrilateral
            FaceNodeTotal = 4;
        }

        
                
        //Particle-Face contact
        
        double distance_point_to_plane        = 0.0;
        double dummy_local_coord_system[3][3] = { {0.0},{0.0},{0.0} }; //MSIMSI comproba que estigui ben inicialitzat. he copiat duna altre lloc
        
        ContactExists = Fast_JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord,  Particle_Coord, Radius, dummy_local_coord_system, distance_point_to_plane);
            
        //if(ContactExists) { ContactType = 1; }
        //The key here is to see that we only need to check for further contact if not having contact with plane, the distance_point_to_plane is lower than the radius. 
        //In this case it might have contact with edges or vertices, otherwise no contact is possible. 
                
        ///Particle-edge contact and Particle-point
        if( (ContactExists == false) && (distance_point_to_plane <= Radius ) )
        {
         
            ContactExists = QuickClosestEdgeVertexDetermination(/*ContactType ,*/Coord, Particle_Coord, facet_size, Radius);
            
        }//no plane contact found

        return ContactExists;
        
    }//EasyIntersection
   
   
   
    static inline bool JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Particle_Coord[3], double rad,
                                                          double LocalCoordSystem[3][3], double Weight[4], double &DistPToB)
    {

      bool If_Contact = false;
      
      GeometryFunctions::Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem);

      ///Cfeng,131007,normal vector should point to particle
//       LocalCoordSystem[0][0] = -LocalCoordSystem[0][0];
//       LocalCoordSystem[0][1] = -LocalCoordSystem[0][1];
//       LocalCoordSystem[0][2] = -LocalCoordSystem[0][2];
// 
//       LocalCoordSystem[1][0] = -LocalCoordSystem[1][0];
//       LocalCoordSystem[1][1] = -LocalCoordSystem[1][1];
//       LocalCoordSystem[1][2] = -LocalCoordSystem[1][2];
// 
//       LocalCoordSystem[2][0] = -LocalCoordSystem[2][0];
//       LocalCoordSystem[2][1] = -LocalCoordSystem[2][1];
//       LocalCoordSystem[2][2] = -LocalCoordSystem[2][2];

      DistPToB = GeometryFunctions::DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

      if(DistPToB < rad )
      {
   
           double IntersectionCoord[3];
           GeometryFunctions::CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);
          
            //if(FaceNodeTotal == 3)
           // {
          
              double TriWeight[3] = {0.0};
              GeometryFunctions::TriAngleWeight(Coord[0], Coord[1], Coord[2], IntersectionCoord, TriWeight);

              if( fabs(TriWeight[0] + TriWeight[1] + TriWeight[2] - 1.0) < 1.0e-15 )
              {
                  Weight[0] = TriWeight[0];
                  Weight[1] = TriWeight[1];
                  Weight[2] = TriWeight[2];

                  If_Contact = true;
              }
              

          //}
          
          
          /*else if(FaceNodeTotal == 4)
          {
            
            
              double FaceArea;
              double TempFace[2];
              GeometryFunctions::TriAngleArea(Coord[0], Coord[1], Coord[2], TempFace[0]);
              GeometryFunctions::TriAngleArea(Coord[2], Coord[3], Coord[0], TempFace[1]);
              FaceArea = TempFace[0] + TempFace[1];

              double AreaComponent[4] = {0.0};
              GeometryFunctions::TriAngleArea(IntersectionCoord, Coord[0], Coord[1], AreaComponent[0]);
              GeometryFunctions::TriAngleArea(IntersectionCoord, Coord[1], Coord[2], AreaComponent[1]);
              GeometryFunctions::TriAngleArea(IntersectionCoord, Coord[2], Coord[3], AreaComponent[2]);
              GeometryFunctions::TriAngleArea(IntersectionCoord, Coord[3], Coord[0], AreaComponent[3]);

              if(fabs( (AreaComponent[0] + AreaComponent[1] + AreaComponent[2] + AreaComponent[3] - FaceArea) / FaceArea) < 1.0e-15)
              {
                
                  If_Contact = true;

                  GeometryFunctions::CalQuadWeightCoefficient(Coord, LocalCoordSystem, IntersectionCoord, Weight);
                  
                  
              }
              
            }*/
        }

        return If_Contact;
          
     } //JudgeIfThisFaceIsContactWithParticle;


    //M.Santasusana Februar 2015     
    static inline bool Fast_JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Particle_Coord[3], double rad,
                                                                double LocalCoordSystem[3][3], double &DistPToB)
    {

        bool If_Contact = false;

        GeometryFunctions::Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem);

        ///Cfeng,131007,normal vector should point to particle
//         LocalCoordSystem[0][0] = -LocalCoordSystem[0][0];
//         LocalCoordSystem[0][1] = -LocalCoordSystem[0][1];
//         LocalCoordSystem[0][2] = -LocalCoordSystem[0][2];
// 
//         LocalCoordSystem[1][0] = -LocalCoordSystem[1][0];
//         LocalCoordSystem[1][1] = -LocalCoordSystem[1][1];
//         LocalCoordSystem[1][2] = -LocalCoordSystem[1][2];
// 
//         LocalCoordSystem[2][0] = -LocalCoordSystem[2][0];
//         LocalCoordSystem[2][1] = -LocalCoordSystem[2][1];
//         LocalCoordSystem[2][2] = -LocalCoordSystem[2][2];

        DistPToB = GeometryFunctions::DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

        if(DistPToB < rad )
        {

            double IntersectionCoord[3];
            GeometryFunctions::CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);
            
            If_Contact = Fast_PointInsideTriangle(Coord[0], Coord[1], Coord[2], IntersectionCoord);                  

            if( (FaceNodeTotal == 4) && (If_Contact==false)) //maybe not found inside first triangle, then check second triangular half of the quadrilateral
            {
            
                bool If_Contact2 = Fast_PointInsideTriangle(Coord[0], Coord[2], Coord[3], IntersectionCoord);                  
                
                if( (If_Contact == true) || (If_Contact2 == true)) { If_Contact = true; }
                    
                    
            }
            
                        
        }

        return If_Contact;
        
    } //JudgeIfThisFaceIsContactWithParticle;

   
    
   
    static inline  bool SameSide(double Coord1[3], double Coord2[3], double ReferenceCoord[3], double JudgeCoord[3]){
        
        double cp1[3]  = {0.0};
        double cp2[3]  = {0.0};
        double b_a[3]  = {0.0};
        double p1_a[3] = {0.0};
        double p2_a[3] = {0.0};
        
        for (int i=0;i<3;i++){
            
            b_a[i]  = Coord2[i]-Coord1[i];
            p1_a[i] = JudgeCoord[i]-Coord1[i];
            p2_a[i] = ReferenceCoord[i]-Coord1[i];
            
        }
        
        GeometryFunctions::CrossProduct(b_a, p1_a, cp1);
        GeometryFunctions::CrossProduct(b_a, p2_a, cp2);
        
        if (GeometryFunctions::DotProduct(cp1, cp2) >= 0) return true;
            
        else return false;
        
    }//SameSide


    static inline  bool Fast_PointInsideTriangle(double Coord1[3], double Coord2[3], double Coord3[3], double JudgeCoord[3]){
        
        if( SameSide(Coord1, Coord2, Coord3, JudgeCoord) == false ) return false;
        if( SameSide(Coord2, Coord3, Coord1, JudgeCoord) == false ) return false;
        if( SameSide(Coord3, Coord1, Coord2, JudgeCoord) == false ) return false;
                               
        return true;    
        
    }
    
    
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
            double c2  =  dist_sq[i];
            
            
            if ( (a2 <= b2 + c2) || (c2 <= a2 + b2) ){  //neglecting obtuse triangles.
                
                TauSqCalculation(tau_sq,base_leng,diag_back,diag_front);
                
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

    
    static inline void TauSqCalculation(double& tau_sq, double d0, double d1, double d2)
    {
      
        tau_sq = 0.25*( d0 + d1 - d2)*( d0 - d1 + d2)*( -d0 + d1 + d2)*( d0 + d1 + d2);
        
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
