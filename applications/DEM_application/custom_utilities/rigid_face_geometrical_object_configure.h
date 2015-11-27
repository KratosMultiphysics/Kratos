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
    
    static inline bool FastIntersection2D(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius)
    {
      //rObj_1 is particle,  and rObj_2 is condition

      const GeometryType& DE_Geom = rObj_1->GetGeometry();
      double Particle_Coord[3]         = {0.0};
      Particle_Coord[0]                = DE_Geom[0].Coordinates()[0];
      Particle_Coord[1]                = DE_Geom[0].Coordinates()[1];
      Particle_Coord[2]                = DE_Geom[0].Coordinates()[2];

      const GeometryType& FE_Geom = rObj_2->GetGeometry();

      std::vector< array_1d<double,3> > Coord(2);

      for (unsigned int i = 0; i<2; i++)
      {
        for (unsigned int j = 0; j<3; j++)
        {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }

      return GeometryFunctions::FastEdgeVertexCheck( Coord[0], Coord[1],  Particle_Coord, Radius );
              
    }//FastIntersection2D
    
    static inline bool FastIntersection3D(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius)
    {
      //rObj_1 is particle,  and rObj_2 is condition
      
      bool ContactExists = false;
      const GeometryType& DE_Geom = rObj_1->GetGeometry();
      double Particle_Coord[3]         = {0.0};
      Particle_Coord[0]                = DE_Geom[0].Coordinates()[0];
      Particle_Coord[1]                = DE_Geom[0].Coordinates()[1];
      Particle_Coord[2]                = DE_Geom[0].Coordinates()[2];

      const GeometryType& FE_Geom = rObj_2->GetGeometry();
      unsigned int FE_size = FE_Geom.size();

      std::vector< array_1d<double,3> > Coord;
      Coord.resize(FE_size, array_1d<double,3>(3,0.0) );

      for (unsigned int i = 0; i<FE_size; i++)
      {
        for (unsigned int j = 0; j<3; j++)
        {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }
      
      double distance_point_to_plane   = 0.0;
      unsigned int current_edge_index  = 0;
      
      ContactExists = GeometryFunctions::FastFacetCheck( Coord,  Particle_Coord, Radius, distance_point_to_plane, current_edge_index);
      
      if(ContactExists == true){return true;}
      
      //The key here is to see that we only need to check for further contact if, when not having contact with plane, the distance_point_to_plane is lower than the radius. 
      //In this case it might have contact with edges or vertices, otherwise no contact is possible. 

      if ( (ContactExists == false) && (distance_point_to_plane < Radius ) )
      {
        bool local_contact_exists = false;
        for (unsigned int e = current_edge_index; e < FE_size; e++ )
        {
          
          local_contact_exists = GeometryFunctions::FastEdgeVertexCheck( Coord[e], Coord[(e+1)%FE_size],  Particle_Coord, Radius );
          if(local_contact_exists) {return true;}
          
        }//for every edge

      }//no plane contact found
      
      return false;    
        
    }//FastIntersection3D
    

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius) //rObj_1 is sphere, rObj_2 is FE
    {

      int facet_size = rObj_2->GetGeometry().size();
      if (facet_size==2)
      {
         return FastIntersection2D(rObj_1, rObj_2, Radius);//, NewContactType);
      }
      else
      {
         return FastIntersection3D(rObj_1, rObj_2, Radius);//, NewContactType);
      }

    }
 
    //******************************************************************************************************************


  //needed for bins
  static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        array_1d<double, 3> center_of_particle1 = rObj_1->GetGeometry().GetPoint(0);
        array_1d<double, 3> center_of_particle2 = rObj_2->GetGeometry().GetPoint(0);


        distance = sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
                        (center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
                        (center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    }
    
    
   static inline bool DistanceHierarchy(SphericParticle* rObj_1, DEMWall* rObj_2, double LocalCoordSystem[3][3], double DistPToB, std::vector<double> Weight, int ContactType)
   {
        std::vector<double>& RF_Param = rObj_1->mNeighbourRigidFacesPram;  //NOTE: THIS member has to be cleared every time step. It is done in: the SphericParticle Function: ComputeNewRigidFaceNeighboursHistoricalData

        int new_ID = rObj_2->Id();

        std::size_t i_current = RF_Param.size();
        std::size_t neighbor_size = i_current / 16;

        bool substitute = false;
        int  position   = i_current;

        for (std::size_t i_old_neigh = 0; i_old_neigh < neighbor_size; i_old_neigh++)
        {
            int i_Old_RF_position = i_old_neigh * 16;

            double Old_Normal_Vector[3] = {0.0};
            Old_Normal_Vector[0] = RF_Param[i_Old_RF_position + 6];
            Old_Normal_Vector[1] = RF_Param[i_Old_RF_position + 7];
            Old_Normal_Vector[2] = RF_Param[i_Old_RF_position + 8];

            double New_dist = DistPToB;
            double Old_dist = RF_Param[i_Old_RF_position + 9];

           double New_projected_on_old = GeometryFunctions::DotProduct(LocalCoordSystem[2], Old_Normal_Vector);
           double New_projected_distance = New_projected_on_old * New_dist;
           double Old_projected_distance = New_projected_on_old * Old_dist;
           //MSI: manera mas simple de demostrarlo
           //if the dot product is positive, we only need to compare the distances,
           //if the dot product is negative, no one cancel no one.
/*
           if( New_projected_on_old > 0.0)
           {
               //hierarchy  
               if(Old_dist < New_dist) //old has hierarchy over new
               {
                 
               }
               else //New has hierarchy over old
               {
                 
               }
           }
           else
           {
             
               //new neighbour
           }
  */         
           if ( ( (New_projected_distance - Old_dist) / fabs(Old_dist) ) > -1.0e-15 ) //old has hierarchy over new
           {
              //DO NOT SAVE NEW NEIGH
              return false;
           }

           if ( ( (Old_projected_distance-New_dist )  / fabs(New_dist) ) > -1.0e-15 )  //new has hierarchy over old
           {
             int old_ID = RF_Param[i_Old_RF_position + 14];
             if (new_ID == old_ID) //SUBSTITUTE
             {
                position = i_Old_RF_position;
                substitute = true;
             }
             else
             {
               RF_Param[i_Old_RF_position + 15] = -1;
             }
           } //new has hierarchy over

         }//Loop over Old Neigh

        std::vector<DEMWall*>& neighbour_rigid_faces = rObj_1->mNeighbourRigidFaces;

        if(substitute)
        {

          neighbour_rigid_faces[position / 16] = rObj_2;

        } //if substitute we wont resize or pushback
        else
        {

          RF_Param.resize(position + 16);
          neighbour_rigid_faces.push_back(rObj_2);
        }

        RF_Param[position + 0]  = LocalCoordSystem[0][0];//Room for improvement: ¿could we just store LocalCoordSystem[2]? Could we regenerate 0 and 1 randomly when we need the coordinate system?
        RF_Param[position + 1]  = LocalCoordSystem[0][1];
        RF_Param[position + 2]  = LocalCoordSystem[0][2];
        RF_Param[position + 3]  = LocalCoordSystem[1][0];
        RF_Param[position + 4]  = LocalCoordSystem[1][1];
        RF_Param[position + 5]  = LocalCoordSystem[1][2];
        RF_Param[position + 6]  = LocalCoordSystem[2][0];
        RF_Param[position + 7]  = LocalCoordSystem[2][1];
        RF_Param[position + 8]  = LocalCoordSystem[2][2];
        RF_Param[position + 9]  = DistPToB;
        RF_Param[position + 10] = Weight[0];
        RF_Param[position + 11] = Weight[1];
        RF_Param[position + 12] = Weight[2];
        RF_Param[position + 13] = Weight[3];
        RF_Param[position + 14] = new_ID;
        RF_Param[position + 15] = ContactType;

        return true;
   }//DistanceHierarchy

    static inline void DoubleHierarchyMethod3D(SphericParticle* rObj_1, DEMWall* rObj_2) 
    {
          
       //rObj_1 is particle,  and rObj_2 is condition
       //TYPE HIERARCHY
       int ContactType          = -1;
       //-1: No contact;
       // 1: Plane;
       // 2: Edge;
       // 3: Point.
       // 4: Edge or Point

       bool ContactExists = false;
       const GeometryType& DE_Geom = rObj_1->GetGeometry();
       double Particle_Coord[3]         = {0.0};
       Particle_Coord[0]                = DE_Geom[0].Coordinates()[0];
       Particle_Coord[1]                = DE_Geom[0].Coordinates()[1];
       Particle_Coord[2]                = DE_Geom[0].Coordinates()[2];
       double Radius                    = DE_Geom[0].FastGetSolutionStepValue(RADIUS);

       const GeometryType& FE_Geom = rObj_2->GetGeometry();
       unsigned int FE_size = FE_Geom.size();

       double local_coord_system[3][3]  = { {0.0},{0.0},{0.0} };
       std::vector<double> Weight(4,0.0);
       std::vector< array_1d<double,3> > Coord;
       Coord.resize(FE_size, array_1d<double,3>(3,0.0) );

       for (unsigned int i = 0; i<FE_size; i++)
       {
         for (unsigned int j = 0; j<3; j++)
         {
             Coord[i][j] = FE_Geom[i].Coordinates()[j];
         }
       }
       
       double distance_point_to_plane   = 0.0;
       unsigned int current_edge_index  = 0;
       ContactExists = GeometryFunctions::FacetCheck(Coord,  Particle_Coord, Radius, local_coord_system,
                                                     distance_point_to_plane, Weight, current_edge_index);

       if (ContactExists == true)
       {
         ContactType = 1;
         ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_plane, Weight, ContactType);
         return;
       }
       //The key here is to see that we only need to check for further contact if, when not having contact with plane, the distance_point_to_plane is lower than the radius.
       //In this case it might have contact with edges or vertices, otherwise no contact is possible.

       //The check should avoid the edges which yielded a OUTSIDE value in the inside-outside check. i.e. we will check from the current edge to the last.

       ///Particle-edge contact and Particle-point
       if ( (ContactExists == false) && (distance_point_to_plane < Radius ) )
       {
          bool local_contact_exists = false;
          for (unsigned int e = current_edge_index; e < FE_size; e++ )
          {

            double eta = 0.5; // dummy initialize
            double distance_point_to_edge = 2*Radius; //dummy big initialization

            local_contact_exists = GeometryFunctions::EdgeCheck( Coord[e], Coord[(e+1)%FE_size],  Particle_Coord, Radius, local_coord_system,
                                                          distance_point_to_edge, eta);
            if (local_contact_exists)
            {
                //save data
                ContactType           = 2;
                Weight[e]             = 1-eta;
                Weight[(e+1)%FE_size] = eta; //the rest remain 0 (valid for triangles and quadrilateral)
                if(FE_size > 4){KRATOS_WATCH("WEIGHTS ALONG EDGE CANT BE CALCULATED WITH SUKUMAR FORMULAE")}
                ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_edge, Weight, ContactType);
                continue; //skip vertex check
            }
            if ( (local_contact_exists == false) && (distance_point_to_edge < Radius ) )
            {
              unsigned int vertex_to_check = -1;
              if(eta<0.0){ vertex_to_check = e;}
              else if(eta>1.0){ vertex_to_check = (e+1)%FE_size;}
              else{continue;}
              double distance_point_to_vertex = 0.0;
              local_contact_exists = GeometryFunctions::VertexCheck(Coord[vertex_to_check], Particle_Coord, Radius, local_coord_system, distance_point_to_vertex);

              if(local_contact_exists)
              {
                ContactType             = 3;
                Weight[vertex_to_check] = 1.0; //the rest weights stay 0.0;
                ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_vertex, Weight, ContactType);
              }

            } // (ContactExists == false) && (distance_point_to_edge < Radius )

          }//for every edge

       }//no plane contact found

    } //DoubleHierarchyMethod3D

    static inline void DoubleHierarchyMethod2D(SphericParticle* rObj_1, DEMWall* rObj_2)
    {
      //rObj_1 is particle,  and rObj_2 is condition

      int ContactType          = -1;
      //-1: No contact;
      // 2: Edge;
      // 3: Vertex.

      bool ContactExists = false;

      const GeometryType& DE_Geom = rObj_1->GetGeometry();
      double Particle_Coord[3]         = {0.0};
      Particle_Coord[0]                = DE_Geom[0].Coordinates()[0];
      Particle_Coord[1]                = DE_Geom[0].Coordinates()[1];
      Particle_Coord[2]                = DE_Geom[0].Coordinates()[2];
      double Radius                    = DE_Geom[0].FastGetSolutionStepValue(RADIUS);

      const GeometryType& FE_Geom = rObj_2->GetGeometry();

      double local_coord_system[3][3]  = { {0.0},{0.0},{0.0} };
      std::vector<double> Weight(2,0.0);
      std::vector< array_1d<double,3> > Coord(2);

      for (unsigned int i = 0; i<2; i++)
      {
        for (unsigned int j = 0; j<3; j++)
        {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }

      double eta = 0.5; // dummy initialize
      double distance_point_to_edge = 2*Radius; //dummy big initialization

      ContactExists = GeometryFunctions::EdgeCheck( Coord[0], Coord[1],  Particle_Coord, Radius, local_coord_system,
                                                        distance_point_to_edge, eta);
      if (ContactExists)
      {
          //save data
          ContactType           = 2;
          Weight[0]             = 1-eta;
          Weight[1]             = eta; //the rest remain 0 (valid for triangles and quadrilateral)

          ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_edge, Weight, ContactType);
          return;
      }
      if ( (ContactExists == false) && (distance_point_to_edge < Radius ) )
      {
        unsigned int vertex_to_check = -1;
        if(eta<0.0){ vertex_to_check = 0;}
        else if(eta>1.0){ vertex_to_check = 1;}
        double distance_point_to_vertex = 0.0;
        ContactExists = GeometryFunctions::VertexCheck( Coord[vertex_to_check], Particle_Coord, Radius, local_coord_system, distance_point_to_vertex);

        if(ContactExists)
        {
          ContactType             = 3;
          Weight[vertex_to_check] = 1.0; //the rest weights stay 0.0;
          ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_vertex, Weight, ContactType);
          return;
        }

      } // (ContactExists == false) && (distance_point_to_edge < Radius )

    }//DoubleHierarchyMethod2D


  

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
