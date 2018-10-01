//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
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
        noalias(rHighPoint) = rObject->GetGeometry()[0];
        noalias(rLowPoint)  = rObject->GetGeometry()[0];

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



        double xyz_min[3] = { 1e20,  1e20,  1e20};
        double xyz_max[3] = {-1e20, -1e20, -1e20};

        for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
        {
                const array_1d<double, 3>& Coord = rObject->GetGeometry()[inode].Coordinates();

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


    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint) {

        double xyz_min[3] = { 1e20,  1e20,  1e20};
        double xyz_max[3] = {-1e20, -1e20, -1e20};

        for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++) {
          const array_1d<double, 3>& Coord = rObject->GetGeometry()[inode].Coordinates();
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
        const array_1d<double, 3>& center_of_particle = rObject->GetGeometry()[0];
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

    static inline bool FastIntersection2D(const GeometryType& DE_Geom, const GeometryType& FE_Geom,  const double& Radius)
    {
      //rObj_1 is particle,  and rObj_2 is condition

      std::vector< array_1d<double,3> > Coord(2);

      for (unsigned int i = 0; i<2; i++) {
        for (unsigned int j = 0; j<3; j++) {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }

      return GeometryFunctions::FastEdgeVertexCheck( Coord[0], Coord[1],  DE_Geom[0].Coordinates(), Radius );
    }//FastIntersection2D

    static inline bool FastIntersection3D(const GeometryType& DE_Geom,const GeometryType& FE_Geom,  const double& Radius)
    {
      //rObj_1 is particle,  and rObj_2 is condition
      bool ContactExists = false;

      unsigned int FE_size = FE_Geom.size();
      std::vector< array_1d<double,3> > Coord;
      Coord.resize(FE_size, array_1d<double,3>(3,0.0) );

      for (unsigned int i = 0; i<FE_size; i++) {
        for (unsigned int j = 0; j<3; j++) {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }

      double distance_point_to_plane   = 0.0;
      unsigned int current_edge_index  = 0;

      ContactExists = GeometryFunctions::FastFacetCheck( Coord,  DE_Geom[0].Coordinates(), Radius, distance_point_to_plane, current_edge_index);

      if (ContactExists){return true;}

      //The key here is to see that we only need to check for further contact if, when not having contact with plane, the distance_point_to_plane is lower than the radius.
      //In this case it might have contact with edges or vertices, otherwise no contact is possible.

      else if (distance_point_to_plane < Radius) {
        bool local_contact_exists = false;
        for (unsigned int e = current_edge_index; e < FE_size; e++ ) {
          local_contact_exists = GeometryFunctions::FastEdgeVertexCheck( Coord[e], Coord[(e+1)%FE_size],  DE_Geom[0].Coordinates(), Radius );
          if (local_contact_exists) {return true;}
        }//for every edge
       }//no plane contact found

      return false;
    }//FastIntersection3D

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2) //rObj_1 is sphere, rObj_2 is FE
    {
        const GeometryType& DE_Geom = rObj_1->GetGeometry();
        const GeometryType& FE_Geom = rObj_2->GetGeometry();
        SphericParticle* p_particle = static_cast<SphericParticle*>(&*rObj_2);
        const double Radius = p_particle->GetSearchRadius();

        int facet_size = FE_Geom.WorkingSpaceDimension();

        if (facet_size==2) {
           return FastIntersection2D(DE_Geom, FE_Geom, Radius);//, NewContactType);
        }
        else {
           return FastIntersection3D(DE_Geom, FE_Geom, Radius);//, NewContactType);
        }
      }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius)
    {
      const GeometryType& DE_Geom = rObj_1->GetGeometry();
      const GeometryType& FE_Geom = rObj_2->GetGeometry();

      int facet_size = FE_Geom.WorkingSpaceDimension();

      if (facet_size==2) {
         return FastIntersection2D(DE_Geom, FE_Geom, Radius);//, NewContactType);
      }
      else {
         return FastIntersection3D(DE_Geom, FE_Geom, Radius);//, NewContactType);
      }
    }

    //Copy of the Intersection function accessible with geometry
    static inline bool FastIntersection(const GeometryType& DE_Geom, const GeometryType& FE_Geom,  const double& Radius) { //rObj_1 is sphere, rObj_2 is FE

      int facet_size = FE_Geom.WorkingSpaceDimension();

      if (facet_size==2) {
         return FastIntersection2D(DE_Geom, FE_Geom, Radius);//, NewContactType);
      }
      else {
         return FastIntersection3D(DE_Geom, FE_Geom, Radius);//, NewContactType);
      }
    }



  //needed for bins
    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance) {
        const array_1d<double, 3>& center_of_particle1 = rObj_1->GetGeometry()[0];
        const array_1d<double, 3>& center_of_particle2 = rObj_2->GetGeometry()[0];

        distance = std::sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
                             (center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
                             (center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    }


   static inline bool DistanceHierarchy(SphericParticle* rObj_1,
                                        DEMWall* rObj_2,
                                        const double LocalCoordSystem[3][3],
                                        const double DistPToB,
                                        std::vector<double> Weight,
                                        int ContactType,
                                        std::vector< double > & Distance_Array,
                                        std::vector< array_1d<double,3> >& Normal_Array,
                                        std::vector< array_1d<double,4> >& Weight_Array,
                                        std::vector<int> & Id_Array,
                                        std::vector<int> & ContactTypes
                                       )
   {

        int new_ID = rObj_2->Id();
        std::size_t i_current = Normal_Array.size();
        std::size_t neighbor_size = Normal_Array.size();

        bool substitute = false;
        int  position   = i_current;

        for (std::size_t i_old_neigh = 0; i_old_neigh < neighbor_size; i_old_neigh++)
        {
            //int i_Old_RF_position = i_old_neigh * 16;

            double Old_Normal_Vector[3] = {0.0};
            Old_Normal_Vector[0] = Normal_Array[i_old_neigh][0];
            Old_Normal_Vector[1] = Normal_Array[i_old_neigh][1];
            Old_Normal_Vector[2] = Normal_Array[i_old_neigh][2];

            double New_Dist = DistPToB;
            double Old_dist = Distance_Array[i_old_neigh];

           double New_projected_on_old = DEM_INNER_PRODUCT_3(LocalCoordSystem[2], Old_Normal_Vector);
           double New_projected_distance = New_projected_on_old * New_Dist;
           double Old_projected_distance = New_projected_on_old * Old_dist;

           if (New_projected_distance - Old_dist > -1.0e-15 * std::abs(Old_dist)) {//old has hierarchy over new  //DO NOT SAVE NEW NEIGH
             return false;
           }

           if (Old_projected_distance - New_Dist > -1.0e-15 * std::abs(New_Dist)) { //new has hierarchy over old

             int old_ID = Id_Array[i_old_neigh];
             if (new_ID == old_ID) {//SUBSTITUTE
                position = i_old_neigh;
                substitute = true;
             }

             else {  //DISABLE THE OLD ONE
               ContactTypes[i_old_neigh] = -1;
             }

           } //new has hierarchy over

         }//Loop over Old Neigh

        std::vector<DEMWall*>& neighbour_rigid_faces = rObj_1->mNeighbourRigidFaces;

        if(!substitute) { //if substitute we wont resize or pushback
          Distance_Array.resize(neighbor_size+1);
          Weight_Array.resize(neighbor_size+1);
          Normal_Array.resize(neighbor_size+1);
          Id_Array.resize(neighbor_size+1);
          ContactTypes.resize(neighbor_size+1);

          neighbour_rigid_faces.push_back(rObj_2);
        }

        Normal_Array[position][0]  = LocalCoordSystem[2][0];
        Normal_Array[position][1]  = LocalCoordSystem[2][1];
        Normal_Array[position][2]  = LocalCoordSystem[2][2];
        Weight_Array[position][0] = Weight[0];
        Weight_Array[position][1] = Weight[1];
        Weight_Array[position][2] = Weight[2];
        Weight_Array[position][3] = Weight[3];
        Distance_Array[position]  = DistPToB;

        Id_Array[position] = new_ID;
        ContactTypes[position] = ContactType;


        return true;
   }//DistanceHierarchy


      static inline void DoubleHierarchyMethod(SphericParticle* rObj_1,
                                               DEMWall* rObj_2,
                                               std::vector< double > & Distance_Array,
                                               std::vector< array_1d<double,3> >& Normal_Array,
                                               std::vector< array_1d<double,4> >& Weight_Array,
                                               std::vector< int >& Id_Array,
                                               std::vector<int> & ContactType_Array
                                               )
      {
        const GeometryType& FE_Geom = rObj_2->GetGeometry();
        unsigned int FE_size = FE_Geom.size();

        if(FE_size==2) {
          DoubleHierarchyMethod2D(rObj_1,rObj_2, Distance_Array, Normal_Array,Weight_Array, Id_Array,ContactType_Array);
        }
        else {
          DoubleHierarchyMethod3D(rObj_1,rObj_2, Distance_Array, Normal_Array,Weight_Array, Id_Array,ContactType_Array);
        }

     }

    static inline void DoubleHierarchyMethod3D(SphericParticle* rObj_1,
                                               DEMWall* rObj_2,
                                               std::vector< double > & Distance_Array,
                                               std::vector< array_1d<double,3> >& Normal_Array,
                                               std::vector< array_1d<double,4> >& Weight_Array,
                                               std::vector< int >& Id_Array,
                                               std::vector< int >& ContactType_Array
                                              )

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

       double Radius                    = rObj_1->GetInteractionRadius();

       const GeometryType& FE_Geom = rObj_2->GetGeometry();
       unsigned int FE_size = FE_Geom.size();

       double local_coord_system[3][3]  = { {0.0},{0.0},{0.0} };
       std::vector<double> Weight(4,0.0);

       double distance_point_to_plane   = 0.0;
       unsigned int current_edge_index  = 0;
       ContactExists = GeometryFunctions::FacetCheck(FE_Geom, DE_Geom[0].Coordinates(), Radius, local_coord_system,
                                                     distance_point_to_plane, Weight, current_edge_index);
       if (ContactExists == true) {
         ContactType = 1;
         ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_plane, Weight, ContactType, Distance_Array, Normal_Array,Weight_Array, Id_Array, ContactType_Array);
         return;
       }
       //The key here is to see that we only need to check for further contact if, when not having contact with plane, the distance_point_to_plane is lower than the radius.
       //In this case it might have contact with edges or vertices, otherwise no contact is possible.

       //The check should avoid the edges which yielded a OUTSIDE value in the inside-outside check. i.e. we will check from the current edge to the last.

       ///Particle-edge contact and Particle-point
       if ( (ContactExists == false) && (distance_point_to_plane < Radius ) ) {

          bool local_contact_exists = false;
          for (unsigned int e = current_edge_index; e < FE_size; e++ ) {
            double eta = 0.5; // dummy initialize
            double distance_point_to_edge = 2.0 * Radius; //dummy big initialization

            local_contact_exists = GeometryFunctions::EdgeCheck( FE_Geom[e], FE_Geom[(e+1)%FE_size],  DE_Geom[0].Coordinates(), Radius, local_coord_system,
                                                          distance_point_to_edge, eta);
            if (local_contact_exists) {
                //save data
                ContactType           = 2;
                Weight[e]             = 1-eta;
                Weight[(e+1)%FE_size] = eta; //the rest remain 0 (valid for triangles and quadrilateral)
                if(FE_size > 4){KRATOS_WATCH("WEIGHTS ALONG EDGE CANT BE CALCULATED WITH SUKUMAR FORMULAE")}
                ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_edge, Weight, ContactType, Distance_Array, Normal_Array,Weight_Array, Id_Array, ContactType_Array);
                continue; //skip vertex check
            }
            if ( (local_contact_exists == false) && (distance_point_to_edge < Radius ) ) {
              unsigned int vertex_to_check = -1;
              if(eta<0.0) {vertex_to_check = e;}
              else if(eta>1.0){ vertex_to_check = (e+1)%FE_size;}
              else {continue;}
              double distance_point_to_vertex = 0.0;
              local_contact_exists = GeometryFunctions::VertexCheck(FE_Geom[vertex_to_check], DE_Geom[0].Coordinates(), Radius, local_coord_system, distance_point_to_vertex);

              if(local_contact_exists) {
                ContactType             = 3;
                Weight[vertex_to_check] = 1.0; //the rest weights stay 0.0;
                ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_vertex, Weight, ContactType, Distance_Array, Normal_Array,Weight_Array, Id_Array, ContactType_Array);
              }
            } // (ContactExists == false) && (distance_point_to_edge < Radius )
          }//for every edge
       }//no plane contact found
       return;
    } //DoubleHierarchyMethod3D


    static inline void DoubleHierarchyMethod2D(SphericParticle* rObj_1,
                                               DEMWall* rObj_2,
                                               std::vector< double > & Distance_Array,
                                               std::vector< array_1d<double,3> >& Normal_Array,
                                               std::vector< array_1d<double,4> >& Weight_Array,
                                               std::vector< int > & Id_Array,
                                               std::vector< int > & ContactType_Array
                                              )
      {
      //rObj_1 is particle,  and rObj_2 is condition
      int ContactType          = -1;
      //-1: No contact;
      // 2: Edge;
      // 3: Vertex.

      bool ContactExists = false;

      const GeometryType& DE_Geom = rObj_1->GetGeometry();

      double Radius                    = rObj_1->GetInteractionRadius();

      const GeometryType& FE_Geom = rObj_2->GetGeometry();

      double local_coord_system[3][3]  = { {0.0},{0.0},{0.0} };
      std::vector<double> Weight(2,0.0);
      std::vector< array_1d<double,3> > Coord(2);

      for (unsigned int i = 0; i<2; i++) {
        for (unsigned int j = 0; j<3; j++) {
            Coord[i][j] = FE_Geom[i].Coordinates()[j];
        }
      }

      double eta = 0.5; // dummy initialize
      double distance_point_to_edge = 2*Radius; //dummy big initialization

      ContactExists = GeometryFunctions::EdgeCheck( Coord[0], Coord[1],  DE_Geom[0].Coordinates(), Radius, local_coord_system,
                                                        distance_point_to_edge, eta);
      if (ContactExists) { //save data
          ContactType           = 2;
          Weight[0]             = 1-eta;
          Weight[1]             = eta; //the rest remain 0 (valid for triangles and quadrilateral)

          ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_edge, Weight, ContactType, Distance_Array, Normal_Array, Weight_Array, Id_Array, ContactType_Array);
          return;
      }

      if ( (ContactExists == false) && (distance_point_to_edge < Radius ) ) {
        unsigned int vertex_to_check = -1;
        if(eta<0.0){ vertex_to_check = 0;}
        else if(eta>1.0){ vertex_to_check = 1;}
        double distance_point_to_vertex = 0.0;
        ContactExists = GeometryFunctions::VertexCheck( Coord[vertex_to_check], DE_Geom[0].Coordinates(), Radius, local_coord_system, distance_point_to_vertex);

        if(ContactExists) {
          ContactType             = 3;
          Weight[vertex_to_check] = 1.0; //the rest weights stay 0.0;
          ContactExists = DistanceHierarchy(rObj_1,rObj_2, local_coord_system, distance_point_to_vertex, Weight, ContactType, Distance_Array, Normal_Array, Weight_Array, Id_Array, ContactType_Array);
          return;
        }
      } // (ContactExists == false) && (distance_point_to_edge < Radius )
      return;
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
