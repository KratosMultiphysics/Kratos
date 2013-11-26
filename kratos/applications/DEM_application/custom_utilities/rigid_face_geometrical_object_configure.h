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
			Coord = rObject->GetGeometry()(inode)->Coordinates();

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
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2,  const double& Radius)
     {
         //Cfeng: rObj_1 is particle,  and rObj_2 is block

         bool If_PB_Contact = false;
		 
	     double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
		 double Weight[4]              = {0.0, 0.0, 0.0, 0.0};
		 double DistPToB               = 0.0;

        double Particle_Coord[3] = {0.0};
        Particle_Coord[0] = rObj_1->GetGeometry()(0)->Coordinates()[0];
        Particle_Coord[1] = rObj_1->GetGeometry()(0)->Coordinates()[1];
        Particle_Coord[2] = rObj_1->GetGeometry()(0)->Coordinates()[2];
		
		
       // double rad = rObj_1->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
	    double rad = Radius;

        double Centroid[3] = {0.0};
        for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++)
        {
            Centroid[0] += rObj_2->GetGeometry()(inode)->Coordinates()[0] / (double)rObj_2->GetGeometry().size();
            Centroid[1] += rObj_2->GetGeometry()(inode)->Coordinates()[1] / (double)rObj_2->GetGeometry().size();
            Centroid[2] += rObj_2->GetGeometry()(inode)->Coordinates()[2] / (double)rObj_2->GetGeometry().size();
        }

		 
         if(rObj_2->GetGeometry().size() == 2)
         {
			double Coord1[3]     = {0.0};
			double Coord2[3]     = {0.0};
		    double tempWeight[2] = {0.0};
			
			Coord1[0] = rObj_2->GetGeometry()(0)->Coordinates()[0];
			Coord1[1] = rObj_2->GetGeometry()(0)->Coordinates()[1];
			Coord1[2] = rObj_2->GetGeometry()(0)->Coordinates()[2];

			Coord2[0] = rObj_2->GetGeometry()(1)->Coordinates()[0];
			Coord2[1] = rObj_2->GetGeometry()(1)->Coordinates()[1];
			Coord2[2] = rObj_2->GetGeometry()(1)->Coordinates()[2];

			//If_PB_Contact = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Particle_Coord, rad);
			
			If_PB_Contact = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Centroid, Particle_Coord, rad,
		                                                                            LocalCoordSystem,  tempWeight, DistPToB);
			
			if(If_PB_Contact == true)
			{
				Weight[0] = tempWeight[0];
				Weight[1] = tempWeight[1];			
			}
			else /// particle point contact search
			{
				for(std::size_t inode1 = 0; inode1 < rObj_2->GetGeometry().size(); inode1++)
				{
					Coord1[0] = rObj_2->GetGeometry()(inode1)->Coordinates()[0];
					Coord1[1] = rObj_2->GetGeometry()(inode1)->Coordinates()[1];
					Coord1[2] = rObj_2->GetGeometry()(inode1)->Coordinates()[2];
					
					If_PB_Contact = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, LocalCoordSystem, DistPToB);
					if(If_PB_Contact == true)
					{
						Weight[inode1] = 1.0;						
						break;
					}
				}				
			}
			
			/////////////////////////////////////
			   
         }        
         else if(rObj_2->GetGeometry().size() > 2)
         {
            double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };

			////Cfeng: Triangle
			int FaceNodeTotal = 3;

			Coord[0][0] = rObj_2->GetGeometry()(0)->Coordinates()[0];
			Coord[0][1] = rObj_2->GetGeometry()(0)->Coordinates()[1];
			Coord[0][2] = rObj_2->GetGeometry()(0)->Coordinates()[2];

			Coord[1][0] = rObj_2->GetGeometry()(1)->Coordinates()[0];
			Coord[1][1] = rObj_2->GetGeometry()(1)->Coordinates()[1];
			Coord[1][2] = rObj_2->GetGeometry()(1)->Coordinates()[2];

			Coord[2][0] = rObj_2->GetGeometry()(2)->Coordinates()[0];
			Coord[2][1] = rObj_2->GetGeometry()(2)->Coordinates()[1];
			Coord[2][2] = rObj_2->GetGeometry()(2)->Coordinates()[2];

			if(rObj_2->GetGeometry().size() == 4)
			{
				Coord[3][0] = rObj_2->GetGeometry()(3)->Coordinates()[0];
				Coord[3][1] = rObj_2->GetGeometry()(3)->Coordinates()[1];
				Coord[3][2] = rObj_2->GetGeometry()(3)->Coordinates()[2];

			  ////Cfeng: Quadral
			  FaceNodeTotal = 4;
			}
		
	        /////Particle-Face contact
		    If_PB_Contact = GeometryFunctions::JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal,Coord,  Centroid, Particle_Coord, rad,
													                           LocalCoordSystem, Weight,  DistPToB);
			///Particle-edge contact
			if(If_PB_Contact == false)
			{
				Weight[0] = Weight[1] = Weight[2] = Weight[3] = 0.0;
				
				for(int inode1 = 0; inode1 < FaceNodeTotal - 1; inode1++)
				{
					int inode2 = inode1 + 1;
					
					double Coord1[3]     = {0.0};
					double Coord2[3]     = {0.0};
					double tempWeight[2] = {0.0};
					
					Coord1[0] = rObj_2->GetGeometry()(inode1)->Coordinates()[0];
					Coord1[1] = rObj_2->GetGeometry()(inode1)->Coordinates()[1];
					Coord1[2] = rObj_2->GetGeometry()(inode1)->Coordinates()[2];

					Coord2[0] = rObj_2->GetGeometry()(inode2)->Coordinates()[0];
					Coord2[1] = rObj_2->GetGeometry()(inode2)->Coordinates()[1];
					Coord2[2] = rObj_2->GetGeometry()(inode2)->Coordinates()[2];
					
					If_PB_Contact = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Centroid, Particle_Coord, rad,
																							LocalCoordSystem,  tempWeight, DistPToB);
					
					if(If_PB_Contact == true)
					{
						int inode3 = (inode2 + 1) % FaceNodeTotal;
						Weight[inode1] = tempWeight[0];
						Weight[inode2] = tempWeight[1];
						Weight[inode3] = 0.0;
						
						if(FaceNodeTotal == 4)
						{
							int inode4 = (inode2 + 2) % FaceNodeTotal;
							Weight[inode4] = 0.0;
						}

						break;
					}
				}
			}
			///////////////////////////////////////////
			
			/////Particle-Point contact
			if(If_PB_Contact == false)
			{
				Weight[0] = Weight[1] = Weight[2] = Weight[3] = 0.0;
				
				for(int inode1 = 0; inode1 < FaceNodeTotal; inode1++)
				{
					double Coord1[3] = {0.0};
					Coord1[0] = rObj_2->GetGeometry()(inode1)->Coordinates()[0];
					Coord1[1] = rObj_2->GetGeometry()(inode1)->Coordinates()[1];
					Coord1[2] = rObj_2->GetGeometry()(inode1)->Coordinates()[2];
					
					If_PB_Contact = GeometryFunctions::JudgeIfThisPointIsContactWithParticle(Coord1, Particle_Coord, rad, LocalCoordSystem, DistPToB);
					if(If_PB_Contact == true)
					{
						Weight[inode1] = 1.0;						
						break;
					}
				}				
			}
			/////////////////////////////////////
         }
		 
		 ////Store the neighbour value, real contact rigidFace
		 if(If_PB_Contact == true)
		 {
			// In geometrical level to search the neighbour, need pointer convert
			Vector & RF_Pram= (boost::dynamic_pointer_cast<Element>(rObj_1))->GetValue(NEIGHBOUR_RIGID_FACES_PRAM);
			  
			 //Vector & RF_Pram= rObj_1->GetValue(NEIGHBOUR_RIGID_FACES_PRAM);
			 std::size_t ino = RF_Pram.size();
			 
			 std::size_t TotalSize = ino / 15;
			 
			 std::size_t isize = 0;
			 
			 for(isize = 0; isize < TotalSize; isize++)
			 {
				 if( static_cast<int> (RF_Pram[isize + 14]) == static_cast<int>(rObj_2->Id()) )
				 {
					 break;
				 }
			 }
			 
			 if(isize == TotalSize)
			 {
				 RF_Pram.resize(ino + 15);
				 RF_Pram[ino + 0]  = LocalCoordSystem[0][0];
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
			 }
			
		 }
		 
        return If_PB_Contact;
      }

    //******************************************************************************************************************
      // it is really used, checked by Cfeng,131004
     static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
     {
         //Cfeng: rObject is block

		
		array_1d<double, 3> Coord;

		double xyz_min[3] = { 1e20,  1e20,  1e20};
		double xyz_max[3] = {-1e20, -1e20, -1e20};
                
		for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
		{
			Coord = rObject->GetGeometry()(inode)->Coordinates();

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
     
    //******************************************************************************************************************

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
#endif	/* rigid_face_geometrical_object_configure.h */
