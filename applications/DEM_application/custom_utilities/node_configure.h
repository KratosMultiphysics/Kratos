//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//



#if !defined(KRATOS_NODE_CONFIGURE_INCLUDED)
#define  KRATOS_NODE_CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>

// Kratos includes
#include "spatial_containers/spatial_search.h"

/* Timer defines */
#include "utilities/timer.h"

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
class NodeConfigure
{

public:
  
    enum { 
        Dimension = TDimension,
        DIMENSION = TDimension,
        MAX_LEVEL = 16,
        MIN_LEVEL = 2
    };

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(NodeConfigure);
    
    typedef SpatialSearch                                           SearchType;

    typedef SearchType::PointType                                   PointType;
    typedef SearchType::NodesContainerType::ContainerType        ContainerType;
    typedef SearchType::NodesContainerType                       NodesContainerType;  // * Comentar -> Create_and_destroy.h
    
    typedef SearchType::NodeType                                 NodeType;
    typedef ContainerType::value_type                               PointerType;
    typedef ContainerType::iterator                                 IteratorType;
    
    typedef SearchType::NodesContainerType::ContainerType        ResultContainerType;
//     typedef SearchType::ResultDistanceType::ContainerType             ResultDistanceType;
    
    typedef ResultContainerType::iterator                           ResultIteratorType;
    typedef std::vector<double>::iterator                           DistanceIteratorType;
    
    typedef ContactPair<PointerType>                                ContactPairType;
    typedef std::vector<ContactPairType>                            ContainerContactType;
    typedef ContainerContactType::iterator                          IteratorContactType;
    typedef ContainerContactType::value_type                        PointerContactType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    NodeConfigure(){};
    virtual ~NodeConfigure(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {    
//         std::cout << "A" << std::endl;
        for(std::size_t i = 0; i < 3; i++)
        {
            rHighPoint[i] = rLowPoint[i]  = (*rObject)[i];
        }
        
        double radius = rObject->FastGetSolutionStepValue(RADIUS);

        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
        }
//         std::cout << "B" << std::endl;
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
//         std::cout << "C" << std::endl;
        for(std::size_t i = 0; i < 3; i++)
        {
            rHighPoint[i] = rLowPoint[i]  = (*rObject)[i];
        }
        
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -Radius;
            rHighPoint[i] += Radius;
        }
//         std::cout << "D" << std::endl;
    }
        
    static inline void CalculateCenter(const PointerType& rObject, PointType& rCenter)
    {
//         std::cout << "E" << std::endl;
        for(std::size_t i = 0; i < 3; i++)
        {
            rCenter[i]  = (*rObject)[i];
        }
//         std::cout << "F" << std::endl;
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
//         std::cout << "G" << std::endl;
        array_1d<double, 3> rObj_2_to_rObj_1;
        
        for(std::size_t i = 0; i < 3; i++)
        {
            rObj_2_to_rObj_1[i]  = (*rObj_1)[i] - (*rObj_2)[i];
        }
        
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);

        const double& radius_1 = rObj_1->FastGetSolutionStepValue(RADIUS);
        const double& radius_2 = rObj_2->FastGetSolutionStepValue(RADIUS);
        double radius_sum      = radius_1 + radius_2;
        bool intersect         = floatle((distance_2 - radius_sum * radius_sum),0);
        
//         std::cout << "H" << std::endl;
        return intersect;
    }


    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double& Radius)
    {
//         std::cout << "I" << std::endl;
        array_1d<double, 3> rObj_2_to_rObj_1;
        
        for(std::size_t i = 0; i < 3; i++)
        {
            rObj_2_to_rObj_1[i]  = (*rObj_1)[i] - (*rObj_2)[i];
        }
        
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);

        const double& radius_1 = Radius;
        const double& radius_2 = rObj_2->FastGetSolutionStepValue(RADIUS);
        double radius_sum      = radius_1 + radius_2;
        bool intersect         = floatle((distance_2 - radius_sum * radius_sum),0);
//         std::cout << "J" << std::endl;
        return intersect;
    }

    //******************************************************************************************************************
    
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
//         std::cout << "K" << std::endl;
//        double separation_from_particle_radius_ratio = 0.1;

        array_1d<double, 3> center_of_particle;
        
        for(std::size_t i = 0; i < 3; i++)
        {
            center_of_particle[i]  = (*rObject)[i];
        }
 
        const double& radius = rObject->FastGetSolutionStepValue(RADIUS);

        bool intersect = (
          floatle(rLowPoint[0]  - radius,center_of_particle[0]) && 
          floatle(rLowPoint[1]  - radius,center_of_particle[1]) && 
          floatle(rLowPoint[2]  - radius,center_of_particle[2]) &&
          floatge(rHighPoint[0] + radius,center_of_particle[0]) && 
          floatge(rHighPoint[1] + radius,center_of_particle[1]) && 
          floatge(rHighPoint[2] + radius,center_of_particle[2]));
//         std::cout << "L" << std::endl;
        return  intersect;
 
    }

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
//       std::cout << "M" << std::endl;
//        double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle;
        
        for(std::size_t i = 0; i < 3; i++)
        {
            center_of_particle[i]  = (*rObject)[i];
        }

        double radius = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        bool intersect = (
          floatle(rLowPoint[0]  - radius,center_of_particle[0]) && 
          floatle(rLowPoint[1]  - radius,center_of_particle[1]) && 
          floatle(rLowPoint[2]  - radius,center_of_particle[2]) &&
          floatge(rHighPoint[0] + radius,center_of_particle[0]) && 
          floatge(rHighPoint[1] + radius,center_of_particle[1]) && 
          floatge(rHighPoint[2] + radius,center_of_particle[2]));
//         std::cout << "N" << std::endl;
        return  intersect;
    }

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
//         std::cout << "O" << std::endl;
        array_1d<double, 3> center_of_particle1;
        array_1d<double, 3> center_of_particle2;
                
        for(std::size_t i = 0; i < 3; i++)
        {
            center_of_particle1[i]  = (*rObj_1)[i];
            center_of_particle2[i]  = (*rObj_2)[i];
        }
// std::cout << "P" << std::endl;
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
    virtual std::string Info() const {return " Spatial Containers Configure for Particles"; }

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
    NodeConfigure& operator=(NodeConfigure const& rOther);

    /// Copy constructor.
    NodeConfigure(NodeConfigure const& rOther);

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
    inline std::istream& operator >> (std::istream& rIStream, NodeConfigure<TDimension> & rThis){
        return rIStream;
        }

    /// output stream function
    template <std::size_t TDimension>
    inline std::ostream& operator << (std::ostream& rOStream, const NodeConfigure<TDimension>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }
        
    ///@}

}   // namespace Kratos.
#endif	/* NODE_CONFIGURE_H */
