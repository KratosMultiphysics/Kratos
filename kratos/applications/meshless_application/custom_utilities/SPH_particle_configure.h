//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine  $
//   Date:                $Date: 2012-05-24 $
//   Revision:            $Revision: 1.0 $
//


#if !defined(KRATOS_SPH_PARTICLE__CONFIGURE_INCLUDED)
#define  KRATOS_SPH_PARTICLE__CONFIGURE_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include "utilities/spatial_containers_configure.h"

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
class SPHParticleConfigure{
public:

    enum { Dimension = TDimension,
           DIMENSION = TDimension,
           MAX_LEVEL = 16,
           MIN_LEVEL = 2
         };

    typedef  Point<Dimension, double>                        PointType;
    typedef  std::vector<double>::iterator                   DistanceIteratorType;
    typedef  std::vector<typename Element::Pointer>          ContainerType; // ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef  ContainerType::value_type                       PointerType;
    typedef  ContainerType::iterator                         IteratorType;
    typedef  ModelPart::ElementsContainerType                ElementsContainerType;     // pointer vector set de elements.
    typedef  ElementsContainerType::ContainerType            ResultContainerType;       // std vector de punters a elements
    typedef  ResultContainerType::value_type                 ResultPointerType;
    typedef  ResultContainerType::iterator                   ResultIteratorType;
    typedef  ContactPair<PointerType>                        ContactPairType;
    typedef  std::vector<ContactPairType>                    ContainerContactType;
    typedef  ContainerContactType::iterator                  IteratorContactType;
    typedef  ContainerContactType::value_type                PointerContactType;

    typedef  std::vector<PointerType>::iterator              PointerTypeIterator;






    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(SPHParticleConfigure);

    ///@}
    ///@name Life Cycle
    ///@{

    SPHParticleConfigure(){};
    virtual ~SPHParticleConfigure(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint){
        
        rHighPoint = rLowPoint  = rObject->GetGeometry().GetPoint(0);


//        double radius = rObject->GetGeometry()(0)->GetSolutionStepValue(RADIUS);

//        for(std::size_t i = 0; i < 3; i++){
//            rLowPoint[i]  += -radius;
//            rHighPoint[i] += radius;
//        }
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius){

        rHighPoint = rLowPoint  = rObject->GetGeometry().GetPoint(0);
        
        for(std::size_t i = 0; i < 3; i++){
            rLowPoint[i]  += -Radius;
            rHighPoint[i] += Radius;
        }

    }

    static inline void CalculateCenter(const PointerType& rObject, PointType& rCenter){

        rCenter  = rObject->GetGeometry().GetPoint(0);
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2){

        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetGeometry().GetPoint(0) - rObj_2->GetGeometry().GetPoint(0);
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        const double& radius_1 = rObj_1->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
        const double& radius_2 = rObj_2->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
        double radius_sum      = radius_1 + radius_2;
        bool intersect         = (distance_2 - radius_sum * radius_sum) <= 0;
        return intersect;
    }


    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double& Radius){
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetGeometry().GetPoint(0) - rObj_2->GetGeometry().GetPoint(0);
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        //         double radius_1 = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        //         double radius_1 = Radius;
        const double& radius_1 = Radius;
        const double& radius_2 = rObj_2->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
        double radius_sum = radius_1 + radius_2;
        bool intersect = (distance_2 - radius_sum * radius_sum) <= 0;
        return intersect;
    }



    //******************************************************************************************************************
    
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {

        //        double separation_from_particle_radius_ratio = 0.1;

        array_1d<double, 3> center_of_particle = rObject->GetGeometry().GetPoint(0);

        const double& radius = rObject->GetGeometry()(0)->GetSolutionStepValue(RADIUS);

        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
                          rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);

        return  intersect;

    }

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
        //        double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle = rObject->GetGeometry().GetPoint(0);

        double radius = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&

                          rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
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
    SPHParticleConfigure& operator=(SPHParticleConfigure const& rOther);

    /// Copy constructor.
    SPHParticleConfigure(SPHParticleConfigure const& rOther);

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
inline std::istream& operator >> (std::istream& rIStream, SPHParticleConfigure<TDimension> & rThis){
    return rIStream;
}

/// output stream function
template <std::size_t TDimension>
inline std::ostream& operator << (std::ostream& rOStream, const SPHParticleConfigure<TDimension>& rThis){
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}   // namespace Kratos.
#endif	/* SPH_PARTICLE_CONFIGURE_H */
