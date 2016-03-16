//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu, Guillermo Casas gcasas@cimne.upc.edu
//



#if !defined(KRATOS_DISCRETE_PARTICLE__CONFIGURE_INCLUDED)
#define  KRATOS_DISCRETE_PARTICLE__CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>

// Kratos includes
#include "spatial_containers/spatial_search.h"
#include "includes/dem_variables.h"
#include "DEM_application_variables.h"
#include "../custom_elements/spheric_particle.h"

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
class DiscreteParticleConfigure
{

public:
  
    enum { 
        Dimension = TDimension,
        DIMENSION = TDimension,
        MAX_LEVEL = 16,
        MIN_LEVEL = 2
    };

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(DiscreteParticleConfigure);
    
    typedef SpatialSearch                                           SearchType;

    typedef SearchType::PointType                                   PointType;
    typedef SearchType::ElementsContainerType::ContainerType        ContainerType;
    typedef SearchType::ElementsContainerType                       ElementsContainerType;
    typedef SearchType::NodesContainerType                          NodesContainerType;
    
    typedef SearchType::ElementType                                 ElementType;
    typedef ContainerType::value_type                               PointerType;
    typedef ContainerType::iterator                                 IteratorType;
    typedef ElementsContainerType::iterator                         ElementIteratorType;
    
    typedef SearchType::ElementsContainerType::ContainerType        ResultContainerType;
//     typedef SearchType::ResultDistanceType::ContainerType             ResultDistanceType;
    
    typedef ResultContainerType::iterator                           ResultIteratorType;
    typedef std::vector<double>::iterator                           DistanceIteratorType;
    
    typedef ContactPair<PointerType>                                ContactPairType;
    typedef std::vector<ContactPairType>                            ContainerContactType;
    typedef ContainerContactType::iterator                          IteratorContactType;
    typedef ContainerContactType::value_type                        PointerContactType;
    
    /////////////////////////////////////////////////////////////////////////////////////////// bins_dynamic_objects.h:279
    
//     typedef SpatialSearch                                                       SearchType;
// 
//     typedef SearchType::PointType                                               PointType;
//     typedef PointerVectorSet<GeometricalObject, IndexedObject>::ContainerType   ContainerType;
//     typedef PointerVectorSet<GeometricalObject, IndexedObject>                  ElementsContainerType;
//     
//     typedef SearchType::ElementType                                             ElementType;
//     typedef ContainerType::value_type                                           PointerType;
//     typedef ContainerType::iterator                                             IteratorType;
//     typedef ElementsContainerType::iterator                                     ElementIteratorType;
//     
//     typedef PointerVectorSet<GeometricalObject, IndexedObject>::ContainerType   ResultContainerType;
// //     typedef SearchType::ResultDistanceType::ContainerType             ResultDistanceType;
//     
//     typedef ResultContainerType::iterator                           ResultIteratorType;
//     typedef std::vector<double>::iterator                           DistanceIteratorType;
//     
//     typedef ContactPair<PointerType>                                ContactPairType;
//     typedef std::vector<ContactPairType>                            ContainerContactType;
//     typedef ContainerContactType::iterator                          IteratorContactType;
//     typedef ContainerContactType::value_type                        PointerContactType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    DiscreteParticleConfigure(){};
    virtual ~DiscreteParticleConfigure(){
	}

    static void SetPeriods(double domain_period_x, double domain_period_y, double domain_period_z)
    {
        mDomainPeriods[0] = domain_period_x;
        mDomainPeriods[1] = domain_period_y;
        mDomainPeriods[2] = domain_period_z;

        if (mDomainPeriods[0] >= 0 && mDomainPeriods[1] >= 0 && mDomainPeriods[2] >= 0){
            mDomainIsPeriodic = true;
        } else {
            mDomainIsPeriodic = false;
        }

    }

    static void GetPeriods(double periods[3])
    {
        periods[0] = mDomainPeriods[0];
        periods[1] = mDomainPeriods[1];
        periods[2] = mDomainPeriods[2];        
    }

    static bool GetDomainPeriodicity()
    {
        return mDomainIsPeriodic;
    }

    static inline void TransformToClosestPeriodicCoordinates(double target[3], double base_coordinates[3])
    {
        for (unsigned int i = 0; i < 3; i++){
            double incr_i = target[i] - base_coordinates[i];

            if (fabs(incr_i) > 0.5 * mDomainPeriods[i]){
                base_coordinates[i] += GetSign(incr_i) * mDomainPeriods[i];
            }
        }

    }
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    { 
        noalias(rHighPoint) = rObject->GetGeometry()[0];
        noalias(rLowPoint)  = rObject->GetGeometry()[0];
        
        SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*rObject);        
        double radius = p_particle->GetSearchRadius();        

        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
        }
    }

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
        
    static inline void CalculateCenter(const PointerType& rObject, PointType& rCenter)
    {
        rCenter  = rObject->GetGeometry()[0];
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        double rObj_2_to_rObj_1[3];
        PeriodicSubtract(rObj_1->GetGeometry()[0], rObj_2->GetGeometry()[0], rObj_2_to_rObj_1);
        
        double distance_2 = DEM_INNER_PRODUCT_3(rObj_2_to_rObj_1, rObj_2_to_rObj_1);

        SphericParticle* p_particle1 = dynamic_cast<SphericParticle*>(&*rObj_1);
        SphericParticle* p_particle2 = dynamic_cast<SphericParticle*>(&*rObj_2);
        double radius_sum      = p_particle1->GetSearchRadius() + p_particle2->GetSearchRadius();
        bool intersect         = floatle((distance_2 - radius_sum * radius_sum),0);
        return intersect;
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double& radius_1)
    {
        double rObj_2_to_rObj_1[3];
        PeriodicSubtract(rObj_1->GetGeometry()[0], rObj_2->GetGeometry()[0], rObj_2_to_rObj_1);

        double distance_2 = DEM_INNER_PRODUCT_3(rObj_2_to_rObj_1, rObj_2_to_rObj_1);

        SphericParticle* p_particle1 = dynamic_cast<SphericParticle*>(&*rObj_1);
        SphericParticle* p_particle2 = dynamic_cast<SphericParticle*>(&*rObj_2);
        double radius_sum      = p_particle1->GetSearchRadius() + p_particle2->GetSearchRadius();
        bool intersect         = floatle((distance_2 - radius_sum * radius_sum),0);
        return intersect;
    }

    //******************************************************************************************************************
    
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        const array_1d<double, 3>& center_of_particle = rObject->GetGeometry()[0];
 
        SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*rObject);
        const double& radius = p_particle->GetSearchRadius();

        bool intersect = (
          floatle(rLowPoint[0]  - radius,center_of_particle[0]) && 
          floatle(rLowPoint[1]  - radius,center_of_particle[1]) && 
          floatle(rLowPoint[2]  - radius,center_of_particle[2]) &&
          floatge(rHighPoint[0] + radius,center_of_particle[0]) && 
          floatge(rHighPoint[1] + radius,center_of_particle[1]) && 
          floatge(rHighPoint[2] + radius,center_of_particle[2]));

        if ((mDomainPeriods[0] > 0.0) && !intersect){
            double closest_representative_coor[3] = {center_of_particle[0], center_of_particle[1], center_of_particle[2]};
            double box_center[3] = {0.5 * (rLowPoint[0] + rHighPoint[0]), 0.5 * (rLowPoint[1] + rHighPoint[1]), 0.5 * (rLowPoint[2] + rHighPoint[2])};
            TransformToClosestPeriodicCoordinates(box_center, closest_representative_coor);
            intersect = (
            floatle(rLowPoint[0]  - radius,closest_representative_coor[0]) &&
            floatle(rLowPoint[1]  - radius,closest_representative_coor[1]) &&
            floatle(rLowPoint[2]  - radius,closest_representative_coor[2]) &&
            floatge(rHighPoint[0] + radius,closest_representative_coor[0]) &&
            floatge(rHighPoint[1] + radius,closest_representative_coor[1]) &&
            floatge(rHighPoint[2] + radius,closest_representative_coor[2]));
            }

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

        if ((mDomainPeriods[0] > 0.0) && !intersect){
            double closest_representative_coor[3] = {center_of_particle[0], center_of_particle[1], center_of_particle[2]};
            double box_center[3] = {0.5 * (rLowPoint[0] + rHighPoint[0]), 0.5 * (rLowPoint[1] + rHighPoint[1]), 0.5 * (rLowPoint[2] + rHighPoint[2])};
            TransformToClosestPeriodicCoordinates(box_center, closest_representative_coor);
            intersect = (
            floatle(rLowPoint[0]  - radius,closest_representative_coor[0]) &&
            floatle(rLowPoint[1]  - radius,closest_representative_coor[1]) &&
            floatle(rLowPoint[2]  - radius,closest_representative_coor[2]) &&
            floatge(rHighPoint[0] + radius,closest_representative_coor[0]) &&
            floatge(rHighPoint[1] + radius,closest_representative_coor[1]) &&
            floatge(rHighPoint[2] + radius,closest_representative_coor[2]));
            }

        return  intersect;
    }

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        double rObj_2_to_rObj_1[3];
        PeriodicSubtract(rObj_1->GetGeometry()[0], rObj_2->GetGeometry()[0], rObj_2_to_rObj_1);
        distance = DEM_MODULUS_3(rObj_2_to_rObj_1);
    }

    static double mDomainPeriods[3];
    static bool mDomainIsPeriodic;

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

    static inline int GetSign(const double value)
    {
        return (0.0 < value) - (value < 0.0);
    }

    static inline void PeriodicSubtract(const array_1d<double, 3>& a, const array_1d<double, 3>& b, double c[3])
    {
        for (unsigned int i = 0; i < 3; i++){
            c[i] = a[i] - b[i];
        }

        if (mDomainPeriods[0] > 0.0){ // Periods have been set (the domain is periodic)
            for (unsigned int i = 0; i < 3; i++){
                if (fabs(c[i]) > 0.5 * mDomainPeriods[i]) c[i] -= GetSign(c[i]) * mDomainPeriods[i];
            }
        }
    }

    static inline bool floateq(double a, double b) {
        return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
    }
    
    static inline bool floatle(double a, double b) {
        return a < b || std::fabs(a - b) < std::numeric_limits<double>::epsilon();
    }
    
    static inline bool floatge(double a, double b) {
        return a > b || std::fabs(a - b) < std::numeric_limits<double>::epsilon();
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
    DiscreteParticleConfigure& operator=(DiscreteParticleConfigure const& rOther);

    /// Copy constructor.
    DiscreteParticleConfigure(DiscreteParticleConfigure const& rOther);

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
    inline std::istream& operator >> (std::istream& rIStream, DiscreteParticleConfigure<TDimension> & rThis){
        return rIStream;
        }

    /// output stream function
    template <std::size_t TDimension>
    inline std::ostream& operator << (std::ostream& rOStream, const DiscreteParticleConfigure<TDimension>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }

    ///@}
    ///
template <std::size_t TDimension>
double DiscreteParticleConfigure<TDimension>::mDomainPeriods[] = {-1.0, -1.0, -1.0};
template <std::size_t TDimension>
bool DiscreteParticleConfigure<TDimension>::mDomainIsPeriodic = false;

}   // namespace Kratos.
#endif	/* DISCRETE_PARTICLE_CONFIGURE_H */
