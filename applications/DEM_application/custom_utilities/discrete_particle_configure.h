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

    DiscreteParticleConfigure(){};
    virtual ~DiscreteParticleConfigure(){
	}

    static void SetDomain(const double domain_min_x, const double domain_min_y, const double domain_min_z,
                          const double domain_max_x, const double domain_max_y, const double domain_max_z)
    {
        mDomainMin[0] = domain_min_x;
        mDomainMin[1] = domain_min_y;
        mDomainMin[2] = domain_min_z;
        mDomainMax[0] = domain_max_x;
        mDomainMax[1] = domain_max_y;
        mDomainMax[2] = domain_max_z;
        SetPeriods(domain_max_x - domain_min_x, domain_max_y - domain_min_y, domain_max_z - domain_min_z);
        mDomainIsPeriodic = (mDomainPeriods[0] >= 0 && mDomainPeriods[1] >= 0 && mDomainPeriods[2] >= 0);

    }

    static void SetPeriods(double domain_period_x, double domain_period_y, double domain_period_z)
    {
        mDomainPeriods[0] = domain_period_x;
        mDomainPeriods[1] = domain_period_y;
        mDomainPeriods[2] = domain_period_z;
    }

    static double* GetMinPoint()
    {
        return mDomainMin;
    }

    static double* GetMaxPoint()
    {
        return mDomainMax;
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

    static inline void TransformToClosestPeriodicCoordinates(const double target[3], double base_coordinates[3])
    {
        TransformToClosestPeriodicCoordinates(target, base_coordinates, mDomainPeriods);
    }

    static inline void TransformToClosestPeriodicCoordinates(const array_1d<double,3>& target, array_1d<double,3>& base_coordinates)
    {
        TransformToClosestPeriodicCoordinates(target, base_coordinates, mDomainPeriods);
    }

    static inline void TransformToClosestPeriodicCoordinates(const double target[3], double base_coordinates[3], const double periods[3])
    {
        for (unsigned int i = 0; i < 3; ++i){
            double incr_i = target[i] - base_coordinates[i];

            if (fabs(incr_i) > 0.5 * periods[i]){
                base_coordinates[i] += GetSign(incr_i) * periods[i];
            }
        }

    }

    static inline void TransformToClosestPeriodicCoordinates(const array_1d<double,3>& target, array_1d<double,3>& base_coordinates, const double periods[3])
    {
        for (unsigned int i = 0; i < 3; ++i){
            double incr_i = target[i] - base_coordinates[i];

            if (fabs(incr_i) > 0.5 * periods[i]){
                base_coordinates[i] += GetSign(incr_i) * periods[i];
            }
        }

    }

    static inline void GetBoxCenter(double box_center[3], const double min_point[3], const double max_point[3])
    {
        for (unsigned int i = 0; i < 3; ++i){
            box_center[i] = 0.5 * (min_point[i] + max_point[i]);

            if (min_point[i] > max_point[i]){ // the box is broken by the boundary in this dimension:  ]    x  * [   (outside),  ]    x    [ * (inside)
                const double& min = mDomainMin[i];
                const double& max = mDomainMax[i];
                box_center[i] += 0.5 * (max - min); // The center of the box and of its complementary are always half a domain period apart: |..]..x..[...*.| (our bet: we assume we should ADD to go from x to *)
                if (box_center[i] > max){ // we made a mistake, we should have gone the opposite way
                    box_center[i] -= max - min;
                }
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

        SphericParticle* p_particle = static_cast<SphericParticle*>(&*rObject);
        double radius = p_particle->GetSearchRadius();

        for(std::size_t i = 0; i < 3; ++i)
        {
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
        }
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
        (void) Radius;
        CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
    }

    static inline void CalculateCenter(const PointerType& rObject, PointType& rCenter)
    {
        rCenter = rObject->GetGeometry()[0];
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        double rObj_2_to_rObj_1[3];
        PointType& point_1 = rObj_1->GetGeometry()[0];
        const double coors_1[3] = {point_1[0], point_1[1], point_1[2]};
        PointType& point_2 = rObj_2->GetGeometry()[0];
        const double coors_2[3] = {point_2[0], point_2[1], point_2[2]};

        PeriodicSubstract(coors_1, coors_2, rObj_2_to_rObj_1);

        const double distance_2 = DEM_INNER_PRODUCT_3(rObj_2_to_rObj_1, rObj_2_to_rObj_1);

        SphericParticle* p_particle1 = static_cast<SphericParticle*>(&*rObj_1);
        SphericParticle* p_particle2 = static_cast<SphericParticle*>(&*rObj_2);
        const double radius_sum      = p_particle1->GetSearchRadius() + p_particle2->GetSearchRadius();
        const bool intersect         = floatle((distance_2 - radius_sum * radius_sum),0);
        return intersect;
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, const double& radius_1)
    {
        return Intersection(rObj_1, rObj_2);
    }

    //******************************************************************************************************************

    static inline bool IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        const array_1d<double, 3>& center_of_particle = rObject->GetGeometry()[0];

        SphericParticle* p_particle = static_cast<SphericParticle*>(&*rObject);
        const double& radius = p_particle->GetSearchRadius();

        bool intersect;

        if (mDomainIsPeriodic){
            double expanded_box_min[3];
            double expanded_box_max[3];

            for (unsigned int i = 0; i < 3; ++i){
                expanded_box_min[i] = rLowPoint[i] - radius;
                expanded_box_max[i] = rHighPoint[i] + radius;
            }

            double box_center[3];
            GetBoxCenter(box_center, expanded_box_min, expanded_box_max);
            //double box_center[3] = {0.5 * (rLowPoint[0] + rHighPoint[0]), 0.5 * (rLowPoint[1] + rHighPoint[1]), 0.5 * (rLowPoint[2] + rHighPoint[2])};
            double representative_center_of_particle[3] = {center_of_particle[0],
                                                           center_of_particle[1],
                                                           center_of_particle[2]};
            TransformToClosestPeriodicCoordinates(box_center, representative_center_of_particle);

            for (unsigned int i = 0; i < 3; ++i){
                const bool is_broken = rLowPoint[i] > rHighPoint[i];

                if (is_broken){ // i.e., we have |  ]  [ | in this direction
                    intersect = floatge(expanded_box_min[i], representative_center_of_particle[i])
                             && floatle(expanded_box_max[i], representative_center_of_particle[i]);
                }

                else { // i.e., we have |  [ ] | in this direction
                    intersect = floatle(expanded_box_min[i], representative_center_of_particle[i])
                             && floatge(expanded_box_max[i], representative_center_of_particle[i]);
                }
            }
        }

        else {
            for (unsigned int i = 0; i < 3; ++i){
                intersect = floatle(rLowPoint[i]  - radius, center_of_particle[i])
                         && floatge(rHighPoint[i] + radius, center_of_particle[i]);
            }
        }

        return  intersect;
    }

    static inline bool IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
        return  IntersectionBox(rObject, rLowPoint, rHighPoint);
    }

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        PointType& point_1 = rObj_1->GetGeometry()[0];
        const double coors_1[3] = {point_1[0], point_1[1], point_1[2]};
        PointType& point_2 = rObj_2->GetGeometry()[0];
        const double coors_2[3] = {point_2[0], point_2[1], point_2[2]};
        double rObj_2_to_rObj_1[3];
        PeriodicSubstract(coors_1, coors_2, rObj_2_to_rObj_1);
        distance = DEM_MODULUS_3(rObj_2_to_rObj_1);
    }

    static inline double GetObjectRadius(const PointerType& rObject, const double& Radius)
    {
        return GetObjectRadius(rObject);
    }

    static inline double GetObjectRadius(const PointerType& rObject)
    {
        SphericParticle* p_particle = static_cast<SphericParticle*>(&*rObject);
        return p_particle->GetSearchRadius();
    }

    static double mDomainPeriods[3];
    static double mDomainMin[3];
    static double mDomainMax[3];
    static bool mDomainIsPeriodic;

    /// Turn back information as a string.
    virtual std::string Info() const {return " Spatial Containers Configure for Discrete Particles"; }

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

    static inline void PeriodicSubstract(const double a[3], const double b[3], double c[3])
    {
        for (unsigned int i = 0; i < 3; ++i){
                c[i] = a[i] - b[i];
        }

        if (mDomainIsPeriodic){ // Periods have been set (the domain is periodic)
            for (unsigned int i = 0; i < 3; ++i){
                if (std::fabs(c[i]) > 0.5 * mDomainPeriods[i]){ // the objects are closer through the boundary
                    c[i] -= GetSign(c[i]) * mDomainPeriods[i];
                }
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
double DiscreteParticleConfigure<TDimension>::mDomainMin[] = {0.0, 0.0, 0.0};
template <std::size_t TDimension>
double DiscreteParticleConfigure<TDimension>::mDomainMax[] = {-1.0, -1.0, -1.0};
template <std::size_t TDimension>
bool DiscreteParticleConfigure<TDimension>::mDomainIsPeriodic = false;

}   // namespace Kratos.
#endif	/* DISCRETE_PARTICLE_CONFIGURE_H */
