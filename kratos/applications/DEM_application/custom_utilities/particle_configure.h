/*
 * File:   spheric_particle_configure.h
 * Author: gcasas
 *
 * Created on 12 de septiembre de 2011, 15:58
 */

#ifndef PARTICLE_CONFIGURE_H
#define	PARTICLE_CONFIGURE_H

// System includes

// Project includes

//#include "custom_utilities/circular_particle.h"
//#include "custom_utilities/circular_particle_hertzian.h"
//#include "custom_utilities/spheric_particle.h"
#include "custom_utilities/spheric_particle_hertzian.h"
//#include "custom_utilities/spheric_rotating_particle.h"
//#include "custom_utilities/spheric_rotating_particle_hertzian.h"

//Database includes


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

  /// Short class definition.
  /** Detail class definition.
  */

template <class T>
class ContactPair{

public:

    T value[2];
    ContactPair(){}

    ContactPair(const T& First,const T& Second){
        value[0]  = First;
        value[1]  = Second;
        }

    ~ContactPair(){}

    T& operator[](std::size_t index){
        return value[index];
        }

    T const& operator[](std::size_t index) const{
        return value[index];
        }

    ContactPair& operator = (ContactPair& Pair){
        value[0] = Pair[0];
        value[1] = Pair[1];
        return *this;
        }

    inline bool operator == (const ContactPair& Pair){
        return  (value[0] == Pair[0]) && (value[1] == Pair[1]) ;
        }

    inline bool  operator == (const ContactPair& Pair) const{
        return ( (value[0] == Pair[0]) && (value[1] == Pair[1]) ) ;
        }

    inline std::ostream& operator << (const ContactPair& Pair){
        std::ostream OStream;
        OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl;
        return OStream;
        }

    inline std::ostream& operator << (const ContactPair& Pair) const{
        std::ostream OStream;
        OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl;
        return OStream;
        }
};

template <class TParticle>
class ParticleConfigure{

public:

      ///@name Type Definitions
      ///@{
    enum {Dimension = 3};
    typedef TParticle                                                   ParticleType;
    typedef Point< 3, double>                                           PointType;
    typedef typename ParticleType::DistanceIteratorType                 DistanceIteratorType;
    typedef typename ParticleType::Pointer                              PointerType;
    typedef typename std::vector<typename ParticleType::Pointer>        ContainerType;
    typedef typename std::vector<PointerType>::iterator                 IteratorType;
    typedef ContainerType                                               ResultContainerType;
    typedef IteratorType                                                ResultIteratorType;

      /// Contact Pairs
    typedef ContactPair<PointerType>                                    ContactPairType;
    typedef  std::vector<ContactPairType>                               ContainerContactType;
    typedef  typename ContainerContactType::iterator                    IteratorContactType;
    typedef  typename ContainerContactType::value_type                  PointerContactType;

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(ParticleConfigure);

      ///@}
      ///@name Life Cycle
      ///@{

    ParticleConfigure(){};
    virtual ~ParticleConfigure(){}

          ///@}
          ///@name Operators
          ///@{


          ///@}
          ///@name Operations
          ///@{

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint){
        rLowPoint = *(rObject->GetPointerToCenterNode());
        rHighPoint = *(rObject->GetPointerToCenterNode());
        double radius = rObject->GetRadius();
        for(std::size_t i = 0; i < 3; i++){
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
            }
        }

    //******************************************************************************************************************

     static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2){
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetPosition() - rObj_2->GetPosition();
        double distance_2 = rObj_2_to_rObj_1[0] * rObj_2_to_rObj_1[0] + rObj_2_to_rObj_1[1] * rObj_2_to_rObj_1[1] + rObj_2_to_rObj_1[2] * rObj_2_to_rObj_1[2];
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_1 = rObj_1->GetRadius();
        double radius_2 = rObj_2->GetRadius();
        double radius_sum = radius_1 + radius_2;
        bool intersect = (distance_2 - radius_sum * radius_sum) <= 0;
        return intersect;
        }

    //******************************************************************************************************************

     static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint){
//        double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle = rObject->GetPosition();
        double radius = rObject->GetRadius();
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
            rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
        return  intersect;
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
    ParticleConfigure& operator=(ParticleConfigure const& rOther);

    /// Copy constructor.
    ParticleConfigure(ParticleConfigure const& rOther);

    ///@}

    }; // Class ParticleConfigure

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template <class TParticle>
    inline std::istream& operator >> (std::istream& rIStream, ParticleConfigure<TParticle> & rThis){
        return rIStream;
        }

    /// output stream function
    template <class TParticle>
    inline std::ostream& operator << (std::ostream& rOStream, const ParticleConfigure<TParticle>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }
    ///@}

}   // namespace Kratos.
#endif	/* PARTICLE_CONFIGURE_H */