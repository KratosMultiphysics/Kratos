/*
 * File:   sparticle_block_configure.h
 * Author: gcasas
 *
 * Created on 12 de septiembre de 2011, 15:58
 */

#ifndef PARTICLE_BLOCK_CONFIGURE_H
#define	PARTICLE_BLOCK_CONFIGURE_H

// System includes

// Project includes
#include "custom_utilities/GeometryFunction.h"
#include "dem_fem__application.h"
//Database includes


namespace Kratos
{
    using namespace GeometryFunction;

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
class ContactPairBlockParticle{

public:

    T value[2];
    ContactPairBlockParticle(){}

    ContactPairBlockParticle(const T& First,const T& Second){
        value[0]  = First;
        value[1]  = Second;
        }

    ~ContactPairBlockParticle(){}

    T& operator[](std::size_t index){
        return value[index];
        }

    T const& operator[](std::size_t index) const{
        return value[index];
        }

    ContactPairBlockParticle& operator = (ContactPairBlockParticle& Pair){
        value[0] = Pair[0];
        value[1] = Pair[1];
        return *this;
        }

    inline bool operator == (const ContactPairBlockParticle& Pair){
        return  (value[0] == Pair[0]) && (value[1] == Pair[1]) ;
        }

    inline bool  operator == (const ContactPairBlockParticle& Pair) const{
        return ( (value[0] == Pair[0]) && (value[1] == Pair[1]) ) ;
        }

    inline std::ostream& operator << (const ContactPairBlockParticle& Pair){
        std::ostream OStream;
        OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl;
        return OStream;
        }

    inline std::ostream& operator << (const ContactPairBlockParticle& Pair) const{
        std::ostream OStream;
        OStream << "Object 1 =  "<<  Pair[0] << "  " << "Object 2 =" << Pair[1] << std::endl;
        return OStream;
        }
};


template <class TParticle>
class ParticleBlockConfigure{

public:

      ///@name Type Definitions
      ///@{
    enum {Dimension = 3};
    typedef TParticle                                                   ParticleType;
    typedef Point< 3, double>                                           PointType;
    //typedef typename ParticleType::DistanceIteratorType                 DistanceIteratorType;
    typedef std::vector<double>::iterator                               DistanceIteratorType;
    typedef typename ParticleType::Pointer                              PointerType;
    typedef typename std::vector<typename ParticleType::Pointer>        ContainerType;
    typedef typename std::vector<PointerType>::iterator                 IteratorType;
    typedef ContainerType                                               ResultContainerType;
    typedef IteratorType                                                ResultIteratorType;


      /// Contact Pairs
    typedef ContactPairBlockParticle<PointerType>                       ContactPairType;
    typedef  std::vector<ContactPairType>                               ContainerContactType;
    typedef  typename ContainerContactType::iterator                    IteratorContactType;
    typedef  typename ContainerContactType::value_type                  PointerContactType;

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(ParticleBlockConfigure);

      ///@}
      ///@name Life Cycle
      ///@{

    ParticleBlockConfigure(){};
    virtual ~ParticleBlockConfigure(){}

          ///@}
          ///@name Operators
          ///@{


          ///@}
          ///@name Operations
          ///@{

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        ///Cfeng:rObject could be particle or block

        double Length = 0.0;

        array_1d<double, 3> Coord;
        array_1d<double, 3> Centroid;

        noalias(Centroid) = ZeroVector(3);

        for(std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
        {
            Coord = rObject->GetGeometry()(inode)->Coordinates();

            Centroid  += Coord / (double)rObject->GetGeometry().size();
        }


        if(rObject->GetGeometry().size() == 1)
        {
            Length = rObject->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        }
        else
        {
             double xyz_min[3] = { 1e20,  1e20,  1e20};
             double xyz_max[3] = {-1e20, -1e20, -1e20};
                
             for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
             {
                Coord = rObject->GetGeometry()(inode)->Coordinates();

                xyz_min[0] = (xyz_min[0] > Coord[0]) ? Coord[0] : xyz_min[0];
                xyz_min[1] = (xyz_min[1] > Coord[1]) ? Coord[1] : xyz_min[1];
                xyz_min[2] = (xyz_min[2] > Coord[2]) ? Coord[2] : xyz_min[2];

                xyz_max[0] = (xyz_max[0] < Coord[0]) ? Coord[0] : xyz_min[0];
                xyz_max[1] = (xyz_max[1] < Coord[1]) ? Coord[1] : xyz_min[1];
                xyz_max[2] = (xyz_max[2] < Coord[2]) ? Coord[2] : xyz_min[2];
             }
             
             double dist1 = xyz_max[0] - xyz_min[0];
             double dist2 = xyz_max[1] - xyz_min[1];
             double dist3 = xyz_max[2] - xyz_min[2];
             
             Length = dist1;            
             Length = (Length < dist2) ? dist2 : Length;
             Length = (Length < dist3) ? dist3 : Length;
             
             Length = Length * 0.5;
        }

   
        rLowPoint  = Centroid;
        rHighPoint = Centroid;
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -Length;
            rHighPoint[i] +=  Length;
        }

    }

    //******************************************************************************************************************

     static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
     {
         //Cfeng: rObj_1 is particle,  and rObj_2 is block

         bool If_PB_Contact = false;

        double Particle_Coord[3] = {0.0};
        Particle_Coord[0] = rObj_1->GetGeometry()(0)->Coordinates()[0];
        Particle_Coord[1] = rObj_1->GetGeometry()(0)->Coordinates()[1];
        Particle_Coord[2] = rObj_1->GetGeometry()(0)->Coordinates()[2];

        double rad = rObj_1->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);

        double Centroid[3] = {0.0};
        for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++)
        {
            Centroid[0] += rObj_2->GetGeometry()(inode)->Coordinates()[0] / (double)rObj_2->GetGeometry().size();
            Centroid[1] += rObj_2->GetGeometry()(inode)->Coordinates()[1] / (double)rObj_2->GetGeometry().size();
            Centroid[2] += rObj_2->GetGeometry()(inode)->Coordinates()[2] / (double)rObj_2->GetGeometry().size();
        }

         std::size_t dim = rObj_2->GetGeometry().WorkingSpaceDimension();
         if(dim == 2)
         {
             double Coord1[3]= {0.0};
             double Coord2[3]= {0.0};

            for (std::size_t iedge = 0; iedge < rObj_2->GetGeometry().Edges().size(); iedge++)
            {
                if( rObj_2->GetValue(IF_BOUNDARY_FACE)[iedge] == 1 )
                {
                    Coord1[0] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[0];
                    Coord1[1] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[1];
                    Coord1[2] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[2];

                    Coord2[0] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[0];
                    Coord2[1] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[1];
                    Coord2[2] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[2];

                    //GeometryFunction::Compute2DimElementEdgeLocalSystem(Coord1, Coord2, Centroid, LocalCoordSystem);
                    If_PB_Contact = GeometryFunction::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Centroid, Particle_Coord, rad);

                    if(If_PB_Contact == true)
                    {
                        break;
                    }
                }
    
            } 
         }        
         else if(dim == 3)
         {
            double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };
            for (std::size_t iface = 0; iface < rObj_2->GetGeometry().Faces().size(); iface++)
            {
                if( rObj_2->GetValue(IF_BOUNDARY_FACE)[iface] == 1 )
                {
                    ////Cfeng: Triangle
                    int FaceNodeTotal = 3;

                    Coord[0][0] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[0];
                    Coord[0][1] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[1];
                    Coord[0][2] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[2];

                    Coord[1][0] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[0];
                    Coord[1][1] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[1];
                    Coord[1][2] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[2];

                    Coord[2][0] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[0];
                    Coord[2][1] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[1];
                    Coord[2][2] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[2];

                    if(rObj_2->GetGeometry().Faces()[iface].size() == 4)
                    {
                        Coord[3][0] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[0];
                        Coord[3][1] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[1];
                        Coord[3][2] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[2];

                      ////Cfeng: Quadral
                      FaceNodeTotal = 4;
                    }

                    //GeometryFunction::Compute2DimElementEdgeLocalSystem(Coord1, Coord2, Centroid, LocalCoordSystem);
                    If_PB_Contact = GeometryFunction::JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord, Centroid, Particle_Coord, rad);

                    if(If_PB_Contact == true)
                    {
                        break;
                    }
                }

            }
             
         }

        return If_PB_Contact;
      }

    //******************************************************************************************************************

     static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
     {
         //Cfeng: rObject is block

        array_1d<double, 3> Centroid;
        array_1d<double, 3> Coord;

        noalias(Centroid) = ZeroVector(3);

        for(std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
        {
            Centroid  += rObject->GetGeometry()(inode)->Coordinates() / (double)rObject->GetGeometry().size();
        }

        double Length = 0.0;

        if(rObject->GetGeometry().size() == 1)
        {
            Length = rObject->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        }
        else
        {
             double xyz_min[3] = { 1e20,  1e20,  1e20};
             double xyz_max[3] = {-1e20, -1e20, -1e20};

             for (std::size_t inode = 0; inode < rObject->GetGeometry().size(); inode++)
             {
                Coord = rObject->GetGeometry()(inode)->Coordinates();

                xyz_min[0] = (xyz_min[0] > Coord[0]) ? Coord[0] : xyz_min[0];
                xyz_min[1] = (xyz_min[1] > Coord[1]) ? Coord[1] : xyz_min[1];
                xyz_min[2] = (xyz_min[2] > Coord[2]) ? Coord[2] : xyz_min[2];

                xyz_max[0] = (xyz_max[0] < Coord[0]) ? Coord[0] : xyz_min[0];
                xyz_max[1] = (xyz_max[1] < Coord[1]) ? Coord[1] : xyz_min[1];
                xyz_max[2] = (xyz_max[2] < Coord[2]) ? Coord[2] : xyz_min[2];
             }

             double dist1 = xyz_max[0] - xyz_min[0];
             double dist2 = xyz_max[1] - xyz_min[1];
             double dist3 = xyz_max[2] - xyz_min[2];

             Length = dist1;
             Length = (Length < dist2) ? dist2 : Length;
             Length = (Length < dist3) ? dist3 : Length;

             Length = Length * 0.5;
        }

        bool intersect = (rLowPoint[0] - Length <= Centroid[0] && rLowPoint[1] - Length <= Centroid[1] && rLowPoint[2] - Length <= Centroid[2] &&
            rHighPoint[0] + Length >= Centroid[0] && rHighPoint[1] + Length >= Centroid[1] && rHighPoint[2] + Length >= Centroid[2]);

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
    ParticleBlockConfigure& operator=(ParticleBlockConfigure const& rOther);

    /// Copy constructor.
    ParticleBlockConfigure(ParticleBlockConfigure const& rOther);

    ///@}

    }; // Class ParticleBlockConfigure

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template <class TParticle>
    inline std::istream& operator >> (std::istream& rIStream, ParticleBlockConfigure<TParticle> & rThis){
        return rIStream;
        }

    /// output stream function
    template <class TParticle>
    inline std::ostream& operator << (std::ostream& rOStream, const ParticleBlockConfigure<TParticle>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }
    ///@}

}   // namespace Kratos.
#endif	/* PARTICLE_CONFIGURE_H */