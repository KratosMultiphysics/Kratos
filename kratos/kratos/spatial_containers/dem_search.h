//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_DEM_SEARCH_H_INCLUDED )
#define  KRATOS_DEM_SEARCH_H_INCLUDED

// include kratos definitions
#include "includes/define.h"

// System includes
#include <string>
#include <iostream>

// External includes
#include "spatial_containers/spatial_search.h"

/* Timer defines */
#include "utilities/timer.h"
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

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

template<std::size_t dim, class T>
class PointDistance2
{
    public:
        inline double operator()( T const& p1, T const& p2 )
        {
            double dist = 0.0;

            double tmp1 = p1[0] - p2[0];
            double tmp2 = p1[1] - p2[1];
            double tmp3 = p1[2] - p2[2];

            dist += tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

//             dist = sqrt(dist) - p2.radius - p1.radius;

            return dist;
        }
};

template<std::size_t Dimension>
class RadiusPoint
{
    public:
      RadiusPoint() {}
      virtual ~RadiusPoint(){}

      void Initialize(SpatialSearch::ElementPointerType baseElem)
      {
          for(std::size_t i = 0; i < Dimension; i++)
              coord[i] = baseElem->GetGeometry()[0][i];

          pNaseElem = baseElem;

//           mRadius = baseElem->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
      }

      void Initialize(SpatialSearch::ElementPointerType baseElem, double Radius)
      {
          for(std::size_t i = 0; i < Dimension; i++)
              coord[i] = baseElem->GetGeometry()[0][i];

          pNaseElem = baseElem;

//           mRadius = Radius;
      }

    public:

      double       mRadius;

      double       coord[Dimension];

      double       & operator[](std::size_t i)       {return coord[i];}
      double const & operator[](std::size_t i) const {return coord[i];}

      SpatialSearch::ElementPointerType pNaseElem;

      void operator=(Point<Dimension> const& Other){
         for(std::size_t i = 0; i < Dimension; i++)
            coord[i] = Other.coord[i];
      }
};

template< std::size_t Dimension >
std::ostream & operator<<( std::ostream& rOut, RadiusPoint<Dimension> & rPoint){
   for(std::size_t i = 0 ; i < Dimension ; i++)
      rOut << rPoint[i] << " ";
   return rOut;
}

template< std::size_t Dimension >
std::istream & operator>>( std::istream& rIn, RadiusPoint<Dimension> & rPoint){
   for(std::size_t i = 0 ; i < Dimension ; i++)
      rIn >> rPoint[i];

   return rIn;
}

template< class TDerived >
class DEMSearch : public SpatialSearch
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(DEMSearch);

      typedef RadiusPoint<Dimension>        PointType;
      typedef PointType*                    PtrPointType;
      typedef std::vector<PtrPointType>*    PointVector;

      PointVector       searchPoints;

      // Tell the compiler which overloaded functions are using so we can
      // avoid "overloaded virtual function" warnings

      using SpatialSearch::SearchElementsInRadiusExclusive;
      using SpatialSearch::SearchElementsInRadiusInclusive;
      using SpatialSearch::SearchNodesInRadiusExclusive;
      using SpatialSearch::SearchNodesInRadiusInclusive;
      using SpatialSearch::SearchConditionsOverElementsInRadiusExclusive;
      using SpatialSearch::SearchConditionsOverElementsInRadiusInclusive;
      using SpatialSearch::SearchElementsOverConditionsInRadiusExclusive;
      using SpatialSearch::SearchElementsOverConditionsInRadiusInclusive;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      DEMSearch(){
        searchPoints = new std::vector<PtrPointType>(0);
      }

      /// Destructor.
      virtual ~DEMSearch(){
        delete searchPoints;
      }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchElementsInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }

      void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchElementsInRadiusInclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }

      void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchElementsInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults);
      }

      void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchElementsInRadiusInclusiveImplementation(StructureElements,InputElements,Radius,rResults);
      }

      void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }

      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusInclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }

      void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults);
      }

      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusInclusiveImplementation(StructureNodes,InputNodes,Radius,rResults);
      }

      void SearchConditionsOverElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultConditionsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchGeometricalInRadiusExclusiveImplementation(StructureElements,InputConditions,Radius,rResults,rResultsDistance);
      }

      void SearchConditionsOverElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultConditionsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchGeometricalInRadiusInclusiveImplementation(StructureElements,InputConditions,Radius,rResults,rResultsDistance);
      }

      void SearchElementsOverConditionsInRadiusExclusive (
          ConditionsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchGeometricalInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }

      void SearchElementsOverConditionsInRadiusInclusive (
          ConditionsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchGeometricalInRadiusInclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
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
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "DemSearch" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DemSearch";}

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
      DEMSearch& operator=(DEMSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      DEMSearch(DEMSearch const& rOther)
      {
          *this = rOther;
      }


      ///@}

    }; // Class DEMSearch

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                     DEMSearch& rThis){return rIStream;}
//
//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                     const DEMSearch& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);
//
//       return rOStream;
//     }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DEM_SEARCH_H_INCLUDED  defined
