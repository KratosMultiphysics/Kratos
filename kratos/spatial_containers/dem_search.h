//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    clabra
//

#pragma once

// include kratos definitions
#include "includes/define.h"

// System includes
#include <string>
#include <iostream>

// External includes
#include "spatial_containers/spatial_search.h"

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

      void Initialize(Element::Pointer baseElem)
      {
          for(std::size_t i = 0; i < Dimension; i++)
              coord[i] = baseElem->GetGeometry()[0][i];

          pNaseElem = baseElem;

//           mRadius = baseElem->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
      }

      void Initialize(Element::Pointer baseElem, double Radius)
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

      Element::Pointer pNaseElem;

      void operator=(Point const& Other){
         for(std::size_t i = 0; i < Dimension; i++)
            coord[i] = Other[i];
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
      DEMSearch(const double domain_min_x = 0.0, const double domain_min_y = 0.0, const double domain_min_z = 0.0,
                const double domain_max_x = -1.0, const double domain_max_y = -1.0, const double domain_max_z = -1.0)
      {
        mDomainMin[0] = domain_min_x;
        mDomainMin[1] = domain_min_y;
        mDomainMin[2] = domain_min_z;
        mDomainMax[0] = domain_max_x;
        mDomainMax[1] = domain_max_y;
        mDomainMax[2] = domain_max_z;
        TDerived::ElementConfigureType::SetDomain(domain_min_x, domain_min_y, domain_min_z, domain_max_x, domain_max_y, domain_max_z);
        TDerived::NodeConfigureType::SetDomain(domain_min_x, domain_min_y, domain_min_z, domain_max_x, domain_max_y, domain_max_z);
        mDomainPeriodicity = TDerived::ElementConfigureType::GetDomainPeriodicity();
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
          VectorDistanceType& rResultsDistance ) override
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
          VectorResultElementsContainerType& rResults ) override
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
          VectorDistanceType& rResultsDistance ) override
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }

      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) override
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusInclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }

      void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults ) override
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults);
      }

      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults ) override
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
          VectorDistanceType& rResultsDistance ) override
      {
          static_cast<TDerived*>(this)->SearchGeometricalInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }

      void SearchElementsOverConditionsInRadiusInclusive (
          ConditionsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance ) override
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
      virtual std::string Info() const override
      {
          std::stringstream buffer;
          buffer << "DemSearch" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "DemSearch";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override {}


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
      bool mDomainPeriodicity;
      array_1d<double, 3> mDomainMin;
      array_1d<double, 3> mDomainMax;

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
