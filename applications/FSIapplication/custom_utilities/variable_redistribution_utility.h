//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Jordi Cotela
//

#if !defined(KRATOS_VARIABLE_REDISTRIBUTION_UTILITY_H_INCLUDED )
#define  KRATOS_VARIABLE_REDISTRIBUTION_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"


namespace Kratos
{
  ///@addtogroup FSIApplication
  ///@{

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

  /// Helper utility to transform between point-wise nodal variables and distributed values. 
  /** The functions are desinged so that both sets of values have the same L2 norm over the
   *  conditions of the provided ModelPart (up to tolerance).
   *  A typical use case is to transform a set of point forces to area-distributed loads
   *  or vice-versa.
   */
  class VariableRedistributionUtility
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of VariableRedistributionUtility
      KRATOS_CLASS_POINTER_DEFINITION(VariableRedistributionUtility);

      ///@}
      ///@name Life Cycle
      ///@{


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      static void ConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< double >& rDistributedVariable,
          const Variable< double >& rPointVariable);

      
      static void ConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< array_1d<double,3> >& rDistributedVariable,
          const Variable< array_1d<double,3> >& rPointVariable);

      static void DistributePointValues(
          ModelPart& rModelPart,
          const Variable< double >& rPointVariable,
          const Variable< double >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

      static void DistributePointValues(
          ModelPart& rModelPart,
          const Variable< array_1d<double,3> >& rPointVariable,
          const Variable< array_1d<double,3> >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

      ///@}
      ///@name Access
      ///@{


      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{


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

      template< class TValueType >
      static void CallSpecializedConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable);

      template< class TValueType >
      static void CallSpecializedDistributePointValues(
          ModelPart& rModelPart,
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
      static void SpecializedConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable);

      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
      static void SpecializedDistributePointValues(
          ModelPart& rModelPart,
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

      static void ComputeNodalSizes(ModelPart& rModelPart);

      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TNumNodes >
      static void ConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>& rMassMatrix);

      template< unsigned int TNumNodes, class TValueType >
      static void UpdateDistributionRHS(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable,
          boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>& rMassMatrix);

      template< class TValueType >
      static double SolveDistributionIteration(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable);

      template< class TValueType >
      static const Variable< TValueType >& GetRHSVariable(const Variable<TValueType>& rVariable);

      template< class TValueType >
      static double AddToNorm(TValueType NodalValue, double NodalSize);
            
      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Default constructor.
      VariableRedistributionUtility();

      /// Assignment operator.
      VariableRedistributionUtility& operator=(VariableRedistributionUtility const& rOther);

      /// Copy constructor.
      VariableRedistributionUtility(VariableRedistributionUtility const& rOther);


      ///@}

    }; // Class VariableRedistributionUtility

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VARIABLE_REDISTRIBUTION_UTILITY_H_INCLUDED  defined
