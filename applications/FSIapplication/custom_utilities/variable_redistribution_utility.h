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

      /// Tranform a variable distributed over the conditions of rModelPart to a set of concentrated nodal values.
      /** The origin and destination values have the same L2 norm over the set of conditions.
       *  A typical use case is to transform a distributed load into an equivalent set of point loads.
       *  Version for scalar magnitudes.
       *  @param rModelPart The model part of the problem
       *  @param rDistributedVariable The variable containing the distributed (origin) values.
       *  @param rPointVariable The variable that will contain the point (destination) values.
       */
      static void ConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< double >& rDistributedVariable,
          const Variable< double >& rPointVariable);

      
      /// Tranform a variable distributed over the conditions of rModelPart to a set of concentrated nodal values.
      /** The origin and destination values have the same L2 norm over the set of conditions.
       *  A typical use case is to transform a distributed load into an equivalent set of point loads.
       *  Version for vector magnitudes.
       *  @param rModelPart The model part of the problem
       *  @param rDistributedVariable The variable containing the distributed (origin) values.
       *  @param rPointVariable The variable that will contain the point (destination) values.
       */
      static void ConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< array_1d<double,3> >& rDistributedVariable,
          const Variable< array_1d<double,3> >& rPointVariable);

      /// Tranform a set of concentrated nodal values to a variable distributed over the conditions of rModelPart.
      /** The origin and destination values have the same L2 norm over the set of conditions.
       *  A typical use case is to transform a set of point loads into an equivalent distributed load.
       *  Version for scalar magnitudes.
       *  @param rModelPart The model part of the problem
       *  @param rPointVariable The variable that will contain the point (origin) values.
       *  @param rDistributedVariable The variable containing the distributed (destination) values.
       *  @param Tolerance Maximum allowed difference (in L2 norm) between origin and destination values.
       *  @param MaximumIterations Maximum number of iterations for the procedure.
       */
      static void DistributePointValues(
          ModelPart& rModelPart,
          const Variable< double >& rPointVariable,
          const Variable< double >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

      /// Tranform a set of concentrated nodal values to a variable distributed over the conditions of rModelPart.
      /** The origin and destination values have the same L2 norm over the set of conditions.
       *  A typical use case is to transform a set of point loads into an equivalent distributed load.
       *  Version for vector magnitudes.
       *  @param rModelPart The model part of the problem
       *  @param rPointVariable The variable that will contain the point (origin) values.
       *  @param rDistributedVariable The variable containing the distributed (destination) values.
       *  @param Tolerance Maximum allowed difference (in L2 norm) between origin and destination values.
       *  @param MaximumIterations Maximum number of iterations for the procedure.
       */
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
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>& rMassMatrix);

      template< class TValueType >
      static double SolveDistributionIteration(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable);

      template< class TValueType >
      static const Variable< TValueType >& GetRHSVariable(const Variable<TValueType>& rVariable);

      template< class TValueType >
      static double AddToNorm(TValueType NodalValue, double NodalSize);

      template< class TValueType >
      static void ThreadsafeAdd(TValueType& rLHS, const TValueType& rRHS);
            
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
