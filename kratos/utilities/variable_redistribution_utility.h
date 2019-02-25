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
  class KRATOS_API(KRATOS_CORE) VariableRedistributionUtility
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

    /**
     * @brief Specialization of ConvertDistributeValuesToPoint according to the variable value type
     *
     * @tparam TValueType variables value type (double or array_1<double,3>)
     * @param rModelPart model part in where the point values accumulation is done
     * @param rDistributedVariable origin distributed variable
     * @param rPointVariable destination point variable
     */
      template< class TValueType >
      static void CallSpecializedConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable);

    /**
     * @brief Specialization of DistributePointValues according to the variable value type
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rPointVariable origin point variable
     * @param rDistributedVariable destination distributed variable
     * @param Tolerance maximum allowed difference (in L2 norm) between origin and destination values.
     * @param MaximumIterations maximum number of iterations for the procedure.
     */
      template< class TValueType >
      static void CallSpecializedDistributePointValues(
          ModelPart& rModelPart,
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

    /**
     * @brief ConvertDistributedValuesToPoint specialization according to geometry family, points number and value type
     *
     * @tparam TFamily geometry family in the model part in where the accumulation is done
     * @tparam TPointNumber points number in the geometries in where the accumulation is done
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the point values accumulation is done
     * @param rDistributedVariable origin distributed variable
     * @param rPointVariable destination point variable
     */
      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
      static void SpecializedConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable);

    /**
     * @brief DummyConvertDistributedValuesToPoint
     * This class does nothing, it is only used in case there is no conditions in the current
     * partition to perform the communication operations that are done in the "standard" case.
     * Otherwise, the MPI synchronism is broken
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the point values accumulation is done
     * @param rDistributedVariable origin distributed variable
     * @param rPointVariable destination point variable
     */
      template< class TValueType >
      static void DummySpecializedConvertDistributedValuesToPoint(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          const Variable< TValueType >& rPointVariable);

    /**
     * @brief DistributePointValues specialization according to geometry family, points number and value type
     *
     * @tparam TFamily geometry family in the model part in where the distribution is done
     * @tparam TPointNumber points number in the geometries in where the distribution is done
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rPointVariable origin point variable
     * @param rDistributedVariable destination distributed variable
     * @param Tolerance maximum allowed difference (in L2 norm) between origin and destination values.
     * @param MaximumIterations maximum number of iterations for the procedure.
     */
      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TPointNumber, class TValueType >
      static void SpecializedDistributePointValues(
          ModelPart& rModelPart,
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

    /**
     * @brief Dummy SpecializedDistributePointValues.
     * This class does nothing, it is only used in case there is no conditions in the current
     * partition to perform the communication operations that are done in the "standard" case.
     * Otherwise, the MPI synchronism is broken
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rDistributedVariable destination distributed variable
     * @param Tolerance maximum allowed difference (in L2 norm) between origin and destination values.
     * @param MaximumIterations maximum number of iterations for the procedure.
     */
      template< class TValueType >
      static void DummySpecializedDistributePointValues(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable,
          double Tolerance,
          unsigned int MaximumIterations);

    /**
     * @brief This function computes the NODAL_MAUX values
     *
     * @param rModelPart model part in where the NODAL_MAUX is computed
     */
      static void ComputeNodalSizes(ModelPart& rModelPart);

    /**
     * @brief Fills the given matrix with the consistent mass matrix values
     *
     * @tparam TFamily geometry family in the model part in where the operation is done
     * @tparam TNumNodes nodes number of the geometries in where the operation is done
     * @param rMassMatrix computed consistent mass matrix that is to be filled
     */
      template< GeometryData::KratosGeometryFamily TFamily, unsigned int TNumNodes >
      static void ConsistentMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix);

    /**
     * @brief Computes the RHS of the distribution problem
     *
     * @tparam TNumNodes nodes number of the geometries in where the operation is done
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rPointVariable origin point variable
     * @param rDistributedVariable destination distributed variable
     * @param rMassMatrix mass matrix that is filled
     */
      template< unsigned int TNumNodes, class TValueType >
      static void UpdateDistributionRHS(
          ModelPart& rModelPart,
          const Variable< TValueType >& rPointVariable,
          const Variable< TValueType >& rDistributedVariable,
          BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix);

    /**
     * @brief Dummy computation of the RHS of the distribution problem
     * This class does nothing, it is only used in case there is no conditions in the current
     * partition to perform the communication operations that are done in the "standard" case
     * Otherwise, the MPI synchronism is broken
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rDistributedVariable destination distributed variable
     */
      template< class TValueType >
      static void DummyUpdateDistributionRHS(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable);

    /**
     * @brief Function that solves the distribution problem. It is called at each iteration.
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rDistributedVariable destination distributed variable
     * @return double accumulated error norm of the distribution
     */
      template< class TValueType >
      static double SolveDistributionIteration(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable);

    /**
     * @brief Dummy function that solves the distribution problem. It is called at each iteration.
     * This class does nothing, it is only used in case there is no conditions in the current
     * partition to perform the communication operations that are done in the "standard" case
     * Otherwise, the MPI synchronism is broken
     * 
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rModelPart model part in where the distribution is done
     * @param rDistributedVariable destination distributed variable
     * @return double accumulated error norm of the distribution
     */
      template< class TValueType >
      static double DummySolveDistributionIteration(
          ModelPart& rModelPart,
          const Variable< TValueType >& rDistributedVariable);

    /**
     * @brief Auxiliar method to retrieve the variable used in the RHS of the distribution problem
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rVariable reference to the RHS variable of the distribution problem
     * @return const Variable< TValueType >& reference to the RHS variable of the distribution problem
     */
      template< class TValueType >
      static const Variable< TValueType >& GetRHSVariable(const Variable<TValueType>& rVariable);

    /**
     * @brief Auxiliar function to compute the error norm according to the variable value type
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param NodalValue nodal value
     * @param NodalSize area associated to the current node
     * @return double product of the nodal value norm times the nodal size
     */
      template< class TValueType >
      static double AddToNorm(TValueType NodalValue, double NodalSize);

    /**
     * @brief Auxiliar function to perform a threadsafe addition
     *
     * @tparam TValueType variables value type (double or array_1d<double,3>)
     * @param rLHS left hand side of the summation (accumulated value)
     * @param rRHS right hand side of the summation (value to be accumulated in the LHS)
     */
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
