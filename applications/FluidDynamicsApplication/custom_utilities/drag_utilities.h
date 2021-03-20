//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_DRAG_UTILITIES_H_INCLUDED )
#define  KRATOS_DRAG_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

  namespace detail
  {

    /** Helper class for looping through elements of a tuple
     *  @param TupleIndex current index in the tuple
     *  @note operations call their matching counterparts of the class with a decremented TupleIndex.
     *  A specialization must be provided for the case where TupleIndex=0, which is the reason for
     *  using classes in the first place.
     */
    template <std::size_t TupleIndex>
    struct RecurseTuple
    {
        /** Helper function for defining operator+= on tuples
         *  @note every type in the tuple must define an operator+=
         */
        template <typename ...TArguments>
        static void InPlaceAdd(std::tuple<TArguments...>& rLeft, const std::tuple<TArguments...>& rRight)
        {
            std::get<TupleIndex>(rLeft) += std::get<TupleIndex>(rRight);
            RecurseTuple<TupleIndex-1>::InPlaceAdd(rLeft, rRight);
        }

        /** Helper function for the atomic addition of tuples
         *  @param rLeft left hand side operand
         *  @param rRight right hand side operand
         *  @note every type in the tuple must have an overload for AtomicAdd
         */
        template <typename ...TArguments>
        static void AtomicAdd(std::tuple<TArguments...>& rLeft, const std::tuple<TArguments...>& rRight)
        {
            Kratos::AtomicAdd(std::get<TupleIndex>(rLeft), std::get<TupleIndex>(rRight));
            RecurseTuple<TupleIndex-1>::AtomicAdd(rLeft, rRight);
        }

        /** Call a set of operations on a component
         *  @param rComponent component (eg. Node<3> or Element) to perform the operation on
         *  @param rObjects output of the corresponding operations
         *  @param rOperations tuple of functors that take a TComponent and return a type that
         *  matches the corresponding element of TObjectTuple 
         */
        template <typename TComponent, typename TObjectTuple, typename TOperationTuple>
        static void OperateOnComponent(TComponent&& rComponent, TObjectTuple& rObjects, TOperationTuple&& rOperations)
        {
            std::get<TupleIndex>(rObjects) = std::get<TupleIndex>(rOperations)(std::forward<TComponent>(rComponent));
            RecurseTuple<TupleIndex-1>::OperateOnComponent(
                std::forward<TComponent>(rComponent),
                rObjects,
                std::forward<TOperationTuple>(rOperations)
            );
        }

        /** Sum all elements in a tuple across MPI processes
         *  @param rModelPart model part for providing the MPI communicator
         *  @param rObjects tuple of objects to be summed
         */
        template <typename TModelPart, typename TObjectTuple>
        static void MPISumAll(TModelPart&& rModelPart, TObjectTuple&& rObjects)
        {
            std::get<TupleIndex>(std::forward<TObjectTuple>(rObjects)) = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(std::get<TupleIndex>(rObjects));
        }
    };

    /// Specialization of tuple-looping class for TupleIndex=0
    /**
     * @note operations are identical to the parent template but do not call their
     * counterparts with TupleIndex-1.
     */
    template <>
    struct RecurseTuple<0>
    {
        template <typename ...TArguments>
        static void InPlaceAdd(std::tuple<TArguments...>& rLeft, const std::tuple<TArguments...>& rRight)
        {
            std::get<0>(rLeft) += std::get<0>(rRight);
        }

        template <typename ...TArguments>
        static void AtomicAdd(std::tuple<TArguments...>& rLeft, const std::tuple<TArguments...>& rRight)
        {
            Kratos::AtomicAdd(std::get<0>(rLeft), std::get<0>(rRight));
        }

        template <typename TComponent, typename TObjectTuple, typename TOperationTuple>
        static void OperateOnComponent(TComponent&& rComponent, TObjectTuple& rObjects, TOperationTuple&& rOperations)
        {
            std::get<0>(rObjects) = std::get<0>(rOperations)(
                std::forward<TComponent>(rComponent)
            );
        }

        template <typename TModelPart, typename TObjectTuple>
        static void MPISumAll(TModelPart&& rModelPart, TObjectTuple&& rObjects)
        {
            std::get<0>(std::forward<TObjectTuple>(rObjects)) = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(std::get<0>(rObjects));
        }
    };
  }

  /** operator+= for tuples to be used in SumReduction::LocalReduce
   *  @param rleft left hand operand
   *  @param rRight right hand operand
   *  @note must be directly in the Kratos namespace, otherwise it won't get picked up in SumReduction::LocalReduce
   */
  template <typename ...TArguments>
  std::tuple<TArguments...>& operator+=(std::tuple<TArguments...>& rLeft, const std::tuple<TArguments...>& rRight)
  {
      detail::RecurseTuple<std::tuple_size<std::tuple<TArguments...>>::value-1>::InPlaceAdd(rLeft, rRight);
      return rLeft;
  }

  /// Atomic addition of tuples to be used in SumReduction::ThreadSafeReduce
  template <typename ...TArguments>
  void AtomicAdd( std::tuple<TArguments...>& target, const std::tuple<TArguments...>& value )
  {
      detail::RecurseTuple<std::tuple_size<std::tuple<TArguments...>>::value-1>::AtomicAdd(target, value);
  }


  ///@addtogroup FluidDynamicsApplication
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

  /// Auxiliary utility to compute the drag force.
  /** For embedded formulations, this utility iterates all the elements of a provided model part. In this iteration
   * calls the calculate method of each element to compute the value of the variable DRAG_FORCE. If the element is split,
   * this method computes the integration of the stress term over the interface. Otherwise, the value is just zero.
   * The obtained values are accumulated to get the total drag force in the model part.
   *
   * Note that if there is more than one embedded object, one just needs to save the surrounding elements to each embedded
   * object in different submodelparts and call this process for each one of that submodelparts.
   *
   * For the body fitted slip case, it integrates the pressure stress term over the given submodelpart conditions (the
   * shear stress term is assumed to be zero).
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DragUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    using GeometryType = Geometry<Node<3>>;
    using IntegrationPointType = IntegrationPoint<3>;
    using IntegrationPointsArrayType = std::vector<IntegrationPointType>;

    template <typename ReturnType>
    using NodalFunctor = std::function<ReturnType(Node<3>&)>;

    template <typename ReturnType>
    using ElementFunctor = std::function<ReturnType(Element&)>;

    using NodalVectorFunctor = NodalFunctor<array_1d<double,3>>;
    using ElementVectorFunctor = ElementFunctor<array_1d<double,3>>;

    /// Pointer definition of DragUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DragUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    DragUtilities();

    /// Destructor.
    ~DragUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /** Perform a set of operations on the specified components in parallel and sum up the results
     *  @param rModelPart model part for providing the MPI communicator and miscellaneous info
     *  @param rComponents a container of nodes or elements to perform the operations on
     *  @param rOperations set of operations to be performed on each component
     *  @return tuple of objects holding the result of each operation, summed up for all components
     */
    template <typename ...TObjects, typename TModelPart, typename TComponentContainer, typename TOperationTuple>
    std::tuple<TObjects...> OperateAndReduceOnComponents(TModelPart&& rModelPart, TComponentContainer&& rComponents, TOperationTuple&& rOperations)
    {
        //using detail::operator+=;
        using TComponent = typename std::remove_reference<TComponentContainer>::type::value_type;
        using TObjectTuple = std::tuple<TObjects...>;

        TObjectTuple reduced_objects = block_for_each<SumReduction<TObjectTuple>>(std::forward<TComponentContainer>(rComponents),
            [&rOperations](TComponent& rComponent) {
                TObjectTuple objects;
                detail::RecurseTuple<std::tuple_size<TObjectTuple>::value-1>::OperateOnComponent(
                    rComponent,
                    objects,
                    std::forward<TOperationTuple>(rOperations)
                );
                return objects;
        });

        // Sum across processes with MPI
        detail::RecurseTuple<std::tuple_size<TObjectTuple>::value-1>::MPISumAll(
            std::forward<TModelPart>(rModelPart),
            reduced_objects
        );
        
        return reduced_objects;
    }

    /**
    * Computes the integral of the pressure stress term normal projection over the conditions
    * of the given modelpart
    * @param rModelPart reference to the model part in where the drag is to be computed
    * @return An array containing the drag force value.
    */
    array_1d<double, 3> CalculateBodyFittedDrag(ModelPart &rModelPart);

    /**
    * Computes the integral of the Cauchy stress term normal projection in the given modelpart elements.
    * @param rModelPart reference to the model part in where the drag is to be computed
    * @return An array containing the drag force value.
    */
    array_1d<double, 3> CalculateEmbeddedDrag(ModelPart &rModelPart);

    /**
    * Calculates the drag force location in embedded formulations
    * @param rModelPart reference to the model part in where the drag force location is to be computed
    * @return An array containing the drag force location coordinates.
    */
    array_1d<double, 3> CalculateEmbeddedDragCenter(const ModelPart &rModelPart);

    /** Perform a set of operations on all nodes of a model and sum up the results
     *  @param rModelPart model part containing the nodes to perform the operations on
     *  @param rOperations set of operations to perform on each node
     *  @note intended for python bindings
     */
    template <typename TOperationTuple, typename ...TObjects>
    std::tuple<TObjects...> SumNodalValues( ModelPart& rModelPart, TOperationTuple&& rOperations )
    {
        return OperateAndReduceOnComponents<TObjects...>(
            rModelPart,
            rModelPart.GetCommunicator().LocalMesh().Nodes(),
            std::forward<TOperationTuple>(rOperations)
        );
    }

    /** Perform a set of operations on all elements of a model and sum up the results
     *  @param rModelPart model part containing the elements to perform the operations on
     *  @param rOperations set of operations to perform on each element
     *  @note intended for python bindings
     */
    template <typename TOperationTuple, typename ...TObjects>
    std::tuple<TObjects...> SumElementValues( ModelPart& rModelPart, TOperationTuple&& rOperations )
    {
        return OperateAndReduceOnComponents<TObjects...>(
            rModelPart,
            rModelPart.Elements(),
            std::forward<TOperationTuple>(rOperations)
        );
    }

    /**
     * Create a functor that computes the moment around a reference point induced by the reactions on a node
     * @param rReference reference point for computing the moment
     * @return lambda function that computes the moment on a node
     */
    NodalVectorFunctor MakeComputeMomentOnNodeFunctor(const array_1d<double,3>& rReference) const;

    ///@}
    ///@name Access
    ///@{

    NodalVectorFunctor GetBodyFittedDragFunctor() const;

    ElementVectorFunctor GetEmbeddedDragFunctor() const;

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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

    ModelPart* mpModelPart;

    NodalVectorFunctor mComputeBodyFittedDragOnNode;

    ElementVectorFunctor mComputeEmbeddedDragOnElement;

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
    DragUtilities& operator=(DragUtilities const& rOther);

    /// Copy constructor.
    DragUtilities(DragUtilities const& rOther);

    ///@}

}; // Class DragUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DragUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DRAG_UTILITIES_H_INCLUDED  defined
