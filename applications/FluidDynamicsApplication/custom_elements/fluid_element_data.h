//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED )
#define  KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"

#include "fluid_dynamics_application_variables.h"


namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
  ///@{

  

  /// Auxiliary class to hold data for elements based on FluidElement
  /** TODO: see how this can work for a generic number of stored arguments.
   */
  template< unsigned int TNumNodes >
  class FluidElementData
    {
    public:
    
      ///@name Public enum
      ///@{

      enum NodalData {
          velocity_x,
          velocity_y,
          velocity_z,
          pressure,
          density,
          viscosity,
          mesh_velocity_x,
          mesh_velocity_y,
          mesh_velocity_z,
          size
      };

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      FluidElementData(Geometry< Node<3> > &rGeom);

      /// Destructor.
      ~FluidElementData();

      ///@}
      ///@name Access
      ///@{

      double GetValue(NodalData Value, unsigned int Node);

      void Evaluate(const Kratos::Vector& rShapeFunctions, NodalData Value, double &rOutput);

      void EvaluateGradient(const Kratos::Matrix& rShapeFunctionGradients, NodalData Value, array_1d<double,3> &rGradient);

      ///@}

    private:
      ///@name Member Variables
      ///@{

      boost::numeric::ublas::bounded_matrix< double, size, TNumNodes > mNodalData;

      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      FluidElementData& operator=(FluidElementData const& rOther);

      /// Copy constructor.
      FluidElementData(FluidElementData const& rOther);


      ///@}

    }; // Class FluidElementData

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLUID_ELEMENT_DATA_H_INCLUDED  defined
