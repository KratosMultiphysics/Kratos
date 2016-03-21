// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                      March 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/vertex_morphing_utilities.h"

// ==============================================================================

namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython()
  {
      using namespace boost::python;

      // ========================================================================
      // Optimization based on Vertex Morphing
      // ========================================================================
      class_<VertexMorphingUtilities, bases<Process> >("VertexMorphingUtilities", init<ModelPart&,
                                                                                       const int,
                                                                                       boost::python::dict,
                                                                                       boost::python::dict,
                                                                                       double,
                                                                                       const int>())

              // ================================================================
              // General geometrical operations
              // ================================================================
              .def("compute_unit_surface_normals", &VertexMorphingUtilities::compute_unit_surface_normals)
              .def("project_grad_on_unit_surface_normal", &VertexMorphingUtilities::project_grad_on_unit_surface_normal)

              // ================================================================
              // For perfoming Vertex Morphing
              // ================================================================
              .def("filter_gradients", &VertexMorphingUtilities::filter_gradients)

              // ================================================================
              // General optimization operations
              // ================================================================
              .def("update_design_variable", &VertexMorphingUtilities::update_design_variable)
              .def("update_shape", &VertexMorphingUtilities::update_shape)

              // ================================================================
              // For running unconstrained descent methods
              // ================================================================
              .def("compute_search_direction_steepest_descent", &VertexMorphingUtilities::compute_search_direction_steepest_descent)

              // ================================================================
              // For running augmented Lagrange method
              // ================================================================
              .def("initialize_augmented_lagrange", &VertexMorphingUtilities::initialize_augmented_lagrange)
              .def("compute_search_direction_augmented_lagrange", &VertexMorphingUtilities::compute_search_direction_augmented_lagrange)
              .def("udpate_augmented_lagrange_parameters", &VertexMorphingUtilities::udpate_augmented_lagrange_parameters)
              .def("get_penalty_fac", &VertexMorphingUtilities::get_penalty_fac)
              .def("get_lambda", &VertexMorphingUtilities::get_lambda)
              .def("get_value_of_augmented_lagrangian", &VertexMorphingUtilities::get_value_of_augmented_lagrangian)

              // ================================================================
              // For running penalized projection method
              // ================================================================
              .def("compute_search_direction_penalized_projection", &VertexMorphingUtilities::compute_search_direction_penalized_projection)
              ;
  }


}  // namespace Python.

} // Namespace Kratos

