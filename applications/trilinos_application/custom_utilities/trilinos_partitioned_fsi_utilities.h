//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Źorrilla
//                   Jordi Cotela
//

#if !defined(KRATOS_TRILINOS_PARTITIONED_FSI_UTILITIES_H )
#define  KRATOS_TRILINOS_PARTITIONED_FSI_UTILITIES_H



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "../FSIapplication/custom_utilities/partitioned_fsi_utilities.hpp"

namespace Kratos
{
  ///@addtogroup TrilinosApplication
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

  /// Trilinos version of the partitioned FSI tools
  /** @see PartitionedFSIUtilities */
  template< class TSpace, unsigned int TDim >
  class TrilinosPartitionedFSIUtilities : public PartitionedFSIUtilities<TSpace,TDim>
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TrilinosPartitionedFSIUtilities
      KRATOS_CLASS_POINTER_DEFINITION(TrilinosPartitionedFSIUtilities);

      typedef typename TSpace::VectorType                     VectorType;
      typedef typename TSpace::MatrixType                     MatrixType;

      typedef typename TSpace::VectorPointerType              VectorPointerType;
      typedef typename TSpace::MatrixPointerType              MatrixPointerType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      TrilinosPartitionedFSIUtilities(){}

      /// Destructor.
      virtual ~TrilinosPartitionedFSIUtilities() override {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      virtual void SetUpInterfaceVector(const ModelPart& rInterfaceModelPart,
                                        VectorPointerType& pInterfaceVector) override
      {

      }

      virtual void ComputeInterfaceVectorResidual(ModelPart& rInterfaceModelPart,
                                                  const Variable<array_1d<double, 3 > >& rOriginalVariable,
                                                  const Variable<array_1d<double, 3 > >& rModifiedVariable,
                                                  VectorType& interface_residual) override
      {
          TSpace::SetToZero(interface_residual);

          // Compute node-by-node residual
          this->ComputeNodeByNodeResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, FSI_INTERFACE_RESIDUAL);

          // Compute consitent residual
          // this->ComputeConsistentResidual(rInterfaceModelPart, rOriginalVariable, rModifiedVariable, FSI_INTERFACE_RESIDUAL);

          // Assemble the final consistent residual values
          #pragma omp parallel for
          for(int k=0; k<static_cast<int>(rInterfaceModelPart.NumberOfNodes()); ++k)
          {
              const ModelPart::NodeIterator it_node = rInterfaceModelPart.NodesBegin()+k;
              const unsigned int base_i = k*TDim;

              const array_1d<double,3>& fsi_res = it_node->FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL);
              for (unsigned int jj=0; jj<TDim; ++jj)
              {
                  //interface_residual[base_i+jj] = fsi_res[jj];
              }
          }

          // Store the L2 norm of the error in the fluid process info
          rInterfaceModelPart.GetProcessInfo().GetValue(FSI_INTERFACE_RESIDUAL_NORM) = TSpace::TwoNorm(interface_residual);

      }

      virtual void UpdateInterfaceValues(ModelPart& rInterfaceModelPart,
                                         const Variable<array_1d<double, 3 > >& rSolutionVariable,
                                         VectorType& rCorrectedGuess) override
      {}

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
      TrilinosPartitionedFSIUtilities& operator=(TrilinosPartitionedFSIUtilities const& rOther){}

      /// Copy constructor.
      TrilinosPartitionedFSIUtilities(TrilinosPartitionedFSIUtilities const& rOther){}


      ///@}

    }; // Class TrilinosPartitionedFSIUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_PARTITIONED_FSI_UTILITIES_H  defined
