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
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

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
      TrilinosPartitionedFSIUtilities(const Epetra_MpiComm& EpetraCommunicator):
        mrEpetraComm(EpetraCommunicator)
      {}

      /// Destructor.
      virtual ~TrilinosPartitionedFSIUtilities() override {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      VectorPointerType SetUpInterfaceVector(ModelPart& rInterfaceModelPart) override
      {
          // Initialize Epetra_Map for the vector
          int NumLocalInterfaceDofs = rInterfaceModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * TDim;
          int NumGlobalInterfaceDofs = NumLocalInterfaceDofs;
          rInterfaceModelPart.GetCommunicator().SumAll(NumGlobalInterfaceDofs);
          int IndexBase = 0; // 0 for C-style vectors, 1 for Fortran numbering
          Epetra_Map InterfaceMap(NumGlobalInterfaceDofs,NumLocalInterfaceDofs,IndexBase,mrEpetraComm);

          // Create new vector using given map
          VectorPointerType p_int_vect = Kratos::make_shared<Epetra_FEVector>(InterfaceMap);

          // Set interface vector to zero
          p_int_vect->PutScalar(0.0);

          return p_int_vect;
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

      virtual void SetLocalValue(VectorType& rVector, int LocalRow, double Value) const override
      {
          // Note: We don't go through TSpace because TrilinosSpace forces a GlobalAssemble after each SetValue
          rVector.ReplaceMyValue(LocalRow,0,Value);
      }

      virtual double GetLocalValue(const VectorType& rVector, int LocalRow) const override
      {
          return rVector[0][LocalRow];
      }


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

      const Epetra_MpiComm& mrEpetraComm;


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
