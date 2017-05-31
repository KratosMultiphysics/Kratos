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

      virtual void SetUpInterfaceVector(ModelPart& rInterfaceModelPart,
                                        VectorPointerType& pInterfaceVector) override
      {
          // Initialize Epetra_Map for the vector
          int NumLocalInterfaceDofs = rInterfaceModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * TDim;
          int NumGlobalInterfaceDofs = NumLocalInterfaceDofs;
          rInterfaceModelPart.GetCommunicator().SumAll(NumGlobalInterfaceDofs);
          int IndexBase = 0; // 0 for C-style vectors, 1 for Fortran numbering
          Epetra_Map InterfaceMap(NumGlobalInterfaceDofs,NumLocalInterfaceDofs,IndexBase,mrEpetraComm);

          // Create new vector using given map
          VectorPointerType Temp = VectorPointerType(new Epetra_FEVector(InterfaceMap));
          pInterfaceVector.swap(Temp);

          // Set interface vector to zero
          pInterfaceVector->PutScalar(0.0);

          std::cout << rInterfaceModelPart.GetCommunicator().MyPID() << ": ";
          KRATOS_WATCH(InterfaceMap.NumGlobalPoints())
          KRATOS_WATCH(InterfaceMap.NumMyPoints())
          KRATOS_WATCH(InterfaceMap.MinAllGID())
          KRATOS_WATCH(InterfaceMap.MaxAllGID())
          KRATOS_WATCH(InterfaceMap.MinMyGID())
          KRATOS_WATCH(InterfaceMap.MaxMyGID())
          KRATOS_WATCH(InterfaceMap.MinLID())
          KRATOS_WATCH(InterfaceMap.MaxLID())
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
          //KRATOS_ERROR_IF_NOT(rVector.Map().MyGID(LocalRow)) << " non-local id: " << LocalRow << ".";
          std::cout << "set: " << LocalRow << " " << Value << std::endl;
          int ierr = rVector.ReplaceMyValue(LocalRow,0,Value);
          KRATOS_ERROR_IF(ierr != 0) << "ReplaceMyValue returns " << ierr << " for row " << LocalRow << " of " << rVector.Map().MaxLID() << std::endl;
          //rVector.GlobalAssemble();
          //KRATOS_WATCH(rVector[0][LocalRow])
      }

      virtual double GetLocalValue(VectorType& rVector, int LocalRow) const override
      {
          std::cout << "get: " << LocalRow << " " << rVector[0][LocalRow] << std::endl;
          //KRATOS_WATCH(rVector[0][LocalRow])
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
