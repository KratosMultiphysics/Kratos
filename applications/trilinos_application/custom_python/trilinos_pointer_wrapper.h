//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					     Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED)
#define KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED

// System includes

// External includes

//Trilinos includes
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_application.h"
#include "trilinos_space.h"

namespace Kratos
{

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;

class AuxiliaryMatrixWrapper
{
  public:
    AuxiliaryMatrixWrapper(TrilinosSparseSpaceType::MatrixPointerType p) : mp(p){};

    TrilinosSparseSpaceType::MatrixPointerType& GetPointer() { return mp; }
    TrilinosSparseSpaceType::MatrixType& GetReference() { return *mp; }

  private:
    TrilinosSparseSpaceType::MatrixPointerType mp;
};

class AuxiliaryVectorWrapper
{
  public:
    AuxiliaryVectorWrapper(TrilinosSparseSpaceType::VectorPointerType p) : mp(p){};

    TrilinosSparseSpaceType::VectorPointerType& GetPointer() { return mp; }
    TrilinosSparseSpaceType::VectorType& GetReference() { return *mp; }

  private:
    TrilinosSparseSpaceType::VectorPointerType mp;
};

} // namespace Kratos.

#endif // KRATOS_TRILINOS_POINTER_WRAPPER_H_INCLUDED  defined