//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


// System includes
#include <pybind11/pybind11.h>

#if defined(KRATOS_PYTHON)
// External includes

//Trilinos includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;

class AuxiliaryMatrixWrapper
{
  public:
    AuxiliaryMatrixWrapper(TrilinosSparseSpaceType::MatrixPointerType p) : mp(p){};

    TrilinosSparseSpaceType::MatrixPointerType &GetPointer() { return mp; }
    TrilinosSparseSpaceType::MatrixType &GetReference() { return *mp; }

  private:
    TrilinosSparseSpaceType::MatrixPointerType mp;
};

class AuxiliaryVectorWrapper
{
  public:
    AuxiliaryVectorWrapper(TrilinosSparseSpaceType::VectorPointerType p) : mp(p){};

    TrilinosSparseSpaceType::VectorPointerType &GetPointer() { return mp; }
    TrilinosSparseSpaceType::VectorType &GetReference() { return *mp; }

  private:
    TrilinosSparseSpaceType::VectorPointerType mp;
};

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined