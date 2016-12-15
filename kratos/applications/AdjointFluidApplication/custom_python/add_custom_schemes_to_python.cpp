// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "custom_python/add_custom_schemes_to_python.h"
#include "custom_strategies/custom_schemes/adjoint_bossak_drag_scheme.h"

namespace Kratos
{

namespace Python
{
  using namespace boost::python;

  void AddCustomSchemesToPython()
  {
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;

    class_< AdjointBossakDragScheme<SparseSpaceType, LocalSpaceType>,
    	    bases<SchemeType>,
    	    boost::noncopyable >
        ( "AdjointBossakDragScheme", init<double>() )
        ;

  }

} // namespace Python

} // namespace Kratos
