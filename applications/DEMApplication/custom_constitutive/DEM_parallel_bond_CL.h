/////////////////////////////////////////////////
// Authors: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: July 2022
/////////////////////////////////////////////////

#if !defined(DEM_PARALLEL_BOND_CL_H_INCLUDE)
#define DEM_PARALLEL_BOND_CL_H_INCLUDE

// Project includes
#include "DEM_continuum_constitutive_law.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_parallel_bond : public DEMContinuumConstitutiveLaw {

        typedef DEMContinuumConstitutiveLaw BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_parallel_bond);

        DEM_parallel_bond() {}

        ~DEM_parallel_bond() {}

    private:

    };
} // namespace Kratos

#endif /*DEM_PARALLEL_BOND_CL_H_INCLUDE defined*/