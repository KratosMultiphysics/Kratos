#include "custom_io/hdf5_model_part_io.h"

#include "custom_io/hdf5_points_data.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"

namespace Kratos
{
namespace HDF5
{
ModelPartIO::ModelPartIO(File::Pointer pFile, std::string const& rPrefix)
: ModelPartIOBase(pFile, rPrefix)
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}



} // namespace HDF5.
} // namespace Kratos.
