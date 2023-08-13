#include "hdf5_file_serial.h"

#include <utility>
#include <vector>
#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "input_output/logger.h"
#include "utilities/data_type_traits.h"
#include "custom_utilities/h5_data_type_traits.h"

namespace Kratos
{
namespace HDF5
{
FileSerial::FileSerial(Parameters& rSettings) : File(rSettings)
{
}

FileSerial::FileSerial(
    const DataCommunicator& rDataCommunicator,
    Parameters Settings)
    : File(rDataCommunicator, Settings)
{
}

FileSerial::FileSerial(FileSerial&& rOther) : File(std::move(rOther))
{
}

FileSerial& FileSerial::operator=(FileSerial&& rOther)
{
    File::operator=(std::move(rOther));
    return *this;
}

} // namespace HDF5.
} // namespace Kratos.