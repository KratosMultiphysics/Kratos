#include "hdf5_io.h"

#include <sstream>

namespace Kratos
{

HDF5IO::HDF5IO():
  IO()
{

}

HDF5IO::~HDF5IO()
{

}

/// Turn back information as a string.
std::string HDF5IO::Info() const
{
  std::stringstream msg;
  msg << "HDF5IO" << std::endl;
  return msg.str();
}

/// Print information about this object.
void HDF5IO::PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.
void HDF5IO::PrintData(std::ostream& rOStream) const
{

}

}
