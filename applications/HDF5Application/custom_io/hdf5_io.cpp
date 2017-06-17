#include "hdf5_io.h"

#include <sstream>

namespace Kratos
{

HDF5IO::HDF5IO(std::string FileName, Flags Options):
  IO(),
  mFileName(FileName)
{

}

HDF5IO::~HDF5IO()
{

}

void HDF5IO::WriteModelPart(ModelPart& rModelPart)
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
