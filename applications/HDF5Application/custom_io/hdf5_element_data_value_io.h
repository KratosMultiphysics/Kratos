//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_HDF5_ELEMENT_DATA_VALUE_IO_H_INCLUDED)
#define KRATOS_HDF5_ELEMENT_DATA_VALUE_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{

class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of element data in HDF5.
class ElementDataValueIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ElementDataValueIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ElementDataValueIO(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    void WriteElementResults(ElementsContainerType const& rElements);

    void ReadElementResults(ElementsContainerType& rElements);

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{
    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<std::string> mVariableNames;
    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // class ElementDataValueIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_ELEMENT_DATA_VALUE_IO_H_INCLUDED defined
