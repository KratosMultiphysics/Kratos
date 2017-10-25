//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED)
#define KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED

// System includes
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/io.h"

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
///@addtogroup HDF5Application
///@{

class HDF5ModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HDF5ModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    bool ReadNodes(NodesContainerType& rNodes) override;

    std::size_t ReadNodesNumber() override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    void ReadElements(NodesContainerType& rNodes,
                      PropertiesContainerType& rProperties,
                      ElementsContainerType& rElements) override;

    std::size_t ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities) override;

    void WriteElements(ElementsContainerType const& rElements) override;

    void ReadConditions(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions) override;

    std::size_t ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities) override;

    void ReadInitialValues(ModelPart& rModelPart) override;

    void ReadInitialValues(NodesContainerType& rNodes,
                           ElementsContainerType& rElements,
                           ConditionsContainerType& rConditions) override;

    void ReadModelPart(ModelPart& rModelPart) override;

    void WriteModelPart(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    unsigned GetPID() const;

    unsigned GetTotalProcesses() const;
    ///@}

private:
    ///@name Member Variables
    ///@{
    HDF5File::Pointer mpFile;
    unsigned m_pid;
    unsigned m_total_processes;
    ///@}

    ///@name Private Operations
    ///@{
        /// Write all nodes into a single array.
        void WriteNodesSerial(NodesContainerType const& rNodes) const;

        /// Write nodes into local and ghost arrays.
        void WriteNodesParallel(NodesContainerType const& rNodes) const;
        
        /// Read nodes from a single array.
        void ReadNodesSerial(NodesContainerType const& rNodes) const;
        
        /// Read nodes from local and ghost arrays.
        void ReadNodesParallel(NodesContainerType const& rNodes) const;
                
    ///@}
};

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED defined
