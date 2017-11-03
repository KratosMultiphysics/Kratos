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
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

/// A class for serial IO of a model part in HDF5.
class ModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ModelPartIO(Parameters& rParams, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    bool ReadNodes(NodesContainerType& rNodes) override;

    std::size_t ReadNodesNumber() override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    void ReadProperties(PropertiesContainerType& rThisProperties) override;

    void WriteProperties(Properties const& rThisProperties) override;

    void WriteProperties(PropertiesContainerType const& rThisProperties) override;

    void ReadElements(NodesContainerType& rNodes,
                      PropertiesContainerType& rProperties,
                      ElementsContainerType& rElements) override;

    void WriteElements(ElementsContainerType const& rElements) override;

    void ReadConditions(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions) override;

    void WriteConditions(ConditionsContainerType const& rConditions) override;

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

    File& GetFile() const;

    virtual void Check();

    ///@}

private:
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<std::string> mElementNames;
    std::vector<const Element*> mElementPointers;
    std::vector<std::string> mConditionNames;
    std::vector<const Condition*> mConditionPointers;

    ///@}

    ///@name Private Operations
    ///@{

    ///@}
};

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED defined
