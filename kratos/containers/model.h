//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
* @class Model
* @ingroup KratosCore
* @brief This class aims to manage different model parts across multi-physics simulations
* @details The class behaves as a manager of the different model parts. It uses unordered_maps of the variables and the model parts for that purpose
* @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) Model final
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the index type
    using IndexType = ModelPart::IndexType;

    /// Pointer definition of Model
    KRATOS_CLASS_POINTER_DEFINITION(Model);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Model(){};

    /// Destructor.
    ~Model()
    {
        mRootModelPartMap.clear();
    }

    Model(const Model&) = delete;

    ///@}
    ///@name Operators
    ///@{

    Model & operator=(const Model&) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method clears the database of modelparts
     * @details Executes a clear on the model part map
     */
    void Reset();

    /**
     * @brief This method creates a new model part contained in the current Model with a given name and buffer size
     * @param ModelPartName The name of the new model part to be created
     * @param NewBufferSize The size of the buffer of the new model part created
     */
    ModelPart& CreateModelPart( const std::string& ModelPartName, IndexType NewBufferSize=1 );

    /**
     * @brief This method deletes a modelpart with a given name
     * @details Raises a warning in case the model part does not exists
     * @param ModelPartName The name of the model part to be removed
     */
    void DeleteModelPart( const std::string& ModelPartName );

    /**
     * @brief This method returns a model part given a certain name
     * @details Iterates over the list of submodelparts of the root model part
     * @param rFullModelPartName The name of the model part to be returned
     * @return Reference to the model part of interest
     */
    ModelPart& GetModelPart(const std::string& rFullModelPartName);

    /**
     * @brief This method returns a model part given a certain name
     * @details Iterates over the list of submodelparts of the root model part
     * @param rFullModelPartName The name of the model part to be returned
     * @return Reference to the model part of interest
     */
    const ModelPart& GetModelPart(const std::string& rFullModelPartName) const;

    /**
     * @brief This method checks if a certain a model part exists given a certain name
     * @details Iterates over the list of submodelparts of the root model part
     * @param rFullModelPartName The name of the model part to be checked
     * @return True if the model part exists, false otherwise
     */
    bool HasModelPart(const std::string& rFullModelPartName) const;

    /**
     * @brief This returns a vector containing a list of model parts names contained on the model
     * @details Iterates over the list of submodelparts of the root model part
     * @return A vector of strings containing the model parts names
     */
    std::vector<std::string> GetModelPartNames() const;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the data value container of the model
     * @details This method is used to access the data value container of the model.
     * @return The data value container of the model
     */
    DataValueContainer& GetDataValueContainer();

    /**
     * @brief This method returns the data value container of the model
     * @details This method is used to access the data value container of the model.
     * @return The data value container of the model
     */
    const DataValueContainer& GetDataValueContainer() const;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /// The map containing the list of model parts
    std::map<std::string, std::unique_ptr<ModelPart>> mRootModelPartMap;

    /// The data value container of the model
    DataValueContainer mDataValueContainer;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method searchs recursively a sub model part in a model part
     * @param rModelPartName The name to be search
     * @param pModelPart Pointer of the model part where search recursively
     * @return The pointer of the model part of interest
     */
    ModelPart* RecursiveSearchByName(const std::string& rModelPartName, ModelPart* pModelPart) const;

    /**
     * @brief This method splits the name of the model part using "." to define the hierarchy
     * @param rFullModelPartName The name with the full hierarchy
     * @return The vector containing each part of the name defining the model part hierarchy
     */
    std::vector<std::string> SplitSubModelPartHierarchy(const std::string& rFullModelPartName) const;

    /**
     * @brief This method creates a new model part contained in the current Model with a given name and buffer size
     * @param ModelPartName The name of the new model part to be created
     * @param NewBufferSize The size of the buffer of the new model part created
     */
    void CreateRootModelPart(const std::string& ModelPartName, ModelPart::IndexType NewBufferSize);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//       Model& operator=(Model const& rOther);

    /// Copy constructor.
//       Model(Model const& rOther);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

    ///@}

}; // Class Model

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                Model& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const Model& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.