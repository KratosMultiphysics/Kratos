// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <filesystem>

// External includes

// Project includes
#include "includes/io.h"


namespace Kratos {

///@addtogroup MedApplication
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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MED_APPLICATION) MedModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MedModelPartIO
    KRATOS_CLASS_POINTER_DEFINITION(MedModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with filename.
    MedModelPartIO(
        const std::filesystem::path& rFileName,
        const Flags Options = IO::READ);

    /// Copy constructor.
    MedModelPartIO(MedModelPartIO const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MedModelPartIO& operator=(MedModelPartIO const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method reads the model part
     * @param rThisModelPart The model part to be read
     */
    void ReadModelPart(ModelPart& rThisModelPart) override;

    /**
     * @brief This method writes the model part
     * @param rThisModelPart The model part to be written
     */
    void WriteModelPart(const ModelPart& rThisModelPart) override;

    /**
     * @brief This method divides a model part into partitions
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    void DivideInputToPartitions(SizeType NumberOfPartitions,
                                 const PartitioningInfo& rPartitioningInfo) override;

    /**
     * @brief This method divides a model part into partitions
     * @param pStreams The stream pointer
     * @param NumberOfPartitions The number of partitions
     * @param rPartitioningInfo Information about partitioning of entities
     */
    void DivideInputToPartitions(Kratos::shared_ptr<std::iostream> * pStreams,
                                SizeType NumberOfPartitions,
                                const PartitioningInfo& rPartitioningInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MedModelPartIO";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MedModelPartIO";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    std::filesystem::path mFileName;

    Flags mOptions;

    class MedFileHandler; // forward declared to avoid "med.h" include in header
    Kratos::shared_ptr<MedFileHandler> mpFileHandler;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}

}; // Class MedModelPartIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MedModelPartIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MedModelPartIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.
