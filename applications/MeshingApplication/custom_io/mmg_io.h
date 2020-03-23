// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MMG_IO)
#define KRATOS_MMG_IO

// System includes

// External includes

// Project includes
#include "includes/io.h"
#include "includes/model_part.h"
#include "custom_utilities/mmg_utilities.h"

// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Index definition
    typedef std::size_t                                               IndexType;

    /// Size definition
    typedef std::size_t                                                SizeType;

    /// Index vector
    typedef std::vector<IndexType>                              IndexVectorType;

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
 * @class MmgIO
 * @ingroup MeshingApplication
 * @brief This class is a IO which uses the MMG library
 * @details This class is an IO tool for MMG .mesh files
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibrary TMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) MmgIO
    : public IO
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MmgIO
    KRATOS_CLASS_POINTER_DEFINITION(MmgIO);

    /// Node containers definition
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    /// Elements containers definition
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Mesh definition
    typedef Mesh<NodeType, Properties, Element, Condition>             MeshType;
    /// Properties container definition
    typedef MeshType::PropertiesContainerType           PropertiesContainerType;
    /// Nodes container definition
    typedef MeshType::NodeConstantIterator                 NodeConstantIterator;
    /// Conditions container definition
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    /// Elements container definition
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;

    /// Conditions array size
    static constexpr SizeType Dimension = (TMMGLibrary == MMGLibrary::MMG2D) ? 2 : 3;

    /// Conditions array size
    static constexpr SizeType ConditionsArraySize = (Dimension == 2) ? 1 : 2;

    /// Elements array size
    static constexpr SizeType ElementsArraySize = (Dimension == 2) ? 1 : 2;

    /// The type of array considered for the tensor
    typedef typename std::conditional<Dimension == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Double vector
    typedef std::vector<double> DoubleVectorType;

    /// Double vector map
    typedef std::unordered_map<DoubleVectorType, IndexType, KeyHasherRange<DoubleVectorType>, KeyComparorRange<DoubleVectorType> > DoubleVectorMapType;

    /// Index vector map
    typedef std::unordered_map<IndexVectorType, IndexType, KeyHasherRange<IndexVectorType>, KeyComparorRange<IndexVectorType> > IndexVectorMapType;

    /// Colors map
    typedef std::unordered_map<IndexType,IndexType> ColorsMapType;

    /// Index pair
    typedef std::pair<IndexType,IndexType> IndexPairType;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /// Constructor with filenames.
    MmgIO(
        std::string const& rFilename,
        Parameters ThisParameters = Parameters(R"({})"),
        const Flags Options = IO::READ | IO::IGNORE_VARIABLES_ERROR.AsFalse() | IO::SKIP_TIMER);

    /// Destructor.
    ~MmgIO() override = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This read the current stream file in order to write a model part
     */
    void ReadModelPart(ModelPart& rModelPart) override;

    /**
     * @brief This writes the current model part info a file
     */
    void WriteModelPart(ModelPart& rModelPart) override;

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
        return "MmgIO";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MmgIO";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    std::string mFilename;                       /// The name of the file
    Parameters mThisParameters;                  /// The parameters (can be used for general pourposes)
    Flags mOptions;                              /// Configuration flags

    MmgUtilities<TMMGLibrary> mMmmgUtilities;    /// The MMG utilities class

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "echo_level"                           : 0
        })" );

        return default_parameters;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

//     /// Assignment operator.
//     MmgIO& operator=(MmgIO const& rOther);

//     /// Copy constructor.
//     MmgIO(MmgIO const& rOther);

    ///@}

};// class MmgIO
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<MMGLibrary TMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  MmgIO<TMMGLibrary>& rThis);

/// output stream function
template<MMGLibrary TMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MmgIO<TMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_MMG_IO defined */
