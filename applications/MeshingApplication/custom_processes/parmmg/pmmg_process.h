// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Marc Nunez
//                   Carlos Roig
//                   Vicente Mataix Ferrandizz
//

#if !defined(KRATOS_PMMG_PROCESS)
#define KRATOS_PMMG_PROCESS

// System includes
#include <unordered_set>
#include <unordered_map>

// External includes

// Project includes
#include "custom_processes/mmg/mmg_process.h"
#include "custom_utilities/parmmg/pmmg_utilities.h"

// NOTE: The following contains the license of the PMMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017- .
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
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
    typedef std::size_t                  IndexType;

    /// Size definition
    typedef std::size_t                   SizeType;

    /// Index vector
    typedef std::vector<IndexType> IndexVectorType;

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
 * @class ParMmgProcess
 * @ingroup MeshingApplication
 * @brief This class is a remesher which uses the PMMG library
 * @details This class is a remesher which uses the PMMG library. The class uses a class for the 2D and 3D cases.
 * The remesher keeps the previous submodelparts and interpolates the nodal values between the old and new mesh
 * @author Marc Nunez
 * @author Carlos Roig
 * @author Vicente Mataix Ferrandiz
 */
template<PMMGLibrary TPMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) ParMmgProcess
    : public MmgProcess<MMGLibrary::MMG3D>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParMmgProcess
    KRATOS_CLASS_POINTER_DEFINITION(ParMmgProcess);

    typedef MmgProcess<MMGLibrary::MMG3D> BaseType;

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Conditions array size
    static constexpr SizeType Dimension = 3;

    /// The type of array considered for the tensor
    typedef typename std::conditional<Dimension == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

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

    /**
     * @brief This is the default constructor, which is used to read the input files
     * @param rThisModelPart The model part
     * @param ThisParameters The parameters
     */
    ParMmgProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ParMmgProcess() override = default;

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

    void operator()();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being execute once before the solution loop but after all of the solvers where built
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function is designed for being execute once before the solution loop but after all of the solvers where built
     */
    void ExecuteBeforeSolutionLoop() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function will be executed at every time step BEFORE  writing the output
     */
    void ExecuteBeforeOutputStep() override;

    /**
     * @brief This function will be executed at every time step AFTER writing the output
     */
    void ExecuteAfterOutputStep() override;

    /**
     * @brief This function is designed for being called at the end of the computations right after reading the model and the groups
     */
    void ExecuteFinalize() override;

    /**
     * @brief This sets the output mesh in a .mdpa format
     */
    void OutputMdpa() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
    */
    const Parameters GetDefaultParameters() const override;

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
        return "ParMmgProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ParMmgProcess";
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

    /**
     * @brief This function generates the mesh MMG5 structure from a Kratos Model Part
     */
    void InitializeMeshData() override;

    /**
     *@brief This function generates the metric MMG5 structure from a Kratos Model Part
     */
    void InitializeSolDataMetric() override;

    /**
     * @brief We execute the MMg library and build the new model part from the old model part
     */
    void ExecuteRemeshing() override;

    /**
     * @brief After we have transfer the information from the previous modelpart we initilize the elements and conditions
     */
    void InitializeElementsAndConditions() override;

    /**
     * @brief It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param PostOutput If the file to save is after or before remeshing
     */
    void SaveSolutionToFile(const bool PostOutput) override;

    /**
     * @brief It frees the memory used during all the process
     */
    void FreeMemory() override;

    /**
     * @brief This function removes the conditions with duplicated geometries
     */
    void ClearConditionsDuplicatedGeometries() override;

    /**
     * @brief This function creates an before/after remesh output file
     * @param rOldModelPart The old model part before remesh
     */
    void CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart) override;

    template<typename TPointerType>
    void SyncMapAcrossRanks(std::unordered_map<IndexType, TPointerType>& rInputMap);


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

    ParMmgUtilities<TPMMGLibrary> mPMmgUtilities;                     /// The PMMG utilities class

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParMmgProcess& operator=(ParMmgProcess const& rOther);

    /// Copy constructor.
    ParMmgProcess(ParMmgProcess const& rOther);

    ///@}

};// class ParMmgProcess
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<PMMGLibrary TPMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  ParMmgProcess<TPMMGLibrary>& rThis);

/// output stream function
template<PMMGLibrary TPMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParMmgProcess<TPMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_PMMG_PROCESS defined */
