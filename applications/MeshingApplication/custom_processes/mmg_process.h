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

#if !defined(KRATOS_MMG_PROCESS)
#define KRATOS_MMG_PROCESS

// System includes
#include <unordered_set>
#include <unordered_map>

// External includes
// The includes related with the MMG library
// #include "mmg/libmmg.h"

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/mmg_utilities.h"
#include "containers/variables_list.h"
#include "meshing_application.h"

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
 * @class MmgProcess
 * @ingroup MeshingApplication
 * @brief This class is a remesher which uses the MMG library
 * @details This class is a remesher which uses the MMG library. The class uses a class for the 2D and 3D cases.
 * The remesher keeps the previous submodelparts and interpolates the nodal values between the old and new mesh
 * @author Vicente Mataix Ferrandiz
 */
template<MMGLibrary TMMGLibrary>
class KRATOS_API(MESHING_APPLICATION) MmgProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MmgProcess
    KRATOS_CLASS_POINTER_DEFINITION(MmgProcess);

    /// Node definition
    typedef Node <3>                                                   NodeType;
    // Geometry definition
    typedef Geometry<NodeType>                                     GeometryType;

    /// Conditions array size
    static constexpr SizeType Dimension = (TMMGLibrary == MMGLibrary::MMG2D) ? 2 : 3;

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
    MmgProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~MmgProcess() override = default;

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
    void OutputMdpa();

    /**
     * @brief Ths function removes superfluous (defined by "not belonging to an element") nodes from the model part
     */
    void CleanSuperfluousNodes();

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
        return "MmgProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MmgProcess";
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

    ModelPart& mrThisModelPart;                                      /// The model part to compute
    Parameters mThisParameters;                                      /// The parameters (can be used for general pourposes)
    NodeType::DofsContainerType mDofs;                               /// Storage for the dof of the node

    MmgUtilities<TMMGLibrary> mMmmgUtilities;                        /// The MMG utilities class

    std::string mFilename;                                           /// I/O file name
    IndexType mEchoLevel;                                            /// The echo level

    FrameworkEulerLagrange mFramework;                               /// The framework

    DiscretizationOption mDiscretization;                            /// The discretization option
    bool mRemoveRegions;                                             /// Cuttig-out specified regions during surface remeshing

    std::unordered_map<IndexType,std::vector<std::string>> mColors;  /// Where the sub model parts IDs are stored

    std::unordered_map<IndexType,Element::Pointer>   mpRefElement;   /// Reference element
    std::unordered_map<IndexType,Condition::Pointer> mpRefCondition; /// Reference condition

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This converts the framework string to an enum
     * @param rString The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */
    static inline FrameworkEulerLagrange ConvertFramework(const std::string& rString)
    {
        if(rString == "Lagrangian" || rString == "LAGRANGIAN")
            return FrameworkEulerLagrange::LAGRANGIAN;
        else if(rString == "Eulerian" || rString == "EULERIAN")
            return FrameworkEulerLagrange::EULERIAN;
        else if(rString == "ALE")
            return FrameworkEulerLagrange::ALE;
        else
            return FrameworkEulerLagrange::EULERIAN;
    }

    /**
     * @brief This converts the discretization string to an enum
     * @param rString The string
     * @return DiscretizationOption: The equivalent enum
     */
    static inline DiscretizationOption ConvertDiscretization(const std::string& rString)
    {
        if(rString == "Lagrangian" || rString == "LAGRANGIAN")
            return DiscretizationOption::LAGRANGIAN;
        else if(rString == "Standard" || rString == "STANDARD")
            return DiscretizationOption::STANDARD;
        else if(rString == "Isosurface" || rString == "ISOSURFACE" || rString == "IsoSurface")
            return DiscretizationOption::ISOSURFACE;
        else
            return DiscretizationOption::STANDARD;
    }

    /**
     * @brief This function generates the mesh MMG5 structure from a Kratos Model Part
     */
    void InitializeMeshData();

    /**
     *@brief This function generates the metric MMG5 structure from a Kratos Model Part
     */
    void InitializeSolDataMetric();

    /**
     *@brief This function generates the MMG5 structure for the distance field from a Kratos Model Part
     */
    void InitializeSolDataDistance();

    /**
     *@brief This function generates the displacement MMG5 structure from a Kratos Model Part
     */
    void InitializeDisplacementData();

    /**
     * @brief We execute the MMg library and build the new model part from the old model part
     */
    void ExecuteRemeshing();

    /**
     * @brief After we have transfer the information from the previous modelpart we initilize the elements and conditions
     */
    void InitializeElementsAndConditions();

    /**
     * @brief It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param PostOutput If the file to save is after or before remeshing
     */
    void SaveSolutionToFile(const bool PostOutput);

    /**
     * @brief It frees the memory used during all the process
     */
    void FreeMemory();

    /**
     * @brief It sets to zero the entity data, using the variables from the orginal model part
     * @param rNewModelPart The new container
     * @param rOldModelPart The old container
     * @tparam TContainerType The container type
     * @todo Interpolate values in the future
     */
    template<class TContainerType>
    void SetToZeroEntityData(
        TContainerType& rNewContainer,
        const TContainerType& rOldContainer
        )
    {
        // Firts we generate the variable list
        std::unordered_set<std::string> list_variables;
        const auto it_begin_old = rOldContainer.begin();
        auto& data = it_begin_old->Data();
        for(auto i = data.begin() ; i != data.end() ; ++i) {
            list_variables.insert((i->first)->Name());
        }

        for (auto& var_name : list_variables) {
            if (KratosComponents<Variable<bool>>::Has(var_name)) {
                const Variable<bool>& r_var = KratosComponents<Variable<bool>>::Get(var_name);
                VariableUtils().SetNonHistoricalVariable(r_var, false, rNewContainer);
            } else if (KratosComponents<Variable<double>>::Has(var_name)) {
                const Variable<double>& r_var = KratosComponents<Variable<double>>::Get(var_name);
                VariableUtils().SetNonHistoricalVariable(r_var, 0.0, rNewContainer);
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(var_name)) {
                const Variable<array_1d<double, 3>>& r_var = KratosComponents<Variable<array_1d<double, 3>>>::Get(var_name);
                const array_1d<double, 3> aux_value = ZeroVector(3);
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
                const Variable<array_1d<double, 4>>& r_var = KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name);
                const array_1d<double, 4> aux_value = ZeroVector(4);
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(var_name)) {
                const Variable<array_1d<double, 6>>& r_var = KratosComponents<Variable<array_1d<double, 6>>>::Get(var_name);
                const array_1d<double, 6> aux_value = ZeroVector(6);
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(var_name)) {
                const Variable<array_1d<double, 9>>& r_var = KratosComponents<Variable<array_1d<double, 9>>>::Get(var_name);
                const array_1d<double, 9> aux_value = ZeroVector(9);
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            } else if (KratosComponents<Variable<Vector>>::Has(var_name)) {
                const Variable<Vector>& r_var = KratosComponents<Variable<Vector>>::Get(var_name);
                Vector aux_value = ZeroVector(it_begin_old->GetValue(r_var).size());
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            } else if (KratosComponents<Variable<Matrix>>::Has(var_name)) {
                const Variable<Matrix>& r_var = KratosComponents<Variable<Matrix>>::Get(var_name);
                const Matrix& ref_matrix = it_begin_old->GetValue(r_var);
                Matrix aux_value = ZeroMatrix(ref_matrix.size1(), ref_matrix.size2());
                VariableUtils().SetNonHistoricalVariable(r_var, aux_value, rNewContainer);
            }
        }
    }

    /**
     * @brief This method collapses the prisms elements into triangles
     */
    void CollapsePrismsToTriangles();

    /**
     * @brief This method extrudes the triangles elements into prisms
     * @param rOldModelPart The old model part
     */
    void ExtrudeTrianglestoPrisms(ModelPart& rOldModelPart);

    /**
     * @brief This function removes the conditions with duplicated geometries
     */
    void ClearConditionsDuplicatedGeometries();

    /**
     * @brief This function creates an before/after remesh output file
     * @param rOldModelPart The old model part before remesh
     */
    void CreateDebugPrePostRemeshOutput(ModelPart& rOldModelPart);

    /**
     * @brief This method is used in order to mark the conditions in a recursive way to avoid remove necessary conditions
     * @param rModelPart The modelpart to be marked
     */
    void MarkConditionsSubmodelParts(ModelPart& rModelPart);

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();

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
    MmgProcess& operator=(MmgProcess const& rOther);

    /// Copy constructor.
    MmgProcess(MmgProcess const& rOther);

    ///@}

};// class MmgProcess
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<MMGLibrary TMMGLibrary>
inline std::istream& operator >> (std::istream& rIStream,
                                  MmgProcess<TMMGLibrary>& rThis);

/// output stream function
template<MMGLibrary TMMGLibrary>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MmgProcess<TMMGLibrary>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}// namespace Kratos.
#endif /* KRATOS_MMG_PROCESS defined */
