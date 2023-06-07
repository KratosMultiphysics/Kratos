// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

/* Geometries defined  */
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @class ShellToSolidShellProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This method transforms triangular and quadrilateral elements into prisms and hexahedra elements
 * @details It used the nodal normal for that pourpose
 * @tparam TNumNodes To distinghuis between triangles and quadrilaterals
 * @author Vicente Mataix Ferrandiz
*/
template<SizeType TNumNodes>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ShellToSolidShellProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShellToSolidShellProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShellToSolidShellProcess);

    /// The index definition
    typedef std::size_t                                     IndexType;

    /// Geometric type definitions
    typedef Node                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;

    /// The definition of the containers
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    ShellToSolidShellProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ShellToSolidShellProcess() override
    = default;

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

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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
        return "ShellToSolidShellProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ShellToSolidShellProcess";
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

    ModelPart& mrThisModelPart;              /// The main model part
    Parameters mThisParameters;              /// The parameters (can be used for general pourposes)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     * @param ReorderAccordingShellConnectivity True if we reorder using the shell connectivity as reference
     */
    void ReorderAllIds(const bool ReorderAccordingShellConnectivity = false);

    /**
     * @brief This is the execution to convert a 2D element into 3D
     */
    void ExecuteExtrusion();

    /**
     * @brief This is the execution to convert a 3D element into 2D
     */
    void ExecuteCollapse();

    /**
     * @brief This replaces the previous geometry
     * @param rGeometryModelPart The previous model part
     * @param rAuxiliaryModelPart The new created model part
     */
    void ReplacePreviousGeometry(
        ModelPart& rGeometryModelPart,
        ModelPart& rAuxiliaryModelPart
        );

    /**
     * @brief If we reassign the constitutive laws
     * @param rGeometryModelPart The previous model part
     * @param rSetIdProperties The set containing the properties ids
     */
    void ReassignConstitutiveLaw(
        ModelPart& rGeometryModelPart,
        std::unordered_set<IndexType>& rSetIdProperties
        );

    /**
     * @brief After we have transfer the information from the previous modelpart we initilize the elements
     */
    void InitializeElements();

    /**
     * @brief After we have created the new elements we export to a *.mdpa file
     */
    void ExportToMDPA();

    /**
     * @brief After we have created the new elements we delete the auxiliary model parts
     */
    void CleanModel();

    /**
     * @brief It computes the mean of the normal on the elements in all the nodes (non historical version)
     */
    inline void ComputeNodesMeanNormalModelPartNonHistorical();

    /**
     * @brief It copies the variable list to the node
     * @param pNodeNew The new node pointer
     * @param pNodeNew The old node pointer
     */
    inline void CopyVariablesList(
        NodeType::Pointer pNodeNew,
        NodeType::Pointer pNodeOld
        );

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
    ShellToSolidShellProcess& operator=(ShellToSolidShellProcess const& rOther) = delete;

    /// Copy constructor.
    //ShellToSolidShellProcess(ShellToSolidShellProcess const& rOther);


    ///@}

}; // Class ShellToSolidShellProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ShellToSolidShellProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ShellToSolidShellProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
