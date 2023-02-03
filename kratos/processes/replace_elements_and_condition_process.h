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
//                   Vicente Mataix Ferrandiz
//                   Philipp Bucher
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "geometries/geometry_data.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class ReplaceElementsAndConditionsProcess
 * @ingroup KratosCore
 * @brief This methods replaces elements and conditions in a model part by a given name
 * @details The submodelparts are later updated
 * @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) ReplaceElementsAndConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReplaceElementsAndConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsAndConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModel The model containing the model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    ReplaceElementsAndConditionsProcess(
        Model& rModel,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModel.GetModelPart(Settings["model_part_name"].GetString())),
            mSettings( Settings)
    {
        KRATOS_TRY

        // Initialize member variables
        InitializeMemberVariables();

        KRATOS_CATCH("")
    }

    /**
     * @brief Default constructor
     * @param rModelPart The model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    ReplaceElementsAndConditionsProcess(
        ModelPart& rModelPart,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModelPart),
            mSettings( Settings)
    {
        KRATOS_TRY

        // Initialize member variables
        InitializeMemberVariables();

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    ReplaceElementsAndConditionsProcess(ReplaceElementsAndConditionsProcess const& rOther) = delete;

    /// Destructor.
    ~ReplaceElementsAndConditionsProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ReplaceElementsAndConditionsProcess& operator=(ReplaceElementsAndConditionsProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override
    {
        return Kratos::make_shared<ReplaceElementsAndConditionsProcess>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the ReplaceElementsAndConditionsProcess algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"({
            "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
            "element_name"    : "PLEASE_CHOOSE_ELEMENT_NAME",
            "condition_name"  : "PLEASE_CHOOSE_CONDITION_NAME"
        } )" );
        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ReplaceElementsAndConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceElementsAndConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;                         /// The main model part where the elements and conditions will be replaced
    Parameters mSettings;                           /// The settings of the problem (names of the conditions and elements)
    std::array<bool, 2> mOnlyOneElementOrCondition; /// This array stores if only one element or condition type is going to be replaced

    /// Equivalent types definition
    const std::unordered_map<GeometryData::KratosGeometryType, std::string> mGeometryTypesToStrings = {
        {GeometryData::KratosGeometryType::Kratos_Point2D, "Point2D"},
        {GeometryData::KratosGeometryType::Kratos_Point3D, "Point3D"},
        {GeometryData::KratosGeometryType::Kratos_Line2D2, "Line2D2"},
        {GeometryData::KratosGeometryType::Kratos_Line2D3, "Line2D3"},
        {GeometryData::KratosGeometryType::Kratos_Line2D4, "Line2D4"},
        {GeometryData::KratosGeometryType::Kratos_Line2D5, "Line2D5"},
        {GeometryData::KratosGeometryType::Kratos_Line3D2, "Line3D2"},
        {GeometryData::KratosGeometryType::Kratos_Line3D3, "Line3D3"},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D3, "Triangle2D3"},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D6, "Triangle2D6"},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D10, "Triangle2D10"},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D15, "Triangle2D15"},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D3, "Triangle3D3"},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D6, "Triangle3D6"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, "Quadrilateral2D4"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, "Quadrilateral2D8"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, "Quadrilateral2D9"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, "Quadrilateral3D4"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, "Quadrilateral3D8"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, "Quadrilateral3D9"},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4, "Tetrahedra3D4"},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10, "Tetrahedra3D10"},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D8, "Hexahedra3D8"},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20, "Hexahedra3D20"},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D27, "Hexahedra3D27"},
        {GeometryData::KratosGeometryType::Kratos_Prism3D6, "Prism3D6"},
        {GeometryData::KratosGeometryType::Kratos_Prism3D15, "Prism3D15"},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D5, "Pyramid3D5"},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D13, "Pyramid3D13"}
    };
    
    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes the member variables
     */
    void InitializeMemberVariables();

    ///@}

}; // Class ReplaceElementsAndConditionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceElementsAndConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceElementsAndConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.