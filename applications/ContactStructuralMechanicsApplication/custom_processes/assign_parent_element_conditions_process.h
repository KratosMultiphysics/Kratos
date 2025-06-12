// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"

namespace Kratos
{
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

/**
 * @class AssignParentElementConditionsProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process initializes the variables related with the ALM
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AssignParentElementConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignParentElementConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignParentElementConditionsProcess);

    /// Definition of the vector indexes considered
    using VectorIndexType = std::vector<std::size_t>;

    /// Definition of the hasher considered
    using VectorIndexHasherType = VectorIndexHasher<VectorIndexType>;

    /// Definition of the key comparor considered
    using VectorIndexComparorType = VectorIndexComparor<VectorIndexType>;

    /// Define the map considered for face ids
    using HashMapVectorIntType = std::unordered_map<VectorIndexType, std::size_t, VectorIndexHasherType, VectorIndexComparorType>;

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    /// Nodes array type definition
    using NodesArrayType = ModelPart::NodesContainerType;

    /// Conditions array type definition
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor that takes the conditions and elements model parts.
     * @param rConditionsModelPart The conditions model part.
     * @param rElementsModelPart The elements model part.
     */
    AssignParentElementConditionsProcess( 
        ModelPart& rConditionsModelPart,
        ModelPart& rElementsModelPart
        )
        : mrConditionsModelPart(rConditionsModelPart),
          mrElementsModelPart(rElementsModelPart)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /**
     * @brief Constructs an AssignParentElementConditionsProcess object.
     * @param rModel The model containing the model parts.
     * @param ThisParameters The parameters defining the process.
     */
    AssignParentElementConditionsProcess(
        Model& rModel,
        Parameters ThisParameters
        ) : mrConditionsModelPart(rModel.GetModelPart(ThisParameters["conditions_model_part_name"].GetString())),
            mrElementsModelPart(rModel.GetModelPart(ThisParameters["elements_model_part_name"].GetString()))
    {
        KRATOS_TRY;

        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mEchoLevel = ThisParameters["echo_level"].GetInt();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~AssignParentElementConditionsProcess() override = default;

    /// Copy constructor.
    AssignParentElementConditionsProcess( AssignParentElementConditionsProcess const& rOther) = delete;

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

    /**
     * @brief Executes the function.
     */
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Executes the function.
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations  right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

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
        return "AssignParentElementConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignParentElementConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrConditionsModelPart; /// The model part of conditions
    ModelPart& mrElementsModelPart;   /// The model part of elements
    int mEchoLevel;                   /// The echo level
    HashMapVectorIntType mMapFaceIds; /// The map containing the face ids of the elements

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
    AssignParentElementConditionsProcess& operator=(AssignParentElementConditionsProcess const& rOther) = delete;

    /// Copy constructor.
    //AssignParentElementConditionsProcess(AssignParentElementConditionsProcess const& rOther);

    ///@}

}; // Class AssignParentElementConditionsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignParentElementConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignParentElementConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.