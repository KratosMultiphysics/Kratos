//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_READ_MATERIALS_H_INCLUDED )
#define  KRATOS_READ_MATERIALS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

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
 * @class ReadMaterialsUtility
 * @ingroup KratosCore
 * @brief Process to read constitutive law and material properties from a json file
 * @details This process reads constitutive law and material properties from a json file
 * and assign them to elements and conditions.
 * The definition includes the creation of subproperties
 * @author Marcelo Raschi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) ReadMaterialsUtility
{
  public:

    ///@name Type Definitions
    ///@{

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the mesh id (always zero)
    static constexpr IndexType mesh_id = 0;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ReadMaterialProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReadMaterialsUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModel The model containing the problem to solve
     */
    ReadMaterialsUtility(Model& rModel) : mrModel(rModel) {};

    /**
     * @brief Constructor reading directly from file, via parameters
     * @param Params The configuration parameters telling where the configuration file can be found
     * @param rModel The model containing the problem to solve
     */
    ReadMaterialsUtility(
        Parameters Params,
        Model& rModel
        );

    /**
     * @brief Constructor reading directly from file, via text
     * @param Params The string telling where the configuration file can be found
     * @param rModel The model containing the problem to solve
     */
    ReadMaterialsUtility(const std::string& rParametersName, Model& rModel);

    /// Destructor.
    virtual ~ReadMaterialsUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This reads the properties from parameters
     * @param MaterialData The configuration parameters defining the properties
     */
    void ReadMaterials(Parameters MaterialData);

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
    std::string Info() const  {
        return "ReadMaterialsUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const  {
        rOStream << "ReadMaterialsUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const  {
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

    /**
     * @brief This method assigns the material parameters to a property from configuration parameters
     * @param MaterialData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the materials are to be assigned
     */
    virtual void AssignMaterialToProperty(
        const Parameters MaterialData,
        Properties& rProperty
        );
    
    /**
     * @brief This method assigns the constitutive law to a property from configuration parameters
     * @param MaterialData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the materials are to be assigned
     */
    virtual void AssignConstitutiveLawToProperty(
        const Parameters MaterialData,
        Properties& rProperty
        );
    
    /**
     * @brief This method assigns the variables to a property from configuration parameters
     * @param MaterialData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the materials are to be assigned
     */
    virtual void AssignVariablesToProperty(
        const Parameters MaterialData,
        Properties& rProperty
        );
    
    /**
     * @brief This method assigns the tables to a property from configuration parameters
     * @param MaterialData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the materials are to be assigned
     */
    virtual void AssignTablesToProperty(
        const Parameters MaterialData,
        Properties& rProperty
        );
    
    /**
     * @brief This method creates an auxiliar Parameters when reading properties in order to avoid error, so these non-registered properties can be processed later
     * @param VariablesParameters The original variable parameters
     * @param PropertyId The current property Id (for a warning)
     * @return The variables filtered if required
     */
    virtual Parameters FilterVariables(
        const Parameters VariablesParameters,
        const IndexType PropertyId = 0
        );
        
    /**
     * @brief Trims out a component name, separating by '."
     * @details Trims out a component name, removing unnecessary module information.
     * For backward compatibility.
     * Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
     * Ex: KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3D -> LinearElastic3D
     * @param rLine Component name in materials json file
     */
    void TrimComponentName(std::string& rLine);
    
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

    Model& mrModel; /// The model containing the model parts

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method creates a list of subproperties and it assigns to the father property from configuration parameters
     * @param rModelPart The currently computed model part
     * @param SubPropertiesData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the subproperties are to be created
     */
    void CreateSubProperties(
        ModelPart& rModelPart,
        const Parameters SubPropertiesData,
        Properties& rProperty
        );

    /**
     * @brief This method assigns the properties to the model parts
     * @param Data The parameters containing all the configurations of the materials
     */
    void AssignPropertyBlock(Parameters Data);

    /**
     * @brief This method gets the properties of the different model parts
     * @param Materials The parameters containing the properties of the materials
     */
    void GetPropertyBlock(Parameters Materials);

    /**
     * @brief Checks if the materials are assigned uniquely to the modelparts
     * @param Materials The parameters containing the properties of the materials
     */
    void CheckUniqueMaterialAssignment(Parameters Materials);

    /**
     * @brief Check if a model part is repeated inside materials
     * @param ModelPartNames vector with the model part names
     */
    void CheckModelPartIsNotRepeated(std::vector<std::string>);


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class ReadMaterialsUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_READ;MATERIALS_H_INCLUDED defined
