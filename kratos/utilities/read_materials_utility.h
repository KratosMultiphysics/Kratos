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
 * @author Marcelo Raschi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) ReadMaterialsUtility
{
  public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    //typedef std::size_t SizeType;

    //typedef ModelPart::NodeType::DofsContainerType DofsContainerType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ReadMaterialProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReadMaterialsUtility);
    
    ///@}
    ///@name Life Cycle
    ///@{

    ReadMaterialsUtility(Parameters Params, Model& rModel);

    ReadMaterialsUtility(const std::string& rParametersName, Model& rModel);

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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
    virtual std::string Info() const  {
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
     * @brief This method creates a property from configuration parameters
     * @param Data The parameters containing all the configurations of the materials
     * @param pNewProperty The pointer to the new property created
     */
    void CreateProperty(
        Parameters Data,
        Properties::Pointer& pNewProperty
        );

    /**
     * @brief This method creates a list of subproperties and it assigns to the father property from configuration parameters
     * @param rModelPart The currently computed model part
     * @param MeshId The mesh id (0 by default, keep it for retrocompatibility)
     * @param Data The parameters containing all the configurations of the materials
     * @param pNewProperty The pointer to the new property created
     */
    void CreateSubProperties(
        ModelPart& rModelPart,
        const IndexType MeshId,
        Parameters Data,
        Properties::Pointer& pNewProperty
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
     * @brief Trims out a component name, separating by '."
     * @details Trims out a component name, removing unnecessary module information.
     * For backward compatibility.
     * Ex: KratosMultiphysics.YOUNG_MODULUS -> YOUNG_MODULUS
     * Ex: KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3D -> LinearElastic3D
     * @param rLine Component name in materials json file
     */
    void TrimComponentName(std::string& rLine);


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
