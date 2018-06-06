//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_READ_MATERIALS_H_INCLUDED )
#define  KRATOS_READ_MATERIALS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/io.h"

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

/// Process to read constitutive law and material properties from a json file
/** This process reads constitutive law and material properties from a json file
 * and assign them to elements and conditions.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ReadMaterialsUtility
{
  public:

    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ReadMaterialProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReadMaterialsUtility);

    //typedef std::size_t SizeType;

    //typedef ModelPart::NodeType::DofsContainerType DofsContainerType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    ReadMaterialsUtility(Parameters& rParameters, Model& rModel);

    ReadMaterialsUtility(const std::string& rParameters, Model& rModel);

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

    Model& mrModel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    void AssignPropertyBlock(Parameters& rData);
    void GetPropertyBlock(Parameters& rMaterials);

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
