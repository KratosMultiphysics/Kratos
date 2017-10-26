//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


#if !defined(KRATOS_SHELL_UTILITIES_HPP_INCLUDED )
#define  KRATOS_SHELL_UTILITIES_HPP_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
   * This class is intended to group functionalities needed for all shells in order to avoid code duplication
  */
  class ShellUtilities
    {
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShellUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ShellUtilities);

    typedef Element::GeometryType GeometryType;

    typedef Properties PropertiesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShellUtilities() = delete;

    /// Destructor.
    ~ShellUtilities() = default;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void CheckVariables();

    static void CheckDofs(GeometryType& rGeom);

    static void CheckProperties(const Element* pTheElement, const ProcessInfo& rCurrentProcessInfo, 
                                const bool IsThickShell = false);


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
    virtual std::string Info() const;

    /// Print information about pTheElement object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    
    static void CheckSpecificProperties(const Element* pTheElement, const PropertiesType & rProps, const bool IsThickShell);


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
    ShellUtilities& operator=(ShellUtilities const& rOther);

    /// Copy constructor.
    ShellUtilities(ShellUtilities const& rOther);


    ///@}

    }; // Class ShellUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    ShellUtilities& rThis);

  /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
// 				    const ShellUtilities& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHELL_UTILITIES_HPP_INCLUDED  defined
