// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_KRATOS_MESHING_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_MESHING_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/element.h"

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

    /**
     * @brief This enum defines the type of MMG library used
     */
    enum class MMGLibrary
    {
        MMG2D = 0,
        MMG3D = 1,
        MMGS  = 2
    };

    /**
     * @brief This enum defines the type of PMMG library used
     */
    enum class PMMGLibrary
    {
        PMMG3D = 1
    };


    /**
     * @brief This enums allows to differentiate the working framework
     */
    enum class FrameworkEulerLagrange {EULERIAN = 0, LAGRANGIAN = 1, ALE = 2};

    /**
     * @brief This enums allows to differentiate the discretization options
     */
    enum class DiscretizationOption {STANDARD = 0, LAGRANGIAN = 1, ISOSURFACE = 2};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MESHING_APPLICATION) KratosMeshingApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosMeshingApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMeshingApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMeshingApplication();

    /// Destructor.
    ~KratosMeshingApplication() override = default;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;



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
        return "KratosMeshingApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in KratosMeshingApplication");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
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

    const Element mTestElement2D;
    const Element mTestElement3D;

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
    KratosMeshingApplication& operator=(KratosMeshingApplication const& rOther);

    /// Copy constructor.
    KratosMeshingApplication(KratosMeshingApplication const& rOther);


    ///@}

}; // Class KratosMeshingApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_MESHING_APPLICATION_H_INCLUDED  defined
