//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-09-28 12:56:44 $
//   Revision:            $Revision: 1.4 $
//
//

#if !defined(KRATOS_FSI_APPLICATION_H_INCLUDED)
#define KRATOS_FSI_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/fsi_variables.h"

namespace Kratos {

///@name Kratos Globals
///@{

class KratosFSIApplication : public KratosApplication {
   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMeshMovingApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFSIApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFSIApplication() : KratosApplication("FSIApplication") {}

    /// Destructor.
    ~KratosFSIApplication() override {}

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
     std::string Info() const override { return "KratosFSIApplication"; }

    /// Print information about this object.
     void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
     void PrintData(std::ostream& rOStream) const override{
        KRATOS_WATCH("in FSIApplication");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
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

    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

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
    KratosFSIApplication& operator=(KratosFSIApplication const& rOther);

    /// Copy constructor.
    KratosFSIApplication(KratosFSIApplication const& rOther);

    ///@}

};  // Class KratosFSIApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif  // KRATOS_FSI_APPLICATION_H_INCLUDED  defined
