//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  G.Casas$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
namespace Kratos
{
class KRATOS_API(PLASMA_DYNAMICS_APPLICATION) KratosPlasmaDynamicsApplication : public KratosApplication
{
public:

    /// Pointer definition of KratosPlasmaDynamicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPlasmaDynamicsApplication);

    /// Default constructor.
    KratosPlasmaDynamicsApplication();

    /// Destructor.
    virtual ~KratosPlasmaDynamicsApplication() {}


    virtual void Register() override;

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "KratosPlasmaDynamicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

protected:

private:    

    /// Assignment operator.
    KratosPlasmaDynamicsApplication& operator=(KratosPlasmaDynamicsApplication const& rOther);

    /// Copy constructor.
    KratosPlasmaDynamicsApplication(KratosPlasmaDynamicsApplication const& rOther);

}; // Class KratosPlasmaDynamicsApplication

}  // namespace Kratos.

#endif // KRATOS_PLASMA_DYNAMICS_APPLICATION_H_INCLUDED  defined


