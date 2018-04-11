
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-04-01 10:25:52 $
//   Revision:            $Revision: 1.7 $
//
//


#if !defined(KRATOS_KRATOS_ULF_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_ULF_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/condition.h"
#include "custom_elements/updated_lagrangian_fluid.h"
#include "custom_elements/updated_lagrangian_fluid3D.h"
#include "custom_elements/updated_lagrangian_fluid_inc.h"
#include "custom_elements/updated_lagrangian_fluid3D_inc.h"
#include "custom_elements/ulf_frac2d.h"
#include "custom_elements/ulf_frac3d.h"
#include "custom_elements/ulf_axisym.h"
//#include "custom_elements/ulf_frac2d_swimming.h"
//#include "custom_elements/ulf_frac3d_swimming.h"
#include "custom_conditions/Point_Neumann3D.h"
#include "custom_conditions/Point_Neumann2D.h"
#include "custom_conditions/Point_Neumann_Axisym.h"
#include "custom_elements/surface_tension.h"
#include "includes/ublas_interface.h"

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KratosULFApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosULFApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosULFApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosULFApplication();

    /// Destructor.
    virtual ~KratosULFApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



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
    virtual std::string Info() const
    {
        return "KratosULFApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in KratosULFApplication");
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
    const UpdatedLagrangianFluid mUpdatedLagrangianFluid2D;
    const UpdatedLagrangianFluid3D mUpdatedLagrangianFluid3D;
    const UpdatedLagrangianFluidInc mUpdatedLagrangianFluid2Dinc;
    const UpdatedLagrangianFluid3Dinc mUpdatedLagrangianFluid3Dinc;
    //
    const UlfFrac2D mUlfFrac2D;
    const UlfFrac3D mUlfFrac3D;
    const UlfAxisym mUlfAxisym;    
    const PointNeumann3D  mPointNeumann3D;
    const PointNeumann2D  mPointNeumann2D;
    const PointNeumannAxisym  mPointNeumannAxisym;
    
       /// 2D instance of the SurfaceTension element
    const SurfaceTension<2> mSurfaceTension2D;
    /// 3D instance of the SurfaceTension element
    const SurfaceTension<3> mSurfaceTension3D;

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
    KratosULFApplication& operator=(KratosULFApplication const& rOther);

    /// Copy constructor.
    KratosULFApplication(KratosULFApplication const& rOther);


    ///@}

}; // Class KratosULFApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_ULF_APPLICATION_H_INCLUDED  defined 
