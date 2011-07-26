//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MixedElement_APPLICATION_H_INCLUDED )
#define  KRATOS_MixedElement_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "custom_elements/sigma_u_element.h"

#include "includes/variables.h"


namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition

    KRATOS_DEFINE_VARIABLE(double, SX)
    KRATOS_DEFINE_VARIABLE(double, SY)
    KRATOS_DEFINE_VARIABLE(double, SZ)
    KRATOS_DEFINE_VARIABLE(double, SXY)
    KRATOS_DEFINE_VARIABLE(double, SXZ)
    KRATOS_DEFINE_VARIABLE(double, SYZ)
    //	KRATOS_DEFINE_VARIABLE(double, AUX_MESH_VAR )
    //	KRATOS_DEFINE_VARIABLE(double, IS_INTERFACE)
    //	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


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
    class KratosMixedElementApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{


        /// Pointer definition of KratosMixedElementApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosMixedElementApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosMixedElementApplication();

        /// Destructor.

        virtual ~KratosMixedElementApplication()
        {
        }


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
            return "KratosMixedElementApplication";
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
            KRATOS_WATCH("in my application");
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData > ().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element > ().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition > ().PrintData(rOStream);
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
        const SigmaUElement   mSigmaUElement2D;
        // 		const Elem3D   mElem3D;


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
        KratosMixedElementApplication & operator=(KratosMixedElementApplication const& rOther);

        /// Copy constructor.
        KratosMixedElementApplication(KratosMixedElementApplication const& rOther);


        ///@}

    }; // Class KratosMixedElementApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos.

#endif // KRATOS_MixedElement_APPLICATION_H_INCLUDED  defined 


