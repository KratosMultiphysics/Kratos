//
//   Project Name:        Kratos
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:54:44 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "structural_application.h"

namespace Kratos
{
    KRATOS_DEFINE_VARIABLE ( Vector, COORDINATES )
    KRATOS_DEFINE_VARIABLE ( Vector, STRESSES )
    KRATOS_DEFINE_VARIABLE ( Vector, FLUID_FLOWS )
    ///@name Kratos Globals
    ///@{

    // Variables definition

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

    class KratosConstitutiveLawsApplication : public KratosApplication
    {

        public:
            ///@name Type Definitions
            ///@{


            /// Pointer definition of KratosExternalSolversApplication
            KRATOS_CLASS_POINTER_DEFINITION ( KratosConstitutiveLawsApplication );

            ///@}
            ///@name Life Cycle
            ///@{

            /// Default constructor.
            KratosConstitutiveLawsApplication() {}

            /// Destructor.
            virtual ~KratosConstitutiveLawsApplication() {}


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
                return "KratosConstitutiveLawsApplication";
            }

            /// Print information about this object.
            virtual void PrintInfo ( std::ostream& rOStream ) const
            {
                rOStream << Info();
                PrintData ( rOStream );
            }

            ///// Print object's data.
            virtual void PrintData ( std::ostream& rOStream ) const
            {
                KRATOS_WATCH ( "in KratosConstitutiveLawsApplication" );
                KRATOS_WATCH ( KratosComponents<VariableData>::GetComponents().size() );
                rOStream << "Variables:" << std::endl;
                KratosComponents<VariableData>().PrintData ( rOStream );
                rOStream << std::endl;
                rOStream << "Elements:" << std::endl;
                KratosComponents<Element>().PrintData ( rOStream );
                rOStream << std::endl;
                rOStream << "Conditions:" << std::endl;
                KratosComponents<Condition>().PrintData ( rOStream );
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
            KratosConstitutiveLawsApplication& operator= ( KratosConstitutiveLawsApplication const& rOther );

            /// Copy constructor.
            KratosConstitutiveLawsApplication ( KratosConstitutiveLawsApplication const& rOther );


            ///@}

    }; // Class KratosConstitutiveLawsApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


}  // namespace Kratos.

#endif // KRATOS_EXTERNAL_CONSTITUTIVE_LAWS_APPLICATION_H_INCLUDED  defined 


