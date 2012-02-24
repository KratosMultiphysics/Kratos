//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: croig $
//   Date:                $Date: 2012-02-10 11:11:35 $
//   Revision:            $Revision: 1.16 $
//
//


#if !defined(KRATOS_KRATOS_MPI_SEARCH_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_MPI_SEARCH_APPLICATION_H_INCLUDED

///@defgroup MPI Search Application



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/condition.h"

namespace Kratos
{
    ///@addtogroup MPISerchApplication
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
     */
    class KratosMPISearchApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{


        /// Pointer definition of KratosMPISolverApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosMPISearchApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosMPISearchApplication();

        /// Destructor.

        virtual ~KratosMPISearchApplication()
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
            return "KratosMPISearchApplication";
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
//             KRATOS_WATCH("in KratosMPISearchApplication");
//             KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
//             rOStream << "Variables:" << std::endl;
//             KratosComponents<VariableData > ().PrintData(rOStream);
//             rOStream << std::endl;
//             rOStream << "Elements:" << std::endl;
//             KratosComponents<Element > ().PrintData(rOStream);
//             rOStream << std::endl;
//             rOStream << "Conditions:" << std::endl;
//             KratosComponents<Condition > ().PrintData(rOStream);
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
//        const Fluid2Dlevelset mFluid2Dlevelset;

        //const ABC2D mABC2D;

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
        KratosMPISearchApplication & operator=(KratosMPISearchApplication const& rOther);

        /// Copy constructor.
        KratosMPISearchApplication(KratosMPISearchApplication const& rOther);


        ///@}

    }; // Class KratosMPISerchApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}

    ///@} KratosMPISearchApplication group

} // namespace Kratos.

#endif // KRATOS_KRATOS_MPI_SEARCH_APPLICATION_H_INCLUDED  defined 


