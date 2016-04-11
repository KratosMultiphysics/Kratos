//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Apr 19, 2012 $
//   Revision:            $Revision: 1.1 $
//
//
//Change log:
//Dec 31, 2012: change application name to phase_field_application

#if !defined(KRATOS_PHASE_FIELD_APPLICATION_H_INCLUDED )
#define  KRATOS_PHASE_FIELD_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/c2c_variables.h"
#include "custom_elements/phase_field_fracture.h"
#include "custom_elements/phase_field_fracture_hybrid.h"
#include "custom_elements/phase_field_kinematic_linear.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{ 

    // Variables definition
    KRATOS_DEFINE_VARIABLE(double, PHASE_FIELD)
    KRATOS_DEFINE_VARIABLE(double, PHASE_FIELD_DUAL_VARIABLE)
    KRATOS_DEFINE_VARIABLE(double, PHASE_FIELD_GRADIENT)
    KRATOS_DEFINE_VARIABLE(double, LENGTH_SCALE)
    KRATOS_DEFINE_VARIABLE(double, REFERENCE_ENERGY_DENSITY)
    KRATOS_DEFINE_VARIABLE(double, ENERGY_FUNCTIONAL_DENSITY)
    KRATOS_DEFINE_VARIABLE(double, ENERGY_FUNCTIONAL)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(INTEGRATION_POINT_GLOBAL_COORDINATES)
    KRATOS_DEFINE_VARIABLE(int, PHASE_FIELD_ORDER)

    // variable imported from structural_application
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESCRIBED_DELTA_DISPLACEMENT)

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
    class KratosPhaseFieldApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{
        

        /// Pointer definition of KratosDiscontinuitiesApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosPhaseFieldApplication);

        ///@}
        ///@name Life Cycle 
        ///@{ 

        /// Default constructor.
        KratosPhaseFieldApplication();

        /// Destructor.
        virtual ~KratosPhaseFieldApplication()
        {}

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
            return "KratosPhaseFieldApplication";
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
            KRATOS_WATCH("in KratosPhaseFieldApplication");
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

        const PhaseFieldFracture mPhaseFieldFracture2D3N;
        const PhaseFieldFracture mPhaseFieldFracture2D4N;
        const PhaseFieldFracture mPhaseFieldFracture2D6N;
        const PhaseFieldFracture mPhaseFieldFracture2D8N;
        const PhaseFieldFracture mPhaseFieldFracture2D9N;
        const PhaseFieldFracture mPhaseFieldFracture3D4N;
        const PhaseFieldFracture mPhaseFieldFracture3D8N;
        const PhaseFieldFracture mPhaseFieldFracture3D10N;
        const PhaseFieldFracture mPhaseFieldFracture3D20N;
        const PhaseFieldFracture mPhaseFieldFracture3D27N;

        const PhaseFieldKinematicLinear mPhaseFieldKinematicLinear3D4N;

        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid2D3N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid2D4N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid2D6N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid2D8N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid2D9N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid3D4N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid3D8N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid3D10N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid3D20N;
        const PhaseFieldFractureHybrid mPhaseFieldFractureHybrid3D27N;

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
        KratosPhaseFieldApplication& operator=(KratosPhaseFieldApplication const& rOther);

        /// Copy constructor.
        KratosPhaseFieldApplication(KratosPhaseFieldApplication const& rOther);


        ///@}    

    }; // Class KratosPhaseFieldApplication 

    ///@} 


    ///@name Type Definitions       
    ///@{ 


    ///@} 
    ///@name Input and output 
    ///@{ 

    ///@} 


}  // namespace Kratos.

#endif // KRATOS_PHASE_FIELD_APPLICATION_H_INCLUDED  defined

