//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_MOR_APPLICATION_H_INCLUDED )
#define  KRATOS_MOR_APPLICATION_H_INCLUDED


// System includes


// External includes
#include "includes/variables.h"

// Project includes
#include "includes/kratos_application.h"

#include "custom_elements/acoustic_element.h"

#include "custom_conditions/acoustic_load_condition.h"
#include "custom_conditions/acoustic_robin_condition.h"
#include "custom_conditions/acoustic_structure_coupling_condition.h"
#include "custom_conditions/displacement_output_condition.h"
#include "custom_conditions/pressure_output_condition.h"

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MOR_APPLICATION) KratosMORApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMORApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMORApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMORApplication();

    /// Destructor.
    ~KratosMORApplication() override {}

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
        return "KratosMORApplication";
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
          KRATOS_WATCH("in my application");
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

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    // Elements
    const AcousticElement mAcousticElement2D2N;
    const AcousticElement mAcousticElement2D4N;
    const AcousticElement mAcousticElement3D4N;
    const AcousticElement mAcousticElement3D8N;

    // Conditions
    const AcousticLoadCondition mAcousticLoadConcition2D2N;
    const AcousticLoadCondition mAcousticLoadConcition3D3N;
    const AcousticLoadCondition mAcousticLoadConcition3D4N;
    const AcousticRobinCondition mAcousticRobinConcition2D2N;
    const AcousticRobinCondition mAcousticRobinConcition3D3N;
    const AcousticRobinCondition mAcousticRobinConcition3D4N;
    const AcousticStructureCouplingCondition<2, false> mAcousticStructureCouplingCondition2D2N;
    const AcousticStructureCouplingCondition<3, false> mAcousticStructureCouplingCondition3D4N;
    const AcousticStructureCouplingCondition<3, false> mAcousticStructureCouplingCondition3D3N;
    const AcousticStructureCouplingCondition<2, true> mAcousticStructureMappingCondition2D2N;
    const AcousticStructureCouplingCondition<3, true> mAcousticStructureMappingCondition3D3N;
    const AcousticStructureCouplingCondition<3, true> mAcousticStructureMappingCondition3D4N;

    const DisplacementOutputCondition mDisplacementOutputCondition3D1N;
    const PressureOutputCondition mPressureOutputCondition3D1N;

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
    KratosMORApplication& operator=(KratosMORApplication const& rOther);

    /// Copy constructor.
    KratosMORApplication(KratosMORApplication const& rOther);


    ///@}

}; // Class KratosMORApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_MOR_APPLICATION_H_INCLUDED  defined
