//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application

#if !defined(KRATOS_IGA_APPLICATION_H_INCLUDED)
#define  KRATOS_IGA_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

//elements
#include "custom_elements/truss_element.h"
#include "custom_elements/truss_embedded_edge_element.h"
#include "custom_elements/iga_membrane_element.h"
#include "custom_elements/shell_3p_element.h"
#include "custom_elements/shell_5p_hierarchic_element.h"
#include "custom_elements/shell_5p_element.h"
#include "custom_elements/laplacian_IGA_element.h"

//conditions
#include "custom_conditions/output_condition.h"
#include "custom_conditions/load_condition.h"
#include "custom_conditions/load_moment_director_5p_condition.h"
#include "custom_conditions/coupling_penalty_condition.h"
#include "custom_conditions/coupling_lagrange_condition.h"
#include "custom_conditions/coupling_nitsche_condition.h"
#include "custom_conditions/support_penalty_condition.h"
#include "custom_conditions/support_lagrange_condition.h"
#include "custom_conditions/support_nitsche_condition.h"
#include "custom_conditions/support_laplacian_condition.h"
#include "custom_conditions/sbm_laplacian_condition.h"

//modelers
#include "custom_modelers/iga_modeler.h"
#include "custom_modelers/iga_thb_modeler.h"
#include "custom_modelers/refinement_modeler.h"
#include "custom_modelers/nurbs_geometry_modeler.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(IGA_APPLICATION) KratosIgaApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosIgaApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosIgaApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosIgaApplication();

    /// Destructor.
    ~KratosIgaApplication() override {}

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information
    std::string Info() const override
    {
        return "KratosIgaApplication";
    }

    /// Print Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in my application");
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

private:

    ///@name Member Variables
    ///@{

    const TrussElement mTrussElement;
    const TrussEmbeddedEdgeElement mTrussEmbeddedEdgeElement;
    const IgaMembraneElement mIgaMembraneElement;
    const Shell3pElement mShell3pElement;
    const Shell5pHierarchicElement mShell5pHierarchicElement;
    const Shell5pElement mShell5pElement;
    const LaplacianIGAElement mLaplacianIGAElement;

    //Conditions
    const OutputCondition mOutputCondition;
    const LoadCondition mLoadCondition;
    const LoadMomentDirector5pCondition mLoadMomentDirector5pCondition;
    const CouplingPenaltyCondition mCouplingPenaltyCondition;
    const CouplingLagrangeCondition mCouplingLagrangeCondition;
    const CouplingNitscheCondition mCouplingNitscheCondition;
    const SupportPenaltyCondition mSupportPenaltyCondition;
    const SupportLagrangeCondition mSupportLagrangeCondition;
    const SupportNitscheCondition mSupportNitscheCondition;
    const SupportLaplacianCondition mSupportLaplacianCondition;
    const SBMLaplacianCondition mSBMLaplacianCondition;

    // Modelers
    const IgaModeler mIgaModeler;
    const IgaThbModeler mIgaThbModeler;
    const RefinementModeler mRefinementModeler;
    const NurbsGeometryModeler mNurbsGeometryModeler;

    ///@}
    ///@name Private methods
    ///@{

    /// Assignment operator.
    KratosIgaApplication& operator=(KratosIgaApplication const& rOther);

    /// Copy constructor.
    KratosIgaApplication(KratosIgaApplication const& rOther);

    ///@}

}; // class KratosIgaApplication

///@}

} // namespace Kratos

#endif // !defined(KRATOS_IGA_APPLICATION_H_INCLUDED)
