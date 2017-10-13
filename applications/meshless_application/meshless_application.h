//
//   Project Name:        Kratos       
//   Last modified by:    $Author: c.karacaova
//   Date:                $Date: 2012-08-11  $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MESHLESS_APPLICATION_H_INCLUDED )
#define  KRATOS_MESHLESS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
//#include "custom_elements/SPHparticle.h"
#include "custom_elements/ULF_LMEparticle.h"
//#include "custom_elements/TLF_LMEparticle.h"
#include "custom_elements/TLF_MLSparticle.h"
#include "custom_elements/ULF_MLSparticle.h"

#include "includes/ublas_interface.h"
#include "meshless_application_variables.h"
#include "custom_conditions/lagrange_multiplier_condition.h"
#include "custom_conditions/lagrange_multiplier_condition_x.h"
#include "custom_conditions/lagrange_multiplier_condition_y.h"
#include "custom_conditions/pointforce2D.h"
//#include "custom_constitutive/isotropic_rankine_damage_2d.h"
//#include "custom_constitutive/brittle_material_2d.h"
//#include "custom_constitutive/plane_strain.h"
//#include "fluency_criteria/fluency_criteria.h"
//#include "fluency_criteria/standard_morh_coulomb_yield_function.h"
//#include "soft_hard_behavior/cohesion.h"
//#include "soft_hard_behavior/cohesion_softening.h"
//#include "soft_hard_behavior/softening_hardening_criteria.h"

namespace Kratos
{

///@name Kratos Globals
///@{

// Variables definition
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
class KratosMeshlessApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosMeshlessApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMeshlessApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMeshlessApplication();

    /// Destructor.
    virtual ~KratosMeshlessApplication(){}


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
        return "KratosMeshlessApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{
    // 		const Elem2D   mElem2D;
    // 		const Elem3D   mElem3D;


//    const SPHparticle<KernelC2,KernelC2,KernelC2> mSPHparticle; // For the SPH particle element
//    const SPHparticle<KernelQuintic,KernelQuintic,KernelQuintic> mSPHparticle; // For the SPH particle element


/*
    const SPHparticle<KernelPoly6,KernelSpiky,KernelMullerViscous> mSPHparticlePoly; // Element 1

    const SPHparticle<KernelSpiky,KernelSpiky,KernelSpiky> mSPHparticlePolyPresSpiky; // Element 2

    const SPHparticle<KernelQuintic,KernelQuintic,KernelQuintic> mSPHparticleQuintic; // Element 3

    const SPHparticle<KernelC2,KernelC2,KernelC2> mSPHparticleC2; // Element 4

    const SPHparticle<KernelQuintic,KernelQuadratic,KernelQuintic> mSPHparticlePolyPresQuad; // Element 5

    const SPHparticle<KernelQuadratic,KernelQuadratic,KernelQuadratic> mSPHparticleGaus; // Element 6

    //const SPHparticle<LinearMLSKernel,LinearMLSKernel,LinearMLSKernel> mSPHparticleGaus;

    //const MLSparticle<LinearMLSKernel> mMLSparticleMls; // Element 7

    const TLF_LMEparticle<LinearLMEKernel> mLinearElement2D3N;  //just for generate mesh
    //const MLSparticle<LinearMLSKernel> mLinearElement2D4N;  //just for generate mesh

    const TLF_LMEparticle<LinearLMEKernel> mTLFLMEparticleLme; // Element 8


    const ULF_LMEparticle<LinearLMEKernel> mULFLMEparticleLme; // Element 9
 */


    const ULF_MLSparticle mLinearElement2D3N;  //just for generate mesh
    const ULF_MLSparticle mULFMLSparticleMls; // Element 9

    const LagrangeMultiplierCondition2D mLagrangeMultiplierCondition2D;
    const LagrangeMultiplierCondition2DX mLagrangeMultiplierCondition2DX;
    const LagrangeMultiplierCondition2DY mLagrangeMultiplierCondition2DY;
    const PointForce2D  mPointForce2D;
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
    KratosMeshlessApplication& operator=(KratosMeshlessApplication const& rOther);

    /// Copy constructor.
    KratosMeshlessApplication(KratosMeshlessApplication const& rOther);


    ///@}

}; // Class KratosMeshlessApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_MESHLESS_APPLICATION_H_INCLUDED  defined
