//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  G.Casas$
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"  //TODO: must be removed eventually
#include "includes/cfd_variables.h"  //TODO: must be removed eventually
#include "includes/legacy_structural_app_vars.h"  //TODO: must be removed eventually
#include "custom_elements/monolithic_dem_coupled.h"
#include "custom_elements/monolithic_dem_coupled_weak.h"
#include "custom_elements/calculate_laplacian_simplex_element.h"
#include "custom_elements/calculate_mat_deriv_simplex_element.h"
#include "custom_elements/calculate_component_gradient_simplex_element.h"
#include "custom_elements/calculate_gradient_Pouliot_2012.h"
#include "custom_elements/calculate_gradient_Pouliot_2012_edge.h"
#include "custom_elements/calculate_velocity_laplacian_component.h"
#include "custom_elements/calculate_velocity_laplacian.h"
#include "custom_elements/shell_rigid.h"
#include "custom_conditions/monolithic_dem_coupled_wall_condition.h"
#include "custom_conditions/calculate_laplacian_simplex_condition.h"

#include "custom_elements/spheric_swimming_particle.h"
#include "../DEMApplication/custom_elements/spheric_particle.h"
#include "../DEMApplication/custom_elements/nanoparticle.h"
#include "../DEMApplication/custom_elements/analytic_spheric_particle.h"

namespace Kratos
{

    #define SWIMMING_COPY_SECOND_TO_FIRST_3(a, b)            a[0]  = b[0]; a[1]  = b[1]; a[2]  = b[2];
    #define SWIMMING_ADD_SECOND_TO_FIRST(a, b)               a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
    #define SWIMMING_SET_COMPONENTS_TO_ZERO_3(a)             a[0]  = 0.0;  a[1]  = 0.0;  a[2]  = 0.0;
    #define SWIMMING_SET_COMPONENTS_TO_ZERO_3x3(a)           a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0; a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 0.0;
    #define SWIMMING_MULTIPLY_BY_SCALAR_3(a, b)              a[0] = b * a[0]; a[1] = b * a[1]; a[2] = b * a[2];
    #define SWIMMING_MODULUS_3(a)                            std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    #define SWIMMING_INNER_PRODUCT_3(a, b)                            (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
    #define SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(a, b, c)    c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
    #define SWIMMING_POW_2(a)                                (a * a)
    #define SWIMMING_POW_3(a)                                (a * a * a)
    #define SWIMMING_POW_4(a)                                (a * a * a * a)
    #define SWIMMING_POW_5(a)                                (a * a * a * a * a)
    #define SWIMMING_POW_6(a)                                (a * a * a * a * a * a)
    #define SWIMMING_POW_7(a)                                (a * a * a * a * a * a * a)

class KRATOS_API(SWIMMING_DEM_APPLICATION) KratosSwimmingDEMApplication : public KratosApplication
{
public:

    /// Pointer definition of KratosSwimmingDEMApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosSwimmingDEMApplication);

    /// Default constructor.
    KratosSwimmingDEMApplication();

    /// Destructor.
    virtual ~KratosSwimmingDEMApplication() {}


    void Register() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosSwimmingDEMApplication";
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
    ///@name Static Member Variables
    ///@{
    /// 2D instance of the MonolithicDEMCoupled element
    const MonolithicDEMCoupled<2> mMonolithicDEMCoupled2D;
    /// 3D instance of the MonolithicDEMCoupled element
    const MonolithicDEMCoupled<3> mMonolithicDEMCoupled3D;

    /// 2D instance of the MonolithicDEMCoupledWeak element
    const MonolithicDEMCoupledWeak<2> mMonolithicDEMCoupledWeak2D;
    /// 3D instance of the MonolithicDEMCoupledWeak element
    const MonolithicDEMCoupledWeak<3> mMonolithicDEMCoupledWeak3D;

    const ComputeLaplacianSimplex<2> mComputeLaplacianSimplex2D;
    const ComputeLaplacianSimplex<3> mComputeLaplacianSimplex3D;

    const ComputeMaterialDerivativeSimplex<2> mComputeMaterialDerivativeSimplex2D;
    const ComputeMaterialDerivativeSimplex<3> mComputeMaterialDerivativeSimplex3D;

    const ComputeComponentGradientSimplex<2> mComputeComponentGradientSimplex2D;
    const ComputeComponentGradientSimplex<3> mComputeComponentGradientSimplex3D;

    const ComputeGradientPouliot2012Edge<2> mComputeGradientPouliot20122DEdge;
    const ComputeGradientPouliot2012Edge<3> mComputeGradientPouliot20123DEdge;

    const ComputeGradientPouliot2012<2> mComputeGradientPouliot20122D;
    const ComputeGradientPouliot2012<3> mComputeGradientPouliot20123D;

    const ComputeVelocityLaplacianComponentSimplex<2> mComputeVelocityLaplacianComponentSimplex2D;
    const ComputeVelocityLaplacianComponentSimplex<3> mComputeVelocityLaplacianComponentSimplex3D;

    const ComputeVelocityLaplacianSimplex<2> mComputeVelocityLaplacianSimplex2D;
    const ComputeVelocityLaplacianSimplex<3> mComputeVelocityLaplacianSimplex3D;

    const  MonolithicDEMCoupledWallCondition<2,2> mMonolithicDEMCoupledWallCondition2D;
    const  MonolithicDEMCoupledWallCondition<3,3> mMonolithicDEMCoupledWallCondition3D;

    const  ComputeLaplacianSimplexCondition<2,2> mComputeLaplacianSimplexCondition2D;
    const  ComputeLaplacianSimplexCondition<3,3> mComputeLaplacianSimplexCondition3D;

    const ShellRigid mRigidShellElement;

    /// swimming derivation of spheric basic DEM element (SphericParticle)
    const SphericSwimmingParticle<SphericParticle> mSphericSwimmingParticle3D;
    const SphericSwimmingParticle<NanoParticle> mSwimmingNanoParticle3D;
    const SphericSwimmingParticle<AnalyticSphericParticle> mSwimmingAnalyticParticle3D;

    /// Assignment operator.
    KratosSwimmingDEMApplication& operator=(KratosSwimmingDEMApplication const& rOther);

    /// Copy constructor.
    KratosSwimmingDEMApplication(KratosSwimmingDEMApplication const& rOther);

}; // Class KratosSwimmingDEMApplication

}  // namespace Kratos.

#endif // KRATOS_SWIMMING_DEM_APPLICATION_H_INCLUDED  defined


