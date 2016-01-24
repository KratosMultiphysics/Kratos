//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes 
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"

#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"
#include "containers/flags.h"

#include "custom_conditions/point_load_condition.hpp"

#include "custom_conditions/line_load_2D_condition.hpp"
#include "custom_conditions/line_normal_load_2D_condition.hpp"
#include "custom_conditions/line_normal_fluid_flux_2D_condition.hpp"
#include "custom_conditions/surface_load_3D_condition.hpp"
#include "custom_conditions/surface_normal_load_3D_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_condition.hpp"

#include "custom_conditions/line_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include "custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.hpp"
#include "custom_conditions/surface_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_load_3D_diff_order_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.hpp"

#include "custom_conditions/line_normal_fluid_flux_2D_FIC_condition.hpp"
#include "custom_conditions/surface_normal_fluid_flux_3D_FIC_condition.hpp"

#include "custom_elements/small_strain_U_Pw_element.hpp"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_elements/small_strain_U_Pw_FIC_element.hpp"

#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_strain_law.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_stress_law.hpp"

namespace Kratos
{

//Define Variables
KRATOS_DEFINE_VARIABLE( double, BETA_NEWMARK )
KRATOS_DEFINE_VARIABLE( double, GAMMA_NEWMARK )
KRATOS_DEFINE_VARIABLE( double, THETA_NEWMARK )

KRATOS_DEFINE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
KRATOS_DEFINE_VARIABLE( std::string, CONSTITUTIVE_LAW_NAME )

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_NORMAL_STRESS )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_TANGENTIAL_STRESS )

KRATOS_DEFINE_VARIABLE( double, DERIVATIVE_WATER_PRESSURE )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_FLUID_PRESSURE )
KRATOS_DEFINE_VARIABLE( double, NORMAL_FLUID_FLUX )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_NORMAL_FLUID_FLUX )

KRATOS_DEFINE_VARIABLE( double, DENSITY_SOLID )
KRATOS_DEFINE_VARIABLE( double, BULK_MODULUS_SOLID )
KRATOS_DEFINE_VARIABLE( double, BULK_MODULUS_FLUID )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_XX )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_YY )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_ZZ )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_XY )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_YZ )
KRATOS_DEFINE_VARIABLE( double, PERMEABILITY_ZX )

KRATOS_DEFINE_VARIABLE( Vector, CAUCHY_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE( Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
KRATOS_DEFINE_VARIABLE( double, VON_MISES_STRESS )

KRATOS_DEFINE_VARIABLE( Vector, FLUID_FLUX )


class KratosPoromechanicsApplication : public KratosApplication
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(KratosPoromechanicsApplication);

    // Default constructor
    KratosPoromechanicsApplication();

    // Destructor
    virtual ~KratosPoromechanicsApplication(){}
    

    virtual void Register();

    // Turn back information as a string
    virtual std::string Info() const
    {
        return "KratosPoromechanicsApplication";
    }

    // Print information about this object
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    // Print object's data
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

private:

// Member Variables

const SmallStrainUPwElement mSmallStrainUPwElement2D3N;
const SmallStrainUPwElement mSmallStrainUPwElement2D4N;
const SmallStrainUPwElement mSmallStrainUPwElement3D4N;
const SmallStrainUPwElement mSmallStrainUPwElement3D8N;

const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D6N;
const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D8N;
const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement2D9N;
const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D10N;
const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D20N;
const SmallStrainUPwDiffOrderElement mSmallStrainUPwDiffOrderElement3D27N;

const SmallStrainUPwFICElement mSmallStrainUPwFICElement2D3N;
const SmallStrainUPwFICElement mSmallStrainUPwFICElement2D4N;
const SmallStrainUPwFICElement mSmallStrainUPwFICElement3D4N;
const SmallStrainUPwFICElement mSmallStrainUPwFICElement3D8N;

const PointLoadCondition mPointLoadCondition2D;
const PointLoadCondition mPointLoadCondition3D;

const LineLoad2DCondition mLineLoadCondition2D2N;
const LineNormalLoad2DCondition mLineNormalLoadCondition2D2N;
const LineNormalFluidFlux2DCondition mLineNormalFluidFluxCondition2D2N;
const SurfaceLoad3DCondition mSurfaceLoadCondition3D3N;
const SurfaceLoad3DCondition mSurfaceLoadCondition3D4N;
const SurfaceNormalLoad3DCondition mSurfaceNormalLoadCondition3D3N;
const SurfaceNormalLoad3DCondition mSurfaceNormalLoadCondition3D4N;
const SurfaceNormalFluidFlux3DCondition mSurfaceNormalFluidFluxCondition3D3N;
const SurfaceNormalFluidFlux3DCondition mSurfaceNormalFluidFluxCondition3D4N;

const LineLoad2DDiffOrderCondition mLineLoadDiffOrderCondition2D3N;
const LineNormalLoad2DDiffOrderCondition mLineNormalLoadDiffOrderCondition2D3N;
const LineNormalFluidFlux2DDiffOrderCondition mLineNormalFluidFluxDiffOrderCondition2D3N;
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D6N; //TODO: the computation of the area in the 3D_6-nodded element needs to be revised
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D8N;
const SurfaceLoad3DDiffOrderCondition mSurfaceLoadDiffOrderCondition3D9N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D6N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D8N;
const SurfaceNormalLoad3DDiffOrderCondition mSurfaceNormalLoadDiffOrderCondition3D9N;
const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D6N;
const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D8N;
const SurfaceNormalFluidFlux3DDiffOrderCondition mSurfaceNormalFluidFluxDiffOrderCondition3D9N;

const LineNormalFluidFlux2DFICCondition mLineNormalFluidFluxFICCondition2D2N;
const SurfaceNormalFluidFlux3DFICCondition mSurfaceNormalFluidFluxFICCondition3D3N;
const SurfaceNormalFluidFlux3DFICCondition mSurfaceNormalFluidFluxFICCondition3D4N;

const LinearElastic3DLaw mLinearElastic3DLaw;
const LinearElastic2DPlaneStrainLaw mLinearElastic2DPlaneStrainLaw;
const LinearElastic2DPlaneStressLaw mLinearElastic2DPlaneStressLaw;

// Assignment operator.
KratosPoromechanicsApplication& operator=(KratosPoromechanicsApplication const& rOther);

// Copy constructor.
KratosPoromechanicsApplication(KratosPoromechanicsApplication const& rOther);

}; // Class KratosPoromechanicsApplication 
}  // namespace Kratos.

#endif // KRATOS_POROMECHANICS_APPLICATION_H_INCLUDED  defined 


