//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 15:41:36 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/eulerian_conv_diff.h"

#include "custom_elements/conv_diff_2d.h"
#include "custom_elements/conv_diff_3d.h"
#include "custom_elements/eulerian_conv_diff.h"
#include "custom_elements/eulerian_diff.h"

#include "custom_elements/laplacian_element.h"
#include "custom_conditions/thermal_face2D.h"
#include "custom_conditions/thermal_face3D.h"
#include "custom_conditions/flux_condition.h"

#include "includes/variables.h"
#include "includes/condition.h"


namespace Kratos
{

///@name Kratos Globals
///@{

// Variables definition
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double,  MELT_TEMPERATURE_1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double,  MELT_TEMPERATURE_2)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double,  BFECC_ERROR)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double,  BFECC_ERROR_1)

KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MEAN_SIZE)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, PROJECTED_SCALAR1)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, DELTA_SCALAR1)//
KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, MEAN_VEL_OVER_ELEM_SIZE)

KRATOS_DEFINE_APPLICATION_VARIABLE( CONVECTION_DIFFUSION_APPLICATION, double, THETA)

KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( CONVECTION_DIFFUSION_APPLICATION, CONVECTION_VELOCITY)

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
class KratosConvectionDiffusionApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosConvectionDiffusionApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConvectionDiffusionApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConvectionDiffusionApplication();

    /// Destructor.
    virtual ~KratosConvectionDiffusionApplication() {}


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
        return "KratosConvectionDiffusionApplication";
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
        KRATOS_WATCH("in KratosConvectionDiffusionApplication");
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

    const EulerianConvectionDiffusionElement<2,3>  mEulerianConvDiff2D;
    const EulerianConvectionDiffusionElement<2,4>  mEulerianConvDiff2D4N;
    const EulerianConvectionDiffusionElement<3,4>  mEulerianConvDiff3D;
    const EulerianConvectionDiffusionElement<3,8>  mEulerianConvDiff3D8N;
    const EulerianDiffusionElement<2,3>  mEulerianDiffusion2D;
    const EulerianDiffusionElement<3,4>  mEulerianDiffusion3D;

    const ConvDiff2D  mConvDiff2D;
    const ConvDiff3D  mConvDiff3D;
    const LaplacianElement mLaplacian2D3N;
    const LaplacianElement mLaplacian3D4N;
    const LaplacianElement mLaplacian3D8N;
    const LaplacianElement mLaplacian3D27N;
    const ThermalFace2D  mThermalFace2D;
    const ThermalFace3D  mThermalFace3D;
    const FluxCondition<2>  mFluxCondition2D2N;
    const FluxCondition<3>  mFluxCondition3D3N;
    const FluxCondition<4>  mFluxCondition3D4N;
   
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
    KratosConvectionDiffusionApplication& operator=(KratosConvectionDiffusionApplication const& rOther);

    /// Copy constructor.
    KratosConvectionDiffusionApplication(KratosConvectionDiffusionApplication const& rOther);


    ///@}

}; // Class KratosConvectionDiffusionApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED  defined 


