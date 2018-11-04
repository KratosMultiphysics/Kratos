// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
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
#include "convection_diffusion_application_variables.h"

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

/** 
 * @class KratosConvectionDiffusionApplication
 * @ingroup KratosConvectionDiffusionApplication
 * @brief The Convection Diffusion Application contains a series of elements and conditions and the corresponding strategies and solvers within Kratos Multiphysics necesaries in order to simulate a convection-diffusion problem
 * @details The application includes tests to check the proper functioning of the application. Features:  
- A set of *Neumann* conditions:
     * Flux conditions
     * Thermal conditions
- Elements:
    * Laplacian element (both 2D/3D)
    * Eulerian convection-diffusion (both 2D/3D)
    * Convection-diffusion (both 2D/3D)
    * Convection-diffusion with change of phase (2D)
- Strategies:
    * Non-linear/linear convection-diffusion strategy
    * Eulerian convection-diffusion strategy
    * Semi-Eulerian convection-diffusion strategy
- Utilities and others:
    * BFECC convection utility
    * BFECC elemental limiter convection utility
    * Convection particle
    * Face-heat utilities
    * Move particle utility
    * Pure convection tools
    * Pure convection (Crank-Nicolson) tools
 * @author Riccardo Rossi
 * @author Pablo Becker
 * @author Jordi Cotela
 * @see KratosApplication
*/
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) KratosConvectionDiffusionApplication : public KratosApplication
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
        return "KratosConvectionDiffusionApplication";
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


