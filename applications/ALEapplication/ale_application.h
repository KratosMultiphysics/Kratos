//
//   Project Name:        Kratos
//   Last Modified by:    $Author: AMini $
//   Date:                $Date: Oct 2014 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_ALE_APPLICATION_H_INCLUDED )
#define  KRATOS_ALE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/laplacian_meshmoving_element_2d.h"
#include "custom_elements/laplacian_meshmoving_element_3d.h"
#include "custom_elements/laplacian_componentwise_meshmoving_element_2d.h"
#include "custom_elements/laplacian_componentwise_meshmoving_element_3d.h"
#include "custom_elements/structural_meshmoving_element.h"
#include "custom_elements/laplacian_componentwise_meshmoving_element_2d_strainbased.h"

//#include "custom_elements/laplacian_componentwise_meshmoving_element_3d_strainbased.h"

#include "custom_elements/structural_meshmoving_element_2d_nonlinear.h"
#include "custom_elements/structural_meshmoving_element_3d_nonlinear.h"


#include "includes/variables.h"


namespace Kratos
{

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
class KratosALEApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosALEApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosALEApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosALEApplication();

    /// Destructor.
    virtual ~KratosALEApplication() {}


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
        return "KratosALEApplication";
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
    const LaplacianMeshMovingElem2D   mLaplacianMeshMovingElem2D;
    const LaplacianMeshMovingElem3D   mLaplacianMeshMovingElem3D;
    const LaplacianMeshMovingElem2D   mLaplacianComponentwiseMeshMovingElem2D;
    const LaplacianMeshMovingElem3D   mLaplacianComponentwiseMeshMovingElem3D;
    const StructuralMeshMovingElement<2> mStructuralMeshMovingElement2D;
    const StructuralMeshMovingElement<3> mStructuralMeshMovingElement3D;
    const LaplacianComponentwiseMeshMovingElem2DStrainbased mLaplacianComponentwiseMeshMovingElem2DStrainbased;
    //const LaplacianComponentwiseMeshMovingElem3DStrainbased mLaplacianComponentwiseMeshMovingElem3DStrainbased;
    const StructuralMeshMovingElem2DNonlin  mStructuralMeshMovingElem2DNonlin;
    const StructuralMeshMovingElem3DNonlin  mStructuralMeshMovingElem3DNonlin;


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
    KratosALEApplication& operator=(KratosALEApplication const& rOther);

    /// Copy constructor.
    KratosALEApplication(KratosALEApplication const& rOther);


    ///@}

}; // Class KratosALEApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_ALE_APPLICATION_H_INCLUDED  defined 


