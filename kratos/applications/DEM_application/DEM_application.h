//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/DEM_FEM_Particle.h"

const long double pi = 3.141592653589793238462643383279;

namespace Kratos
{
    /* Define In Global variables.h
        KRATOS_DEFINE_VARIABLE(double,  DELTA_TIME);
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_ROTATE_SPRING_FAILURE_TYPE)
        typedef vector<array_1d<double,3> > VectorArray3Double;
        KRATOS_DEFINE_VARIABLE( VectorArray3Double, PARTICLE_ROTATE_SPRING_MOMENT )
     */

        KRATOS_DEFINE_VARIABLE (int, NEIGH_INITIALIZED)
        KRATOS_DEFINE_VARIABLE (int, TOTAL_CONTACTS )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)

        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
        KRATOS_DEFINE_VARIABLE(WeakPointerVector<Element >,     NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

  
        KRATOS_DEFINE_VARIABLE(int, VIRTUAL_MASS_OPTION)
        KRATOS_DEFINE_VARIABLE(double, NODAL_MASS_COEFF)

        KRATOS_DEFINE_VARIABLE(double, HISTORICAL_MIN_K)
        KRATOS_DEFINE_VARIABLE(int, CRITICAL_TIME_OPTION)

        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
        KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

        KRATOS_DEFINE_VARIABLE(Vector,     INITIAL_AXES_TRACKING)
        KRATOS_DEFINE_VARIABLE(int,     plot_OPTIONS)

 
        KRATOS_DEFINE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
        KRATOS_DEFINE_VARIABLE(Vector, IF_BOUNDARY_FACE)



class KratosDEMApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosDEMSpheresApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosDEMApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosDEMApplication();

    /// Destructor.
    virtual ~KratosDEMApplication() {}


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
        return "KratosDEMApplication";
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
        //KRATOS_WATCH("in KratosDEMApplication");
        //KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
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
    const SphericParticle mSphericParticle2D;
    const SphericParticle mSphericParticle3D;

    const DEM_FEM_Particle mDEM_FEM_Particle2D;
    const DEM_FEM_Particle mDEM_FEM_Particle3D;


    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{
// 		const Elem2D   mElem2D;
// 		const Elem3D   mElem3D;


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
    KratosDEMApplication& operator=(KratosDEMApplication const& rOther);

    /// Copy constructor.
    KratosDEMApplication(KratosDEMApplication const& rOther);


    ///@}

}; // Class KratosDEMApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_DEM_APPLICATION_H_INCLUDED  defined 


