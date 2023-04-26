//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SHALLOW_WATER_APPLICATION_H_INCLUDED
#define KRATOS_SHALLOW_WATER_APPLICATION_H_INCLUDED

///@defgroup ShallowWaterApplication Kratos Shallow Water Application
///@brief Basic set of tools to solve the shallow water equations.
/// The Shallow Water Application implements a basic set of tools to
/// solve shallow water problems. This applications contains a basic FEM
/// implementation of common techniques using eulerian schemes.


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Shallow water includes
#include "custom_elements/wave_element.h"
#include "custom_elements/primitive_element.h"
#include "custom_elements/crank_nicolson_wave_element.h"
#include "custom_elements/boussinesq_element.h"
#include "custom_elements/conservative_element.h"
#include "custom_elements/conservative_element_rv.h"
#include "custom_elements/conservative_element_fc.h"
#include "custom_conditions/wave_condition.h"
#include "custom_conditions/primitive_condition.h"
#include "custom_conditions/boussinesq_condition.h"
#include "custom_conditions/conservative_condition.h"
#include "custom_modelers/mesh_moving_modeler.h"


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

    /// Short class definition.
    /** Detail class definition.
    */
    class KRATOS_API(SHALLOW_WATER_APPLICATION) KratosShallowWaterApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{


        /// Pointer definition of KratosShallowWaterApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosShallowWaterApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosShallowWaterApplication();

        /// Destructor.
        virtual ~KratosShallowWaterApplication(){}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Register() override;

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
        virtual std::string Info() const override
        {
            return "KratosShallowWaterApplication";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
            KRATOS_WATCH("in Shallow Water Application");
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

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{

        // Elements
        const WaveElement<3> mWaveElement2D3N;
        const WaveElement<6> mWaveElement2D6N;
        const WaveElement<4> mWaveElement2D4N;
        const WaveElement<8> mWaveElement2D8N;
        const WaveElement<9> mWaveElement2D9N;
        const WaveElement<3> mPrimitiveElement2D3N;
        const WaveElement<4> mPrimitiveElement2D4N;
        const CrankNicolsonWaveElement<3> mCrankNicolsonWaveElement2D3N;
        const BoussinesqElement<3> mBoussinesqElement2D3N;
        const BoussinesqElement<4> mBoussinesqElement2D4N;
        const ConservativeElement<3> mConservativeElementGJ2D3N;
        const ConservativeElementRV<3> mConservativeElementRV2D3N;
        const ConservativeElementFC<3> mConservativeElementFC2D3N;
        // Conditions
        const WaveCondition<2> mWaveCondition2D2N;
        const WaveCondition<3> mWaveCondition2D3N;
        const PrimitiveCondition<2> mPrimitiveCondition2D2N;
        const BoussinesqCondition<2> mBoussinesqCondition2D2N;
        const ConservativeCondition<2> mConservativeCondition2D2N;

        // Modelers
        const MeshMovingModeler mMeshMovingModeler;

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
        KratosShallowWaterApplication& operator=(KratosShallowWaterApplication const& rOther);

        /// Copy constructor.
        KratosShallowWaterApplication(KratosShallowWaterApplication const& rOther);


        ///@}

    }; // Class KratosShallowWaterApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


}  // namespace Kratos.

#endif // KRATOS_SHALLOW_WATER_APPLICATION_H_INCLUDED  defined
