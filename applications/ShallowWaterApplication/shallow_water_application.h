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
/// implementation of common techniques using both explicit pfem2 and
/// eulerian shemes.


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Shallow water includes
#include "custom_elements/swe.h"
#include "custom_elements/wave_element.h"
#include "custom_elements/crank_nicolson_wave_element.h"
#include "custom_elements/conservative_element.h"
#include "custom_elements/shallow_water_2d_3.h"
#include "custom_conditions/wave_condition.h"
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
        const SWE<3, Eulerian> mSWE2D3N;
        const SWE<4, Eulerian> mSWE2D4N;
        const SWE<3, PFEM2> mLagrangianSWE2D3N;
        const SWE<4, PFEM2> mLagrangianSWE2D4N;
        const WaveElement<3> mWaveElement2D3N;
        const WaveElement<6> mWaveElement2D6N;
        const WaveElement<4> mWaveElement2D4N;
        const WaveElement<8> mWaveElement2D8N;
        const WaveElement<9> mWaveElement2D9N;
        const CrankNicolsonWaveElement<3> mCrankNicolsonWaveElement2D3N;
        const ConservativeElement<3> mConservativeElement2D3N;
        const ShallowWater2D3 mShallowWater2D3N;
        // Conditions
        const WaveCondition<2> mWaveCondition2D2N;
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
