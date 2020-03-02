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

#if !defined(KRATOS_SHALLOW_WATER_APPLICATION_H_INCLUDED )
#define  KRATOS_SHALLOW_WATER_APPLICATION_H_INCLUDED

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
#include "custom_elements/shallow_element.h"
#include "custom_elements/rv_swe.h"
#include "custom_elements/cv_swe.h"
#include "custom_elements/swe.h"
#include "custom_elements/conserved_element.h"
#include "custom_conditions/nothing_condition.hpp"


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

        // Elements
        const ShallowElement mShallowElement2D3N;
        const RV_SWE<3, Eulerian> mRVSWE2D3N;
        const RV_SWE<4, Eulerian> mRVSWE2D4N;
        const RV_SWE<3, PFEM2> mPFEM2RVSWE2D3N;
        const RV_SWE<4, PFEM2> mPFEM2RVSWE2D4N;
        const CV_SWE<3, Eulerian> mCVSWE2D3N;
        const CV_SWE<4, Eulerian> mCVSWE2D4N;
        const CV_SWE<3, PFEM2> mPFEM2CVSWE2D3N;
        const CV_SWE<4, PFEM2> mPFEM2CVSWE2D4N;
        const SWE<3, Eulerian> mSWE2D3N;
        const SWE<4, Eulerian> mSWE2D4N;
        const SWE<3, PFEM2> mLagrangianSWE2D3N;
        const SWE<4, PFEM2> mLagrangianSWE2D4N;
        const ConservedElement<3> mConservedElement2D3N;
        const ConservedElement<4> mConservedElement2D4N;
        // Condition
        const NothingCondition<2> mNothingCondition2D2N;


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
