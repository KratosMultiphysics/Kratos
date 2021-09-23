//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Eduard Gómez
//

#if !defined(KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED )
#define  KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "containers/model.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
  ///@{

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

  /// 
  /** This process applies diferent boundary conditions accoring to the mach regime.
   * 
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ApplyCompressibleInlet: public Process
  {
  public:
        ///@name Type Definitions
        ///@{
        //
        typedef ModelPart::NodeType NodeType;
        
        /**
         * @brief This class validates and manages a variable to fix
         * and the value to fix it to.
         */
        class BoundaryConditionUtility
        {
        public:
            BoundaryConditionUtility(const std::string & variable_name,
                                     const double Value,
                                     const double Start,
                                     const double End);
            /**
             * @brief Checks that the curent time is within the interval, and decides accordingly what function 
             * to use during Enforce. This avoids checking the time at every node.
             */
            void ActivateIfInsideInterval(const double time);

            /**
             * @brief Enforces (or not, depending on ActivateIfInsideInterval) the boundary condition.
             */
            void Enforce(NodeType & rNode) const;
        private:
            const Variable<double> * mpVariable;
            const double mValue; // Value to enforce
            const double mStart; // Start of the active interval
            const double mEnd;   // End of the active interval
            
            /// This function will be called when the time is within the specified interval
            static void EnforceActive(const BoundaryConditionUtility & rUtility, NodeType & rNode);

            /// This function will be called when the  time is outside the specified interval
            static void EnforcePassive(const BoundaryConditionUtility & rUtility, NodeType & rNode);

            /// This object points to the proper Enforce(Active|Passive) function
            std::function<decltype(EnforceActive)> mEnforceInternal = EnforcePassive;
        };

        /// Pointer definition of ApplyCompressibleInlet
        KRATOS_CLASS_POINTER_DEFINITION(ApplyCompressibleInlet);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor
        ApplyCompressibleInlet(Model& rModel, Parameters& rParameters);

        /// Destructor.
        ~ApplyCompressibleInlet() override;


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void ExecuteInitializeSolutionStep() override;

        const Parameters GetDefaultParameters() const override;

        bool ValidateParameters(Parameters Params) const;

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
        std::string Info() const override;

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override;

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override;


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

        ModelPart * mpModelPart;
        std::vector<BoundaryConditionUtility> mSubsonicBCs = {};
        std::vector<BoundaryConditionUtility> mSupersonicBCs = {};

        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{
        
        /**
         * @brief Reads the JSON parameters for a variable and appends the resulting
         * BoundaryConditionUtility to the provided list
         */
        void ReadBoundaryCondition(std::vector<BoundaryConditionUtility> & rBCList, Parameters Parameters);

        /**
         * @brief Performs the conversion [0, 1, 2] -> ['X', 'Y', 'Z'] at compile-time
         */
        static constexpr char IndexToAxis(const unsigned int i);

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
        ApplyCompressibleInlet& operator=(ApplyCompressibleInlet const& rOther);

        /// Copy constructor.
        ApplyCompressibleInlet(ApplyCompressibleInlet const& rOther);


        ///@}

        }; // Class ApplyCompressibleInlet

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyCompressibleInlet& rThis);

    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED  defined
