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

#if !defined(KRATOS_APPLY_MACH_DEPENDENT_BOUNDARY_CONDITIONS_H_INCLUDED )
#define  KRATOS_APPLY_MACH_DEPENDENT_BOUNDARY_CONDITIONS_H_INCLUDED



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
#include "utilities/interval_utility.h"
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

  /** 
   * @ brief This process applies diferent boundary conditions accoring to the
   * mach regime.
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ApplyMachDependentBoundaryConditions: public Process
  {
  public:
        ///@name Type Definitions
        ///@{
        //
        typedef ModelPart::NodeType NodeType;
        
        /**
         * @brief This class validates and manages a variable to fix
         * and the value to fix it to.
         * 
         * Method ActivateIfInsideTimeInterval must be called at the begining
         * of the time step to see if the BC is active or not. This will cause
         * the following behaviour in the other two public methods:
         * - Enforce:
         *     - If t inside interval  -> Fix dof
         *     - If t outside interval -> Free dof
         * - Release
         *     - Always -> Free dof
         * 
         * This allows checking the time only once at the beginning of the
         * time-step rather than at every node.
         */
        class BoundaryConditionUtility
        {
        public:
            BoundaryConditionUtility(const std::string & variable_name,
                                     const double Value,
                                     const IntervalUtility & rIntervalUtility);
            /**
             * @brief Checks that the current time is within the interval, and
             * decides accordingly what function to bind to mEnforceInternal
             * (either FreeDof or FixDof).
             */
            void ActivateIfInsideTimeInterval(const double time);

            /**
             * @brief Calls mEnforceInternal. This will either fix or free the
             * Dof depending on the result of ActivateIfInsideTimeInterval.
             */
            void Enforce(NodeType & rNode) const;

            /// @brief Frees the Dof.
            void Release(NodeType & rNode) const;

        private:
            const Variable<double> * mpVariable;
            const double mValue; // Value to enforce
            IntervalUtility mInterval;
            
            /* @brief Fixes the Dof managed by this class and sets the value
             * to the one specified by mValue
             */
            static void FixDof(const BoundaryConditionUtility & rUtility, NodeType & rNode);

            /// @brief Frees the Dof managed by this class
            static void FreeDof(const BoundaryConditionUtility & rUtility, NodeType & rNode);

            /** @brief This object points to the proper (Fix|Free)Dof function
             * according to the time interval.
             */
            std::function<decltype(FixDof)> mEnforceInternal = FreeDof;
        };

        /// Pointer definition of ApplyMachDependentBoundaryConditions
        KRATOS_CLASS_POINTER_DEFINITION(ApplyMachDependentBoundaryConditions);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor
        ApplyMachDependentBoundaryConditions(Model& rModel, Parameters Parameters);


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void ExecuteInitializeSolutionStep() override;

        const Parameters GetDefaultParameters() const override;

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
         * @brief Reads the JSON parameters for a variable and appends the 
         * resulting BoundaryConditionUtility to the provided list
         */
        void ReadBoundaryCondition(std::vector<BoundaryConditionUtility> & rBCList, Parameters Parameters);

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
        ApplyMachDependentBoundaryConditions& operator=(ApplyMachDependentBoundaryConditions const& rOther);

        /// Copy constructor.
        ApplyMachDependentBoundaryConditions(ApplyMachDependentBoundaryConditions const& rOther);


        ///@}

        }; // Class ApplyMachDependentBoundaryConditions

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyMachDependentBoundaryConditions& rThis);

    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_MACH_DEPENDENT_BOUNDARY_CONDITIONS_H_INCLUDED  defined
