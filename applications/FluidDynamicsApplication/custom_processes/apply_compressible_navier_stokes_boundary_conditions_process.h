//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Eduard Gómez
//

#if !defined(KRATOS_APPLY_COMPRESSIBLE_NAVIER_STOKES_BOUNDARY_CONDITIONS_PROCESS_INCLUDED )
#define  KRATOS_APPLY_COMPRESSIBLE_NAVIER_STOKES_BOUNDARY_CONDITIONS_PROCESS_INCLUDED


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
   * @brief This process applies different boundary conditions according to the
   * mach regime.
   * Parameters:
   * @param model_part_name The model part onto which to enforce the boundary
   *    conditions.
   * @param flow_direction_variable A variable to use to know direction of the
   *    flow. Its magnitude will me ignored. must be historical.
   * @param refresh_normals_every_time_step Whether to recompute the normal
   *    direction every time-step.Leave at false unless your mesh or geometry
   *    is expected to change.
   * @param subsonic_boundary_conditions A list of variables to fix when the
   *    mach projected onto the normal of the boundary is below 1.
   * @param supersonic_boundary_conditions  A list of variables to fix when
   *    the mach projected onto the normal of the boundary is above 1.
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ApplyCompressibleNavierStokesBoundaryConditionsProcess: public Process
  {
  public:
        ///@name Type Definitions
        ///@{
        //
        typedef ModelPart::NodeType NodeType;
        typedef Variable<array_1d<double, 3>> VectorVariable;
        
        /**
         * @brief This class validates and manages a variable to fix
         * and the value to fix it to.
         * 
         * Method ActivateIfInsideTimeInterval must be called at the beginning
         * of the time step to see if the BC is active or not. This will cause
         * the following behaviour in the public method:
         * - Enforce:
         *     - If t inside interval  -> Fix dof
         *     - If t outside interval -> Do nothing
         * 
         * This allows for checking the time only once at the beginning of the
         * time-step rather than at every node.
         */
        class BoundaryConditionUtility
        {
        public:
            BoundaryConditionUtility(
                const std::string& rVariableName,
                const double Value,
                const IntervalUtility& rIntervalUtility);
            /**
             * @brief Checks that the current time is within the interval, and
             * decides accordingly what function to bind to mEnforceInternal
             * (either FreeDof or DoNothing).
             */
            void ActivateIfInsideTimeInterval(const double Time);

            /**
             * @brief Calls mEnforceInternal. This will either fix the Dof
             * or do nothing depending on the result of 
             * ActivateIfInsideTimeInterval.
             */
            void Enforce(NodeType& rNode) const;

            /// @brief Returns a reference to the stored variable
            const Variable<double> & GetVariable() const;

        private:
            const Variable<double> * mpVariable;    // Variable to fix
            const double mValue;                    // Value to enforce
            IntervalUtility mInterval;              // Interval to enforce in
            
            /** 
             * @brief Fixes the Dof managed by this class and sets the value
             * to the one specified by mValue
             */
            static void FixDof(const BoundaryConditionUtility& rUtility, NodeType& rNode);

            /** 
             * @brief This object points to the proper (FixDof|nullptr)
             * function according to the time interval.
             */
            decltype(FixDof) * mEnforceInternal = nullptr;
        };

        /// Pointer definition of ApplyCompressibleNavierStokesBoundaryConditionsProcess
        KRATOS_CLASS_POINTER_DEFINITION(ApplyCompressibleNavierStokesBoundaryConditionsProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor
        ApplyCompressibleNavierStokesBoundaryConditionsProcess(Model& rModel, Parameters Parameters);


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{
        
        void ExecuteInitialize() override;

        void ExecuteInitializeSolutionStep() override;

        void ExecuteFinalizeSolutionStep() override;

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

        ModelPart * mpModelPart;
        VectorVariable const * mpFlowDirectionVariable;
        bool mRefreshNormalsEveryTimeStep;
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
        ApplyCompressibleNavierStokesBoundaryConditionsProcess& operator=(ApplyCompressibleNavierStokesBoundaryConditionsProcess const& rOther);

        /// Copy constructor.
        ApplyCompressibleNavierStokesBoundaryConditionsProcess(ApplyCompressibleNavierStokesBoundaryConditionsProcess const& rOther);


        ///@}

        }; // Class ApplyCompressibleNavierStokesBoundaryConditionsProcess

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyCompressibleNavierStokesBoundaryConditionsProcess& rThis);

    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_COMPRESSIBLE_NAVIER_STOKES_BOUNDARY_CONDITIONS_PROCESS_INCLUDED  defined
