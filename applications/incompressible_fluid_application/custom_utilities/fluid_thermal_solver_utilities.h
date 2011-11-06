//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_FLUID_THERMAL_SOLVER_UTILITIES_H_INCLUDED )
#define  KRATOS_FLUID_THERMAL_SOLVER_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "modeler/duplicate_mesh_modeler.h"


namespace Kratos {
    ///@addtogroup ApplicationNameApplication
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

    /// Short class definition.

    /** Detail class definition.
     */
    class FluidThermalSolverUtilities {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of FluidThermalSolverUtilities
        KRATOS_CLASS_POINTER_DEFINITION(FluidThermalSolverUtilities);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        FluidThermalSolverUtilities(ModelPart& rFluidModelPart, ModelPart& rThermalModelPart) : mrFluidModelPart(rFluidModelPart), mrThermalModelPart(rThermalModelPart) {
            // Check all variables
            // Generate Copy of Elements
            DuplicateMeshModeler copy_modeler(mrFluidModelPart);

            copy_modeler.GenerateMesh(mrThermalModelPart, KratosComponents<Element>::Get("ConvDiff3D"), KratosComponents<Condition>::Get("ThermalFace3D"));


            mSkinNodes.clear();
            // Compute list of skin nodes
            for(ModelPart::ConditionsContainerType::iterator i_condition = mrThermalModelPart.ConditionsBegin(); i_condition != mrThermalModelPart.ConditionsEnd(); i_condition++)
            {
                Condition::GeometryType& geometry = i_condition->GetGeometry();
                for(std::size_t i = 0 ; i < geometry.size() ; i++)
                    mSkinNodes.push_back(geometry(i));
            }

            mSkinNodes.Unique();


            Check();


        }


        /// Destructor.

        virtual ~FluidThermalSolverUtilities() {
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        int Check()
        {
            //veryfying that the model part has all the variables needed
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(DENSITY) == false)
                KRATOS_ERROR(std::logic_error, "Add DENSITY variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_ERROR(std::logic_error, "Add TEMPERATURE variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(CONDUCTIVITY) == false)
                KRATOS_ERROR(std::logic_error, "Add CONDUCTIVITY variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(HEAT_FLUX) == false)
                KRATOS_ERROR(std::logic_error, "Add HEAT_FLUX variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(FACE_HEAT_FLUX) == false)
                KRATOS_ERROR(std::logic_error, "Add FACE_HEAT_FLUX variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_ERROR(std::logic_error, "Add MESH_VELOCITY variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR(std::logic_error, "Add VELOCITY variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(CONVECTION_COEFFICIENT) == false)
                KRATOS_ERROR(std::logic_error, "Add CONVECTION_COEFFICIENT variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_ERROR(std::logic_error, "Add DISTANCE variable!!!!!! ERROR", "");
            if (mrThermalModelPart.NodesBegin()->SolutionStepsDataHas(VISCOSITY) == false)
                KRATOS_ERROR(std::logic_error, "Add VISCOSITY variable!!!!!! ERROR", "");
            if(mrThermalModelPart.ElementsBegin()->GetGeometry().Dimension() != 3)
                KRATOS_ERROR(std::invalid_argument, "This algorithm only works with 3D geometries", "")

		  return 0;

        }

//        void AddVariablesToThermalModelPart()
//        {
//            mrThermalModelPart.AddNodalSolutionStepVariable(DENSITY);
//            mrThermalModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
//            mrThermalModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
//            mrThermalModelPart.AddNodalSolutionStepVariable(HEAT_FLUX);
//            mrThermalModelPart.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
//            mrThermalModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
//            mrThermalModelPart.AddNodalSolutionStepVariable(VELOCITY);
//            mrThermalModelPart.AddNodalSolutionStepVariable(CONVECTION_COEFFICIENT);
//        }

        /// Copy constructor.

        FluidThermalSolverUtilities(FluidThermalSolverUtilities const& rOther)
        : mrFluidModelPart(rOther.mrFluidModelPart), mrThermalModelPart(rOther.mrThermalModelPart) {
        }

        void ProjectFromThermalToFluid() {
            // Transfer viscosity
            for(size_t i=0; i<mrFluidModelPart.Nodes().size(); i++)
            {
                (mrFluidModelPart.NodesBegin() + i)->FastGetSolutionStepValue(VISCOSITY)
                        = (mrThermalModelPart.NodesBegin() + i)->FastGetSolutionStepValue(VISCOSITY);
            }
        }

        void ProjectFromFluidToThermal() {
            // Transfer velocity and distance
            for(size_t i=0; i<mrThermalModelPart.Nodes().size(); i++)
            {
                (mrThermalModelPart.NodesBegin() + i)->FastGetSolutionStepValue(DISTANCE)
                        = (mrFluidModelPart.NodesBegin() + i)->FastGetSolutionStepValue(DISTANCE);
                (mrThermalModelPart.NodesBegin() + i)->FastGetSolutionStepValue(VELOCITY)
                        = (mrFluidModelPart.NodesBegin() + i)->FastGetSolutionStepValue(VELOCITY);
            }

        }


        void ApplyTables() {

        }

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

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "FluidThermalSolverUtilities";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
            rOStream << "FluidThermalSolverUtilities";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
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

        ModelPart& mrFluidModelPart;
        ModelPart& mrThermalModelPart;
        
        ModelPart::NodesContainerType mSkinNodes;


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

        FluidThermalSolverUtilities & operator=(FluidThermalSolverUtilities const& rOther) {
            return *this;
        }


        ///@}

    }; // Class FluidThermalSolverUtilities

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
            FluidThermalSolverUtilities& rThis) {
        return rIStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const FluidThermalSolverUtilities& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FLUID_THERMAL_SOLVER_UTILITIES_H_INCLUDED  defined


