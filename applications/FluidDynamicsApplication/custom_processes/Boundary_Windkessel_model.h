//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//                   Jordi Cotela
//

#if !defined(KRATOS_WINDKESSEL_H_INCLUDED )
#define  KRATOS_WINDKESSEL_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"


// Application includes
#include "fluid_dynamics_application_variables.h"

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

/// An implementation of the Windkessel model for boundary condition of incompressible flows.
/** Detail class definition.
 */

class WindkesselModel : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WindkesselModel
    KRATOS_CLASS_POINTER_DEFINITION(WindkesselModel);

    ///@}
    ///@name Life Cycle
    ///@{

    WindkesselModel(

            ModelPart& ThisModelPart
            //typename TLinearSolver::Pointer pLinearSolver,
            //unsigned int DomainSize,
            //double NonLinearTol,
            //unsigned int MaxIter,
            //bool ReformDofSet,
            //unsigned int TimeOrder
            )
            : mr_model_part(ThisModelPart)
            //mmax_it(MaxIter), mtime_order(TimeOrder),madapt_for_fractional_step(false)
            //, mdomain_size(DomainSize), mtol(NonLinearTol), mmax_it(MaxIter), mtime_order(TimeOrder),madapt_for_fractional_step(false)
        {
         }

    /// Destructor.

    ~WindkesselModel() override
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Solve an iteration of the turbulent viscosity
    void Execute() override
    {
        KRATOS_TRY
        double Flow_total=0.0;
        double Resistance = 5e9;
        if ((Resistance <=0.0)){
             KRATOS_THROW_ERROR(std::logic_error, "Resistance must be higher than zero:", Resistance);
            }
        ModelPart::ConditionsContainerType rConditions = mr_model_part.Conditions();
        //1_Compute Flow in the outlew area
        for(ModelPart::ConditionsContainerType::iterator i = rConditions.begin(); i!=rConditions.end(); i++){
            //KRATOS_WATCH(i->GetProperties().Id())
            if(i->GetProperties().Id() == 1001){
                //i->NodesContainerType rNodes = mr_model_part.Nodes();
                double Area = i->GetGeometry().Area();
                //double Velocity_total_X=0.0;
                //double Velocity_total_Y=0.0;
                //double Velocity_total_Z=0.0;
                double total_velocity=0.0;
                double velocity =0.0;
                double Flow_total_node=0.0;
                for (unsigned int j = 0 ; j < i->GetGeometry().PointsNumber(); j++){
                    const double Velocity_node_X= i->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X);
                    const double Velocity_node_Y= i->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y);
                    const double Velocity_node_Z= i->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z);
                    velocity=sqrt(pow(Velocity_node_X,2.0) + pow(Velocity_node_Y,2.0) + pow(Velocity_node_Z,2.0));
                    total_velocity += velocity;
//                  Velocity_total_X +=Velocity_node_X;
//                  Velocity_total_Y +=Velocity_node_Y;
//                  Velocity_total_Z +=Velocity_node_Z;
                }
                //Flow_total=((Velocity_total_X+Velocity_total_Y+Velocity_total_Z)/3);

                Flow_total_node=total_velocity/3;
                Flow_total+= Flow_total_node*Area;
                //KRATOS_WATCH(Area)
            }
        }
        //2_Compute the Pressure according to the Windkessel Model
        double Update_Pressure = Flow_total * Resistance;
        //KRATOS_WATCH(Flow_total)
        //3_Assign Pressure over the outlet area
        ModelPart::NodesContainerType rNodes = mr_model_part.Nodes();
        for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
        {
            const int Flag = i->FastGetSolutionStepValue(FLAG_VARIABLE);
            if(Flag == 1001)
            {
                i->FastGetSolutionStepValue(PRESSURE) = Update_Pressure;
                //KRATOS_WATCH(Update_Pressure)
            }
        }
        KRATOS_WATCH("Windkessel_Pressure")
        KRATOS_CATCH("");
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

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "WindkesselModel";
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WindkesselModel";
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
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

    ModelPart& mr_model_part;
    //ModelPart mspalart_model_part;






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

    WindkesselModel & operator=(WindkesselModel const& rOther)
    {
        return *this;
    }

    /// Copy constructor.

    WindkesselModel(WindkesselModel const& rOther)
        : mr_model_part(rOther.mr_model_part)
    {
    }


    ///@}

}; // Class WindkesselModel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
                                  WindkesselModel& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const WindkesselModel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_WINDKESSEL_H_INCLUDED defined


