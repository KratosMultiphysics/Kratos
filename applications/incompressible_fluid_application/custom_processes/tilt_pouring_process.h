//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_TILT_POURING_PROCESS_H_INCLUDED )
#define  KRATOS_TILT_POURING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"


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

/// This Process takes a mesh and extract out the individual bodies by analysing the connectivity.
/** .
*/
class TiltPouringProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TiltPouringProcess
    KRATOS_CLASS_POINTER_DEFINITION(TiltPouringProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
	/// @param: InitialAngle in degree
	/// @param: RotationAxis: 1 for x, 2 for y, 3 for z axis
    TiltPouringProcess(ModelPart& rModelPart, double InitialAngle, double EndTime, Point<3> const& RotationPoint, int RotationAxis, array_1d<double,3> const& Gravity) 
		: mrModelPart(rModelPart) , mInitialAngle(InitialAngle), mEndTime(EndTime), mRotationPoint(RotationPoint), mRotationAxis(RotationAxis), mGravity(Gravity)
    {
    }

    /// Copy constructor.
    TiltPouringProcess(TiltPouringProcess const& rOther)
		: mrModelPart(rOther.mrModelPart) , mInitialAngle(rOther.mInitialAngle), mEndTime(rOther.mEndTime), mRotationPoint(rOther.mRotationPoint), mRotationAxis(rOther.mRotationAxis), mGravity(rOther.mGravity)
    {
    }

    /// Destructor.
    virtual ~TiltPouringProcess() {}


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void ExecuteInitializeSolutionStep()
    {
		double current_time = mrModelPart.GetProcessInfo()[TIME];
		double alpha = (mEndTime - current_time) * mInitialAngle * 3.1415 / (mEndTime * 180.0);

		Matrix rotation_matrix(3,3);
		CalculateRotationMatrix(rotation_matrix, mRotationAxis, alpha);

		array_1d<double,3> new_gravity = prod(trans(rotation_matrix), mGravity);
		KRATOS_WATCH(mGravity);
		KRATOS_WATCH(new_gravity);
		ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
			i_node->GetSolutionStepValue(BODY_FORCE) = new_gravity;

			array_1d<double,3> offset = i_node->Coordinates();
			offset -= mRotationPoint;

			array_1d<double,3> new_coordinates = prod(rotation_matrix, offset);
			new_coordinates += mRotationPoint;

			array_1d<double,3> displacement = new_coordinates - i_node->Coordinates();
			i_node->GetSolutionStepValue(DISPLACEMENT) = displacement;

        }

	}
 

	void CalculateRotationMatrix(Matrix& RotationMatrix, int RotationAxis, double Alpha)
	{
		if(RotationAxis == 1) // x
		{
			RotationMatrix(0,0) = 1.00;   RotationMatrix(0,1) = 0.00;         RotationMatrix(0,2) = 0.00;    
			RotationMatrix(1,0) = 0.00;   RotationMatrix(1,1) = cos(Alpha);  RotationMatrix(1,2) = -sin(Alpha); 
			RotationMatrix(2,0) = 0.00;   RotationMatrix(2,1) = sin(Alpha);  RotationMatrix(2,2) = cos(Alpha);  
		}
		else if(RotationAxis == 2) // y
		{
			RotationMatrix(0,0) = cos(Alpha);   RotationMatrix(0,1) = 0.00;  RotationMatrix(0,2) = sin(Alpha);    
			RotationMatrix(1,0) = 0.00;			RotationMatrix(1,1) = 1.00;  RotationMatrix(1,2) = 0.00; 
			RotationMatrix(2,0) = -sin(Alpha);	RotationMatrix(2,1) = 0.00;  RotationMatrix(2,2) = cos(Alpha);  
		}
		else if(RotationAxis == 3) // z
		{
			RotationMatrix(0,0) = cos(Alpha);   RotationMatrix(0,1) = -sin(Alpha);  RotationMatrix(0,2) = 0.00;    
			RotationMatrix(1,0) = sin(Alpha);	RotationMatrix(1,1) = cos(Alpha);   RotationMatrix(1,2) = 0.00; 
			RotationMatrix(2,0) = 0.00;			RotationMatrix(2,1) = 0.00;			RotationMatrix(2,2) = 1.00;  
		}
		else // error
		{
			KRATOS_ERROR(std::invalid_argument, "Invalid rotation axis. The valid axis are 1 for x, 2 for y, 3 for z axis and given axis is: ", RotationAxis);
		}
	}

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TiltPouringProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& mrModelPart;
	double mInitialAngle;
	double mEndTime;
	Point<3> mRotationPoint;
	int mRotationAxis;
	array_1d<double,3> mGravity;

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
    TiltPouringProcess& operator=(TiltPouringProcess const& rOther)
    {
        return *this;
    }


    ///@}

}; // Class TiltPouringProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TiltPouringProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TiltPouringProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TILT_POURING_PROCESS_H_INCLUDED  defined 


