//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//


#if !defined(KRATOS_ROTATION_OPERATION_H_INCLUDED )
#define  KRATOS_ROTATION_OPERATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/table.h"
#include "includes/mesh.h"



namespace Kratos
{

///@name Kratos Classes
///@{

/// The Rotation Operation.
/** This Operation is a derived class from the process.h
 *  
*/

class RotationOperation : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(RotationOperation);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RotationOperation(ModelPart& model_part, DenseVector<int> group_ids,DenseVector<int> table_ids,unsigned int echo_level=0):
	Process(),mr_model_part(model_part),mgroup_ids(group_ids),mtable_ids(table_ids)
     {
	mecho_level=echo_level;
     }

    /// Destructor.
    ~RotationOperation() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the Process algorithms.
    /*
    virtual void Execute() 
    {

    }


    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }
*/

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
				KRATOS_TRY
		if ((mr_model_part.NumberOfTables())==0)
			            KRATOS_THROW_ERROR(std::logic_error, "Tables of the modelpart are empty", "");
		if (mgroup_ids.size()==0)
			KRATOS_THROW_ERROR(std::logic_error, "No groups to rotate", "");
		if (mtable_ids.size()<6)
			KRATOS_THROW_ERROR(std::logic_error, "Table's Vector too small!. Must be at least of size 6 for the 3 rotations + 3 reference(center) coordiantes", "");

		ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double time = CurrentProcessInfo[TIME];
		//double delta_t = CurrentProcessInfo[DELTA_TIME];	
			
		//tables are used in the following way
		//table(table_ids[0]) is the x displacement, table(table_ids[1]) is the y_displacement and table(table_ids[2]) is the z_displacement
		//table0(time,displacement(time)
		Table<double,double>& RotationTableA = mr_model_part.GetTable(mtable_ids[0]);
		//KRATOS_WATCH(mtable_ids[0]);
		Table<double,double>& RotationTableB = mr_model_part.GetTable(mtable_ids[1]);
		//KRATOS_WATCH(mtable_ids[1]);
		Table<double,double>& RotationTableC = mr_model_part.GetTable(mtable_ids[2]);

		Table<double,double>& ReferenceTableX = mr_model_part.GetTable(mtable_ids[3]);
		Table<double,double>& ReferenceTableY = mr_model_part.GetTable(mtable_ids[4]);
		Table<double,double>& ReferenceTableZ = mr_model_part.GetTable(mtable_ids[5]);
		//KRATOS_WATCH(mtable_ids[2]);
		array_1d<double,3> rotation = ZeroVector(3); //the 3 angles of the rotation
		rotation(0)=RotationTableA(time);
		rotation(1)=RotationTableB(time);
		rotation(2)=RotationTableC(time);
		array_1d<double,3> reference_point = ZeroVector(3); //the reference center of coordinates for the rotation.
		reference_point(0)=ReferenceTableX(time);
		reference_point(1)=ReferenceTableY(time);
		reference_point(2)=ReferenceTableZ(time);

		BoundedMatrix<double, 3, 3 > rotation_matrix;
		const double c1=cos(rotation(0));
		const double c2=cos(rotation(1));
		const double c3=cos(rotation(2));
		const double s1=sin(rotation(0));
		const double s2=sin(rotation(1));
		const double s3=sin(rotation(2));

		rotation_matrix(0,0)=c2*c3;     rotation_matrix(0,1)=c1*s3+s1*s2*c3;      rotation_matrix(0,2)=s1*s3-c1*s2*c3;
		
		rotation_matrix(1,0)=-c2*s3;    rotation_matrix(1,1)=c1*c3-s1*s2*s3;      rotation_matrix(1,2)=s1*c3+c1*s2*s3;

		rotation_matrix(2,0)=s2;	rotation_matrix(2,1)=-s1*c2;  		   rotation_matrix(2,2)=c1*c2;




		
		for (unsigned int mesh_index=0;mesh_index<mgroup_ids.size();mesh_index++) //we loop around the desired groups
		{
			const int mesh_id=mgroup_ids[mesh_index];
			ModelPart::MeshType& current_mesh = mr_model_part.GetMesh(mesh_id);
			ModelPart::NodesContainerType::iterator inodebegin = current_mesh.NodesBegin();
			//ModelPart::NodesContainerType::iterator inodeend = mgroup_container(mesh_id).NodesEnd();
			#pragma omp parallel for
			for(int ii=0; ii< static_cast<int>(current_mesh.Nodes().size()); ii++)
			{
				ModelPart::NodesContainerType::iterator pnode = inodebegin+ii;
				//pnode->Coordinates()=pnode->X0()+translation(0);
				const array_1d<double,3> relative_position = pnode->GetInitialPosition().Coordinates() - reference_point;
				const array_1d<double,3> new_position = prod(rotation_matrix,relative_position) + reference_point ;
								
				if (pnode->SolutionStepsDataHas(DISPLACEMENT_X)) //
					pnode->FastGetSolutionStepValue(DISPLACEMENT) = new_position - pnode->GetInitialPosition().Coordinates();
				pnode->Coordinates()  = new_position;
			}
		}
		// = (i->FastGetSolutionStepValue(PRESS_PROJ_NO_RO));
        KRATOS_CATCH("")
    }
/*
    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
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
    virtual std::string Info() const
    {
        return "Process";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Process";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
*/

    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ModelPart& mr_model_part;
   // ModelPart::MeshType& mgroup_container;
    DenseVector<int> mgroup_ids;
    DenseVector<int> mtable_ids;
    unsigned int mecho_level;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TranslationOperation& operator=(Process const& rOther);

    /// Copy constructor.
    //Process(Process const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/*
/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Process& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Process& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
* */
///@}


}  // namespace Kratos.

#endif // KRATOS_ROTATION_OPERATION_H_INCLUDED  defined 


