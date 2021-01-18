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

#if !defined(KRATOS_TRANSLATION_OPERATION_H_INCLUDED )
#define  KRATOS_TRANSLATION_OPERATION_H_INCLUDED



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

/// The Transalation Operation.
/** This Operation is a derived class from the process.h
 *  
*/

class TranslationOperation : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(TranslationOperation);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TranslationOperation(ModelPart& model_part, DenseVector<int> group_ids,DenseVector<int> table_ids,unsigned int echo_level=0):
	Process(),mr_model_part(model_part),mgroup_ids(group_ids),mtable_ids(table_ids)
     {
	mecho_level=echo_level;
     }

    /// Destructor.
    ~TranslationOperation() override {}


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
		KRATOS_THROW_ERROR(std::logic_error, "No groups to translate", "");
	if (mtable_ids.size()<3)
		KRATOS_THROW_ERROR(std::logic_error, "Table's Vector too small!. Must be at least of size 3 for the 3 displacements", "");
	ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	double time = CurrentProcessInfo[TIME];

	//tables are used in the following way
	//table(table_ids[0]) is the x displacement, table(table_ids[1]) is the y_displacement and table(table_ids[2]) is the z_displacement
	//table0(time,displacement(time)
    
	Table<double,double>& TranslationTableX = mr_model_part.GetTable(mtable_ids[0]);
	Table<double,double>& TranslationTableY = mr_model_part.GetTable(mtable_ids[1]);
	Table<double,double>& TranslationTableZ = mr_model_part.GetTable(mtable_ids[2]);

	array_1d<double,3> translation = ZeroVector(3);
	translation(0)=TranslationTableX(time);
	translation(1)=TranslationTableY(time);
	translation(2)=TranslationTableZ(time);
	
	
	for (unsigned int mesh_index=0;mesh_index<mgroup_ids.size();mesh_index++) //we loop around the desired groups
	{
      
		const int mesh_id=mgroup_ids[mesh_index];
		ModelPart::MeshType& current_mesh = mr_model_part.GetMesh(mesh_id);
		ModelPart::NodesContainerType::iterator inodebegin = current_mesh.NodesBegin();

		#pragma omp parallel 
        
        #pragma omp for
		for(int ii=0; ii< static_cast<int>(current_mesh.Nodes().size()); ii++)
		{
                    
			ModelPart::NodesContainerType::iterator pnode = inodebegin+ii;
            
			pnode->Coordinates() = pnode->GetInitialPosition().Coordinates()+translation;
            
			if (pnode->SolutionStepsDataHas(DISPLACEMENT_X))
            {

				pnode->FastGetSolutionStepValue(DISPLACEMENT) = translation;
            }

		}
	}
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

#endif // KRATOS_TRANSLATION_OPERATION_H_INCLUDED  defined 


