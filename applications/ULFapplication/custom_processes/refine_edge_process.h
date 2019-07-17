//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov

#if !defined(KRATOS_REFINE_EDGE_PROCESS_INCLUDED )
#define  KRATOS_REFINE_EDGE_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "custom_utilities/geometry_utilities2D.h"



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
	Update the PRESSURE_FORCE on the nodes


*/

class RefineEdgeProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(RefineEdgeProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineEdgeProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    ~RefineEdgeProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{




    ///@}
    ///@name Operations
    ///@{

    void RefineFreeSurfEdge2D(const double & factor)
    {
        KRATOS_TRY

        int step_data_size = mr_model_part.GetNodalSolutionStepDataSize();
        
	//REFINING EDGES FOR GLASS SIMULATION PROCESSES
	//NOTE: INBLOW SURFACE IS MARKED WITH FLAG_VARIABLE=1 and OUTER SURFACE IS MARKED WITH IS_INTERFACE=1.. BOTH ARE IS_FREE_SURFACE
	int id = (mr_model_part.Nodes().end() - 1)->Id() + 1;	
	//KRATOS_WATCH(id)	
	for(ModelPart::ConditionsContainerType::iterator ic = mr_model_part.ConditionsBegin() ; ic != mr_model_part.ConditionsEnd() ; ic++)
	{

		if (ic->GetGeometry().size()==2)
		{
		unsigned int n_fs=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE);			
		n_fs+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE);

		unsigned int n_str=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);			
		n_str+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);

								
		//THIS REFINES THE NODES OF INTERBAL ELEMENTS OF THE SURFACE WHERE THE INBLOW IS: FLAG_VAR=1
		if (n_fs==ic->GetGeometry().size() ) //|| (n_fs==1.0 && n_str==1.0 ))
				{
				double x0=ic->GetGeometry()[0].X(); 	double y0=ic->GetGeometry()[0].Y(); 	
				double x1=ic->GetGeometry()[1].X(); 	double y1=ic->GetGeometry()[1].Y(); 	


				double edge=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));

				double nodal_h=ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
				nodal_h+=ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
				
				nodal_h*=0.5;
				//if the edge of the facet (condition) is too long, we split it into two by adding a node in the middle
				

				Node<3>::DofsContainerType& reference_dofs = (mr_model_part.NodesBegin())->GetDofs();
				

				if (edge>factor*nodal_h)
					{
					id++;

					double x=0.5*(x0+x1); 	double y=0.5*(y0+y1);	double z=0.0;
					Node<3>::Pointer pnode = mr_model_part.CreateNewNode(id,x,y,z);

					//putting the new node also in an auxiliary list
					//KRATOS_WATCH("adding nodes AT THE EDGE==================================")
																							
					
					//std::cout << "new node id = " << pnode->Id() << std::endl;
					//generating the dofs
					for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
						{
						Node<3>::DofType& rDof = *iii;
						Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
						
						(p_new_dof)->FreeDof();
						}
					Geometry<Node<3> >& geom = ic->GetGeometry();
					InterpolateOnEdge(geom, step_data_size, pnode);
					const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
					pnode->X0() = pnode->X() - disp[0];
					pnode->Y0() = pnode->Y() - disp[1];
					pnode->Z0() = 0.0;
					//KRATOS_WATCH("Added node at the EDGE=================================================================")
					}

				

				}
			}
		}
        //int new_last_id = (mr_model_part.Nodes().end() - 1)->Id() + 1;	
	//KRATOS_WATCH(new_last_id)
        KRATOS_CATCH("")
    }

	void InterpolateOnEdge( Geometry<Node<3> >& geom, unsigned int step_data_size, Node<3>::Pointer pnode)
		{
		unsigned int buffer_size = pnode->GetBufferSize();
		//KRATOS_WATCH(buffer_size)

		for(unsigned int step = 0; step<buffer_size; step++)
		{	
				
			//getting the data of the solution step
			double* step_data = (pnode)->SolutionStepData().Data(step);		
			
											
			double* node0_data = geom[0].SolutionStepData().Data(step);
			double* node1_data = geom[1].SolutionStepData().Data(step);						
				
			//copying this data in the position of the vector we are interested in
			for(unsigned int j= 0; j<step_data_size; j++)
				{ 

					step_data[j] = 0.5*(node0_data[j] + node1_data[j]);
								

				}						
			}
			//now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
			//we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
			//not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way 		
			//as the normal variables, like Velocity etc		
			
			//pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
			//pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=1.0;			

			pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)=0.0;
			pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=1.0;
			pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;
			

			pnode->Set(TO_ERASE,false);
			
			if (pnode->FastGetSolutionStepValue(IS_INTERFACE)>0.99)
				pnode->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
			else 
				pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;

							

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
        return "RefineEdgeProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RefineEdgeProcess";
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
    ModelPart& mr_model_part;

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
//		RefineEdgeProcess& operator=(RefineEdgeProcess const& rOther);

    /// Copy constructor.
//		RefineEdgeProcess(RefineEdgeProcess const& rOther);


    ///@}

}; // Class RefineEdgeProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RefineEdgeProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineEdgeProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_EDGE_PROCESS_INCLUDED  defined 


