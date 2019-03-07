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



#if !defined(KRATOS_INITIALIZE_NODAL_AREA_PROCESS_H_INCLUDED )
#define  KRATOS_INITIALIZE_NODAL_AREA_PROCESS_H_INCLUDED



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
#include "utilities/geometry_utilities.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;


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
class InitializeNodalAreaProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InitializeNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(InitializeNodalAreaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    InitializeNodalAreaProcess(ModelPart& model_part, unsigned int domain_size)
        : mr_model_part(model_part), mdomain_size(domain_size)
    {
    }

    /// Destructor.
    ~InitializeNodalAreaProcess() override
    {
    }


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

    void Execute() override
    {
        KRATOS_TRY

        std::cout<<"                        INITIALIZE NODAL AREA PROCESS"<<std::endl;
        //set to zero the nodal area
        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();
                in!=mr_model_part.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
            in->FastGetSolutionStepValue(FLUID_FRACTION_RATE) = 0.0;


            if(in->Is(RIGID) || in->Is(ISOLATED)){
           	    WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);
                unsigned int numberOfNeighbourElements= neighb_elems.size();
                if(numberOfNeighbourElements==0 || in->Is(ISOLATED)){
                    // std::cout<<"node "<<in->Id()<<" is an isolated node "<<in->Coordinates()<<std::endl;
                   in->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
                   in->FastGetSolutionStepValue(FLUID_FRACTION) = 1.0;
                   in->FastGetSolutionStepValue(FLUID_FRACTION_OLD) = 1.0;
                   in->FastGetSolutionStepValue(FLUID_FRACTION_RATE) = 0.0;
                }  

            }
        }

        if(mdomain_size == 2)
        {
            double area = 0.0;
            for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin();
                    i!=mr_model_part.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();

                area = GeometryUtils::CalculateVolume2D(geom);
                area *= 0.333333333333333333333333333;


                geom[0].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[1].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[2].FastGetSolutionStepValue(NODAL_AREA) += area;
            }
        }
        else if(mdomain_size == 3)
        {
            for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin();
                    i!=mr_model_part.ElementsEnd(); i++)
            {
                double vol;
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();

                vol = GeometryUtils::CalculateVolume3D(geom);
                vol *= 0.25;

                geom[0].FastGetSolutionStepValue(NODAL_AREA) += vol;
                geom[1].FastGetSolutionStepValue(NODAL_AREA) += vol;
                geom[2].FastGetSolutionStepValue(NODAL_AREA) += vol;
                geom[3].FastGetSolutionStepValue(NODAL_AREA) += vol;
            }
        }

        mr_model_part.GetCommunicator().AssembleCurrentData(NODAL_AREA);



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
        return "InitializeNodalAreaProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InitializeNodalAreaProcess";
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
    unsigned int mdomain_size;


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
    InitializeNodalAreaProcess& operator=(InitializeNodalAreaProcess const& rOther);

    /// Copy constructor.
    //InitializeNodalAreaProcess(InitializeNodalAreaProcess const& rOther);


    ///@}

}; // Class InitializeNodalAreaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InitializeNodalAreaProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InitializeNodalAreaProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INITIALIZE_NODAL_AREA_PROCESS_H_INCLUDED  defined 


