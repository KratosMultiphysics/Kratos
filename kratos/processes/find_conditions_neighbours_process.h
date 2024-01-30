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


#if !defined(KRATOS_FIND_CONDITIONS_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_CONDITIONS_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ConditionsContainerType ConditionsContainerType;


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
class FindConditionsNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindConditionsNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindConditionsNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindConditionsNeighboursProcess(ModelPart& model_part, int TDim, unsigned int avg_conds = 10)
        : mr_model_part(model_part)
    {
        mavg_conds = avg_conds;
        mTDim=TDim;
// 	mavg_nodes = avg_nodes;
    }

    /// Destructor.
    ~FindConditionsNeighboursProcess() override
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
        NodesContainerType& rNodes = mr_model_part.Nodes();
        ConditionsContainerType& rConds = mr_model_part.Conditions();

        //first of all the neighbour nodes and conditions array are initialized to the guessed size
        //and empties the old entries
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            (in->GetValue(NEIGHBOUR_CONDITIONS)).reserve(mavg_conds);
            GlobalPointersVector<Condition >& rC = in->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end() );
        }
        for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ic++)
        {
            (ic->GetValue(NEIGHBOUR_CONDITIONS)).reserve(3);
            GlobalPointersVector<Condition >& rC = ic->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end() );
        }

        //add the neighbour conditions to all the nodes in the mesh
        for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ic++)
        {
            Condition::GeometryType& pGeom = ic->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                //KRATOS_WATCH( pGeom[i] );
                (pGeom[i].GetValue(NEIGHBOUR_CONDITIONS)).push_back( Condition::WeakPointer( *(ic.base()) ) );
                //KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_CONDITIONS)).size() );
            }
        }

        //adding the neighbouring conditions to the condition
        //loop over faces
        if (mTDim==3)
        {
            for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ic++)
            {
                //face nodes
                Geometry<Node >& geom = (ic)->GetGeometry();
                //vector of the 3 faces around the given face
                (ic->GetValue(NEIGHBOUR_CONDITIONS)).resize(3);
                GlobalPointersVector< Condition >& neighb_faces = ic->GetValue(NEIGHBOUR_CONDITIONS);
                //neighb_face is the vector containing pointers to the three faces around ic
                //neighb_face[0] = neighbour face over edge 1-2 of element ic;
                //neighb_face[1] = neighbour face over edge 2-0 of element ic;
                //neighb_face[2] = neighbour face over edge 0-1 of element ic;
                neighb_faces(0) = CheckForNeighbourFaces(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_CONDITIONS), ic->Id());
                neighb_faces(1) = CheckForNeighbourFaces(geom[2].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_CONDITIONS), ic->Id());
                neighb_faces(2) = CheckForNeighbourFaces(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_CONDITIONS), ic->Id());
            }
        }
    }


    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            GlobalPointersVector<Condition >& rC = in->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end());
        }
        ConditionsContainerType& rConds = mr_model_part.Conditions();
        for(ConditionsContainerType::iterator ic = rConds.begin(); ic!=rConds.end(); ic++)
        {
            GlobalPointersVector<Condition >& rC = ic->GetValue(NEIGHBOUR_CONDITIONS);
            rC.erase(rC.begin(),rC.end());
        }

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
        return "FindConditionsNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindConditionsNeighboursProcess";
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
    unsigned int mavg_conds;
    int mTDim;
// 	unsigned int mavg_nodes;


    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    template< class TDataType > void  AddUniqueWeakPointer
    (GlobalPointersVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename GlobalPointersVector< TDataType >::iterator i = v.begin();
        typename GlobalPointersVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    Condition::WeakPointer CheckForNeighbourFaces (unsigned int Id_1, unsigned int Id_2, GlobalPointersVector< Condition >& neighbour_face, unsigned int face)
    {
        //look for the faces around node Id_1
        for( GlobalPointersVector< Condition >::iterator i =neighbour_face.begin(); i != neighbour_face.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node >& neigh_face_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_face_geometry.size(); node_i++)
            {
                if (neigh_face_geometry[node_i].Id() == Id_2)
                {
                    if(i->Id() != face)
                    {
                        return *(i.base());
                    }
                }
            }
        }
        return Condition::WeakPointer();
    }
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
    FindConditionsNeighboursProcess& operator=(FindConditionsNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindConditionsNeighboursProcess(FindConditionsNeighboursProcess const& rOther);


    ///@}

}; // Class FindConditionsNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindConditionsNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindConditionsNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_CONDITIONS_NEIGHBOURS_PROCESS_H_INCLUDED  defined 


