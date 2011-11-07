//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_METIS_SCALAR_REORDER_PROCESS_H_INCLUDED )
#define  KRATOS_METIS_SCALAR_REORDER_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include <parmetis.h>
#include "includes/model_part.h"

extern "C"
{
    //extern void METIS_PartMeshDual(int*, int*, idxtype*, int*, int*, int*, int*, idxtype*, idxtype*);
    void METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *);
};

namespace Kratos
{
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

    /// This class reorders the nodes using metis (Fill in reordering)

    /** Detail class definition.
     */
    class MetisScalarReorder : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MetisScalarReorder
        KRATOS_CLASS_POINTER_DEFINITION(MetisScalarReorder);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        MetisScalarReorder(ModelPart& rModelPart) : mrModelPart(rModelPart)
        {
        }

        /// Destructor.

        virtual ~MetisScalarReorder()
        {
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        virtual void Execute()
        {
            KRATOS_TRY

            if (NEIGHBOUR_NODES.Key() == 0)
                KRATOS_ERROR(std::logic_error, "Metis Application not registered correctly. NEIGBHBOUR NODES has 0 key!", "");

            //first of all compute the size of the graph, and ensure that the nodes are ordered consecutively 1...N (with no gaps)
            //we should verify here that the neighbours are calculated
            int N = mrModelPart.Nodes().size();
            int total_neighbs = 0;
            int id = 1;

            for (ModelPart::NodesContainerType::iterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); it++)
            {
                WeakPointerVector<Node < 3 > >& nodal_neighb = it->GetValue(NEIGHBOUR_NODES);
                int nn = (nodal_neighb).size();
                if (nn == 0)
                    KRATOS_ERROR(std::logic_error, "isolated node found, or neighbours not calculated, ID= ", it->Id());
                //                    std::cout << "isolated node found, or neighbours not calculated, ID= " << it->Id() << std::endl;
                total_neighbs += nn ;

                it->SetId(id++);
            }
            KRATOS_WATCH(total_neighbs);
            KRATOS_WATCH(N);

            //allocate data needed for Metis
            idxtype* adjncy = new idxtype[total_neighbs];
            idxtype* xadj = new idxtype[N + 1];
            idxtype* perm = new idxtype[N];
            idxtype* iperm = new idxtype[N];
            int size8 = 8;
            int* options = new int[size8];

            //fill the graph
            int counter = 0;
            int row_counter = 0;
            for (ModelPart::NodesContainerType::iterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); it++)
            {
                WeakPointerVector<Node < 3 > >& nodal_neighb = it->GetValue(NEIGHBOUR_NODES);

                xadj[row_counter++] = counter;

                for (WeakPointerVector<Node < 3 > >::iterator iii = nodal_neighb.begin(); iii != nodal_neighb.end(); iii++)
                {
                    if (iii->Id() != it->Id())
                        adjncy[counter++] = iii->Id() - 1;
                }

            }

            xadj[row_counter++] = counter; //finalize the array with the N+1 entry

                    //call metis
                    int numflag = 0;
            options[0] = 0;
            METIS_NodeND(&N, xadj, adjncy, &numflag, options, perm, iperm);


            //perform the permutation
            counter = 0;
            for (ModelPart::NodesContainerType::iterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); it++)
            {
                it->SetId(iperm[counter++] + 1);
            }
            mrModelPart.Nodes().Sort();

                    //deallocate mem
                    delete [] xadj;
            delete [] adjncy;
            delete [] perm;
            delete [] iperm;
            delete [] options;



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

        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "MetisScalarReorder";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "MetisScalarReorder";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
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

        ModelPart& mrModelPart;


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

        MetisScalarReorder & operator=(MetisScalarReorder const& rOther)
        {
            this->mrModelPart = rOther.mrModelPart;
            return *this;
        }

        /// Copy constructor.
        //
        //        MetisScalarReorder(MetisScalarReorder const& rOther)
        //        {
        //        }


        ///@}

    }; // Class MetisScalarReorder

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
            MetisScalarReorder& rThis)
    {
        return rIStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const MetisScalarReorder& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_METIS_SCALAR_REORDER_PROCESS_H_INCLUDED  defined


