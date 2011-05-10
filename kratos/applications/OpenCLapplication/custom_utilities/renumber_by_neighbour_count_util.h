//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_RENUMBER_BY_NEIGHBOUR_COUNT_UTIL_H_INCLUDED )
#define  KRATOS_RENUMBER_BY_NEIGHBOUR_COUNT_UTIL_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include "boost/smart_ptr.hpp"


// External includes 


// Project includes
#include "includes/define.h"


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

    /// This class renumbers the nodes grouping them according to the number of their neighbours

    /** Detail class definition.
     */
    class RenumberByNeighbourCountUtil
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of RenumberByNeighbourCountUtil
        KRATOS_CLASS_POINTER_DEFINITION(RenumberByNeighbourCountUtil);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        RenumberByNeighbourCountUtil()
        {
        }

        /// Destructor.

        virtual ~RenumberByNeighbourCountUtil()
        {
        }

        void Renumber(ModelPart::NodesContainerType& rNodes)
        {
            KRATOS_TRY

            //firs of all verify that the neighbours are calculated

            //now reorder the nodes list by a user defined ordering
            std::sort(rNodes.ptr_begin(), rNodes.ptr_end(), &Kratos::RenumberByNeighbourCountUtil::cmp );

            //set the ide by size
            int i = 1;
            for (ModelPart::NodesContainerType::iterator it = rNodes.begin(); it != rNodes.end(); it++)
            {
//                std::cout << it->GetValue(NEIGHBOUR_NODES).size() << std::endl;
                it->SetId(i++);
            }


            KRATOS_CATCH("");
        }
        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{


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
            buffer << "RenumberByNeighbourCountUtil";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "RenumberByNeighbourCountUtil";
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


        ///@}
        ///@name Private Operators
        ///@{

        //        struct compare_by_neighb_count
        //        {
        //
        //            bool operator ()(const Node<3>::Pointer it1, const Node<3>::Pointer it2) const
        //            {
        //                if (it1->GetValue(NEIGHBOUR_NODES).size() < it2->GetValue(NEIGHBOUR_NODES).size()) return true;
        //                return false;
        //            }
        ////            bool operator ()(ModelPart::NodesContainerType::iterator it1, ModelPart::NodesContainerType::iterator it2) const
        ////            {
        ////                if (it1->GetValue(NEIGHBOUR_NODES).size() < it2->GetValue(NEIGHBOUR_NODES).size()) return true;
        ////                return false;
        ////            }
        //
        //        };

        //        class compare_by_neighb_count : public std::binary_function< boost::shared_ptr<Node<3> >, boost::shared_ptr<Node<3> >, std::less<typename TGetKeyType::result_type>::result_type >
        //        {
        //     public:
        //       bool operator()(boost::shared_ptr<Node<3> > a, boost::shared_ptr<Node<3> > b) const
        //       {return ((*a).Id() < (*b).Id());}
        //     };

        static bool cmp(boost::shared_ptr<Node < 3 > > a, boost::shared_ptr<Node < 3 > > b)
        {
            return a->GetValue(NEIGHBOUR_NODES).size() < b->GetValue(NEIGHBOUR_NODES).size() ;
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

        RenumberByNeighbourCountUtil & operator=(RenumberByNeighbourCountUtil const& rOther)
        {
            return *this;
        }

        /// Copy constructor.

        RenumberByNeighbourCountUtil(RenumberByNeighbourCountUtil const& rOther)
        {
        }


        ///@}

    }; // Class RenumberByNeighbourCountUtil

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
            RenumberByNeighbourCountUtil& rThis)
    {
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const RenumberByNeighbourCountUtil& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RENUMBER_BY_NEIGHBOUR_COUNT_UTIL_H_INCLUDED  defined


