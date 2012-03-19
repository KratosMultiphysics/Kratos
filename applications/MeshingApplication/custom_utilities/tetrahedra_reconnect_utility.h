//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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

/// Short class definition.
/** Detail class definition.
*/
class TetrahedraReconnectUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TetrahedraReconnectUtility
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedraReconnectUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TetrahedraReconnectUtility() {}

    /// Destructor.
    virtual ~TetrahedraReconnectUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    /// this function evaluates the quality of all the tetrahedras within a given model_part
    void EvaluateQuality(ModelPart& r_model_part)
    {
        //loop on nodes
        for (ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            std::cout << it->Id() << std::endl;
        }

        for (ModelPart::ElementsContainerType::iterator el_it=r_model_part.ElementsBegin(); el_it!=r_model_part.ElementsEnd(); el_it++)
        {
            std::cout << el_it->GetGeometry().Area() << std::endl;
			
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			for(unsigned int i=0; i<geom.size(); i++)
			{
				const double x = geom[i].X();
				const double Y = geom[i].Y();
				const double z = geom[i].Z();
				geom[i].Id();
				array_1d<double,3>& coords = geom[i].Coordinates();
			}
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "TetrahedraReconnectUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TetrahedraReconnectUtility";
    }

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
    TetrahedraReconnectUtility& operator=(TetrahedraReconnectUtility const& rOther) {}

    /// Copy constructor.
    TetrahedraReconnectUtility(TetrahedraReconnectUtility const& rOther) {}


    ///@}

}; // Class TetrahedraReconnectUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TetrahedraReconnectUtility& rThis) {}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TetrahedraReconnectUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_RECONNECT_H_INCLUDED  defined 


