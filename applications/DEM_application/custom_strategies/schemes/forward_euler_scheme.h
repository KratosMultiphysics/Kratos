//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2012-04-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//



#if !defined(KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED )
#define  KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cfloat>

// External includes 


// Project includes
#include "integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "DEM_application.h"


namespace Kratos {

    class ForwardEulerScheme : public IntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of ForwardEulerScheme
        KRATOS_CLASS_POINTER_DEFINITION(ForwardEulerScheme);

        /// Default constructor.

        ForwardEulerScheme() {
        }

        /// Destructor.

        virtual ~ForwardEulerScheme() {
        }

        void Calculate(ModelPart& model_part);
        void CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes);
        void CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes);
        void CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part);

        /// Turn back information as a string.

        virtual std::string Info() const {
            std::stringstream buffer;
            buffer << "ForwardEulerScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
            rOStream << "ForwardEulerScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
        }


    protected:


    private:


        /// Assignment operator.

        ForwardEulerScheme& operator=(ForwardEulerScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        ForwardEulerScheme(ForwardEulerScheme const& rOther) {
            *this = rOther;
        }


        ///@}    

    }; // Class ForwardEulerScheme 

    ///@} 

    ///@name Type Definitions       
    ///@{ 


    ///@} 
    ///@name Input and output 
    ///@{ 


    /// input stream function

    inline std::istream& operator>>(std::istream& rIStream,
            ForwardEulerScheme& rThis) {
        return rIStream;
    }

    /// output stream function

    inline std::ostream& operator<<(std::ostream& rOStream,
            const ForwardEulerScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED  defined 
