//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_CHECK_CONSTITUTIVELAW_UTILITY_H_INCLUDED )
#define  KRATOS_CHECK_CONSTITUTIVELAW_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/constitutive_law.h"
#include "geometries/triangle_2d_3.h"

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
    class CheckConstitutiveLawUtility
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of CheckConstitutiveLawUtility
        KRATOS_CLASS_POINTER_DEFINITION(CheckConstitutiveLawUtility);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        CheckConstitutiveLawUtility(Properties::Pointer pprop, ProcessInfo& r_process_info, double characteristic_size)
        {
            mp_property = pprop;
//            mr_model_part.pGetProperties(property_index,0);
            mp_constitutive_law = (mp_property->GetValue(CONSTITUTIVE_LAW))->Clone();
            
            //create new nodes
            Node < 3 > ::Pointer pnode1 = Node < 3 > ::Pointer(new Node<3>(1, 0.0, 0.0, 0.0 ));
            Node < 3 > ::Pointer pnode2 = Node < 3 > ::Pointer(new Node<3>(2, characteristic_size, 0.0, 0.0 ));
            Node < 3 > ::Pointer pnode3 = Node < 3 > ::Pointer(new Node<3>(3, characteristic_size, characteristic_size, 0.0 ));

            //create new geometry
            mpgeom = Geometry<Node < 3 > >::Pointer(new Triangle2D3<Node < 3 > >(pnode1,pnode2,pnode3));
            
            mp_constitutive_law->Check(*mp_property,*mpgeom,r_process_info);
        }

        /// Destructor.

        virtual ~CheckConstitutiveLawUtility()
        {
        }


        ///@}
        ///@name Operators
        ///@{

        void CheckConstitutiveLaw(
                Vector& strain,
                Vector& stress,
                Matrix& C,
                ProcessInfo& rCurrentProcessInfo,
                bool CalculateStresses,
                int CalculateTangent,
                bool SaveInternalVariables)
        {
            KRATOS_TRY
            
            unsigned int strain_size = mp_constitutive_law->GetStrainSize();

            Vector N(strain_size,0.0);
            Matrix F(2,2);

            if(strain.size()!=strain_size)
                KRATOS_ERROR(std::logic_error,"wrong size of strain ... expected = ",strain_size)
            if(stress.size()!=strain_size)
                KRATOS_ERROR(std::logic_error,"wrong size of stress ... expected = ",strain_size)
            if(C.size2()!=strain_size)
                KRATOS_ERROR(std::logic_error,"wrong size of C ... expected = ",strain_size)

            //compute continuous stress
            mp_constitutive_law->CalculateMaterialResponse(strain,F,stress,C,rCurrentProcessInfo,*mp_property,*mpgeom,N,CalculateStresses,CalculateTangent,SaveInternalVariables);


            KRATOS_CATCH("")
        }

        double GetScalarVar(Variable<double>& rVariable)
        {
            double value = 0.0;
            double tmp = mp_constitutive_law->GetValue(rVariable,value);
            return tmp;
        }


        void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
        {
            unsigned int strain_size = mp_constitutive_law->GetStrainSize();
            Vector N(strain_size,0.0);
            mp_constitutive_law->InitializeSolutionStep(*mp_property,*mpgeom,N,rCurrentProcessInfo);
        }

        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
        {
            unsigned int strain_size = mp_constitutive_law->GetStrainSize();
            Vector N(strain_size,0.0);
            mp_constitutive_law->FinalizeSolutionStep(*mp_property,*mpgeom,N,rCurrentProcessInfo);
        }
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
            buffer << "CheckConstitutiveLawUtility";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "CheckConstitutiveLawUtility";
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
        ConstitutiveLaw::Pointer mp_constitutive_law;
        Geometry<Node < 3 > >::Pointer mpgeom;
        Properties::Pointer mp_property;
        
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
        //      CheckConstitutiveLawUtility& operator=(CheckConstitutiveLawUtility const& rOther){}

        /// Copy constructor.

//        CheckConstitutiveLawUtility(CheckConstitutiveLawUtility const& rOther)
//        {
//        }


        ///@}

    }; // Class CheckConstitutiveLawUtility

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
            CheckConstitutiveLawUtility& rThis)
    {
        return rIStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const CheckConstitutiveLawUtility& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CHECK_CONSTITUTIVELAW_UTILITY_H_INCLUDED  defined


