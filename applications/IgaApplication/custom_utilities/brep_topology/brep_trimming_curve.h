#if !defined(KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED )
#define  KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED


// System includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"


namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    class BrepTrimmingCurve
    {
    public:
        ///@name Type Definitions
        ///@{
        //typedef std::vector<array_1d<double, 4>> ControlPointVector;

        /// Pointer definition of KratosNurbsBrepApplication
        //KRATOS_CLASS_POINTER_DEFINITION(BrepTrimmingCurve);

        ///@}
        ///@name Life Cycle 
        ///@{ 
        // Utilities
        /* Returns trimming curve index
        */
        unsigned int& GetTrimIndex();

        //TODO: you need to give reading access to your internals through the Calculate function
        /// Constructor.
        //TODO: pass by reference not by value
        BrepTrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
            unsigned int p,
            Vector& active_range);

        /// Destructor.
        virtual ~BrepTrimmingCurve() {};

        ///@} 
    protected:

    private:
        ///@name Member Variables
        ///@{ 
        unsigned int m_trim_index;
        bool m_curve_direction;
        Vector m_knot_vector_u;
        unsigned int m_p;
        //std::vector<array_1d<double, 4>> m_control_points;
        Vector m_active_range;
        ///@}
        ///@name Un accessible methods 
        ///@{



        ///@}

    }; // Class BrepTrimmingCurve 

}  // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED  defined