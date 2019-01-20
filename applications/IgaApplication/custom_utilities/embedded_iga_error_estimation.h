#if !defined(KRATOS_EMBEDDED_IGA_ERROR_ESTIMATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_ERROR_ESTIMATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "iga_application_variables.h"


namespace Kratos
{

    class EmbeddedIgaErrorEstimation
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaErrorEstimation);

        ///@}
        ///@name functions
        ///@{
        void PrintTriangleGaussPoints(); 
        ///@}
        ///@name Life Cycle
        ///@{
        /// Constructor.
        EmbeddedIgaErrorEstimation(std::vector<Matrix> rTriangles); 

        /// Destructor.
        virtual ~EmbeddedIgaErrorEstimation() 
        {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        
        std::vector<Matrix> mTriangles;
        

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    };

}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_IGA_ERROR_ESTIMATION_H_INCLUDED defined


