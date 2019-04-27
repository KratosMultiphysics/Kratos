#if !defined(KRATOS_EMBEDDED_IGA_ERROR_ESTIMATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_ERROR_ESTIMATION_H_INCLUDED

// System includes
#include <cmath>
// External includes

// Project includes
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "iga_application_variables.h"
#include "custom_utilities/embedded_iga/embedded_iga_mapper.h"

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
    static void InsertGaussPointsExactSurface(
        const BrepFace& rFaceGeometry,
        const std::vector<Matrix>& rTriangulation_uv,
        std::vector<Matrix>& rGaussPoints_xyz);
    
    static void InsertGaussPointsApproxSurface(
        const BrepFace& rFaceGeometry,
        const std::vector<Matrix>& rTriangulation_uv,
        std::vector<Matrix>& rGaussPoints_xyz);
    
    static void GetError(
        const std::vector<Matrix>& rGaussPointsExact, 
        const std::vector<Matrix>& rGaussPointsApprox, 
        Vector& rError);
        
    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor.
    EmbeddedIgaErrorEstimation(); 

    /// Destructor.
    virtual ~EmbeddedIgaErrorEstimation() 
    {};

    ///@}
protected:

private:
    ///@name Member Variables
    ///@{

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


