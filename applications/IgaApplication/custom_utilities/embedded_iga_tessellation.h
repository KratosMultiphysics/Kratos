#if !defined(KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED

// System includes

// External includes
#include "anurbs.h"

// Project includes
#include "iga_application_variables.h"
#include "nurbs_brep_modeler.h"

namespace Kratos
{
    class EmbeddedIgaTessellation
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaTessellation);

        ///@}
        ///@name functions
        ///@{
        
        void CreateTessellationCurve(
            std::vector<array_1d<double, 3> >& rPolygon);
        
        void CreateTessellationParameterCurve(
            std::vector<array_1d<double, 3> >& rPolygon);
        ///@}
        ///@name Life Cycle
        ///@{
        /// Constructor.
        EmbeddedIgaTessellation(std::vector<BrepModel>&  rBrepModelVector);

        /// Destructor.
        virtual ~EmbeddedIgaTessellation()
        {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        std::vector<BrepModel>  mBrepModelVector; 
        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    };

}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED defined


