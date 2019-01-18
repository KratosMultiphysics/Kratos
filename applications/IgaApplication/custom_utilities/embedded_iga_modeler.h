#if !defined(KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED

// System includes

// External includes
#include "anurbs.h"
#include "containers/model.h"
// #include "triangle.h"

// Project includes
#include "iga_application_variables.h"
#include "nurbs_brep_modeler.h"
#include "includes/model_part.h"

namespace Kratos
{
    struct DPState 
    {
        bool visible;
        double weight;
        long bestvertex;
    };
    
    struct Diagonal {
        long index1;
        long index2;
    };

    class EmbeddedIgaModeler : public NurbsBrepModeler
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaModeler);

        ///@}
        ///@name functions
        ///@{
        
        void CreateTessellationCurve(ANurbs::Pointer<ANurbs::CurveTessellation3D>& rTessellation); 
        void CreateTessellationParameterCurve(std::vector<array_1d<double, 3> >& rPolygon); 
        void CreateElements2D(ModelPart& rSkinModelPart);
        bool CreateElements3D(const std::vector<array_1d<double,3> >& rpolygon, std::vector<Matrix>& rtriangles); 

        std::vector<Matrix> Triangulate(); 
        
        bool IsConvex(
            const array_1d<double, 3>& p1, const array_1d<double, 3>& p2, 
            const array_1d<double, 3>& p3); 
        
        bool InCone(
            array_1d<double, 3> &p1, array_1d<double, 3> &p2,
            array_1d<double, 3> &p3, array_1d<double, 3> &p) ; 
        
        bool Intersects(
            array_1d<double, 3>& p11, array_1d<double, 3>& p12,
            array_1d<double, 3>& p21, array_1d<double, 3>& p22); 

        double Distance(array_1d<double, 3> p1, array_1d<double, 3> p2); 
        double GetAreaOfTriangle(const Matrix& triangle);


        std::vector<double> PrintNodesX();
        std::vector<double> PrintNodesY();
        std::vector<double> PrintNodesX3D();
        std::vector<double> PrintNodesY3D();
        std::vector<double> TessellationX(); 
        std::vector<double> TessellationY(); 

        ///@}
        ///@name Life Cycle
        ///@{
        /// Constructor.
        EmbeddedIgaModeler(ModelPart& rModelPart);

        /// Destructor.
        virtual ~EmbeddedIgaModeler() override
        {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        ModelPart&                 m_model_part;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    };

}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_IGA_H_INCLUDED defined


