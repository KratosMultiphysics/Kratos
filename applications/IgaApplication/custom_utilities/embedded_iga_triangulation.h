#if !defined(KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED

// System includes

// External includes
#include "anurbs.h"

// Project includes
#include "iga_application_variables.h"


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

    class EmbeddedIgaTriangulation
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaTriangulation);

        ///@}
        ///@name functions
        ///@{

        bool CreateTrianglesEmpire(
            const std::vector<array_1d<double,3> >& rpolygon, 
            std::vector<Matrix>& rtriangles); 

        bool IsConvex(
            const array_1d<double, 3>& p1, const array_1d<double, 3>& p2, 
            const array_1d<double, 3>& p3); 

        bool InCone(
            array_1d<double, 3> &p1, array_1d<double, 3> &p2,
            array_1d<double, 3> &p3, array_1d<double, 3> &p);

        bool Intersects(
            array_1d<double, 3>& p11, array_1d<double, 3>& p12,
            array_1d<double, 3>& p21, array_1d<double, 3>& p22);  

        double Distance(array_1d<double, 3> p1, array_1d<double, 3> p2); 

        double GetAreaOfTriangle(const Matrix& triangle);

        void PrintTriangleGaussPoints(); 

        ///@}
        ///@name Life Cycle
        ///@{
        /// Constructor.
        EmbeddedIgaTriangulation()
        {};

        /// Destructor.
        virtual ~EmbeddedIgaTriangulation() 
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
#endif // KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED defined


