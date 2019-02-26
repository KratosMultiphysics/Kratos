#if !defined(KRATOS_EMBEDDED_IGA_MAPPER_H_INCLUDED)
#define  KRATOS_EMBEDDED_IGA_MAPPER_H_INCLUDED



// System includes

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/brep_topology/brep_face.h"

namespace Kratos
{
class EmbeddedIgaMapper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaMapper);

    ///@}
    ///@name functions
    ///@{
    
    static void MapCartesianSpace(
        const BrepFace& rFaceGeometry,
        const std::vector<Matrix>& rPoints_uv,
              std::vector<Matrix>& rPoints_xyz) 
    {
        /**
         * This static class maps points from the parametric space into the cartesian space.
         * The number of points per row (element of the std::vector) is dependent on the
         * number of input points.
        */
        const auto number_triangles = rPoints_uv.size(); 
        const auto surface_geometry = rFaceGeometry.GetSurface();
        
        rPoints_xyz.resize(
            rPoints_uv.size(), ZeroMatrix(rPoints_uv[0].size1(),3));


        for (unsigned int tri_i = 0; tri_i < number_triangles; ++tri_i)
        {    
            for (unsigned int point_i = 0; point_i < 3; ++point_i)
            {
                auto point_xyz = surface_geometry->PointAt(
                        rPoints_uv[tri_i](point_i,0), rPoints_uv[tri_i](point_i,1));
                
                for (unsigned int j = 0; j < 3; ++j)    
                {
                    rPoints_xyz[tri_i](point_i, j) = point_xyz[j]; 
                }
            }
        }      
    }; 

    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor.
    EmbeddedIgaMapper();

    /// Destructor.
    virtual ~EmbeddedIgaMapper()
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
#endif // KRATOS_EMBEDDED_IGA_MAPPER_H_INCLUDED defined