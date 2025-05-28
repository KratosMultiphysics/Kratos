//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#if !defined(KRATOS_BREP_SBM_UTILITIES_H_INCLUDED)
#define KRATOS_BREP_SBM_UTILITIES_H_INCLUDED

// Std includes
#include <list>

// System includes
#include "includes/define.h"

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"

#include "includes/node.h"

#include "includes/model_part.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    //using namespace ClipperLib;

    class KRATOS_API(KRATOS_CORE) BrepSBMUtilities
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef std::size_t IndexType;
        typedef std::size_t SizeType;

        typedef IntegrationPoint<3> IntegrationPointType;
        typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;


        //template<class TBrepLoopType, class TPointType>
        static void CreateBrepSurfaceSBMIntegrationPoints(
            IntegrationPointsArrayType& rIntegrationPoints,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            ModelPart& rSurrogateModelPart_inner, 
            ModelPart& rSurrogateModelPart_outer,
            IntegrationInfo& rIntegrationInfo);

        static void CreateBrepSurfaceSBMExternalIntegrationPoints(
            IntegrationPointsArrayType& rIntegrationPoints,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            ModelPart& rSurrogateModelPart_inner, 
            ModelPart& rSurrogateModelPart_outer,
            IntegrationInfo& rIntegrationInfo);

        
        static void CreateBrepVolumeSBMIntegrationPoints(
            IntegrationPointsArrayType& rIntegrationPoints,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            const std::vector<double>& rSpansW,
            ModelPart& rSurrogateModelPart_inner, 
            ModelPart& rSurrogateModelPart_outer,
            IntegrationInfo& rIntegrationInfo);
    
    private:
    
        static int FindKnotSpans1D(
            const std::vector<double>& rSpans, const double coord);

        ///@}
    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_UTILITIES_H_INCLUDED defined
