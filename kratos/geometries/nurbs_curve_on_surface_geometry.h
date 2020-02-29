//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//                   Tobias Teschemacher
//                   Andreas Apostolatos
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos {

template <int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
class NurbsCurveOnSurfaceGeometry : public Geometry<typename TSurfaceContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<typename TSurfaceContainerPointType::value_type> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef NurbsSurfaceGeometry<3, TSurfaceContainerPointType> NurbsSurfaceType;
    typedef NurbsCurveGeometry<2, TCurveContainerPointType> NurbsCurveType;

    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of NurbsCurveOnSurfaceGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsCurveOnSurfaceGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    NurbsCurveOnSurfaceGeometry(
        typename NurbsSurfaceType::Pointer pSurface,
        typename NurbsCurveType::Pointer pCurve)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsSurface(pSurface)
        , mpNurbsCurve(pCurve)
    {
    }

    /// Default constructor
    NurbsCurveOnSurfaceGeometry()
        : BaseType(PointsArrayType(), &msGeometryData)
    {};

    /// Copy constructor
    NurbsCurveOnSurfaceGeometry(NurbsCurveOnSurfaceGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Copy constructor, with different point type.
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType> NurbsCurveOnSurfaceGeometry(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TOtherSurfaceContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mpNurbsSurface(rOther.mpNurbsSurface)
        , mpNurbsCurve(rOther.mpNurbsCurve)
    {
    }

    /// Destructor
    ~NurbsCurveOnSurfaceGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    NurbsCurveOnSurfaceGeometry& operator=(const NurbsCurveOnSurfaceGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }

    /**
     * @brief Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherCurveContainerPointType, class TOtherSurfaceContainerPointType>
    NurbsCurveOnSurfaceGeometry& operator=(
        NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TOtherCurveContainerPointType, TOtherSurfaceContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsSurface = rOther.mpNurbsSurface;
        mpNurbsCurve = rOther.mpNurbsCurve;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /*typename BaseType::Pointer Create(
        TSurfaceContainerPointType const& ThisPoints) const override
    {
        return Kratos::make_shared<NurbsCurveOnSurfaceGeometry>(ThisPoints);
    }*/

    ///@}
    ///@name Operation within Global Space
    ///@{

    /*
    * @brief This method maps from dimension space to working space.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        // Compute the coordinates of the embedded curve in the parametric space of the surface
        CoordinatesArrayType result_local = mpNurbsCurve->GlobalCoordinates(rResult, rLocalCoordinates);
        
        // Compute and return the coordinates of the surface in the geometric space
        return mpNurbsSurface->GlobalCoordinates(rResult, result_local);
    }

    /** 
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * From ANurbs library (https://github.com/oberbichler/ANurbs)
    * @param LocalCoordinates The local coordinates in dimension space
    * @param Derivative Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rCoordinates,
        const SizeType DerivativeOrder) const override
    {
        // Check size of output
        if (rGlobalSpaceDerivatives.size() != DerivativeOrder + 1) {
            rGlobalSpaceDerivatives.resize(DerivativeOrder + 1);
        }

        // Compute the gradients of the embedded curve in the parametric space of the surface
        std::vector<array_1d<double, 3>> curve_derivatives;
        mpNurbsCurve->GlobalSpaceDerivatives(curve_derivatives, rCoordinates, DerivativeOrder);
        
        // Compute the gradients of the surface in the geometric space
        array_1d<double, 3> surface_coordinates =  ZeroVector(3);
        surface_coordinates[0] = curve_derivatives[0][0];
        surface_coordinates[1] = curve_derivatives[0][1];
        std::vector<array_1d<double, 3>> surface_derivatives;
        mpNurbsSurface->GlobalSpaceDerivatives(surface_derivatives, surface_coordinates, DerivativeOrder);

        std::function<array_1d<double, 3>(int, int, int)> c;
        c = [&](int DerivativeOrder, int i, int j) -> array_1d<double, 3> {
            if (DerivativeOrder > 0) {
                array_1d<double, 3> result = ZeroVector(3);

                for (int a = 1; a <= DerivativeOrder; a++) {
                    result += (
                        c(DerivativeOrder - a, i + 1, j) * curve_derivatives[a][0] +
                        c(DerivativeOrder - a, i, j + 1) * curve_derivatives[a][1]
                        ) * NurbsUtilities::GetBinomCoefficient(DerivativeOrder - 1, a - 1);
                }

                return result;
            }
            else {
                const int index = NurbsSurfaceShapeFunction::IndexOfShapeFunctionRow(i, j);
                return surface_derivatives[index];
            }
        };
        for (SizeType i = 0; i <= DerivativeOrder; i++) {
            rGlobalSpaceDerivatives[i] = c(i, 0, 0);
        }
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    typename NurbsSurfaceType::Pointer mpNurbsSurface;
    typename NurbsCurveType::Pointer mpNurbsCurve;

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("pNurbsSurface", mpNurbsSurface);
        rSerializer.save("pNurbsCurve", mpNurbsCurve);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("pNurbsSurface", mpNurbsSurface);
        rSerializer.load("pNurbsCurve", mpNurbsCurve);
    }

    ///@}

}; // class NurbsCurveOnSurfaceGeometry

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryData NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TCurveContainerPointType, class TSurfaceContainerPointType>
const GeometryDimension NurbsCurveOnSurfaceGeometry<TWorkingSpaceDimension, TCurveContainerPointType, TSurfaceContainerPointType>::msGeometryDimension(
    1, TWorkingSpaceDimension, 1);

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_ON_SURFACE_H_INCLUDED defined