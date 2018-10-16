//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                    
//



#if !defined(KRATOS_BOUNDING_BOX_UTILITIES_INCLUDED )
#define  KRATOS_BOUNDING_BOX_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <set>
#include <time.h>



#ifdef _OPENMP
#include <omp.h>
#endif

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/mesh.h"

#include "geometries/geometry.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"
#include "utilities/geometry_utilities.h"
#include "utilities/timer.h"

namespace Kratos
{

///******************************************************************************************************************
///******************************************************************************************************************

template <std::size_t TDimension>
class BoxFunction
{
public:
    template< class TPointType, class TPointerType>
    void operator ()(TPointerType& rObject, TPointType& rLowPoint, TPointType& rHighPoint)
    {
        rHighPoint = rObject.GetGeometry().GetPoint(0);
        rLowPoint  = rObject.GetGeometry().GetPoint(0);
        for (unsigned int point = 0; point<rObject.GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<TDimension; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject.GetGeometry().GetPoint(point)[i] ) ?  rObject.GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject.GetGeometry().GetPoint(point)[i] ) ?  rObject.GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }
};



///******************************************************************************************************************
///******************************************************************************************************************

class TriBoxOverlapFunction
{
public:
    template< class TPointType, class TPointerType>
    bool operator ()(TPointerType& rObject,  const TPointType& rLowPoint, const TPointType& rHighPoint)
    {
        return rObject.GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }
};


///******************************************************************************************************************
///******************************************************************************************************************

class TDistanceFunction
{
public:
    template<class TPointerType>
    bool operator ()(TPointerType& rObj_1, TPointerType& rObj_2)
    {
        Element::GeometryType& geom_1 = rObj_1.GetGeometry();
        Element::GeometryType& geom_2 = rObj_2.GetGeometry();
        return  geom_1.HasIntersection(geom_2);
    }
};

template<class TPointType, std::size_t TDimension>
class Segment : public Point
{
    enum {Dimension = TDimension };
    typedef Point    PointType;
    typedef array_1d<double, Dimension> VectorType;
public:

    Segment(const PointType& rPoint1, const PointType& rPoint2) :
        mPoint1(rPoint1), mPoint2(rPoint2)
    {
    }

    PointType Center()
    {
        return 0.50 * (mPoint1 + mPoint2);
    }

    double Length()
    {
        return norm_2(Direction());
    }

    VectorType Direction()
    {
        return mPoint2 - mPoint1;
    }

    double Extent()
    {
        return    0.50 * Length();
    }

    void Normalize(const double epsilon = 1E-9)
    {
        const double length = Length();
        VectorType result   = Direction();
        if (length > epsilon)
        {
            const double invLength = 1.00 / length;
            for(std::size_t i = 0; i<Dimension; i++)
                result[i] =  result * invLength;
        }
        else
        {
            for(std::size_t i = 0; i<Dimension; i++)
                result[i] =  0.00;
        }

    }


    PointType  mPoint1;
    PointType  mPoint2;



};


///******************************************************************************************************************
///******************************************************************************************************************


class BoundingBoxUtilities
{
public:


    const static std::size_t dimension = 2;
    typedef Point                                            PointType;
    typedef Segment<PointType, 2>                            SegmentType;
    typedef SegmentType*                                     SegmentPointer;
    typedef std::vector<SegmentType>                         ContainerSegmentType;
    typedef ModelPart::ElementsContainerType::ContainerType  ContainerType;
    typedef ContainerType::value_type                        PointerType;
    typedef ContainerType::iterator                          IteratorType;
    typedef SpatialContainersConfigure<dimension>            Configure2D;
    typedef Cell<Configure2D>                                CellType;
    typedef std::vector<CellType>                            CellContainerType;
    typedef CellContainerType::iterator                      CellContainerIterator;
    typedef std::vector<PointerType>::iterator               PointerTypeIterator;
    typedef ContactPair<PointerType>                         ContactPairType;
    typedef std::vector<ContactPairType>                     ContainerContactPair;
    typedef ContainerContactPair::iterator                   IteratorContainerContactPair;
    typedef ContainerContactPair::value_type                 PointerContainerContactPair;
    typedef Element::GeometryType                            GeomType;
    typedef Node<3>                                          NodeType;

    BoundingBoxUtilities(ModelPart& model_part, const unsigned int& dimension) : mr_model_part(model_part), mrdimension(dimension)
    {
    }

    virtual ~BoundingBoxUtilities() {}

    void Test()
    {
    }
};




}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


