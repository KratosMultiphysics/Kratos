//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "utilities/divide_geometry.h"

namespace Kratos
{
    /// IndexedPoint class implementation
    /// Constructors
    IndexedPoint::IndexedPoint()
        : Point() , IndexedObject(0) {};

    IndexedPoint::IndexedPoint(const unsigned int Id)
        : Point() , IndexedObject(Id) {};

    IndexedPoint::IndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id)
        : Point(rCoords) , IndexedObject(Id) {};

    /// Destructor
    IndexedPoint::~IndexedPoint() {};

    /// Turn back information as a string.
    std::string IndexedPoint::Info() const {
        std::stringstream info_string;
        info_string << "Indexed point class created as a combination of the Point and IndexedObject classes.\n";
        info_string << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
        return info_string.str();
    };

    /// Print information about this object.
    void IndexedPoint::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Indexed point class created as a combination of the Point and IndexedObject classes.\n";
        rOStream << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
    };

    /// Print object's data.
    void IndexedPoint::PrintData(std::ostream& rOStream) const {
        rOStream << "Indexed point object:\n";
        rOStream << "\tIndex value: " << this->Id() << std::endl;

        const array_1d<double, 3> point_coords = this->Coordinates();
        std::stringstream coordinates_buffer;
        coordinates_buffer << "\tCoordinates: ( " << point_coords(0) << " , " << point_coords(1) << " , " << point_coords(2) << " )";
        rOStream << coordinates_buffer.str() << std::endl;
    };

    /// DivideGeometry implementation
    /// Default constructor
    template<class TPointType>
    DivideGeometry<TPointType>::DivideGeometry(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        mrInputGeometry(rInputGeometry),
        mrNodalDistances(rNodalDistances) {
        this->IsSplit(); // Fast operation to state if the element is split or not.
    };

    /// Destructor
    template<class TPointType>
    DivideGeometry<TPointType>::~DivideGeometry() {};

    /// Turn back information as a string.
    template<class TPointType>
    std::string DivideGeometry<TPointType>::Info() const {
        return "Base class for geometries splitting operations.";
    };

    /// Print information about this object.
    template<class TPointType>
    void DivideGeometry<TPointType>::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations.";
    };

    /// Print object's data.
    template<class TPointType>
    void DivideGeometry<TPointType>::PrintData(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations constructed with:\n";
        rOStream << "   Geometry type: " << mrInputGeometry.Info() << "\n";
        std::stringstream distances_buffer;
        std::ostringstream stm;
        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            stm << mrNodalDistances(i);
            distances_buffer << stm.str()<< " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    template<class TPointType>
    Geometry < TPointType > DivideGeometry<TPointType>::GetInputGeometry() const {
        return mrInputGeometry;
    };

    template<class TPointType>
    Vector DivideGeometry<TPointType>::GetNodalDistances() const {
        return mrNodalDistances;
    };

    template<class TPointType>
    void DivideGeometry<TPointType>::IsSplit() {
        unsigned int n_pos = 0 , n_neg = 0;

        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            if (mrNodalDistances(i) < 0.0) {
                n_neg++;
            } else if (mrNodalDistances(i) > 0.0) {
                n_pos++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0)) {
            mIsSplit = true;
        } else {
            mIsSplit = false;
        }
    };

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetPositiveSubdivisions() const
    {
        return mPositiveSubdivisions;
    }

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetNegativeSubdivisions() const
    {
        return mNegativeSubdivisions;
    }

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetPositiveInterfaces() const
    {
        return mPositiveInterfaces;
    }

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetNegativeInterfaces() const
    {
        return mNegativeInterfaces;
    }

    template<class TPointType>
    std::vector<unsigned int> DivideGeometry<TPointType>::GetPositiveInterfacesParentIds() const
    {
        return mPositiveInterfacesParentIds;
    }

    template<class TPointType>
    std::vector<unsigned int> DivideGeometry<TPointType>::GetNegativeInterfacesParentIds() const
    {
        return mNegativeInterfacesParentIds;
    }

    template<class TPointType>
    std::vector<unsigned int> DivideGeometry<TPointType>::GetContactInterface() const
    {
        return mContactInterface;
    }
    
    template<class TPointType>
    std::vector<unsigned int> DivideGeometry<TPointType>::GetContactEdge() const
    {
        return mContactEdge;
    }

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointGeometryPointerType> DivideGeometry<TPointType>::GetContactLine() const
    {
        return mContactLine;
    }

    template<class TPointType>
    std::vector<unsigned int> DivideGeometry<TPointType>::GetContactFace() const
    {
        return mContactFace;
    }

    template<class TPointType>
    std::vector<typename DivideGeometry<TPointType>::IndexedPointPointerType> DivideGeometry<TPointType>::GetContactPoint() const
    {
        return mContactPoint;
    }

    template class DivideGeometry<Node>;
    template class DivideGeometry<IndexedPoint>;
};
