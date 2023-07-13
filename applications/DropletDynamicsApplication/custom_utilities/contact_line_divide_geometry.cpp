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
#include "contact_line_divide_geometry.h"

namespace Kratos
{
    /// ContactLineIndexedPoint class implementation
    /// Constructors
    ContactLineIndexedPoint::ContactLineIndexedPoint()
        : Point() , IndexedObject(0) {};

    ContactLineIndexedPoint::ContactLineIndexedPoint(const unsigned int Id)
        : Point() , IndexedObject(Id) {};

    ContactLineIndexedPoint::ContactLineIndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id)
        : Point(rCoords) , IndexedObject(Id) {};

    /// Destructor
    ContactLineIndexedPoint::~ContactLineIndexedPoint() {};

    /// Turn back information as a string.
    std::string ContactLineIndexedPoint::Info() const {
        std::stringstream info_string;
        info_string << "Indexed point class created as a combination of the Point and IndexedObject classes.\n";
        info_string << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
        return info_string.str();
    };

    /// Print information about this object.
    void ContactLineIndexedPoint::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Indexed point class created as a combination of the Point and IndexedObject classes.\n";
        rOStream << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
    };

    /// Print object's data.
    void ContactLineIndexedPoint::PrintData(std::ostream& rOStream) const {
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
    ContactLineDivideGeometry<TPointType>::ContactLineDivideGeometry(const GeometryType& rInputGeometry, const Vector& rNodalDistances) : DivideGeometry<TPointType>(rInputGeometry, rNodalDistances),
        mrInputGeometry(rInputGeometry),
        mrNodalDistances(rNodalDistances) {
        this->IsSplit(); // Fast operation to state if the element is split or not.
    };

    /// Destructor
    template<class TPointType>
    ContactLineDivideGeometry<TPointType>::~ContactLineDivideGeometry() {};

    /// Turn back information as a string.
    template<class TPointType>
    std::string ContactLineDivideGeometry<TPointType>::Info() const {
        return "Base class for geometries splitting operations.";
    };

    /// Print information about this object.
    template<class TPointType>
    void ContactLineDivideGeometry<TPointType>::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations.";
    };

    /// Print object's data.
    template<class TPointType>
    void ContactLineDivideGeometry<TPointType>::PrintData(std::ostream& rOStream) const {
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
    Geometry<TPointType> ContactLineDivideGeometry<TPointType>::GetInputGeometry() const {
        return mrInputGeometry;
    };

    template<class TPointType>
    Vector ContactLineDivideGeometry<TPointType>::GetNodalDistances() const {
        return mrNodalDistances;
    };

    template<class TPointType>
    bool ContactLineDivideGeometry<TPointType>::IsSplit() {
        unsigned int n_pos = 0 , n_neg = 0;

        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            if (mrNodalDistances(i) < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0)) {
            mIsSplit = true;
        } else {
            mIsSplit = false;
        }
        return mIsSplit;
    };

    template<class TPointType>
    std::vector<typename ContactLineDivideGeometry<TPointType>::IndexedPointGeometryPointerType> ContactLineDivideGeometry<TPointType>::GetPositiveSubdivisions() const
    {
        return mPositiveSubdivisions;
    }

   template<class TPointType>
    std::vector<typename ContactLineDivideGeometry<TPointType>::IndexedPointGeometryPointerType> ContactLineDivideGeometry<TPointType>::GetNegativeSubdivisions() const
    {
        return mNegativeSubdivisions;
    }

    template<class TPointType>
    std::vector<typename ContactLineDivideGeometry<TPointType>::IndexedPointGeometryPointerType> ContactLineDivideGeometry<TPointType>::GetPositiveInterfaces() const
    {
        return mPositiveInterfaces;
    }

    template<class TPointType>
    std::vector<typename ContactLineDivideGeometry<TPointType>::IndexedPointGeometryPointerType> ContactLineDivideGeometry<TPointType>::GetNegativeInterfaces() const
    {
        return mNegativeInterfaces;
    }

    template<class TPointType>
    std::vector<unsigned int> ContactLineDivideGeometry<TPointType>::GetPositiveInterfacesParentIds() const
    {
        return mPositiveInterfacesParentIds;
    }

    template<class TPointType>
    std::vector<unsigned int> ContactLineDivideGeometry<TPointType>::GetNegativeInterfacesParentIds() const
    {
        return mNegativeInterfacesParentIds;
    }

    template class ContactLineDivideGeometry<Node>;
    template class ContactLineDivideGeometry<IndexedPoint>;
    
};
