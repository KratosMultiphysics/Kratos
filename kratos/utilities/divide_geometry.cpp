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
    DivideGeometry::DivideGeometry(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        mrInputGeometry(rInputGeometry),
        mrNodalDistances(rNodalDistances) {
        this->IsSplit(); // Fast operation to state if the element is split or not.
    };

    /// Destructor
    DivideGeometry::~DivideGeometry() {};

    /// Turn back information as a string.
    std::string DivideGeometry::Info() const {
        return "Base class for geometries splitting operations.";
    };

    /// Print information about this object.
    void DivideGeometry::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations.";
    };

    /// Print object's data.
    void DivideGeometry::PrintData(std::ostream& rOStream) const {
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

    DivideGeometry::GeometryType DivideGeometry::GetInputGeometry() const {
        return mrInputGeometry;
    };

    Vector DivideGeometry::GetNodalDistances() const {
        return mrNodalDistances;
    };

    void DivideGeometry::IsSplit() {
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
    };

    std::vector<DivideGeometry::IndexedPointGeometryPointerType> DivideGeometry::GetPositiveSubdivisions() const
    {
        return mPositiveSubdivisions;
    }

    std::vector<DivideGeometry::IndexedPointGeometryPointerType> DivideGeometry::GetNegativeSubdivisions() const
    {
        return mNegativeSubdivisions;
    }

    std::vector<DivideGeometry::IndexedPointGeometryPointerType> DivideGeometry::GetPositiveInterfaces() const
    {
        return mPositiveInterfaces;
    }

    std::vector<DivideGeometry::IndexedPointGeometryPointerType> DivideGeometry::GetNegativeInterfaces() const
    {
        return mNegativeInterfaces;
    }

    std::vector<unsigned int> DivideGeometry::GetPositiveInterfacesParentIds() const
    {
        return mPositiveInterfacesParentIds;
    }

    std::vector<unsigned int> DivideGeometry::GetNegativeInterfacesParentIds() const
    {
        return mNegativeInterfacesParentIds;
    }
};
