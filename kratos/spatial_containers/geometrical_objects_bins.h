//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


#if !defined(KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <unordered_set>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/bounding_box.h"
#include "geometries/point.h"
#include "includes/geometrical_object.h"
#include "includes/global_pointer.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{


template <typename TObjectType>
class SpatialSearchResult {
	using TPointerType = GlobalPointer<TObjectType>;
	TPointerType mpObject;
	double mDistance2;
	bool mIsObjectFound;
	bool mIsDistanceCalculated;

public:
	SpatialSearchResult() : mpObject(nullptr), mDistance2(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {}
	SpatialSearchResult(TObjectType* pObject) : mpObject(pObject), mDistance2(0.00), mIsObjectFound(false), mIsDistanceCalculated(false) {
		if (mpObject.get() != nullptr)
			mIsObjectFound = true;
	}

	SpatialSearchResult(SpatialSearchResult const& /* Other */) = default;

	SpatialSearchResult(SpatialSearchResult&& /* Other */) = default;

	TPointerType Get() { return mpObject; }
	TPointerType const Get() const { return mpObject; }
	void Set(TObjectType* pObject) {
		mpObject = pObject;
		mIsObjectFound = true;
	}
	bool IsObjectFound() const { return mIsObjectFound; }

	double GetDistance2() const { return mDistance2; }
	void SetDistance2(double TheDistance2) {
		mDistance2 = TheDistance2;
		mIsDistanceCalculated = true;
	}
	bool IsDistanceCalculated() const { return mIsDistanceCalculated; }

	void Reset() {
		mpObject = mpObject(nullptr);
		mDistance2 = 0.00;
		mIsObjectFound = false;
		mIsDistanceCalculated = false;
	}

        SpatialSearchResult& operator=(SpatialSearchResult const& /*Other*/) = default;

};

/// A bins container for 3 dimensional geometries.
/** Detail class definition.
*/
class GeometricalObjectsBins
{
    static constexpr int Dimension = 3;
    static constexpr double Tolerance = 1e-12;
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObjectsBins
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObjectsBins);

    using CellType = std::vector<GeometricalObject*>;
    using ResultType = SpatialSearchResult<GeometricalObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    GeometricalObjectsBins() = delete;

    template<typename TIteratorType>
	GeometricalObjectsBins(TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd) {
        std::size_t number_of_objects = std::distance(GeometricalObjectsBegin, GeometricalObjectsEnd);
        if(number_of_objects > 0){
            mBoundingBox.Set(GeometricalObjectsBegin->GetGeometry().begin(), GeometricalObjectsBegin->GetGeometry().end());
            for(TIteratorType i_object = GeometricalObjectsBegin ; i_object != GeometricalObjectsEnd ; i_object++){
                mBoundingBox.Extend(i_object->GetGeometry().begin() , i_object->GetGeometry().end());
            }
	    }
        mBoundingBox.Extend(Tolerance);
		CalculateCellSize(number_of_objects);
        mCells.resize(GetTotalNumberOfCells());
        AddObjectsToCells(GeometricalObjectsBegin, GeometricalObjectsEnd);
    }


    /// Destructor.
    virtual ~GeometricalObjectsBins(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    const BoundingBox<Point>& GetBoundingBox() const {
        return mBoundingBox;
    }

    const array_1d<double, 3>& GetCellSizes(){
        return mCellSizes;
    }

    const array_1d<std::size_t, 3>& GetNumberOfCells(){
        return mNumberOfCells;
    }

    const std::size_t GetTotalNumberOfCells(){
        return mNumberOfCells[0] * mNumberOfCells[1] * mNumberOfCells[2];
    }

    CellType& GetCell(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mNumberOfCells[0] + K * mNumberOfCells[1] * mNumberOfCells[0];
        return mCells[index];
    }

    BoundingBox<Point> GetCellBoundingBox(std::size_t I, std::size_t J, std::size_t K){
        BoundingBox<Point> result;
        
        result.GetMinPoint()[0] = mBoundingBox.GetMinPoint()[0] + I * mCellSizes[0];
        result.GetMinPoint()[1] = mBoundingBox.GetMinPoint()[1] + J * mCellSizes[1];
        result.GetMinPoint()[2] = mBoundingBox.GetMinPoint()[2] + K * mCellSizes[2];
         
        result.GetMaxPoint()[0] = mBoundingBox.GetMinPoint()[0] + (I + 1) * mCellSizes[0];
        result.GetMaxPoint()[1] = mBoundingBox.GetMinPoint()[1] + (J + 1) * mCellSizes[1];
        result.GetMaxPoint()[2] = mBoundingBox.GetMinPoint()[2] + (K + 1) * mCellSizes[2];

        return result;
    }

    template<typename TPointType>
    void SearchInRadius(TPointType const& ThePoint, double Radius, std::vector<ResultType>& rResults) {
        std::unordered_set<GeometricalObject*> results;

        array_1d< std::size_t, 3 > min_position;
        array_1d< std::size_t, 3 > max_position;

        for(int i = 0; i < 3; i++ ) {
            min_position[ i ] = CalculatePosition( ThePoint[i] - Radius, i );
            max_position[ i ] = CalculatePosition( ThePoint[i] + Radius, i ) + 1;
        }
        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    auto& cell = GetCell(i,j,k);
                    SearchInRadiusInCell(cell, ThePoint, Radius, results);
                }
            }
        }

        rResults.clear();
        for(auto p_object : results){
            rResults.push_back(ResultType(p_object));
        }
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GeometricalObjectsBins" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GeometricalObjectsBins";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
	
    BoundingBox<Point> mBoundingBox;
	array_1d<std::size_t, Dimension> mNumberOfCells;
    array_1d<double, 3>  mCellSizes;
	array_1d<double, 3>  mInverseOfCellSize;
    std::vector<CellType> mCells;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateCellSize(std::size_t NumberOfPoints) {
		std::size_t avarage_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(NumberOfPoints), 1.00 / Dimension));
		std::array<double, 3> lengths;
		double avarage_length = 0.00;
		for (int i = 0; i < Dimension; i++) {
			lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
			avarage_length += lengths[i];
		}
		avarage_length *= 1.00 / 3.00;

		if (avarage_length < std::numeric_limits<double>::epsilon()) {
			mNumberOfCells = ScalarVector(3, 1);
			return;
		}

		for (int i = 0; i < Dimension; i++) {
			mNumberOfCells[i] = static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1;
			if (mNumberOfCells[i] > 1)
				mCellSizes[i] = lengths[i] / mNumberOfCells[i];
			else
				mCellSizes[i] = avarage_length;

			mInverseOfCellSize[i] = 1.00 / mCellSizes[i];
		}

	}


    template<typename TIteratorType>
	void AddObjectsToCells(TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd) {
        for(auto i_geometrical_object = GeometricalObjectsBegin ; i_geometrical_object != GeometricalObjectsEnd ; i_geometrical_object++){
            array_1d< std::size_t, 3 > min_position;
            array_1d< std::size_t, 3 > max_position;
            CalculateMinMaxPositions(i_geometrical_object->GetGeometry(), min_position, max_position);
            for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                        auto cell_bounding_box = GetCellBoundingBox(i,j,k);
                        if(IsIntersected(i_geometrical_object->GetGeometry(), cell_bounding_box, Tolerance)){
                            GetCell(i,j,k).push_back(&(*i_geometrical_object));
                        }
                    }
                }
            }
        }        
    }


    template<typename TGeometryType>
    void CalculateMinMaxPositions(TGeometryType const& TheGeometry, array_1d< std::size_t, 3 >& MinPosition, array_1d< std::size_t, 3 >& MaxPosition){
        if(TheGeometry.empty())
            return;

        BoundingBox<Point> bounding_box(TheGeometry.begin(), TheGeometry.end());

        for(int i = 0; i < 3; i++ ) {
            MinPosition[ i ] = CalculatePosition( bounding_box.GetMinPoint()[i], i );
            MaxPosition[ i ] = CalculatePosition( bounding_box.GetMaxPoint()[i], i ) + 1;
        }
    }


    void CalculateMinMaxPositions(BoundingBox<Point> TheBox, array_1d< std::size_t, 3 >& MinPosition, array_1d< std::size_t, 3 >& MaxPosition){
        for(int i = 0; i < 3; i++ ) {
            MinPosition[ i ] = CalculatePosition( TheBox.GetMinPoint()[i], i );
            MaxPosition[ i ] = CalculatePosition( TheBox.GetMaxPoint()[i], i ) + 1;
        }
    }

    std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const {
        auto distance = Coordinate - mBoundingBox.GetMinPoint()[ ThisDimension ];
        distance = ( distance < 0.00 ) ? 0.00 : distance;
        std::size_t position =
            static_cast< std::size_t >( distance * mInverseOfCellSize[ ThisDimension ] );
        std::size_t result= ( position > mNumberOfCells[ ThisDimension ] - 1 )
                                ? mNumberOfCells[ ThisDimension ] - 1
                                : position;
        return result;
    }
    
    template<typename TGeometryType>
    static inline bool IsIntersected(TGeometryType& TheGeometry, BoundingBox<Point> Box, const double tolerance)
    {
        Point rLowPointTolerance;
        Point rHighPointTolerance;
        
        for(std::size_t i = 0; i<3; i++)
        {
            rLowPointTolerance[i]  =  Box.GetMinPoint()[i] - tolerance;
            rHighPointTolerance[i] =  Box.GetMaxPoint()[i] + tolerance;
        }
        
        return  TheGeometry.HasIntersection(rLowPointTolerance,rHighPointTolerance);
    }

    template<typename TPointType>
    void SearchInRadiusInCell(CellType const& TheCell, TPointType const& ThePoint, double Radius, std::unordered_set<GeometricalObject*>& rResults) {
        for(auto p_geometrical_object : TheCell){  
            auto& geometry = p_geometrical_object->GetGeometry();
            // TODO: Change this to new Distance method of the geometry to be more general
            double distance = GeometryUtils::PointDistanceToTriangle3D(
			geometry[0],
			geometry[1],
			geometry[2],
			ThePoint);
            if((Radius + Tolerance) > distance){
                rResults.insert(p_geometrical_object);
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator deleted.
    GeometricalObjectsBins& operator=(GeometricalObjectsBins const& rOther) = delete;

    /// Copy constructor deleted.
    GeometricalObjectsBins(GeometricalObjectsBins const& rOther) = delete;

    ///@}

}; // Class GeometricalObjectsBins

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                GeometricalObjectsBins& rThis){
                    return rIStream;
                }

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const GeometricalObjectsBins& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED  defined


