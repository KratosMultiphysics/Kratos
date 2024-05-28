// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef ELEMENT_CONTAINER_INCLUDE_H
#define ELEMENT_CONTAINER_INCLUDE_H

//// STL includes
#include <cstring>
#include <sstream>
#include <stdexcept>
//// Project includes
#include "containers/element.hpp"
#include "includes/parameters.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  ElementContainer
 * @author Manuel Messmer
 * @brief  Stores elements in vector and provides fast access via Id map.
 * @note Only active elements/knot spans are stored.
 * @todo Refactor. Store elements as unique_ptr
*/
template<typename TElementType>
class ElementContainer {

public:
    ///@name Type Defintitions
    ///@{
    typedef TElementType ElementType;
    typedef typename ElementType::IntegrationPointType IntegrationPointType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;

    typedef Unique<ElementType> ElementPtrType;
    typedef std::vector<ElementPtrType> ElementVectorPtrType;
    typedef std::vector<IntegrationPointType> IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef std::unordered_map<IndexType, IndexType> ElementIdMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ElementContainer(const Parameters& rParameters){
        mNumberOfElements = rParameters.NumberOfElements();
        mLastElementId = 0;
    }

    // Delete copy constructor
    ElementContainer(ElementContainer const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    const ElementType& GetElement(std::size_t id) const{
        return *pGetElement(id);
    }

    const ElementType* pGetElement(std::size_t id) const {
        auto found_key = mElementIdMap.find(id);
        if( found_key == mElementIdMap.end() )
            QuESo_ERROR << "ID does not exist.\n";
        return mElements[found_key->second].get();
    }

    ElementType* pGetElement(std::size_t id, bool& found){
        auto found_key = mElementIdMap.find(id);
        found = false;
        if( found_key != mElementIdMap.end() ){
            found = true;
            return mElements[found_key->second].get();
        }
        return nullptr;
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> begin() {
        return dereference_iterator(mElements.begin());
    }


    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> begin() const {
        return dereference_iterator(mElements.begin());
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> end() {
        return dereference_iterator(mElements.end());
    }

    DereferenceIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> end() const {
        return dereference_iterator(mElements.end());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> begin_to_ptr() {
        return raw_pointer_iterator(mElements.begin());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> begin_to_ptr() const {
        return raw_pointer_iterator(mElements.begin());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::iterator> end_to_ptr() {
        return raw_pointer_iterator(mElements.end());
    }

    RawPointerIterator<typename std::vector<std::unique_ptr<ElementType>>::const_iterator> end_to_ptr() const {
        return raw_pointer_iterator(mElements.end());
    }

    const ElementVectorPtrType& GetElements() const{
        return mElements;
    }

    void AddElement(ElementPtrType& rElement){
        const int current_id = rElement->GetId();
        auto found_key = mElementIdMap.find(current_id);
        if( found_key == mElementIdMap.end() ){
            // critical section
            if( rElement->GetId() > static_cast<IndexType>(mLastElementId) ){
                mLastElementId = rElement->GetId();
            }
            mElementIdMap.insert(std::pair<IndexType, IndexType>(rElement->GetId(), mElements.size()));
            mElements.push_back(std::move(rElement));
        }
        else {
            QuESo_ERROR << "ID already exists.\n";
        }
    }

    std::size_t size() const {
        return mElements.size();
    }

    void reserve(std::size_t new_capacity){
        mElements.reserve(new_capacity);
        mElementIdMap.reserve(new_capacity);
    }

    ElementType* pGetNextElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = id + 1;
        int next_index = id + 1;
        auto indices = GetMatrixIndicesFromVectorIndex(id);
        if( indices[0] == mNumberOfElements[0]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);
        if( found == false){  // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementType* pGetNextElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        // Make sure current element exists
        // TODO:: if id >= mLastElement error
        local_end = false;
        int next_index = GetNextIndexY(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[1] == mNumberOfElements[1]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetNextElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;

        int next_index = GetNextIndexZ(id, local_end);
        next_id = next_index;
        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[2] == mNumberOfElements[2]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetPreviousElementInX(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        local_end = false;
        next_id = id - 1;
        int next_index = id - 1;
        auto indices = GetMatrixIndicesFromVectorIndex(id);
        if( indices[0] == mNumberOfElements[0]-1) {
            local_end = true;
        }

        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }
        return found_element;
    }

    ElementType* pGetPreviousElementInY(std::size_t id, std::size_t& next_id, bool& found, bool& local_end) {
        // Make sure current element exists
        // TODO:: if id >= mLastElement error
        local_end = false;
        int next_index = GetPreviousIndexY(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[1] == mNumberOfElements[1]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    ElementType* pGetPreviousElementInZ(std::size_t id, std::size_t& next_id, bool& found, bool& local_end){
        local_end = false;

        int next_index = GetPreviousIndexZ(id, local_end);
        next_id = next_index;

        auto indices = GetMatrixIndicesFromVectorIndex(next_id-1); // Matrix starts with 0. Here Id's start with 1.
        if( indices[2] == mNumberOfElements[2]-1) {
            local_end = true;
        }
        auto found_element = pGetElement(next_index, found);

        if( found == false){                            // Element is not found
            local_end = true;
        }

        return found_element;
    }

    bool IsLast(std::size_t id, std::size_t direction ){
        auto indices = GetMatrixIndicesFromVectorIndex(id-1);

        switch( direction )
        {
        case 0: // Forward x
            return (indices[0] == (mNumberOfElements[0]-1));
        case 1: // Backward X
            return (indices[0] == 0);
        case 2: // Forward Y
            return (indices[1] == (mNumberOfElements[1]-1));
        case 3: // Backward Y
            return (indices[1] == 0);
        case 4: // Forward Z
            return (indices[2] == (mNumberOfElements[2]-1));
        case 5: // Backward Z
            return (indices[2] == 0);
        default:
            QuESo_ERROR << "There are only 6 different directions! \n";
        }
    }

    const IntegrationPointVectorPtrType pGetPoints(const char* type) const {
        IntegrationPointVectorPtrType points = MakeUnique<IntegrationPointVectorType>();
        const auto begin_el_itr_ptr = this->begin();
        for( IndexType i = 0; i < this->size(); ++i){
            const auto& el_ptr = *(begin_el_itr_ptr + i);
            IntegrationPointVectorType points_tmp;
            if( std::strcmp(type,"Trimmed") == 0 ){
                if( el_ptr->IsTrimmed() )
                    points_tmp = el_ptr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"Inside") == 0 ){
                if( !el_ptr->IsTrimmed() )
                    points_tmp = el_ptr->GetIntegrationPoints();
            }
            else if( std::strcmp(type,"All") == 0 ){
                points_tmp = el_ptr->GetIntegrationPoints();
            }
            else {
                QuESo_ERROR << "Given type '" << type << "' not available.\n";
            }
            points->insert(points->end(), points_tmp.begin(), points_tmp.end());
        }
        return points;
    }

    double GetVolumeOfAllIPs(){
        double volume = 0.0;
        const auto el_it_ptr_begin = this->begin();
        #pragma omp parallel for reduction(+ : volume)
        for( int i = 0; i < static_cast<int>(this->size()); ++i ){
            const auto& el_ptr = (*(el_it_ptr_begin + i));
            const double det_j = el_ptr->DetJ();
            const auto& r_points = el_ptr->GetIntegrationPoints();
            for( const auto& r_point : r_points ){
                volume += r_point.Weight()*det_j;
            }
        }

        return volume;
    }

private:

    IndexType GetNextIndexX(IndexType i){
        return i + 1;
    }

    IndexType GetNextIndexY(IndexType i, bool& local_end) const {
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] < mNumberOfElements[1]-1) {
            indices[1] += 1;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[1] = 0;
        }
        else {
            indices[2] += 1;
            indices[1] = 0;
            indices[0] = 0;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetNextIndexZ(IndexType i, bool& local_end) const {
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] < mNumberOfElements[2]-1) {
            indices[2] += 1;
        }
        else if( indices[0] < mNumberOfElements[0]-1){
            indices[0] += 1;
            indices[2] = 0;
        }
        else {
            indices[1] += 1;
            indices[2] = 0;
            indices[0] = 0;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexX(IndexType i){
        return i - 1;
    }

    IndexType GetPreviousIndexY(IndexType i, bool& local_end) const {
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[1] > 0) {
            indices[1] -= 1;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[1] = mNumberOfElements[1]-1;
        }
        else {
            indices[2] -= 1;
            indices[1] = mNumberOfElements[1]-1;
            indices[0] = mNumberOfElements[0]-1;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    IndexType GetPreviousIndexZ(IndexType i, bool& local_end) const {
        auto indices = GetMatrixIndicesFromVectorIndex(i-1);
        if( indices[2] > 0 ){
            indices[2] -= 1;
        }
        else if( indices[0] > 0){
            indices[0] -= 1;
            indices[2] = mNumberOfElements[2]-1;;
        }
        else {
            indices[1] -= 1;
            indices[2] = mNumberOfElements[2]-1;;
            indices[0] = mNumberOfElements[0]-1;;
        }

        IndexType new_index = GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);

        return new_index+1;
    }

    inline std::array<IndexType,3> GetMatrixIndicesFromVectorIndex(
        const IndexType Index) const noexcept
    {
        std::array<IndexType,3> result;
        const IndexType index_in_row_column_plane = Index % (mNumberOfElements[0]*mNumberOfElements[1]);
        result[0] = index_in_row_column_plane % mNumberOfElements[0]; // row
        result[1] = index_in_row_column_plane / mNumberOfElements[0]; // column
        result[2] = Index / (mNumberOfElements[0]*mNumberOfElements[1]);   // depth

        return result;
    }

    inline IndexType GetVectorIndexFromMatrixIndices(
        const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex) const noexcept
    {
        return DepthIndex * (mNumberOfElements[1]*mNumberOfElements[0]) + ColumnIndex * mNumberOfElements[0] + RowIndex;
    }

    ///@}
    ///@name Private member variables
    ///@{

    int mLastElementId;
    ElementVectorPtrType mElements{};
    ElementIdMapType mElementIdMap{};
    Vector3i mNumberOfElements{};
    ///@}
}; // End class Element container
///@} // End QuESo classes

} // End namespace queso
#endif // ELEMENT_CONTAINER_INCLUDE_H