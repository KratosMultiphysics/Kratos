
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Abel Coll
//                   Pooyan Dadvand
//

#if !defined(KRATOS_OCTREE_H_INCLUDED )
#define  KRATOS_OCTREE_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "includes/node.h"

#include "octree_binary_cell.h"

#define KRATOS_WATCH_3(name) std::cout << #name << " : " << name[0] << ", " << name[1] << ", " << name[2] << std::endl;

namespace Kratos {
    ///@addtogroup ApplicationNameApplication
    ///@{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.

    /** Detail class definition.
     */
    template <class TCellType>
    class OctreeBinary {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Octree
        //KRATOS_CLASS_POINTER_DEFINITION(OctreeBinary);

        typedef TCellType cell_type;

        typedef typename cell_type::key_type key_type;

        typedef typename cell_type::configuration_type configuration_type;

        typedef double coordinate_type;

        typedef Node<3> NodeType;

        typedef Geometry<NodeType> GeometryType;

        enum {
            CHILDREN_NUMBER = cell_type::CHILDREN_NUMBER,
            DIMENSION = cell_type::DIMENSION,
            MAX_LEVEL = cell_type::MAX_LEVEL,
            ROOT_LEVEL = cell_type::ROOT_LEVEL,
            MIN_LEVEL = cell_type::MIN_LEVEL // must be greater iqual to 2
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        OctreeBinary() : root_(new cell_type), number_of_cells_(CHILDREN_NUMBER + 1), number_of_leaves_(1), levels_(0) {

            for(int i = 0 ; i < DIMENSION ; i++)
            {
                mScaleFactor[i] = 1.00;
                mOffset[i] = 0.00;
            }
        }

        OctreeBinary(const double*  NewScaleFactor, const double* NewOffset) : root_(new cell_type), number_of_cells_(CHILDREN_NUMBER + 1), number_of_leaves_(1), levels_(0) {
            for(int i = 0 ; i < DIMENSION ; i++)
            {
                mScaleFactor[i] = NewScaleFactor[i];
                mOffset[i] = NewOffset[i];
            }
        }

        /// Destructor.

        virtual ~OctreeBinary() {
            delete root_;
        }

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void SetBoundingBox(const coordinate_type * Low, const coordinate_type * High)
        {
            for(int i = 0 ; i < DIMENSION ; i++)
            {
                mScaleFactor[i] = 1/(High[i] - Low[i]);
                mOffset[i] = -Low[i];
            }
        }

        double CalcSizeNormalized(const cell_type* cell) const {
            const double scale = 1.00 / (1 << ROOT_LEVEL);

            return (1 << cell->GetLevel()) * scale; // I'm doing it in this way to avoid division
        }

        double CalcMinCellNormalizedSize() const{
          const double scale = 1.00 / (1 << ROOT_LEVEL);
          return (1 << MIN_LEVEL) * scale; // I'm doing it in this way to avoid division
        }

        void NormalizeCoordinates(coordinate_type* Coordinates) const
        {
            for(int i = 0 ; i < DIMENSION ; i++)
            {
                Coordinates[i] += mOffset[i];
                Coordinates[i] *= mScaleFactor[i];
            }
        }

        void NormalizeCoordinates(const coordinate_type * Coordinates, coordinate_type * NormalizedCoordinates) const
        {
            for(int i = 0 ; i < 3 ; i++)
            {
                NormalizedCoordinates[i] = Coordinates[i] + mOffset[i];
                NormalizedCoordinates[i] *= mScaleFactor[i];
            }
        }

        void CalculateCoordinateNormalized(const key_type key, coordinate_type& NormalizedCoordinate) const {
            const double scale = 1.00 / (1 << ROOT_LEVEL);
            NormalizedCoordinate = static_cast<double>(key * scale);
        }

        void CalculateCoordinatesNormalized(const key_type* keys, coordinate_type * NormalizedCoordinates) const {
            const double scale = 1.00 / (1 << ROOT_LEVEL);
            for(int i = 0 ; i < DIMENSION ; i++)
                NormalizedCoordinates[i] = static_cast<double>(keys[i] * scale);
        }

        void CalculateCoordinates(key_type* keys, coordinate_type * ResultCoordinates) const {
            const double scale = 1.00 / (1 << ROOT_LEVEL);
            for(int i = 0 ; i < DIMENSION ; i++)
                ResultCoordinates[i] = (static_cast<double>(keys[i] * scale) / mScaleFactor[i]) - mOffset[i];
        }

        void ScaleBackToOriginalCoordinate(coordinate_type * ThisCoordinates) const
        {
            for(int i = 0 ; i < DIMENSION ; i++)
            {
                ThisCoordinates[i] /=  mScaleFactor[i];
                ThisCoordinates[i] -= mOffset[i];
            }

        }

        void ScaleBackToOriginalCoordinate(const coordinate_type * ThisCoordinates, coordinate_type * ResultCoordinates) const
        {
            for(int i = 0 ; i < DIMENSION ; i++)
            {
                ResultCoordinates[i] = ThisCoordinates[i] /  mScaleFactor[i];
                ResultCoordinates[i] -= mOffset[i];
            }

        }

        //pooyan. uncomment this when needed
        //key_type CalcKey(coordinate_type coordinate) const {
        //  //pooyan. NEED TO SCALE COORDINATE!!!
        //  return CalcKeyNormalized(coordinate);
        //}

        key_type CalcKeyNormalized(coordinate_type coordinate) const {
          assert(coordinate>=0.); assert(coordinate<=1.);
          return static_cast<key_type> ((1 << ROOT_LEVEL) * coordinate);
        }

        void InsertNormalized(coordinate_type* point) {
            key_type x_key = CalcKeyNormalized(point[0]);
            key_type y_key = CalcKeyNormalized(point[1]);
            key_type z_key = CalcKeyNormalized(point[2]);

            cell_type* cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > MIN_LEVEL; i--) {
                if (cell->IsLeaf()) {
                    SubdivideCell(cell);
                }
                cell = cell->pGetChild(x_key, y_key, z_key);
            }

        }

        void Insert(coordinate_type* point) {
            coordinate_type normalized_point[3];
            NormalizeCoordinates(point, normalized_point);
            key_type x_key = CalcKeyNormalized(normalized_point[0]);
            key_type y_key = CalcKeyNormalized(normalized_point[1]);
            key_type z_key = CalcKeyNormalized(normalized_point[2]);

            cell_type* cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > MIN_LEVEL; i--) {
                if (cell->IsLeaf()) {
                    SubdivideCell(cell);
                }
                cell = cell->pGetChild(x_key, y_key, z_key);
            }

        }

        bool CheckConstrain2To1() const{
          //POOYAN. This function must return true if the octree is balanced (constrained2:1) and false if not.
          return true;
        }

        void Constrain2To1() {
          std::vector<cell_type*> leaves;
          std::vector<cell_type*> next_leaves;

          //when the function will be at upper level (in mesher instead of octree) this vector (leaves) should be passed and copied instead of recomputed
          GetAllLeavesVector(leaves);

          for (char i_level = MIN_LEVEL; i_level < ROOT_LEVEL - 1; i_level++) {
            for (std::size_t i_cell = 0; i_cell < leaves.size(); i_cell++) {
              cell_type* current_leaf = leaves[i_cell];
              if (current_leaf->GetLevel() == i_level) {
                key_type neighbour_key[3];
                //18 is the number of neighbours counting faces and edges of the cell
                for (int i_direction = 0; i_direction < 18; i_direction++) {
                  if (current_leaf->GetNeighbourKey(i_direction, neighbour_key)) {
                    cell_type* neighbour_cell = pGetCell(neighbour_key);
                    if (neighbour_cell->GetLevel() > i_level + 1) {
                      cell_type* temp_neighbour_cell = neighbour_cell;
                      for (char j_level = neighbour_cell->GetLevel(); j_level > i_level + 1; j_level--) {
                        SubdivideCell(temp_neighbour_cell);
                        temp_neighbour_cell->TransferObjectsToSonsNormalized();
                        std::size_t child_index = temp_neighbour_cell->GetChildIndex(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
                        for (std::size_t j = 0; j < CHILDREN_NUMBER; j++) {
                          if (j != child_index) {
                            next_leaves.push_back(temp_neighbour_cell->GetChildren() + j);
                          }
                        }
                        temp_neighbour_cell = temp_neighbour_cell->GetChildren() + child_index;
                        if (j_level == neighbour_cell->GetLevel() - 1) // the last loop we add all the child as leaf
                          next_leaves.push_back(temp_neighbour_cell);
                      }
                    }
                  }
                }
              } else if (current_leaf->IsLeaf()) { // becuase it may be divided
                next_leaves.push_back(current_leaf);
              }
            }
            leaves.swap(next_leaves);
            next_leaves.clear();
          }
        }

        void Constrain2To1New() {
            std::vector<cell_type*> leaves;
            std::vector<cell_type*> next_leaves;

            GetAllLeavesVector(leaves);

            for (char i_level = MIN_LEVEL; i_level < ROOT_LEVEL - 1; i_level++) {
                for (int i_direction = 0; i_direction < 18; i_direction++) {
                    for (std::size_t i_cell = 0; i_cell < leaves.size(); i_cell++) {
                        cell_type* current_leaf = leaves[i_cell];
                        if (current_leaf->GetLevel() == i_level) {
                            key_type neighbour_key[3];
                            if (current_leaf->GetNeighbourKey(i_direction, neighbour_key)) {
                                cell_type* neighbour_cell = pGetCell(neighbour_key);
                                if (neighbour_cell->GetLevel() > i_level + 1) {
                                    cell_type* temp_neighbour_cell = neighbour_cell;
                                    for (char j_level = neighbour_cell->GetLevel(); j_level > i_level + 1; j_level--) {
                                        SubdivideCell(temp_neighbour_cell);

                                        std::size_t child_index = temp_neighbour_cell->GetChildIndex(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
                                        for (std::size_t j = 0; j < CHILDREN_NUMBER; j++) {
                                            if (j != child_index) {
                                                next_leaves.push_back(temp_neighbour_cell->GetChildren() + j);
                                            }
                                        }
                                        temp_neighbour_cell = temp_neighbour_cell->GetChildren() + child_index;
                                        if (j_level == neighbour_cell->GetLevel() - 1) // the last loop we add all the child as leaf
                                            next_leaves.push_back(temp_neighbour_cell);
                                    }
                                }
                            }
                        } else if (i_direction == 0) {
                            if (current_leaf->IsLeaf()) { // becuase it may be divided
                                next_leaves.push_back(current_leaf);
                            }
                        }
                    }
                }
                leaves.swap(next_leaves);
                next_leaves.clear();
                KRATOS_WATCH(leaves.size())
            }
            KRATOS_WATCH(number_of_leaves_);
        }

        void GetLeavesInBoundingBoxNormalized(const double* coord1, const double* coord2,
          std::vector<cell_type*>& leaves) const
        {
            key_type min_x_key = CalcKeyNormalized(coord1[0]);
            key_type min_y_key = CalcKeyNormalized(coord1[1]);
            key_type min_z_key = CalcKeyNormalized(coord1[2]);

            key_type max_x_key = CalcKeyNormalized(coord2[0]);
            key_type max_y_key = CalcKeyNormalized(coord2[1]);
            key_type max_z_key = CalcKeyNormalized(coord2[2]);

            key_type delta_x = min_x_key^max_x_key;
            key_type delta_y = min_y_key^max_y_key;
            key_type delta_z = min_z_key^max_z_key;

            // finding the level of the cell containing the entire region
            std::size_t min_level_1 = ROOT_LEVEL;
            std::size_t min_level_2 = ROOT_LEVEL;
            std::size_t min_level = ROOT_LEVEL;

            const std::size_t one = 1;
            while (!(delta_x & (one << min_level_1)) && (min_level_1 > MIN_LEVEL)) min_level_1--;
            while (!(delta_y & (one << min_level_2)) && (min_level_2 > min_level_1)) min_level_2--;
            while (!(delta_z & (one << min_level)) && (min_level > min_level_2)) min_level--;
            min_level++;

            cell_type* range_cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > min_level; i--) {
              if (range_cell->IsLeaf()) {
                break;
              }
              range_cell = range_cell->pGetChild(min_x_key, min_y_key, min_z_key);

            }
            // Now we have the cell (or leaf) containing the entire range and from now on we have to gather the leaves
            std::vector<cell_type*> cells_stack;
            cells_stack.push_back(range_cell);
            while (!cells_stack.empty()) {
              cell_type* cell = cells_stack.back();
              cells_stack.pop_back();
              if (cell->HasChildren()) {
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++){
                  //abel. to be optimized
                  cell_type* child=cell->pGetChild(i);

                  double low[3];
                  double high[3];
                  child->GetMinPoint(low);
                  child->GetMaxPoint(high);
                  if (Collides(coord1, coord2, low, high))
                    cells_stack.push_back(cell->pGetChild(i));
                }
              } else
                leaves.push_back(cell);
            }


            return;
        }


        bool Collides(const double* Low1, const double* High1, const double* Low2, const double* High2)
        {
            return (Low1[0] <= High2[0]) &&
                   (Low1[1] <= High2[1]) &&
                   (Low1[2] <= High2[2]) &&
                   (Low2[0] <= High1[0]) &&
                   (Low2[1] <= High1[1]) &&
                   (Low2[2] <= High1[2]);
        }


        int GetAllLeavesVector(std::vector<cell_type*>& all_leaves) const {
          std::vector<cell_type*> cells_stack;
          cells_stack.push_back(root_);
          while (!cells_stack.empty()) {
            cell_type* cell = cells_stack.back();
            cells_stack.pop_back();
            if (cell->HasChildren()) {
              for (std::size_t i = 0; i < CHILDREN_NUMBER; i++)
                        cells_stack.push_back(cell->pGetChild(i));
                } else
                    all_leaves.push_back(cell);
            }

            return 0;
        }

        cell_type * pGetCellNormalized(const coordinate_type * point) const {
            key_type keys[3];
            keys[0] = CalcKeyNormalized(point[0]);
            keys[1] = CalcKeyNormalized(point[1]);
            keys[2] = CalcKeyNormalized(point[2]);

            return pGetCell(keys);
        }

        cell_type * pGetCell(key_type * keys) const {
            cell_type* cell = root_;

            for (std::size_t i = 0; i < ROOT_LEVEL; i++) {
                if (cell->IsLeaf()) {
                    return cell;
                }
                cell = cell->pGetChild(keys[0], keys[1], keys[2]);
            }
            return cell;
        }

        cell_type * pGetCell(key_type* keys, std::size_t level) const {
            cell_type* cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > level; i--) {
                if (cell->IsLeaf()) {
                    return cell;
                }
                cell = cell->pGetChild(keys[0], keys[1], keys[2]);
            }
            return cell;
        }

        cell_type * pGetLeftCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetLeftKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetLeftCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetLeftKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetRightCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetRightKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetRightCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetRightKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetBackCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetBackKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetBackCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetBackKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetFrontCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetFrontKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetFrontCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetFrontKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetTopCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetTopKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetTopCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetTopKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetBottomCell(const cell_type * p_cell) {
          key_type keys[3];
          if (p_cell->GetBottomKey(keys)) {
            return pGetCell(keys);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetBottomCell(cell_type* p_cell, std::size_t level) {
          key_type keys[3];
          if (p_cell->GetBottomKey(keys)) {
            return pGetCell(keys, level);
          }
          return NULL; // no neighbour
        }

        cell_type * pGetNeighbourCell(const cell_type* p_cell, std::size_t direction) {
          key_type keys[3];

            if (p_cell->GetNeighbourKey(direction, keys)) {
                return pGetCell(keys);
            }
            return NULL; // no neighbour
        }

        cell_type * pGetNeighbourCell(cell_type* p_cell, std::size_t position, std::size_t direction) {
            key_type keys[3];

            if (p_cell->GetNeighbourKey(position, direction, keys)) {
                return pGetCell(keys);
            }
            return NULL; // no neighbour
        }

        int SubdivideCell(cell_type* p_cell) {
          number_of_cells_ += CHILDREN_NUMBER;
          number_of_leaves_ += CHILDREN_NUMBER - 1;

          return p_cell->SubdivideCell();
        }


        int SubvidiveUntilSizeNormalized(double* coord, const double desired_size){
            key_type x_key = CalcKeyNormalized(coord[0]);
            key_type y_key = CalcKeyNormalized(coord[1]);
            key_type z_key = CalcKeyNormalized(coord[2]);

            cell_type* cell = root_;


            // I'm assuming that the desired size is also normalized to the [0,1] span
            std::size_t scaled_size = std::size_t(desired_size * (1 << ROOT_LEVEL));

            if(scaled_size < (1 << MIN_LEVEL))
                scaled_size = (1 << MIN_LEVEL);

            for (std::size_t i = ROOT_LEVEL; (std::size_t(1) << i) > scaled_size; i--) {
                if (cell->IsLeaf()) {
                    SubdivideCell(cell);
                }
                cell = cell->pGetChild(x_key, y_key, z_key);
            }

          return 0;
        }

        int RefineWithUniformSizeNormalized(const double uniform_size){
            const double min_size = double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL);
            double cell_size = uniform_size;
            if(cell_size < min_size)
                cell_size = min_size;

            std::vector<cell_type*> cells_stack;
            cells_stack.push_back(root_);

            while (!cells_stack.empty()) {
              cell_type* cell = cells_stack.back();
              cells_stack.pop_back();
              if(CalcSizeNormalized(cell) > cell_size){
                if (cell->IsLeaf()) {
                  SubdivideCell(cell);
                }
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++)
                  cells_stack.push_back(cell->pGetChild(i));
              }
            }


          return 0;
        }

        void InsertNormalized(typename cell_type::pointer_type object){

            const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL) ; // 0.1% of the min size

            double min_coord[3]={0.00, 0.00, 0.00};
            double max_coord[3]={0.00, 0.00, 0.00};



            //
            // to be added to configure

            configuration_type::GetBoundingBox(object, min_coord,  max_coord);
            //KRATOS_THROW_ERROR(std::logic_error,"To be added to configure", "")

            key_type min_x_key = CalcKeyNormalized(min_coord[0]);
            key_type min_y_key = CalcKeyNormalized(min_coord[1]);
            key_type min_z_key = CalcKeyNormalized(min_coord[2]);

            key_type max_x_key = CalcKeyNormalized(max_coord[0]);
            key_type max_y_key = CalcKeyNormalized(max_coord[1]);
            key_type max_z_key = CalcKeyNormalized(max_coord[2]);

            key_type delta_x = min_x_key^max_x_key;
            key_type delta_y = min_y_key^max_y_key;
            key_type delta_z = min_z_key^max_z_key;

            // finding the level of the cell containing the entire region
            std::size_t min_level_1 = ROOT_LEVEL;
            std::size_t min_level_2 = ROOT_LEVEL;
            std::size_t min_level = ROOT_LEVEL;

            const std::size_t one = 1;
            while (!(delta_x & (one << min_level_1)) && (min_level_1 > MIN_LEVEL)) min_level_1--;
            while (!(delta_y & (one << min_level_2)) && (min_level_2 > min_level_1)) min_level_2--;
            while (!(delta_z & (one << min_level)) && (min_level > min_level_2)) min_level--;
            min_level++;

            cell_type* range_cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > min_level ; i--) {
                if (range_cell->IsLeaf()) {
                    //SubdivideCell(range_cell);
                  break;
                }
                range_cell = range_cell->pGetChild(min_x_key, min_y_key, min_z_key);
            }

            // Now we have the cell (or leaf) containing the entire range and from now on we have to intersect the object with all childs
            std::vector<cell_type*> cells_stack;
            cells_stack.push_back(range_cell);
            while (!cells_stack.empty()) {
              cell_type* cell = cells_stack.back();
              cells_stack.pop_back();
              if (cell->HasChildren()) {
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++){
                  //abel. to be optimized
                  cell_type* child=cell->pGetChild(i);
                  double low[3];
                  double high[3];
                  child->GetMinPointNormalized(low);
                  child->GetMaxPointNormalized(high);
                  if (Collides(min_coord, max_coord, low, high))
                    cells_stack.push_back(child);
                }
              } else{
                // we are in a leaf and we can check if it intersects with the object
                double cell_min_point[3];
                double cell_max_point[3];

                cell->GetMinPointNormalized(cell_min_point);
                cell->GetMaxPointNormalized(cell_max_point);

                // I have to put this in no normailzed part. Pooyan.
//                ScaleBackToOriginalCoordinate(cell_min_point);
//                ScaleBackToOriginalCoordinate(cell_max_point);

                const int is_intersected = /*configuration_type::*/IsIntersected(object,tolerance, cell_min_point, cell_max_point);
                if(is_intersected)
                  cell->Insert(object);

              }
            }
        }

        void Insert(typename cell_type::pointer_type object){

            const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL) ; // 0.1% of the min size

            double min_coord[3]={0.00, 0.00, 0.00};
            double max_coord[3]={0.00, 0.00, 0.00};

            //
            // to be added to configure

            configuration_type::GetBoundingBox(object, min_coord,  max_coord);
            //KRATOS_THROW_ERROR(std::logic_error,"To be added to configure", "")
            NormalizeCoordinates(min_coord);
            NormalizeCoordinates(max_coord);

            key_type min_x_key = CalcKeyNormalized(min_coord[0]);
            key_type min_y_key = CalcKeyNormalized(min_coord[1]);
            key_type min_z_key = CalcKeyNormalized(min_coord[2]);

            key_type max_x_key = CalcKeyNormalized(max_coord[0]);
            key_type max_y_key = CalcKeyNormalized(max_coord[1]);
            key_type max_z_key = CalcKeyNormalized(max_coord[2]);

            key_type delta_x = min_x_key^max_x_key;
            key_type delta_y = min_y_key^max_y_key;
            key_type delta_z = min_z_key^max_z_key;

            // finding the level of the cell containing the entire region
            std::size_t min_level_1 = ROOT_LEVEL;
            std::size_t min_level_2 = ROOT_LEVEL;
            std::size_t min_level = ROOT_LEVEL;

            const std::size_t one = 1;
            while (!(delta_x & (one << min_level_1)) && (min_level_1 > MIN_LEVEL)) min_level_1--;
            while (!(delta_y & (one << min_level_2)) && (min_level_2 > min_level_1)) min_level_2--;
            while (!(delta_z & (one << min_level)) && (min_level > min_level_2)) min_level--;
            min_level++;

            cell_type* range_cell = root_;

            for (std::size_t i = ROOT_LEVEL; i > min_level ; i--) {
                if (range_cell->IsLeaf()) {
                    //SubdivideCell(range_cell);
                  break;
                }
                range_cell = range_cell->pGetChild(min_x_key, min_y_key, min_z_key);

            }

            // Now we have the cell (or leaf) containing the entire range and from now on we have to intersect the object with all childs
            std::vector<cell_type*> cells_stack;
            cells_stack.push_back(range_cell);
            while (!cells_stack.empty()) {
              cell_type* cell = cells_stack.back();
              cells_stack.pop_back();
              if (cell->HasChildren()) {
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++){
                  //abel. to be optimized
                  cell_type* child=cell->pGetChild(i);
                  double low[3];
                  double high[3];
                  child->GetMinPointNormalized(low);
                  child->GetMaxPointNormalized(high);
                  if (Collides(min_coord, max_coord, low, high))
                    cells_stack.push_back(child);
                }
              } else{
                // we are in a leaf and we can check if it intersects with the object
                double cell_min_point[3];
                double cell_max_point[3];

                cell->GetMinPointNormalized(cell_min_point);
                cell->GetMaxPointNormalized(cell_max_point);

                // I have to put this in no normailzed part. Pooyan.
                ScaleBackToOriginalCoordinate(cell_min_point);
                ScaleBackToOriginalCoordinate(cell_max_point);

                const int is_intersected = /*configuration_type::*/IsIntersected(object,tolerance, cell_min_point, cell_max_point);
                if(is_intersected)
                  cell->Insert(object);

              }
            }


//            std::cout << "min_coord : [" << min_coord[0] << "," << min_coord[1] << "," << min_coord[2] << std::endl;
//            std::cout << "max_coord : [" << max_coord[0] << "," << max_coord[1] << "," << max_coord[2] << std::endl;



        }

    /**
    * @brief This method fills the intersected leaves
    * @param pObject Pointer to the object of interest containing the geometry to check the intersection
    * @param rLeaves Leaves to be filled for search
    * @param ToleranceCoefficient The proportion of the minimal size to compute the tolerance (0.1% of the min size by default)
    */
    void GetIntersectedLeaves(
        typename cell_type::pointer_type pObject,
        std::vector<cell_type*>& rLeaves,
        const double ToleranceCoefficient = 0.001  // 0.1% of the min size
        )
    {
        const double tolerance = ToleranceCoefficient * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL);

        double min_coord[3]={0.00, 0.00, 0.00};
        double max_coord[3]={0.00, 0.00, 0.00};

        // TODO: To be added to configure

        configuration_type::GetBoundingBox(pObject, min_coord,  max_coord);
        NormalizeCoordinates(min_coord);
        NormalizeCoordinates(max_coord);

        key_type min_x_key = CalcKeyNormalized(min_coord[0]);
        key_type min_y_key = CalcKeyNormalized(min_coord[1]);
        key_type min_z_key = CalcKeyNormalized(min_coord[2]);

        key_type max_x_key = CalcKeyNormalized(max_coord[0]);
        key_type max_y_key = CalcKeyNormalized(max_coord[1]);
        key_type max_z_key = CalcKeyNormalized(max_coord[2]);

        key_type delta_x = min_x_key^max_x_key;
        key_type delta_y = min_y_key^max_y_key;
        key_type delta_z = min_z_key^max_z_key;

        // Finding the level of the cell containing the entire region
        std::size_t min_level_1 = ROOT_LEVEL;
        std::size_t min_level_2 = ROOT_LEVEL;
        std::size_t min_level = ROOT_LEVEL;

        const std::size_t one = 1;
        while (!(delta_x & (one << min_level_1)) && (min_level_1 > MIN_LEVEL)) min_level_1--;
        while (!(delta_y & (one << min_level_2)) && (min_level_2 > min_level_1)) min_level_2--;
        while (!(delta_z & (one << min_level)) && (min_level > min_level_2)) min_level--;
        min_level++;

        cell_type* range_cell = root_;

        for (std::size_t i = ROOT_LEVEL; i > min_level ; i--) {
            if (range_cell->IsLeaf()) {
                //SubdivideCell(range_cell);
                break;
            }
            range_cell = range_cell->pGetChild(min_x_key, min_y_key, min_z_key);
        }

        // Now we have the cell (or leaf) containing the entire range and from now on we have to intersect the pObject with all childs
        std::vector<cell_type*> cells_stack;
        cells_stack.push_back(range_cell);
        while (!cells_stack.empty()) {
            cell_type* cell = cells_stack.back();
            cells_stack.pop_back();
            if (cell->HasChildren()) {
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++) {
                    // TODO: Abel. to be optimized
                    cell_type* child=cell->pGetChild(i);
                    double low[3];
                    double high[3];
                    child->GetMinPointNormalized(low);
                    child->GetMaxPointNormalized(high);
                    if (Collides(min_coord, max_coord, low, high)) {
                        cells_stack.push_back(child);
                    }
                }
            } else {
                // we are in a leaf and we can check if it intersects with the pObject
                double cell_min_point[3];
                double cell_max_point[3];

                cell->GetMinPointNormalized(cell_min_point);
                cell->GetMaxPointNormalized(cell_max_point);

                // I have to put this in no normailzed part. Pooyan.
                ScaleBackToOriginalCoordinate(cell_min_point);
                ScaleBackToOriginalCoordinate(cell_max_point);

                const bool is_intersected = IsIntersected(pObject,tolerance, cell_min_point, cell_max_point);

                if(is_intersected) {
                    rLeaves.push_back(cell);
                }
            }
        }
    }


    inline bool  IsIntersected(typename cell_type::pointer_type rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)
    {
        Point low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
        Point high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);


        return HasIntersection(rObject->GetGeometry(), low_point, high_point);
    }

    /**
     * @brief Detect if  triangle and box are intersected
     */
    virtual bool HasIntersection(
        GeometryType& rGeometry,
        const Point& rLowPoint,
        const Point& rHighPoint
        )
    {
        return rGeometry.HasIntersection(rLowPoint, rHighPoint);
    }

        cell_type* pGetCellContainRegion(key_type min_x_key, key_type min_y_key, key_type min_z_key,
                                         key_type max_x_key, key_type max_y_key, key_type max_z_key)
        {
            cell_type* cell = root_;

        }


        ///@}
        ///@name Access
        ///@{

        cell_type * pGetRoot() {
            return root_;
        }


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        void PrintGiDMesh(std::ostream & rOStream) const {
            std::vector<cell_type*> leaves;

            GetAllLeavesVector(leaves);

            std::cout << "writing " << leaves.size() << " leaves" << std::endl;
            rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
            rOStream << "# color 96 96 96" << std::endl;
            rOStream << "Coordinates" << std::endl;
            rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

            std::size_t node_index = 1;
            for (std::size_t i = 0; i < leaves.size(); i++) {
                cell_type* cell = leaves[i];
                double min_point[3];
                cell->GetMinPoint(min_point);

                double cell_size = cell->CalcSize();

                for (std::size_t j = 0; j < 2; j++)
                    for (std::size_t k = 0; k < 2; k++)
                        for (std::size_t h = 0; h < 2; h++) {
                            rOStream << node_index++ << "  " << min_point[0] + j * cell_size << "  " << min_point[1] + k * cell_size << "  " << min_point[2] + h * cell_size << std::endl;
                        }
            }

            rOStream << "end coordinates" << std::endl;
            rOStream << "Elements" << std::endl;
            rOStream << "# element node_1 node_2 node_3 material_number" << std::endl;

            for (std::size_t i = 0; i < leaves.size(); i++) {
                if ((leaves[i]->pGetData()))
                    rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << leaves[i]->GetLevel() + 100 << std::endl;
                else
                    rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << int(leaves[i]->GetLevel()) << std::endl;

            }
            rOStream << "end elements" << std::endl;

        }

        double GetCoordinateNormalized(key_type key) const {
          const double scale = 1.00 / (1 << ROOT_LEVEL);

          return static_cast<double>(key * scale);
        }

        void PrintGiDMeshNew(std::ostream & rOStream) const {
            std::vector<cell_type*> leaves;

            GetAllLeavesVector(leaves);

            std::cout << "writing " << leaves.size() << " leaves" << std::endl;
            rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
            rOStream << "# color 96 96 96" << std::endl;
            rOStream << "Coordinates" << std::endl;
            rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;
            std::size_t node_number = 0;
            //            std::size_t node_index = 1;
            for (std::size_t i = 0; i < leaves.size(); i++) {
                cell_type* leaf = leaves[i];
                for (std::size_t i_point = 0; i_point < 8; i_point++) {
                    std::size_t node_id = (*(leaf->pGetData()))[i_point]->Id();
                    if (node_id > node_number) {
                        key_type point_key[3];
                        leaf->GetKey(i_point, point_key);
                        double point_coordinate[3];

                        for (std::size_t j = 0; j < DIMENSION; j++) {
                            point_coordinate[j] = leaf->GetCoordinate(point_key[j]);
                        }
                        rOStream << node_id << "  " << point_coordinate[0] << "  " << point_coordinate[1] << "  " << point_coordinate[2] << std::endl;
                        node_number++;
                    }
                }
            }

            rOStream << "end coordinates" << std::endl;
            rOStream << "Elements" << std::endl;
            rOStream << "# element node_1 node_2 node_3 material_number" << std::endl;

            for (std::size_t i = 0; i < leaves.size(); i++) {
                cell_type* leaf = leaves[i];
                rOStream << i + 1 << "  ";
                for (std::size_t i_point = 0; i_point < 8; i_point++)
                    rOStream << (*(leaf->pGetData()))[i_point]->Id() << "  ";

                rOStream << std::endl;
            }
            rOStream << "end elements" << std::endl;

        }

        /// Turn back information as a string.

        virtual std::string Info() const {
            return "Octree";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream & rOStream) const {
            rOStream << Info();
        }

        /// Print object's data.

        virtual void PrintData(std::ostream & rOStream) const {
            rOStream << "Number of cells  : " << number_of_cells_ << std::endl;
            rOStream << "Number of leaves : " << number_of_leaves_ << std::endl;
            //rOStream << *root_;
        }


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

        cell_type* root_;

        std::size_t number_of_cells_;
        std::size_t number_of_leaves_;
        std::size_t levels_;

        coordinate_type mOffset[3];
        coordinate_type mScaleFactor[3];


        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.

        //OctreeBinary & operator=(OctreeBinary const& rOther) {
        //   return *this;
        //}

        /// Copy constructor.

        //OctreeBinary(OctreeBinary const& rOther) {

        //}


        ///@}

    }; // Class Octree

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    //    inline std::istream & operator >>(std::istream& rIStream,
    //            Octree& rThis);

    /// output stream function

    template <class TCellType>
    inline std::ostream & operator <<(std::ostream& rOStream,
    const OctreeBinary<TCellType>& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_OCTREE_H_INCLUDED  defined


