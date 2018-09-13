//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: miguel $
//   Date:                $Date: 2014-03-27 14:21:15 $
//   Revision:            $Revision: 1.35 $
//
//


#pragma once


// System includes
#include <string>
#include <iostream> 
#include <assert.h>

//temporary

// External includes 

#define KRATOS_INDEPENDENT

// Project includes
#ifdef KRATOS_INDEPENDENT
#else
#include "includes/define.h"
#endif

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
        //KRATOS_CLASS_POINTER_DEFINITION(Octree_Pooyan);

        typedef TCellType cell_type;

        typedef typename cell_type::key_type key_type;

        typedef typename cell_type::configuration_type configuration_type;

        typedef double coordinate_type;

        enum {
            CHILDREN_NUMBER = cell_type::CHILDREN_NUMBER,
            DIMENSION = cell_type::DIMENSION,
            MAX_LEVEL = cell_type::MAX_LEVEL,
            ROOT_LEVEL = cell_type::ROOT_LEVEL,
            MIN_LEVEL = cell_type::MIN_LEVEL
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        OctreeBinary() : root_(new cell_type), number_of_cells_(CHILDREN_NUMBER + 1), number_of_leaves_(1), levels_(0) {
        }

        /// Destructor.

        virtual ~OctreeBinary() {
            delete root_;
        }

        double CalcSizeNormalized(const cell_type* cell) const {
            const double scale = 1.00 / (1 << ROOT_LEVEL);

            return (1 << cell->GetLevel()) * scale; // I'm doing it in this way to avoid division
        }

        double CalcMinCellNormalizedSize() const{          
          const double scale = 1.00 / (1 << ROOT_LEVEL);    
          return (1 << MIN_LEVEL) * scale; // I'm doing it in this way to avoid division
        }
        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

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
                        //                                        for (std::size_t j = 0; j < CHILDREN_NUMBER; j++) {
                        //                                            next_leaves.push_back(temp_neighbour_cell->GetChildren() + j);
                        //                                        }
                        //                                        temp_neighbour_cell = temp_neighbour_cell->pGetChild(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
                        std::size_t child_index = temp_neighbour_cell->GetChildIndex(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
                        for (std::size_t j = 0; j < CHILDREN_NUMBER; j++) {
                          //                                            *((temp_neighbour_cell->GetChildren() + j)->pGetDataPointer()) = TCellType::configuration_type::AllocateData();
                          //                                            *((temp_neighbour_cell->GetChildren() + j)->pGetData()) = 1;
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
            //                KRATOS_WATCH(leaves.size())
          }
#ifdef KRATOS_INDEPENDENT
#else
          KRATOS_WATCH(number_of_leaves_);
#endif
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

                                        //                                        for (std::size_t j = 0; j < CHILDREN_NUMBER; j++) {
                                        //                                            next_leaves.push_back(temp_neighbour_cell->GetChildren() + j);
                                        //                                        }
                                        //                                        temp_neighbour_cell = temp_neighbour_cell->pGetChild(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
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
#ifdef KRATOS_INDEPENDENT
#else
                KRATOS_WATCH(leaves.size())
#endif
            }
#ifdef KRATOS_INDEPENDENT
#else
            KRATOS_WATCH(number_of_leaves_);
#endif
        }

        bool BoundingBoxesCollide(const double min_cornerA[3],const double max_cornerA[3],
          const double min_cornerB[3],const double max_cornerB[3],const double tolerance) const {
            for (int idim=0;idim<3;idim++){
              if(min_cornerA[idim]>max_cornerB[idim]+tolerance) return false;
              if(max_cornerA[idim]<min_cornerB[idim]-tolerance) return false;
            }    
            return true;
        }

        void GetLeavesInBoundingBoxNormalized(const double* coord1, const double* coord2,
          std::vector<cell_type*>& leaves) const
        {        
	  // const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL) ; // 0.1% of the min size

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
			cells_stack.reserve(200);
            cells_stack.push_back(range_cell);
            while (!cells_stack.empty()) {
              cell_type* cell = cells_stack.back();
              cells_stack.pop_back();
              if (cell->HasChildren()) {
                for (std::size_t i = 0; i < CHILDREN_NUMBER; i++){
                  //abel. to be optimized
                  cell_type* child=cell->pGetChild(i);                 
                  double min_coord[3];                     
                  child->GetMinPointNormalized(min_coord);                               
                  double max_coord[3];   
                  child->GetMaxPointNormalized(max_coord);                 
                  const double tol=0.;
                  if (BoundingBoxesCollide(coord1,coord2,min_coord,max_coord,tol)){                    
                    cells_stack.push_back(cell->pGetChild(i));
                  }
                }
              } else
                leaves.push_back(cell);
            }


            return;
        }

        int GetAllLeavesVector(std::vector<cell_type*>& all_leaves) const {
          std::vector<cell_type*> cells_stack;
          cells_stack.push_back(root_);
          while (!cells_stack.empty()) {
            cell_type* cell = cells_stack.back();
            cells_stack.pop_back();
            if (cell->HasChildren()) {
              for (std::size_t i = 0; i < CHILDREN_NUMBER; i++){
                cells_stack.push_back(cell->pGetChild(i));
              }
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
#ifdef _DEBUG
          //abel
          cell_type* ret=NULL;
          switch (direction){
          case 0: ret = pGetLeftCell(p_cell);
            break;
          case 1: ret = pGetRightCell(p_cell);
            break;
          case 2:  ret = pGetBackCell(p_cell);
            break;
          case 3:  ret = pGetFrontCell(p_cell);
            break;
          case 4:  ret = pGetBottomCell(p_cell);
            break;
          case 5: ret = pGetTopCell(p_cell);
            break;
          default: assert(0); break;
          }
          //return ret;
          if (p_cell->GetNeighbourKey(direction, keys)) {
            assert(ret==pGetCell(keys));
          }
          else assert(!ret); // no neighbour
#endif
            if (p_cell->GetNeighbourKey(direction, keys)) {
                return pGetCell(keys);
            }
            return NULL; // no neighbour
        }

        cell_type * pGetNeighbourCell(cell_type* p_cell, std::size_t position, std::size_t direction) {
            //                 KRATOS_WATCH(position);
            //                KRATOS_WATCH(direction);
            //                KRATOS_WATCH(int(p_cell->GetLevel()));
            key_type keys[3];
            //                p_cell->GetMinKey(keys[0], keys[1], keys[2]);
            //                std::cout << "min_keys: " << keys[0] << ", " << keys[1] << ", " << keys[2] << std::endl;
            //                p_cell->GetKey(position, keys);
            //                std::cout << "pos_keys: " << keys[0] << ", " << keys[1] << ", " << keys[2] << std::endl;
            if (p_cell->GetNeighbourKey(position, direction, keys)) {

                //                std::cout << "new_keys: " << keys[0] << ", " << keys[1] << ", " << keys[2] << std::endl;
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
               
        void InsertNormalized(typename cell_type::pointer_type object)
        {
            const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL) ; // 0.1% of the min size
            //const double tolerance = EPSILON4;
            const double min_coord_object[3]={0.00, 0.00, 0.00};
            const double max_coord_object[3]={0.00, 0.00, 0.00};

            configuration_type::GetBoundingBox(object, min_coord_object,  max_coord_object);

            //BoundingBox bbox_to_check;                     
            //bbox_to_check.Join(min_coord_object);
            //bbox_to_check.Join(max_coord_object);

            key_type min_x_key = CalcKeyNormalized(min_coord_object[0]);
            key_type min_y_key = CalcKeyNormalized(min_coord_object[1]);
            key_type min_z_key = CalcKeyNormalized(min_coord_object[2]);

            key_type max_x_key = CalcKeyNormalized(max_coord_object[0]);
            key_type max_y_key = CalcKeyNormalized(max_coord_object[1]);
            key_type max_z_key = CalcKeyNormalized(max_coord_object[2]);

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
#ifdef _DEBUG
            {
              //abel
              double cell_min_point[3];
              double cell_max_point[3];

              range_cell->GetMinPointNormalized(cell_min_point);
              range_cell->GetMaxPointNormalized(cell_max_point);
              for(int i = 0 ; i < 3 ; i++)
              {
                assert(cell_min_point[i] <= min_coord_object[i]);
                assert(cell_max_point[i] >= max_coord_object[i]);
              }
              //assert(configuration_type::IsIntersected(object,tolerance, cell_min_point, cell_max_point);
            }
#endif

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
                  double min_coord_cell[3];                     
                  child->GetMinPointNormalized(min_coord_cell);                                
                  double max_coord_cell[3];
                  child->GetMaxPointNormalized(max_coord_cell);                  
                  if (BoundingBoxesCollide(min_coord_object,max_coord_object,
                    min_coord_cell,max_coord_cell,tolerance)){                     
                      cells_stack.push_back(child);
                  }
                }
              } else{
                // we are in a leaf and we can check if it intersects with the object                  
                double cell_min_point[3];
                double cell_max_point[3];

                cell->GetMinPointNormalized(cell_min_point);
                cell->GetMaxPointNormalized(cell_max_point);

                //POOYAN. NEED TO SCALE BACK TO ORIGINAL COORDINATE IN NO NORMALIZED FUNCTION
               
                const int is_intersected = configuration_type::IsIntersected(object,tolerance, cell_min_point, cell_max_point);
                if(is_intersected)
                  cell->Insert(object);

              }
            }


//            std::cout << "min_coord_object : [" << min_coord_object[0] << "," << min_coord_object[1] << "," << min_coord_object[2] << std::endl;
//            std::cout << "max_coord_object : [" << max_coord_object[0] << "," << max_coord_object[1] << "," << max_coord_object[2] << std::endl;



        }



        #ifdef KRATOS_INDEPENDENT

#else


        void GetIntersectedLeaves(typename cell_type::pointer_type object, std::vector<cell_type*>& leaves){

            const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL) ; // 0.1% of the min size

            double min_coord[3]={0.00, 0.00, 0.00};
            double max_coord[3]={0.00, 0.00, 0.00};



            //
            // to be added to configure

            configuration_type::GetBoundingBox(object, min_coord,  max_coord);
            //KRATOS_ERROR(std::logic_error,"To be added to configure", "")

            key_type min_x_key = CalcKey(min_coord[0]);
            key_type min_y_key = CalcKey(min_coord[1]);
            key_type min_z_key = CalcKey(min_coord[2]);

            key_type max_x_key = CalcKey(max_coord[0]);
            key_type max_y_key = CalcKey(max_coord[1]);
            key_type max_z_key = CalcKey(max_coord[2]);

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
#ifdef _DEBUG
            {
              //abel
              double cell_min_point[3];
              double cell_max_point[3];

              range_cell->GetMinPoint(cell_min_point);
              range_cell->GetMaxPoint(cell_max_point);
              for(int i = 0 ; i < 3 ; i++)
              {
                assert(cell_min_point[i] <= min_coord[i]);
                assert(cell_max_point[i] >= max_coord[i]);
              }
              //assert(configuration_type::IsIntersected(object,tolerance, cell_min_point, cell_max_point);
            }
#endif

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
                  child->GetMinPoint(low);
                  child->GetMaxPoint(high);
//                  KRATOS_WATCH_3(low);
//                  KRATOS_WATCH_3(high);
//                  KRATOS_WATCH_3(min_coord);
//                  KRATOS_WATCH_3(max_coord);
//                  KRATOS_WATCH(Collides(min_coord, max_coord, low, high));
                  if (Collides(min_coord, max_coord, low, high))
                    cells_stack.push_back(child);
                }
              } else{
                // we are in a leaf and we can check if it intersects with the object
                double cell_min_point[3];
                double cell_max_point[3];

                cell->GetMinPoint(cell_min_point);
                cell->GetMaxPoint(cell_max_point);

//                KRATOS_WATCH(object->GetGeometry());
//                KRATOS_WATCH(tolerance);
//                std::cout << "cell_min_point : " << cell_min_point[0] << "," << cell_min_point[1] << "," << cell_min_point[2] << std::endl;
//                std::cout << "cell_max_point : " << cell_max_point[0] << "," << cell_max_point[1] << "," << cell_max_point[2] << std::endl;

                const int is_intersected = /*configuration_type::*/IsIntersected(object,tolerance, cell_min_point, cell_max_point);

                if(is_intersected)
                  leaves.push_back(cell);

              }
            }


//            std::cout << "min_coord : [" << min_coord[0] << "," << min_coord[1] << "," << min_coord[2] << std::endl;
//            std::cout << "max_coord : [" << max_coord[0] << "," << max_coord[1] << "," << max_coord[2] << std::endl;



        }






  //////////////////////////////////////////////////////////




       inline bool  IsIntersected(const Element::Pointer rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)
    {
        Point low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
        Point high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);


        return HasIntersection(rObject->GetGeometry(), low_point, high_point);
    }




    /// detect if  triangle and box are intersected
    virtual bool HasIntersection(Element::GeometryType& geom_1, const Point& rLowPoint, const Point& rHighPoint )
    {
//        const BaseType& geom_1 = rGeometry;

        Point boxcenter;
        Point boxhalfsize;

        boxcenter[0]   = 0.50 * ( rLowPoint[0] + rHighPoint[0] );
        boxcenter[1]   = 0.50 * ( rLowPoint[1] + rHighPoint[1] );
        boxcenter[2]   = 0.50 * ( rLowPoint[2] + rHighPoint[2] );

        boxhalfsize[0] = 0.50 * ( rHighPoint[0] - rLowPoint[0] );

        boxhalfsize[1] = 0.50 * ( rHighPoint[1] - rLowPoint[1] );
        boxhalfsize[2] = 0.50 * ( rHighPoint[2] - rLowPoint[2] );

        std::size_t size = geom_1.size();

        std::vector<Point > triverts;

        triverts.resize( size );

        for ( unsigned int i = 0; i < size; i++ )
        {
            triverts[i] =  geom_1.GetPoint( i );
        }

//        std::cout << "triangle: " << triverts[0] << ", " << triverts[1] << ", " << triverts[2];

        return TriBoxOverlap( boxcenter, boxhalfsize, triverts );

    }

        inline bool TriBoxOverlap( Point& boxcenter, Point& boxhalfsize, std::vector< Point >& triverts )
    {

        /*    use separating axis theorem to test overlap between triangle and box */
        /*    need to test for overlap in these directions: */
        /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
        /*       we do not even need to test these) */
        /*    2) normal of the triangle */
        /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
        /*       this gives 3x3=9 more tests */

        double min, max, d, p0, p1, p2, rad, fex, fey, fez;
        array_1d<double, 3 > v0, v1, v2;
        array_1d<double, 3 > axis;
        array_1d<double, 3 > normal, e0, e1 , e2;
//
//     /* This is the fastest branch on Sun */
//     /* move everything so that the boxcenter is in (0,0,0) */
        noalias( v0 ) = triverts[0] - boxcenter;
        noalias( v1 ) = triverts[1] - boxcenter;
        noalias( v2 ) = triverts[2] - boxcenter;
//
//     /* compute triangle edges */
        noalias( e0 ) = v1 - v0;    /* tri edge 0 */
        noalias( e1 ) = v2 - v1;    /* tri edge 1 */
        noalias( e2 ) = v0 - v2;    /* tri edge 2 */
//
//     /* Bullet 3:  */
//     /*  test the 9 tests first (this was faster) */
        fex = fabs( e0[0] );
        fey = fabs( e0[1] );
        fez = fabs( e0[2] );
        //AXISTEST_X01(e0[2], e0[1], fez, fey);

        if ( AxisTest_X01( e0[2], e0[1], fez, fey, p0, p2, min, max, rad, v0, v2, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Y02(e0[2], e0[0], fez, fex);
        if ( AxisTest_Y02( e0[2], e0[0], fez, fex, p0, p2, min, max, rad, v0, v2, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Z12(e0[1], e0[0], fey, fex);
        if ( AxisTest_Z12( e0[1], e0[0], fey, fex, p1, p2, min, max, rad, v1, v2, boxhalfsize ) == 0 ) return false;



        fex = fabs( e1[0] );

        fey = fabs( e1[1] );

        fez = fabs( e1[2] );

        //AXISTEST_X01(e1[2], e1[1], fez, fey);
        if ( AxisTest_X01( e1[2], e1[1], fez, fey, p0, p2, min, max, rad, v0, v2, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Y02(e1[2], e1[0], fez, fex);
        if ( AxisTest_Y02( e1[2], e1[0], fez, fex, p0, p2, min, max, rad, v0, v2, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Z0(e1[1], e1[0], fey, fex);
        if ( AxisTest_Z0( e1[1], e1[0], fey, fex, p0,  p1, min, max, rad, v0, v1, boxhalfsize ) == 0 ) return false;


        fex = fabs( e2[0] );

        fey = fabs( e2[1] );

        fez = fabs( e2[2] );

        //AXISTEST_X2(e2[2], e2[1], fez, fey);
        if ( AxisTest_X2( e2[2], e2[1], fez, fey, p0, p1, min, max, rad, v0, v1, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Y1(e2[2], e2[0], fez, fex);
        if ( AxisTest_Y1( e2[2], e2[0], fez, fex, p0, p1, min, max, rad, v0, v1, boxhalfsize ) == 0 ) return false;

        //AXISTEST_Z12(e2[1], e2[0], fey, fex);
        if ( AxisTest_Z12( e2[1], e2[0], fey, fex, p1, p2, min, max, rad, v1, v2, boxhalfsize ) == 0 ) return false;


        /* Bullet 1: */
        /*  first test overlap in the {x,y,z}-directions */
        /*  find min, max of the triangle each direction, and test for overlap in */
        /*  that direction -- this is equivalent to testing a minimal AABB around */
        /*  the triangle against the AABB */

        /* test in X-direction */
        FindMinMax( v0[0], v1[0], v2[0], min, max );

        if ( min > boxhalfsize[0] || max < -boxhalfsize[0] ) return false;

        /* test in Y-direction */
        FindMinMax( v0[1], v1[1], v2[1], min, max );

        if ( min > boxhalfsize[1] || max < -boxhalfsize[1] ) return false;

        /* test in Z-direction */
        FindMinMax( v0[2], v1[2], v2[2], min, max );

        if ( min > boxhalfsize[2] || max < -boxhalfsize[2] ) return false;

        /* Bullet 2: */
        /*  test if the box intersects the plane of the triangle */
        /*  compute plane equation of triangle: normal*x+d=0 */
        MathUtils<double>::CrossProduct( normal, e0, e1 );

        d = -inner_prod( normal, v0 );  /* plane eq: normal.x+d=0 */

        if ( !planeBoxOverlap( normal, d, boxhalfsize ) ) return false;

        return true;   /* box and triangle overlaps */
    }

   inline void FindMinMax( const double& x0,
                            const double& x1,
                            const double& x2,
                            double& min,
                            double& max )
    {
        min = max = x0;

        if ( x1 < min ) min = x1;

        if ( x1 > max ) max = x1;

        if ( x2 < min ) min = x2;

        if ( x2 > max ) max = x2;
    }


    inline  bool planeBoxOverlap( const array_1d<double, 3 >& normal,  const double& d, const array_1d<double, 3 >& maxbox )
    {
        int q;
        array_1d<double, 3 >  vmin, vmax;

        for ( q = 0; q <= 2; q++ )
        {
            if ( normal[q] > 0.0f )
            {
                vmin[q] = -maxbox[q];
                vmax[q] = maxbox[q];
            }
            else
            {
                vmin[q] = maxbox[q];
                vmax[q] = -maxbox[q];
            }
        }

        if ( inner_prod( normal, vmin ) + d > 0.0f ) return false;

        if ( inner_prod( normal, vmax ) + d >= 0.0f ) return true;

        return false;
    }

    /*======================== X-tests ========================*/
    inline unsigned int  AxisTest_X01( double& a,   double& b,
                                       double& fa,  double& fb,
                                       double& p0,  double& p2,
                                       double& min, double& max, double& rad,
                                       array_1d<double, 3 >& v0,
                                       array_1d<double, 3 >& v2,
                                       Point& boxhalfsize
                                     )
    {
        p0 = a * v0[1] - b * v0[2];
        p2 = a * v2[1] - b * v2[2];

        if ( p0 < p2 )
        {
            min = p0;
            max = p2;
        }
        else
        {
            min = p2;
            max = p0;
        }

        rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }

    inline unsigned int  AxisTest_X2( double& a,   double& b,
                                      double& fa,  double& fb,
                                      double& p0,  double& p1,
                                      double& min, double& max, double& rad,
                                      array_1d<double, 3 >& v0,
                                      array_1d<double, 3 >& v1,
                                      Point& boxhalfsize
                                    )
    {
        p0 = a * v0[1] - b * v0[2];
        p1 = a * v1[1] - b * v1[2];

        if ( p0 < p1 )
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }

        rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }

    /*======================== Y-tests ========================*/
    inline unsigned int  AxisTest_Y02( double& a, double& b,
                                       double& fa,  double& fb,
                                       double& p0,  double& p2,
                                       double& min, double& max, double& rad,
                                       array_1d<double, 3 >& v0,
                                       array_1d<double, 3 >& v2,
                                       Point& boxhalfsize
                                     )
    {

        p0 = -a * v0[0] + b * v0[2];
        p2 = -a * v2[0] + b * v2[2];

        if ( p0 < p2 )
        {
            min = p0;
            max = p2;
        }
        else
        {
            min = p2;
            max = p0;
        }

        rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }

    inline unsigned int  AxisTest_Y1( double& a,   double& b,
                                      double& fa,  double& fb,
                                      double& p0,  double& p1,
                                      double& min, double& max, double& rad,
                                      array_1d<double, 3 >& v0,
                                      array_1d<double, 3 >& v1,
                                      Point& boxhalfsize
                                    )

    {
        p0 = -a * v0[0] + b * v0[2];
        p1 = -a * v1[0] + b * v1[2];

        if ( p0 < p1 )
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }

        rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }

    /*======================== Z-tests ========================*/

    inline unsigned int  AxisTest_Z12( double& a, double& b,
                                       double& fa,  double& fb,
                                       double& p1,  double& p2,
                                       double& min, double& max, double& rad,
                                       array_1d<double, 3 >& v1,
                                       array_1d<double, 3 >& v2,
                                       Point& boxhalfsize
                                     )
    {
        p1 = a * v1[0] - b * v1[1];
        p2 = a * v2[0] - b * v2[1];

        if ( p2 < p1 )
        {
            min = p2;
            max = p1;
        }
        else
        {
            min = p1;
            max = p2;
        }

        rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }

    inline unsigned int  AxisTest_Z0( double& a, double& b,
                                      double& fa,  double& fb,
                                      double& p0,  double& p1,
                                      double& min, double& max, double& rad,
                                      array_1d<double, 3 >& v0,
                                      array_1d<double, 3 >& v1,
                                      Point& boxhalfsize
                                    )
    {
        p0 = a * v0[0] - b * v0[1];
        p1 = a * v1[0] - b * v1[1];

        if ( p0 < p1 )
        {
            min = p0;
            max = p1;
        }
        else
        {
            min = p1;
            max = p0;
        }

        rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];

        if ( min > rad || max < -rad ) return 0;
        else return 1;
    }















///////////////////////////////////////////



#endif // KRATOS_INDEPENDENT


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
                double cell_size = cell->GetSize();

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

          return (static_cast<double>(key) * scale);
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
                    std::size_t node_id = (*(leaf->pGetData()))[i_point];
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
                    rOStream << (*(leaf->pGetData()))[i_point] << "  ";

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

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
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

        OctreeBinary & operator=(OctreeBinary const& rOther) {
            return *this;
        }

        /// Copy constructor.

        OctreeBinary(OctreeBinary const& rOther) {

        }


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




