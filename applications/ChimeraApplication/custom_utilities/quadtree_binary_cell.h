//
//   Project Name:        Kratos
//   Last Modified by:    $Author: abel $
//   Date:                $Date: 2013-05-27 16:48:42 $
//   Revision:            $Revision: 1.23 $
//
//

#if !defined(KRATOS_QUADTREE_BINARY_CELL_H_INCLUDED)
#define KRATOS_QUADTREE_BINARY_CELL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#ifdef KRATOS_INDEPENDENT
#else
#include "includes/define.h"
#endif

namespace Kratos
{
///@addtogroup SpatialContainer
///@{

///@}
///@name Kratos Classes
///@{

/// This class represents a cell in an quadtree to be used with quadtree class

/** This class represents a cell in an quadtree and holds its level, min_key
     *  ,the children of the cell and pointer to a data class defined by
     *  configuration class.
     *  The level_ start from ROOT_LEVEL and each children has 1 level less than
     *  its parents but more than MIN_LEVEL.
     *  The keys are the binary bisection pattern in each dimension for the
     *  min point of the cell.
     *  The children are stored in in an array of 8 in following order
     *
     *                    
     *
     *   
     *                 Front
     *            ________________
     *           |        |       | 
     *           |    2   |   3   |
     *    Left   |_______ |_______| Right
     *           |        |       |
     *           |    0   |   1   |
     *           |________|_______|   
     *  y   
     *  |                Back  
     *  |_ __x
     *
    */
template <class TConfiguration>
class QuadtreeBinaryCell
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QuadtreeBinaryCell
    //KRATOS_CLASS_POINTER_DEFINITION(QuadtreeBinaryCell);

    typedef typename TConfiguration::data_type data_type;

    typedef TConfiguration configuration_type;

    typedef typename TConfiguration::pointer_type pointer_type;

    typedef std::vector<pointer_type> object_container_type;

    typedef std::size_t key_type;

    enum
    {
        CHILDREN_NUMBER = 4,
        DIMENSION = TConfiguration::DIMENSION,
        MAX_LEVEL = TConfiguration::MAX_LEVEL,
        ROOT_LEVEL = MAX_LEVEL - 1,
        MIN_LEVEL = TConfiguration::MIN_LEVEL
    };

    enum
    {
        LEFT = 0,
        RIGHT = 1,

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    QuadtreeBinaryCell(char Level = ROOT_LEVEL) : level_(Level), children_(NULL), data_(NULL)
    {
        for (std::size_t i = 0; i < DIMENSION; i++)
            min_key_[i] = 0;
    }

    /// Destructor.

    virtual ~QuadtreeBinaryCell()
    {

        if (data_)
            configuration_type::DeleteData(data_);
        delete[] children_;
    }

    void DeleteChildren()
    {
        delete[] children_;
        children_ = NULL;
    }
    void DeleteData()
    {
        if (data_)
        {
            configuration_type::DeleteData(data_);
            data_ = NULL;
        }
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::size_t GetChildIndex(key_type x_key, key_type y_key) const
    {
        char next_level = (char)(level_ - 1);
        key_type level_bit = 1 << next_level;
        return (((x_key & level_bit) >> next_level) + (((y_key & level_bit) >> next_level) << 1));
    }

    int SubdivideCell()
    {
        if (level_ == 0)
            return 1;
        if (children_)
            return 1;

        children_ = new QuadtreeBinaryCell[CHILDREN_NUMBER];

        char next_level = (char)(level_ - 1);

        for (std::size_t i = 0; i < CHILDREN_NUMBER; i++)
        {
            children_[i].SetMinKey(min_key_[0] | ((i & 1) << next_level), min_key_[1] | (((i & 2) >> 1)) << next_level);
            children_[i].SetLevel(next_level);
        }

        return 0; // Zero says no error!
    }

    void GetMinPointNormalized(double *min_point) const
    {
        for (std::size_t i = 0; i < DIMENSION; i++)
        {
            min_point[i] = GetCoordinateNormalized(min_key_[i]);
        }
    }

    void GetMaxPointNormalized(double *max_point) const
    {
        double size = CalcSizeNormalized();
        for (std::size_t i = 0; i < DIMENSION; i++)
        {
            max_point[i] = GetCoordinateNormalized(min_key_[i]) + size;
        }
    }

    int GetLeftKey(key_type *keys) const
    {

        if (min_key_[0] >= 1)
        {
            keys[0] = min_key_[0] - 1;
            keys[1] = min_key_[1];

            return 1; // There is a neighbour
        }
        return 0; // There is no neighbour
    }

    int GetRightKey(key_type *keys) const
    {
        if (min_key_[0] < static_cast<key_type>((1 << ROOT_LEVEL) - (1 << level_)))
        {
            keys[0] = min_key_[0] + (static_cast<key_type>(1) << level_);
            keys[1] = min_key_[1];

            return 1; // There is a neighbour
        }
        return 0; // There is no neighbour
    }

    int GetBackKey(key_type *keys) const
    {
        //
        if (min_key_[1] >= 1)
        {
            keys[0] = min_key_[0];
            keys[1] = min_key_[1] - 1;

            return 1; // There is a neighbour
        }
        return 0; // There is no neighbour
    }

    int GetFrontKey(key_type *keys) const
    {
        if (min_key_[1] < static_cast<key_type>((1 << ROOT_LEVEL) - (1 << level_)))
        {
            keys[0] = min_key_[0];
            keys[1] = min_key_[1] + (static_cast<key_type>(1) << level_);

            return 1; // There is a neighbour
        }
        return 0; // There is no neighbour
    }

    int GetKey(std::size_t position, key_type *keys) const
    {
        //in total, there are a maximum of 9 nodes. The only difference is that here
        //the central node of the quad (cell) is just after the lineal nodes.
        // the node index              :  1  2  3  4  5  6  7  8  9
        const std::size_t x_position[] = {0, 2, 2, 0, 1, 1, 2, 1, 0};
        const std::size_t y_position[] = {0, 0, 2, 2, 1, 0, 1, 2, 1};

        keys[0] = min_key_[0] + ((x_position[position]) << (level_ - 1));
        keys[1] = min_key_[1] + ((y_position[position]) << (level_ - 1));

        return 1;
    }

    int GetNeighbourKey(std::size_t direction, key_type *keys) const
    {

        //the key returned is inside the cell (minkey +1), to ensure that the corresponding
        //cell get in pGetCell is the right one.

        assert(direction < 4);
        const std::size_t x_offset[] = {0, 2, 1, 1};
        const std::size_t y_offset[] = {1, 1, 0, 2};

        const std::size_t x_coef[] = {0, 1, 0, 0};
        const std::size_t y_coef[] = {0, 0, 0, 1};

        std::size_t size = (1 << level_);

        // here i'm adding 2 to each dimension and it won't overflow
        keys[0] = min_key_[0] + x_offset[direction] + x_coef[direction] * size;
        keys[1] = min_key_[1] + y_offset[direction] + y_coef[direction] * size;

        for (int i = 0; i < DIMENSION; i++)
        {
            if (keys[i] == 0)
                return 0; // There is no neighbour
            else
                (keys[i])--;

            //now we are in correct location and we can check if the right nieigbours exist
            if (keys[i] > static_cast<key_type>(1 << ROOT_LEVEL))
                return 0; // There is no neighbour
        }
        return 1; // There is a neighbour
    }

    int GetNeighbourKey(std::size_t position, std::size_t direction, key_type *keys) const
    {

        GetKey(position, keys);

        // here i'm adding 2 to each dimension and it won't overflow
        keys[0] += (direction & 1) << 1;
        keys[1] += (direction & 2);

        for (int i = 0; i < DIMENSION; i++)
        {
            if (keys[i] == 0)
                return 0; // There is no neighbour
            else
                (keys[i])--;

            //now we are in correct location and we can check if the right nieigbours exist
            if (keys[i] > static_cast<key_type>(1 << ROOT_LEVEL))
                return 0; // There is no neighbour
        }
        return 1; // There is a neighbour
    }

    std::size_t GetLocalPosition(key_type *keys)
    {
        key_type position[2];
        //in total, there are a maximum of 9 nodes (as the Quadratic9 Hexa in GiD). The only difference is that here
        //the central node of the cell is just after the lineal nodes.
        // the node index              : 0  1  2  3  4  5  6  7  8
        const std::size_t local_index[] = {0, 5, 1, 8, 4, 6, 3, 7, 2};

        for (std::size_t i = 0; i < DIMENSION; i++)
        {
            position[i] = (keys[i] - min_key_[i]) >> (level_ - 1);
        }
        std::size_t index = position[0] + position[1] * 3;
        assert(index <= 8);
        //            if(local_index[index] > 26)
        //            {
        //                KRATOS_WATCH(int(level_));
        //                KRATOS_WATCH(index);
        //                KRATOS_WATCH(min_key_[0]);
        //                KRATOS_WATCH(min_key_[1]);
        //                KRATOS_WATCH(min_key_[2]);
        //                KRATOS_WATCH(keys[0]);
        //                KRATOS_WATCH(keys[1]);
        //                KRATOS_WATCH(keys[2]);
        //                KRATOS_WATCH(position[0]);
        //                KRATOS_WATCH(position[1]);
        //                KRATOS_WATCH(position[2]);
        //                KRATOS_WATCH(local_index[index]);
        //            }
        return local_index[index];
    }

    void Insert(pointer_type object)
    {
        objects_.push_back(object);
    }

    void TransferObjectsToSonsNormalized()
    {

        if (!objects_.size())
            return;
        assert(this->HasChildren());

        const double tolerance = 0.001 * double(1 << MIN_LEVEL) / double(1 << ROOT_LEVEL); // 0.1% of the min size
        double min_coord[3] = {0.00, 0.00, 0.00};
        double max_coord[3] = {0.00, 0.00, 0.00};

        for (std::size_t i = 0; i < CHILDREN_NUMBER; i++)
        {
            QuadtreeBinaryCell *son = pGetChild(i);
            if (son->HasChildren())
            {
                son->TransferObjectsToSonsNormalized();
                continue;
            }
            son->GetMinPointNormalized(min_coord);
            son->GetMaxPointNormalized(max_coord);
            pointer_type object;
            for (int j = 0; j < (int)objects_.size(); j++)
            {
                object = objects_[j];

                const int is_intersected = configuration_type::IsIntersected(object, tolerance, min_coord, max_coord);
                if (is_intersected)
                    son->Insert(object);
            }
        }

        //clear the memory of objects_ (now the objects are transfered to children)
        object_container_type temp;
        objects_.swap(temp);
    }

    ///@}
    ///@name Access
    ///@{

    unsigned char GetLevel() const
    {
        return level_;
    }

    char SetLevel(char level)
    {
        level_ = level;
        return level_;
    }

    void GetMinKey(key_type &min_key_x, key_type &min_key_y) const
    {
        min_key_x = min_key_[0];
        min_key_y = min_key_[1];
    }

    void SetMinKey(key_type min_key_x, key_type min_key_y)
    {
        min_key_[0] = min_key_x;
        min_key_[1] = min_key_y;
    }

    QuadtreeBinaryCell &rGetChild(std::size_t pos) const
    {
        return children_[pos];
    }

    QuadtreeBinaryCell *pGetChild(std::size_t pos) const
    {
        return children_ + pos;
    }

    QuadtreeBinaryCell *pGetChild(key_type x_key, key_type y_key) const
    {
        return pGetChild(GetChildIndex(x_key, y_key));
    }

    QuadtreeBinaryCell *GetChildren()
    {
        return children_;
    }

    QuadtreeBinaryCell const *GetChildren() const
    {
        return children_;
    }

    data_type *pGetData() const
    {
        return data_;
    }

    data_type **pGetDataPointer()
    {
        return &data_;
    }

    const std::vector<pointer_type> *pGetObjects() const
    {
        return &objects_;
    }

    std::vector<pointer_type> *pGetObjects()
    {
        return &objects_;
    }

    void EmptyObjects()
    {
        object_container_type tmp;
        tmp.swap(objects_);
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool IsLeaf() const
    {
        return (children_ == NULL);
    }

    bool HasChildren() const
    {
        return (children_ != NULL);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        return "QuadtreeBinaryCell";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream &rOStream) const
    {
        for (char i = ROOT_LEVEL; i > level_; i--)
        {
            rOStream << "  ";
        }
        rOStream << Info() << " at level " << static_cast<int>(level_);
    }

    /// Print object's data.

    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream << "(" << GetCoordinateNormalized(min_key_[0]) << "," << GetCoordinateNormalized(min_key_[1]) << ","
                 << "),";
        rOStream << "(" << GetCoordinateNormalized(min_key_[0]) + CalcSizeNormalized() << "," << GetCoordinateNormalized(min_key_[1]) + CalcSizeNormalized() << ")" << std::endl;

        for (std::size_t i = 0; i < CHILDREN_NUMBER; i++)
        {
            if (children_)
            {
                for (char j = ROOT_LEVEL + 1; j > level_; j--)
                {
                    rOStream << "  ";
                }

                rOStream << "child #" << i;

                children_[i].PrintData(rOStream);
            }
        }
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

    char level_;
    key_type min_key_[DIMENSION];
    QuadtreeBinaryCell *children_;
    data_type *data_;
    object_container_type objects_;

    double CalcSizeNormalized() const
    {
        const double scale = 1.00 / (1 << ROOT_LEVEL);

        return (1 << level_) * scale; // I'm doing it in this way to avoid division
    }

    double GetCoordinateNormalized(key_type key) const
    {
        const double scale = 1.00 / (1 << ROOT_LEVEL);

        return static_cast<double>(key * scale);
    }

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

    QuadtreeBinaryCell &operator=(QuadtreeBinaryCell const &rOther)
    {
        return *this;
    }

    /// Copy constructor.

    QuadtreeBinaryCell(QuadtreeBinaryCell const &rOther)
    {
    }

    ///@}

}; // Class QuadtreeBinaryCell

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

//    /// input stream function
//    inline std::istream & operator >>(std::istream& rIStream,
//            QuadtreeBinaryCell& rThis);

/// output stream function
template <class TConfiguration>
inline std::ostream &operator<<(std::ostream &rOStream,
                                const QuadtreeBinaryCell<TConfiguration> &rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_Quadtree_BYNARY_CELL_H_INCLUDED  defined
