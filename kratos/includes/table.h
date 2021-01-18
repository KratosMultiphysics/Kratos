//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_TABLE_H_INCLUDED )
#define  KRATOS_TABLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "input_output/logger.h"
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{
///@addtogroup KratosCore
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
    
/** 
 * @class Table
 * @ingroup KratosCore
 * @brief This class represents the value of its variable depending to other variable.
 * @details Table class stores the value of its second variable respect to the value of its first variable.
 * It also provides a double to double table with piecewise linear interpolator/extrapolator for getting intermediate values.
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @tparam TArgumentType The type of argument considered
 * @tparam TResultType The type of result obtained
 * @tparam TResultsColumns The number of columns considered
 */
template<class TArgumentType, class TResultType = TArgumentType, std::size_t TResultsColumns = 1>
class Table
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Table
    KRATOS_CLASS_POINTER_DEFINITION(Table);

    typedef TArgumentType argument_type; // To be STL conformance.
    typedef TResultType result_type; // To be STL conformance.

    typedef std::array<TResultType, TResultsColumns>  result_row_type;

    typedef std::pair<argument_type, result_row_type> RecordType;

    typedef std::vector<RecordType> TableContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Table() : mData()
    {
    }

    /// Destructor.
    virtual ~Table()
    {
    }

    /// Assignment operator.
    Table& operator=(Table const& rOther)
    {
        mData = rOther.mData;
        return *this;
    }

    ///@}
    ///@name Operators
    ///@{

    // This operator gives the first column result for the nearest argument found in table
    result_type const& operator()(argument_type const& X) const
    {
        return GetNearestRow(X)[0];
    }

    // This operator gives the first column result for the nearest argument found in table
    result_type& operator()(argument_type const& X)
    {
        return GetNearestRow(X)[0];
    }

    // This operator gives the result in the Jth column for the nearest argument found in table
    result_type const& operator()(argument_type const& X, std::size_t J) const
    {
        return GetNearestRow(X)[J];
    }

    // This operator gives the result in the Jth column for the nearest argument found in table
    result_type& operator()(argument_type const& X, std::size_t J)
    {
        return GetNearestRow(X)[J];
    }

    // This operator gives the row for the nearest value to argument found in table
    result_row_type const & operator[](argument_type const& X) const
    {
        return GetNearestRow(X);
    }

    // This operator gives the row for the nearest value to argument found in table
    result_row_type & operator[](argument_type& X)
    {
        return GetNearestRow(X);
    }

    ///@}
    ///@name Operations
    ///@{


    // Get the nesrest value for the given argument
    result_type& GetNearestRow(argument_type const& X)
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second;

        if(X <= mData[0].first)
            return mData[0].second;

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second : mData[i].second;

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return mData[size-1].second;
    }

    // Get the nesrest value for the given argument
    result_type const& GetNearestRow(argument_type const& X)  const
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second;

        if(X <= mData[0].first)
            return mData[0].second;

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second : mData[i].second;

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return mData[size-1].second;
    }

    // inserts a row in a sorted position where Xi-1 < X < Xi+1 and fills the first column with Y
    void insert(argument_type const& X, result_type const& Y)
    {
        result_row_type a = {{Y}};
        insert(X,a);
    }

    // inserts a row in a sorted position where Xi-1 < X < Xi+1 and fills the first column with Y
    // assumes that Y has [] operator with TResultsColumns element
    template<class TArrayType>
    void insert(argument_type const& X, TArrayType const& Y)
    {
        result_row_type a;
        for(std::size_t i = 0 ; i < TResultsColumns ; i++)
            a[i] = Y[i];
        insert(X,a);
    }

    // inserts a row in a sorted position where Xi-1 < X < Xi+1
    void insert(argument_type const& X, result_row_type const& Y)
    {
        std::size_t size = mData.size();

        if(size == 0)
            mData.push_back(RecordType(X,Y));

        else if(X <= mData[0].first)
            mData.insert(mData.begin(), RecordType(X,Y));
        else if(X <= mData.back().first)
            mData.push_back(RecordType(X,Y));
        else
            for(std::size_t i = 1 ; i < size ; i++)
                if((X > mData[i-1].first) && (X <= mData[i].first))
                {
                    mData.insert(mData.begin() + i, RecordType(X,Y));
                    break;
                }

    }


    // assumes that the X is the greater than the last argument and put the row at the end.
    // faster than insert.
    void PushBack(argument_type const& X, result_type const& Y)
    {
        result_row_type a = {{Y}};
        mData.push_back(RecordType(X,a));
    }

    // assumes that the X is the greater than the last argument and put the row at the end.
    // assumes that Y has [] operator with TResultsColumns element
    // faster than insert.
    template<class TArrayType>
    void PushBack(argument_type const& X, TArrayType const& Y)
    {
        result_row_type a;
        for(std::size_t i = 0 ; i < TResultsColumns ; i++)
            a[i] = Y[i];
        mData.push_back(RecordType(X,a));
    }

    // assumes that the X is the greater than the last argument and put the row at the end.
    // faster than insert.
    template<class TArrayType>
    void PushBack(argument_type const& X, result_row_type const& Y)
    {
        mData.push_back(RecordType(X,Y));
    }
    
    /**
     * @brief This method clears databse
     */
    void Clear()
    {
        mData.clear();
    }

    ///@}
    ///@name Access
    ///@{


    TableContainerType& Data()
    {
        return mData;
    }

    TableContainerType const& Data() const
    {
        return mData;
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
        return "Table";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    TableContainerType mData;

    ///@}
    ///@name Private Operators
    ///@{



    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        std::size_t  local_size = mData.size();

        rSerializer.save("size", local_size);

        for(auto i_row = mData.begin() ; i_row != mData.end() ; i_row++){
            rSerializer.save("Argument", i_row->first);
            for(auto j = i_row->second.begin() ; j != i_row->second.end(); j++)
                rSerializer.save("Column", j);
        }
    }

    virtual void load(Serializer& rSerializer)
    {
        std::size_t local_size;

        rSerializer.load("size", local_size);

        mData.resize(local_size);

        for(auto i_row = mData.begin() ; i_row != mData.end() ; i_row++){
            rSerializer.load("Argument", i_row->first);
            for(auto j = i_row->second.begin() ; j != i_row->second.end() ; j++)
                rSerializer.load("Column", j);
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

    /// Copy constructor.
    Table(Table const& rOther);


    ///@}

}; // Class Table

template<>
class Table<double, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Table
    KRATOS_CLASS_POINTER_DEFINITION(Table);

    typedef double argument_type; // To be STL conformance.
    typedef double result_type; // To be STL conformance.

    typedef std::array<result_type, 1>  result_row_type;

    typedef std::pair<argument_type, result_row_type> RecordType;

    typedef std::vector<RecordType> TableContainerType;

    typedef Variable<argument_type> XVariableType;
    typedef Variable<result_type> YVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Table() : mData()
    {
    }


    /// Copy constructor.
    Table(Table const& rOther): mData(rOther.mData)
    {

    }

    /// Matrix constructor. the template parameter must have (i,j) access operator and  size1 methods defined.
    template<class TMatrixType>
    Table(TMatrixType const& ThisMatrix): mData()
    {
        for(unsigned int i = 0 ; i < ThisMatrix.size1() ; i++)
            PushBack(ThisMatrix(i,0), ThisMatrix(i,1));
    }

    /// Destructor.
    virtual ~Table()
    {
    }


    /// Assignment operator.
    Table& operator=(Table const& rOther)
    {
        mData = rOther.mData;
        return *this;
    }


    ///@}
    ///@name Operators
    ///@{

    // I want to put operator(i,j) for accessing, operator(i) for first column and operator[i] for getting the complete row

    // This operator calculates the piecewise linear interpolation for
    // given argument
    result_type operator()(argument_type const& X) const
    {
        return GetValue(X);
    }

    // This operator gives the result for the nearest value to argument found in table
    result_type const & operator[](argument_type const& X) const
    {
        return GetNearestValue(X);
    }

    // This operator gives the result for the nearest value to argument found in table
    result_type & operator[](argument_type& X)
    {
        return GetNearestValue(X);
    }

    ///@}
    ///@name Operations
    ///@{

    // Get the value for the given argument using piecewise linear
    result_type GetValue(argument_type const& X) const
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second[0];

        result_type result;
        if(X <= mData[0].first)
            return Interpolate(X, mData[0].first, mData[0].second[0], mData[1].first, mData[1].second[0], result);

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return Interpolate(X, mData[i-1].first, mData[i-1].second[0], mData[i].first, mData[i].second[0], result);

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return Interpolate(X, mData[size-2].first, mData[size-2].second[0], mData[size-1].first, mData[size-1].second[0], result);
    }

    // Get the nesrest value for the given argument
    result_row_type& GetNearestRow(argument_type const& X)
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second;

        if(X <= mData[0].first)
            return mData[0].second;

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second : mData[i].second;

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return mData[size-1].second;
    }

    // Get the nesrest value for the given argument
    result_type const& GetNearestValue(argument_type const& X)  const
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second[0];

        if(X <= mData[0].first)
            return mData[0].second[0];

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second[0] : mData[i].second[0];

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return mData[size-1].second[0];
    }

    // Get the nesrest value for the given argument
    result_type & GetNearestValue(argument_type const& X)
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return mData.begin()->second[0];

        if(X <= mData[0].first)
            return mData[0].second[0];

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second[0] : mData[i].second[0];

        // now the x is outside the table and we hae to extrapolate it using last two records of table.
        return mData[size-1].second[0];
    }

    result_type& Interpolate(argument_type const& X, argument_type const& X1, result_type const& Y1, argument_type const& X2, result_type const& Y2, result_type& Result) const
    {
        const double epsilon = 1e-12;

        double dx = X2 - X1;
        result_type dy = Y2 - Y1;

        double scale = 0.00;

        if (dx > epsilon)
            scale = (X - X1) / dx;

        Result = Y1 + dy * scale;

        return Result;

    }

    // inserts a row in a sorted position where Xi-1 < X < Xi+1 and fills the first column with Y
    void insert(argument_type const& X, result_type const& Y)
    {
        result_row_type a = {{Y}};
        insert(X,a);
    }


    // inserts a row in a sorted position where Xi-1 < X < Xi+1
    void insert(argument_type const& X, result_row_type const& Y)
    {
        std::size_t size = mData.size();

        if(size == 0)
            mData.push_back(RecordType(X,Y));
        else if(X <= mData[0].first)
            mData.insert(mData.begin(), RecordType(X,Y));
        else if(X > mData.back().first)
            mData.push_back(RecordType(X,Y));
        else
            for(std::size_t i = 1 ; i < size ; i++)
                if((X > mData[i-1].first) && (X <= mData[i].first))
                {
                    mData.insert(mData.begin() + i, RecordType(X,Y));
                    break;
                }
    }

    // assumes that the X is the greater than the last argument and put the row at the end.
    // faster than insert.
    void PushBack(argument_type const& X, result_type const& Y)
    {
        result_row_type a = {{Y}};
        mData.push_back(RecordType(X,a));
    }

     // Get the derivative for the given argument using piecewise linear
    result_type GetDerivative(argument_type const& X) const
    {
        std::size_t size = mData.size();

        KRATOS_ERROR_IF(size == 0) << "Get value from empty table" << std::endl;

        if(size==1) // constant table. Returning the only value we have.
            return 0.0;

        result_type result;
        if(X <= mData[0].first)
            //return Interpolate(X, mData[0].first, mData[0].second[0], mData[1].first, mData[1].second[0], result);
            return 0.0;

        for(std::size_t i = 1 ; i < size ; i++)
            if(X <= mData[i].first)
                return InterpolateDerivative( mData[i-1].first, mData[i-1].second[0], mData[i].first, mData[i].second[0], result);

        // If it lies outside the table values we will return 0.0.
        return 0.0;
    }
     result_type& InterpolateDerivative( argument_type const& X1, result_type const& Y1, argument_type const& X2, result_type const& Y2, result_type& Result) const
    {
        const double epsilon = 1e-12;
        argument_type dx = X2 - X1;
        result_type dy = Y2 - Y1;
        if (dx < epsilon)
        {
            dx=epsilon;
            KRATOS_WARNING("") 
            << "*******************************************\n"
            << "*** ATTENTION: SMALL dX WHEN COMPUTING  ***\n"
            << "*** DERIVATIVE FROM TABLE. SET TO 1E-12 ***\n"
            << "*******************************************" <<std::endl;
        }
        Result= dy/dx;
        return Result;
    }
    
    /**
     * @brief This method clears databse
     */
    void Clear()
    {
        mData.clear();
    }
    
    ///@}
    ///@name Access
    ///@{


    TableContainerType& Data()
    {
        return mData;
    }

    TableContainerType const& Data() const
    {
        return mData;
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
        return "Piecewise Linear Table";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for(std::size_t i = 0 ; i < mData.size() ; i++)
            rOStream << mData[i].first << "\t\t" << mData[i].second[0] << std::endl;
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

    TableContainerType mData;

    ///@}
    ///@name Private Operators
    ///@{



    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        std::size_t  local_size = mData.size();

        rSerializer.save("size", local_size);

        for(auto i_row = mData.begin() ; i_row != mData.end() ; i_row++){
            rSerializer.save("Argument", i_row->first);
            for(auto j = i_row->second.begin() ; j != i_row->second.end(); j++)
                rSerializer.save("Column", *j);
        }
    }

    void load(Serializer& rSerializer)
    {
        std::size_t local_size;

        rSerializer.load("size", local_size);

        mData.resize(local_size);

        for(auto i_row = mData.begin() ; i_row != mData.end() ; i_row++){
            rSerializer.load("Argument", i_row->first);
            for(auto j = i_row->second.begin() ; j != i_row->second.end() ; j++)
                rSerializer.load("Column", *j);
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


    ///@}

}; // Class Table

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TArgumentType, class TResultType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Table<TArgumentType, TResultType>& rThis);

/// output stream function
template<class TArgumentType, class TResultType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Table<TArgumentType, TResultType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TABLE_H_INCLUDED  defined 


