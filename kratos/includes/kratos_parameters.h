//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_KRATOS_PARAMETERS_H_INCLUDED )
#define  KRATOS_KRATOS_PARAMETERS_H_INCLUDED

// System includes

// External includes
#include "json/json_fwd.hpp" // Import forward declaration nlohmann json library

// Project includes
#include "includes/serializer.h"
#include "includes/ublas_interface.h"

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
 * @class Parameters
 * @ingroup KratosCore
 * @brief This class provides to Kratos a data structure for I/O based on the standard of JSON
 * @details In computing, JavaScript Object Notation or JSON is an open-standard file format that uses human-readable text to transmit data objects consisting of attribute–value pairs and array data types (or any other serializable value). It is a very common data format used for asynchronous browser–server communication, including as a replacement for XML in some AJAX-style systems. More info: https://json.org/
 * This class uses nlohmann JSON header only library
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) Parameters
{
private:
    ///@name Nested clases
    ///@{
    /**
     * @class iterator_adaptor
     * @ingroup KratosCore
     * @brief This nested class can be used to adapt a Parameter iterator
     * @author Riccardo Rossi
     */
    class KRATOS_API(KRATOS_CORE) iterator_adaptor
        : public std::iterator<std::forward_iterator_tag, Parameters>
    {
        ///@name Type Definitions
        ///@{

        using value_iterator = nlohmann::detail::iter_impl<nlohmann::json>; /// Iterator definition

        ///@}
        ///@name Member Variables
        ///@{

        std::size_t mDistance = 0;                       /// The iterator distance
        nlohmann::json& mrValue;                         /// The original container
        std::unique_ptr<Parameters> mpParameters;        /// The unique pointer to the base Parameter

        ///@}
    public:
        ///@name Life Cycle
        ///@{

        /**
         * @brief Default constructor (iterator + root Parameter)
         * @param itValue The iterator to adapt
         * @param pRoot The root Parameter pointer
         */
        iterator_adaptor(value_iterator itValue, nlohmann::json* pValue,  Kratos::shared_ptr<nlohmann::json> pRoot);

        /**
         * @brief Default constructor (just iterator)
         * @param itValue The iterator to adapt
         */
        iterator_adaptor(const iterator_adaptor& itValue);

        ///@}
        ///@name Operators
        ///@{

        /**
         * @brief operator ++
         * @details This adds one to the current iterator
         * @return The next iterator
         */
        iterator_adaptor& operator++();

        /**
         * @brief operator ++int
         * @details This adds N to the current iterator
         * @param int N increment of iterations
         * @return The +N iterator
         */
        iterator_adaptor operator++(int);

        /**
         * @brief operator ==
         * @details This operator check if the iterator is equal to another given iterator
         * @return True if equal, false otherwise
         */
        bool operator==(const iterator_adaptor& rhs) const;

        /**
         * @brief operator !=
         * @details This operator check if the iterator is not equal to another given iterator
         * @return True if not equal, false otherwise
         */
        bool operator!=(const iterator_adaptor& rhs) const;

        /**
         * @brief operator*
         * @details This operator returns the pointer of a given iterator
         * @return The Pointer of the given iterator
         */
        Parameters& operator*() const;

        /**
         * @brief operator ->
         * @details This operator acces to the pointer of the Parameter
         * @return The pointer of the parameter
         */
        Parameters* operator->() const;

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief This method returs the current iterator
         * @return The current iterator
         */
        inline value_iterator GetCurrentIterator() const;

        /**
         * @brief This method returns the key of the current Parameter iterator
         * @return The key (name) of the Parameter iterator
         */
        const std::string name();

        ///@}
    };

    /**
     * @class const_iterator_adaptor
     * @ingroup KratosCore
     * @brief This nested class can be used to adapt a Parameter constant iterator
     * @author Riccardo Rossi
     */
    class KRATOS_API(KRATOS_CORE) const_iterator_adaptor
        : public std::iterator<std::forward_iterator_tag, Parameters>
    {
        ///@name Type Definitions
        ///@{

        using value_iterator = nlohmann::detail::iter_impl<const nlohmann::json>; /// Iterator definition

        ///@}
        ///@name Member Variables
        ///@{

        std::size_t mDistance = 0;                       /// The iterator distance
        nlohmann::json& mrValue;                         /// The original container
        std::unique_ptr<Parameters> mpParameters;        /// The unique pointer to the base Parameter

        ///@}
    public:
        ///@name Life Cycle
        ///@{

        /**
         * @brief Default constructor (constant iterator + root Parameter)
         * @param itValue The iterator to adapt
         * @param pRoot The root Parameter pointer
         */
        const_iterator_adaptor(value_iterator itValue, nlohmann::json* pValue,  Kratos::shared_ptr<nlohmann::json> pRoot);

        /**
         * @brief Default constructor (just constant iterator)
         * @param itValue The iterator to adapt
         * @todo Use copy constructor in the following method
         */
        const_iterator_adaptor(const const_iterator_adaptor& itValue);

        ///@}
        ///@name Operators
        ///@{

        /**
         * @brief operator ++
         * @details This adds one to the current iterator
         * @return The next iterator (const)
         */
        const_iterator_adaptor& operator++();

        /**
         * @brief operator ++int
         * @details This adds N to the current iterator
         * @param int N increment of iterations
         * @return The +N iterator (const)
         */
        const_iterator_adaptor operator++(int);

        /**
         * @brief operator ==
         * @details This operator check if the iterator is equal to another given iterator
         * @return True if equal, false otherwise
         */
        bool operator==(const const_iterator_adaptor& rhs) const;

        /**
         * @brief operator !=
         * @details This operator check if the iterator is not equal to another given iterator
         * @return True if not equal, false otherwise
         */
        bool operator!=(const const_iterator_adaptor& rhs) const;

        /**
         * @brief operator*
         * @details This operator returns the pointer of a given iterator
         * @return The Pointer of the given iterator
         */
        const Parameters& operator*() const;

        /**
         * @brief operator ->
         * @details This operator acces to the pointer of the Parameter
         * @return The pointer of the parameter
         */
        const Parameters* operator->() const;

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief This method returs the current iterator
         * @return The current iterator
         */
        inline value_iterator GetCurrentIterator() const;

        /**
         * @brief This method returns the key of the current Parameter iterator
         * @return The key (name) of the Parameter iterator
         */
        const std::string name();

        ///@}
    };

    ///@}

public:
    ///@name Type Definitions
    ///@{

    /// Index definition
    typedef std::size_t IndexType;

    /// Size definition
    typedef std::size_t SizeType;

    /// Pointer definition of MmgProcess
    KRATOS_CLASS_POINTER_DEFINITION(Parameters);

    /// Definiton of the iterators
    using iterator = iterator_adaptor;
    using const_iterator = const_iterator_adaptor;

    /// Iterators from nlohmann::json
    typedef nlohmann::detail::iter_impl<nlohmann::json> json_iterator;
    typedef nlohmann::detail::iter_impl<const nlohmann::json> json_const_iterator;
    typedef nlohmann::detail::iteration_proxy<json_iterator> json_iteration_proxy;
    typedef nlohmann::detail::iteration_proxy<json_const_iterator> json_const_iteration_proxy;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @brief It assigns null pointers to the member variables
     */
    Parameters();

    /**
     * @brief String constructor. It takes a string as input, which parses into a nlohmann::json class
     * @param rJsonString The string to be parsed into a nlohmann::json class
     */
    Parameters(const std::string& rJsonString);

    /// Copy constructor.
    Parameters(Parameters const& rOther);

    /// Move constructor.
    Parameters(Parameters&& rOther);

    /// Destructor.
    virtual ~Parameters()
    {
//         delete[] mpValue;
//         mpRoot = nullptr;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Parameters& operator=(Parameters const& rOther);

    /**
     * @brief This metrod returns the Parameter corresponding to a given key
     * @param rEntry The key identifier of the parameter
     * @return The desired Parameter
     */
    Parameters operator[](const std::string& rEntry);

    /**
     * @brief This method allows to acces to an array item with the operator []
     * @param Index The index of the term of interest
     * @return The desired Parameter
     */
    Parameters operator[](const IndexType Index);

    /**
     * @brief This is the move operator
     * @param rOther The othe parameter to compute the move
     */
    Parameters& operator=(Parameters&& rOther);

    ///@}
    ///@name Operations
    ///@{

    //generates a clone of the current document
    Parameters Clone();

    /**
     * @brief This method returns a string with the corresponding text to the equivalent *.json file
     * @return The corresponding text
     */
    const std::string WriteJsonString() const;

    /**
     * @brief This method returns a string with the corresponding text to the equivalent *.json file (this version is prettier, and considers tabulations)
     * @return The corresponding text
     */
    const std::string PrettyPrintJsonString() const;

    /**
     * @brief This method returns the Parameter corresponding to a certain entry
     * @param rEntry The key identifier of the parameter
     * @return The corresponding parameter
     */
    Parameters GetValue(const std::string& rEntry);

    /**
     * @brief This method sets an existing parameter with a given parameter
     * @param rEntry The key identifier of the parameter
     * @param rOtherValue The value to set
     */
    void SetValue(
        const std::string& rEntry,
        const Parameters& rOtherValue
        );

    /**
     * @brief This method sets a non-existing parameter with a given parameter
     * @param rEntry The key identifier of the parameter
     * @param rOtherValue The value to set
     */
    void AddValue(
        const std::string& rEntry,
        const Parameters& rOtherValue
        );

    /**
     * @brief This method adds an empty parameter
     * @param rEntry The key identifier of the parameter
     */
    Parameters AddEmptyValue(const std::string& rEntry);

    /**
     * @brief This method removes an entry of the Parameters given a certain key
     * @param rEntry The key identifier of the parameter
     * @return False if failed, true otherwise
     */
    bool RemoveValue(const std::string& rEntry);

    /**
     * @brief This method returns the items of the current parameter
     * @return The items of the current Parameter
     */
    json_iteration_proxy items() noexcept;

    /**
     * @brief This method returns the items of the current parameter (const)
     * @return The items of the current Parameter (const)
     */
    json_const_iteration_proxy items() const noexcept;

    /**
     * @brief This method checks if the Parameter contains a certain entry
     * @param rEntry The key identifier of the parameter
     * @return True if it contains, false otherwise
     */
    bool Has(const std::string& rEntry) const;

    /**
     * @brief This method checks if the parameter is a null
     * @return True if it is null, false otherwise
     */
    bool IsNull() const;

    /**
     * @brief This method checks if the parameter is a number
     * @return True if it is a number, false otherwise
     */
    bool IsNumber() const;

    /**
     * @brief This method checks if the parameter is a double
     * @return True if it is a double, false otherwise
     */
    bool IsDouble() const;

    /**
     * @brief This method checks if the parameter is a integer
     * @return True if it is a integer, false otherwise
     */
    bool IsInt() const;

    /**
     * @brief This method checks if the parameter is a boolean
     * @return True if it is a boolean, false otherwise
     */
    bool IsBool() const;

    /**
     * @brief This method checks if the parameter is a string
     * @return True if it is a string, false otherwise
     */
    bool IsString() const;

    /**
     * @brief This method checks if the parameter is an array
     * @return True if it is an array, false otherwise
     */
    bool IsArray() const;

    /**
     * @brief This method checks if the parameter is a vector
     * @return True if it is a vector, false otherwise
     */
    bool IsVector() const;

    /**
     * @brief This method checks if the parameter is a matrix
     * @return True if it is a matrix, false otherwise
     */
    bool IsMatrix() const;

    /**
     * @brief This method checks if the parameter is a subparameter
     * @return True if it is a suparameter, false otherwise
     */
    bool IsSubParameter() const;

    /**
     * @brief This method returns the double contained in the current Parameter
     * @return The double value
     */
    double GetDouble() const;

    /**
     * @brief This method returns the integer contained in the current Parameter
     * @return The integer value
     */
    int GetInt() const;

    /**
     * @brief This method returns the boolean contained in the current Parameter
     * @return The boolean value
     */
    bool GetBool() const;

    /**
     * @brief This method returns the string contained in the current Parameter
     * @return The string value
     */
    std::string GetString() const;

    /**
     * @brief This method returns the vector contained in the current Parameter
     * @return The vector value
     */
    Vector GetVector() const;

    /**
     * @brief This method returns the matrix contained in the current Parameter
     * @return The matrix value
     */
    Matrix GetMatrix() const;

    /**
     * @brief This method sets the double contained in the current Parameter
     * @param Value The double value
     */
    void SetDouble(const double Value);

    /**
     * @brief This method sets the integer contained in the current Parameter
     * @param Value The integer value
     */
    void SetInt(const int Value);

    /**
     * @brief This method sets the bool contained in the current Parameter
     * @param Value The bool value
     */
    void SetBool(const bool Value);

    /**
     * @brief This method sets the string contained in the current Parameter
     * @param rValue The string value
     */
    void SetString(const std::string& rValue);

    /**
     * @brief This method sets the vector contained in the current Parameter
     * @param rValue The vector value
     */
    void SetVector(const Vector& rValue);

    /**
     * @brief This method sets the matrix contained in the current Parameter
     * @param Value The matrix value
     */
    void SetMatrix(const Matrix& rValue);

    /**
     * @brief This returns the begin iterator
     * @return The begin iterator
     */
    iterator begin();

    /**
     * @brief This returns the end iterator
     * @return The end iterator
     */
    iterator end();

    /**
     * @brief This returns the constant begin iterator
     * @return The constant begin iterator
     */
    const_iterator begin() const;

    /**
     * @brief This returns the constant end iterator
     * @return The constant end iterator
     */
    const_iterator end() const;

    /**
     * @brief This method returns the total size of the current array parameter
     * @return The size of the current array parameter
     */
    SizeType size() const;

    /**
     * @brief This method does a swap between two parameters
     * @param rOther The othe parameter to compute the swap
     */
    void swap(Parameters& rOther) noexcept;

    /**
     * @brief This method resets the whole parameter (it assigns an empty parameter)
     */
    void Reset() noexcept;

    /**
     * @brief This method returns an array item given an index
     * @param Index The index of the parameter to obtain
     * @return The parameter corresponding to the given index
     */
    Parameters GetArrayItem(const IndexType Index);

    /**
     * @brief This method sets an array item given an index
     * @param Index The index of the parameter to set
     * @param rOtherArrayItem The parameter corresponding to the given index
     */
    void SetArrayItem(
        const IndexType Index,
        const Parameters& rOtherArrayItem
        );

    /**
     * @brief This method add a new entry with no value assigned
     * @param rEntry The key identifier of the parameter
     */
    void AddEmptyArray(const std::string& rEntry);

    /**
     * @brief This method appends into an array a double value
     * @param Value The double value to append
     */
    void Append(const double Value);

    /**
     * @brief This method appends into an array a integer value
     * @param Value The integer value to append
     */
    void Append(const int Value);

    /**
     * @brief This method appends into an array a boolean value
     * @param Value The boolean value to append
     */
    void Append(const bool Value);

    /**
     * @brief This method appends into an array a string value
     * @param rValue The string value to append
     */
    void Append(const std::string& rValue);

    /**
     * @brief This method appends into an array a vector value
     * @param rValue The vector value to append
     */
    void Append(const Vector& rValue);

    /**
     * @brief This method appends into an array a matrix value
     * @param rValue The matrix value to append
     */
    void Append(const Matrix& rValue);

    /**
     * @brief This method appends into an array a Parameter value
     * @param rValue The Parameter value to append
     */
    void Append(const Parameters& rValue);

    /**
     * @brief This method looks in a recursive way in the json structure
     * @param rBaseValue The value where to find
     * @param rValueToFind The value to look
     */
    void RecursivelyFindValue(
        const nlohmann::json& rBaseValue,
        const nlohmann::json& rValueToFind
        ) const;

    /**
     * @brief Checks if the names and values are the same, no importance to the order.
     * @details Lists have to be ordered, though! Take into account that in Kratos some physical vectors are represented with a list.
     * @param rParameters The parameters to be checked
     * @return True if it has, false othersise
     */
    bool IsEquivalentTo(Parameters& rParameters);

    /**
     * @brief Checks if the names and the type of values are the same, no importance to the order.
     * @details Lists have to be ordered, though! Take into account that in Kratos some physical vectors are represented with a list.
     * @param rParameters The parameters to be checked
     * @return True if it has, false othersise
     */
    bool HasSameKeysAndTypeOfValuesAs(Parameters& rParameters);

    /**
     * @brief This function is designed to verify that the parameters under testing match the form prescribed by the rDefaultParameters.
     * @details If the parameters contain values that do not appear in the rDefaultParameters, an error is thrown, whereas if a parameter is found in the rDefaultParameters but not in the Parameters been tested, it is copied to the parameters.
     * This version of the function only walks one level, without descending in the branches
     * @param rDefaultParameters Parameters of reference which we use to check
     */
    void ValidateAndAssignDefaults(Parameters& rDefaultParameters);
    /**
     * @brief This function is designed to verify that the parameters under testing match the form prescribed by the defaults.
     * @details If the parameters contain values that do not appear in the defaults, an error is thrown, whereas if a parameter is found in the defaults but not in the Parameters been tested, it is copied to the parameters.
     * This version walks and validates the entire json tree below the point at which the function is called
     * @param rDefaultParameters Parameters of reference which we use to check
     */
    void RecursivelyValidateAndAssignDefaults(Parameters& rDefaultParameters);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return this->PrettyPrintJsonString();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Parameters Object " << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Parameters Object " << Info();
    };

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

    friend class Serializer;

    void save(Serializer& rSerializer) const;

    void load(Serializer& rSerializer);

    ///@}
    ///@name Member Variables
    ///@{

    nlohmann::json* mpValue;                   /// This is where the json is actually stored
    Kratos::shared_ptr<nlohmann::json> mpRoot; /// This is a shared pointer to the root structure (this is what allows us to acces in a tree structure to the JSON database)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Direct constructor. It takes as parameters the "member" variables of the Parameters class
     * @param pValue The nlohmann::json class raw pointer
     * @param pRoot A shared pointer to a nlohmann::json class
     * @warning Please DO NOT use this constructor. It assumes nlohmann::json and hence it should be considered as an implementation detail
     */
    Parameters(nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot);

    /**
     * @brief Direct constructor. It takes as parameters the "member" variables of the Parameters class
     * @param itValue The nlohmann::json class iterator
     * @param pRoot A shared pointer to a nlohmann::json class
     * @warning Please DO NOT use this constructor. It assumes nlohmann::json and hence it should be considered as an implementation detail
     */
    Parameters(json_iterator itValue, nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot);

    /**
     * @brief Direct constructor. It takes as parameters the "member" variables of the Parameters class
     * @param itValue The nlohmann::json class iterator
     * @param pRoot A shared pointer to a nlohmann::json class
     * @warning Please DO NOT use this constructor. It assumes nlohmann::json and hence it should be considered as an implementation detail
     */
    Parameters(json_const_iterator itValue, nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot);

    /**
     * @brief This method is created in order to access from the iterators to the database
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    nlohmann::json* GetUnderlyingStorage();

    /**
     * @brief This method is created in order to access from the iterators to the database
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    nlohmann::json* GetUnderlyingStorage() const;

    /**
     * @brief This method is created in order to set the database
     * @param pNewValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    void SetUnderlyingSotrage(nlohmann::json* pNewValue);

    /**
     * @brief This method is created in order to access from the iterators to the database
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    Kratos::shared_ptr<nlohmann::json> GetUnderlyingRootStorage();

    /**
     * @brief This method is created in order to access from the iterators to the database
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    Kratos::shared_ptr<nlohmann::json> GetUnderlyingRootStorage() const;

    /**
     * @brief This method is created in order to set the database
     * @param pNewValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    void SetUnderlyingRootStorage(Kratos::shared_ptr<nlohmann::json> pNewValue);
    /**
     * @brief This method sets the database from other Parameters
     * @param rOtherValue The database to copy
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    void InternalSetValue(const Parameters& rOtherValue);

}; // Parameters class

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Parameters& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Parameters& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_KRATOS_PARAMETERS_H_INCLUDED  defined
