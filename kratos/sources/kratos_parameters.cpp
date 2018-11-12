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

// System includes

// External includes
#include "json/json.hpp" // Import nlohmann json library

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/define.h"
#include "input_output/logger.h"

namespace Kratos
{

Parameters::iterator_adaptor::iterator_adaptor(Parameters::iterator_adaptor::value_iterator itValue, nlohmann::json* pValue,  Kratos::shared_ptr<nlohmann::json> pRoot) :mDistance(std::distance(pValue->begin(), itValue)), mrValue(*pValue), mpParameters(new Parameters(itValue, pValue, pRoot)) {}

/***********************************************************************************/
/***********************************************************************************/

Parameters::iterator_adaptor::iterator_adaptor(const Parameters::iterator_adaptor& itValue) : mDistance(itValue.mDistance), mrValue(itValue.mrValue),  mpParameters(new Parameters(itValue->GetUnderlyingStorage(), itValue->GetUnderlyingRootStorage())) {}

/***********************************************************************************/
/***********************************************************************************/

Parameters::iterator_adaptor& Parameters::iterator_adaptor::operator++()
{
    ++mDistance;
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::iterator_adaptor Parameters::iterator_adaptor::operator++(int)
{
    iterator_adaptor tmp(*this);
    operator++();
    return tmp;
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::iterator_adaptor::operator==(const Parameters::iterator_adaptor& rhs) const
{
    return mDistance == rhs.mDistance;
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::iterator_adaptor::operator!=(const Parameters::iterator_adaptor& rhs) const
{
    return mDistance != rhs.mDistance;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters& Parameters::iterator_adaptor::operator*() const
{
    auto it = GetCurrentIterator();
    if (it != mrValue.end()) {
        mpParameters->mpValue = &(*it);
    }
    return *mpParameters;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters* Parameters::iterator_adaptor::operator->() const
{
    auto it = GetCurrentIterator();
    if (it != mrValue.end()) {
        mpParameters->mpValue = &(*it);
    }
    return mpParameters.get();
}

/***********************************************************************************/
/***********************************************************************************/

inline Parameters::iterator_adaptor::value_iterator Parameters::iterator_adaptor::GetCurrentIterator() const
{
    auto it = mrValue.begin();
    for (std::size_t i = 0; i < mDistance; ++i)
        it++;
    return it;
}

/***********************************************************************************/
/***********************************************************************************/

const std::string Parameters::iterator_adaptor::name()
{
    auto it = GetCurrentIterator();
    return it.key();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator_adaptor::const_iterator_adaptor(value_iterator itValue, nlohmann::json* pValue,  Kratos::shared_ptr<nlohmann::json> pRoot) : mDistance(std::distance(pValue->cbegin(), itValue)), mrValue(*pValue), mpParameters(new Parameters(itValue, pValue, pRoot)) {}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator_adaptor::const_iterator_adaptor(const const_iterator_adaptor& itValue) : mDistance(itValue.mDistance), mrValue(itValue.mrValue), mpParameters(new Parameters(itValue->GetUnderlyingStorage(), itValue->GetUnderlyingRootStorage()))  {}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator_adaptor& Parameters::const_iterator_adaptor::operator++()
{
    ++mDistance;
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator_adaptor Parameters::const_iterator_adaptor::operator++(int)
{
    const_iterator_adaptor tmp(*this);
    operator++();
    return tmp;
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::const_iterator_adaptor::operator==(const Parameters::const_iterator_adaptor& rhs) const
{
    return mDistance == rhs.mDistance;
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::const_iterator_adaptor::operator!=(const Parameters::const_iterator_adaptor& rhs) const
{
    return mDistance != rhs.mDistance;
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters& Parameters::const_iterator_adaptor::operator*() const
{
    auto it = GetCurrentIterator();
    if (it != mrValue.cend())
        mpParameters->mpValue = const_cast<nlohmann::json*>(&(*it));
    return *mpParameters;
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters* Parameters::const_iterator_adaptor::operator->() const
{
    auto it = GetCurrentIterator();
    if (it != mrValue.cend())
        mpParameters->mpValue = const_cast<nlohmann::json*>(&(*it));
    return mpParameters.get();
}

/***********************************************************************************/
/***********************************************************************************/

inline Parameters::const_iterator_adaptor::value_iterator Parameters::const_iterator_adaptor::GetCurrentIterator() const
{
    auto it = mrValue.cbegin();
    for (std::size_t i = 0; i < mDistance; ++i)
        it++;
    return it;
}

/***********************************************************************************/
/***********************************************************************************/

const std::string Parameters::const_iterator_adaptor::name()
{
    auto it = GetCurrentIterator();
    return it.key();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters()
{
    mpRoot = nullptr;
    mpValue = nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters(const std::string& rJsonString)
{
    mpRoot = Kratos::make_shared<nlohmann::json>(nlohmann::json::parse( rJsonString ));
    mpValue = mpRoot.get();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters(Parameters const& rOther)
{
    //TODO: verify if mpValue is not null and eventually destruct correctly the data
    mpRoot = rOther.mpRoot;
    mpValue = rOther.mpValue;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters& Parameters::operator=(Parameters const& rOther)
{
    if(mpRoot.get() ==  mpValue || mpRoot == nullptr) {
        mpRoot = Kratos::make_shared<nlohmann::json>(nlohmann::json::parse(rOther.WriteJsonString()));
        mpValue = mpRoot.get();
    } else {
        *mpValue = nlohmann::json( nlohmann::json::parse( rOther.WriteJsonString() ) );
        // note that mpRoot is unchanged
    }

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::operator[](const std::string& rEntry)
{
    return this->GetValue(rEntry);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::operator[](const IndexType Index)
{
    return this->GetArrayItem(Index);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::Clone()
{
    //TODO: make a clone
    //TODO: find a better way to make the copy
    return Parameters(mpValue->dump());                     //new json(*mpValue));
}

/***********************************************************************************/
/***********************************************************************************/

const std::string Parameters::WriteJsonString() const
{
    return mpValue->dump();
}

/***********************************************************************************/
/***********************************************************************************/

const std::string Parameters::PrettyPrintJsonString() const
{
    return mpValue->dump(4);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::GetValue(const std::string& rEntry)
{
    auto j = mpValue->find(rEntry);
    KRATOS_ERROR_IF(j == mpValue->end()) << "Getting a value that does not exist. entry string : " << rEntry << std::endl;
    return Parameters(&(*j), mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetValue(
    const std::string& rEntry,
    const Parameters& rOtherValue
    )
{
    KRATOS_ERROR_IF(mpValue->find(rEntry) == mpValue->end()) << "Value must exist to be set. Use AddValue instead" << std::endl;
    (*mpValue)[rEntry] = *(rOtherValue.mpValue);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::AddValue(
    const std::string& rEntry,
    const Parameters& rOtherValue
    )
{
    if(mpValue->find(rEntry) == mpValue->end()) {
        (*mpValue)[rEntry] = *(rOtherValue.mpValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::AddEmptyValue(const std::string& rEntry)
{
    if(this->Has(rEntry) == false) {
        return Parameters(&(*mpValue)[rEntry], mpRoot);
    }
    return this->GetValue(rEntry);
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::RemoveValue(const std::string& rEntry)
{
    return static_cast<bool>(mpValue->erase(rEntry));
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::json_iteration_proxy Parameters::items() noexcept
{
    return mpValue->items();
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::json_const_iteration_proxy Parameters::items() const noexcept
{
    return json_const_iteration_proxy(*mpValue);
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::Has(const std::string& rEntry) const
{
    return mpValue->find(rEntry) != mpValue->end();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsNull() const
{
    return mpValue->is_null();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsNumber() const
{
    return mpValue->is_number();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsDouble() const
{
    return mpValue->is_number_float();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsInt() const
{
    return mpValue->is_number_integer();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsBool() const
{
    return mpValue->is_boolean();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsString() const
{
    return mpValue->is_string();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsArray() const
{
    return mpValue->is_array();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsVector() const
{
    if (!mpValue->is_array())
        return false;

    auto& r_array = (*mpValue);
    for (IndexType i = 0; i < mpValue->size(); ++i) {
        if (!r_array[i].is_number())
            return false;
    }
    return true; // All entries are numbers or Vector is empty
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsMatrix() const
{
    if (!mpValue->is_array()) // mpValue != [ ... ]
        return false;

    const SizeType nrows = mpValue->size();
    if (nrows == 0) // mpValue is an empty array/vector => "[]"
        return false;

    for (IndexType i = 0; i < nrows; ++i) {
        auto& row_i = (*mpValue)[i];
        if (!row_i.is_array())
            return false;

        IndexType ncols = row_i.size();
        if (ncols != (*mpValue)[0].size()) // Compare number of columns to first row
            return false;                  // Number of columns is not consistent

        for (IndexType j = 0; j < ncols; ++j) { // Check all values in column
            if (!row_i[j].is_number())
            return false;
        }
    }

    return true; // All entries are numbers or Matrix is empty ([[]] or
                    // [[],[],[],...])
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsSubParameter() const
{
    return mpValue->is_object();
}

/***********************************************************************************/
/***********************************************************************************/

double Parameters::GetDouble() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_number()) << "Argument must be a number" << std::endl;
    return mpValue->get<double>();
}

/***********************************************************************************/
/***********************************************************************************/

int Parameters::GetInt() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_number()) << "Argument must be a number" << std::endl;
    return mpValue->get<int>();
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::GetBool() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_boolean()) << "Argument must be a bool" << std::endl;
    return mpValue->get<bool>();
}

/***********************************************************************************/
/***********************************************************************************/

std::string Parameters::GetString() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_string()) << "Argument must be a string" << std::endl;
    return mpValue->get<std::string>();
}

/***********************************************************************************/
/***********************************************************************************/

Vector Parameters::GetVector() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "Argument must be a Vector (a json list)" << std::endl;

    const SizeType size = mpValue->size();

    Vector aux_V(size);

    for (IndexType i = 0; i < size; ++i) {
        KRATOS_ERROR_IF_NOT((*mpValue)[i].is_number()) << "Entry " << i << " is not a number!" << std::endl;
        aux_V(i) = (*mpValue)[i].get<double>();
    }

    return aux_V;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix Parameters::GetMatrix() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "Argument must be a Matrix (a json list of lists)" << std::endl;

    const SizeType nrows = mpValue->size();
    KRATOS_ERROR_IF(nrows == 0) << "Argument must be a Matrix (a json list of lists)" << std::endl;

    IndexType ncols = 0;
    if ((*mpValue)[0].is_array())
        ncols = (*mpValue)[0].size();

    Matrix aux_A(nrows, ncols);

    for (IndexType i = 0; i < nrows; ++i) {
        auto &row_i = (*mpValue)[i];
        KRATOS_ERROR_IF_NOT(row_i.is_array()) << "Not an array on row " << i << std::endl;
        KRATOS_ERROR_IF_NOT(row_i.size() == ncols) << "Wrong size of row " << i << std::endl;
        for (IndexType j = 0; j < ncols; ++j) {
            KRATOS_ERROR_IF_NOT((row_i)[j].is_number()) << "Entry (" << i << "," << j << ") is not a number!" << std::endl;
            aux_A(i, j) = (row_i)[j].get<double>();
        }
    }

    return aux_A;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetDouble(const double Value)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_number()) << "Parameter must be a number" << std::endl;
    *mpValue=Value;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetInt(const int Value)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_number()) << "Parameter must be a number" << std::endl;
    *mpValue=Value;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetBool(const bool Value)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_boolean()) << "Parameter must be a bool" << std::endl;
    *mpValue=Value;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetString(const std::string& rValue)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_string()) << "Parameter must be a string" << std::endl;
    *mpValue=rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetVector(const Vector& rValue)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "Parameter must be a Vector (a json list)" << std::endl;

    const SizeType size = rValue.size();

    nlohmann::json j_array(0.0, size);
    (*mpValue) = j_array;

    for (IndexType i = 0; i < size; ++i) {
        (*mpValue)[i] = rValue[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetMatrix(const Matrix& rValue)
{
//     KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "Parameter must be a Matrix (a json list of lists)" << std::endl;

    const SizeType nrows = rValue.size1();
    const SizeType ncols = rValue.size2();

    nlohmann::json j_col_array(0.0, ncols);
    nlohmann::json j_row_array(0.0, nrows);
    (*mpValue) = j_row_array;

    for (IndexType i = 0; i < nrows; ++i) {
        (*mpValue)[i] = j_col_array;

        for (IndexType j = 0; j < ncols; ++j) {
            (*mpValue)[i][j] = rValue(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::iterator Parameters::begin()
{
    return iterator(mpValue->begin(), mpValue, mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::iterator Parameters::end()
{
    return iterator(mpValue->end(), mpValue, mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator Parameters::begin() const
{
    return const_iterator(mpValue->cbegin(), mpValue, mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::const_iterator Parameters::end() const
{
    return const_iterator(mpValue->cend(), mpValue, mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::SizeType Parameters::size() const
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array())  << "Size can only be queried if the value if of Array type" << std::endl;
    return mpValue->size();
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::swap(Parameters& rOther) noexcept
{
    std::swap(mpValue, rOther.mpValue);
    std::swap(mpRoot, rOther.mpRoot);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Reset() noexcept
{
    Parameters p;
    swap(p);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters Parameters::GetArrayItem(const IndexType Index)
{
    if(mpValue->is_array() == false) {
        KRATOS_ERROR << "GetArrayItem only makes sense if the value if of Array type" << std::endl;
    } else {
        KRATOS_ERROR_IF(Index >= mpValue->size()) << "Index exceeds array size. Index value is : " << Index << std::endl;
        return Parameters(&((*mpValue)[Index]), mpRoot);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetArrayItem(
    const IndexType Index,
    const Parameters& rOtherArrayItem
    )
{
    if(mpValue->is_array() == false) {
        KRATOS_ERROR << "SetArrayItem only makes sense if the value if of Array type" << std::endl;
    } else {
        KRATOS_ERROR_IF(Index >= mpValue->size()) << "Index exceeds array size. Index value is : " << Index <<std::endl;
        (*mpValue)[Index] = *rOtherArrayItem.mpValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::AddEmptyArray(const std::string& rEntry)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    if(mpValue->find(rEntry) == mpValue->end()) {
        nlohmann::json j_array(nlohmann::json::value_t::array);
        (*mpValue)[rEntry] = j_array;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const double Value)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    nlohmann::json j_number_float(nlohmann::json::value_t::number_float);
    j_number_float = Value;
    mpValue->push_back(j_number_float);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const int Value)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    nlohmann::json j_number_integer(nlohmann::json::value_t::number_integer);
    j_number_integer = Value;
    mpValue->push_back(j_number_integer);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const bool Value)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    nlohmann::json j_boolean(nlohmann::json::value_t::boolean);
    j_boolean = Value;
    mpValue->push_back(j_boolean);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const std::string& rValue)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    nlohmann::json j_string(nlohmann::json::value_t::string);
    j_string = rValue;
    mpValue->push_back(j_string);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const Vector& rValue)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    const SizeType size = rValue.size();

    nlohmann::json j_array(0.0, size);

    for (IndexType i = 0; i < size; ++i) {
        j_array = rValue[i];
    }

    mpValue->push_back(j_array);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const Matrix& rValue)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    const SizeType nrows = rValue.size1();
    const SizeType ncols = rValue.size2();

    nlohmann::json j_col_array(0.0, ncols);
    nlohmann::json j_row_array(0.0, nrows);

    for (IndexType i = 0; i < nrows; ++i) {
        for (IndexType j = 0; j < ncols; ++j) {
            j_col_array[j] = rValue(i, j);
        }

        j_row_array[i] = j_col_array;
    }

    mpValue->push_back(j_row_array);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::Append(const Parameters& rValue)
{
    KRATOS_ERROR_IF_NOT(mpValue->is_array()) << "It must be an Array parameter to append" << std::endl;
    nlohmann::json j_object = nlohmann::json( nlohmann::json::parse( rValue.WriteJsonString() ) );
    mpValue->push_back(j_object);
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::RecursivelyFindValue(
    const nlohmann::json& rBaseValue,
    const nlohmann::json& rValueToFind
    ) const
{
    for (auto itr = rBaseValue.begin(); itr != rBaseValue.end(); ++itr) {
        if (&(itr.value()) == &rValueToFind) {
            KRATOS_INFO("Parameters") << "Base = " << PrettyPrintJsonString() << std::endl
                        << "Problematic var name " << itr.key() << " value " << itr.value() << std::endl;
        } else {
            if (itr->is_object()) RecursivelyFindValue(itr.value(), rValueToFind);
            //TODO: it could be an array
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::IsEquivalentTo(Parameters& rParameters)
{
    for (auto itr = this->mpValue->begin(); itr != this->mpValue->end(); ++itr) {
        const std::string& item_name = itr.key();

        bool found = false;

        for (auto itr_ref = rParameters.mpValue->begin(); itr_ref != rParameters.mpValue->end(); ++itr_ref) {
            if (item_name == itr_ref.key()) {
                found = true;
                Parameters subobject = (*this)[item_name];
                Parameters reference_subobject = rParameters[item_name];

                if (itr->is_object()) {
                    if (!subobject.IsEquivalentTo(reference_subobject))
                        return false;
                } else {
                    if (itr.value() != itr_ref.value())
                        return false;
                }
                break;
            }
        }

        if (!found)
            return false;
    }

    // Reverse check: the rParameters can contain fields that are missing in the object
    for (auto itr = rParameters.mpValue->begin();  itr != rParameters.mpValue->end(); ++itr) {
        const std::string& item_name = itr.key();

        bool found = false;

        for (auto itr_ref = this->mpValue->begin(); itr_ref != this->mpValue->end(); ++itr_ref) {
            if (item_name == itr_ref.key()) {
                found = true;
                // No need to check the values here, if they were found in the previous loop, values were checked there
                break;
            }
        }

        if (!found)
            return false;
    }

    return true;
}

/***********************************************************************************/
/***********************************************************************************/

bool Parameters::HasSameKeysAndTypeOfValuesAs(Parameters& rParameters)
{
    for (auto itr = this->mpValue->begin(); itr != this->mpValue->end(); ++itr) {
        const std::string& item_name = itr.key();

        bool found = false;

        for (auto itr_ref = rParameters.mpValue->begin(); itr_ref != rParameters.mpValue->end(); ++itr_ref) {
            if (item_name == itr_ref.key()) {
                found = true;
                Parameters subobject = (*this)[item_name];
                Parameters reference_subobject = rParameters[item_name];

                if (itr->is_object()) {
                    if (!subobject.HasSameKeysAndTypeOfValuesAs(reference_subobject))
                        return false;
                } else {
                    if (itr.value().type() != itr_ref.value().type()) {
                        return false;
                    }
                }
                break;
            }
        }

        if (!found)
            return false;
    }

    // Reverse check: the rParameters can contain fields that are missing in the object
    for (auto itr = rParameters.mpValue->begin(); itr != rParameters.mpValue->end(); ++itr) {
        const std::string& item_name = itr.key();

        bool found = false;

        for (auto itr_ref =  this->mpValue->begin(); itr_ref != this->mpValue->end(); ++itr_ref) {
            if (item_name == itr_ref.key()) {
                found = true;
                // No need to check the types here, if they were found in the previous loop, types were checked there
                break;
            }
        }

        if (!found)
            return false;
    }

    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::ValidateAndAssignDefaults(Parameters& rDefaultParameters)
{
    KRATOS_TRY

    // First verifies that all the entries in the current parameters have a correspondance in the rDefaultParameters.
    // If it is not the case throw an error
    for (auto itr = this->mpValue->begin(); itr != this->mpValue->end(); ++itr) {
        const std::string& r_item_name = itr.key();
        if(!rDefaultParameters.Has(r_item_name) ) {
            std::stringstream msg;
            msg << "The item with name \"" << r_item_name << "\" is present in this Parameters but NOT in the default values" << std::endl;
            msg << "Hence Validation fails" << std::endl;
            msg << "Parameters being validated are : " << std::endl;
            msg << this->PrettyPrintJsonString() << std::endl;
            msg << "Defaults against which the current parameters are validated are :" << std::endl;
            msg << rDefaultParameters.PrettyPrintJsonString() << std::endl;
            KRATOS_ERROR << msg.str() << std::endl;
        }

        bool type_coincides = false;
        auto value_defaults = (rDefaultParameters[r_item_name]).GetUnderlyingStorage();
        if(itr->is_number() && value_defaults->is_number()) type_coincides = true;
//             if(itr->is_number_integer() && value_defaults->is_number_integer()) type_coincides = true;
//             if(itr->is_number_float() && value_defaults->is_number_float()) type_coincides = true;
        if(itr->is_boolean() && value_defaults->is_boolean()) type_coincides = true;
        if(itr->is_null() && value_defaults->is_null()) type_coincides = true;
        if(itr->is_array() && value_defaults->is_array()) type_coincides = true;
        if(itr->is_string() && value_defaults->is_string()) type_coincides = true;
        if(itr->is_object() && value_defaults->is_object()) type_coincides = true;

        if(type_coincides == false) {
            std::stringstream msg;
            msg << "******************************************************************************************************" << std::endl;
            msg << "The item with name :\"" << r_item_name << "\" does not have the same type as the corresponding one in the default values" << std::endl;
            msg << "******************************************************************************************************" << std::endl;
            msg << "Parameters being validated are : " << std::endl;
            msg << this->PrettyPrintJsonString() << std::endl;
            msg << "Defaults against which the current parameters are validated are :" << std::endl;
            msg << rDefaultParameters.PrettyPrintJsonString() << std::endl;
            KRATOS_ERROR << msg.str() << std::endl;
        }

    }

    // Now iterate over all the rDefaultParameters. In the case a default value is not assigned in the current Parameters add an item copying its value
    if (rDefaultParameters.IsSubParameter()) {
        for (auto itr = rDefaultParameters.mpValue->begin(); itr != rDefaultParameters.mpValue->end(); ++itr) {
            const std::string& r_item_name = itr.key();
            if(mpValue->find(r_item_name) == mpValue->end()) {
                (*mpValue)[r_item_name] = itr.value();
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::RecursivelyValidateAndAssignDefaults(Parameters& rDefaultParameters)
{
    KRATOS_TRY

    // First verifies that all the entries in the current parameters have a correspondance in the rDefaultParameters.
    // If it is not the case throw an error
    for (auto itr = this->mpValue->cbegin(); itr != this->mpValue->cend(); ++itr) {
        const std::string& r_item_name = itr.key();

        if(!rDefaultParameters.Has(r_item_name) ) {
            std::stringstream msg;
            msg << "The item with name \"" << r_item_name << "\" is present in this Parameters but NOT in the default values" << std::endl;
            msg << "Hence Validation fails" << std::endl;
            msg << "Parameters being validated are : " << std::endl;
            msg << this->PrettyPrintJsonString() << std::endl;
            msg << "Defaults against which the current parameters are validated are :" << std::endl;
            msg << rDefaultParameters.PrettyPrintJsonString() << std::endl;
            KRATOS_ERROR << msg.str() << std::endl;
        }

        bool type_coincides = false;
        auto value_defaults = (rDefaultParameters[r_item_name]).GetUnderlyingStorage();
        if(itr->is_number() && value_defaults->is_number()) type_coincides = true;
//             if(itr->is_number_integer() && value_defaults->is_number_integer()) type_coincides = true;
//             if(itr->is_number_float() && value_defaults->is_number_float()) type_coincides = true;
        if(itr->is_boolean() && value_defaults->is_boolean()) type_coincides = true;
        if(itr->is_null() && value_defaults->is_null()) type_coincides = true;
        if(itr->is_array() && value_defaults->is_array()) type_coincides = true;
        if(itr->is_string() && value_defaults->is_string()) type_coincides = true;
        if(itr->is_object() && value_defaults->is_object()) type_coincides = true;

        if(type_coincides == false) {
            std::stringstream msg;
            msg << "The item with name :\"" << r_item_name << "\" does not have the same type as the corresponding one in the default values" << std::endl;
            msg << "Parameters being validated are : " << std::endl;
            msg << this->PrettyPrintJsonString() << std::endl;
            msg << "Defaults against which the current parameters are validated are :" << std::endl;
            msg << rDefaultParameters.PrettyPrintJsonString() << std::endl;
            KRATOS_ERROR << msg.str() << std::endl;
        }

        // Now walk the tree recursively
        if(itr->is_object()) {
            Parameters subobject = (*this)[r_item_name];
            Parameters defaults_subobject = rDefaultParameters[r_item_name];
            subobject.RecursivelyValidateAndAssignDefaults(defaults_subobject);
        }
    }

    // Now iterate over all the rDefaultParameters. In the case a default value is not assigned in the current Parameters add an item copying its value
    if (rDefaultParameters.IsSubParameter()) {
        for (auto itr = rDefaultParameters.mpValue->begin(); itr != rDefaultParameters.mpValue->end(); ++itr) {
            const std::string& r_item_name = itr.key();

            if(mpValue->find(r_item_name) == mpValue->end()) {
                (*mpValue)[r_item_name] = itr.value();
            }

            // Now walk the tree recursively
            if(itr->is_object()) {
                Parameters subobject = (*this)[r_item_name];
                Parameters defaults_subobject = rDefaultParameters[r_item_name];

                subobject.RecursivelyValidateAndAssignDefaults(defaults_subobject);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::save(Serializer& rSerializer) const
{
    rSerializer.save("Data", this->WriteJsonString());
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::load(Serializer& rSerializer)
{
    std::string parameters_data;
    rSerializer.load("Data", parameters_data);
    *this = Parameters(parameters_data);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters(nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot)
    : mpValue(pValue),
      mpRoot(pRoot)
{}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters(json_iterator itValue, nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot)
    : mpValue(nullptr),
        mpRoot(pRoot)
{
    if (itValue != pValue->end())
        mpValue = &(*itValue);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters::Parameters(json_const_iterator itValue, nlohmann::json* pValue, Kratos::shared_ptr<nlohmann::json> pRoot)
    : mpValue(nullptr),
        mpRoot(pRoot)
{
    if (itValue != pValue->cend())
        mpValue = const_cast<nlohmann::json*>(&(*itValue));
}

/***********************************************************************************/
/***********************************************************************************/

nlohmann::json* Parameters::GetUnderlyingStorage()
{
    return mpValue;
}

/***********************************************************************************/
/***********************************************************************************/

nlohmann::json* Parameters::GetUnderlyingStorage() const
{
    return mpValue;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetUnderlyingSotrage(nlohmann::json* pNewValue)
{
    mpValue = pNewValue;
}

/***********************************************************************************/
/***********************************************************************************/

Kratos::shared_ptr<nlohmann::json> Parameters::GetUnderlyingRootStorage()
{
    return mpRoot;
}

/***********************************************************************************/
/***********************************************************************************/

Kratos::shared_ptr<nlohmann::json> Parameters::GetUnderlyingRootStorage() const
{
    return mpRoot;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::SetUnderlyingRootStorage(Kratos::shared_ptr<nlohmann::json> pNewValue)
{
    mpRoot = pNewValue;
}

/***********************************************************************************/
/***********************************************************************************/

void Parameters::InternalSetValue(const Parameters& rOtherValue)
{
    delete[] mpValue;
    mpValue = new nlohmann::json( nlohmann::json::parse( rOtherValue.WriteJsonString()));
}

}  // namespace Kratos.
