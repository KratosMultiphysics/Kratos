//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_KRATOS_PARAMETERS_H_INCLUDED)
#define KRATOS_KRATOS_PARAMETERS_H_INCLUDED

// System includes

#include <iostream>
#include <sstream>
#include <string>
#include <utility>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"

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
class Parameters {
private:
  ///@name Nested clases
  ///@{
  class iterator_adaptor
      : public std::iterator<std::forward_iterator_tag, Parameters> {
    using value_iterator = rapidjson::Value::MemberIterator;
    value_iterator mValueIterator;
    std::unique_ptr<Parameters> mpParameters;

  public:
    iterator_adaptor(value_iterator it,
                     Kratos::shared_ptr<rapidjson::Document> pdoc)
        : mValueIterator(it), mpParameters(new Parameters(&it->value, pdoc)) {}
    iterator_adaptor(const iterator_adaptor &it)
        : mValueIterator(it.mValueIterator),
          mpParameters(new Parameters(*(it.mpParameters))) {}
    iterator_adaptor &operator++() {
      mValueIterator++;
      return *this;
    }
    iterator_adaptor operator++(int) {
      iterator_adaptor tmp(*this);
      operator++();
      return tmp;
    }
    bool operator==(const iterator_adaptor &rhs) const {
      return mValueIterator == rhs.mValueIterator;
    }
    bool operator!=(const iterator_adaptor &rhs) const {
      return mValueIterator != rhs.mValueIterator;
    }
    Parameters &operator*() const {
      mpParameters->SetUnderlyingSotrage(&mValueIterator->value);
      return *mpParameters;
    }
    Parameters *operator->() const {
      mpParameters->SetUnderlyingSotrage(&mValueIterator->value);
      return mpParameters.get();
    }
    value_iterator &base() { return mValueIterator; }
    value_iterator const &base() const { return mValueIterator; }
    std::string name() { return mValueIterator->name.GetString(); }
  };

  class const_iterator_adaptor
      : public std::iterator<std::forward_iterator_tag, Parameters> {
    using value_iterator = rapidjson::Value::ConstMemberIterator;
    value_iterator mValueIterator;
    std::unique_ptr<Parameters> mpParameters;

  public:
    const_iterator_adaptor(value_iterator it,
                           Kratos::shared_ptr<rapidjson::Document> pdoc)
        : mValueIterator(it),
          mpParameters(new Parameters(
              const_cast<rapidjson::Value *>(&it->value), pdoc)) {}
    const_iterator_adaptor(const const_iterator_adaptor &it)
        : mValueIterator(it.mValueIterator),
          mpParameters(new Parameters(*(it.mpParameters))) {}
    const_iterator_adaptor &operator++() {
      mValueIterator++;
      return *this;
    }
    const_iterator_adaptor operator++(int) {
      const_iterator_adaptor tmp(*this);
      operator++();
      return tmp;
    }
    bool operator==(const const_iterator_adaptor &rhs) const {
      return mValueIterator == rhs.mValueIterator;
    }
    bool operator!=(const const_iterator_adaptor &rhs) const {
      return mValueIterator != rhs.mValueIterator;
    }
    const Parameters &operator*() const {
      mpParameters->SetUnderlyingSotrage(
          const_cast<rapidjson::Value *>(&mValueIterator->value));
      return *mpParameters;
    }
    const Parameters *operator->() const {
      mpParameters->SetUnderlyingSotrage(
          const_cast<rapidjson::Value *>(&mValueIterator->value));
      return mpParameters.get();
    }
    value_iterator &base() { return mValueIterator; }
    value_iterator const &base() const { return mValueIterator; }
    std::string name() { return mValueIterator->name.GetString(); }
  };

  ///@}

public:
  KRATOS_CLASS_POINTER_DEFINITION(Parameters);

  using iterator = iterator_adaptor;
  using const_iterator = const_iterator_adaptor;

  Parameters(const std::string& json_string = "{}") {

    mpdoc = Kratos::make_shared<rapidjson::Document>();
    rapidjson::ParseResult ok = mpdoc->Parse<0>(json_string.c_str());

    if (!ok) {
      std::stringstream msg;
      msg << rapidjson::GetParseError_En(ok.Code())
          << " offset of the error from the beginning of the string = "
          << ok.Offset() << std::endl;
      msg << "a much more explicative error message can be obtained by "
             "analysing the input string with an online analyzer such for "
             "example json lint"
          << std::endl;
      msg << "the value of the string that was attempted to parse is :"
          << std::endl
          << std::endl;
      msg << json_string;
      KRATOS_ERROR << "error found in parsing the json_string, the value of "
                      "the json string was: \n"
                   << msg.str() << std::endl;
    }

    mpvalue = (mpdoc.get());
  }

  /// Assignment operator.
  Parameters &operator=(Parameters const &rOther) {
    mpvalue->CopyFrom(*(rOther.GetUnderlyingStorage()), mpdoc->GetAllocator());

    return *this;
  }

  Parameters& operator=(Parameters&& rOther) {
      Reset();
      swap(rOther);
      return *this;
  }

  /// Copy constructor.
  Parameters(Parameters const &rOther) {
    mpdoc = rOther.mpdoc;
    mpvalue = rOther.mpvalue;
  }

  Parameters(Parameters&& rOther) {
      Reset();
      swap(rOther);
  }

  // generates a clone of the current document
  Parameters Clone() {
    auto pnew_cloned_doc =
        Kratos::make_shared<rapidjson::Document>();
    rapidjson::ParseResult ok =
        pnew_cloned_doc->Parse<0>(this->WriteJsonString().c_str());
    if (!ok) {
      std::stringstream msg;
      msg << rapidjson::GetParseError_En(ok.Code())
          << " offset of the error from the beginning of the string = "
          << ok.Offset() << std::endl;
      msg << "a much more explicative error message can be obtained by "
             "analysing the input string "
          << std::endl;
      msg << "with an online analyzer such for example json lint" << std::endl;
      msg << "the value of the string that was attempted to parse is :"
          << std::endl
          << std::endl;
      msg << this->WriteJsonString();
      KRATOS_ERROR << "error found in parsing the json_string, the value of "
                      "the json string was: \n"
                   << msg.str() << std::endl;
    }
    return Parameters(pnew_cloned_doc.get(), pnew_cloned_doc);
  }

  /// Destructor.
  virtual ~Parameters() {}

  const std::string WriteJsonString() const {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    mpvalue->Accept(writer);
    return buffer.GetString();
  }

  const std::string PrettyPrintJsonString() const {
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    mpvalue->Accept(writer);
    return buffer.GetString();
  }

  //*******************************************************************************************************
  Parameters GetValue(const std::string& entry) {
    KRATOS_ERROR_IF_NOT(this->Has(entry))
        << "--------- ERROR : --------- getting a value that does not exist. "
           "entry string : "
        << entry << std::endl;
    rapidjson::Value *pvalue = &((*mpvalue)[entry.c_str()]);

    return Parameters(pvalue, mpdoc);
  }
  Parameters operator[](const std::string& entry) {
    return this->GetValue(entry);
  }
  void SetValue(const std::string& entry, const Parameters &other_value) {
    KRATOS_ERROR_IF_NOT(this->Has(entry))
        << "value must exist to be set. Use AddValue instead" << std::endl;
    Parameters tmp(&(*mpvalue)[entry.c_str()], mpdoc);
    tmp.InternalSetValue(other_value);
  }

  void AddValue(const std::string& entry, const Parameters &other_value) {
    if (this->Has(entry) == false) {
      rapidjson::Value tmp;
      tmp.CopyFrom(*(other_value.GetUnderlyingStorage()),
                   mpdoc->GetAllocator()); // this will be moved away
      rapidjson::Value name(entry.c_str(),
                            mpdoc->GetAllocator()); // rhis will be moved away
      this->mpvalue->AddMember(name, tmp, mpdoc->GetAllocator());
    }
  }
  Parameters AddEmptyValue(const std::string& entry) {
    if (this->Has(entry) == false) {
      rapidjson::Value tmp;
      rapidjson::Value name(entry.c_str(),
                            mpdoc->GetAllocator()); // rhis will be moved away
      this->mpvalue->AddMember(name, tmp, mpdoc->GetAllocator());
    }
    return this->GetValue(entry);
  }

  bool RemoveValue(const std::string& entry) {
    return mpvalue->RemoveMember(entry.c_str());
  }

  //*******************************************************************************************************
  bool Has(const std::string& entry) const {
    return mpvalue->HasMember(entry.c_str());
  }

  bool IsNull() const { return mpvalue->IsNull(); }
  bool IsNumber() const { return mpvalue->IsNumber(); }
  bool IsDouble() const { return mpvalue->IsDouble(); }
  bool IsInt() const { return mpvalue->IsInt(); }
  bool IsBool() const { return mpvalue->IsBool(); }
  bool IsString() const { return mpvalue->IsString(); }
  bool IsArray() const { return mpvalue->IsArray(); }
  bool IsVector() const {
    if (!mpvalue->IsArray())
      return false;

    for (unsigned int i = 0; i < mpvalue->Size(); ++i) {
      if (!(*mpvalue)[i].IsNumber())
        return false;
    }
    return true; // All entries are numbers or Vector is empty
  }
  bool IsMatrix() const {
    if (!mpvalue->IsArray()) // mpvalue != [ ... ]
      return false;

    unsigned int nrows = mpvalue->Size();
    if (nrows == 0) // mpvalue is an empty array/vector => "[]"
      return false;

    for (unsigned int i = 0; i < nrows; ++i) {
      auto &row_i = (*mpvalue)[i];
      if (!row_i.IsArray())
        return false;

      unsigned int ncols = row_i.Size();
      if (ncols !=
          (*mpvalue)[0].Size()) // compare number of columns to first row
        return false;           // number of columns is not consistent

      for (unsigned int j = 0; j < ncols; ++j) // Check all values in column
      {
        if (!row_i[j].IsNumber())
          return false;
      }
    }
    return true; // All entries are numbers or Matrix is empty ([[]] or
                 // [[],[],[],...])
  }
  bool IsSubParameter() const { return mpvalue->IsObject(); }

  double GetDouble() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsNumber()) << "argument must be a number"
                                             << std::endl;
    return mpvalue->GetDouble();
  }
  int GetInt() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsNumber()) << "argument must be a number"
                                             << std::endl;
    return mpvalue->GetInt();
  }
  bool GetBool() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsBool()) << "argument must be a bool"
                                           << std::endl;
    return mpvalue->GetBool();
  }
  std::string GetString() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsString()) << "argument must be a string"
                                             << std::endl;
    return mpvalue->GetString();
  }
  Vector GetVector() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "argument must be a Vector (a json list)" << std::endl;

    const unsigned int size = mpvalue->Size();

    Vector V(size);

    for (unsigned int i = 0; i < size; ++i) {
      KRATOS_ERROR_IF_NOT((*mpvalue)[i].IsNumber())
          << "Entry " << i << " is not a number!" << std::endl;
      V(i) = (*mpvalue)[i].GetDouble();
    }

    return V;
  }
  Matrix GetMatrix() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "argument must be a Matrix (a json list of lists)" << std::endl;

    const unsigned int nrows = mpvalue->Size();
    KRATOS_ERROR_IF(nrows == 0)
        << "argument must be a Matrix (a json list of lists)" << std::endl;

    unsigned int ncols = 0;
    if ((*mpvalue)[0].IsArray())
      ncols = (*mpvalue)[0].Size();

    Matrix A(nrows, ncols);

    for (unsigned int i = 0; i < nrows; ++i) {
      auto &row_i = (*mpvalue)[i];
      KRATOS_ERROR_IF_NOT(row_i.IsArray()) << "not an array on row " << i
                                           << std::endl;
      KRATOS_ERROR_IF_NOT(row_i.Size() == ncols) << "wrong size of row " << i
                                                 << std::endl;
      for (unsigned int j = 0; j < ncols; ++j) {
        KRATOS_ERROR_IF_NOT((row_i)[j].IsNumber())
            << "Entry (" << i << "," << j << ") is not a number!" << std::endl;
        A(i, j) = (row_i)[j].GetDouble();
      }
    }

    return A;
  }

  void SetDouble(const double value) { mpvalue->SetDouble(value); }
  void SetInt(const int value) { mpvalue->SetInt(value); }
  void SetBool(const bool value) { mpvalue->SetBool(value); }
  void SetString(const std::string& value) {
    rapidjson::Value tmp(value.c_str(), mpdoc->GetAllocator());
    *mpvalue = tmp;
  }
  void SetVector(const Vector &vec) {
    const unsigned int size = vec.size();

    mpvalue->SetArray();
    mpvalue->Reserve(size, mpdoc->GetAllocator());

    for (unsigned int i = 0; i < size; ++i) {
      mpvalue->PushBack(vec[i], mpdoc->GetAllocator());
    }
  }
  void SetMatrix(const Matrix &mat) {
    const unsigned int nrows = mat.size1();
    const unsigned int ncols = mat.size2();

    mpvalue->SetArray();
    mpvalue->Reserve(nrows, mpdoc->GetAllocator());

    for (unsigned int i = 0; i < nrows; ++i) {
      mpvalue->PushBack(0, mpdoc->GetAllocator()); // Pushing back a default
                                                   // element to allocate memory
      (*mpvalue)[i].SetArray(); // change that default element to an array
      (*mpvalue)[i].Reserve(ncols, mpdoc->GetAllocator());

      for (unsigned int j = 0; j < ncols; ++j) {
        (*mpvalue)[i].PushBack(mat(i, j), mpdoc->GetAllocator());
      }
    }
  }

  iterator begin() { return iterator(this->mpvalue->MemberBegin(), mpdoc); }

  iterator end() { return iterator(this->mpvalue->MemberEnd(), mpdoc); }

  const_iterator begin() const {
    return const_iterator(this->mpvalue->MemberBegin(), mpdoc);
  }

  const_iterator end() const {
    return const_iterator(this->mpvalue->MemberEnd(), mpdoc);
  }

  //*******************************************************************************************************
  // methods for array
  unsigned int size() const {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "size can only be queried if the value if of Array type"
        << std::endl;
    return mpvalue->Size();
  }

  void swap(Parameters& rOther) noexcept
  {
      std::swap(mpdoc, rOther.mpdoc);
      std::swap(mpvalue, rOther.mpvalue);
  }

  void Reset() noexcept
  {
      Parameters p;
      swap(p);
  }

  Parameters GetArrayItem(unsigned int index) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "GetArrayItem only makes sense if the value if of Array type"
        << std::endl;
    KRATOS_ERROR_IF(index >= mpvalue->Size())
        << "index exceeds array size. Index value is : " << index << std::endl;
    return Parameters(&(*mpvalue)[index], mpdoc);
  }

  void SetArrayItem(unsigned int index, const Parameters &other_array_item) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "SetArrayItem only makes sense if the value if of Array type"
        << std::endl;
    KRATOS_ERROR_IF(index >= mpdoc->Size())
        << "index exceeds array size. Index value is : " << index << std::endl;

#if RAPIDJSON_HAS_CXX11_RVALUE_REFS
    (*mpvalue)[index] = rapidjson::Value(
        *other_array_item.GetUnderlyingStorage(), mpdoc->GetAllocator());
#else
    (*mpvalue)[index].CopyFrom(*other_array_item.GetUnderlyingStorage(),
                               mpdoc->GetAllocator());
#endif
  }

  void AddEmptyArray(const std::string &entry) {
    if (this->Has(entry) == false) {
      rapidjson::Value tmp;
      tmp.SetArray();
      rapidjson::Value name(entry.c_str(),
                            mpdoc->GetAllocator()); // rhis will be moved away
      this->mpvalue->AddMember(name, tmp, mpdoc->GetAllocator());
    }
  }

  void Append(const double value) {
    rapidjson::Value tmp_value;
    tmp_value.SetDouble(value);
    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }
  void Append(const int value) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value;
    tmp_value.SetInt(value);
    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  void Append(const bool value) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value;
    tmp_value.SetBool(value);
    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  void Append(const std::string& value) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value(value.c_str(), mpdoc->GetAllocator());
    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  void Append(const Vector &vec) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value;

    const unsigned int size = vec.size();

    tmp_value.SetArray();
    tmp_value.Reserve(size, mpdoc->GetAllocator());

    for (unsigned int i = 0; i < size; ++i) {
      tmp_value.PushBack(vec[i], mpdoc->GetAllocator());
    }

    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  void Append(const Matrix &mat) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value;

    const unsigned int nrows = mat.size1();
    const unsigned int ncols = mat.size2();

    tmp_value.SetArray();
    tmp_value.Reserve(nrows, mpdoc->GetAllocator());

    for (unsigned int i = 0; i < nrows; ++i) {
      tmp_value.PushBack(0, mpdoc->GetAllocator()); // Pushing back a default
                                                    // element to allocate
                                                    // memory
      tmp_value[i].SetArray(); // change that default element to an array
      tmp_value[i].Reserve(ncols, mpdoc->GetAllocator());

      for (unsigned int j = 0; j < ncols; ++j) {
        tmp_value[i].PushBack(mat(i, j), mpdoc->GetAllocator());
      }
    }

    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  void Append(const Parameters &object) {
    KRATOS_ERROR_IF_NOT(mpvalue->IsArray())
        << "it must be an Array parameter to append" << std::endl;
    rapidjson::Value tmp_value;
    tmp_value.CopyFrom(*(object.GetUnderlyingStorage()), mpdoc->GetAllocator());
    mpvalue->PushBack(tmp_value, mpdoc->GetAllocator());
  }

  Parameters operator[](unsigned int index) {
    return this->GetArrayItem(index);
  }

  /**This function is designed to verify that the parameters under testing match
   * the
   * form prescribed by the defaults.
   * If the parameters contain values that do not appear in the defaults, an
   * error is thrown,
   * whereas if a parameter is found in the defaults but not in the Parameters
   * been tested,
   * it is copied to the parameters.
   *
   * this version of the function only walks one level, without descending in
   * the branches
   */
  void ValidateAndAssignDefaults(Parameters &defaults) {
    KRATOS_TRY

    // first verifies that all the enries in the current parameters
    // have a correspondance in the defaults.
    // if it is not the case throw an error
    for (rapidjson::Value::ConstMemberIterator itr =
             this->mpvalue->MemberBegin();
         itr != this->mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      if (!defaults.Has(item_name)) {
        std::stringstream msg;
        msg << "***************************************************************"
               "***************************************"
            << std::endl;
        msg << "the item with name \"" << item_name
            << "\" is present in this Parameters but NOT in the default values"
            << std::endl;
        msg << "***************************************************************"
               "***************************************"
            << std::endl;
        msg << "hence Validation fails" << std::endl;
        msg << "parameters being validated are : " << std::endl;
        msg << this->PrettyPrintJsonString() << std::endl;
        msg << "defaults against which the current parameters are validated "
               "are :"
            << std::endl;
        msg << defaults.PrettyPrintJsonString() << std::endl;
        KRATOS_ERROR << msg.str() << std::endl;
      }

      bool type_coincides = false;
      rapidjson::Value *value_defaults =
          (defaults[item_name.c_str()]).GetUnderlyingStorage();
      if (itr->value.IsInt() && value_defaults->IsNumber())
        type_coincides = true;
      if (itr->value.IsBool() && value_defaults->IsBool())
        type_coincides = true;
      if (itr->value.IsDouble() && value_defaults->IsDouble())
        type_coincides = true;
      if (itr->value.IsArray() && value_defaults->IsArray())
        type_coincides = true;
      if (itr->value.IsString() && value_defaults->IsString())
        type_coincides = true;
      if (itr->value.IsObject() && value_defaults->IsObject())
        type_coincides = true;
      if (itr->value.IsNull() && value_defaults->IsNull())
        type_coincides = true;

      if (type_coincides == false) {
        std::stringstream msg;
        msg << "***************************************************************"
               "***************************************"
            << std::endl;
        msg << "the item with name :\"" << item_name
            << "\" does not have the same type as the corresponding one in the "
               "default values"
            << std::endl;
        msg << "***************************************************************"
               "***************************************"
            << std::endl;
        msg << "parameters being validated are : " << std::endl;
        msg << this->PrettyPrintJsonString() << std::endl;
        msg << "defaults against which the current parameters are validated "
               "are :"
            << std::endl;
        msg << defaults.PrettyPrintJsonString() << std::endl;
        KRATOS_ERROR << msg.str() << std::endl;
      }
    }

    // now iterate over all the defaults. In the case a default value is not
    // assigned in the current Parameters
    // add an item copying its value
    if (defaults.IsSubParameter()) {
      for (rapidjson::Value::MemberIterator itr =
               defaults.mpvalue->MemberBegin();
           itr != defaults.mpvalue->MemberEnd(); ++itr) {
        std::string item_name = itr->name.GetString();
        if (!this->Has(item_name)) {
          rapidjson::Value *pvalue = &itr->value;

          this->AddValue(item_name, Parameters(pvalue, defaults.mpdoc));
        }
      }
    }

    KRATOS_CATCH("")
  }

  /**This function is designed to verify that the parameters under testing match
   * the
   * form prescribed by the defaults.
   * If the parameters contain values that do not appear in the defaults, an
   * error is thrown,
   * whereas if a parameter is found in the defaults but not in the Parameters
   * been tested,
   * it is copied to the parameters.
   *
   * this version walks and validates the entire json tree below
   * the point at which the function is called
  */
  void RecursivelyValidateAndAssignDefaults(Parameters &defaults) {
    KRATOS_TRY

    // first verifies that all the enries in the current parameters
    // have a correspondance in the defaults.
    // if it is not the case throw an error
    for (rapidjson::Value::ConstMemberIterator itr =
             this->mpvalue->MemberBegin();
         itr != this->mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      if (!defaults.Has(item_name)) {
        std::stringstream msg;
        msg << "the item with name \"" << item_name
            << "\" is present in this Parameters but NOT in the default values"
            << std::endl;
        msg << "hence Validation fails" << std::endl;
        msg << "parameters being validated are : " << std::endl;
        msg << this->PrettyPrintJsonString() << std::endl;
        msg << "defaults against which the current parameters are validated "
               "are :"
            << std::endl;
        msg << defaults.PrettyPrintJsonString() << std::endl;
        KRATOS_ERROR << msg.str() << std::endl;
      }

      bool type_coincides = false;
      rapidjson::Value *value_defaults =
          (defaults[item_name.c_str()]).GetUnderlyingStorage();
      if (itr->value.IsInt() && value_defaults->IsInt())
        type_coincides = true;
      if (itr->value.IsBool() && value_defaults->IsBool())
        type_coincides = true;
      if (itr->value.IsDouble() && value_defaults->IsDouble())
        type_coincides = true;
      if (itr->value.IsArray() && value_defaults->IsArray())
        type_coincides = true;
      if (itr->value.IsString() && value_defaults->IsString())
        type_coincides = true;
      if (itr->value.IsObject() && value_defaults->IsObject())
        type_coincides = true;
      if (itr->value.IsNull() && value_defaults->IsNull())
        type_coincides = true;

      if (type_coincides == false) {
        std::stringstream msg;
        msg << "the item with name :\"" << item_name
            << "\" does not have the same type as the corresponding one in the "
               "default values"
            << std::endl;
        msg << "parameters being validated are : " << std::endl;
        msg << this->PrettyPrintJsonString() << std::endl;
        msg << "defaults against which the current parameters are validated "
               "are :"
            << std::endl;
        msg << defaults.PrettyPrintJsonString() << std::endl;
        KRATOS_ERROR << msg.str() << std::endl;
      }
      // now walk the tree recursively
      if (itr->value.IsObject()) {
        Parameters subobject = (*this)[item_name];
        Parameters defaults_subobject = defaults[item_name];
        subobject.RecursivelyValidateAndAssignDefaults(defaults_subobject);
      }
    }

    // now iterate over all the defaults. In the case a default value is not
    // assigned in the current Parameters
    // add an item copying its value
    if (defaults.IsSubParameter()) {
      for (rapidjson::Value::MemberIterator itr =
               defaults.mpvalue->MemberBegin();
           itr != defaults.mpvalue->MemberEnd(); ++itr) {
        std::string item_name = itr->name.GetString();
        if (!this->Has(item_name)) {
          rapidjson::Value *pvalue = &itr->value;

          this->AddValue(item_name, Parameters(pvalue, defaults.mpdoc));
        }

        // now walk the tree recursively
        if (itr->value.IsObject()) {
          Parameters subobject = (*this)[item_name];
          Parameters defaults_subobject = defaults[item_name];
          subobject.ValidateAndAssignDefaults(defaults_subobject);
        }
      }
    }

    KRATOS_CATCH("")
  }

  void RecursivelyFindValue(const rapidjson::Value &rbase_value,
                            const rapidjson::Value &rvalue_to_find) const {
    for (rapidjson::Value::ConstMemberIterator itr = rbase_value.MemberBegin();
         itr != rbase_value.MemberEnd(); ++itr) {
      if (&(itr->value) == &rvalue_to_find) {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        mpvalue->Accept(writer);
        std::cout << "base = " << buffer.GetString() << std::endl;
        std::cout << "problematic var name " << itr->name.GetString()
                  << " value " << itr->value.GetString() << std::endl;
      } else {
        if (itr->value.IsObject())
          RecursivelyFindValue(itr->value, rvalue_to_find);
        // TODO: it could be an array
      }
    }
  }

  bool IsEquivalentTo(Parameters &reference) {
    // Checks if the names and values are the same, no importance to the order.
    // Lists have to be ordered, though! Take into account that in Kratos some
    // physical vectors are represented with a list.

    for (rapidjson::Value::ConstMemberIterator itr =
             this->mpvalue->MemberBegin();
         itr != this->mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      bool found = false;

      for (rapidjson::Value::ConstMemberIterator itr_ref =
               reference.mpvalue->MemberBegin();
           itr_ref != reference.mpvalue->MemberEnd(); ++itr_ref) {
        if (item_name == itr_ref->name.GetString()) {
          found = true;
          Parameters subobject = (*this)[item_name];
          Parameters reference_subobject = reference[item_name];

          if (itr->value.IsObject()) {
            if (!subobject.IsEquivalentTo(reference_subobject))
              return false;
          } else {
            if (itr->value != itr_ref->value)
              return false;
          }
          break;
        }
      }

      if (!found)
        return false;
    }

    // reverse check: the reference can contain fields that are missing in the
    // object
    for (rapidjson::Value::ConstMemberIterator itr =
             reference.mpvalue->MemberBegin();
         itr != reference.mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      bool found = false;

      for (rapidjson::Value::ConstMemberIterator itr_ref =
               this->mpvalue->MemberBegin();
           itr_ref != this->mpvalue->MemberEnd(); ++itr_ref) {
        if (item_name == itr_ref->name.GetString()) {
          found = true;
          // no need to check the values here, if they were found in the
          // previous loop, values were checked there
          break;
        }
      }

      if (!found)
        return false;
    }

    return true;
  }

  bool HasSameKeysAndTypeOfValuesAs(Parameters &reference) {
    // Checks if the names and the type of values are the same, no importance to
    // the order.
    // Lists have to be ordered, though! Take into account that in Kratos some
    // physical vectors are represented with a list.

    for (rapidjson::Value::ConstMemberIterator itr =
             this->mpvalue->MemberBegin();
         itr != this->mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      bool found = false;

      for (rapidjson::Value::ConstMemberIterator itr_ref =
               reference.mpvalue->MemberBegin();
           itr_ref != reference.mpvalue->MemberEnd(); ++itr_ref) {
        if (item_name == itr_ref->name.GetString()) {
          found = true;
          Parameters subobject = (*this)[item_name];
          Parameters reference_subobject = reference[item_name];

          if (itr->value.IsObject()) {
            if (!subobject.HasSameKeysAndTypeOfValuesAs(reference_subobject))
              return false;
          } else {
            if (itr->value.GetType() != itr_ref->value.GetType()) {
              /*std::stringstream msg;
              msg << "The item with name :\"" << item_name << "\" does not have
              the same type as the corresponding one in the default values" <<
              std::endl;
              msg << "The Parameters being validated are : " << std::endl;
              msg << this->PrettyPrintJsonString() << std::endl;
              msg << "Defaults against which the current parameters are
              validated are :" << std::endl;
              msg << reference.PrettyPrintJsonString() << std::endl;
              KRATOS_THROW_ERROR(std::invalid_argument,"",msg.str());*/
              return false;
            }
          }
          break;
        }
      }

      if (!found)
        return false;
    }

    // reverse check: the reference can contain fields that are missing in the
    // object
    for (rapidjson::Value::ConstMemberIterator itr =
             reference.mpvalue->MemberBegin();
         itr != reference.mpvalue->MemberEnd(); ++itr) {
      std::string item_name = itr->name.GetString();

      bool found = false;

      for (rapidjson::Value::ConstMemberIterator itr_ref =
               this->mpvalue->MemberBegin();
           itr_ref != this->mpvalue->MemberEnd(); ++itr_ref) {
        if (item_name == itr_ref->name.GetString()) {
          found = true;
          // no need to check the types here, if they were found in the previous
          // loop, types were checked there
          break;
        }
      }

      if (!found)
        return false;
    }

    return true;
  }

  /// Turn back information as a string.
  virtual std::string Info() const { return this->PrettyPrintJsonString(); }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const {
    rOStream << "Parameters Object " << Info();
  }

  /// Print object's data.
  virtual void PrintData(std::ostream &rOStream) const {};

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const 
    {
        rSerializer.save("Data", this->WriteJsonString());
    }

    void load(Serializer& rSerializer) 
    {
        std::string parameters_data;
        rSerializer.load("Data", parameters_data);
        *this = Parameters(parameters_data);
    }


  // ATTENTION: please DO NOT use this constructor. It assumes rapidjson and
  // hence it should be considered as an implementation detail
  Parameters(rapidjson::Value *pvalue,
             Kratos::shared_ptr<rapidjson::Document> pdoc)
      : mpvalue(pvalue), mpdoc(pdoc) {}

  rapidjson::Value *mpvalue = nullptr;
  Kratos::shared_ptr<rapidjson::Document> mpdoc;

  // ATTENTION: please DO NOT use this method. It is a low level accessor, and
  // may change in the future
  rapidjson::Value *GetUnderlyingStorage() const { return mpvalue; }

  void SetUnderlyingSotrage(rapidjson::Value *pNewValue) {
    mpvalue = pNewValue;
  }

  void InternalSetValue(const Parameters &other_value) {
    mpvalue->CopyFrom(*(other_value.GetUnderlyingStorage()),
                      mpdoc->GetAllocator());
  }
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream, Parameters &rThis) {
  return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const Parameters &rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_KRATOS_PARAMETERS_H_INCLUDED  defined
