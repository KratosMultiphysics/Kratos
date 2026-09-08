//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <numeric>
#include <sstream>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"

// Include baes h
#include "combined_tensor_adaptor.h"

namespace Kratos {

namespace {
[[nodiscard]] DenseVector<unsigned int> GetModifiedTensorShape(
    const IndexType NumberOfDimensions,
    const DenseVector<unsigned int>& rCurrentTensorShape)
{
    auto modified_current_tensor_shape = DenseVector<unsigned int>(NumberOfDimensions, 1);
    for (IndexType i_dim = 0; i_dim < rCurrentTensorShape.size(); ++i_dim) {
        modified_current_tensor_shape[i_dim] = rCurrentTensorShape[i_dim];
    }
    return modified_current_tensor_shape;
}
} // namespace

template<class TDataType>
CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor(
    const TensorAdaptorVectorType& rTensorAdaptorVector,
    const unsigned int Axis,
    const bool PerformCollectDataRecursively,
    const bool PerformStoreDataRecursively,
    const bool Copy)
    : mPerformCollectDataRecursively(PerformCollectDataRecursively),
      mPerformStoreDataRecursively(PerformStoreDataRecursively),
      mAxis(Axis)
{
    KRATOS_TRY

    for (const auto& p_tensor_adaptor : rTensorAdaptorVector) {
        if (Copy) {
            mTensorAdaptors.push_back(p_tensor_adaptor->Clone());
        } else {
            mTensorAdaptors.push_back(p_tensor_adaptor);
        }
    }

    // It is not allowed to create combined tensor adaptors with empty list of tensor adaptors vector
    KRATOS_ERROR_IF(rTensorAdaptorVector.empty())
        << "Creating a CombinedTensorAdaptor with empty list of tensor adaptors is not allowed.";

    // looking for the tensor with the highest number of dimensions. This is done to
    // append all the other tensor dimensions with 1 at the end for the combination.
    DenseVector<unsigned int> tensor_shape{};
    for (const auto& p_tensor_adaptor : mTensorAdaptors) {
        const auto& r_current_shape = p_tensor_adaptor->Shape();
        if (tensor_shape.size() < r_current_shape.size()) {
            tensor_shape = r_current_shape;
        }
    }

    KRATOS_ERROR_IF(mAxis >= static_cast<int>(tensor_shape.size()))
        << "Axis must be lesser than the number of dimensions [ axis = " << mAxis
        << "tensor shape = " << tensor_shape << " ].\n";

    tensor_shape[Axis] = 0;

    // first check whether all the other axes in the dimensions except for the specified @p Axis
    // is having same number of components.
    std::for_each(rTensorAdaptorVector.begin(), rTensorAdaptorVector.end(), [&rTensorAdaptorVector, &tensor_shape, Axis](const typename BaseType::Pointer& pTensorAdaptor) {
        auto r_current_tensor_shape = pTensorAdaptor->Shape();

        KRATOS_ERROR_IF_NOT(r_current_tensor_shape.size() <= tensor_shape.size())
            << "The number of dimensions mismatch. Each tensors number of dimensions must be equivalent or smaller than the first tensor's number of dimensions."
            << "\n\ttensor_1 = " << rTensorAdaptorVector.front()
            << "\n\ttensor_2 = " << *pTensorAdaptor << " ].\n";

        // modifying the number of dimensions by adding 1 at the end.
        auto modified_current_tensor_shape = GetModifiedTensorShape(tensor_shape.size(), r_current_tensor_shape);

        for (IndexType i_dim = 0; i_dim < tensor_shape.size(); ++i_dim) {
            KRATOS_ERROR_IF(i_dim != Axis && modified_current_tensor_shape[i_dim] != tensor_shape[i_dim])
                << "Number of components in each dimension should be same in all the tensors except in the axis dimension [ checking dimension = "
                << i_dim << ", axis dimension = " << Axis
                << ", modified tensor shape of tensor_2 = " << modified_current_tensor_shape << " ]"
                << "\n\ttensor_1 = " << rTensorAdaptorVector.front()
                << "\n\ttensor_2 = " << *pTensorAdaptor << " ].\n";
        }

        tensor_shape[Axis] += modified_current_tensor_shape[Axis];
    });

    // now create the storage
    // here we create a storage with a nullptr for the container, because
    // the combined tensor adaptor will not be able to represent the underlying
    // container without incurring heavy additional cost.
    this->mpStorage = Kratos::make_shared<typename BaseType::Storage>(tensor_shape);

    this->mpContainer = std::nullopt;

    KRATOS_CATCH("");
}

template<class TDataType>
CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor(
    const TensorAdaptorVectorType& rTensorAdaptorVector,
    const bool PerformCollectDataRecursively,
    const bool PerformStoreDataRecursively,
    const bool Copy)
    : mPerformCollectDataRecursively(PerformCollectDataRecursively),
      mPerformStoreDataRecursively(PerformStoreDataRecursively),
      mAxis(-1) // represent that this combined tensor adaptor is used with ravel
{
    KRATOS_TRY

    for (const auto& p_tensor_adaptor : rTensorAdaptorVector) {
        if (Copy) {
            mTensorAdaptors.push_back(p_tensor_adaptor->Clone());
        } else {
            mTensorAdaptors.push_back(p_tensor_adaptor);
        }
    }

    // this combined tensor will have everything raveled. Therefore, it will only have one dimension.
    DenseVector<unsigned int> tensor_shape(1);

    tensor_shape[0] = std::accumulate(rTensorAdaptorVector.begin(), rTensorAdaptorVector.end(), 0u, [](const auto Value, const auto& pTensorAdaptor) {
        return Value + pTensorAdaptor->Size();
    });

    this->mpStorage = Kratos::make_shared<typename BaseType::Storage>(tensor_shape);

    this->mpContainer = std::nullopt;

    KRATOS_CATCH("");
}

template<class TDataType>
CombinedTensorAdaptor<TDataType>::CombinedTensorAdaptor(
    const CombinedTensorAdaptor& rOther,
    const bool PerformCollectDataRecursively,
    const bool PerformStoreDataRecursively,
    const bool Copy)
    : BaseType(rOther, Copy),
      mPerformCollectDataRecursively(PerformCollectDataRecursively),
      mPerformStoreDataRecursively(PerformStoreDataRecursively),
      mAxis(rOther.mAxis)
{
    KRATOS_TRY

    for (const auto& p_tensor_adaptor : rOther.mTensorAdaptors) {
        if (Copy) {
            mTensorAdaptors.push_back(p_tensor_adaptor->Clone());
        } else {
            mTensorAdaptors.push_back(p_tensor_adaptor);
        }
    }

    KRATOS_CATCH("");
}

template<class TDataType>
typename TensorAdaptor<TDataType>::Pointer CombinedTensorAdaptor<TDataType>::Clone() const
{
    return Kratos::make_shared<CombinedTensorAdaptor<TDataType>>(*this, this->mPerformCollectDataRecursively, this->mPerformStoreDataRecursively);
}

template<class TDataType>
void CombinedTensorAdaptor<TDataType>::Check() const
{
    KRATOS_TRY

    std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [](auto& pTensorAdaptor) {
        pTensorAdaptor->Check();
    });

    // now check the size compatibility of the combined tensor adaptor storage.
    if (mAxis == -1) {
        // this is the case everything is raveled.
        const auto ravelled_size = std::accumulate(
                                                    mTensorAdaptors.begin(),
                                                    mTensorAdaptors.end(),
                                                    IndexType{},
                                                    [](const auto& rValue, const auto& pTensorAdaptor) { return rValue + pTensorAdaptor->Size(); }
                                        );

        KRATOS_ERROR_IF_NOT(this->Size() == ravelled_size)
            << "Size mismatch [ tensor adaptor size = " << this->Size()
            << ", ravelled size of combined tensor adaptors = " << ravelled_size
            << ", combined tensor adaptor = " << *this << " ].\n";

    } else {
        const auto& combined_shape = this->Shape();
        IndexType axis_components = 0;
        for (const auto& p_sub_tensor_adaptor : mTensorAdaptors) {
            const auto& orig_sub_shape = p_sub_tensor_adaptor->Shape();
            DenseVector<unsigned int> sub_shape(combined_shape.size(), 1);
            std::copy(orig_sub_shape.begin(), orig_sub_shape.end(), sub_shape.begin());

            KRATOS_ERROR_IF_NOT(combined_shape.size() == sub_shape.size())
                << "Number of dimensions mismatch [ combined tensor adaptor shape = " << combined_shape << ", "
                << "sub tensor adaptor shape = " << sub_shape << ", combined tensor adaptor = " << *this << ", "
                << "sub tensor adaptor = " << *p_sub_tensor_adaptor << " ].\n";

            for (int i = 0; i < static_cast<int>(combined_shape.size()); ++i) {
                KRATOS_ERROR_IF(i != mAxis && combined_shape[i] != sub_shape[i])
                    << "Number of components in dimension = " << i << " mismatch [ combined tensor adaptor shape = " << combined_shape << ", "
                    << "sub tensor adaptor shape = " << sub_shape << ", combined tensor adaptor = " << *this << ", "
                    << "sub tensor adaptor = " << *p_sub_tensor_adaptor << " ].\n";
            }
            axis_components += sub_shape[mAxis];
        }

        KRATOS_ERROR_IF_NOT(axis_components == combined_shape[mAxis])
            << "Number of components in combining dimension = " << mAxis << " mismatch [ combined tensor adaptor shape = " << combined_shape << ", "
            << "combined tensor adaptor = " << *this << " ].\n";
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void CombinedTensorAdaptor<TDataType>::CollectData()
{
    KRATOS_TRY

    // if recursively collection is chose, first collect data to the sub TensorAdaptors
    if (mPerformCollectDataRecursively) {
        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [](auto& pTensorAdaptor) {
            pTensorAdaptor->CollectData();
        });
    }

    auto combined_span = this->ViewData();
    IndexType offset = 0;

    // now collect data to the CombinedTensorAdaptor
    if (mAxis == 0 || mAxis == -1) {
        // for axis 0, or ravel we just need to append items.
        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [&combined_span, &offset](const auto& pTensorAdaptor) {
            const auto& current_span = pTensorAdaptor->ViewData();
            IndexPartition<IndexType>(current_span.size()).for_each([&combined_span, &current_span, &offset](const auto Index) {
                combined_span[offset + Index] = current_span[Index];
            });
            offset += current_span.size();
        });
    } else {
        // this needs to append according to the axis
        const auto& tensor_shape = this->Shape();
        const auto skip_size = std::accumulate(tensor_shape.begin() + mAxis, tensor_shape.end(), 1u, std::multiplies<unsigned int>{});

        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [this, &tensor_shape, &combined_span, &offset, skip_size](const auto& pTensorAdaptor) {
            const auto& current_span = pTensorAdaptor->ViewData();
            const auto& current_tensor_shape = GetModifiedTensorShape(tensor_shape.size(), pTensorAdaptor->Shape());
            const auto axis_append_size = std::accumulate(current_tensor_shape.begin() + this->mAxis, current_tensor_shape.end(), 1u, std::multiplies<unsigned int>{});
            IndexPartition<IndexType>(current_span.size()).for_each([&combined_span, &current_span, offset, skip_size, axis_append_size](const auto Index) {
                combined_span[offset + skip_size * (Index / axis_append_size) + (Index % axis_append_size)] = current_span[Index];
            });
            offset += axis_append_size;
        });
    }

    KRATOS_CATCH("");
}

template<class TDataType>
void CombinedTensorAdaptor<TDataType>::StoreData()
{
    KRATOS_TRY

    const auto& combined_span = this->ViewData();
    IndexType offset = 0;

    // now store data to the sub tensor adaptors
    if (mAxis == 0 || mAxis == -1) {
        // for axis 0, or ravel we just need to append items.
        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [&combined_span, &offset](auto& pTensorAdaptor) {
            const auto& current_span = pTensorAdaptor->ViewData();
            IndexPartition<IndexType>(current_span.size()).for_each([&combined_span, &current_span, &offset](const auto Index) {
                current_span[Index] = combined_span[offset + Index];
            });
            offset += current_span.size();
        });
    } else {
        // this needs to append according to the axis
        const auto& tensor_shape = this->Shape();
        const auto skip_size = std::accumulate(tensor_shape.begin() + mAxis, tensor_shape.end(), 1u, std::multiplies<unsigned int>{});

        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [this, &combined_span, &offset, skip_size](auto& pTensorAdaptor) {
            const auto current_span = pTensorAdaptor->ViewData();
            const auto& current_tensor_shape = pTensorAdaptor->Shape();
            const auto axis_append_size = std::accumulate(current_tensor_shape.begin() + this->mAxis, current_tensor_shape.end(), 1u, std::multiplies<unsigned int>{});
            IndexPartition<IndexType>(current_span.size()).for_each([&combined_span, &current_span, offset, skip_size, axis_append_size](const auto Index) {
                current_span[Index] = combined_span[offset + skip_size * (Index / axis_append_size) + (Index % axis_append_size)];
            });
            offset += axis_append_size;
        });
    }

    // if recursively storing is chosen, lastly call store data in the sub TensorAdaptors
    if (mPerformStoreDataRecursively) {
        std::for_each(mTensorAdaptors.begin(), mTensorAdaptors.end(), [](auto& pTensorAdaptor) {
            pTensorAdaptor->StoreData();
        });
    }

    KRATOS_CATCH("");
}

template<class TDataType>
typename CombinedTensorAdaptor<TDataType>::TensorAdaptorVectorType CombinedTensorAdaptor<TDataType>::GetTensorAdaptors() const
{
    return mTensorAdaptors;
}

template<class TDataType>
std::string CombinedTensorAdaptor<TDataType>::Info() const
{
    std::stringstream info;
    info << "CombinedTensorAdaptor: "
         << "collect data recursively = " << (mPerformCollectDataRecursively ? "true" : "false") << ", "
         << "store data recursively = " << (mPerformStoreDataRecursively ? "true" : "false") << ", ";

    if (mAxis == -1) {
        info << "Ravel = true";
    } else {
        info << "Axis = " << mAxis;
    }
    info << ", Shape = " << this->Shape() << " with " << mTensorAdaptors.size() << "  TensorAdaptors \n[";
    for (IndexType i = 0; i < mTensorAdaptors.size(); ++i) {
        info << "\n\t" << mTensorAdaptors[i]->Info();
    }
    info << "\n]";
    return info.str();
}

// template instantiations
template class CombinedTensorAdaptor<bool>;
template class CombinedTensorAdaptor<int>;
template class CombinedTensorAdaptor<double>;

} // namespace Kratos