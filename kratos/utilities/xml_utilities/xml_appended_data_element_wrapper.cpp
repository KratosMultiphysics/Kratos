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

// External includes
#include "zlib.h"

// Project includes
#include "utilities/data_type_traits.h"
#include "xml_elements_array.h"
#include "xml_utils.h"

// Include base h
#include "xml_appended_data_element_wrapper.h"

namespace Kratos {

XmlAppendedDataElementWrapper::XmlAppendedDataElementWrapper(const XmlOutputType OutputType)
    : BaseType("AppendedData"),
      mOutputType(OutputType)
{
    this->AddAttribute("encoding", "raw");
}

XmlElement::Pointer XmlAppendedDataElementWrapper::Get(
    const std::string& rDataArrayName,
    NDDataPointerType pNDData)
{
    KRATOS_TRY

    auto p_element = BaseType::Get(rDataArrayName, pNDData);
    p_element->AddAttribute("offset", std::to_string(mData.size()));
    p_element->AddAttribute("format", "appended");

    std::visit([this, &rDataArrayName](auto p_nd_data) {
        switch (this->mOutputType) {
            case RAW: {
                const auto number_of_values = p_nd_data->Size();
                std::uint64_t data_size = 0;
                if (number_of_values > 0) {
                    data_size = number_of_values * sizeof(*p_nd_data->ViewData().data());
                }
                XmlUtilities::AppendData(this->mData, p_nd_data->ViewData().data(), number_of_values, &data_size, 1u);
                break;
            }
            case COMPRESSED_RAW: {
                using nd_data_type = typename BareType<decltype(*p_nd_data)>::DataType;

                const auto uncompressed_data = reinterpret_cast<Bytef const *>(p_nd_data->ViewData().data());
                const auto uncompressed_length = sizeof(nd_data_type) * p_nd_data->Size();

                std::vector<std::uint64_t> header(3, 0);
                // storage of total data
                std::vector<Bytef> compressed_total_data;

                if (uncompressed_length > 0) {
                    // here we need to chunk data to 32 kB
                    const int chunk_size = 32 * 1024;
                    const IndexType number_of_chunks = ( uncompressed_length - 1 ) / chunk_size + 1;

                    // temporary storage for compressed chunks
                    uLongf compressed_length = compressBound(chunk_size);
                    std::vector<Bytef> compressed_data(compressed_length);

                    int offset = 0;
                    int remainder = uncompressed_length;

                    for (IndexType i_chunk = 0; i_chunk < number_of_chunks; ++i_chunk) {

                        const auto current_chunk_size = std::min(chunk_size, remainder);
                        uLongf current_compressed_length = compressed_length;

                        // now compress the data using zlib
                        const auto res = compress(compressed_data.data(), &current_compressed_length, uncompressed_data + offset, current_chunk_size);
                        KRATOS_ERROR_IF_NOT(res == Z_OK)
                            << "Zlib compression failed compressing \"" << rDataArrayName << "\" with data " << *p_nd_data << " [ ZLib error code = " << res << " ].\n";

                        // insert compressed data to the total compressed data
                        compressed_total_data.insert(compressed_total_data.end(), compressed_data.data(), compressed_data.data() + current_compressed_length);
                        header.push_back(current_compressed_length);

                        offset += current_chunk_size;
                        remainder -= current_chunk_size;
                    }

                    header[0] = header.size() - 3;                  // number of chunks
                    header[1] = chunk_size;                         // chunk size
                    header[2] = uncompressed_length % chunk_size;   // remainder
                    XmlUtilities::AppendData(this->mData, compressed_total_data.data(), compressed_total_data.size(), header.data(), header.size());
                    break;
                } else {
                    XmlUtilities::AppendData(this->mData, compressed_total_data.data(), 0, header.data(), header.size());
                    break;
                }
            }
        }
    }, pNDData);

    return p_element;

    KRATOS_CATCH("");
}

void XmlAppendedDataElementWrapper::Write(
    std::ostream& rOStream,
    const IndexType Level) const
{
    KRATOS_TRY

    if (mData.empty()) {
        this->WriteEmptyElementTag(rOStream, Level);
    } else {
        this->WriteElementTagStart(rOStream, Level);

        const std::string tabbing(Level * 3, ' ');

        rOStream << tabbing << "  _";

        rOStream.write(mData.data(), mData.size());

        rOStream << "\n";

        this->WriteElementTagEnd(rOStream, Level);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos