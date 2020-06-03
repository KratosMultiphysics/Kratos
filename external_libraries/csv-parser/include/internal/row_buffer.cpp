/** @file
 *  Defines an object which can store CSV data in
 *  continuous regions of memory
 */

#include "row_buffer.hpp"

namespace csv {
    namespace internals {
        //////////////
        // ColNames //
        //////////////
        CSV_INLINE std::vector<std::string> ColNames::get_col_names() const {
            return this->col_names;
        }

        CSV_INLINE void ColNames::set_col_names(const std::vector<std::string>& cnames) {
            this->col_names = cnames;

            for (size_t i = 0; i < cnames.size(); i++) {
                this->col_pos[cnames[i]] = i;
            }
        }

        CSV_INLINE int ColNames::index_of(csv::string_view col_name) const {
            auto pos = this->col_pos.find(col_name.data());
            if (pos != this->col_pos.end())
                return (int)pos->second;

            return CSV_NOT_FOUND;
        }

        CSV_INLINE size_t ColNames::size() const {
            return this->col_names.size();
        }

        CSV_INLINE RowData RawRowBuffer::get_row() {
            return {
                this->get_row_string(),
                this->get_splits()
            };
        }

        /** Get the current row in the buffer
         *  @note Has the side effect of updating the current end pointer
         */
        CSV_INLINE std::pair<size_t, size_t> RawRowBuffer::get_row_string() {
            auto ret = std::make_pair(
                this->current_end, // Beginning of string
                (this->buffer.size() - this->current_end) // Count
            );

            this->current_end = this->buffer.size();
            return ret;
        }

        CSV_INLINE ColumnPositions RawRowBuffer::get_splits()
        {
            const size_t head_idx = this->current_split_idx,
                new_split_idx = this->split_buffer.size();
            StrBufferPos n_cols = (new_split_idx - head_idx > 0) ?
                (StrBufferPos)(new_split_idx - head_idx + 1): 0;

            this->current_split_idx = new_split_idx;
            return ColumnPositions(head_idx, n_cols);
        }

        CSV_INLINE size_t RawRowBuffer::size() const {
            return this->buffer.size() - this->current_end;
        }

        CSV_INLINE size_t RawRowBuffer::splits_size() const {
            return this->split_buffer.size() - this->current_split_idx;
        }
        
        HEDLEY_WARN_UNUSED_RESULT CSV_INLINE
        BufferPtr RawRowBuffer::reset() const {
            // Save current row in progress
            auto new_buff = BufferPtr(new RawRowBuffer());

            // Save text
            new_buff->buffer = this->buffer.substr(
                this->current_end,   // Position
                (this->buffer.size() - this->current_end) // Count
            );

            // Save split buffer in progress
            for (size_t i = this->current_split_idx; i < this->split_buffer.size(); i++) {
                new_buff->split_buffer.push_back(this->split_buffer[i]);
            }

            new_buff->col_names = this->col_names;

            // No need to remove unnecessary bits from this buffer
            // (memory savings would be marginal anyways)
            return new_buff;
        }
    }
}