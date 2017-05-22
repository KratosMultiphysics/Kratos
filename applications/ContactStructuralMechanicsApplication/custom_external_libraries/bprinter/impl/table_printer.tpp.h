#include <boost/spirit/include/karma.hpp>
namespace karma = boost::spirit::karma;

namespace bprinter{
template<typename T> void TablePrinter::OutputDecimalNumber(T input){
  *out_stream_ << karma::format(
                 karma::maxwidth(column_widths_.at(j_))[
                   karma::right_align(column_widths_.at(j_))[
                     karma::double_
                   ]
                 ], input
               );

  if (j_ == get_num_columns()-1){
    *out_stream_ << "|\n";
    i_ = i_ + 1;
    j_ = 0;
  } else {
    *out_stream_ << separator_;
    j_ = j_ + 1;
  }
}
}
