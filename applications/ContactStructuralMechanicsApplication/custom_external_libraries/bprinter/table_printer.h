#ifndef BPRINTER_TABLE_PRINTER_H_
#define BPRINTER_TABLE_PRINTER_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

namespace bprinter {
class endl{};
/** \class TablePrinter

  Print a pretty table into your output of choice.

  Usage:
    TablePrinter tp(&std::cout);
    tp.AddColumn("Name", 25);
    tp.AddColumn("Age", 3);
    tp.AddColumn("Position", 30);

    tp.PrintHeader();
    tp << "Dat Chu" << 25 << "Research Assistant";
    tp << "John Doe" << 26 << "Professional Anonymity";
    tp << "Jane Doe" << tp.SkipToNextLine();
    tp << "Tom Doe" << 7 << "Student";
    tp.PrintFooter();

  \todo Add support for padding in each table cell
  */
class TablePrinter{
public:
  TablePrinter(std::ostream * output, const std::string & separator = "|");
  ~TablePrinter();

  int get_num_columns() const;
  int get_table_width() const;
  void set_separator(const std::string & separator);
  void set_flush_left();
  void set_flush_right();

  void AddColumn(const std::string & header_name, int column_width);
  void PrintHeader();
  void PrintFooter();

  TablePrinter& operator<<(endl input){
    while (j_ != 0){
      *this << "";
    }
    return *this;
  }

  // Can we merge these?
  TablePrinter& operator<<(float input);
  TablePrinter& operator<<(double input);

  template<typename T> TablePrinter& operator<<(T input){
    if (j_ == 0)
      *out_stream_ << "|";

    if(flush_left_)
      *out_stream_ << std::left;
    else
      *out_stream_ << std::right; 

    // Leave 3 extra space: One for negative sign, one for zero, one for decimal
    *out_stream_ << std::setw(column_widths_.at(j_))
                 << input;

    if (j_ == get_num_columns()-1){
      *out_stream_ << "|\n";
      i_ = i_ + 1;
      j_ = 0;
    } else {
      *out_stream_ << separator_;
      j_ = j_ + 1;
    }

    return *this;
  }

private:
  void PrintHorizontalLine();

  template<typename T> void OutputDecimalNumber(T input);

  std::ostream * out_stream_;
  std::vector<std::string> column_headers_;
  std::vector<int> column_widths_;
  std::string separator_;

  int i_; // index of current row
  int j_; // index of current column

  int table_width_;
  bool flush_left_;
};

}

#include "impl/table_printer.tpp.h"
#endif
