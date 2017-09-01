#ifndef BPRINTER_TABLE_PRINTER_H_
#define BPRINTER_TABLE_PRINTER_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#if defined(USE_BOOST_KARMA)
#include <boost/spirit/include/karma.hpp>
namespace karma = boost::spirit::karma;
#endif //USE_BOOST_KARMA

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
  TablePrinter(std::ostream * output, const std::string & separator = "|", const bool use_bool_font = true);
  ~TablePrinter();

  unsigned int get_num_columns() const;
  unsigned int get_table_width() const;
  void set_separator(const std::string & separator);
  void set_flush_left();
  void set_flush_right();

  void AddColumn(const std::string & header_name, int column_width);
  void PrintHeader();
  void PrintFooter();

  TablePrinter& operator<<(endl input)
  {
    while (j_ != 0)
    {
        *this << "";
    }
    
    return *this;
  }

  // Can we merge these?
  TablePrinter& operator<<(float input);
  TablePrinter& operator<<(double input);

  template<typename TClass> 
  TablePrinter& operator<<(TClass input)
  {
      if (j_ == 0)
      {
          *out_stream_ << "|";
      }
      
      if(flush_left_)
      {
          *out_stream_ << std::left;
      }
      else
      {
          *out_stream_ << std::right; 
      }

      // Leave 3 extra space: One for negative sign, one for zero, one for decimal
      *out_stream_ << std::setw(column_widths_.at(j_)) << input;

      if (j_ == get_num_columns()-1)
      {
          *out_stream_ << "|\n";
          i_ = i_ + 1;
          j_ = 0;
      } 
      else 
      {
          *out_stream_ << separator_;
          j_ = j_ + 1;
      }
      
      return *this;
  }

private:
    void PrintHorizontalLine();
    
    #if defined(USE_BOOST_KARMA)
    template<typename TClass> 
    void OutputDecimalNumber(TClass input)
    {
        *out_stream_ << karma::format(
                    karma::maxwidth(column_widths_.at(j_))[
                    karma::right_align(column_widths_.at(j_))[
                        karma::double_
                    ]
                    ], input
                );
        
        if (j_ == get_num_columns()-1)
        {
            *out_stream_ << "|\n";
            i_ = i_ + 1;
            j_ = 0;
        } 
        else 
        {
            *out_stream_ << separator_;
            j_ = j_ + 1;
        }
    }
    #else
    template<typename TClass> 
    void OutputDecimalNumber(TClass input)
    {
        // If we cannot handle this number, indicate so
        if (input < 10*(column_widths_.at(j_)-1) || input > 10*column_widths_.at(j_))
        {
            std::stringstream string_out;
            string_out 
//             << std::setiosflags(std::ios::fixed)
            << std::setiosflags(std::ios::scientific)
//             << std::setprecision(column_widths_.at(j_))
            << std::setprecision(3)
            << std::uppercase
            << std::setw(column_widths_.at(j_))
            << input;

            std::string string_to_print = string_out.str();
//             std::string string_rep_of_number = string_out.str();
// 
//             string_rep_of_number[column_widths_.at(j_)-1] = '*';
//             std::string string_to_print = string_rep_of_number.substr(0, column_widths_.at(j_));
            *out_stream_ << string_to_print;
        } 
        else 
        {
//             // Determine what precision we need
//             int precision = column_widths_.at(j_) - 1; // leave room for the decimal point
//             if (input < 0)
//             {
//                 --precision; // leave room for the minus sign
//             }
// 
//             // Leave room for digits before the decimal?
//             if (input < -1 || input > 1)
//             {
//                 unsigned int num_digits_before_decimal = 1 + (int)log10(std::abs(input));
//                 precision -= num_digits_before_decimal;
//             }
//             else
//             {
//                 precision --; // e.g. 0.12345 or -0.1234
//             }
// 
//             if (precision < 0)
//             {
//                 precision = 0; // don't go negative with precision
//             }

            *out_stream_ 
//             << std::setiosflags(std::ios::fixed)
            << std::setiosflags(std::ios::scientific)
//             << std::setprecision(precision)
            << std::setprecision(3)
            << std::uppercase
            << std::setw(column_widths_.at(j_))
            << input;
        }

        if (j_ == get_num_columns()-1)
        {
            *out_stream_ << "|\n";
            i_ = i_ + 1;
            j_ = 0;
        } 
        else 
        {
            *out_stream_ << separator_;
            j_ = j_ + 1;
        }
    }
    #endif //USE_BOOST_KARMA
    
    std::ostream * out_stream_;
    std::vector<std::string> column_headers_;
    std::vector<int> column_widths_;
    std::string separator_;

    unsigned int i_; // index of current row
    unsigned int j_; // index of current column

    unsigned int table_width_;
    bool flush_left_;
    bool bold_font_;
};

}
#endif
