#ifndef VIENNACL_SCHEDULER_STATEMENT_HPP
#define VIENNACL_SCHEDULER_STATEMENT_HPP

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/scheduler/forwards.h
    @brief Provides the datastructures for dealing with a single statement such as 'x = y + z;'
*/

#include "viennacl/forwards.h"

#include <vector>

namespace viennacl
{
  namespace scheduler
  {

    /** @brief Exception for the case the scheduler is unable to deal with the operation */
    class statement_not_supported_exception : public std::exception
    {
    public:
      statement_not_supported_exception() : message_() {}
      statement_not_supported_exception(std::string message) : message_("ViennaCL: Internal error: The scheduler encountered a problem with the operation provided: " + message) {}

      virtual const char* what() const throw() { return message_.c_str(); }

      virtual ~statement_not_supported_exception() throw() {}
    private:
      std::string message_;
    };


    /** @brief Optimization enum for grouping operations into unary or binary operations. Just for optimization of lookups. */
    enum operation_node_type_family
    {
      OPERATION_INVALID_TYPE_FAMILY = 0,

      // unary or binary expression
      OPERATION_UNARY_TYPE_FAMILY,
      OPERATION_BINARY_TYPE_FAMILY
    };

    /** @brief Enumeration for identifying the possible operations */
    enum operation_node_type
    {
      OPERATION_INVALID_TYPE = 0,

      // unary expression
      OPERATION_UNARY_ABS_TYPE,
      OPERATION_UNARY_ACOS_TYPE,
      OPERATION_UNARY_ASIN_TYPE,
      OPERATION_UNARY_ATAN_TYPE,
      OPERATION_UNARY_CEIL_TYPE,
      OPERATION_UNARY_COS_TYPE,
      OPERATION_UNARY_COSH_TYPE,
      OPERATION_UNARY_EXP_TYPE,
      OPERATION_UNARY_FABS_TYPE,
      OPERATION_UNARY_FLOOR_TYPE,
      OPERATION_UNARY_LOG_TYPE,
      OPERATION_UNARY_LOG10_TYPE,
      OPERATION_UNARY_SIN_TYPE,
      OPERATION_UNARY_SINH_TYPE,
      OPERATION_UNARY_SQRT_TYPE,
      OPERATION_UNARY_TAN_TYPE,
      OPERATION_UNARY_TANH_TYPE,
      OPERATION_UNARY_TRANS_TYPE,
      OPERATION_UNARY_NORM_1_TYPE,
      OPERATION_UNARY_NORM_2_TYPE,
      OPERATION_UNARY_NORM_INF_TYPE,

      // binary expression
      OPERATION_BINARY_ACCESS_TYPE,
      OPERATION_BINARY_ASSIGN_TYPE,
      OPERATION_BINARY_INPLACE_ADD_TYPE,
      OPERATION_BINARY_INPLACE_SUB_TYPE,
      OPERATION_BINARY_ADD_TYPE,
      OPERATION_BINARY_SUB_TYPE,
      OPERATION_BINARY_MAT_VEC_PROD_TYPE,
      OPERATION_BINARY_MAT_MAT_PROD_TYPE,
      OPERATION_BINARY_MULT_TYPE,    // scalar times vector/matrix
      OPERATION_BINARY_DIV_TYPE,     // vector/matrix divided by scalar
      OPERATION_BINARY_ELEMENT_PROD_TYPE,
      OPERATION_BINARY_ELEMENT_DIV_TYPE,
      OPERATION_BINARY_INNER_PROD_TYPE
    };



    namespace result_of
    {
      /** @brief Helper metafunction for obtaining the operation ID as well as the operation family for unary and binary operations on vectors and matrices. */
      template <typename T>
      struct op_type_info
      {
        typedef typename T::ERROR_UNKNOWN_OP_TYPE   error_type;
      };

      /** \cond */

      // unary operations
      template <> struct op_type_info<op_element_unary<op_abs>   > { enum { id = OPERATION_UNARY_ABS_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_acos>  > { enum { id = OPERATION_UNARY_ACOS_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_asin>  > { enum { id = OPERATION_UNARY_ASIN_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_atan>  > { enum { id = OPERATION_UNARY_ATAN_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_ceil>  > { enum { id = OPERATION_UNARY_CEIL_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_cos>   > { enum { id = OPERATION_UNARY_COS_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_cosh>  > { enum { id = OPERATION_UNARY_COSH_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_exp>   > { enum { id = OPERATION_UNARY_EXP_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_fabs>  > { enum { id = OPERATION_UNARY_FABS_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_floor> > { enum { id = OPERATION_UNARY_FLOOR_TYPE, family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_log>   > { enum { id = OPERATION_UNARY_LOG_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_log10> > { enum { id = OPERATION_UNARY_LOG10_TYPE, family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_sin>   > { enum { id = OPERATION_UNARY_SIN_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_sinh>  > { enum { id = OPERATION_UNARY_SINH_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_sqrt>  > { enum { id = OPERATION_UNARY_SQRT_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_tan>   > { enum { id = OPERATION_UNARY_TAN_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_unary<op_tanh>  > { enum { id = OPERATION_UNARY_TANH_TYPE,  family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_norm_1                  > { enum { id = OPERATION_UNARY_NORM_1_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_norm_2                  > { enum { id = OPERATION_UNARY_NORM_2_TYPE,   family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_norm_inf                > { enum { id = OPERATION_UNARY_NORM_INF_TYPE, family = OPERATION_UNARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_trans                   > { enum { id = OPERATION_UNARY_TRANS_TYPE, family = OPERATION_UNARY_TYPE_FAMILY }; };

      // binary operations
      template <> struct op_type_info<op_assign>                   { enum { id = OPERATION_BINARY_ASSIGN_TYPE,       family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_inplace_add>              { enum { id = OPERATION_BINARY_INPLACE_ADD_TYPE,  family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_inplace_sub>              { enum { id = OPERATION_BINARY_INPLACE_SUB_TYPE,  family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_add>                      { enum { id = OPERATION_BINARY_ADD_TYPE,          family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_sub>                      { enum { id = OPERATION_BINARY_SUB_TYPE,          family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_prod>                     { enum { id = OPERATION_BINARY_MAT_VEC_PROD_TYPE, family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_mat_mat_prod>             { enum { id = OPERATION_BINARY_MAT_MAT_PROD_TYPE, family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_mult>                     { enum { id = OPERATION_BINARY_MULT_TYPE,         family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_div>                      { enum { id = OPERATION_BINARY_DIV_TYPE,          family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_binary<op_prod> > { enum { id = OPERATION_BINARY_ELEMENT_PROD_TYPE, family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_element_binary<op_div>  > { enum { id = OPERATION_BINARY_ELEMENT_DIV_TYPE,  family = OPERATION_BINARY_TYPE_FAMILY }; };
      template <> struct op_type_info<op_inner_prod>               { enum { id = OPERATION_BINARY_INNER_PROD_TYPE,   family = OPERATION_BINARY_TYPE_FAMILY }; };

      /** \endcond */
    } // namespace result_of





    /** @brief Groups the type of a node in the statement tree. Used for faster dispatching */
    enum statement_node_type_family
    {
      INVALID_TYPE_FAMILY = 0,

      // LHS or RHS are again an expression:
      COMPOSITE_OPERATION_FAMILY,

      // device scalars:
      SCALAR_TYPE_FAMILY,

      // vector:
      VECTOR_TYPE_FAMILY,

      // matrices:
      MATRIX_TYPE_FAMILY
    };

    /** @brief Encodes the type of a node in the statement tree. */
    enum statement_node_subtype
    {
      INVALID_SUBTYPE = 0, //when type is COMPOSITE_OPERATION_FAMILY

      HOST_SCALAR_TYPE,
      DEVICE_SCALAR_TYPE,

      DENSE_VECTOR_TYPE,
      IMPLICIT_VECTOR_TYPE,

      DENSE_ROW_MATRIX_TYPE,
      DENSE_COL_MATRIX_TYPE,
      IMPLICIT_MATRIX_TYPE,

      COMPRESSED_MATRIX_TYPE,
      COORDINATE_MATRIX_TYPE,
      ELL_MATRIX_TYPE,
      HYB_MATRIX_TYPE

      // other matrix types to be added here
    };

    /** @brief Encodes the type of a node in the statement tree. */
    enum statement_node_numeric_type
    {
      INVALID_NUMERIC_TYPE = 0, //when type is COMPOSITE_OPERATION_FAMILY

      CHAR_TYPE,
      UCHAR_TYPE,
      SHORT_TYPE,
      USHORT_TYPE,
      INT_TYPE,
      UINT_TYPE,
      LONG_TYPE,
      ULONG_TYPE,
      HALF_TYPE,
      FLOAT_TYPE,
      DOUBLE_TYPE
    };


    namespace result_of
    {
      ///////////// numeric type ID deduction /////////////

      /** @brief Helper metafunction for obtaining the runtime type ID for a numerical type */
      template <typename T>
      struct numeric_type_id {};

      /** \cond */

      template <> struct numeric_type_id<char>           { enum { value = CHAR_TYPE   }; };
      template <> struct numeric_type_id<unsigned char>  { enum { value = UCHAR_TYPE  }; };
      template <> struct numeric_type_id<short>          { enum { value = SHORT_TYPE  }; };
      template <> struct numeric_type_id<unsigned short> { enum { value = USHORT_TYPE }; };
      template <> struct numeric_type_id<int>            { enum { value = INT_TYPE    }; };
      template <> struct numeric_type_id<unsigned int>   { enum { value = UINT_TYPE   }; };
      template <> struct numeric_type_id<long>           { enum { value = LONG_TYPE   }; };
      template <> struct numeric_type_id<unsigned long>  { enum { value = ULONG_TYPE  }; };
      template <> struct numeric_type_id<float>          { enum { value = FLOAT_TYPE  }; };
      template <> struct numeric_type_id<double>         { enum { value = DOUBLE_TYPE }; };

      /** \endcond */

      ///////////// matrix layout ID deduction /////////////

      /** @brief Helper metafunction for obtaining the memory layout (row-/column-major) for a matrix. */
      template <typename F>
      struct layout_type_id {};

      /** \cond */

      template <> struct layout_type_id<viennacl::column_major> { enum { value = DENSE_COL_MATRIX_TYPE }; };
      template <> struct layout_type_id<viennacl::row_major   > { enum { value = DENSE_ROW_MATRIX_TYPE }; };

      /** \endcond */
    }



    /** @brief A class representing the 'data' for the LHS or RHS operand of the respective node.
      *
      * If it represents a compound expression, the union holds the array index within the respective statement array.
      * If it represents a object with data (vector, matrix, etc.) it holds the respective pointer (scalar, vector, matrix) or value (host scalar)
      *
      * The member 'type_family' is an optimization for quickly retrieving the 'type', which denotes the currently 'active' member in the union
      */
    struct lhs_rhs_element
    {
      statement_node_type_family   type_family;
      statement_node_subtype       subtype;
      statement_node_numeric_type  numeric_type;

      union
      {
        /////// Case 1: Node is another compound expression:
        vcl_size_t        node_index;

        /////// Case 2: Node is a leaf, hence carries an operand:

        // host scalars:
        char               host_char;
        unsigned char      host_uchar;
        short              host_short;
        unsigned short     host_ushort;
        int                host_int;
        unsigned int       host_uint;
        long               host_long;
        unsigned long      host_ulong;
        float              host_float;
        double             host_double;

        // Note: ViennaCL types have potentially expensive copy-CTORs, hence using pointers:

        // scalars:
        //viennacl::scalar<char>             *scalar_char;
        //viennacl::scalar<unsigned char>    *scalar_uchar;
        //viennacl::scalar<short>            *scalar_short;
        //viennacl::scalar<unsigned short>   *scalar_ushort;
        //viennacl::scalar<int>              *scalar_int;
        //viennacl::scalar<unsigned int>     *scalar_uint;
        //viennacl::scalar<long>             *scalar_long;
        //viennacl::scalar<unsigned long>    *scalar_ulong;
        viennacl::scalar<float>            *scalar_float;
        viennacl::scalar<double>           *scalar_double;

        // vectors:
        //viennacl::vector_base<char>             *vector_char;
        //viennacl::vector_base<unsigned char>    *vector_uchar;
        //viennacl::vector_base<short>            *vector_short;
        //viennacl::vector_base<unsigned short>   *vector_ushort;
        //viennacl::vector_base<int>              *vector_int;
        //viennacl::vector_base<unsigned int>     *vector_uint;
        //viennacl::vector_base<long>             *vector_long;
        //viennacl::vector_base<unsigned long>    *vector_ulong;
        viennacl::vector_base<float>            *vector_float;
        viennacl::vector_base<double>           *vector_double;

        // implicit vectors:
        //viennacl::implicit_vector_base<char>             *implicit_vector_char;
        //viennacl::implicit_vector_base<unsigned char>    *implicit_vector_uchar;
        //viennacl::implicit_vector_base<short>            *implicit_vector_short;
        //viennacl::implicit_vector_base<unsigned short>   *implicit_vector_ushort;
        //viennacl::implicit_vector_base<int>              *implicit_vector_int;
        //viennacl::implicit_vector_base<unsigned int>     *implicit_vector_uint;
        //viennacl::implicit_vector_base<long>             *implicit_vector_long;
        //viennacl::implicit_vector_base<unsigned long>    *implicit_vector_ulong;
        viennacl::implicit_vector_base<float>            *implicit_vector_float;
        viennacl::implicit_vector_base<double>           *implicit_vector_double;

        // row-major matrices:
        //viennacl::matrix_base<char>             *matrix_row_char;
        //viennacl::matrix_base<unsigned char>    *matrix_row_uchar;
        //viennacl::matrix_base<short>            *matrix_row_short;
        //viennacl::matrix_base<unsigned short>   *matrix_row_ushort;
        //viennacl::matrix_base<int>              *matrix_row_int;
        //viennacl::matrix_base<unsigned int>     *matrix_row_uint;
        //viennacl::matrix_base<long>             *matrix_row_long;
        //viennacl::matrix_base<unsigned long>    *matrix_row_ulong;
        viennacl::matrix_base<float>            *matrix_row_float;
        viennacl::matrix_base<double>           *matrix_row_double;

        // column-major matrices:
        //viennacl::matrix_base<char,           viennacl::column_major>    *matrix_col_char;
        //viennacl::matrix_base<unsigned char,  viennacl::column_major>    *matrix_col_uchar;
        //viennacl::matrix_base<short,          viennacl::column_major>    *matrix_col_short;
        //viennacl::matrix_base<unsigned short, viennacl::column_major>    *matrix_col_ushort;
        //viennacl::matrix_base<int,            viennacl::column_major>    *matrix_col_int;
        //viennacl::matrix_base<unsigned int,   viennacl::column_major>    *matrix_col_uint;
        //viennacl::matrix_base<long,           viennacl::column_major>    *matrix_col_long;
        //viennacl::matrix_base<unsigned long,  viennacl::column_major>    *matrix_col_ulong;
        viennacl::matrix_base<float,          viennacl::column_major>    *matrix_col_float;
        viennacl::matrix_base<double,         viennacl::column_major>    *matrix_col_double;

        //viennacl::implicit_matrix_base<char>             *implicit_matrix_char;
        //viennacl::implicit_matrix_base<unsigned char>    *implicit_matrix_uchar;
        //viennacl::implicit_matrix_base<short>            *implicit_matrix_short;
        //viennacl::implicit_matrix_base<unsigned short>   *implicit_matrix_ushort;
        //viennacl::implicit_matrix_base<int>              *implicit_matrix_int;
        //viennacl::implicit_matrix_base<unsigned int>     *implicit_matrix_uint;
        //viennacl::implicit_matrix_base<long>             *implicit_matrix_long;
        //viennacl::implicit_matrix_base<unsigned long>    *implicit_matrix_ulong;
        viennacl::implicit_matrix_base<float>            *implicit_matrix_float;
        viennacl::implicit_matrix_base<double>           *implicit_matrix_double;

        //viennacl::compressed_matrix<float>    *compressed_matrix_char;
        //viennacl::compressed_matrix<double>   *compressed_matrix_uchar;
        //viennacl::compressed_matrix<float>    *compressed_matrix_short;
        //viennacl::compressed_matrix<double>   *compressed_matrix_ushort;
        //viennacl::compressed_matrix<float>    *compressed_matrix_int;
        //viennacl::compressed_matrix<double>   *compressed_matrix_uint;
        //viennacl::compressed_matrix<float>    *compressed_matrix_long;
        //viennacl::compressed_matrix<double>   *compressed_matrix_ulong;
        viennacl::compressed_matrix<float>    *compressed_matrix_float;
        viennacl::compressed_matrix<double>   *compressed_matrix_double;

        //viennacl::coordinate_matrix<float>    *coordinate_matrix_char;
        //viennacl::coordinate_matrix<double>   *coordinate_matrix_uchar;
        //viennacl::coordinate_matrix<float>    *coordinate_matrix_short;
        //viennacl::coordinate_matrix<double>   *coordinate_matrix_ushort;
        //viennacl::coordinate_matrix<float>    *coordinate_matrix_int;
        //viennacl::coordinate_matrix<double>   *coordinate_matrix_uint;
        //viennacl::coordinate_matrix<float>    *coordinate_matrix_long;
        //viennacl::coordinate_matrix<double>   *coordinate_matrix_ulong;
        viennacl::coordinate_matrix<float>    *coordinate_matrix_float;
        viennacl::coordinate_matrix<double>   *coordinate_matrix_double;

        //viennacl::ell_matrix<float>    *ell_matrix_char;
        //viennacl::ell_matrix<double>   *ell_matrix_uchar;
        //viennacl::ell_matrix<float>    *ell_matrix_short;
        //viennacl::ell_matrix<double>   *ell_matrix_ushort;
        //viennacl::ell_matrix<float>    *ell_matrix_int;
        //viennacl::ell_matrix<double>   *ell_matrix_uint;
        //viennacl::ell_matrix<float>    *ell_matrix_long;
        //viennacl::ell_matrix<double>   *ell_matrix_ulong;
        viennacl::ell_matrix<float>    *ell_matrix_float;
        viennacl::ell_matrix<double>   *ell_matrix_double;

        //viennacl::hyb_matrix<float>    *hyb_matrix_char;
        //viennacl::hyb_matrix<double>   *hyb_matrix_uchar;
        //viennacl::hyb_matrix<float>    *hyb_matrix_short;
        //viennacl::hyb_matrix<double>   *hyb_matrix_ushort;
        //viennacl::hyb_matrix<float>    *hyb_matrix_int;
        //viennacl::hyb_matrix<double>   *hyb_matrix_uint;
        //viennacl::hyb_matrix<float>    *hyb_matrix_long;
        //viennacl::hyb_matrix<double>   *hyb_matrix_ulong;
        viennacl::hyb_matrix<float>    *hyb_matrix_float;
        viennacl::hyb_matrix<double>   *hyb_matrix_double;
      };
    };


    /** @brief Struct for holding the type family as well as the type of an operation (could be addition, subtraction, norm, etc.) */
    struct op_element
    {
      operation_node_type_family   type_family;
      operation_node_type          type;
    };

    /** @brief Main datastructure for an node in the statement tree */
    struct statement_node
    {
      lhs_rhs_element    lhs;
      op_element         op;
      lhs_rhs_element    rhs;
    };

    namespace result_of
    {

      /** @brief Helper metafunction for obtaining the number of nodes of an expression template tree. */
      template <class T> struct num_nodes { enum { value = 0 }; };
      /** \cond */
      template <class LHS, class OP, class RHS> struct num_nodes<       vector_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      template <class LHS, class OP, class RHS> struct num_nodes< const vector_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      template <class LHS, class OP, class RHS> struct num_nodes<       matrix_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      template <class LHS, class OP, class RHS> struct num_nodes< const matrix_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      template <class LHS, class OP, class RHS> struct num_nodes<       scalar_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      template <class LHS, class OP, class RHS> struct num_nodes< const scalar_expression<LHS, RHS, OP> > { enum { value = 1 + num_nodes<LHS>::value + num_nodes<RHS>::value }; };
      /** \endcond */

    }

    /** \brief The main class for representing a statement such as x = inner_prod(y,z); at runtime.
      *
      * This is the equivalent to an expression template tree, but entirely built at runtime in order to perform really cool stuff such as kernel fusion.
      */
    class statement
    {
      public:
        typedef statement_node              value_type;
        typedef viennacl::vcl_size_t        size_type;
        typedef std::vector<value_type>     container_type;

        statement(container_type const & custom_array) : array_(custom_array) {}

        /** @brief Generate the runtime statement from an expression template.
          *
          * Constructing a runtime statement from expression templates makes perfect sense, because this way only a single allocation is needed when creating the statement. */
        template <typename LHS, typename OP, typename RHS>
        statement(LHS & lhs, OP const &, RHS const & rhs) : array_(1 + result_of::num_nodes<RHS>::value)
        {
          // set OP:
          array_[0].op.type_family = operation_node_type_family(result_of::op_type_info<OP>::family);
          array_[0].op.type        = operation_node_type(result_of::op_type_info<OP>::id);

          // set LHS:
          add_lhs(0, 1, lhs);

          // set RHS:
          add_rhs(0, 1, rhs);
        }

        container_type const & array() const { return array_; }

        size_type root() const { return 0; }

      private:

        ///////////// Scalar node helper ////////////////

        // TODO: add integer vector overloads here
        void assign_element(lhs_rhs_element & elem, viennacl::scalar<float>  const & t) { elem.scalar_float  = const_cast<viennacl::scalar<float> *>(&t); }
        void assign_element(lhs_rhs_element & elem, viennacl::scalar<double> const & t) { elem.scalar_double = const_cast<viennacl::scalar<double> *>(&t); }

        ///////////// Vector node helper ////////////////
        // TODO: add integer vector overloads here
        void assign_element(lhs_rhs_element & elem, viennacl::vector_base<float>  const & t) { elem.vector_float  = const_cast<viennacl::vector_base<float> *>(&t); }
        void assign_element(lhs_rhs_element & elem, viennacl::vector_base<double> const & t) { elem.vector_double = const_cast<viennacl::vector_base<double> *>(&t); }

        ///////////// Matrix node helper ////////////////
        // TODO: add integer matrix overloads here
        void assign_element(lhs_rhs_element & elem, viennacl::matrix_base<float,  viennacl::column_major> const & t) { elem.matrix_col_float  = const_cast<viennacl::matrix_base<float,  viennacl::column_major> *>(&t); }
        void assign_element(lhs_rhs_element & elem, viennacl::matrix_base<float,  viennacl::row_major>    const & t) { elem.matrix_row_float  = const_cast<viennacl::matrix_base<float,  viennacl::row_major>    *>(&t); }
        void assign_element(lhs_rhs_element & elem, viennacl::matrix_base<double, viennacl::column_major> const & t) { elem.matrix_col_double = const_cast<viennacl::matrix_base<double, viennacl::column_major> *>(&t); }
        void assign_element(lhs_rhs_element & elem, viennacl::matrix_base<double, viennacl::row_major>    const & t) { elem.matrix_row_double = const_cast<viennacl::matrix_base<double, viennacl::row_major>    *>(&t); }

        void assign_element(lhs_rhs_element & elem, viennacl::compressed_matrix<float>  const & m) { elem.compressed_matrix_float  = const_cast<viennacl::compressed_matrix<float>  *>(&m); }
        void assign_element(lhs_rhs_element & elem, viennacl::compressed_matrix<double> const & m) { elem.compressed_matrix_double = const_cast<viennacl::compressed_matrix<double> *>(&m); }

        void assign_element(lhs_rhs_element & elem, viennacl::coordinate_matrix<float>  const & m) { elem.coordinate_matrix_float  = const_cast<viennacl::coordinate_matrix<float>  *>(&m); }
        void assign_element(lhs_rhs_element & elem, viennacl::coordinate_matrix<double> const & m) { elem.coordinate_matrix_double = const_cast<viennacl::coordinate_matrix<double> *>(&m); }

        void assign_element(lhs_rhs_element & elem, viennacl::ell_matrix<float>  const & m) { elem.ell_matrix_float  = const_cast<viennacl::ell_matrix<float>  *>(&m); }
        void assign_element(lhs_rhs_element & elem, viennacl::ell_matrix<double> const & m) { elem.ell_matrix_double = const_cast<viennacl::ell_matrix<double> *>(&m); }

        void assign_element(lhs_rhs_element & elem, viennacl::hyb_matrix<float>  const & m) { elem.hyb_matrix_float  = const_cast<viennacl::hyb_matrix<float>  *>(&m); }
        void assign_element(lhs_rhs_element & elem, viennacl::hyb_matrix<double> const & m) { elem.hyb_matrix_double = const_cast<viennacl::hyb_matrix<double> *>(&m); }

        //////////// Tree leaves (terminals) ////////////////////

        vcl_size_t add_element(vcl_size_t       next_free,
                                lhs_rhs_element & elem,
                                float const &     t)
        {
          elem.type_family  = SCALAR_TYPE_FAMILY;
          elem.subtype      = HOST_SCALAR_TYPE;
          elem.numeric_type = FLOAT_TYPE;
          elem.host_float   = t;
          return next_free;
        }

        vcl_size_t add_element(vcl_size_t       next_free,
                                lhs_rhs_element & elem,
                                double const &    t)
        {
          elem.type_family  = SCALAR_TYPE_FAMILY;
          elem.subtype      = HOST_SCALAR_TYPE;
          elem.numeric_type = DOUBLE_TYPE;
          elem.host_double  = t;
          return next_free;
        }

        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::scalar<T> const & t)
        {
          elem.type_family  = SCALAR_TYPE_FAMILY;
          elem.subtype      = DEVICE_SCALAR_TYPE;
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }


        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::vector_base<T> const & t)
        {
          elem.type_family           = VECTOR_TYPE_FAMILY;
          elem.subtype               = DENSE_VECTOR_TYPE;
          elem.numeric_type          = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }

        template <typename T, typename F>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::matrix_base<T, F> const & t)
        {
          elem.type_family  = MATRIX_TYPE_FAMILY;
          elem.subtype      = statement_node_subtype(result_of::layout_type_id<F>::value);
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }

        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::compressed_matrix<T> const & t)
        {
          elem.type_family  = MATRIX_TYPE_FAMILY;
          elem.subtype      = COMPRESSED_MATRIX_TYPE;
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }

        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::coordinate_matrix<T> const & t)
        {
          elem.type_family  = MATRIX_TYPE_FAMILY;
          elem.subtype      = COORDINATE_MATRIX_TYPE;
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }

        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::ell_matrix<T> const & t)
        {
          elem.type_family  = MATRIX_TYPE_FAMILY;
          elem.subtype      = ELL_MATRIX_TYPE;
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }

        template <typename T>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element            & elem,
                                viennacl::hyb_matrix<T> const & t)
        {
          elem.type_family  = MATRIX_TYPE_FAMILY;
          elem.subtype      = HYB_MATRIX_TYPE;
          elem.numeric_type = statement_node_numeric_type(result_of::numeric_type_id<T>::value);
          assign_element(elem, t);
          return next_free;
        }


        //////////// Tree nodes (non-terminals) ////////////////////

        template <typename LHS, typename RHS, typename OP>
        vcl_size_t add_element(vcl_size_t       next_free,
                                lhs_rhs_element & elem,
                                viennacl::scalar_expression<LHS, RHS, OP> const & t)
        {
          elem.type_family  = COMPOSITE_OPERATION_FAMILY;
          elem.subtype      = INVALID_SUBTYPE;
          elem.numeric_type = INVALID_NUMERIC_TYPE;
          elem.node_index   = next_free;
          return add_node(next_free, next_free + 1, t);
        }

        template <typename LHS, typename RHS, typename OP>
        vcl_size_t add_element(vcl_size_t       next_free,
                                lhs_rhs_element & elem,
                                viennacl::vector_expression<LHS, RHS, OP> const & t)
        {
          elem.type_family  = COMPOSITE_OPERATION_FAMILY;
          elem.subtype      = INVALID_SUBTYPE;
          elem.numeric_type = INVALID_NUMERIC_TYPE;
          elem.node_index   = next_free;
          return add_node(next_free, next_free + 1, t);
        }

        template <typename LHS, typename RHS, typename OP>
        vcl_size_t add_element(vcl_size_t next_free,
                                lhs_rhs_element & elem,
                                viennacl::matrix_expression<LHS, RHS, OP> const & t)
        {
          elem.type_family   = COMPOSITE_OPERATION_FAMILY;
          elem.subtype      = INVALID_SUBTYPE;
          elem.numeric_type = INVALID_NUMERIC_TYPE;
          elem.node_index    = next_free;
          return add_node(next_free, next_free + 1, t);
        }


        //////////// Helper routines ////////////////////


        template <typename T>
        vcl_size_t add_lhs(vcl_size_t current_index, vcl_size_t next_free, T const & t)
        {
          return add_element(next_free, array_[current_index].lhs, t);
        }

        template <typename T>
        vcl_size_t add_rhs(vcl_size_t current_index, vcl_size_t next_free, T const & t)
        {
          return add_element(next_free, array_[current_index].rhs, t);
        }

        //////////// Internal interfaces ////////////////////

        template <template <typename, typename, typename> class ExpressionT, typename LHS, typename RHS, typename OP>
        vcl_size_t add_node(vcl_size_t current_index, vcl_size_t next_free, ExpressionT<LHS, RHS, OP> const & proxy)
        {
          // set OP:
          array_[current_index].op.type_family = operation_node_type_family(result_of::op_type_info<OP>::family);
          array_[current_index].op.type        = operation_node_type(result_of::op_type_info<OP>::id);

          // set LHS and RHS:
          if (array_[current_index].op.type_family == OPERATION_UNARY_TYPE_FAMILY)
          {
            // unary expression: set rhs to invalid:
            array_[current_index].rhs.type_family  = INVALID_TYPE_FAMILY;
            array_[current_index].rhs.subtype      = INVALID_SUBTYPE;
            array_[current_index].rhs.numeric_type = INVALID_NUMERIC_TYPE;
            return add_lhs(current_index, next_free, proxy.lhs());
          }

          return add_rhs(current_index, add_lhs(current_index, next_free, proxy.lhs()), proxy.rhs());

        }

        container_type   array_;
    };

    namespace detail
    {
      /** @brief Deals with x = RHS where RHS is an expression and x is either a scalar, a vector, or a matrix */
      inline void execute_composite(statement const & /* s */, statement_node const & /* root_node */);
    }

  } // namespace scheduler

} // namespace viennacl

#endif

