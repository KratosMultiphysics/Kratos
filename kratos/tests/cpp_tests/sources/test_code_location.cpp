//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//
	           
// System includes


// External includes 


// Project includes
#include "testing/testing.h"



namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(CodeLocation, KratosCoreFastSuite)
		{
			CodeLocation location1(R"string([C:\kratos\kratos\sources\model_part_io.cpp)string",
				R"string(iterators::indirect_iterator< std::_Vector_iterator< std::_Vector_val<struct std::_Simple_types< shared_ptr< Element> > > >, struct iterators::use_default, struct iterators::use_default, struct iterators::use_default, struct iterators::use_default>  ModelPartIO::FindKey< PointerVectorSet< Element, IndexedObject, struct std::less<unsigned __int64>, struct std::equal_to<unsigned __int64>, shared_ptr< Element>, std::vector< shared_ptr< Element>, std::allocator< shared_ptr< Element> > > >, unsigned __int64>(PointerVectorSet< Element, IndexedObject, struct std::less<unsigned __int64>, struct std::equal_to<unsigned __int64>, shared_ptr< Element>, std::vector< shared_ptr< Element>, std::allocator< shared_ptr< Element> > > > &, unsigned __int64, std::basic_string<char, struct std::char_traits<char>, std::allocator<char> >))string", 3552);
			KRATOS_CHECK_LESS(location1.CleanFunctionName().size(), 240); // the original size was 826
			KRATOS_CHECK_C_STRING_EQUAL(location1.CleanFileName().c_str(), "kratos/sources/model_part_io.cpp");

			std::string functions = R"string(void  ResidualBasedBlockBuilderAndSolver< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, LinearSolver< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, Reorderer< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > > > > >::SystemSolveWithPhysics( ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> > &, ublas::vector<double, ublas::unbounded_array<double> > &, ublas::vector<double, ublas::unbounded_array<double> > &, ModelPart &)
void  ResidualBasedBlockBuilderAndSolver< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, LinearSolver< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, Reorderer< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > > > > >::BuildAndSolve( shared_ptr< Scheme< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > > > >, ModelPart &, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> > &, ublas::vector<double, ublas::unbounded_array<double> > &, ublas::vector<double, ublas::unbounded_array<double> > &) [ C:\KratosSource\kratos\solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h , Line 610 ] 
 
double  ResidualBasedLinearStrategy< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, LinearSolver< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, Reorderer< UblasSpace<double, ublas::compressed_matrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>,0, ublas::unbounded_array<unsigned __int64, std::allocator<unsigned __int64> >, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > >, UblasSpace<double, ublas::DenseMatrix<double,struct ublas::basic_row_major<unsigned __int64,__int64>, ublas::unbounded_array<double> >, ublas::vector<double, ublas::unbounded_array<double> > > > > >::Solve(void) [ C:\KratosSource\kratos\solving_strategies/strategies/residualbased_linear_strategy.h , Line 465 ] 

)string";
			CodeLocation location2(R"string(C:\KratosSource\kratos\solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h)string", functions, 559);
			KRATOS_CHECK_LESS(location2.CleanFunctionName().size(), 1000); // the original size was 5858
			//std::cout << location2.GetFunctionName() << std::endl;
			std::cout << "size: " << location2.GetFunctionName().size() << std::endl;
			//std::cout << " ERROR -------------------------------------------------------------------------------------------------------------------------" << std::endl;
			//std::cout << location2.CleanFunctionName() << std::endl;
			//std::cout << " -------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			//std::cout << "size: " << location2.CleanFunctionName().size() << std::endl;
		}

		KRATOS_TEST_CASE_IN_SUITE(CodeLocationTemplateParameterReduction, KratosCoreFastSuite)
		{
			CodeLocation location1("", "ublas::vector<double, ublas::unbounded_array<double> >", 0);
			KRATOS_CHECK_C_STRING_EQUAL(location1.CleanFunctionName().c_str(), "Vector");
			CodeLocation location2("", "iterators::indirect_iterator< std::_Vector_iterator< std::_Vector_val<struct std::_Simple_types< shared_ptr< Element> > > >, struct iterators::use_default, struct iterators::use_default, struct iterators::use_default, struct iterators::use_default>", 0);
			KRATOS_CHECK_C_STRING_EQUAL(location2.CleanFunctionName().c_str(), "iterators::indirect_iterator< _Vector_iterator< _Vector_val<struct _Simple_types< shared_ptr< Element> > > >,...>");

		}

	}
}  // namespace Kratos.


