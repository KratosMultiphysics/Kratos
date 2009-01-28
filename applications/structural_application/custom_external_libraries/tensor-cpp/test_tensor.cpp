/**
 * Tests of the class Tensor.
 */

#include <iostream>
#include <vector>

#include "tensor.h"

using namespace std;


void test_vector();
void test_matrix();
void test_addition();
void test_substraction();
void test_product();
void test_contraction();
void test_swap();
void test_brackets();
void test_base_conversion();


int main()
{
    test_vector();
    test_matrix();
    test_addition();
    test_substraction();
    test_product();
    test_contraction();
    test_swap();
    test_brackets();

    test_base_conversion();

    return 0;
}




/*
 * Tensor as a vector.
 */
void test_vector()
{
    cout << "Testing tensor as a vector...\n";

    Tensor<double> v(1, 3);  // rank=1, dimension=3

    v(0) = 1.0;
    v(1) = 2.0;
    v(2) = 3.0;

    cout << "Rank:      " << v.rank()      << endl;
    cout << "Dimension: " << v.dimension() << endl;
    cout << "Value:\n"    << v             << endl;

    cout << "--------" << endl;
}




/*
 * Tensor as a matrix.
 */
void test_matrix()
{
    cout << "Testing tensor as a matrix...\n";

    Tensor<double> m(2, 4);  // rank=2, dimension=4
    
    m(0, 0) = 1;   m(0, 1) = 2;   m(0, 2) = 3;   m(0, 3) = 4;
    m(1, 0) = 5;   m(1, 1) = 6;   m(1, 2) = 7;   m(1, 3) = 8;
    m(2, 0) = 9;   m(2, 1) = 10;  m(2, 2) = 11;  m(2, 3) = 12;
    m(3, 0) = 13;  m(3, 1) = 14;  m(3, 2) = 15;  m(3, 3) = 16;

    cout << "Rank:      " << m.rank()      << endl;
    cout << "Dimension: " << m.dimension() << endl;
    cout << "Value:\n"    << m             << endl;

    cout << "----------" << endl;
}




/*
 * Addition.
 */
void test_addition()
{
    cout << "Testing tensor addition, rank 1...\n";

    Tensor<double> v1(1, 3);  // rank=1, dimension=3

    v1(0) = 1.0;  v1(1) = 2.0;  v1(2) = 3.0;

    Tensor<double> v2(1, 3);  // rank=1, dimension=3

    v2(0) = 4.0;  v2(1) = 5.0;  v2(2) = 6.0;

    cout << v1 << "\n + \n" << v2 << "\n = \n" << v1 + v2 << endl;
    
    cout << "----------" << endl;

    /*
     * More complicated.
     */
    cout << "Testing tensor addition, rank 2...\n";
    Tensor<double> m1(2, 3);
    Tensor<double> m2(2, 3);

    m1(0, 0) = 1;   m1(0, 1) = 2;   m1(0, 2) = 3;
    m1(1, 0) = 4;   m1(1, 1) = 5;   m1(1, 2) = 6;
    m1(2, 0) = 7;   m1(2, 1) = 8;   m1(2, 2) = 9;

    m2(0, 0) = 10;  m2(0, 1) = 11;  m2(0, 2) = 12;
    m2(1, 0) = 13;  m2(1, 1) = 14;  m2(1, 2) = 15;
    m2(2, 0) = 16;  m2(2, 1) = 17;  m2(2, 2) = 18;

    cout << m1 << "\n + \n" << m2 << "\n = \n" << m1 + m2 << endl;

    cout << "----------" << endl;
}


/*
 * Substraction
 */
void test_substraction()
{
    cout << "Testing tensor substraction rank 1...\n";

    Tensor<double> v1(1, 3);  // rank=1, dimension=3

    v1(0) = 1.0;  v1(1) = 2.0;  v1(2) = 3.0;

    Tensor<double> v2(1, 3);  // rank=1, dimension=3

    v2(0) = 4.0;  v2(1) = 5.0;  v2(2) = 6.0;

    cout << v1 << "\n - \n" << v2 << "\n = \n" << v1 - v2 << endl;
    
    cout << "----------" << endl;

    /*
     * More complicated.
     */
    cout << "Testing tensor substraction, rank 2...\n";
    Tensor<double> m1(2, 3);
    Tensor<double> m2(2, 3);

    m1(0, 0) = 1;   m1(0, 1) = 2;   m1(0, 2) = 3;
    m1(1, 0) = 4;   m1(1, 1) = 5;   m1(1, 2) = 6;
    m1(2, 0) = 7;   m1(2, 1) = 8;   m1(2, 2) = 9;

    m2(0, 0) = 10;  m2(0, 1) = 11;  m2(0, 2) = 12;
    m2(1, 0) = 13;  m2(1, 1) = 14;  m2(1, 2) = 15;
    m2(2, 0) = 16;  m2(2, 1) = 17;  m2(2, 2) = 18;

    cout << m1 << "\n - \n" << m2 << "\n = \n" << m1 - m2 << endl;

    cout << "----------" << endl;
}



/*
 * Tensor product.
 */
void test_product()
{
    cout << "Testing tensor product, vector x vector...\n";

    Tensor<double> v1(1, 2);
    Tensor<double> v2(1, 2);

    v1(0) = 1;
    v1(1) = 2;

    v2(0) = 3;
    v2(1) = 4;

    Tensor<double> p = v1 * v2;

    cout << v1 << "\n X \n" << v2 << "\n = \n" << p << endl;
    cout << "That was, a rank " << v1.rank() << " tensor times a rank "
	 << v2.rank() << " tensor gives a rank " << p.rank() << " tensor\n";

    cout << "--------"  << endl;

    /*
     * More complicated.
     */
    cout << "Testing tensor product, matrix x vector...\n";

    Tensor<double> m(2, 2);
    Tensor<double> v(1, 2);

    m(0, 0) = 1;  m(0, 1) = 2;
    m(1, 0) = 3;  m(1, 1) = 4;
    
    v(0) = 5;
    v(1) = 6;

    Tensor<double> p2 = m * v;

    cout << m << "\n X \n" << v << "\n = \n" << p2 << endl;
    cout << "That was, a rank " << m.rank() << " tensor times a rank "
	 << v.rank() << " tensor gives a rank " << p2.rank() << " tensor\n";

    cout << "--------" << endl;
}




/*
 * Tensor contraction.
 */
void test_contraction()
{
    cout << "Testing index contraction, rank 2 tensor...\n";

    Tensor<double> m(2, 3);
    
    m(0, 0) = 1;   m(0, 1) = 2;   m(0, 2) = 3;
    m(1, 0) = 4;   m(1, 1) = 5;   m(1, 2) = 6;
    m(2, 0) = 7;   m(2, 1) = 8;   m(2, 2) = 9;

    cout << "Contraction of indices 0, 1 of tensor:\n" << m << "\n = \n"
	 << contract(m, 0, 1) << endl;
    cout << "That was, its trace in matrix language.\n";

    cout << "--------" << endl;

    /*
     * More complicated.
     */
    cout << "Testing index contraction, rank 3 tensor, last indices...\n";

    Tensor<double> t(3, 2);

    t(0, 0, 0) = 1;  t(0, 0, 1) = 2;
    t(0, 1, 0) = 3;  t(0, 1, 1) = 4;

    t(1, 0, 0) = 5;  t(1, 0, 1) = 6;
    t(1, 1, 0) = 7;  t(1, 1, 1) = 8;

    Tensor<double> c = contract(t, 1, 2);

    cout << "Contraction of indices 1, 2 of tensor:\n" << t
	 << "\n = \n" << c << endl;

    cout << "-----------" << endl;

    cout << "Testing matrix product...\n";

    Tensor<double> m1(2, 3);
    Tensor<double> m2(2, 3);

    m1(0, 0) = 1;   m1(0, 1) = 2;   m1(0, 2) = 3;
    m1(1, 0) = 4;   m1(1, 1) = 5;   m1(1, 2) = 6;
    m1(2, 0) = 7;   m1(2, 1) = 8;   m1(2, 2) = 9;

    m2(0, 0) = 10;  m2(0, 1) = 11;  m2(0, 2) = 12;
    m2(1, 0) = 13;  m2(1, 1) = 14;  m2(1, 2) = 15;
    m2(2, 0) = 16;  m2(2, 1) = 17;  m2(2, 2) = 18;

    cout << m1 << "\n * \n" << m2 << "\n = \n";
    cout << contract(m1 * m2, 1, 2) << endl;

    cout << "-----------" << endl;
}


/*
 * Test swapping indices.
 */
void test_swap()
{
    cout << "Testing indices swapping...\n";

    Tensor<double> m(2, 3);

    m(0, 0) = 1;   m(0, 1) = 2;   m(0, 2) = 3;
    m(1, 0) = 4;   m(1, 1) = 5;   m(1, 2) = 6;
    m(2, 0) = 7;   m(2, 1) = 8;   m(2, 2) = 9;

    cout << "swapping indices 0,1 of tensor:\n" << m
	 << "\n = \n" << swap(m, 0, 1) << endl;

    cout << "-----------" << endl;

    /*
     * More complicated.
     */
    cout << "Testing indices swapping, rank 3 tensor...\n";

    Tensor<double> t(3, 2);

    t(0, 0, 0) = 1;  t(0, 0, 1) = 2;
    t(0, 1, 0) = 3;  t(0, 1, 1) = 4;

    t(1, 0, 0) = 5;  t(1, 0, 1) = 6;
    t(1, 1, 0) = 7;  t(1, 1, 1) = 8;

    Tensor<double> s = swap(t, 1, 2);

    cout << "Swap of indices 1, 2 of tensor:\n" << t
	 << "\n = \n" << s << endl;

    cout << "-----------" << endl;
}



/*
 * Access to elements through brackets.
 */
void test_brackets()
{
    Tensor<double> t(3, 2);

    t(0, 0, 0) = 1;  t(0, 0, 1) = 2;
    t(0, 1, 0) = 3;  t(0, 1, 1) = 4;

    t(1, 0, 0) = 5;  t(1, 0, 1) = 6;
    t(1, 1, 0) = 7;  t(1, 1, 1) = 8;

    cout << "Reading tensor with operator (int, ...):\n";

    for (int i1 = 0; i1 < t.dimension(); i1++) {
	for (int i2 = 0; i2 < t.dimension(); i2++) {
	    for (int i3 = 0; i3 < t.dimension(); i3++)
		cout << t(i1, i2, i3) << "\t";

	    cout << endl;
	}
	cout << endl;
    }

    cout << "-----------" << endl;

    Tensor<double> s(0, 3);

    s() = 17;

    cout << "Scalar: " << s() << endl;

    cout << "-----------" << endl;
}



/* --------------------------- */

/*
 * Tests auxiliary functions for index contraction.
 */
void test_base_conversion()
{
    const int NDIGITS = 5;
    const int BASE    = 8;
    const int NUMBER  = 623;

    cout << "Testing base change...\n";

    vector<long> digit = decimal2base(NDIGITS, BASE, NUMBER);

    cout << NUMBER << " in base " << BASE << " = ";
    for (int i = (int) digit.size() - 1; i >= 0; i--) {
	cout << digit[i] << " ";
    }
    
    cout << " ---> " << base2decimal(BASE, digit) << endl;

    cout << "-----------" << endl;
}
