//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(MATH_HELPERS_H_INCLUDED)
#define MATH_HELPERS_H_INCLUDED

#include <iostream>
#include <iomanip>

const double Pi = 3.1415926535897932384626433832795;

const double CoeffToRadians = 3.1415926535897932384626433832795 / 180.0;

const double CoeffToDegrees = 180.0 / 3.1415926535897932384626433832795;


namespace Kratos
{


	class MathHelpers
	{

	private:

		MathHelpers() {}
		~MathHelpers() {}
		MathHelpers(const MathHelpers & other);
		MathHelpers & operator = (const MathHelpers & other);

	public:


	public:

		inline static double ToRadians( double degrees ) { return degrees * CoeffToRadians; }

		inline static double ToDegrees( double radians ) { return radians * CoeffToDegrees; }

	public:

		template< class TVector, class TFormat >
		inline static std::string VectorToString(const TVector & v, std::streamsize prec = 4, TFormat fmt = std::scientific)
		{
			std::stringstream ss;
            for(size_t i = 0; i < v.size(); i++)
				ss << fmt << std::setprecision(prec) << v[i] << ", ";
			ss << std::endl;
			return ss.str();
		}

		template< class TMatrix, class TFormat >
		inline static std::string MatrixToString(const TMatrix & m, std::streamsize prec = 4, TFormat fmt = std::scientific)
		{
			std::stringstream ss;
            for(size_t i = 0; i< m.size1(); i++){
                for(size_t j=0; j<m.size2(); j++) {
					ss << fmt << std::setprecision(prec) << m(i,j) << ", ";
				}
				ss << std::endl;
			}
			ss << std::endl;
			return ss.str();
		}

	public:

		template< class TVector3 >
		inline static double LengthSquaredOfVector3(const TVector3 & V)
		{
			return V(0) * V(0) + V(1) * V(1) + V(2) * V(2);
		}

		template< class TVector3 >
		inline static double NormalizeVector3(TVector3 & V)
		{
			double vn( std::sqrt( V(0) * V(0) + V(1) * V(1) + V(2) * V(2) ) );
			V *= (1.0 / vn);
			return vn;
		}

		template< class TVector3 >
		inline static bool TryNormalizeVector3(TVector3 & V)
		{
			double vn( V(0)*V(0) + V(1)*V(1) + V(2)*V(2) );
			if(vn == 0.0) return false;
			if(vn == 1.0) return true;
			vn = std::sqrt(vn);
			V(0) /= vn;
			V(1) /= vn;
			V(2) /= vn;
			return true;
		}

		template< class TVec1, class TVec2, class TVec3>
		inline static void Cross3(const TVec1& a, const TVec2& b, TVec3& c)
		{
			c(0) = a(1)*b(2) - a(2)*b(1);
			c(1) = a(2)*b(0) - a(0)*b(2);
			c(2) = a(0)*b(1) - a(1)*b(0);
		}

	public:

		/**
		Returns true if the two vectors are oriented counter-clock-wise.
		The TVector3 type is the template parameter. It should act like a C-Style vector with 3 components.
		Here it is assumed that the 2 input vectors are already normalized!
		*/
		template<class TVector3>
		inline static bool AreVectorsCounterClockWise(const TVector3& a, const TVector3& b)
		{
			return (a(1)*b(2)-b(1)*a(2) - (a(0)*b(2)-b(0)*a(2)) + a(0)*b(1)-b(0)*a(1)) > 0.0;
		}

		/**
		Calculates the angle(in radians) between the two vectors.
		The sign of the angle is choosen based on the orientation of the two vectors using the
		function AreVectorsCounterClockWise(a,b).
		The TVector3 type is the template parameter. It should act like a C-Style vector with 3 components.
		Here it is assumed that the 2 input vectors are already normalized!
		*/
		template<class TVector3>
		inline static double CalculateAngleBetweenTwoVectors(const TVector3& a, const TVector3& b)
		{
			double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
			if(a_dot_b > 1.0)
				a_dot_b = 1.0;
			else if(a_dot_b < -1.0)
				a_dot_b = -1.0;
			double angle = std::acos( a_dot_b );
			if(!AreVectorsCounterClockWise(a, b))
				angle = -angle;
			return angle;
		}

	public:

		template< class TVec, class TMat>
		inline static void Spin(const TVec & V, TMat & S)
		{
			S(0, 0) =	0.00;		S(0, 1) = - V(2);		S(0, 2) =   V(1);
			S(1, 0) =	V(2);		S(1, 1) =   0.00;		S(1, 2) = - V(0);
			S(2, 0) = - V(1);		S(2, 1) =   V(0);		S(2, 2) =   0.00;
		}

		template< class TVec, class TMat>
		inline static void Spin_AtRow(const TVec & V, TMat & S, size_t row_index)
		{
			size_t i0 = row_index;
			size_t i1 = 1 + row_index;
			size_t i2 = 2 + row_index;
			double v0 = V(i0);
			double v1 = V(i1);
			double v2 = V(i2);
			S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
			S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
			S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
		}

		template< class TVec, class TMat>
		inline static void Spin_AtRow(const TVec & V, TMat & S, size_t vector_index, size_t matrix_row_index)
		{
			size_t i0 = matrix_row_index;
			size_t i1 = 1 + matrix_row_index;
			size_t i2 = 2 + matrix_row_index;
			double v0 = V(vector_index);
			double v1 = V(vector_index + 1);
			double v2 = V(vector_index + 2);
			S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
			S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
			S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
		}

		template< class TVec, class TMat>
		inline static void Spin(const TVec & V, TMat & S, double mult)
		{
			S(0, 0) =	0.00;			S(0, 1) = - mult * V(2);	S(0, 2) =   mult * V(1);
			S(1, 0) =	mult * V(2);	S(1, 1) =   0.00;			S(1, 2) = - mult * V(0);
			S(2, 0) = - mult * V(1);	S(2, 1) =   mult * V(0);	S(2, 2) =   0.00;
		}

		template< class TVec, class TMat>
		inline static void Spin_AtRow(const TVec & V, TMat & S, double mult, size_t row_index)
		{
			size_t i0 = row_index;
			size_t i1 = 1 + row_index;
			size_t i2 = 2 + row_index;
			double v0 = mult * V(i0);
			double v1 = mult * V(i1);
			double v2 = mult * V(i2);
			S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
			S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
			S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
		}

		template< class TVec, class TMat>
		inline static void Spin_AtRow(const TVec & V, TMat & S, double mult, size_t vector_index, size_t matrix_row_index)
		{
			size_t i0 = matrix_row_index;
			size_t i1 = 1 + matrix_row_index;
			size_t i2 = 2 + matrix_row_index;
			double v0 = mult * V(vector_index);
			double v1 = mult * V(vector_index + 1);
			double v2 = mult * V(vector_index + 2);
			S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
			S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
			S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
		}

	public:

		template<class TMat>
		inline static double Determinant3(const TMat & a)
		{
			double c1 =  a(1,1)*a(2,2) - a(1,2)*a(2,1);
			double c2 = -a(1,0)*a(2,2) + a(1,2)*a(2,0);
			double c3 =  a(1,0)*a(2,1) - a(1,1)*a(2,0);
			return a(0,0)*c1 + a(0,1)*c2 + a(0,2)*c3;
		}

		template<class TMat1, class TMat2>
		inline static void InvertMatrix3(const TMat1 & a, TMat2 & b)
		{
			b(0,0) =  a(1,1)*a(2,2) - a(1,2)*a(2,1);
			b(1,0) = -a(1,0)*a(2,2) + a(1,2)*a(2,0);
			b(2,0) =  a(1,0)*a(2,1) - a(1,1)*a(2,0);

			b(0,1) = -a(0,1)*a(2,2) + a(0,2)*a(2,1);
			b(1,1) =  a(0,0)*a(2,2) - a(0,2)*a(2,0);
			b(2,1) = -a(0,0)*a(2,1) + a(0,1)*a(2,0);

			b(0,2) =  a(0,1)*a(1,2) - a(0,2)*a(1,1);
			b(1,2) = -a(0,0)*a(1,2) + a(0,2)*a(1,0);
			b(2,2) =  a(0,0)*a(1,1) - a(0,1)*a(1,0);

			double det = a(0,0)*b(0,0) + a(0,1)*b(1,0) + a(0,2)*b(2,0);

			b /= det;
		}

		template<class TMat>
		inline static void EigenValuesOfMatrix3(const TMat & m, double & e1, double & e2, double & e3)
		{
			double m12 = m(0,1);
			double m13 = m(0,2);
			double m23 = m(1,2);
			double p = m12*m12 + m13*m13 + m23*m23;
			if(p == 0.0) 
			{
				e1 = m(0,0);
				e2 = m(1,1);
				e3 = m(2,2);
			}
			else
			{
				double m11 = m(0, 0);
				double m22 = m(1, 1);
				double m33 = m(2, 2);
				double q = (m11 + m22 + m33) / 3.0;
				m11 -= q;
				m22 -= q;
				m33 -= q;
				p = m11*m11 + m22*m22 + m33*m33 + 2.0 * p;
				p = std::sqrt(p / 6.0);
				Matrix B = (1.0 / p) * (m - q * IdentityMatrix(3, 3));
				double r = Determinant3( B ) / 2.0;

				double phi;
				if(r <= -1.0) 
					phi = Pi / 3.0;
			    else if (r >= 1.0)
					phi = 0.0;
				else
					phi = std::acos(r) / 3.0;

				// the eigenvalues satisfy eig3 <= eig2 <= eig1
				e1 = q + 2.0 * p * std::cos(phi);
				e3 = q + 2.0 * p * std::cos(phi + Pi * (2.0/3.0));
				e2 = 3.0 * q - e1 - e3; //since trace(m) = e1 + e2 + e3
			}
		}

	};


}


#endif // MATH_HELPERS_H_INCLUDED
