//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Massimo Petracca
//                    
//

#if !defined(QUATERNION_H_INCLUDED)
#define QUATERNION_H_INCLUDED

#include "includes/global_variables.h"
#include "includes/serializer.h"

namespace Kratos
{

    /** \brief Quaternion
	* A simple class that implements the main features of quaternion algebra
	*/
	template<class T>
	class Quaternion
	{

	public:
		
		///@name Life Cycle
		///@{
		
		/**
		Creates a Zero Quaternion.
		*/
		Quaternion()
			: mX(0.0)
			, mY(0.0)
			, mZ(0.0)
			, mW(0.0)
		{
		}

		/**
		Creates a Quaternion from its coefficients.
		@param w w coefficient
		@param x x coefficient
		@param y y coefficient
		@param z z coefficient
		*/
		Quaternion(T w, T x, T y, T z)
			: mX(x)
			, mY(y)
			, mZ(z)
			, mW(w)
		{
		}
		
		/**
		Creates a Quaternion from another Quaternion.
		@param other the other Quaternion
		*/
		Quaternion(const Quaternion& other)
			: mX(other.mX)
			, mY(other.mY)
			, mZ(other.mZ)
			, mW(other.mW)
		{
		}

		///@}
                
                /// Destructor.
                virtual ~Quaternion(){};
		
	public:
		
		///@name Operators
		///@{
		
		/**
		Copies a Quaternion.
		@param other the other Quaternion
		*/
		Quaternion& operator= (const Quaternion& other)
		{
			if(this != &other) {
				mX = other.mX;
				mY = other.mY;
				mZ = other.mZ;
				mW = other.mW;
			}
			return *this;
		}
                
                
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 

    template <typename TExpressionType, std::size_t TCategory>
    Quaternion& operator = (AMatrix::MatrixExpression<TExpressionType, TCategory> const& Other)
                {
                    mX = Other.expression()[0];
                    mY = Other.expression()[1];
                    mZ = Other.expression()[2];
                    mW = Other.expression()[3];
                    
                    return *this;
                }
#else
                template<class AE>
                BOOST_UBLAS_INLINE
                Quaternion& operator = (const boost::numeric::ublas::vector_expression<AE> &ae)
                {
                    mX = ae()(0);
                    mY = ae()(1);
                    mZ = ae()(2);
                    mW = ae()(3);
                    
                    return *this;
                }
#endif // ifdef KRATOS_USE_AMATRIX
		
		///@}
		
	public:

		///@name Access
		///@{
		
		/**
		Returns the X coefficient of this quaternion.
		@return the X coefficient of this quaternion.
		*/
		inline const T x()const { return mX; }
		
		/**
		Returns the Y coefficient of this quaternion.
		@return the Y coefficient of this quaternion.
		*/
		inline const T y()const { return mY; }
		
		/**
		Returns the Z coefficient of this quaternion.
		@return the Z coefficient of this quaternion.
		*/
		inline const T z()const { return mZ; }
		
		/**
		Returns the W coefficient of this quaternion.
		@return the W coefficient of this quaternion.
		*/
		inline const T w()const { return mW; }
		
		///@}
		
	public:

		///@name Operations
		///@{
		
		/**
		Returns the squared norm of this quaternion.
		x*x + y*y + z*z + w*w
		@return the squared norm of this quaternion.
		*/
		inline const T squaredNorm()const
		{
			return mX*mX + mY*mY + mZ*mZ + mW*mW;
		}
		
		/**
		Returns the norm of this quaternion.
		sqrt(x*x + y*y + z*z + w*w)
		@return the norm of this quaternion.
		*/
		inline const T norm()const
		{
			return std::sqrt( this->squaredNorm() );
		}
		
		/**
		Makes this Quaternion a Unit Quaternion.
		If this Quaternion is already normalized this is a no-op
		*/
		inline void normalize()
		{
			T n = this->squaredNorm();
			if(n > 0.0 && n != 1.0) {
				n = std::sqrt(n);
				mX /= n;
				mY /= n;
				mZ /= n;
				mW /= n;
			}
		}
		
		/**
	    Returns the Conjugate of this Quaternion, which represents the opposite rotation
		@return the Conjugate of this Quaternion
		*/
		inline Quaternion conjugate()const
		{
			return Quaternion(mW, -mX, -mY, -mZ);
		}
		
		/**
		Constructs a Rotation Matrix from this Quaternion.
		The rotation matrix type is the template argument, no check is made on the type of
		this matrix.
		These assumptions are made:
		The matrix should provide an indexed access like m(i, j) where i and j are indices
		from 0 to 2.
		This means that the input matrix is a C-Style 3x3 Matrix.
		All the 9 coefficients are properly set so there's no need to set the matrix to Zero
		before calling this function.
		@param R the output rotation matrix
		*/
		template<class TMatrix3x3>
		inline void ToRotationMatrix(TMatrix3x3& R)const
		{
		        if( R.size1()!=3 || R.size2()!=3 )
		          R.resize(3,3,false);
		  
			R(0, 0) = 2.0 * ( mW*mW + mX*mX - 0.5 );
			R(0, 1) = 2.0 * ( mX*mY - mW*mZ );
			R(0, 2) = 2.0 * ( mX*mZ + mW*mY );

			R(1, 0) = 2.0 * ( mY*mX + mW*mZ );
			R(1, 1) = 2.0 * ( mW*mW + mY*mY - 0.5 );
			R(1, 2) = 2.0 * ( mY*mZ - mX*mW );

			R(2, 0) = 2.0 * ( mZ*mX - mW*mY );
			R(2, 1) = 2.0 * ( mZ*mY + mW*mX );
			R(2, 2) = 2.0 * ( mW*mW + mZ*mZ - 0.5 );
		}
		
		/**
		Constructs the Euler Angles from this Quaternion.
		Euler Angles expresed in Z(-X)Z sequence as in GiD
		@param EA the output rotation matrix
		*/
		inline void ToEulerAngles(array_1d<double, 3>& EA)const
		{
			double test = mW * mW - mX * mX - mY * mY + mZ * mZ;                    
			if (test > (1.0 - 1.0e-6)) { // singularity at north pole                     
				EA[0] = -atan2 (2 * mZ * mW, (mW * mW - mZ * mZ));
				EA[1] = 0.0;
				EA[2] = 0.0;
			}
			else if (test < (-1.0 +  1.0e-6)) { // singularity at south pole
				EA[0] = atan2 (2 * mZ * mW, (mW * mW - mZ * mZ));
				EA[1] = Globals::Pi;
				EA[2] = 0.0;
			}
			else {                    
				EA[0] = atan2((mX * mZ + mY * mW), -(mY * mZ - mX * mW));
				EA[1] = -acos (test);
				EA[2] = atan2((mX * mZ - mY * mW), (mY * mZ + mX * mW));
			}
		}

		/**
		Extracts the Rotation Vector from this Quaternion
		@param rx the output x component if the rotation vector
		@param ry the output y component if the rotation vector
		@param rz the output z component if the rotation vector
		*/
		inline void ToRotationVector(T& rx, T& ry, T& rz)const
		{
			T xx, yy, zz, ww;
			
			if(mW < 0.0) {
				xx = -mX;
				yy = -mY;
				zz = -mZ;
				ww = -mW;
			}
			else {
				xx = mX;
				yy = mY;
				zz = mZ;
				ww = mW;
			}

			T vNorm = xx*xx + yy*yy + zz*zz;
			if(vNorm == 0.0) {
				rx = 0.0;
				ry = 0.0;
				rz = 0.0;
				return;
			}

			if(vNorm != 1.0)
				vNorm = std::sqrt(vNorm);

			T mult = (vNorm < ww) ? (2.0 / vNorm * std::asin(vNorm)) : (2.0 / vNorm * std::acos(ww));

			rx = xx * mult;
			ry = yy * mult;
			rz = zz * mult;
		}
		
		/**
		Extracts the Rotation Vector from this Quaternion
		The vector type is the template parameter. No check is made on this type.
		The following assumptions are made:
		The vector type should provide indexing like vector(i) where i goes from 0 to 2.
		(i.e. a C-Style vector of size 3)
		@param v the output rotation vector
		*/
		template<class TVector3>
		inline void ToRotationVector(TVector3& v)const
		{
		        if( v.size()!=3 )
			  v.resize(3);

			this->ToRotationVector(v(0), v(1), v(2));
		}

		/**
		Rotates a vector using this quaternion.
		Note: this is faster than constructing the rotation matrix and perform the matrix
		multiplication for a single vector.
		The vector type is the template parameter. No check is made on this type.
		The following assumptions are made:
		The vector type should provide indexing like vector(i) where i goes from 0 to 2.
		(i.e. a C-Style vector of size 3)
		@param a the input source vector
		@param b the output rotated vector
		*/
		template<class TVector3_A, class TVector3_B>
		inline void RotateVector3(const TVector3_A& a, TVector3_B& b)const
		{
			// b = 2.0 * cross( this->VectorialPart, a )
			b(0) = 2.0 * (mY * a(2) - mZ * a(1));
			b(1) = 2.0 * (mZ * a(0) - mX * a(2));
			b(2) = 2.0 * (mX * a(1) - mY * a(0));

			// c = cross( this->VectorialPart, b )
			T c0 = mY * b(2) - mZ * b(1);
			T c1 = mZ * b(0) - mX * b(2);
			T c2 = mX * b(1) - mY * b(0);

			// set results
			b(0) = a(0) + b(0)*mW + c0;
			b(1) = a(1) + b(1)*mW + c1;
			b(2) = a(2) + b(2)*mW + c2;
		}

		/**
		Rotates a vector using this quaternion.
		Note: this is faster than constructing the rotation matrix and perform the matrix
		multiplication for a single vector.
		The vector type is the template parameter. No check is made on this type.
		The following assumptions are made:
		The vector type should provide indexing like vector(i) where i goes from 0 to 2.
		(i.e. a C-Style vector of size 3)
		@param a the input source vector - rotated on exit
		*/
		template<class TVector3>
		inline void RotateVector3(TVector3& a)const
		{
			// b = 2.0 * cross( this->VectorialPart, a )
			T b0 = 2.0 * (mY * a(2) - mZ * a(1));
			T b1 = 2.0 * (mZ * a(0) - mX * a(2));
			T b2 = 2.0 * (mX * a(1) - mY * a(0));

			// c = cross( this->VectorialPart, b )
			T c0 = mY * b2 - mZ * b1;
			T c1 = mZ * b0 - mX * b2;
			T c2 = mX * b1 - mY * b0;

			// set results
			a(0) += b0*mW + c0;
			a(1) += b1*mW + c1;
			a(2) += b2*mW + c2;
		}

		///@}
		
	public:
		
		///@name Static Operations
		///@{
		
		/**
		Returns the Identity Quaternion (i.e. a Quaternion that represents a Zero rotation)
		@return the Identity Quaternion
		*/
		static inline Quaternion Identity()
		{
			return Quaternion(1.0, 0.0, 0.0, 0.0);
		}
		
		/**
		Returns a Quaternion that represents a rotation of an angle 'radians' around the axis (x, y, z)
		@param x the x component of the rotation axis
		@param y the y component of the rotation axis
		@param z the z component of the rotation axis
		@param radians the rotation angle in radians
		@return a Quaternion that represents a rotation of an angle 'radians' around the axis (x, y, z)
		*/
		static inline Quaternion FromAxisAngle(T x, T y, T z, T radians)
		{
			T sqnorm = x*x + y*y + z*z;
			if (sqnorm == 0.0)
				return Quaternion::Identity();

			if(sqnorm != 1.0) {
				T norm = std::sqrt(sqnorm);
				x /= norm;
				y /= norm;
				z /= norm;
			}
			
			T halfAngle = radians * 0.5;

			T s = std::sin(halfAngle);
			T q0 = std::cos(halfAngle);

			Quaternion result(q0, s*x, s*y, s*z);
			result.normalize();

			return result;
		}

		/**
		Returns a Quaternion from a rotation vector
		@param rx the x component of the source rotation vector
		@param ry the y component of the source rotation vector
		@param rz the z component of the source rotation vector
		@return a Quaternion from a rotation vector
		*/
		static inline Quaternion FromRotationVector(T rx, T ry, T rz)
		{
			T rModulus = rx*rx + ry*ry + rz*rz;
			if(rModulus == 0.0)
				return Quaternion::Identity();

			if(rModulus != 1.0) {
				rModulus = std::sqrt(rModulus);
				rx /= rModulus;
				ry /= rModulus;
				rz /= rModulus;
			}

			T halfAngle = rModulus * 0.5;

			T q0 = std::cos(halfAngle);
			T s = std::sin(halfAngle);

			return Quaternion(q0, rx*s, ry*s, rz*s);
		}
		
		/**
		Returns a Quaternion from a rotation vector.
		The vector type is the template parameter. No check is made on this type.
		The following assumptions are made:
		The vector type should provide indexing like vector(i) where i goes from 0 to 2.
		(i.e. a C-Style vector of size 3)
		@param v the source rotation vector
		@return a Quaternion from a rotation vector
		*/
		template<class TVector3>
		static inline Quaternion FromRotationVector(const TVector3& v)
		{
			return Quaternion::FromRotationVector(v(0), v(1), v(2));
		}

		/**
		Returns a Quaternion from a Rotation Matrix.
		The rotation matrix type is the template argument, no check is made on the type of
		this matrix.
		These assumptions are made:
		The matrix should provide an indexed access like m(i, j) where i and j are indices
		from 0 to 2.
		This means that the input matrix is a C-Style 3x3 Matrix.
		@param m the source rotation matrix
		@return a Quaternion from a Rotation Matrix
		*/
		template<class TMatrix3x3>
		static inline Quaternion FromRotationMatrix(const TMatrix3x3 & m)
		{
			T xx = m(0, 0);
			T yy = m(1, 1);
			T zz = m(2, 2);
			T tr = xx + yy + zz;
			Quaternion Q;
			if ((tr > xx) && (tr > yy) && (tr > zz))
			{
				T S = std::sqrt(tr + 1.0) * 2.0;
				Q = Quaternion(
					0.25 * S,
					(m(2, 1) - m(1, 2)) / S,
					(m(0, 2) - m(2, 0)) / S,
					(m(1, 0) - m(0, 1)) / S
					);
			}
			else if ((xx > yy) && (xx > zz))
			{
				T S = std::sqrt(1.0 + xx - yy - zz) * 2.0;
				Q = Quaternion(
					(m(2, 1) - m(1, 2)) / S,
					0.25 * S,
					(m(0, 1) + m(1, 0)) / S,
					(m(0, 2) + m(2, 0)) / S
					);
			}
			else if (yy > zz)
			{
				T S = std::sqrt(1.0 + yy - xx - zz) * 2.0;
				Q = Quaternion(
					(m(0, 2) - m(2, 0)) / S,
					(m(0, 1) + m(1, 0)) / S,
					0.25 * S,
					(m(1, 2) + m(2, 1)) / S
					);
			}
			else
			{
				T S = std::sqrt(1.0 + zz - xx - yy) * 2.0;
				Q =  Quaternion(
					(m(1, 0) - m(0, 1)) / S,
					(m(0, 2) + m(2, 0)) / S,
					(m(1, 2) + m(2, 1)) / S,
					0.25 * S
					);
			}


			Q.normalize();
			return Q;
		}
		
		/**
		Returns a Quaternion from Euler Angles.
		Euler Angles expresed in Z(-X)Z sequence as in GiD
		@param EA the source rotation Euler Angles
		@return a Quaternion from a Euler Angles
		*/
		static inline Quaternion FromEulerAngles(const array_1d<double, 3> & EA)
		{
			Quaternion Q;
			double c2 = cos(-EA[1]/2); double c1p3 = cos((EA[0]+EA[2])/2); double c1m3 = cos((EA[0]-EA[2])/2);
			double s2 = sin(-EA[1]/2); double s1p3 = sin((EA[0]+EA[2])/2); double s1m3 = sin((EA[0]-EA[2])/2);
			
			Q = Quaternion(	
					c2 * c1p3,
					s2 * c1m3,
					s2 * s1m3,
					c2 * s1p3);        
			Q.normalize();
			return Q;
		}
                
		///@}
                virtual void PrintInfo(std::ostream& rOStream) const {
                    rOStream << Info();
                }

                /// Print object's data.
                virtual void PrintData(std::ostream& rOStream) const {
                    rOStream << std::endl << this->mX <<"  " << this->mY << "  " << this->mZ << "  " <<this->mW<< std::endl;                                
                }
		
	private:

		///@name Member Variables
		///@{
		
		T mX;
		T mY;
		T mZ;
		T mW;
		
		///@}

		///@name Serialization
		///@{

		friend class Serializer;

		void save(Serializer& rSerializer) const
		{
			rSerializer.save("X", mX);
			rSerializer.save("Y", mY);
			rSerializer.save("Z", mZ);
			rSerializer.save("W", mW);
		}

		void load(Serializer& rSerializer)
		{
			rSerializer.load("X", mX);
			rSerializer.load("Y", mY);
			rSerializer.load("Z", mZ);
			rSerializer.load("W", mW);
		}		
                
                virtual std::string Info() const {
                    std::stringstream buffer;
                    buffer << "Quaternion " ;
                    return buffer.str();
                }                                
		
	};

	/**
	Performs a Quaternion-Quaternion product and returns the concatenation
	of the two input quaternions (i.e. the compound rotation)
	@param a the first quaternion
	@param b the second quaternion
	@return the compound quaternion
	*/
	template<class T>
	inline Quaternion<T> operator* (const Quaternion<T>& a, const Quaternion<T>& b)
	{
		return Quaternion<T>(
			a.w() * b.w() - a.x() * b.x() - a.y() * b.y() - a.z() * b.z(),
			a.w() * b.x() + a.x() * b.w() + a.y() * b.z() - a.z() * b.y(),
			a.w() * b.y() + a.y() * b.w() + a.z() * b.x() - a.x() * b.z(),
			a.w() * b.z() + a.z() * b.w() + a.x() * b.y() - a.y() * b.x()
		);
	}    
        
        
        template<class T>
        inline std::istream& operator >> (std::istream& rIStream, Quaternion<T>& rThis);

        template<class T>
        inline std::ostream& operator << (std::ostream& rOStream, const Quaternion<T>& rThis)
        {
            rThis.PrintInfo(rOStream);
            rOStream << " : ";
            rThis.PrintData(rOStream);

            return rOStream;
        }
        
} // namespace Kratos

#endif // QUATERNION_H_INCLUDED
