//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_utilities/solid_mechanics_math_utilities.hpp


#if !defined(KRATOS_PARTICLE_MECHANICS_MATH_UTILITIES_H_INCLUDED)
#define      KRATOS_PARTICLE_MECHANICS_MATH_UTILITIES_H_INCLUDED


#ifdef FIND_MAX
#undef FIND_MAX
#endif

#define FIND_MAX(a, b) ((a)>(b)?(a):(b))

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "particle_mechanics_application_variables.h"

#if !defined(INITIAL_CURRENT)
#define INITIAL_CURRENT
    enum Configuration {Initial = 0, Current = 1};
#endif

namespace Kratos
{
template<class TDataType>
class ParticleMechanicsMathUtilities
{
public:
    /**
     * @name type definitions
     * @{
     */
    typedef Matrix MatrixType;

    typedef Vector VectorType;

    typedef unsigned int IndexType;

    typedef unsigned int SizeType;

    typedef MathUtils<TDataType> MathUtilsType;

	typedef Node<3> NodeType;

	typedef Geometry< Node<3> > GeometryType;

    typedef DenseVector<Vector> Second_Order_Tensor;

    typedef DenseVector<Second_Order_Tensor> Third_Order_Tensor;

    typedef DenseVector<DenseVector<Matrix> > Fourth_Order_Tensor;

    typedef DenseMatrix<Second_Order_Tensor> Matrix_Second_Tensor;

    /**
     * This function returns rotation matrix from given normal vector and dimension
     * @param rRotationMatrix: Rotation Matrix
     * @param rNormalVector: Normal Vector at the material point condition
     * @param Dimension: given dimension
     */
    static inline void GetRotationMatrix(
        MatrixType& rRotationMatrix,
        const VectorType& rNormalVector,
        const unsigned int Dimension
        )
    {
        if(Dimension == 2){
            if (rRotationMatrix.size1() != 2 && rRotationMatrix.size2() != 2)
                rRotationMatrix.resize(2,2,false);
            noalias(rRotationMatrix) = IdentityMatrix(Dimension);

            double aux = rNormalVector[0]*rNormalVector[0] + rNormalVector[1]*rNormalVector[1];
            aux = std::sqrt(aux);
            if (std::abs(aux) < std::numeric_limits<double>::epsilon()) aux = std::numeric_limits<double>::epsilon();

            rRotationMatrix(0,0) =  rNormalVector[0]/aux;
            rRotationMatrix(0,1) =  rNormalVector[1]/aux;
            rRotationMatrix(1,0) = -rNormalVector[1]/aux;
            rRotationMatrix(1,1) =  rNormalVector[0]/aux;
        }
        else if (Dimension == 3){
            if (rRotationMatrix.size1() != 3 && rRotationMatrix.size2() != 3)
                rRotationMatrix.resize(2,2,false);
            noalias(rRotationMatrix) = IdentityMatrix(Dimension);

            double aux = rNormalVector[0]*rNormalVector[0] + rNormalVector[1]*rNormalVector[1] + rNormalVector[2]*rNormalVector[2];
            aux = std::sqrt(aux);
            if (std::abs(aux) < std::numeric_limits<double>::epsilon()) aux = std::numeric_limits<double>::epsilon();

            rRotationMatrix(0,0) = rNormalVector[0]/aux;
            rRotationMatrix(0,1) = rNormalVector[1]/aux;
            rRotationMatrix(0,2) = rNormalVector[2]/aux;

            // Define the new coordinate system, where the first vector is aligned with the normal
            // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
            Vector rT1 = ZeroVector(3);
            rT1(0) = 1.0;
            rT1(1) = 0.0;
            rT1(2) = 0.0;
            double dot = rRotationMatrix(0,0); //this->Dot(rN,rT1);

            // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
            // If this is the case, repeat the procedure using (0,1,0)
            if ( std::abs(dot) > 0.99 )
            {
                rT1(0) = 0.0;
                rT1(1) = 1.0;
                rT1(2) = 0.0;

                dot = rRotationMatrix(0,1); //this->Dot(rN,rT1);
            }

            // Calculate projection and normalize
            rT1[0] -= dot*rRotationMatrix(0,0);
            rT1[1] -= dot*rRotationMatrix(0,1);
            rT1[2] -= dot*rRotationMatrix(0,2);
            ParticleMechanicsMathUtilities<double>::Normalize(rT1);

            rRotationMatrix(1,0) = rT1[0];
            rRotationMatrix(1,0) = rT1[1];
            rRotationMatrix(1,0) = rT1[2];

            // The third base component is choosen as N x T1, which is normalized by construction
            rRotationMatrix(2,0) = rRotationMatrix(0,1)*rT1[2] - rRotationMatrix(0,2)*rT1[1];
            rRotationMatrix(2,1) = rRotationMatrix(0,2)*rT1[0] - rRotationMatrix(0,0)*rT1[2];
            rRotationMatrix(2,2) = rRotationMatrix(0,0)*rT1[1] - rRotationMatrix(0,1)*rT1[0];
        }
        else{
            KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
        }
    }

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    static inline double CalculateRadius(
        const Matrix& rN,
        GeometryType& Geom,
        const Configuration ThisConfiguration = Current
        )
    {
        double radius = 0.0;

        for (unsigned int iNode = 0; iNode < Geom.size(); iNode++)
        {
            // Displacement from the reference to the current configuration
            if (ThisConfiguration == Current)
            {
                const array_1d<double, 3 >& delta_displacement = Geom[iNode].FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3 >& reference_position = Geom[iNode].Coordinates();

                const array_1d<double, 3 > current_position = reference_position + delta_displacement;
                radius += current_position[0] * rN(0, iNode);
            }
            else
            {
                const array_1d<double, 3 >& reference_position = Geom[iNode].Coordinates();
                radius += reference_position[0] * rN(0, iNode);
            }
        }

        return radius;
    }

    /**
     * Calculates the radius of axisymmetry for a point
     * @param Geom: The geometry studied
     * @return The radius of axisymmetry
     */

    static inline double CalculateRadiusPoint(
        GeometryType& Geom,
        const Configuration ThisConfiguration = Current
        )
    {
        // Displacement from the reference to the current configuration
        if (ThisConfiguration == Current)
        {
            const array_1d<double, 3 >& delta_displacement = Geom[0].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3 >& reference_position = Geom[0].Coordinates();

            const array_1d<double, 3 > current_position = reference_position + delta_displacement;
            return current_position[0];
        }
        else
        {
            const array_1d<double, 3 >& reference_position = Geom[0].Coordinates();
            return reference_position[0];
        }
    }

    /**
     * @brief Calculates the QR Factorization of given square matrix A=QR.
     *        The Factorization is performed using the householder algorithm.
     * @param[in]  A The given square matrix the factorization is to be calculated.
     * @param[out] Q The result matrix Q.
     * @param[out] R The result matrix R.
     */
    static inline void QRFactorization(const MatrixType& A, MatrixType& Q, MatrixType& R)
    {

        if(A.size1()!= A.size2())
            KRATOS_ERROR<<" GIVEN MATRIX IS NOT A SQUARE MATRIX: QRFactorization calculation"<<std::endl;

        // QR Factorization with Householder-Algorithm
        unsigned int dimension= A.size1();

        Vector y(dimension);

        Vector w(dimension);

        R.resize(dimension,dimension,false);

        noalias(R)=ZeroMatrix(dimension, dimension);

        Q.resize(dimension,dimension,false);

        noalias(Q)=ZeroMatrix(dimension, dimension);

        Matrix Help(A.size1(),A.size2());
	    noalias(Help) = A;

        Matrix unity(dimension,dimension);
	    noalias(unity) = ZeroMatrix(dimension, dimension);

        for(unsigned int j=0; j<dimension; j++)
            unity(j,j)=1.0;

        std::vector<Matrix> HelpQ(dimension-1);

        std::vector<Matrix> HelpR(dimension-1);

        for(unsigned int i=0; i< dimension-1; i++)
        {
            HelpQ[i].resize(dimension,dimension,false);
            HelpR[i].resize(dimension,dimension,false);
            noalias(HelpQ[i])= unity;
            noalias(HelpR[i])= ZeroMatrix(dimension, dimension);
        }

        for(unsigned int iteration=0; iteration< dimension-1; iteration++)
        {
            // Vector y
            for(unsigned int i=iteration; i<dimension; i++)
                y[i]= Help(i,iteration);


            // Helpvalue l
            double normy=0.0;

            for(unsigned int i=iteration; i<dimension; i++)
                normy += y[i]*y[i];

            normy= std::sqrt(normy);

            double l= std::sqrt((normy*(normy+std::abs(y(iteration))))/2);

            double k=0.0;

            if(y[iteration] !=0)
                k= - y(iteration)/std::abs(y(iteration))*normy;
            else
                k= -normy;

            for(unsigned int i=iteration; i<dimension; i++)
            {
                double e=0;

                if(i==iteration)
                    e=1;

                w[i]= 1/(2*l)*(y[i]-k*e);
            }

            for(unsigned int i=iteration; i<dimension; i++)
                for(unsigned int j=iteration; j<dimension; j++)
                    HelpQ[iteration](i,j)= unity(i,j)- 2*w[i]*w[j];


            for(unsigned int i=iteration; i<dimension; i++)
                for(unsigned int j=iteration; j<dimension; j++)
                    for(unsigned int k=iteration; k<dimension; k++)
                        HelpR[iteration](i,j)+= HelpQ[iteration](i,k)*Help(k,j);

            Help= HelpR[iteration];

        }

        // Assembling R
        for(unsigned int k=0; k<dimension-1; k++)
        {
            for(unsigned int i=k; i<dimension; i++)
                for(unsigned int j=k; j<dimension; j++)
                    R(i,j) =HelpR[k](i,j);

        }

        for(unsigned int k=1; k<dimension-1; k++)
        {
            for(unsigned int i=0; i<dimension; i++)
                for(unsigned int j=0; j<dimension; j++)
                    for(unsigned int l=0; l<dimension; l++)
                        Q(i,j)+= HelpQ[(k-1)](i,l)*HelpQ[k](l,j);
            noalias(HelpQ[k])=Q;
        }
        if(dimension-1==1)
            noalias(Q)=HelpQ[0];

    }


    /**
     * @brief Calculates Eigenvalues of given square matrix A.
     *        The QR Algorithm with shifts is used.
     * @param[in] A the given square matrix the eigenvalues are to be calculated.
     * @param[in] rTolerance convergence criteria.
     * @param[in] rZeroTolerance number treated as zero.
     * @return Vector of eigenvalues.
     * @warning Only valid for 2*2 and 3*3 Matrices yet.
     */
    static inline Vector EigenValues(const Matrix& A, const double rTolerance = 1e-9, const double rZeroTolerance = 1e-9)
    {
        unsigned int dimension= A.size1();

        Matrix Convergence(2,dimension);
	    noalias(Convergence) = ZeroMatrix(2,dimension);

        double delta = 0.0;

        double abs   = 0.0;

        Vector Result(dimension);
	    noalias(Result) = ZeroVector(dimension);

        Matrix HelpA(dimension,dimension);
	    noalias(HelpA) = ZeroMatrix(dimension, dimension);

        Matrix HelpQ(dimension,dimension);
	    noalias(HelpQ) = ZeroMatrix(dimension, dimension);

        Matrix HelpR(dimension,dimension);
	    noalias(HelpR) = ZeroMatrix(dimension, dimension);

        HelpA=A;

        bool is_converged=false;
        unsigned int max_iters = 10;
        unsigned int iter = 0;

        while(is_converged == false && iter<max_iters )
        {
            double shift= HelpA((dimension-1),(dimension-1));

            for(unsigned int i=0; i<dimension; i++)
            {
                HelpA(i,i) = HelpA(i,i)- shift;
            }

            ParticleMechanicsMathUtilities<double>::QRFactorization(HelpA, HelpQ, HelpR);

            HelpA= ZeroMatrix(dimension, dimension);

            for(unsigned int i=0; i<dimension; i++)
            {
                HelpA(i,i) += shift;
                for(unsigned int j=0; j< dimension; j++)
                {
                    for(unsigned int k=0; k< dimension; k++)
                    {
                        HelpA(i,j) += HelpR(i,k)*HelpQ(k,j);
                    }
                }
            }

            delta= 0.0;

            abs = 0.0;

            for(unsigned int i=0; i<dimension; i++)
            {
                Convergence(0,i)=Convergence(1,i);
                Convergence(1,i)=HelpA(i,i);
                delta+= (Convergence(1,i)-Convergence(0,i))*(Convergence(1,i)-Convergence(0,i));
                abs+=(Convergence(1,i))*(Convergence(1,i));
            }

            delta= std::sqrt(delta);

            abs=std::sqrt(abs);

            if(abs < rZeroTolerance)
                abs=1.0;

            if(delta < rZeroTolerance || (delta/abs) < rTolerance)
                is_converged=true;

            iter++;
        }


        for(unsigned int i=0; i<dimension; i++)
        {
            Result[i]= HelpA(i,i);

            if(std::abs(Result[i]) <rZeroTolerance)
                Result[i]=0.0;
        }

        return Result;
    }


    /**
     * @brief Calculates the Eigenvalues using a direct method.
     * @param A the given square matrix the eigenvalues are to be calculated.
     * @return Vector of eigenvalues.
     * @warning only valid symmetric 3*3 Matrices.
     */
    static inline Vector EigenValuesDirectMethod(const Matrix& A)
    {
        // Given a real symmetric 3x3 matrix A, compute the eigenvalues
        unsigned int dimension= A.size1();
        Vector Result(dimension);
	    noalias(Result) = ZeroVector(dimension);

        const double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        if (p1 == 0)
        {//A is diagonal.
            Result[0] = A(0,0);
            Result[1] = A(1,1);
            Result[2] = A(2,2);
            return Result;
        }

        const double q = (A(0,0) + A(1,1) + A(2,2)) / 3.0;
        const double p2 = (A(0,0) - q) * (A(0,0) - q) + (A(1,1) - q) * (A(1,1) - q) + (A(2,2) - q) * (A(2,2) - q) + 2.0 * p1;
        const double p = std::sqrt(p2 / 6.0);

        Matrix B(3,3);
        const double inv_p = 1.0/p;

        // B = (1 / p) * (A - q * I)  where  I is the identity matrix
        B(0,0) = inv_p * (A(0,0) - q);
        B(1,1) = inv_p * (A(1,1) - q);
        B(2,2) = inv_p * (A(2,2) - q);
        B(0,1) = inv_p * A(0,1);
        B(1,0) = inv_p * A(1,0);
        B(0,2) = inv_p * A(0,2);
        B(2,0) = inv_p * A(2,0);
        B(1,2) = inv_p * A(1,2);
        B(2,1) = inv_p * A(2,1);

        //r = det(B) / 2
        double r = 0.5 * ( B(0,0)*B(1,1)*B(2,2) + B(0,1)*B(1,2)*B(2,0) + B(1,0)*B(2,1)*B(0,2) - B(2,0)*B(1,1)*B(0,2) - B(1,0)*B(0,1)*B(2,2) - B(0,0)*B(2,1)*B(1,2) );

        // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range.
        double phi = 0.0;
        if (r <= -1) { phi = Globals::Pi / 3.0; }
        else if (r >= 1) { phi = 0.0; }
        else { phi = acos(r) / 3.0;}

        // the eigenvalues satisfy eig3 <= eig2 <= eig1
        Result[0] = q + 2.0 * p * cos(phi);
        Result[2] = q + 2.0 * p * cos(phi + (2.0 * Globals::Pi / 3.0));
        Result[1] = 3.0 * q - Result[0] - Result[2];     //% since trace(A) = eig1 + eig2 + eig3

        return Result;
    }


    /**
     * @brief Calculates the Eigenvectors and Eigenvalues of given symmetric matrix A.
     *        The eigenvectors and eigenvalues are calculated using the iterative
     *        Gauss-Seidel-method.
     * @param[in] A the given symmetric matrix where the eigenvectors have to be calculated.
     * @param[out] rEigenVectors Where the eigenvectors will be stored.
     * @param[out] rEigenValues Where the eigenvalues will be stored.
     * @param[in] rZeroTolerance The largest value considered to be zero.
     * @param[in] rMaxIteration Maximum number of iteration allowed.
     */
    static inline void EigenVectors(const MatrixType& A,
				    MatrixType& rEigenVectors,
				    VectorType& rEigenValues,
				    const double rZeroTolerance = 1e-9,
				    const unsigned int rMaxIteration = 10)
    {
        if(A.size1()!=3 || A.size2()!=3)
            KRATOS_ERROR<<" GIVEN MATRIX IS NOT 3x3: Eigenvectors calculation"<<std::endl;

        Matrix Help= A;

        for(unsigned int i=0; i<3; i++)
            for(unsigned int j=0; j<3; j++)
                Help(i,j)= Help(i,j);


        rEigenVectors.resize(Help.size1(),Help.size2(),false);

        rEigenValues.resize(Help.size1(),false);

        Matrix HelpDummy(Help.size1(),Help.size2());

        bool is_converged = false;

        Matrix unity(Help.size1(),Help.size2());
	    noalias(unity) = ZeroMatrix(Help.size1(),Help.size2());

        for(unsigned int i=0; i< Help.size1(); i++)
            unity(i,i)= 1.0;

        Matrix V= unity;

        Matrix VDummy(Help.size1(),Help.size2());

        Matrix Rotation(Help.size1(),Help.size2());


        for(unsigned int iterations=0; iterations<rMaxIteration; iterations++)
        {

            is_converged= true;

            double a= 0.0;

            unsigned int index1= 0;

            unsigned int index2= 1;

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=(i+1); j< Help.size2(); j++)
                {
                    if((std::abs(Help(i,j)) > a ) && (std::abs(Help(i,j)) > rZeroTolerance))
                    {
                        a= std::abs(Help(i,j));

                        index1= i;
                        index2= j;

                        is_converged= false;
                    }
                }
            }

            if(is_converged)
                break;

            // Calculation of Rotationangle

            double gamma= (Help(index2,index2)-Help(index1,index1))/(2*Help(index1,index2));

            double u=1.0;

            if(std::abs(gamma) > rZeroTolerance && std::abs(gamma)< (1/rZeroTolerance))
            {
                u= gamma/std::abs(gamma)*1.0/(std::abs(gamma)+std::sqrt(1.0+gamma*gamma));
            }
            else
            {
                if  (std::abs(gamma)>= (1.0/rZeroTolerance))
                    u= 0.5/gamma;
            }

            double c= 1.0/(std::sqrt(1.0+u*u));

            double s= c*u;

            double teta= s/(1.0+c);

            // Rotation of the Matrix
            HelpDummy= Help;

            HelpDummy(index2,index2)= Help(index2,index2)+u*Help(index1,index2);
            HelpDummy(index1,index1)= Help(index1,index1)-u*Help(index1,index2);
            HelpDummy(index1,index2)= 0.0;
            HelpDummy(index2,index1)= 0.0;

            for(unsigned int i=0; i<Help.size1(); i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    HelpDummy(index2,i)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));
                    HelpDummy(i,index2)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));

                    HelpDummy(index1,i)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                    HelpDummy(i,index1)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                }
            }


            Help= HelpDummy;

            //Calculation of the eigenvectors V
            Rotation =unity;
            Rotation(index2,index1)=-s;
            Rotation(index1,index2)=s;
            Rotation(index1,index1)=c;
            Rotation(index2,index2)=c;

            noalias(VDummy) = ZeroMatrix(Help.size1(), Help.size2());

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=0; j< Help.size1(); j++)
                {
                    for(unsigned int k=0; k< Help.size1(); k++)
                    {
                        VDummy(i,j) += V(i,k)*Rotation(k,j);
                    }
                }
            }
            V= VDummy;

        }

        if(!(is_converged))
        {
            KRATOS_WARNING("ParticleMechanicsMathUtilities")<<"########################################################"<<std::endl;
            KRATOS_WARNING("ParticleMechanicsMathUtilities")<<"rMaxIteration exceed in Jacobi-Seidel-Iteration (eigenvectors)"<<std::endl;
            KRATOS_WARNING("ParticleMechanicsMathUtilities")<<"########################################################"<<std::endl;
        }

        for(unsigned int i=0; i< Help.size1(); i++)
        {
            for(unsigned int j=0; j< Help.size1(); j++)
            {
                rEigenVectors(i,j)= V(j,i);
            }
        }

        for(unsigned int i=0; i<Help.size1(); i++)
            rEigenValues[i]= Help(i,i);

    }


    /**
     * @brief Normalises a vector. Vector is scaled by \f$ V_{norm} = \frac{V}{|V|} \f$
     * @param[in/out] v Vector to be normalized
     */
    template< class TVectorType >
    static inline void Normalize( TVectorType& v )
    {
        if (MathUtilsType::Norm( v ) > std::numeric_limits<double>::epsilon())
            v *= 1.0/(MathUtilsType::Norm( v ));
    }


    /**
    * @brief Builds the norm of a given second order tensor
    * @param[in] rTensor the given second order tensor
    * @return The norm of the given tensor
    */
    static double NormTensor(Matrix& rTensor)
    {
        double result=0.0;
        for(unsigned int i=0; i< rTensor.size1(); i++)
            for(unsigned int j=0; j< rTensor.size2(); j++)
                result += rTensor(i,j)*rTensor(i,j);

        result = std::sqrt(result);

        return result;
    }


    /**
    * @brief Transforms a given fourth order tensor to a corresponing Matrix
    * @param[in] rTensor the given symmetric second order tensor
    * @param[out] rVector the vector
    */
    static inline void TensorToMatrix(const Fourth_Order_Tensor& rTensor, Matrix& rMatrix)
    {
        if (rTensor[0].size()== 3)
        {
            if(rMatrix.size1()!=6 || rMatrix.size2()!=6)
                rMatrix.resize(6,6,false);

            rMatrix(0,0) = rTensor[0][0](0,0);
            rMatrix(0,1) = rTensor[0][0](1,1);
            rMatrix(0,2) = rTensor[0][0](2,2);
            rMatrix(0,3) = rTensor[0][0](0,1);
            rMatrix(0,4) = rTensor[0][0](0,2);
            rMatrix(0,5) = rTensor[0][0](1,2);

            rMatrix(1,0) = rTensor[1][1](0,0);
            rMatrix(1,1) = rTensor[1][1](1,1);
            rMatrix(1,2) = rTensor[1][1](2,2);
            rMatrix(1,3) = rTensor[1][1](0,1);
            rMatrix(1,4) = rTensor[1][1](0,2);
            rMatrix(1,5) = rTensor[1][1](1,2);

            rMatrix(2,0) = rTensor[2][2](0,0);
            rMatrix(2,1) = rTensor[2][2](1,1);
            rMatrix(2,2) = rTensor[2][2](2,2);
            rMatrix(2,3) = rTensor[2][2](0,1);
            rMatrix(2,4) = rTensor[2][2](0,2);
            rMatrix(2,5) = rTensor[2][2](1,2);

            rMatrix(3,0) = rTensor[0][1](0,0);
            rMatrix(3,1) = rTensor[0][1](1,1);
            rMatrix(3,2) = rTensor[0][1](2,2);
            rMatrix(3,3) = rTensor[0][1](0,1);
            rMatrix(3,4) = rTensor[0][1](0,2);
            rMatrix(3,5) = rTensor[0][1](1,2);

            rMatrix(4,0) = rTensor[0][2](0,0);
            rMatrix(4,1) = rTensor[0][2](1,1);
            rMatrix(4,2) = rTensor[0][2](2,2);
            rMatrix(4,3) = rTensor[0][2](0,1);
            rMatrix(4,4) = rTensor[0][2](0,2);
            rMatrix(4,5) = rTensor[0][2](1,2);

            rMatrix(5,0) = rTensor[1][2](0,0);
            rMatrix(5,1) = rTensor[1][2](1,1);
            rMatrix(5,2) = rTensor[1][2](2,2);
            rMatrix(5,3) = rTensor[1][2](0,1);
            rMatrix(5,4) = rTensor[1][2](0,2);
            rMatrix(5,5) = rTensor[1][2](1,2);
        }
        else
        {
            if(rMatrix.size1()!=3 || rMatrix.size2()!=3)
                rMatrix.resize(3,3,false);

            rMatrix(0,0) = rTensor[0][0](0,0);
            rMatrix(0,1) = rTensor[0][0](1,1);
            rMatrix(0,2) = rTensor[0][0](0,1);
            rMatrix(1,0) = rTensor[1][1](0,0);
            rMatrix(1,1) = rTensor[1][1](1,1);
            rMatrix(1,2) = rTensor[1][1](0,1);
            rMatrix(2,0) = rTensor[0][1](0,0);
            rMatrix(2,1) = rTensor[0][1](1,1);
            rMatrix(2,2) = rTensor[0][1](0,1);

        }

    }


    /**
    * @brief Transforms a given 6*6 Matrix to a corresponing 4th order tensor
    * @param[in] rTensor the given Matrix
    * @param[out] rTensor the rTensor
    */
    static void MatrixToTensor(const MatrixType& A, std::vector<std::vector<Matrix> >& rTensor)
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;

        rTensor.resize(3);

        for(unsigned int i=0; i<3; i++)
        {
            rTensor[i].resize(3);
            for(unsigned int j=0; j<3; j++)
            {
                rTensor[i][j].resize(3,3,false);
                noalias(rTensor[i][j])= ZeroMatrix(3,3);
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                    {
                        if(i==j) help1= i;
                        else
                        {
                            if((i==0 && j==1) || (i==1 && j==0)) help1= 3;
                            if((i==1 && j==2) || (i==2 && j==1)) help1= 4;
                            if((i==2 && j==0) || (i==0 && j==2)) help1= 5;
                        }
                        if(k==l)
                        {
                            help2= k;
                            coeff=1.0;
                        }
                        else
                        {
                            coeff=0.5;
                            if((k==0 && l==1) || (k==1 && l==0)) help2= 3;
                            if((k==1 && l==2) || (k==2 && l==1)) help2= 4;
                            if((k==2 && l==0) || (k==0 && l==2)) help2= 5;
                        }

                        rTensor[i][j](k,l)= A(help1,help2)*coeff;
                    }
            }
        }


    }


    /**
    * @brief Transforms a given 6*6 Matrix to a corresponing 4th order tensor
    * @param[in] rTensor the given Matrix
    * @param[out] rTensor the rTensor
    */
    static void MatrixToTensor(const MatrixType& A,array_1d<double, 81>& rTensor)
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;
        std::fill(rTensor.begin(), rTensor.end(), 0.0);
        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                    {
                        if(i==j) help1= i;
                        else
                        {
                            if((i==0 && j==1) || (i==1 && j==0)) help1= 3;
                            if((i==1 && j==2) || (i==2 && j==1)) help1= 4;
                            if((i==2 && j==0) || (i==0 && j==2)) help1= 5;
                        }
                        if(k==l)
                        {
                            help2= k;
                            coeff=1.0;
                        }
                        else
                        {
                            coeff=0.5;
                            if((k==0 && l==1) || (k==1 && l==0)) help2= 3;
                            if((k==1 && l==2) || (k==2 && l==1)) help2= 4;
                            if((k==2 && l==0) || (k==0 && l==2)) help2= 5;
                        }

                        rTensor[i*27+j*9+k*3+l]= A(help1,help2)*coeff;
                    }
            }
        }
    }


    /**
    * @brief Transforms a given 4th order tensor to a corresponing 6*6 Matrix
    * @param[in] rTensor the given Tensor
    * @param[out] rMatrix the Matrix
    */
    static void TensorToMatrix(const std::vector<std::vector<Matrix> >& rTensor, Matrix& rMatrix)
    {
        int help1 = 0;
        int help2 = 0;
        int help3 = 0;
        int help4 = 0;
        double coeff = 1.0;

        if(rMatrix.size1()!=6 || rMatrix.size2()!=6)
            rMatrix.resize(6,6,false);

        for(unsigned int i=0; i<6; i++)
            for(unsigned int j=0; j<6; j++)
            {
                if(i<3)
                {
                    help1= i;
                    help2= i;
                }
                else
                {
                    if(i==3)
                    {
                        help1= 0;
                        help2= 1;
                    }
                    if(i==4)
                    {
                        help1= 1;
                        help2= 2;
                    }
                    if(i==5)
                    {
                        help1= 2;
                        help2= 0;
                    }
                }

                if(j<3)
                {
                    help3= j;
                    help4= j;
                    coeff= 1.0;
                }
                else
                {
                    if(j==3)
                    {
                        help3= 0;
                        help4= 1;
                    }
                    if(j==4)
                    {
                        help3= 1;
                        help4= 2;
                    }
                    if(j==5)
                    {
                        help3= 2;
                        help4= 0;
                    }
                    coeff= 2.0;
                }

                rMatrix(i,j)= rTensor[help1][help2](help3,help4)*coeff;
            }

    }


    /**
     * @brief Transforms a given 4th order Tensor to a corresponing 6*6 Matrix
     * @param[in] rTensor the given Tensor
     * @param[out] Matrix the Matrix
     */
    static void TensorToMatrix( const array_1d<double, 81>& rTensor, Matrix& rMatrix )
    {
        if(rMatrix.size1()!=6 || rMatrix.size2()!=6)
            rMatrix.resize(6,6,false);

        rMatrix(0,0) = rTensor[0];
        rMatrix(0,1) = rTensor[4];
        rMatrix(0,2) = rTensor[8];
        rMatrix(0,3) = 2.0*rTensor[1];
        rMatrix(0,4) = 2.0*rTensor[5];
        rMatrix(0,5) = 2.0*rTensor[6];

        rMatrix(1,0) = rTensor[36];
        rMatrix(1,1) = rTensor[40];
        rMatrix(1,2) = rTensor[44];
        rMatrix(1,3) = 2.0*rTensor[37];
        rMatrix(1,4) = 0.0*rTensor[41];
        rMatrix(1,5) = 0.0*rTensor[42];

        rMatrix(2,0) = rTensor[72];
        rMatrix(2,1) = rTensor[76];
        rMatrix(2,2) = rTensor[80];
        rMatrix(2,3) = 2.0*rTensor[73];
        rMatrix(2,4) = 2.0*rTensor[77];
        rMatrix(2,5) = 2.0*rTensor[78];

        rMatrix(3,0) = rTensor[9];
        rMatrix(3,1) = rTensor[13];
        rMatrix(3,2) = rTensor[18];
        rMatrix(3,3) = 2.0*rTensor[10];
        rMatrix(3,4) = 2.0*rTensor[14];
        rMatrix(3,5) = 2.0*rTensor[15];

        rMatrix(4,0) = rTensor[45];
        rMatrix(4,1) = rTensor[49];
        rMatrix(4,2) = rTensor[53];
        rMatrix(4,3) = 2.0*rTensor[46];
        rMatrix(4,4) = 0.0*rTensor[50];
        rMatrix(4,5) = 0.0*rTensor[51];

        rMatrix(5,0) = rTensor[54];
        rMatrix(5,1) = rTensor[58];
        rMatrix(5,2) = rTensor[62];
        rMatrix(5,3) = 2.0*rTensor[55];
        rMatrix(5,4) = 2.0*rTensor[59];
        rMatrix(5,5) = 2.0*rTensor[60];

    }


private:
};// class ParticleMechanicsMathUtilities
}
#endif /* KRATOS_PARTICLE_MECHANICS_MATH_UTILITIES_H_INCLUDED defined */
