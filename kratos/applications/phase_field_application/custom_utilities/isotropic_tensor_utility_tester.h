//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 11 Sep 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_ISOTROPIC_TENSOR_UTILITY_TESTER_H_INCLUDED )
#define  KRATOS_ISOTROPIC_TENSOR_UTILITY_TESTER_H_INCLUDED



// System includes
#include <iostream>
#include <string>
#include <vector>
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "isotropic_tensor_utility.h"


namespace Kratos
{


class IsotropicTensorUtilityTester
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicTensorUtilityTester);

    IsotropicTensorUtilityTester()
    {
    }

    virtual ~IsotropicTensorUtilityTester()
    {
    }

    template<std::size_t TDimension>
    class TestFunction : public SingleIsotropicTensorFunction<TDimension>
    {
        virtual double Value(const double& x, const double& y)
        {
            return x;
        }

        virtual double Value(const double& x, const double& y, const double& z)
        {
            return x;
        }

        virtual void Derivative(const double& x, const double& y, Vector& Y)
        {
            Y[0] = 1.0;
            Y[1] = 0.0;
        }

        virtual void Derivative(const double& x, const double& y, const double& z, Vector& Y)
        {
            Y[0] = 1.0;
            Y[1] = 0.0;
            Y[2] = 0.0;
        }
    };

    template<std::size_t TDimension>
    class TestFunction2 : public SingleIsotropicTensorFunction<TDimension>
    {
        virtual double Value(const double& x, const double& y)
        {
            if(x > 0)
                return x;
            else
                return 0.0;
        }

        virtual double Value(const double& x, const double& y, const double& z)
        {
            if(x > 0)
                return x;
            else
                return 0.0;
        }

        virtual void Derivative(const double& x, const double& y, Vector& Y)
        {
            if(x > 0)
                Y[0] = 1.0;
            else
                Y[0] = 0.0;
            Y[1] = 0.0;
        }

        virtual void Derivative(const double& x, const double& y, const double& z, Vector& Y)
        {
            if(x > 0)
                Y[0] = 1.0;
            else
                Y[0] = 0.0;
            Y[1] = 0.0;
            Y[2] = 0.0;
        }
    };

    void Test1()
    {
        Matrix X(2, 2);
        X(0, 0) = 1.0;
        X(0, 1) = 2.0;
        X(1, 0) = 2.0;
        X(1, 1) = 3.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<2>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0] + E[1])
        KRATOS_WATCH(e[0]*E[0] + e[1]*E[1])
        KRATOS_WATCH(IsotropicTensorUtility<2>::TwoDotProduct(E[0], E[0]))
        KRATOS_WATCH(IsotropicTensorUtility<2>::TwoDotProduct(E[1], E[1]))
        KRATOS_WATCH(IsotropicTensorUtility<2>::TwoDotProduct(E[0], E[1]))
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction<2> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<2>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<2>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<2>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

    void Test1_1()
    {
        Matrix X(2, 2);
        X(0, 0) = 0.0;
        X(0, 1) = 0.0;
        X(1, 0) = 0.0;
        X(1, 1) = 0.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<2>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0])
        KRATOS_WATCH(e[0]*E[0])
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction<2> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<2>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<2>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<2>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

    void Test1_2()
    {
        Matrix X(2, 2);
        X(0, 0) = 0.0;
        X(0, 1) = 0.0;
        X(1, 0) = 0.0;
        X(1, 1) = 0.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<2>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0])
        KRATOS_WATCH(e[0]*E[0])
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction2<2> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<2>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<2>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<2>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<2>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

    void Test2()
    {
        Matrix X(3, 3);
        X(0, 0) = 1.0;
        X(0, 1) = 2.0;
        X(0, 2) = 3.0;
        X(1, 0) = 2.0;
        X(1, 1) = 4.0;
        X(1, 2) = 5.0;
        X(2, 0) = 3.0;
        X(2, 1) = 5.0;
        X(2, 2) = 6.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<3>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0] + E[1] + E[2])
        KRATOS_WATCH(e[0]*E[0] + e[1]*E[1] + e[2]*E[2])
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[0], E[0]))
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[1], E[1]))
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[2], E[2]))
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[0], E[1]))
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[1], E[2]))
        KRATOS_WATCH(IsotropicTensorUtility<3>::TwoDotProduct(E[0], E[2]))
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction<3> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<3>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<3>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<3>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<3>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

    void Test2_1()
    {
        Matrix X(3, 3);
        X(0, 0) = 0.0;
        X(0, 1) = 0.0;
        X(0, 2) = 0.0;
        X(1, 0) = 0.0;
        X(1, 1) = 0.0;
        X(1, 2) = 0.0;
        X(2, 0) = 0.0;
        X(2, 1) = 0.0;
        X(2, 2) = 0.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<3>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0])
        KRATOS_WATCH(e[0]*E[0])
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction<3> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<3>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<3>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<3>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<3>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

    void Test2_2()
    {
        Matrix X(3, 3); // (-0.229217,8.58688e-17,-2.0383e-17),(8.58688e-17,0,0),(-2.0383e-17,0,0)
        X(0, 0) = -0.229217;
        X(0, 1) = 8.58688e-17;
        X(0, 2) = -2.0383e-17;
        X(1, 0) = 8.58688e-17;
        X(1, 1) = 0.0;
        X(1, 2) = 0.0;
        X(2, 0) = -2.0383e-17;
        X(2, 1) = 0.0;
        X(2, 2) = 0.0;

        Vector e;
        std::vector<Matrix> E;
        SDParameters decomp_params = IsotropicTensorUtility<3>::SpectralDecomposition(X, e, E);

        KRATOS_WATCH(e)

        for(unsigned int i = 0; i < E.size(); ++i)
            KRATOS_WATCH(E[i])

        KRATOS_WATCH(E[0])
        KRATOS_WATCH(e[0]*E[0])
        std::cout << "---------------------------------" << std::endl;

        /////////////////////////////////////////////////////////////
        TestFunction<3> Phi;
        Matrix V;
        Matrix Ce;
        IsotropicTensorUtility<3>::Fourth_Order_Tensor D;
        IsotropicTensorUtility<3>::Value(X, e, E, decomp_params, Phi, V);
        IsotropicTensorUtility<3>::Derivative(X, e, E, decomp_params, Phi, D);
        IsotropicTensorUtility<3>::FourthOrderTensorToMatrix(D, Ce);
        KRATOS_WATCH(V)
        KRATOS_WATCH(D)
        KRATOS_WATCH(Ce)
    }

private:

};

}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_TENSOR_UTILITY_H_INCLUDED  defined 

