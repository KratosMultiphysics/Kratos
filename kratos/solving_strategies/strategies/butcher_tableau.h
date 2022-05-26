//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

#if !defined(KRATOS_RUNGE_KUTTA_BUTCHER_TABLEAU_H)
#define KRATOS_RUNGE_KUTTA_BUTCHER_TABLEAU_H

/* System includes */

/* External includes */

/* Project includes */
#include "containers/array_1d.h"

namespace Kratos
{

///@name Kratos Globals
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{


///@}
///@name  Functions
///@{


///@}
///@name Kratos Classes
///@{


/** @brief Butcher tableau for Runge-Kutta method.
 *
 * Contains all info necessary of a particular RK method.
 * It specifies the coefficients of the particular Runge-Kutta method.
 *
 * @tparam: Derived:       A derived class that must contain the methods specified below.
 * @tparam: TOrder:        The order of integration.
 * @tparam: TSubstepCount: The number of sub-steps.
 *
 * Child class must provide methods:
 * - static const MatrixType GenerateRKMatrix()     -> Provides coefficients a_ij
 * - static const VectorType GenerateWeights()      -> Provides coefficients b_i
 * - static const VectorType GenerateThetasVector() -> Provides coefficients c_i
 * - static std::string Name()
 */
template<typename Derived, unsigned int TOrder, unsigned int TSubstepCount>
class ButcherTableau
{

public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double, TSubstepCount> VectorType;

    /* Using the following constructs allows us to multiply parts of vectors with parts of matrices
     * while avoiding BOOST's size checks. This is useful to skip multiplications by zero, since
     * for all explicit runge-kutta methods a_ij = 0 for i>j
     */

    typedef array_1d<double,TSubstepCount-1> RowType;
    typedef std::array<RowType, TSubstepCount-1> MatrixType;

    struct ArraySlice
    {
        typename RowType::const_iterator begin;
        typename RowType::const_iterator end;
    };

    static constexpr unsigned int Order() {return TOrder;}
    static constexpr unsigned int SubstepCount() {return TSubstepCount; }

    ///@}
    ///@name Life Cycle
    ///@{


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * The runge kutta must perform for all substeps 1...N-1:
     *
     *  du^(i) = dt * A_ij*k_j
     *
     * This method return the coefficients A_i[1...i]. The rest of coefficents
     * A_i[i+1...n] are skipped. This is they are always zero for explicit
     * Runge-Kutta.
     *
     * @param SubstepIndex: The i in the formula (the row of the matrix). Note that it counts from 1 to n.
     * @param rK: The k in the formula (the reaction)
     *
     * @return: A struct with iterators pointing to A_i1 (.begin) and A_ii (.end)
     *          intended to be used with std::inner_product.
     */

    ArraySlice GetMatrixRow(const unsigned int SubStepIndex) const
    {
        return ArraySlice{
            mA[SubStepIndex - 1].begin(),
            mA[SubStepIndex - 1].begin() + SubStepIndex
        };
    }

    constexpr const VectorType& GetWeights() const
    {
        return mB;
    }

    constexpr double GetIntegrationTheta(const unsigned int SubStepIndex) const
    {
        return mC(SubStepIndex - 1);
    }

    ///@}
    ///@name Input and output
    ///@{

    static std::string Name()
    {
        return Derived::Name();
    }

    virtual std::string Info() const = 0;

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    const MatrixType mA = Derived::GenerateRKMatrix();     // Runge-Kutta matrix
    const VectorType mB = Derived::GenerateWeights();      // Weights vector
    const VectorType mC = Derived::GenerateThetasVector(); // Nodes vector

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

   ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

private:

    ///@}
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{
};

///@}

///@name Type Definitions
///@{



class ButcherTableauForwardEuler : public ButcherTableau<ButcherTableauForwardEuler, 1, 1>
{
public:
    typedef ButcherTableau<ButcherTableauForwardEuler, 1, 1> BaseType;

    static const BaseType::MatrixType GenerateRKMatrix()
    {
        return {};
    }

    static const BaseType::VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0;
        return B;
    }

    static const BaseType::VectorType GenerateThetasVector()
    {
        VectorType C;
        C(0) = 0.0;
        return C;
    }

    static std::string Name()
    {
        return "butcher_tableau_forward_euler";
    }

    std::string Info() const override
    {
        return "ButcherTableauForwardEuler";
    }
};


class ButcherTableauMidPointMethod : public ButcherTableau<ButcherTableauMidPointMethod, 2, 2>
{
public:
    typedef ButcherTableau<ButcherTableauMidPointMethod, 2, 2> BaseType;

    static const BaseType::MatrixType GenerateRKMatrix()
    {
        MatrixType A;
        A[0][0] = 0.5;
        return A;
    }

    static const BaseType::VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 0.0;
        B(1) = 1.0;
        return B;
    }

    static const BaseType::VectorType GenerateThetasVector()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 0.5;
        return C;
    }

    static std::string Name()
    {
        return "butcher_tableau_midpoint_method";
    }

    std::string Info() const override
    {
        return "ButcherTableauMidPointMethod";
    }
};

/** @brief Explicit total variation diminishing 3rd order Runge-Kutta
 *
 * @details Implementation of Proposition 3.2 of:
 * Gottlieb, Sigal, and Chi-Wang Shu. "Total variation diminishing Runge-Kutta schemes."
 * Mathematics of computation 67.221 (1998): 73-85.
 */
class ButcherTableauRK3TVD : public ButcherTableau<ButcherTableauRK3TVD, 3, 3>
{
public:
    typedef ButcherTableau<ButcherTableauRK3TVD, 3, 3> BaseType;

    static const BaseType::MatrixType GenerateRKMatrix()
    {
        MatrixType A;
        A[0][0] = 1;
        A[1][0] = 0.25;
        A[0][1] = 0.0;
        A[1][1] = 0.25;
        return A;
    }

    static const BaseType::VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0 / 6.0;
        B(1) = 1.0 / 6.0;
        B(2) = 2.0 / 3.0;
        return B;
    }

    static const BaseType::VectorType GenerateThetasVector()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 1.0;
        C(2) = 0.5;
        return C;
    }

    static std::string Name()
    {
        return "butcher_tableau_RK3TVD";
    }

    std::string Info() const override
    {
        return "ButcherTableauRK3TVD";
    }
};


class ButcherTableauRK4 : public ButcherTableau<ButcherTableauRK4, 4, 4>
{
public:
    typedef ButcherTableau<ButcherTableauRK4, 4, 4> BaseType;
    static const BaseType::MatrixType GenerateRKMatrix()
    {
        MatrixType A;
        std::fill(begin(A), end(A), ZeroVector(SubstepCount()-1));
        A[0][0] = 0.5;
        A[1][1] = 0.5;
        A[2][2] = 1.0;
        return A;
    }

    static const BaseType::VectorType GenerateWeights()
    {
        VectorType B;
        B(0) = 1.0 / 6.0;
        B(1) = 1.0 / 3.0;
        B(2) = 1.0 / 3.0;
        B(3) = 1.0 / 6.0;
        return B;
    }

    static const BaseType::VectorType GenerateThetasVector()
    {
        VectorType C;
        C(0) = 0.0;
        C(1) = 0.5;
        C(2) = 0.5;
        C(3) = 1.0;
        return C;
    }

    static std::string Name()
    {
        return "butcher_tableau_RK4";
    }

    std::string Info() const override
    {
        return "ButcherTableauRK4";
    }
};

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RUNGE_KUTTA_BUTCHER_TABLEAU_H  defined */
