#include "test_symbolic_lib.h"
#include <iomanip>

static constexpr std::size_t nnodes = 4;
static constexpr std::size_t dim = 2;
static constexpr std::size_t blocksize = dim + 2;

using ShapeFun = Vector<nnodes>;
using ShapeFunGrad = Matrix<nnodes, dim>;
using ElementData = ElementDataT<nnodes, dim, blocksize>;

int test_substitute_rho_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_mom_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_tot_ener_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_lhs_2D_OSS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_rhs_2D_OSS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_lhs_2D_ASGS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);
int test_substitute_rhs_2D_ASGS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data);

int main()
{
    int result = 0;

    const auto data = RankineHugoniotQuadData();
    ShapeFunGrad DN_DX;
    ShapeFun N;

    double xi = -sqrt(3)/3;
    double eta = sqrt(3)/3;

    QuadShapeFunctions(N, DN_DX, xi, eta);

    result += test_substitute_rho_proj_2D(N, DN_DX, data);
    result += test_substitute_mom_proj_2D(N, DN_DX, data);
    result += test_substitute_tot_ener_proj_2D(N, DN_DX, data);
    result += test_substitute_lhs_2D_OSS(N, DN_DX, data);
    result += test_substitute_rhs_2D_OSS(N, DN_DX, data);
    result += test_substitute_lhs_2D_ASGS(N, DN_DX, data);
    result += test_substitute_rhs_2D_ASGS(N, DN_DX, data);
    
    return result;
}

int test_substitute_rho_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<nnodes> rho_proj;

    // ------------------------- substitute_rho_proj_2D -------------------------
    //substitute_rho_proj_2D

    Vector<nnodes> expected;
    expected[0] = 0.1725866667;
    expected[1] = 0.04624445796;
    expected[2] = 0.1725866667;
    expected[3] = 0.6441022087;

    std::cout << std::setprecision(10) << "rho_proj:\n" << rho_proj << std::endl;

    return expected.assert_almost_equal(rho_proj);
}


int test_substitute_mom_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<dim * nnodes> mom_proj;
// ------------------------- substitute_mom_proj_2D -------------------------
//substitute_mom_proj_2D
    std::cout << std::setprecision(10) << "mom_proj:\n" << mom_proj << std::endl;
    return 0;
}

int test_substitute_tot_ener_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<nnodes> tot_ener_proj;
// ------------------------- substitute_tot_ener_proj_2D -------------------------
//substitute_tot_ener_proj_2D
    std::cout << std::setprecision(10) << "tot_ener_proj:\n" << tot_ener_proj << std::endl;
    return 0;
}

int test_substitute_lhs_2D_OSS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
// ------------------------- substitute_lhs_2D_OSS -------------------------
//substitute_lhs_2D_OSS
    return 0;
}

int test_substitute_rhs_2D_OSS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<blocksize*nnodes> rRightHandSideBoundedVector;
    const double stab_c1 = 12;
    const double stab_c2 = 2;
    const double stab_c3 = 1;
// ------------------------- substitute_rhs_2D_OSS -------------------------
//substitute_rhs_2D_OSS
    std::cout << std::setprecision(10) << "rRightHandSideBoundedVector (OSS):\n" << rRightHandSideBoundedVector << std::endl;
    return 0;
}

int test_substitute_lhs_2D_ASGS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
// ------------------------- substitute_lhs_2D_ASGS -------------------------
//substitute_lhs_2D_ASGS
    return 0;
}

int test_substitute_rhs_2D_ASGS(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<blocksize*nnodes> rRightHandSideBoundedVector;
    const double stab_c1 = 12;
    const double stab_c2 = 2;
    const double stab_c3 = 1;
// ------------------------- substitute_rhs_2D_ASGS -------------------------
//substitute_rhs_2D_ASGS
    std::cout << std::setprecision(10) << "rRightHandSideBoundedVector (ASGS):\n" << rRightHandSideBoundedVector << std::endl;
    return 0;
}
