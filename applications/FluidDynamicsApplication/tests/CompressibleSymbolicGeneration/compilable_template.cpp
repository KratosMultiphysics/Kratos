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
    expected[0] = 0.122021666664901829;
    expected[1] = 0.0326956070419600897;
    expected[2] = 0.122021666664901829;
    expected[3] = 0.455391059617647198;

    std::cout << "test_substitute_rho_proj_2D:";
    int test = expected.assert_almost_equal(rho_proj);
    if(test == 0)
    {
        std::cout << " OK";
    }
    std::cout << std::endl;
    
    return test;
}


int test_substitute_mom_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<dim * nnodes> mom_proj;

// ------------------------- substitute_mom_proj_2D -------------------------
//substitute_mom_proj_2D

    Vector<dim * nnodes> expected;
    expected[0] = -4394.0308473424238;
    expected[1] = 4394.31412179187009;
    expected[2] = -1177.37701706284429;
    expected[3] = 1177.45292022280978;
    expected[4] = -4394.03084734242384;
    expected[5] = 4394.31412179187009;
    expected[6] = -16398.7463723068468;
    expected[7] = 16399.803566944669;
    
    std::cout << "test_substitute_mom_proj_2D:";
    int test = expected.assert_almost_equal(mom_proj);
    if(test == 0)
    {
        std::cout << " OK";
    }
    std::cout << std::endl;
    
    return test;
}

int test_substitute_tot_ener_proj_2D(ShapeFun const& N, ShapeFunGrad const& DN_DX, ElementData const& data)
{
    Vector<nnodes> tot_ener_proj;

// ------------------------- substitute_tot_ener_proj_2D -------------------------
//substitute_tot_ener_proj_2D
    
    Vector<nnodes> expected;
    expected[0] = 0.0944634901631313112;
    expected[1] = 0.0253114159034363399;
    expected[2] = 0.0944634901631313112;
    expected[3] = 0.352542544749088915;

    std::cout << "test_substitute_tot_ener_proj_2D:";
    int test = expected.assert_almost_equal(tot_ener_proj);
    if(test == 0)
    {
        std::cout << " OK";
    }
    std::cout << std::endl;
    
    return test;
}

int test_substitute_lhs_2D_OSS(ShapeFun const&, ShapeFunGrad const&, ElementData const&)
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

    Vector<blocksize*nnodes> expected;
    expected[0] = 0.05903093967;
    expected[1] = 0.04992055143;
    expected[2] = -0.0305187576;
    expected[3] = 0.08084174258;
    expected[4] = 0.08200115221;
    expected[5] = 0.03082145388;
    expected[6] = -0.005004484361;
    expected[7] = 0.1144782285;
    expected[8] = 0.1782081785;
    expected[9] = 0.1268689995;
    expected[10] = -0.03037472342;
    expected[11] = 0.2502758974;
    expected[12] = -0.2567402704;
    expected[13] = 0.230500113;
    expected[14] = -0.1575540898;
    expected[15] = -0.3587281757;

    std::cout << "test_substitute_rhs_2D_OSS:";
    int test = expected.assert_almost_equal(rRightHandSideBoundedVector);
    if(test == 0)
    {
        std::cout << " OK";
    }
    std::cout << std::endl;
    
    return test;
}

int test_substitute_lhs_2D_ASGS(ShapeFun const&, ShapeFunGrad const&, ElementData const&)
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

    Vector<blocksize*nnodes> expected;
    expected[0] = 0.05903093967;
    expected[1] = 0.04992055143;
    expected[2] = -0.0305187576;
    expected[3] = 0.08084174258;
    expected[4] = 0.08200115221;
    expected[5] = 0.03082145388;
    expected[6] = -0.005004484361;
    expected[7] = 0.1144782285;
    expected[8] = 0.1782081785;
    expected[9] = 0.1268689995;
    expected[10] = -0.03037472342;
    expected[11] = 0.2502758974;
    expected[12] = -0.2567402704;
    expected[13] = 0.230500113;
    expected[14] = -0.1575540898;
    expected[15] = -0.3587281757;

    std::cout << "test_substitute_rhs_2D_ASGS:";
    int test = expected.assert_almost_equal(rRightHandSideBoundedVector);
    if(test == 0)
    {
        std::cout << " OK";
    }
    std::cout << std::endl;
    
    return test;
}
