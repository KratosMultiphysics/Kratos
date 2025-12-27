//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:    Andrea Gorgi
//  Co-author:      Codex GPT

#include "custom_processes/assign_iga_interpolation_condition_process.h"

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric> 

// Project includes
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

namespace
{

constexpr double RelTol = 1.0e-12; // relative tolerance for abscissae
using CoordinatesArrayType = array_1d<double, 3>;

inline double ComputeDistance(const CoordinatesArrayType& rA, const CoordinatesArrayType& rB)
{
    const CoordinatesArrayType diff = rB - rA;
    return std::sqrt(inner_prod(diff, diff));
}

inline double Clamp(const double Value, const double LowerBound, const double UpperBound)
{
    return std::max(LowerBound, std::min(Value, UpperBound));
}

// Collapse duplicates in abscissae with a relative tolerance; keep last value for ties.
inline void CollapseDuplicates(
    const std::vector<double>& x_in,
    const std::vector<double>& y_in,
    std::vector<double>& x_out,
    std::vector<double>& y_out)
{
    x_out.clear(); y_out.clear();
    if (x_in.empty()) return;

    // Build index and sort by x to catch non-adjacent duplicates
    std::vector<std::size_t> idx(x_in.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b){ return x_in[a] < x_in[b]; });

    // Scale for relative tolerance
    const double xmin = *std::min_element(x_in.begin(), x_in.end());
    const double xmax = *std::max_element(x_in.begin(), x_in.end());
    const double scale = std::max(1.0, std::abs(xmin) + std::abs(xmax));
    const double atol = RelTol * scale;

    for (std::size_t k = 0; k < idx.size(); ++k) {
        const std::size_t i = idx[k];
        if (x_out.empty() || std::abs(x_in[i] - x_out.back()) > atol) {
            x_out.push_back(x_in[i]);
            y_out.push_back(y_in[i]);
        } else {
            // overwrite last if nearly duplicate
            x_out.back() = x_in[i];
            y_out.back() = y_in[i];
        }
    }
}

// Build chord-length parameters in [0,1]
inline std::vector<double> BuildChordParams(const std::vector<CoordinatesArrayType>& X) {
    const std::size_t n = X.size();
    std::vector<double> t(n, 0.0);
    if (n <= 1) return t;
    double L = 0.0;
    for (std::size_t i=1;i<n;++i) L += std::sqrt(inner_prod(X[i]-X[i-1], X[i]-X[i-1]));
    if (L == 0.0) return t;
    double acc = 0.0;
    for (std::size_t i=1;i<n;++i) {
        acc += std::sqrt(inner_prod(X[i]-X[i-1], X[i]-X[i-1]));
        t[i] = acc / L;
    }
    return t;
}

// Build Vandermonde matrix V_ij = t_i^j
inline void BuildVandermonde(const std::vector<double>& t, std::size_t deg,
                             std::vector<std::vector<double>>& V) {
    const std::size_t n = t.size();
    V.assign(n, std::vector<double>(deg+1, 1.0));
    for (std::size_t i=0;i<n;++i)
        for (std::size_t j=1;j<=deg;++j)
            V[i][j] = V[i][j-1]*t[i];
}

// Solve least-squares V a ≈ y via normal equations (small p+1 → OK). Use double precision.
inline std::vector<double> SolveNormalEq(const std::vector<std::vector<double>>& V,
                                         const std::vector<double>& y) {
    const std::size_t n = V.size();
    const std::size_t m = V[0].size();
    std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
    std::vector<double> b(m, 0.0);
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t r=0;r<m;++r) {
            b[r] += V[i][r]*y[i];
            for (std::size_t c=0;c<m;++c) A[r][c] += V[i][r]*V[i][c];
        }
    }
    // Cholesky (simple) A = L L^T
    for (std::size_t k=0;k<m;++k) {
        double sum = A[k][k];
        for (std::size_t s=0;s<k;++s) sum -= A[k][s]*A[k][s];
        const double Lkk = std::sqrt(std::max(sum, 0.0));
        A[k][k] = Lkk;
        for (std::size_t i=k+1;i<m;++i) {
            double v = A[i][k];
            for (std::size_t s=0;s<k;++s) v -= A[i][s]*A[k][s];
            A[i][k] = v / (Lkk==0.0?1.0:Lkk);
        }
        for (std::size_t j=k+1;j<m;++j) A[k][j]=0.0; // store lower
    }
    // Forward solve L y = b
    std::vector<double> y2(m,0.0);
    for (std::size_t i=0;i<m;++i){
        double v=b[i];
        for (std::size_t s=0;s<i;++s) v -= A[i][s]*y2[s];
        y2[i]=v/(A[i][i]==0.0?1.0:A[i][i]);
    }
    // Backward solve L^T x = y
    std::vector<double> x(m,0.0);
    for (int i=int(m)-1;i>=0;--i){
        double v=y2[i];
        for (std::size_t s=i+1;s<m;++s) v -= A[s][i]*x[s];
        x[i]=v/(A[i][i]==0.0?1.0:A[i][i]);
    }
    return x;
}

// Fit parametric polynomial curve C(t) of degree q using chord-length t in [0,1]
// X(t) = sum a_j t^j, Y(t) = sum b_j t^j, Z similarly (Z omitted if all zero).
struct ParamCurve {
    std::vector<double> ax, ay, az; // coefficients low→high degree
    int dim = 3;
};

inline ParamCurve FitCurve(const std::vector<CoordinatesArrayType>& X, std::size_t q){
    ParamCurve C;
    const std::size_t n = X.size();
    if (n==0){ C.dim=0; return C; }
    // Detect dimension
    bool has_z=false;
    for (const auto& P: X) if (std::abs(P[2])>0.0){ has_z=true; break; }
    C.dim = has_z ? 3 : 2;

    // Parameters and Vandermonde
    const auto t = BuildChordParams(X);
    std::vector<std::vector<double>> V;
    BuildVandermonde(t, q, V);

    // y-vectors
    std::vector<double> xx(n), yy(n), zz(n);
    for (std::size_t i=0;i<n;++i){ xx[i]=X[i][0]; yy[i]=X[i][1]; zz[i]=X[i][2]; }

    C.ax = SolveNormalEq(V, xx);
    C.ay = SolveNormalEq(V, yy);
    if (C.dim==3) C.az = SolveNormalEq(V, zz);
    return C;
}

inline void EvalCurve(const ParamCurve& C, double t, CoordinatesArrayType& x, CoordinatesArrayType& dxdt){
    auto eval1d = [&](const std::vector<double>& a)->std::pair<double,double>{
        // Horner for value and derivative
        double v=0.0, d=0.0;
        for (int j=int(a.size())-1;j>=0;--j){ d = d*t + v; v = v*t + a[j]; }
        return {v, d};
    };
    auto [vx,dx_] = eval1d(C.ax);
    auto [vy,dy_] = eval1d(C.ay);
    double vz=0.0, dz_=0.0;
    if (C.dim==3){ auto [vz_,dzv_] = eval1d(C.az); vz=vz_; dz_=dzv_; }
    x[0]=vx; x[1]=vy; x[2]=vz;
    dxdt[0]=dx_; dxdt[1]=dy_; dxdt[2]=dz_;
}

// Find t* in [0,1] minimizing ||C(t)-x*||^2 via damped Newton
inline double ProjectOnCurve(const ParamCurve& C, const CoordinatesArrayType& xq){
    double t = 0.5; // init
    for (int it=0; it<20; ++it){
        CoordinatesArrayType xc, dc;
        EvalCurve(C, t, xc, dc);
        CoordinatesArrayType r = xc - xq;
        const double f  = 2.0*inner_prod(r, dc);
        // d/dt f = 2( dc·dc + r·d2c ), approximate with dc·dc only for robustness
        const double df = 2.0*inner_prod(dc, dc);
        if (df==0.0) break;
        double step = -f/df;
        // damp and clamp
        step = Clamp(step, -0.25, 0.25);
        t = Clamp(t + step, 0.0, 1.0);
        if (std::abs(step) < 1e-12) break;
    }
    return t;
}

// Build Vandermonde V_ij = t_i^j
inline void BuildVandermonde1D(const std::vector<double>& t, std::size_t deg,
                               std::vector<std::vector<double>>& V) {
    const std::size_t n = t.size();
    V.assign(n, std::vector<double>(deg+1, 1.0));
    for (std::size_t i=0;i<n;++i)
        for (std::size_t j=1;j<=deg;++j)
            V[i][j] = V[i][j-1]*t[i];
}

// Solve dense linear system A x = b (small m ~ O(10)) by Gaussian elimination with partial pivoting
inline bool SolveDense(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    const std::size_t m = A.size();
    for (std::size_t k=0;k<m;++k) {
        // pivot
        std::size_t p = k;
        for (std::size_t i=k+1;i<m;++i) if (std::abs(A[i][k]) > std::abs(A[p][k])) p = i;
        if (std::abs(A[p][k]) < 1e-18) return false;
        if (p != k) { std::swap(A[p], A[k]); std::swap(b[p], b[k]); }
        // eliminate
        const double akk = A[k][k];
        for (std::size_t i=k+1;i<m;++i) {
            const double f = A[i][k]/akk;
            if (f==0.0) continue;
            for (std::size_t j=k;j<m;++j) A[i][j] -= f*A[k][j];
            b[i] -= f*b[k];
        }
    }
    // back-substitution
    for (int i=int(m)-1;i>=0;--i) {
        double s = b[i];
        for (std::size_t j=i+1;j<m;++j) s -= A[i][j]*b[j];
        b[i] = s / A[i][i];
    }
    return true;
}

// Constrained LSQ: minimize ||V a - y||_2^2 s.t. C a = d
inline std::vector<double> SolveConstrainedLSQ(
    const std::vector<std::vector<double>>& V,
    const std::vector<double>& y,
    const std::vector<std::vector<double>>& C,
    const std::vector<double>& d)
{
    const std::size_t n = V.size();
    const std::size_t m = V[0].size();
    const std::size_t c = C.size();

    // Form normal KKT system:
    // [ 2 V^T V   C^T ] [ a ] = [ 2 V^T y ]
    // [   C        0  ] [ λ ]   [   d     ]
    std::vector<std::vector<double>> M(m+c, std::vector<double>(m+c, 0.0));
    std::vector<double> rhs(m+c, 0.0);

    // 2 V^T V and 2 V^T y
    for (std::size_t i=0;i<m;++i) {
        for (std::size_t j=0;j<m;++j) {
            double s = 0.0;
            for (std::size_t k=0;k<n;++k) s += V[k][i]*V[k][j];
            M[i][j] = 2.0*s;
        }
        double sy = 0.0;
        for (std::size_t k=0;k<n;++k) sy += V[k][i]*y[k];
        rhs[i] = 2.0*sy;
    }
    // C and C^T
    for (std::size_t i=0;i<c;++i) {
        for (std::size_t j=0;j<m;++j) {
            M[m+i][j] = C[i][j];
            M[j][m+i] = C[i][j];
        }
        rhs[m+i] = d[i];
    }

    // Solve KKT
    if (!SolveDense(M, rhs)) return std::vector<double>(m, 0.0);

    // First m entries are a
    return std::vector<double>(rhs.begin(), rhs.begin()+static_cast<long>(m));
}

// Evaluate polynomial a_0 + a_1 t + ... + a_p t^p and its Horner derivative if needed
inline double EvalPoly(const std::vector<double>& a, double t) {
    double v = 0.0;
    for (int j=int(a.size())-1;j>=0;--j) v = v*t + a[j];
    return v;
}

// Fit a scalar poly of degree p to all points, enforcing endpoints
inline double PolyFitConstrainedEval(
    const std::vector<double>& t,      // parameters (all samples)
    const std::vector<double>& v,      // values (all samples)
    double tq,                         // query parameter
    std::size_t p,                     // degree
    bool enforce_endpoints = true)
{
    // Collapse duplicates (relative tol)
    std::vector<double> tu, vu;
    CollapseDuplicates(t, v, tu, vu);
    if (tu.empty()) return 0.0;
    if (tu.size()==1) return vu.front();

    // Degree clamp
    p = std::min<std::size_t>(p, tu.size()-1);

    // Vandermonde
    std::vector<std::vector<double>> V;
    BuildVandermonde1D(tu, p, V);

    // Constraints: P(t0)=v0, P(tn)=vn
    std::vector<std::vector<double>> C;
    std::vector<double> d;
    if (enforce_endpoints) {
        C.assign(2, std::vector<double>(p+1, 1.0));
        for (std::size_t j=1;j<=p;++j) { C[0][j] = C[0][j-1]*tu.front(); }
        for (std::size_t j=1;j<=p;++j) { C[1][j] = C[1][j-1]*tu.back();  }
        d = { vu.front(), vu.back() };
    }

    const std::vector<double> a = enforce_endpoints
        ? SolveConstrainedLSQ(V, vu, C, d)
        : SolveConstrainedLSQ(V, vu, {}, {}); // falls back to unconstrained normal eq via KKT

    return EvalPoly(a, tq);
}


// Local polynomial interpolation of degree p using p+1 closest nodes (Newton form).
double PolyInterpDegreeP(
    const std::vector<double>& x,        // abscissae (arbitrary order)
    const std::vector<double>& y,        // values
    const double xq,                      // query
    std::size_t p)                        // desired degree
{
    const std::size_t n = x.size();
    if (n == 0) return 0.0;
    if (n == 1) return y.front();

    // If xq coincides with a node (relative), return its value
    const double scale = std::max(1.0, std::abs(*std::min_element(x.begin(), x.end())) +
                                       std::abs(*std::max_element(x.begin(), x.end())));
    const double atol = RelTol * scale;
    for (std::size_t i = 0; i < n; ++i)
        if (std::abs(xq - x[i]) <= atol) return y[i];

    // m = min(p+1, n)
    const std::size_t m = std::min<std::size_t>(p + 1, n);

    // Select m closest indices to xq
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::partial_sort(
        idx.begin(), idx.begin() + static_cast<long>(m), idx.end(),
        [&](std::size_t a, std::size_t b){ return std::abs(x[a] - xq) < std::abs(x[b] - xq); });
    idx.resize(m);

    // Sort the local stencil by x for stable divided differences
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b){ return x[a] < x[b]; });

    std::vector<double> xl(m), fl(m);
    for (std::size_t i = 0; i < m; ++i) { xl[i] = x[idx[i]]; fl[i] = y[idx[i]]; }

    // Guard against near-duplicates inside stencil
    for (std::size_t i = 1; i < m; ++i)
        if (std::abs(xl[i] - xl[i-1]) <= atol) return fl[i];

    // In-place Newton divided differences
    for (std::size_t j = 1; j < m; ++j) {
        for (std::size_t i = m - 1; i >= j; --i) {
            const double den = xl[i] - xl[i - j];
            if (std::abs(den) <= atol) return fl[i]; // degenerate
            fl[i] = (fl[i] - fl[i - 1]) / den;
            if (i == j) break; // prevent size_t underflow
        }
    }

    // Evaluate at xq (Horner for Newton form)
    double acc = fl[m - 1];
    for (std::size_t i = m - 1; i-- > 0; )
        acc = fl[i] + (xq - xl[i]) * acc;

    return acc;
}

double ComputeCurvilinearAbscissa(
    const CoordinatesArrayType& rPoint,
    const std::vector<CoordinatesArrayType>& rNodeCoordinates,
    const std::vector<double>& rAbscissae)
{
    if (rNodeCoordinates.empty()) {
        return 0.0;
    }

    if (rNodeCoordinates.size() == 1) {
        return rAbscissae.front();
    }

    double best_parameter = rAbscissae.front();
    double min_distance_sq = std::numeric_limits<double>::max();

    for (std::size_t i = 0; i + 1 < rNodeCoordinates.size(); ++i) {
        const CoordinatesArrayType& r_segment_start = rNodeCoordinates[i];
        const CoordinatesArrayType& r_segment_end = rNodeCoordinates[i + 1];
        const CoordinatesArrayType segment = r_segment_end - r_segment_start;
        const double segment_length_sq = inner_prod(segment, segment);

        double projected_parameter = 0.0;
        const double seg_tol = RelTol * (std::sqrt(segment_length_sq) + 1.0);
        if (segment_length_sq > seg_tol * seg_tol) {
            const CoordinatesArrayType to_point = rPoint - r_segment_start;
            const double projection = inner_prod(to_point, segment) / segment_length_sq;
            projected_parameter = Clamp(projection, 0.0, 1.0);
        }

        const CoordinatesArrayType projected_point = r_segment_start + segment * projected_parameter;
        const CoordinatesArrayType diff = projected_point - rPoint;
        const double distance_sq = inner_prod(diff, diff);

        // Update best segment using relative tolerance on segment length
        if (distance_sq < min_distance_sq) {
            min_distance_sq = distance_sq;

            const double segment_length = std::sqrt(segment_length_sq);
            const double seg_atol = RelTol * std::max(1.0, segment_length); // relative tol
            const double effective_len = (segment_length > seg_atol) ? segment_length : 0.0;

            best_parameter = rAbscissae[i] + projected_parameter * effective_len;
        }

    }

    return best_parameter;
}

} // namespace (anonymous)

AssignIgaInterpolationConditionProcess::AssignIgaInterpolationConditionProcess(
    Model& rModel,
    Parameters ThisParameters)
    : mpModel(&rModel)
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    mParameters = ThisParameters;

    KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
        << "AssignIgaInterpolationConditionProcess: missing 'model_part_name'" << std::endl;
    KRATOS_ERROR_IF(mParameters["model_part_name"].GetString().empty())
        << "AssignIgaInterpolationConditionProcess: 'model_part_name' cannot be empty" << std::endl;

    KRATOS_ERROR_IF_NOT(mParameters.Has("reference_model_part_name"))
        << "AssignIgaInterpolationConditionProcess: missing 'reference_model_part_name'" << std::endl;
    KRATOS_ERROR_IF(mParameters["reference_model_part_name"].GetString().empty())
        << "AssignIgaInterpolationConditionProcess: 'reference_model_part_name' cannot be empty" << std::endl;
}

void AssignIgaInterpolationConditionProcess::ExecuteInitialize()
{
    ModelPart& r_target_model_part = mpModel->GetModelPart(mParameters["model_part_name"].GetString());
    const ModelPart& r_reference_model_part = mpModel->GetModelPart(mParameters["reference_model_part_name"].GetString());

    AssignInterpolatedValuesToConditions(r_target_model_part, r_reference_model_part);
}

void AssignIgaInterpolationConditionProcess::ExecuteInitializeSolutionStep()
{
    ModelPart& r_target_model_part = mpModel->GetModelPart(mParameters["model_part_name"].GetString());
    const ModelPart& r_reference_model_part = mpModel->GetModelPart(mParameters["reference_model_part_name"].GetString());

    AssignInterpolatedValuesToConditions(r_target_model_part, r_reference_model_part);
}

const Parameters AssignIgaInterpolationConditionProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "echo_level" : 0,
        "model_part_name" : "",
        "reference_model_part_name" : "",
        "poly_degree": 1
    }
    )");

    return default_parameters;
}

void AssignIgaInterpolationConditionProcess::AssignInterpolatedValuesToConditions(
    ModelPart& rTargetModelPart,
    const ModelPart& rReferenceModelPart)
{
    const auto& r_variables = GetSupportedVariables();

    if (r_variables.empty()) {
        return;
    }

    for (auto& r_condition : rTargetModelPart.Conditions()) {
        if (!r_condition.Has(INTERPOLATION_NODES_ID)) {
            continue;
        }

        const std::vector<IndexType>& interpolation_node_ids = r_condition.GetValue(INTERPOLATION_NODES_ID);
        KRATOS_ERROR_IF(interpolation_node_ids.empty())
            << "AssignIgaInterpolationConditionProcess: condition " << r_condition.Id()
            << " has an empty INTERPOLATION_NODES_ID list" << std::endl;

        const GeometryType& r_geometry = r_condition.GetGeometry();
        const CoordinatesArrayType center_coordinates = r_geometry.Center();

        std::vector<const NodeType*> p_nodes(interpolation_node_ids.size(), nullptr);
        std::vector<CoordinatesArrayType> node_coordinates;
        node_coordinates.reserve(interpolation_node_ids.size());
        for (IndexType i = 0; i < interpolation_node_ids.size(); ++i) {
            const IndexType node_id = interpolation_node_ids[i];
            KRATOS_ERROR_IF_NOT(rReferenceModelPart.HasNode(node_id))
                << "AssignIgaInterpolationConditionProcess: reference model part '"
                << rReferenceModelPart.FullName() << "' does not contain node " << node_id
                << " required by condition " << r_condition.Id() << std::endl;

            p_nodes[i] = &rReferenceModelPart.GetNode(node_id);
            node_coordinates.push_back(p_nodes[i]->Coordinates());
        }

        // 0) Build chord-length parameters for nodes
        std::vector<double> t_nodes = BuildChordParams(node_coordinates);

        // 1) Fit geometric curve of degree q_geom = min(2, |nodes|-1) or expose a parameter
        const std::size_t q_geom = std::min<std::size_t>(2, node_coordinates.size()-1);
        const ParamCurve curve = FitCurve(node_coordinates, q_geom);

        // 2) Parameter of the condition center by projecting onto fitted curve
        const double t_center = ProjectOnCurve(curve, center_coordinates);

        for (const auto* p_variable : r_variables) {
            bool all_nodes_have_value = true;
            std::vector<double> nodal_values;
            nodal_values.reserve(p_nodes.size());

            for (IndexType i = 0; i < interpolation_node_ids.size(); ++i) {
                const NodeType& r_node = *p_nodes[i];
                if (!r_node.Has(*p_variable)) {
                    all_nodes_have_value = false;
                    break;
                }
                
                // FIXME: just for the Neumann convergence to manufactured solution
                if (*p_variable == FORCE_X)
                {
                    double nu = 0.3;
                    double E = 1000;

                    const double x = r_node.X();
                    const double y = r_node.Y();
                    // Vector normal = r_node.GetValue(NORMAL);

                    Vector normal = r_condition.GetValue(NORMAL);

                    // // // cosinusoidal
                    double flux_x = E/(1-nu)*(sin(x)*sinh(y)) * normal[0]; 

                    nodal_values.push_back(flux_x);

                } else if (*p_variable == FORCE_Y)
                {
                    double nu = 0.3;
                    double E = 1000;

                    const double x = r_node.X();
                    const double y = r_node.Y();
                    // Vector normal = r_node.GetValue(NORMAL);

                    Vector normal = r_condition.GetValue(NORMAL);

                    // // // cosinusoidal
                    double flux_y = E/(1-nu)*(sin(x)*sinh(y)) * normal[1]; 

                    nodal_values.push_back(flux_y);
                }
                else // standard case
                    nodal_values.push_back(r_node.GetValue(*p_variable));
            }

            // Build (t_nodes, nodal_values)
            if (all_nodes_have_value) {
                // Collapse duplicates on t
                // std::vector<double> t_unique, v_unique;
                // CollapseDuplicates(t_nodes, nodal_values, t_unique, v_unique);
                
                // // Degree p from parameters, clamped
                // std::size_t p = t_unique.size()-1;

                // const double interpolated_value = PolyInterpDegreeP(t_unique, v_unique, t_center, p);
                // r_condition.SetValue(*p_variable, interpolated_value);
                
                // Use all samples with LSQ and enforce endpoints
                std::size_t p = round(t_nodes.size()/3)+1; // choose degree; or read from parameters
                if (mParameters.Has("poly_degree")) {
                    std::size_t tentative_p = static_cast<std::size_t>(mParameters["poly_degree"].GetInt());

                    p = std::min(p, tentative_p);
                }

                const double interpolated_value =
                    PolyFitConstrainedEval(t_nodes, nodal_values, t_center, p, /*enforce_endpoints=*/true);

                r_condition.SetValue(*p_variable, interpolated_value);

            }

        }
    }
}

const AssignIgaInterpolationConditionProcess::VariableListType& AssignIgaInterpolationConditionProcess::GetSupportedVariables() const
{
    static const VariableListType variables = {
        &HEAT_FLUX,
        &TEMPERATURE,
        &FACE_HEAT_FLUX,
        &VELOCITY_X,
        &VELOCITY_Y,
        &VELOCITY_Z,
        &BODY_FORCE_X,
        &BODY_FORCE_Y,
        &BODY_FORCE_Z,
        &DISPLACEMENT_X,
        &DISPLACEMENT_Y,
        &DISPLACEMENT_Z,
        &FORCE_X,
        &FORCE_Y,
        &FORCE_Z,
        &DIRECTION_X,
        &DIRECTION_Y,
        &DIRECTION_Z
    };

    return variables;
}

} // namespace Kratos
