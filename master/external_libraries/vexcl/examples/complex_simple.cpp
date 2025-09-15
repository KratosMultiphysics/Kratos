#include <iostream>
#include <complex>
#include <vector>

#include <vexcl/vexcl.hpp>

//---------------------------------------------------------------------------
// Complex multiplication:
//---------------------------------------------------------------------------
VEX_FUNCTION(cl_double2, cmul, (cl_double2, a)(cl_double2, b),
    double2 r = {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
    return r;
);

//---------------------------------------------------------------------------
// Complex division:
//---------------------------------------------------------------------------
VEX_FUNCTION(cl_double2, cdiv, (cl_double2, a)(cl_double2, b),
    double d = b.x * b.x + b.y * b.y;
    double2 r = {(a.x * b.x + a.y * b.y) / d, (a.y * b.x - a.x * b.y) / d};
    return r;
);

int main() {
    vex::Context ctx(vex::Filter::DoublePrecision);
    std::cout << ctx << std::endl;

    const int n = 16;

    // Host side input vectors:
    std::vector<std::complex<double>> x(n), y(n);
    for(int i = 0; i < n; ++i) {
        x[i] = {(double) i,  (double) n-i};
        y[i] = {(double) n-i,(double) i};
    }

    // Transfer host-side complex numbers into device-side cl_double2 vectors:
    vex::vector<cl_double2> X(ctx, n, reinterpret_cast<cl_double2*>(x.data()));
    vex::vector<cl_double2> Y(ctx, n, reinterpret_cast<cl_double2*>(y.data()));

    // Do device-side multiplication and division:
    vex::vector<cl_double2> T(ctx, n);

    T = cmul(X, Y);

    std::vector<std::complex<double>> TT(n);
    vex::copy(T.begin(), T.end(), reinterpret_cast<cl_double2*>(TT.data()));

    for(int i=0; i<n; i++) {
        std::cout << "X * Y = " << x[i] << " * " << y[i] << " = " << TT[i] << "i" << std::endl;
    }

    T = cdiv(X, Y);
    vex::copy(T.begin(), T.end(), reinterpret_cast<cl_double2*>(TT.data()));
    
    for(int i=0; i<n; i++) {
        std::cout << "X / Y = " << x[i] << " / " << y[i] << " = " << TT[i] << std::endl;
    }

}
