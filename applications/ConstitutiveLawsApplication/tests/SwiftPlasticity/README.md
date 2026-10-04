# Small-strain Swift J2 plasticity

`Materials.json` is a minimal example using MPa for stress. Import
`KratosMultiphysics.ConstitutiveLawsApplication` before reading materials and
adjust `model_part_name` and `properties_id` to the analysis. Use a 3D
small-displacement element. Supply density separately if the analysis needs it.
The numbers illustrate the interface; they are not a calibrated steel model.

The law is `SmallStrainSwiftJ2Plasticity3DLaw`. Its inputs are `YOUNG_MODULUS`,
`POISSON_RATIO`, `SWIFT_COEFFICIENT` (K), `SWIFT_INITIAL_STRAIN` (epsilon_0),
and `SWIFT_HARDENING_EXPONENT` (n). The latter two are dimensionless. E, K,
epsilon_0 and n must be finite and positive, and -1 < nu < 0.5. Powers and
elastic moduli must be representable in double precision. There is no input
initial yield stress: it is exactly K*epsilon_0^n (398.1071705534972 MPa here).

## Integration and tangent

History comprises the engineering plastic strain vector and accumulated
equivalent plastic strain p. Stress order is [xx, yy, zz, xy, yz, xz]; strain
order is [xx, yy, zz, 2xy, 2yz, 2xz]. Calculate evaluates from committed
history; Finalize recomputes and commits the accepted strain. GetValue exposes
committed `ACCUMULATED_PLASTIC_STRAIN` and `PLASTIC_STRAIN_VECTOR`.
CalculateValue for `YIELD_STRESS` evaluates the Swift equation at committed p.
Cloning and serialization preserve both history variables.

With r = ||s_trial||, N = s_trial/r and c = sqrt(2/3), the plastic increment is
Delta epsilon_p = gamma*N and Delta p = c*gamma. This matches the multiplier
normalization of the existing SmallStrainJ2Plasticity3D. For a yield function
q-sigma_y, its associated multiplier would instead be Delta p.

The scalar equation is

    R(gamma) = r - 2G*gamma - c*K*(epsilon_0+p_old+c*gamma)^n = 0
    R'(gamma) = -2G - c^2*H(p_old+c*gamma)
    H(p) = n*K*(epsilon_0+p)^(n-1).

A safeguarded Newton iteration maintains a bracket from zero to
(r-c*sigma_y(p_old))/(2G), using bisection when Newton leaves the bracket.
Both yield stress and H are updated at every iterate. The stress residual
tolerance is 1e-12*r, with at most 100 iterations and an explicit error on
failure. A positive trial residual selects the plastic branch; a positive
yield stress keeps r away from zero on that branch.

For fixed committed history, dr = 2G*N:d_epsilon and
dgamma = dr/(2G+c^2*H). Differentiating s = (1-2G*gamma/r)*s_trial gives

    a = 1 - 2G*gamma/r
    b = 2G*(2G/(2G+c^2*H) - 2G*gamma/r)
    C_alg = K_bulk*I (x) I + 2G*a*P_dev - b*N (x) N.

H is evaluated at the converged p. In the engineering Voigt matrix,
2G*P_dev = C_elastic-K_bulk*I (x) I: its shear diagonal is G. The outer
product N_i*N_j has no extra shear factors because N:d_epsilon is the
ordinary dot product of the stress-like N vector with engineering strain.
The elastic branch returns C_elastic. The derivative is one-sided at the
elastic/plastic switch; a single central derivative is not defined there.

## Material-point verification

`tests/cpp_tests/test_small_strain_swift_j2_plasticity_3d.cpp` belongs to
KratosConstitutiveLawsFastSuite; select it with `--gtest_filter=*SwiftJ2*`.
It covers the equation, initial yield, elasticity, transition/pure shear,
uniaxial tension/compression, unloading, strongly nonlinear hardening,
accumulated plastic strain on nonproportional
paths, three tangent regimes, cloning, polymorphic restart, material JSON
registration, invalid parameters, output flags, and strain derived from F.

Uniaxial stress is obtained by Newton-solving the two transverse strains
until the transverse stress norm is below 1e-8 MPa. The independent reference
bisects E*(epsilon_x-p) = K*(epsilon_0+p)^n. Exponents 0.2, 1 and 2 exercise
concave, linear and convex hardening. Equivalent plastic strain is checked
from accumulated sqrt(2/3*Delta epsilon_p:Delta epsilon_p), with half weights
for engineering shear components.

Central differences use strain perturbations of 1e-8 without finalization,
always from identical committed history. The just-plastic point is 0.1%
beyond first yield so neither perturbation crosses the branch switch. The
absolute stiffness tolerance is 3e-7*E, allowing integration-residual/h and
floating-point error; stress and uniaxial reference checks use tighter
absolute tolerances. Tangent checks also assert symmetry and unchanged
committed history.

Only infinitesimal 3D plasticity is supported. The F fallback uses sym(F-I),
not finite-strain kinematics. There is no damage, fracture, rate dependence,
temperature dependence, kinematic hardening, or plane-stress reduction.
