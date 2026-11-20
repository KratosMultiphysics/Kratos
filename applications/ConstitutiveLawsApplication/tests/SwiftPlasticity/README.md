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

The public law is a thin wrapper around `GenericSmallStrainIsotropicPlasticity`
with `GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<
VonMisesPlasticPotential<6>>, SwiftHardeningLaw>`. There is no separate Swift
J2 integration implementation. The stateless `SwiftHardeningLaw` policy supplies
only the threshold and its derivative. Another monotone accumulated-strain
hardening policy, such as Voce, could reuse this integration and history without
copying the J2 algorithm; Voce is not implemented here.

The generic law stores accumulated equivalent plastic strain alongside plastic
strain, threshold and plastic dissipation. For each accepted increment,

    Delta p = sqrt(2/3 * Delta epsilon_p : Delta epsilon_p)
    p_new = p_old + Delta p.

This is path length, not the norm of total plastic strain or its stress
projection. Stress order is [xx, yy, zz, xy, yz, xz]; strain order is
[xx, yy, zz, 2xy, 2yz, 2xz]. The tensor conversion used for the norm halves
engineering shear. Calculate uses local copies of committed history; Finalize
recomputes and commits the accepted increment. GetValue exposes committed
`ACCUMULATED_PLASTIC_STRAIN` and `PLASTIC_STRAIN_VECTOR`. CalculateValue for
`YIELD_STRESS` evaluates the Swift equation at committed p. Clone, copy and
serialization preserve the history.

Swift deliberately does not support `INTERNAL_VARIABLES`: the legacy generic
vector omits accumulated p. Its getter and setter reject that variable. For
explicit state transfer, read and set both `PLASTIC_STRAIN_VECTOR` and
`ACCUMULATED_PLASTIC_STRAIN`; cloning and serialization also preserve both.
Existing generic laws retain their legacy `INTERNAL_VARIABLES` layout.

The generic von Mises potential returns the engineering flow vector

    g = 3/(2*q) * [s_xx, s_yy, s_zz, 2*s_xy, 2*s_yz, 2*s_xz]
    Delta epsilon_p = lambda*g.

In tensor notation g:g = 3/2, hence dp/dlambda = sqrt(2/3*g:g) = 1.
The scalar consistency equation in the generic integrator is therefore

    R(lambda) = q_trial - 3G*lambda - sigma_y(p_old+lambda) = 0
    R'(lambda) = -3G - H(p_old+lambda)
    H(p) = n*K*(epsilon_0+p)^(n-1).

The generic plastic denominator is 1/(f:C:g + H*dp/dlambda), with f=g.
A safeguarded Newton iteration brackets the root between zero and
(q_trial-sigma_y(p_old))/(3G); it bisects if Newton leaves the bracket.
Threshold and derivative are evaluated at every current iterate. Newton
corrections are never added as physical path increments. The stress tolerance
is 1e-12*q_trial and the default iteration limit is 100, overridable with
`MAX_NUMBER_NL_CL_ITERATIONS`. Non-convergence throws an informative error.

For N = s_trial/||s_trial|| and a = 3G*lambda/q_trial, the exact derivative
of the discrete return mapping is

    C_alg = K_bulk*I (x) I + 2G*(1-a)*P_dev
            - 2G*(3G/(3G+H)-a)*N (x) N.

H is evaluated at converged p. In engineering Voigt notation P_dev has shear
diagonal 1/2. The N outer product needs no extra shear factors because
N:d_epsilon is the ordinary dot product with engineering strain. This is the
analytical algorithmic tangent, including the derivative of the radial
direction; the continuum rank-one elastoplastic tangent alone is insufficient.
The elastic branch returns C_elastic. A central derivative is not defined at
the elastic/plastic switch itself.

Existing generic laws using the default (void) policy retain their dissipation
curves, projected `EQUIVALENT_PLASTIC_STRAIN`, local iteration and tangent
selection. Their 3D small-strain response also records accumulated p, without
feeding it back into their constitutive equations. Swift needs no fracture
energy or dissipation normalization: `PLASTIC_DISSIPATION` remains zero for
this policy. Accumulated history output is not advertised for reduced 2D or
finite-strain laws, which require separate kinematic derivations.

The public law name, Python constructor, properties and material JSON remain
unchanged. The serialized layout now includes the generic base state; restart
files written by the former dedicated Swift law or the previous supported 3D
small-strain generic layout must be regenerated. Unsupported reduced 2D and
finite-strain laws retain their previous layouts: they do not serialize p.
New-format Swift restart continuation and previous-layout compatibility for
representative plane-strain and finite-strain laws are tested.

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

The `*GenericPlasticityAccumulatedHistory*` tests additionally exercise the
shared history with both Swift and the existing perfect-plasticity J2 law:
proportional loading, unloading/reloading, repeated trials, nonproportional
path length versus total-strain norm/projection, cloning and restart.

Only infinitesimal 3D plasticity is supported. The F fallback uses sym(F-I),
not finite-strain kinematics. There is no damage, fracture, rate dependence,
temperature dependence, kinematic hardening, or plane-stress reduction.
