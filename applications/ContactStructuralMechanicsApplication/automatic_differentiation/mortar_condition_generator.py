"""Shared machinery of the mortar contact condition generators.

The five generated condition families (ALM frictionless, ALM frictionless-components,
ALM frictional, penalty frictionless and penalty frictional) share everything except
the functional of each active-set branch. This module owns the shared part:

* the symbolic unknowns, test functions, mortar operators and derived kinematic
  quantities (``SymbolSet``), including the "AD exceptions" of the thesis (Appendix C):
  ``D``, ``M`` and, with normal variation, ``n`` are undefined functions of the DoFs whose
  derivatives are supplied at run time by ``DerivativesUtilities``;
* the differentiation of the functional (RHS = derivative w.r.t. the test functions,
  LHS = minus the derivative of the RHS w.r.t. the DoFs), the replacement of the
  DoF-dependent nodes by plain symbols, the common-subexpression elimination and the
  C++ printing (through ``custom_sympy_fe_utilities``);
* the C++ layout of every family (``FamilySpec``): preambles of ``CalculateLocalLHS`` and of
  the static ``StaticCalculateLocalRHS``, the per-node active-set branching and the
  substitution into the ``*_template.cpp`` file.

The RHS does not depend on the derivatives of the normal, so it is generated only once per
geometry (``TNormalVariation = false``); the ``true`` specialisation forwards to it. The LHS
is generated for both. Everything is assembled in memory and written once.

The functionals themselves (the physics) live in the family scripts / notebooks next to
this module: ``functional(symbols, node, branch) -> sympy expression``.
"""

# System imports
import importlib.util
import os
import re
import time
from dataclasses import dataclass

# External imports
import sympy


def _ImportCustomSympyUtilities():
    """Import ``custom_sympy_fe_utilities`` from the source tree next to this module
    (preferred, so that the generators always use the version they are committed with),
    falling back to the installed ``KratosMultiphysics`` package."""
    source_tree_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python_scripts", "custom_sympy_fe_utilities.py")
    source_tree_file = os.path.abspath(source_tree_file)
    if os.path.isfile(source_tree_file):
        spec = importlib.util.spec_from_file_location("custom_sympy_fe_utilities", source_tree_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    from KratosMultiphysics.ContactStructuralMechanicsApplication import custom_sympy_fe_utilities
    return custom_sympy_fe_utilities


csu = _ImportCustomSympyUtilities()

# (dim, number of slave nodes, number of master nodes) of every registered specialisation
DEFAULT_COMBINATIONS = ((2, 2, 2), (3, 3, 3), (3, 4, 4), (3, 3, 4), (3, 4, 3))

# Active-set branch layouts. Each entry is (branch identifier passed to the functional,
# C++ text opening the branch, extra indentation of the block). The layout also carries the
# text closing the last branch.
BRANCH_LAYOUTS = {
    "active_inactive": {
        "branches": (
            ("inactive", "    if (r_geometry[{node}].IsNot(ACTIVE)) {{ // INACTIVE\n", "    "),
            ("active", "    }} else {{ // ACTIVE\n", "    "),
        ),
        "closing": "    }\n",
    },
    "slip_stick": {
        "branches": (
            ("inactive", "    if (r_geometry[{node}].IsNot(ACTIVE)) {{ // INACTIVE\n", "    "),
            ("slip", "    }} else if (r_geometry[{node}].Is(SLIP)) {{ // ACTIVE-SLIP\n", "    "),
            ("stick", "    }} else {{ // ACTIVE-STICK\n", "    "),
        ),
        "closing": "    }\n",
    },
    "slip_stick_objective": {
        "branches": (
            ("inactive", "    if (r_geometry[{node}].IsNot(ACTIVE)) {{ // INACTIVE\n", "    "),
            ("slip_objective", "    }} else if (r_geometry[{node}].Is(SLIP)) {{ // ACTIVE-SLIP\n        if (is_objetive) {{ // OBJECTIVE-SLIP\n", "        "),
            ("slip_non_objective", "        }} else {{ // NONOBJECTIVE-SLIP\n", "        "),
            ("stick_objective", "        }}\n    }} else {{ // ACTIVE-STICK\n        if (is_objetive) {{ // OBJECTIVE-STICK\n", "        "),
            ("stick_non_objective", "        }} else {{ // NONOBJECTIVE-STICK\n", "        "),
        ),
        "closing": "        }\n    }\n",
    },
}

_BANNER = "/***********************************************************************************/\n/***********************************************************************************/\n"

_LHS_SIGNATURE = (
    "template<>\n"
    "void {class_name}<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::CalculateLocalLHS(\n"
    "    Matrix& rLocalLHS,\n"
    "    const MortarConditionMatrices& rMortarConditionMatrices,\n"
    "    const DerivativeDataType& rDerivativeData,\n"
    "    const IndexType rActiveInactive,\n"
    "    const ProcessInfo& rCurrentProcessInfo\n"
    "    )\n"
    "{{\n"
    "    // Initialize\n"
    "    for (std::size_t i = 0; i < MatrixSize; ++i)\n"
    "        for (std::size_t j = 0; j < MatrixSize; ++j)\n"
    "            rLocalLHS(i, j) = 0.0;\n"
    "\n"
    "    // The geometry of the condition\n"
    "    const GeometryType& r_geometry = this->GetParentGeometry();\n"
    "\n"
)

_RHS_PARAMETERS_FRICTIONAL = (
    "    PairedCondition* pCondition,\n"
    "    const MortarBaseConditionMatrices& rPreviousMortarOperators,\n"
    "    const array_1d<double, TNumNodes>& mu,\n"
    "    Vector& rLocalRHS,\n"
)

_RHS_PARAMETERS_FRICTIONLESS = (
    "    PairedCondition* pCondition,\n"
    "    Vector& rLocalRHS,\n"
)

_RHS_SIGNATURE = (
    "template<>\n"
    "void {class_name}<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::StaticCalculateLocalRHS(\n"
    "{parameters}"
    "    const MortarConditionMatrices& rMortarConditionMatrices,\n"
    "    const DerivativeDataType& rDerivativeData,\n"
    "    const IndexType rActiveInactive,\n"
    "    const ProcessInfo& rCurrentProcessInfo\n"
    "    )\n"
    "{{\n"
)

_RHS_BODY_START = (
    "    // Initialize\n"
    "    for (std::size_t i = 0; i < MatrixSize; ++i)\n"
    "        rLocalRHS[i] = 0.0;\n"
    "\n"
    "    // The geometry of the condition\n"
    "    const GeometryType& r_geometry = pCondition->GetParentGeometry();\n"
    "\n"
)

_RHS_FORWARDER_FRICTIONAL = (
    "    {class_name}<TDim,TNumNodes, false, TNumNodesMaster>::StaticCalculateLocalRHS(\n"
    "      pCondition,\n"
    "      rPreviousMortarOperators,\n"
    "      mu,\n"
    "      rLocalRHS,\n"
    "      rMortarConditionMatrices,\n"
    "      rDerivativeData,\n"
    "      rActiveInactive,\n"
    "      rCurrentProcessInfo\n"
    "      );\n"
    "}}\n"
)

_RHS_FORWARDER_FRICTIONLESS = (
    "    {class_name}<TDim,TNumNodes, false, TNumNodesMaster>::StaticCalculateLocalRHS(\n"
    "      pCondition,\n"
    "      rLocalRHS,\n"
    "      rMortarConditionMatrices,\n"
    "      rDerivativeData,\n"
    "      rActiveInactive,\n"
    "      rCurrentProcessInfo\n"
    "      );\n"
    "}}\n"
)

_DELTA_NORMAL_LINE = "    const array_1d<BoundedMatrix<double, TNumNodes, TDim>, SIZEDERIVATIVES1>& DeltaNormalSlave = rDerivativeData.DeltaNormalSlave;\n\n"


@dataclass
class FamilySpec:
    """Static description of one generated condition family.

    Attributes:
    name -- Short name (folder of the generator).
    class_name -- C++ class template name.
    template_file -- ``*_template.cpp`` with the ``// replace_lhs`` / ``// replace_rhs`` markers.
    output_file -- Name of the generated ``.cpp`` (written into ``custom_conditions/``).
    lm_kind -- Lagrange multiplier DoFs: ``"scalar"`` (``LAGRANGE_MULTIPLIER_CONTACT_PRESSURE``),
               ``"vector"`` (``VECTOR_LAGRANGE_MULTIPLIER``) or ``"none"`` (penalty).
    frictional -- Whether the family carries the frictional quantities (tangent, previous
                  operators, friction coefficient) and the static RHS receives them.
    branch_layout -- Key of ``BRANCH_LAYOUTS``.
    preamble_values -- C++ declarations shared by the LHS and the RHS bodies, written with the
                       placeholders ``TDim``, ``TNumNodes``, ``TNumNodesMaster`` and
                       ``{previous_operators}`` (``mPreviousMortarOperators`` in the member LHS,
                       ``rPreviousMortarOperators`` in the static RHS).
    preamble_lhs_only -- Declarations that exist only in the LHS (operator derivatives, friction coefficient).
    preamble_common_tail -- Declarations emitted after the LHS-only block in both bodies
                            (the objective/non-objective switch), with ``{condition}`` being
                            ``this`` or ``pCondition``.
    """
    name: str
    class_name: str
    template_file: str
    output_file: str
    lm_kind: str
    frictional: bool
    branch_layout: str
    preamble_values: str
    preamble_lhs_only: str = ""
    preamble_common_tail: str = ""

    def NumberOfDofs(self, dim, nnodes, nnodes_master):
        """Size of the local system."""
        number_dof = dim * (nnodes + nnodes_master)
        if self.lm_kind == "scalar":
            number_dof += nnodes
        elif self.lm_kind == "vector":
            number_dof += dim * nnodes
        elif self.lm_kind != "none":
            raise ValueError("Unknown lm_kind '{}'".format(self.lm_kind))
        return number_dof

    @property
    def branches(self):
        return [branch[0] for branch in BRANCH_LAYOUTS[self.branch_layout]["branches"]]


class SymbolSet:
    """All the symbols and derived quantities of one geometry combination.

    Notation (thesis chapter 4): ``x1 = X1 + u1`` and ``x2 = X2 + u2`` are the current slave and
    master coordinates; ``x1old``/``x2old`` those of the previous step (``u1old`` being the offset from
    ``X1`` to the previous position). ``D``/``M`` are the mortar operators, ``DOperatorold``/``MOperatorold``
    those of the previous step (constants). ``w1``, ``w2``, ``wLM`` are the test functions.

    Sign conventions:
    * the weighted gap ``NormalGap = -n.(D x1 - M x2)`` is positive when the bodies are apart;
    * every ``<quantity>w`` (``NormalwGap``, ``TangentwSlip*``) is **minus** the variation of the
      quantity in the direction of the test functions (``Xw = -dX``), so that the virtual work of
      a traction ``t`` reads ``t . Xw`` and its linearisation is consistent; ``NormalwGap`` is
      ``+n.(D w1 - M w2)`` because ``NormalGap`` carries a minus sign;
    * both tangential slips approximate the same relative tangential motion (they give
      ``-D_jj * delta`` for a slave sliding by ``delta`` over a fixed master, like ``WEIGHTED_SLIP``):
      ``TangentSlipObjective = tau[(D - Dold) x1 - (M - Mold) x2]`` and
      ``TangentSlipNonObjective = -tau[D (x1 - x1old) - M (x2 - x2old)]`` (thesis eqs. 4.65-4.69).
    """

    def __init__(self, dim, nnodes, nnodes_master, normal_variation, lm_kind, frictional):
        self.dim = dim
        self.nnodes = nnodes
        self.nnodes_master = nnodes_master
        self.normal_variation = normal_variation
        self.lm_kind = lm_kind
        self.frictional = frictional

        # Unknowns
        self.u1 = csu.DefineMatrix("u1", nnodes, dim)          # current displacement of the slave nodes
        self.u2 = csu.DefineMatrix("u2", nnodes_master, dim)   # current displacement of the master nodes
        if lm_kind == "vector":
            self.LM = csu.DefineMatrix("LM", nnodes, dim)      # vector Lagrange multiplier
        elif lm_kind == "scalar":
            self.LMNormal = csu.DefineVector("LMNormal", nnodes)  # normal (scalar) Lagrange multiplier
        if frictional:
            self.u1old = csu.DefineMatrix("u1old", nnodes, dim)
            self.u2old = csu.DefineMatrix("u2old", nnodes_master, dim)

        # Test functions
        self.w1 = csu.DefineMatrix("w1", nnodes, dim)
        self.w2 = csu.DefineMatrix("w2", nnodes_master, dim)
        if lm_kind == "vector":
            self.wLM = csu.DefineMatrix("wLM", nnodes, dim)
        elif lm_kind == "scalar":
            self.wLMNormal = csu.DefineVector("wLMNormal", nnodes)

        # Geometry, normal and tangent
        self.X1 = csu.DefineMatrix("X1", nnodes, dim)
        self.X2 = csu.DefineMatrix("X2", nnodes_master, dim)
        self.NormalSlave = csu.DefineMatrix("NormalSlave", nnodes, dim)
        if frictional:
            self.TangentSlave = csu.DefineMatrix("TangentSlave", nnodes, dim)

        # Mortar operators
        self.DOperator = csu.DefineMatrix("DOperator", nnodes, nnodes)
        self.MOperator = csu.DefineMatrix("MOperator", nnodes, nnodes_master)
        if frictional:
            self.DOperatorold = csu.DefineMatrix("DOperatorold", nnodes, nnodes)
            self.MOperatorold = csu.DefineMatrix("MOperatorold", nnodes, nnodes_master)

        # Parameters
        self.DynamicFactor = csu.DefineVector("DynamicFactor", nnodes)
        self.PenaltyParameter = csu.DefineVector("PenaltyParameter", nnodes)
        self.ScaleFactor = sympy.Symbol("ScaleFactor", positive=True)
        if frictional:
            self.mu = csu.DefineVector("mu", nnodes)
            self.TangentFactor = sympy.Symbol("TangentFactor", positive=True)
            self.delta_time = sympy.Symbol("delta_time", positive=True)

        # DoF lists (the ordering of u1_var + u2_var is the one of DeltaDOperator[k] / DeltaMOperator[k])
        self.u1_var = []
        self.u2_var = []
        csu.CreateVariableMatrixList(self.u1_var, self.u1)
        csu.CreateVariableMatrixList(self.u2_var, self.u2)
        self.u12_var = self.u1_var + self.u2_var

        # AD exceptions: D, M (and n with normal variation) depend on the DoFs
        self.dependencies = {}
        if normal_variation:
            self.NormalSlave = csu.DefineDofDependencyMatrix(self.NormalSlave, self.u1_var)
            self.dependencies["NormalSlave"] = (self.NormalSlave, self.u1_var)
        self.DOperator = csu.DefineDofDependencyMatrix(self.DOperator, self.u12_var)
        self.MOperator = csu.DefineDofDependencyMatrix(self.MOperator, self.u12_var)
        self.dependencies["DOperator"] = (self.DOperator, self.u12_var)
        self.dependencies["MOperator"] = (self.MOperator, self.u12_var)
        self.replacement = csu.BuildDependencyReplacement(self.dependencies)

        # Normal and tangential components of the multiplier
        if lm_kind == "vector":
            self.LMNormal = csu.DefineVector("LMNormal", nnodes)
            self.wLMNormal = csu.DefineVector("wLMNormal", nnodes)
            self.LMTangent = csu.DefineMatrix("LMTangent", nnodes, dim)
            self.wLMTangent = csu.DefineMatrix("wLMTangent", nnodes, dim)
            for node in range(nnodes):
                self.LMNormal[node] = self.LM.row(node).dot(self.NormalSlave.row(node))
                self.wLMNormal[node] = self.wLM.row(node).dot(self.NormalSlave.row(node))
                for idim in range(dim):
                    self.LMTangent[node, idim] = self.LM[node, idim] - self.LMNormal[node] * self.NormalSlave[node, idim]
                    self.wLMTangent[node, idim] = self.wLM[node, idim] - self.wLMNormal[node] * self.NormalSlave[node, idim]

        # Kinematics
        self.x1 = self.X1 + self.u1
        self.x2 = self.X2 + self.u2
        self.Dx1Mx2 = self.DOperator * self.x1 - self.MOperator * self.x2
        self.Dw1Mw2 = self.DOperator * self.w1 - self.MOperator * self.w2
        self.NormalGap = csu.DefineVector("NormalGap", nnodes)
        self.NormalwGap = csu.DefineVector("NormalwGap", nnodes)
        for node in range(nnodes):
            self.NormalGap[node] = - self.Dx1Mx2.row(node).dot(self.NormalSlave.row(node))
            self.NormalwGap[node] = self.Dw1Mw2.row(node).dot(self.NormalSlave.row(node))

        if frictional:
            self.x1old = self.X1 + self.u1old
            self.x2old = self.X2 + self.u2old
            self.DDeltax1MDeltax2 = self.DOperator * (self.x1 - self.x1old) - self.MOperator * (self.x2 - self.x2old)
            self.DeltaDx1DeltaMx2 = (self.DOperator - self.DOperatorold) * self.x1 - (self.MOperator - self.MOperatorold) * self.x2
            self.DeltaDw1DeltaMw2 = (self.DOperator - self.DOperatorold) * self.w1 - (self.MOperator - self.MOperatorold) * self.w2
            self.TangentSlipNonObjective = csu.DefineMatrix("TangentSlipNonObjective", nnodes, dim)
            self.TangentwSlipNonObjective = csu.DefineMatrix("TangentwSlipNonObjective", nnodes, dim)
            self.TangentSlipObjective = csu.DefineMatrix("TangentSlipObjective", nnodes, dim)
            self.TangentwSlipObjective = csu.DefineMatrix("TangentwSlipObjective", nnodes, dim)
            for node in range(nnodes):
                normal = self.NormalSlave.row(node)
                # Estimations of the time derivative of the relative motion (the delta_time cancels out)
                gap_time_derivative_non_objective = - self.DDeltax1MDeltax2.row(node) / self.delta_time
                gap_time_derivative_non_objective_w = self.Dw1Mw2.row(node) / self.delta_time
                gap_time_derivative_objective = self.DeltaDx1DeltaMx2.row(node) / self.delta_time
                gap_time_derivative_objective_w = - self.DeltaDw1DeltaMw2.row(node) / self.delta_time
                # Tangential projection (I - n x n)
                aux_slip_non_objective = self.delta_time * (gap_time_derivative_non_objective - gap_time_derivative_non_objective.dot(normal) * normal)
                aux_wslip_non_objective = self.delta_time * (gap_time_derivative_non_objective_w - gap_time_derivative_non_objective_w.dot(normal) * normal)
                aux_slip_objective = self.delta_time * (gap_time_derivative_objective - gap_time_derivative_objective.dot(normal) * normal)
                aux_wslip_objective = self.delta_time * (gap_time_derivative_objective_w - gap_time_derivative_objective_w.dot(normal) * normal)
                for idim in range(dim):
                    self.TangentSlipNonObjective[node, idim] = aux_slip_non_objective[idim]
                    self.TangentwSlipNonObjective[node, idim] = aux_wslip_non_objective[idim]
                    self.TangentSlipObjective[node, idim] = aux_slip_objective[idim]
                    self.TangentwSlipObjective[node, idim] = aux_wslip_objective[idim]

        # DoFs and test functions, ordered as [master displacements, slave displacements, multipliers]
        dofs = []
        testfunc = []
        for i in range(nnodes_master):
            for k in range(dim):
                dofs.append(self.u2[i, k])
                testfunc.append(self.w2[i, k])
        for i in range(nnodes):
            for k in range(dim):
                dofs.append(self.u1[i, k])
                testfunc.append(self.w1[i, k])
        if lm_kind == "vector":
            for i in range(nnodes):
                for k in range(dim):
                    dofs.append(self.LM[i, k])
                    testfunc.append(self.wLM[i, k])
        elif lm_kind == "scalar":
            for i in range(nnodes):
                dofs.append(self.LMNormal[i])
                testfunc.append(self.wLMNormal[i])
        self.dofs = sympy.Matrix(dofs)
        self.testfunc = sympy.Matrix(testfunc)
        self.number_dof = len(dofs)


def _Indent(block, spaces):
    """Indent every non-empty line of ``block`` by ``spaces``."""
    return "".join((spaces + line if line.strip() else line) for line in block.splitlines(True))


def _StripTrailingWhitespace(text):
    return "\n".join(line.rstrip() for line in text.split("\n"))


def _SubstituteTemplateParameters(text, dim, nnodes, nnodes_master, normal_variation, matrix_size):
    """Replace the placeholders of the emitted bodies (never applied to the template file itself)."""
    replacements = (
        ("TNumNodesMaster", str(nnodes_master)),
        ("TNumNodes", str(nnodes)),
        ("TDim", str(dim)),
        ("TNormalVariation", "true" if normal_variation else "false"),
        ("MatrixSize", str(matrix_size)),
        ("SIZEDERIVATIVES1", str(nnodes * dim)),
        ("SIZEDERIVATIVES2", str((nnodes + nnodes_master) * dim)),
    )
    for placeholder, value in replacements:
        text = re.sub(r"\b" + placeholder + r"\b", value, text)
    return text


def _ComputeBranch(symbols, functional, node, branch, do_simplifications, mode, log):
    """Differentiate the functional of one branch and return the C++ text of its LHS and RHS blocks."""
    t0 = time.time()
    rv = sympy.Matrix([[functional(symbols, node, branch)]])
    if do_simplifications:
        rv[0, 0] = sympy.simplify(rv[0, 0])
    rhs, lhs = csu.Compute_RHS_and_LHS(rv, symbols.testfunc, symbols.dofs, False)
    lhs = csu.ReplaceDependenciesBySymbols(lhs, symbols.replacement)
    rhs = csu.ReplaceDependenciesBySymbols(rhs, symbols.replacement)
    lhs_out = csu.OutputMatrix_CollectingFactorsNonZero(lhs, "lhs", mode, 1)
    rhs_out = csu.OutputVector_CollectingFactorsNonZero(rhs, "rhs", mode, 1)
    log("    node {} {:<20s} LHS {}x{} ({} non-zero), RHS {} ({} non-zero) [{:.1f} s]".format(
        node, branch, lhs.shape[0], lhs.shape[1], lhs_out.count("lhs("), rhs.shape[0], rhs_out.count("rhs["), time.time() - t0))
    lhs_out = lhs_out.replace("lhs(", "rLocalLHS(")
    rhs_out = rhs_out.replace("rhs[", "rLocalRHS[")
    return lhs_out, rhs_out


def _AssembleBranches(spec, nnodes, blocks):
    """Wrap the per-node, per-branch blocks into the C++ active-set dispatch."""
    layout = BRANCH_LAYOUTS[spec.branch_layout]
    text = ""
    for node in range(nnodes):
        text += "\n    // NODE {}\n".format(node)
        for branch, opening, indentation in layout["branches"]:
            text += opening.format(node=node)
            text += _Indent(blocks[(node, branch)], indentation)
        text += layout["closing"]
    return text


def GenerateSpecialisation(spec, functional, dim, nnodes, nnodes_master, normal_variation, do_simplifications=False, mode="c", log=print):
    """Generate the ``CalculateLocalLHS`` and ``StaticCalculateLocalRHS`` bodies of one specialisation.

    Returns ``(lhs_text, rhs_text)``. For ``normal_variation=True`` the RHS text is the forwarder
    to the ``false`` specialisation and no RHS is differentiated.
    """
    symbols = SymbolSet(dim, nnodes, nnodes_master, normal_variation, spec.lm_kind, spec.frictional)
    matrix_size = spec.NumberOfDofs(dim, nnodes, nnodes_master)
    if symbols.number_dof != matrix_size:
        raise RuntimeError("Inconsistent number of DoFs: {} vs {}".format(symbols.number_dof, matrix_size))
    log("  {} {}D {}N{}N normal variation {}: {} DoFs".format(spec.class_name, dim, nnodes, nnodes_master, normal_variation, matrix_size))

    lhs_blocks = {}
    rhs_blocks = {}
    for node in range(nnodes):
        for branch in spec.branches:
            lhs_blocks[(node, branch)], rhs_blocks[(node, branch)] = _ComputeBranch(symbols, functional, node, branch, do_simplifications, mode, log)

    # LHS
    previous_operators = "mPreviousMortarOperators"
    lhs_text = _BANNER + "\n" + _LHS_SIGNATURE.format(class_name=spec.class_name)
    lhs_text += spec.preamble_values.replace("{previous_operators}", previous_operators)
    lhs_text += spec.preamble_lhs_only
    lhs_text += spec.preamble_common_tail.replace("{condition}", "this")
    if normal_variation:
        lhs_text += _DELTA_NORMAL_LINE
    lhs_text += _AssembleBranches(spec, nnodes, lhs_blocks)
    lhs_text += "}\n"

    # RHS
    parameters = _RHS_PARAMETERS_FRICTIONAL if spec.frictional else _RHS_PARAMETERS_FRICTIONLESS
    rhs_text = _BANNER + "\n" + _RHS_SIGNATURE.format(class_name=spec.class_name, parameters=parameters)
    if normal_variation:
        rhs_text += (_RHS_FORWARDER_FRICTIONAL if spec.frictional else _RHS_FORWARDER_FRICTIONLESS).format(class_name=spec.class_name)
    else:
        rhs_text += _RHS_BODY_START
        rhs_text += spec.preamble_values.replace("{previous_operators}", "rPreviousMortarOperators")
        rhs_text += spec.preamble_common_tail.replace("{condition}", "pCondition")
        rhs_text += _AssembleBranches(spec, nnodes, rhs_blocks)
        rhs_text += "}\n"

    lhs_text = _SubstituteTemplateParameters(lhs_text, dim, nnodes, nnodes_master, normal_variation, matrix_size)
    rhs_text = _SubstituteTemplateParameters(rhs_text, dim, nnodes, nnodes_master, normal_variation, matrix_size)
    return _StripTrailingWhitespace(lhs_text), _StripTrailingWhitespace(rhs_text)


def _JoinBodies(texts):
    """Concatenate the bodies, separated by a blank line and the banner; the template already
    carries the banner before the marker, so the leading one is dropped."""
    code = "\n".join(text.strip("\n") + "\n" for text in texts)
    if code.startswith(_BANNER):
        code = code[len(_BANNER):].lstrip("\n")
    return code


def CheckGeneratedCode(text):
    """Raise if the generated code still contains symbolic leftovers."""
    for pattern in (r"Derivative\(", r"//subsvar_", r"Not supported", r"(?<![A-Za-z0-9_])[A-Za-z]\w*?_\d+_\d+\b", r"\bTNumNodes\b", r"\bTDim\b", r"\bMatrixSize\b"):
        match = re.search(pattern, text)
        if match:
            line = text[:match.start()].count("\n") + 1
            raise RuntimeError("Generated code contains '{}' at line {}".format(match.group(0), line))


def Generate(spec, functional, template_dir, output_dir, combinations=DEFAULT_COMBINATIONS, normal_variations=(False, True), do_simplifications=False, log=print):
    """Generate the whole ``.cpp`` file of a family.

    Keyword arguments:
    spec -- The ``FamilySpec``.
    functional -- ``functional(symbols, node, branch)`` returning the Galerkin functional of the branch.
    template_dir -- Folder containing ``spec.template_file``.
    output_dir -- Folder where ``spec.output_file`` is written.
    combinations -- ``(dim, nnodes, nnodes_master)`` tuples.
    normal_variations -- Values of ``TNormalVariation`` to generate.
    do_simplifications -- Apply ``sympy.simplify`` to the functionals (slow).
    log -- Progress callback.

    Returns the path of the written file.
    """
    t0 = time.time()
    template_path = os.path.join(template_dir, spec.template_file)
    with open(template_path, "r") as template_file:  # universal newlines: the template may be CRLF
        template = template_file.read()
    for marker in ("// replace_lhs", "// replace_rhs"):
        if template.count(marker) != 1:
            raise RuntimeError("Template {} must contain exactly one '{}' marker".format(template_path, marker))

    lhs_texts = []
    rhs_texts = []
    for normal_variation in normal_variations:
        for dim, nnodes, nnodes_master in combinations:
            lhs_text, rhs_text = GenerateSpecialisation(spec, functional, dim, nnodes, nnodes_master, normal_variation, do_simplifications, "c", log)
            lhs_texts.append(lhs_text)
            rhs_texts.append(rhs_text)

    lhs_code = _JoinBodies(lhs_texts)
    rhs_code = _JoinBodies(rhs_texts)
    CheckGeneratedCode(lhs_code)
    CheckGeneratedCode(rhs_code)

    output = template.replace("// replace_lhs", lhs_code).replace("// replace_rhs", rhs_code)
    output = _StripTrailingWhitespace(output)
    output_path = os.path.join(output_dir, spec.output_file)
    with open(output_path, "w", newline="\n") as output_file:
        output_file.write(output)
    log("Written {} ({} lines) in {:.0f} s".format(output_path, output.count("\n"), time.time() - t0))
    return output_path


def DefaultDirectories(family_folder):
    """``(template_dir, output_dir)`` for a generator living in ``automatic_differentiation/<family>``."""
    family_folder = os.path.abspath(family_folder)
    output_dir = os.path.abspath(os.path.join(family_folder, "..", "..", "custom_conditions"))
    return family_folder, output_dir


##############################################################################
# Family specifications
##############################################################################

_ALM_FRICTIONAL_VALUES = (
    "    // Initialize values\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2old = rDerivativeData.u2old;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n"
    "\n"
    "    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);\n"
    "\n"
    "    // The normal and tangent vectors\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim> TangentSlave = MortarUtilities::ComputeTangentMatrix<TNumNodes,TDim>(r_geometry);\n"
    "\n"
    "    // The ALM parameters\n"
    "    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n"
    "    const double ScaleFactor = rDerivativeData.ScaleFactor;\n"
    "    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n"
    "    const double TangentFactor = rDerivativeData.TangentFactor;\n"
    "\n"
    "    // Mortar operators\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperatorold = {previous_operators}.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperatorold = {previous_operators}.DOperator;\n"
    "\n"
)

_ALM_FRICTIONAL_LHS_ONLY = (
    "    // Mortar operators derivatives\n"
    "    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodesMaster>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n"
    "    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n"
    "\n"
    "    // We get the friction coefficient\n"
    "    const array_1d<double, TNumNodes> mu = GetFrictionCoefficient();\n"
    "\n"
)

_OBJECTIVE_SWITCH = (
    "//    // The delta time\n"
    "//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];\n"
    "\n"
    "    const double OperatorThreshold = rCurrentProcessInfo[OPERATOR_THRESHOLD];\n"
    "    const double norm_delta_M = norm_frobenius(MOperator - MOperatorold);\n"
    "    const double norm_delta_D = norm_frobenius(DOperator - DOperatorold);\n"
    "    const bool is_objetive = (norm_delta_D > OperatorThreshold && norm_delta_M > OperatorThreshold) ? true : false;\n"
    "    {condition}->Set(MODIFIED, !is_objetive);\n"
    "\n"
)

ALM_FRICTIONAL = FamilySpec(
    name="ALM_frictional_mortar_condition",
    class_name="AugmentedLagrangianMethodFrictionalMortarContactCondition",
    template_file="ALM_frictional_mortar_contact_condition_template.cpp",
    output_file="ALM_frictional_mortar_contact_condition.cpp",
    lm_kind="vector",
    frictional=True,
    branch_layout="slip_stick_objective",
    preamble_values=_ALM_FRICTIONAL_VALUES,
    preamble_lhs_only=_ALM_FRICTIONAL_LHS_ONLY,
    preamble_common_tail=_OBJECTIVE_SWITCH,
)

_PENALTY_FRICTIONAL_VALUES = (
    "    // Initialize values\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n"
    "//    const BoundedMatrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n"
    "//    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2old = rDerivativeData.u2old;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n"
    "\n"
    "    // The normal and tangent vectors\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim> TangentSlave = ComputeTangentMatrixSlip(r_geometry);\n"
    "\n"
    "    // The penalty parameters\n"
    "    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n"
    "    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n"
    "    const double TangentFactor = rDerivativeData.TangentFactor;\n"
    "\n"
    "    // Mortar operators\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperatorold = {previous_operators}.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperatorold = {previous_operators}.DOperator;\n"
    "\n"
)

_PENALTY_FRICTIONAL_TAIL = (
    "//    // The delta time\n"
    "//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];\n"
    "\n"
)

PENALTY_FRICTIONAL = FamilySpec(
    name="penalty_frictional_mortar_condition",
    class_name="PenaltyMethodFrictionalMortarContactCondition",
    template_file="penalty_frictional_mortar_contact_condition_template.cpp",
    output_file="penalty_frictional_mortar_contact_condition.cpp",
    lm_kind="none",
    frictional=True,
    branch_layout="slip_stick",
    preamble_values=_PENALTY_FRICTIONAL_VALUES,
    preamble_lhs_only=_ALM_FRICTIONAL_LHS_ONLY,
    preamble_common_tail=_PENALTY_FRICTIONAL_TAIL,
)

_FRICTIONLESS_LHS_ONLY = (
    "    // Mortar operators derivatives\n"
    "    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodesMaster>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n"
    "    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n"
    "\n"
)

_ALM_FRICTIONLESS_VALUES = (
    "    // Initialize values\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n"
    "\n"
    "    const array_1d<double, TNumNodes> LMNormal = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0);\n"
    "\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n"
    "\n"
    "    // The ALM parameters\n"
    "    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n"
    "    const double ScaleFactor = rDerivativeData.ScaleFactor;\n"
    "    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n"
    "\n"
    "    // Mortar operators\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"
    "\n"
)

ALM_FRICTIONLESS = FamilySpec(
    name="ALM_frictionless_mortar_condition",
    class_name="AugmentedLagrangianMethodFrictionlessMortarContactCondition",
    template_file="ALM_frictionless_mortar_contact_condition_template.cpp",
    output_file="ALM_frictionless_mortar_contact_condition.cpp",
    lm_kind="scalar",
    frictional=False,
    branch_layout="active_inactive",
    preamble_values=_ALM_FRICTIONLESS_VALUES,
    preamble_lhs_only=_FRICTIONLESS_LHS_ONLY,
)

_ALM_FRICTIONLESS_COMPONENTS_VALUES = (
    "    // Initialize values\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n"
    "\n"
    "    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(r_geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);\n"
    "\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n"
    "\n"
    "    // The ALM parameters\n"
    "    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n"
    "    const double ScaleFactor = rDerivativeData.ScaleFactor;\n"
    "    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n"
    "\n"
    "    // Mortar operators\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"
    "\n"
)

ALM_FRICTIONLESS_COMPONENTS = FamilySpec(
    name="ALM_frictionless_components_mortar_condition",
    class_name="AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition",
    template_file="ALM_frictionless_components_mortar_contact_condition_template.cpp",
    output_file="ALM_frictionless_components_mortar_contact_condition.cpp",
    lm_kind="vector",
    frictional=False,
    branch_layout="active_inactive",
    preamble_values=_ALM_FRICTIONLESS_COMPONENTS_VALUES,
    preamble_lhs_only=_FRICTIONLESS_LHS_ONLY,
)

_PENALTY_FRICTIONLESS_VALUES = (
    "    // Initialize values\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n"
    "    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n"
    "\n"
    "    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n"
    "\n"
    "    // The Penalty parameters\n"
    "    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n"
    "    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n"
    "\n"
    "    // Mortar operators\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n"
    "    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"
    "\n"
)

PENALTY_FRICTIONLESS = FamilySpec(
    name="penalty_frictionless_mortar_condition",
    class_name="PenaltyMethodFrictionlessMortarContactCondition",
    template_file="penalty_frictionless_mortar_contact_condition_template.cpp",
    output_file="penalty_frictionless_mortar_contact_condition.cpp",
    lm_kind="none",
    frictional=False,
    branch_layout="active_inactive",
    preamble_values=_PENALTY_FRICTIONLESS_VALUES,
    preamble_lhs_only=_FRICTIONLESS_LHS_ONLY,
)

FAMILIES = {spec.name: spec for spec in (ALM_FRICTIONLESS, ALM_FRICTIONLESS_COMPONENTS, ALM_FRICTIONAL, PENALTY_FRICTIONLESS, PENALTY_FRICTIONAL)}


def ParseCommandLine(argv=None):
    """Common command line of the thin generator scripts.

    ``--combinations 2,2,2 3,3,3`` restricts the geometries, ``--normal-variation false|true|both``
    the ``TNormalVariation`` values, ``--output-dir`` overrides ``custom_conditions/``.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Generate the local system of a mortar contact condition family with sympy.")
    parser.add_argument("--combinations", nargs="*", default=None, help="(dim,nnodes,nnodes_master) triplets, e.g. 2,2,2 3,3,4")
    parser.add_argument("--normal-variation", choices=("false", "true", "both"), default="both")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--simplify", action="store_true", help="Apply sympy.simplify to the functionals (slow)")
    args = parser.parse_args(argv)
    combinations = DEFAULT_COMBINATIONS
    if args.combinations:
        combinations = tuple(tuple(int(value) for value in triplet.split(",")) for triplet in args.combinations)
    normal_variations = {"false": (False,), "true": (True,), "both": (False, True)}[args.normal_variation]
    return combinations, normal_variations, args.output_dir, args.simplify


def Main(spec, functional, family_folder, argv=None):
    """Entry point shared by the thin generator scripts."""
    combinations, normal_variations, output_dir, simplify = ParseCommandLine(argv)
    template_dir, default_output_dir = DefaultDirectories(family_folder)
    return Generate(spec, functional, template_dir, output_dir or default_output_dir, combinations, normal_variations, simplify)
