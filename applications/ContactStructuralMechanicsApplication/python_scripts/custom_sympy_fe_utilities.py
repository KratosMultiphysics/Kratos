"""Symbolic helpers for the code generators of ``automatic_differentiation/``.

This module extends the core ``KratosMultiphysics.sympy_fe_utilities`` with the
pieces that the mortar contact conditions need and that the core does not offer:

* **DoF-dependency injection** (the "AD exceptions" of the thesis, Appendix C).
  The mortar operators ``D``, ``M`` and, when the normal variation is considered,
  the slave normal ``n`` are not expressed in terms of the displacements. They are
  declared as *undefined functions* of the DoFs so that ``sympy.diff`` applies the
  chain rule and produces unevaluated ``Derivative`` nodes. Their values and their
  derivatives are computed at run time by ``DerivativesUtilities``.
* **Replacement of those nodes by plain symbols** (``BuildDependencyReplacement``)
  before the common-subexpression elimination, so that the C++ printer only sees
  ordinary symbols: ``DOperator_0_1(u1_0_0, ...)`` becomes ``DOperator_0_1`` and
  ``Derivative(DOperator_0_1(...), u1_0_0)`` becomes ``DeltaDOperator_0_0_1``.
* **C++ output with collected factors** that skips the zero entries and accumulates
  with ``+=`` (``Output*_CollectingFactorsNonZero``), plus the dense ``=`` variants
  kept for the legacy mesh-tying generator.

The module runs on any modern sympy (tested with 1.14) and does not need a compiled
Kratos: the core utilities are imported from the ``KratosMultiphysics`` package when it
is available and from ``kratos/python_scripts`` of the source tree otherwise.
"""

# System imports
import importlib.util
import os
import re

# External imports
import sympy
from sympy.core.function import AppliedUndef, Derivative


def _ImportCoreSympyUtilities():
    """Import the core ``sympy_fe_utilities`` from the package or from the source tree."""
    try:
        import KratosMultiphysics.sympy_fe_utilities as core_utilities
        return core_utilities
    except ImportError:
        pass

    this_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(this_dir, "..", "..", "..", "kratos", "python_scripts", "sympy_fe_utilities.py"),  # source tree
        os.path.join(this_dir, "..", "sympy_fe_utilities.py"),  # installed package layout
    ]
    for candidate in candidates:
        candidate = os.path.abspath(candidate)
        if os.path.isfile(candidate):
            spec = importlib.util.spec_from_file_location("sympy_fe_utilities", candidate)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise ImportError("Cannot import KratosMultiphysics.sympy_fe_utilities nor locate kratos/python_scripts/sympy_fe_utilities.py")


sympy_fe_utilities = _ImportCoreSympyUtilities()

# Regular expressions turning the underscore-indexed symbol names into C++ accessors.
# Anchored on identifier boundaries so that e.g. ``array_1d`` or ``delta_time`` are untouched.
_IDENTIFIER_START = r"(?<![A-Za-z0-9_])"
_THREE_INDICES = re.compile(_IDENTIFIER_START + r"([A-Za-z]\w*?)_(\d+)_(\d+)_(\d+)\b")
_TWO_INDICES = re.compile(_IDENTIFIER_START + r"([A-Za-z]\w*?)_(\d+)_(\d+)\b")
_ONE_INDEX = re.compile(_IDENTIFIER_START + r"([A-Za-z]\w*?)_(\d+)\b")

##############################################################################
# Symbol definition
##############################################################################

def DefineMatrix(name, m, n, mode="Symbol"):
    """Define a symbolic ``m x n`` matrix whose entries are the symbols ``name_i_j``.

    Keyword arguments:
    name -- Name of the variables.
    m -- Number of rows.
    n -- Number of columns.
    mode -- Kept for backwards compatibility and ignored: entries are always plain
            symbols. Use ``DefineDofDependencyMatrix`` to make them depend on the DoFs.
    """
    return sympy.Matrix(m, n, lambda i, j: sympy.Symbol(name + "_%d_%d" % (i, j)))

def DefineSymmetricMatrix(name, m, n=-1, mode="Symbol"):
    """Define a symbolic symmetric matrix (``name_i_j`` with ``j >= i``).

    Keyword arguments:
    name -- Name of the variables.
    m -- Number of rows.
    n -- Number of columns (defaults to ``m``).
    mode -- Ignored (see ``DefineMatrix``).
    """
    if n < 0:
        n = m
    tmp = DefineMatrix(name, m, n)
    for i in range(tmp.shape[0]):
        for j in range(i + 1, tmp.shape[1]):
            tmp[j, i] = tmp[i, j]
    return tmp

def DefineVector(name, m, mode="Symbol"):
    """Define a symbolic column vector whose entries are the symbols ``name_i``.

    Keyword arguments:
    name -- Name of the variables.
    m -- Number of components.
    mode -- Ignored (see ``DefineMatrix``).
    """
    return sympy.Matrix(m, 1, lambda i, j: sympy.Symbol(name + "_%d" % i))

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False):
    """Define shape functions and derivatives (delegates to the core utilities)."""
    return sympy_fe_utilities.DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

def DefineCustomShapeFunctions(nnodes, dim, name):
    """Define shape functions ``name`` and derivatives ``Dname`` with a custom name."""
    DN = DefineMatrix("D" + name, nnodes, dim)
    N = DefineVector(name, nnodes)
    return N, DN

def DefineJacobian(J, DN, x):
    """Fill the Jacobian ``J`` of the mapping defined by the nodal coordinates ``x``.

    Keyword arguments:
    J -- The Jacobian matrix (modified in place and returned).
    DN -- The shape function derivatives.
    x -- The nodal coordinates.
    """
    nnodes, dim = x.shape
    localdim = dim - 1

    if dim == 2:
        if nnodes == 2:
            J[0, 0] = 0.5 * (x[1, 0] - x[0, 0])
            J[1, 0] = 0.5 * (x[1, 1] - x[0, 1])
    else:
        if nnodes == 3:
            J[0, 0] = - (x[0, 0] + x[1, 0])
            J[1, 0] = - (x[0, 1] + x[1, 1])
            J[2, 0] = - (x[0, 2] + x[1, 2])
            J[0, 1] = - (x[0, 0] + x[2, 0])
            J[1, 1] = - (x[0, 1] + x[2, 1])
            J[2, 1] = - (x[0, 2] + x[2, 2])
        else:
            for i in range(dim):
                for j in range(localdim):
                    J[i, j] = 0
            for i in range(nnodes):
                for k in range(dim):
                    for m in range(localdim):
                        J[k, m] += x[i, k] * DN[i, m]

    return J

##############################################################################
# Variable lists and DoF dependency (the "AD exceptions")
##############################################################################

def CreateVariableMatrixList(variable_list, variable_matrix):
    """Append the entries of ``variable_matrix`` (row by row) to ``variable_list``."""
    nnodes, dim = variable_matrix.shape
    for i in range(nnodes):
        for k in range(dim):
            variable_list.append(variable_matrix[i, k])

def CreateVariableVectorList(variable_list, variable_vector):
    """Append the entries of ``variable_vector`` to ``variable_list``."""
    for i in range(variable_vector.shape[0]):
        variable_list.append(variable_vector[i])

def DefineDofDependencyScalar(scalar, variable_list):
    """Turn the symbol ``scalar`` into an undefined function of ``variable_list``.

    ``DOperator_0_1`` becomes ``DOperator_0_1(u1_0_0, u1_0_1, ...)``, so that
    ``sympy.diff`` produces ``Derivative(DOperator_0_1(...), u1_0_0)`` instead of zero.
    This is the sympy >= 1.3 equivalent of calling a ``Symbol`` (removed in 1.3).
    """
    if isinstance(scalar, AppliedUndef):
        return scalar
    return sympy.Function(str(scalar))(*variable_list)

def DefineDofDependencyVector(vector, variable_list):
    """Inject the dependency on ``variable_list`` into every entry of ``vector``."""
    for i in range(vector.shape[0]):
        vector[i, 0] = DefineDofDependencyScalar(vector[i, 0], variable_list)
    return vector

def DefineDofDependencyMatrix(matrix, variable_list):
    """Inject the dependency on ``variable_list`` into every entry of ``matrix``."""
    for i in range(matrix.shape[0]):
        for k in range(matrix.shape[1]):
            matrix[i, k] = DefineDofDependencyScalar(matrix[i, k], variable_list)
    return matrix

def BuildDependencyReplacement(dependencies, derivative_prefix="Delta"):
    """Build the ``xreplace`` dictionary mapping the DoF-dependent functions and
    their first derivatives to plain symbols.

    ``F_i_j(dofs...)`` is mapped to ``Symbol("F_i_j")`` and
    ``Derivative(F_i_j(dofs...), dof_k)`` to ``Symbol("DeltaF_k_i_j")`` (``k`` being the
    position of the DoF in the dependency list, which must match the ordering used by
    ``DerivativesUtilities`` to fill ``DeltaDOperator[k]``, ``DeltaMOperator[k]`` and
    ``DeltaNormalSlave[k]``). The index rewriting of the output functions then prints
    ``F(i,j)`` and ``DeltaF[k](i,j)``.

    Keyword arguments:
    dependencies -- ``{name: (matrix_or_vector_or_scalar, dof_list)}``
    derivative_prefix -- Prefix of the derivative symbols (``Delta`` by default).
    """
    replacement = {}
    for name, (variable, dof_list) in dependencies.items():
        entries = list(variable) if isinstance(variable, sympy.MatrixBase) else [variable]
        for entry in entries:
            if not isinstance(entry, AppliedUndef):
                raise TypeError("Entry {} of {} is not DoF-dependent; call DefineDofDependency* first".format(entry, name))
            function_name = entry.func.__name__
            if not function_name.startswith(name):
                raise ValueError("Function {} does not belong to the dependency {}".format(function_name, name))
            indices = function_name[len(name):]  # "_i_j", "_i" or ""
            replacement[entry] = sympy.Symbol(function_name)
            for k, dof in enumerate(dof_list):
                replacement[Derivative(entry, dof)] = sympy.Symbol(derivative_prefix + name + "_" + str(k) + indices)
    return replacement

def ReplaceDependenciesBySymbols(expression, replacement):
    """Apply ``replacement`` (see ``BuildDependencyReplacement``) to ``expression`` in a
    single ``xreplace`` pass and check that no DoF-dependent node survives.

    A single pass is required: the ``Derivative`` keys must be matched before the
    functions they contain, otherwise ``Derivative(F_i_j, dof)`` would be left behind.
    """
    result = expression.xreplace(replacement)
    if result.has(Derivative):
        raise RuntimeError("Unreplaced Derivative nodes: {}".format(sorted(str(d) for d in result.atoms(Derivative))[:5]))
    leftover = result.atoms(AppliedUndef)
    if leftover:
        raise RuntimeError("Unreplaced DoF-dependent functions: {}".format(sorted(str(f) for f in leftover)[:5]))
    return result

##############################################################################
# Delegations to the core utilities
##############################################################################

def StrainToVoigt(M):
    """Transform the strain matrix to Voigt notation."""
    return sympy_fe_utilities.StrainToVoigt(M)

def MatrixB(DN):
    """Define the deformation matrix B."""
    return sympy_fe_utilities.MatrixB(DN)

def grad_sym_voigtform(DN, x):
    """Define a symmetric gradient in Voigt form."""
    return sympy_fe_utilities.grad_sym_voigtform(DN, x)

def grad(DN, x):
    """Define a gradient."""
    return sympy_fe_utilities.grad(DN, x)

def DfjDxi(DN, f):
    """Gradient returning ``D(i,j) = D(fj)/D(xi)``."""
    return sympy_fe_utilities.DfjDxi(DN, f)

def DfiDxj(DN, f):
    """Gradient returning ``D(i,j) = D(fi)/D(xj)``."""
    return sympy_fe_utilities.DfiDxj(DN, f)

def div(DN, x):
    """Define the divergence."""
    return sympy_fe_utilities.div(DN, x)

def SubstituteMatrixValue(where_to_substitute, what_to_substitute, substituted_value):
    """Substitute values into a matrix."""
    return sympy_fe_utilities.SubstituteMatrixValue(where_to_substitute, what_to_substitute, substituted_value)

def SubstituteScalarValue(where_to_substitute, what_to_substitute, substituted_value):
    """Substitute values into a scalar."""
    return sympy_fe_utilities.SubstituteScalarValue(where_to_substitute, what_to_substitute, substituted_value)

def Compute_RHS(functional, testfunc, do_simplifications=False):
    """Compute the RHS vector ``r_a = d(functional)/d(w_a)``."""
    return sympy_fe_utilities.Compute_RHS(functional, testfunc, do_simplifications)

def Compute_LHS(rhs, testfunc, dofs, do_simplifications=False):
    """Compute the LHS matrix ``K_ab = -d(r_a)/d(u_b)``."""
    return sympy_fe_utilities.Compute_LHS(rhs, testfunc, dofs, do_simplifications)

def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications=False):
    """Compute the RHS vector and the LHS matrix from the functional."""
    return sympy_fe_utilities.Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications)

##############################################################################
# Output
##############################################################################

def ReplaceIndices(code, mode="c"):
    """Rewrite the underscore-indexed names of ``code`` as accessors.

    ``F_k_i_j`` becomes ``F[k](i,j)`` (C) or ``F[k][i,j]`` (Python),
    ``F_i_j`` becomes ``F(i,j)`` (C) or ``F[i,j]`` (Python) and ``F_i`` becomes ``F[i]``.
    """
    if mode == "c":
        code = _THREE_INDICES.sub(r"\1[\2](\3,\4)", code)
        code = _TWO_INDICES.sub(r"\1(\2,\3)", code)
    else:
        code = _THREE_INDICES.sub(r"\1[\2][\3,\4]", code)
        code = _TWO_INDICES.sub(r"\1[\2,\3]", code)
    code = _ONE_INDEX.sub(r"\1[\2]", code)
    return code

def SubstituteIndex(outstring, mode="python", max_index=30):
    """Backwards-compatible alias of ``ReplaceIndices`` (``max_index`` is ignored)."""
    return ReplaceIndices(outstring, mode)

def _CodeGen(expression, mode):
    """Print ``expression`` in the requested language, with C++ math functions."""
    if mode == "c":
        code = sympy.ccode(expression)  # strict: any leftover Derivative/AppliedUndef raises
        return code.replace("pow(", "std::pow(").replace("sqrt(", "std::sqrt(")
    elif mode == "python":
        return sympy.pycode(expression)
    else:
        raise ValueError("Unknown output mode '{}'".format(mode))

def _Assignment(name, mode, indices, assignment_op, expression, initial_spaces):
    """One output line ``name(i,j) op expression;``."""
    if mode == "c":
        accessor = "(" + ",".join(str(i) for i in indices) + ")" if len(indices) > 1 else "[" + str(indices[0]) + "]"
        return initial_spaces + name + accessor + assignment_op + _CodeGen(expression, mode) + ";\n"
    accessor = "[" + ",".join(str(i) for i in indices) + "]"
    return initial_spaces + name + accessor + assignment_op + _CodeGen(expression, mode) + "\n"

def OutputVector(rhs, name, mode="python", initial_tabs=1, max_index=None, assignment_op="=", skip_zeros=False):
    """Convert into text the vector ``rhs`` (one assignment per entry).

    Keyword arguments:
    rhs -- The vector.
    name -- The name of the variable.
    mode -- The output language ("c" or "python").
    initial_tabs -- The number of tabulations considered.
    max_index -- Ignored (kept for backwards compatibility).
    assignment_op -- The assignment operator ("=" or "+=").
    skip_zeros -- Do not print the entries that are identically zero.
    """
    initial_spaces = "    " * initial_tabs
    outstring = ""
    for i in range(rhs.shape[0]):
        if skip_zeros and rhs[i, 0] == 0:
            continue
        outstring += _Assignment(name, mode, (i,), assignment_op, rhs[i, 0], initial_spaces)
    return ReplaceIndices(outstring, mode)

def OutputMatrix(lhs, name, mode="python", initial_tabs=1, max_index=None, assignment_op="=", skip_zeros=False):
    """Convert into text the matrix ``lhs`` (one assignment per entry).

    Keyword arguments:
    lhs -- The matrix.
    name -- The name of the variable.
    mode -- The output language ("c" or "python").
    initial_tabs -- The number of tabulations considered.
    max_index -- Ignored (kept for backwards compatibility).
    assignment_op -- The assignment operator ("=" or "+=").
    skip_zeros -- Do not print the entries that are identically zero.
    """
    initial_spaces = "    " * initial_tabs
    outstring = ""
    for i in range(lhs.shape[0]):
        for j in range(lhs.shape[1]):
            if skip_zeros and lhs[i, j] == 0:
                continue
            outstring += _Assignment(name, mode, (i, j), assignment_op, lhs[i, j], initial_spaces)
    return ReplaceIndices(outstring, mode)

def OutputVectorNonZero(rhs, name, mode="python", initial_tabs=1, max_index=None):
    """Convert into text the non-zero entries of the vector ``rhs``, accumulating with ``+=``."""
    return OutputVector(rhs, name, mode, initial_tabs, max_index, "+=", True)

def OutputMatrixNonZero(lhs, name, mode="python", initial_tabs=1, max_index=None):
    """Convert into text the non-zero entries of the matrix ``lhs``, accumulating with ``+=``."""
    return OutputMatrix(lhs, name, mode, initial_tabs, max_index, "+=", True)

def OutputSymbolicVariable(var, mode="python", varname="", initial_tabs=1, max_index=None):
    """Convert into text the expression ``var`` (one line, with the line terminator)."""
    initial_spaces = "    " * initial_tabs
    if mode == "c":
        return initial_spaces + _CodeGen(var, mode) + ";\n"
    return initial_spaces + _CodeGen(var, mode) + "\n"

def _CollectFactors(A, name, mode, initial_tabs, optimizations):
    """Common-subexpression elimination of ``A``; returns the reduced ``A`` and the
    declarations of the factors ``c<name><n>``."""
    symbol_name = "c" + name
    factors, collected = sympy.cse(A, sympy.numbered_symbols(symbol_name), optimizations)
    A = collected[0]
    initial_spaces = "    " * initial_tabs
    declaration = "const double " if mode == "c" else ""
    coefficients = ""
    for varname, value in factors:
        # NOTE: the expression is indented by one level after "=", reproducing the historical layout of the generated files
        coefficients += initial_spaces + declaration + str(varname) + " = " + OutputSymbolicVariable(value, mode, str(varname), 1)
    return A, ReplaceIndices(coefficients, mode)

def OutputMatrix_CollectingFactors(A, name, mode, initial_tabs=1, max_index=None, optimizations="basic"):
    """Convert into text the matrix ``A`` after collecting its common factors (dense, ``=``)."""
    A, coefficients = _CollectFactors(A, name, mode, initial_tabs, optimizations)
    return coefficients + "\n" + OutputMatrix(A, name, mode, initial_tabs, max_index)

def OutputVector_CollectingFactors(A, name, mode, initial_tabs=1, max_index=None, optimizations="basic"):
    """Convert into text the vector ``A`` after collecting its common factors (dense, ``=``)."""
    A, coefficients = _CollectFactors(A, name, mode, initial_tabs, optimizations)
    return coefficients + "\n" + OutputVector(A, name, mode, initial_tabs, max_index)

def OutputMatrix_CollectingFactorsNonZero(A, name, mode, initial_tabs=1, max_index=None, optimizations="basic"):
    """Convert into text the non-zero entries of ``A`` after collecting its common factors (``+=``)."""
    A, coefficients = _CollectFactors(A, name, mode, initial_tabs, optimizations)
    return coefficients + "\n" + OutputMatrixNonZero(A, name, mode, initial_tabs, max_index)

def OutputVector_CollectingFactorsNonZero(A, name, mode, initial_tabs=1, max_index=None, optimizations="basic"):
    """Convert into text the non-zero entries of ``A`` after collecting its common factors (``+=``)."""
    A, coefficients = _CollectFactors(A, name, mode, initial_tabs, optimizations)
    return coefficients + "\n" + OutputVectorNonZero(A, name, mode, initial_tabs, max_index)
