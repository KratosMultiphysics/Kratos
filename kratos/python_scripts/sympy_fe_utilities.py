import re
import sympy

def DefineMatrix(name, m, n):
    """
    This method defines a symbolic matrix.

    Keyword arguments:
    - name -- Name of variables.
    - m -- Number of rows.
    - n -- Number of columns.
    """
    return sympy.Matrix(m, n, lambda i, j: sympy.var("{name}_{i}_{j}".format(name=name, i=i, j=j)))

def DefineSymmetricMatrix(name, m, n):
    """
    This method defines a symbolic symmetric matrix.

    Keyword arguments:
    - name -- Name of variables.
    - m -- Number of rows.
    - n -- Number of columns.
    """
    return sympy.Matrix(m, n, lambda i, j:
        sympy.var("{name}_{i}_{j}".format(name=name, i=min(i,j), j=max(i,j))))

def DefineVector( name, m):
    """
    This method defines a symbolic vector.

    Keyword arguments:
    - name -- Name of variables.
    - m -- Number of components.
    """
    return sympy.Matrix(m, 1, lambda i,_: sympy.var("{name}_{i}".format(name=name, i=i)))

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False):
    """
    This method defines shape functions and derivatives.
    Note that partition of unity is imposed
    the name HAS TO BE --> N and DN

    Keyword arguments:
    - nnodes -- Number of nodes
    - dim -- Dimension of the space
    - impose_partion_of_unity -- Impose the partition of unity
    """
    DN = DefineMatrix('DN', nnodes, dim)
    N = DefineVector('N', nnodes)

    #impose partition of unity
    if impose_partion_of_unity:
        N[nnodes-1] = 1
        for i in range(nnodes-1):
            N[nnodes-1] -= N[i]

        DN[nnodes-1,:] = -DN[0,:]
        for i in range(1,nnodes-1):
            DN[nnodes-1,:] -= DN[i,:]

    return N, DN

def StrainToVoigt(M):
    """
    This method transform the strains matrix to Voigt notation.

    Keyword arguments:
    - M -- The strain matrix
    """
    if M.shape[0] == 2:
        vm = sympy.Matrix(3, 1, lambda _: 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = 2.0*M[0,1]
    elif M.shape[0] == 3:
        raise NotImplementedError()
    return vm

def MatrixB(DN):
    """
    This method defines the deformation matrix B.

    Keyword arguments:
    - DN -- The shape function derivatives
    """
    dim = DN.shape[1]
    if dim == 2:
        strain_size = 3
        nnodes = DN.shape[0]
        B = sympy.Matrix(sympy.zeros(strain_size, nnodes*dim))
        for i in range(nnodes):
            for _ in range(dim):
                B[0, i*dim] = DN[i,0]
                B[0, i*dim+1] = 0
                B[1, i*dim] = 0
                B[1, i*dim+1] = DN[i,1]
                B[2, i*dim] = DN[i,1]
                B[2, i*dim+1] = DN[i,0]
    elif dim == 3:
        strain_size = 6
        nnodes = DN.shape[0]
        B = sympy.Matrix(sympy.zeros(strain_size, nnodes*dim))
        for i in range(nnodes):
            B[0, i*3 ] = DN[i, 0]
            B[1, i*3 + 1] = DN[i, 1]
            B[2, i*3 + 2] = DN[i, 2]
            B[3, i*3] = DN[i, 1]
            B[3, i*3 + 1] = DN[i, 0]
            B[4, i*3 + 1] = DN[i, 2]
            B[4, i*3 + 2] = DN[i, 1]
            B[5, i*3] = DN[i, 2]
            B[5, i*3 + 2] = DN[i, 0]
    else:
        print("dimension asked in Matrix B is ", dim)
        raise ValueError("wrong dimension")
    return B

def grad_sym_voigtform(DN, x):
    """
    This method defines a symmetric gradient.

    Keyword arguments:
    - DN -- The shape function derivatives
    - x -- The variable to compute the gradient
    """
    [nnodes, dim] = x.shape

    B = MatrixB(DN)

    # Put the x components one after the other in a vector
    xvec = sympy.Matrix(sympy.zeros(B.shape[1], 1))
    for i in range(nnodes):
        for k in range(dim):
            xvec[i*dim+k] = x[i,k]

    return sympy.simplify(B*xvec)

def DfjDxi(DN,f):
    """
    This method defines a gradient. Returns a matrix D such that D(i,j) = D(fj)/D(xi)

    This is the standard in fluid dynamics, that is:
        D(f1)/D(x1) D(f2)/D(x1) D(f3)/D(x1)
        D(f1)/D(x2) D(f2)/D(x2) D(f3)/D(x2)
        D(f1)/D(x3) D(f2)/D(x3) D(f3)/D(x3)

    Keyword arguments:
    - DN -- The shape function derivatives
    - f-- The variable to compute the gradient
    """
    return sympy.simplify(DN.transpose()*f)

def DfiDxj(DN,f):
    """
    This method defines a gradient This returns a matrix D such that D(i,j) = D(fi)/D(xj).

    This is the standard in structural mechanics, that is:
        D(f1)/D(x1) D(f1)/D(x2) D(f1)/D(x3)
        D(f2)/D(x1) D(f2)/D(x2) D(f2)/D(x3)
        D(f3)/D(x1) D(f3)/D(x2) D(f3)/D(x3)

    Keyword arguments:
    - DN -- The shape function derivatives
    - f -- The variable to compute the gradient
    """
    return (DfjDxi(DN,f)).transpose()

def div(DN,x):
    """
    This method defines the divergence.

    Keyword arguments:
    - DN -- The shape function derivatives
    - x -- The variable to compute the gradient
    """
    if DN.shape != x.shape:
        raise ValueError("shapes are not compatible")

    div_x = 0
    for i in range(DN.shape[0]):
        for k in range(DN.shape[1]):
            div_x += DN[i,k]*x[i,k]

    return sympy.Matrix([sympy.simplify(div_x)])

def SubstituteMatrixValue(where_to_substitute, what_to_substitute, substituted_value):
    """
    This method substitutes values into a matrix.

    Keyword arguments:
    - where_to_substitute -- Coordinates where to substitute
    - what_to_substitute -- Components to substitute
    - substituted_value -- Variable to substitute
    """
    for lll in range(where_to_substitute.shape[0]):
        for kkk in range(where_to_substitute.shape[1]):
            tmp = where_to_substitute[lll, kkk]
            for i in range(what_to_substitute.shape[0]):
                for j in range(what_to_substitute.shape[1]):
                    tmp = tmp.subs(what_to_substitute[i,j], substituted_value[i,j])

            where_to_substitute[lll, kkk] = tmp

    return where_to_substitute

def SubstituteScalarValue(where_to_substitute, what_to_substitute, substituted_value):
    """
    This method substitutes values into a scalar.

    Keyword arguments:
    - where_to_substitute -- Coordinates where to substitute
    - what_to_substitute -- Components to substitute
    - substituted_value -- Variable to substitute
    """
    for lll in range(where_to_substitute.shape[0]):
        tmp  = where_to_substitute[lll]
        tmp = tmp.subs( what_to_substitute, substituted_value)
        where_to_substitute[lll] = tmp
    return where_to_substitute

def Compute_RHS(functional, testfunc, do_simplifications=False):
    """
    This computes the RHS vector.

    Keyword arguments:
    - functional -- The functional to derivate
    - testfunc -- The test functions
    - do_simplifications -- If apply simplifications
    """
    rhs = sympy.Matrix(sympy.zeros(testfunc.shape[0],1))
    for i in range(testfunc.shape[0]):
        rhs[i] = sympy.diff(functional[0,0], testfunc[i])

        if do_simplifications:
            rhs[i] = sympy.simplify(rhs[i])

    return rhs

def Compute_LHS(rhs, testfunc, dofs, do_simplifications=False):
    """
    This computes the LHS matrix.

    Keyword arguments:
    - rhs -- The RHS vector
    - testfunc -- The test functions
    - dofs -- The dofs vectors
    - do_simplifications -- If apply simplifications
    """
    lhs = sympy.Matrix(sympy.zeros(testfunc.shape[0],dofs.shape[0]) )
    for i in range(lhs.shape[0]):
        for j in range(lhs.shape[1]):
            lhs[i,j] = -sympy.diff(rhs[i,0], dofs[j,0])

            if do_simplifications:
                lhs[i,j] = sympy.simplify(lhs[i,j])

    return lhs

def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications=False):
    """
    This computes the LHS matrix and the RHS vector.

    Keyword arguments:
    - functional -- The functional to derivate
    - testfunc -- The test functions
    - dofs -- The dofs vectors
    - do_simplifications -- If apply simplifications
    """
    rhs = Compute_RHS(functional, testfunc, do_simplifications)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    return rhs,lhs

# Output functions
def _Indentation(indentation_level):
    """Returns the indentation string."""
    return "    " * indentation_level

def _CodeGen(language, value):
    return  {
        "c"     : sympy.ccode,
        "python": sympy.pycode
    }[language](value)

def _VariableDeclaration(language, variable_name, variable_expression):
    """"
    Returns the variable declaration, without indentation nor suffix.

    The expression must have been turned into code already.
    """
    return  {
        "c"     : "const double {name} = {expr}",
        "python": "{name} = {expr}"
    }[language].format(name=variable_name, expr=variable_expression)

def _Suffix(language):
    """Returns the endline suffix."""
    return  {
        "c"     : ";\n",
        "python": "\n"
    }[language]

def _ReplaceIndices(language, expression):
    """Replaces array access with underscored variable.

    For matrices: `variable_3_7` becomes `variable[3,7]`
    For vectors:  `variable_3`   becomes `variable[3]`

    The accessor for matrices is chosen according to the language (`[]` vs. `()`)
    """
    #Matrices
    pattern = r"_(\d+)_(\d+)"
    replacement = r"[\1,\2]" if language == 'python' else r"(\1,\2)"
    expression = re.sub(pattern, replacement, expression)

    # Vectors
    pattern = r"_(\d+)"
    replacement = r"[\1]"
    expression = re.sub(pattern, replacement, expression)

    return expression

def OutputScalar(scalar_expression, name, language, indentation_level=0, replace_indices=True, assignment_op="="):
    """
    This function generates code to assign to a (pre-declared) scalar

    Keyword arguments:
    - scalar_expression -- A scalar
    - name -- The name of the variables
    - language -- The language of output
    - indentation_level -- The number of tabulations considered
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    prefix = _Indentation(indentation_level)
    suffix = _Suffix(language)

    fmt = prefix + "{var}{op}{expr}" + suffix

    expression = _CodeGen(language, scalar_expression)
    outstring = fmt.format(var=name, op=assignment_op, expr=expression)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def OutputVector(vector_expression, name, language, indentation_level=0, replace_indices=True, assignment_op="="):
    """
    This function generates code to fill a (pre-declared) vector.

    Keyword arguments:
    - rhs -- The RHS vector
    - name -- The name of the variables
    - language -- The language of output
    - indentation_level -- The number of tabulations considered
    - replace_indices -- Set to `True` to replace matrix_i_j with matrix[i,j] (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    prefix = _Indentation(indentation_level)
    suffix = _Suffix(language)
    fmt = prefix + "{var}[{i}]{op}{expr}" + suffix

    outstring = str("")
    for i in range(vector_expression.shape[0]):
        expression = _CodeGen(language, vector_expression[i,0])
        outstring += fmt.format(var=name, i=i, op=assignment_op, expr=expression)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring


def OutputMatrix(matrix_expression, name, language, indentation_level=0, replace_indices=True, assignment_op="="):
    """
    This function generates code to fill a (pre-declared) matrix.

    Keyword arguments:
    - matrix_expression -- The matrix
    - name -- The name of the variables
    - language -- The language of output
    - indentation_level -- The number of tabulations considered
    - replace_indices -- Set to `True` to replace `matrix_i_j` with `matrix[i,j]` (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    prefix = _Indentation(indentation_level)
    suffix = _Suffix(language)

    fmt = prefix \
          + ("{var}[{i},{j}]{op}{expr}" if language == "python" else "{var}({i},{j}){op}{expr}") \
          + suffix

    outstring = str("")
    for i in range(matrix_expression.shape[0]):
        for j in range(matrix_expression.shape[1]):
            expression = _CodeGen(language, matrix_expression[i,j])
            outstring += fmt.format(var=name, i=i, j=j, op=assignment_op, expr=expression)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def OutputSymbolicVariable(expression, language, replace_indices=True):
    """
    This function generates code from an expression..

    Keyword arguments:
    - expression -- The expression to geneate code from
    - language -- The language of output
    - indentation_level -- The number of tabulations considered
    - max_index -- The maximum index
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    """
    outstring = _CodeGen(language, expression) + _Suffix(language)

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def OutputSymbolicVariableDeclaration(expression, name, language, indentation_level=0, replace_indices=True):
    """
    This function generates code to declare and assign an expression, such as:

    ```C++
        const double variable = expression;

    ```

    Keyword arguments:
    - expression -- The variable to define symbolic
    - language -- The language of output
    - name -- The name of the variables
    - indentation_level -- The number of tabulations considered
    - max_index -- DEPRECATED The maximum index
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    """
    prefix = _Indentation(indentation_level)
    value = _CodeGen(language, expression)
    expr = _VariableDeclaration(language, name, value)
    suffix = _Suffix(language)

    outstring = prefix + expr + suffix

    if replace_indices:
        outstring = _ReplaceIndices(language, outstring)

    return outstring

def _AuxiliaryOutputCollectionFactors(A, name, language, indentation_level, optimizations, replace_indices, assignment_op, output_func):
    """
    This method collects the constants of the replacement for matrices, vectors and scalars.

    Keyword arguments:
    - A -- The  factors
    - name -- The name of the constant
    - language -- The language of replacement
    - indentation_level -- The depth of the indentation (4 spaces per level)
    - optimizations -- The level of optimizations
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    - assignment_op -- The assignment operation
    - output_func -- The output function. Must have the same signature as OutputMatrix and OutputVector
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A, sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] #overwrite lhs with the one with the collected components

    Acoefficient_str = str("")
    for factor in A_factors:
        varname = str(factor[0])
        value = factor[1]
        Acoefficient_str += OutputSymbolicVariableDeclaration(value, varname, language, indentation_level, replace_indices)

    A_out = Acoefficient_str + output_func(A, name, language, indentation_level, replace_indices, assignment_op)
    return A_out


def OutputMatrix_CollectingFactors(A, name, language, indentation_level=0, max_index=None, optimizations='basic', replace_indices=True, assignment_op="="):
    """
    This method collects the constants of the replacement for matrices.

    Keyword arguments:
    - A -- The  factors
    - name -- The name of the constant
    - language -- The language of replacement
    - indentation_level -- The depth of the indentation (4 spaces per level)
    - max_index -- DEPRECATED The max number of indexes
    - optimizations -- The level of optimizations
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputMatrix_CollectingFactors")

    return _AuxiliaryOutputCollectionFactors(A, name, language, indentation_level, optimizations, replace_indices, assignment_op, OutputMatrix)


def OutputVector_CollectingFactors(A, name, language, indentation_level=0, max_index=None, optimizations='basic', replace_indices=True, assignment_op="="):
    """
    This method collects the constants of the replacement for vectors.

    Keyword arguments:
    - A -- The  factors
    - name -- The name of the constant
    - language -- The language of replacement
    - indentation_level -- The depth of the indentation (4 spaces per level)
    - max_index -- DEPRECATED The max number of indexes
    - optimizations -- The level of optimizations
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    if max_index is not None:
        print("Warning: max_index parameter is deprecated in OutputVector_CollectingFactors")

    return _AuxiliaryOutputCollectionFactors(A, name, language, indentation_level, optimizations, replace_indices, assignment_op, OutputVector)


def OutputScalar_CollectingFactors(A, name, language, indentation_level=0, optimizations='basic', replace_indices=True, assignment_op="="):
    """
    This method collects the constants of the replacement for vectors.

    Keyword arguments:
    - A -- The  factors
    - name -- The name of the constant
    - language -- The language of replacement
    - indentation_level -- The depth of the indentation (4 spaces per level)
    - optimizations -- The level of optimizations
    - replace_indices -- Set to `True` to replace matrix[i,j] with matrix_i_j (And similarly for vectors)
    - assignment_op -- The assignment operation
    """
    return _AuxiliaryOutputCollectionFactors(A, name, language, indentation_level, optimizations, replace_indices, assignment_op, OutputScalar)
