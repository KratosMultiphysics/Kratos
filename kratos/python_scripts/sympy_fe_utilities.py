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

def DefineSymmetricFourthOrderTensor(name, m, n, o, p):
    """
    This method defines a symbolic symmetric 4th order tensor.

    Keyword arguments:
    - name -- Name of variables.
    - m -- 1st dimension.
    - n -- 2nd dimension.
    - o -- 3rd dimension.
    - p -- 4th dimension.
    """
    if m != n:
        raise ValueError("Provided sizes do not respect first symmetry.")
    if o != p:
        raise ValueError("Provided sizes do not respect second symmetry.")

    tensor = sympy.MutableDenseNDimArray(sympy.zeros(m**4),shape=(m,n,o,p))
    for i in range(m):
        for j in range(n):
            for k in range(o):
                for l in range(p):
                    tensor[i,j,k,l] = sympy.var("{name}_{m}_{n}_{o}_{p}".format(name=name, m=min(i,j), n=max(i,j), o=min(k,l), p=max(k,l)))

    return tensor

def DefineVector( name, m):
    """
    This method defines a symbolic vector.

    Keyword arguments:
    - name -- Name of variables.
    - m -- Number of components.
    """
    return sympy.Matrix(m, 1, lambda i,_: sympy.var("{name}_{i}".format(name=name, i=i)))

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False, shape_functions_name='N', first_derivatives_name='DN', second_derivatives_name=None):
    """
    This method defines shape functions and derivatives.
    Second order derivatives can be optionally defined as well.
    Note that partition of unity can be imposed.

    Keyword arguments:
    - nnodes -- Number of nodes
    - dim -- Dimension of the space
    - impose_partion_of_unity -- Impose the partition of unity
    - shape_functions_name -- Name for the shape functions symbols
    - first_derivatives_name -- Name for the shape functions first derivatives symbols
    - second_derivatives_name -- Name for the shape functions second derivatives symbols
    """
    DN = DefineMatrix(first_derivatives_name, nnodes, dim)
    N = DefineVector(shape_functions_name, nnodes)

    # Impose partition of unity
    if impose_partion_of_unity:
        N[nnodes-1] = 1
        for i in range(nnodes-1):
            N[nnodes-1] -= N[i]

        DN[nnodes-1,:] = -DN[0,:]
        for i in range(1,nnodes-1):
            DN[nnodes-1,:] -= DN[i,:]

    if second_derivatives_name:
        if not impose_partion_of_unity:
            DDN = sympy.Matrix(1, nnodes, lambda _, j : (sympy.Matrix(dim, dim, lambda m, n : sympy.var(f"{second_derivatives_name}_{j}_{m}_{n}"))))
        else:
            raise Exception("Partition of unity imposition is not implemented for shape functions second derivatives.")

    if second_derivatives_name:
        return N, DN, DDN
    else:
        return N, DN

def StrainToMatrix(V):
    """
    This method transfoms a Voigt strain vector to matrix.

    Keyword arguments:
    - V -- The strain vector in Voigt notation
    """
    if V.shape[1] != 1:
        raise Exception(f"Provided array is not a vector. Shape is {V.shape}.")

    if V.shape[0] == 3:
        M = sympy.Matrix(2, 2, lambda i,j: 0.0)
        M[0,0] = V[0]
        M[0,1] = 0.5*V[2]
        M[1,0] = 0.5*V[2]
        M[1,1] = V[1]
    elif V.shape[0] == 6:
        M = sympy.Matrix(3, 3, lambda i,j : 0.0)
        M[0,0] = V[0]
        M[0,1] = 0.5*V[3]
        M[0,2] = 0.5*V[5]
        M[1,0] = 0.5*V[3]
        M[1,1] = V[1]
        M[1,2] = 0.5*V[4]
        M[2,0] = 0.5*V[5]
        M[2,1] = 0.5*V[4]
        M[2,2] = V[2]

    return M

def StrainToVoigt(M):
    """
    This method transform the strains matrix to Voigt notation.

    Keyword arguments:
    - M -- The strain matrix
    """
    if M.shape[0] == 2:
        vm = sympy.Matrix(3, 1, lambda i,j: 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = 2.0*M[0,1]
    elif M.shape[0] == 3:
        vm = sympy.Matrix(6, 1, lambda i,j : 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = M[2,2]
        vm[3,0] = 2.0*M[0,1]
        vm[4,0] = 2.0*M[1,2]
        vm[5,0] = 2.0*M[0,2]

    return vm

def VoigtToMatrix(V):
    """
    This method transforms a Voigt vector to matrix.

    Keyword arguments:
    - V -- The vector in Voigt notation
    """
    if V.shape[1] != 1:
        raise Exception(f"Provided array is not a vector. Shape is {V.shape}.")

    if V.shape[0] == 3:
        M = sympy.Matrix(2, 2, lambda i,j: 0.0)
        M[0,0] = V[0]
        M[0,1] = V[2]
        M[1,0] = V[2]
        M[1,1] = V[1]
    elif V.shape[0] == 6:
        M = sympy.Matrix(3, 3, lambda i,j : 0.0)
        M[0,0] = V[0]
        M[0,1] = V[3]
        M[0,2] = V[5]
        M[1,0] = V[3]
        M[1,1] = V[1]
        M[1,2] = V[4]
        M[2,0] = V[5]
        M[2,1] = V[4]
        M[2,2] = V[2]

    return M

def MatrixToVoigt(M):
    """
    This method transform a symmetric matrix to Voigt notation.

    Keyword arguments:
    - M -- The input matrix
    """
    if M.shape[0] == 2:
        vm = sympy.Matrix(3, 1, lambda i,j : 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = M[0,1]
    elif M.shape[0] == 3:
        vm = sympy.Matrix(6, 1, lambda i,j : 0.0)
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = M[2,2]
        vm[3,0] = M[0,1]
        vm[4,0] = M[1,2]
        vm[5,0] = M[0,2]

    return vm

def ConvertTensorToVoigtMatrix(C):
    """
    This method converts a 4th order tensor given to a matrix in Voigt notation.

    Keyword arguments:
    - C -- The 4th order tensor to be converted to Voigt notation
    """
    # Get input Voigt matrix strain size
    dim = C.shape[0]
    if dim != C.shape[1] or dim != C.shape[2] or dim != C.shape[3]:
        raise ValueError(f"Input 4th order tensor is not valid. Shape is ({dim},{C.shape[1]},{C.shape[2]},{C.shape[3]})")

    # Set the tensor to Voigt indices conversion dictionary from the input strain size
    conversion = _GetTensorToVoigtConversionIndices(dim)

    # Set and fill the Voigt notation matrix
    strain_size = 3 if dim == 2 else 6
    C_mat = sympy.zeros(strain_size,strain_size)
    for i in range(strain_size):
        index_1 = conversion[i]
        for j in range(strain_size):
            index_2 = conversion[j]
            C_mat[i,j] = C[index_1[0],index_1[1],index_2[0],index_2[1]]

    return C_mat

def ConvertVoigtMatrixToTensor(C):
    """
    This method converts a matrix given in Voigt notation in a 4th order tensor.

    Keyword arguments:
    - C -- The matrix in Voigt notation to be converted to 4th order tensor
    """
    # Get input Voigt matrix strain size
    strain_size = C.shape[0]
    if strain_size != C.shape[1]:
        raise ValueError("Input Voigt matrix is not square. Shape is ({},{})".format(strain_size, C.shape[1]))

    # Set the Voigt to tensor indices conversion dictionary from the input strain size
    conversion = _GetVoigtToTensorConversionIndices(strain_size)
    dim = 2 if strain_size == 3 else 3

    # Set and fill the symmetric fourth order tensor
    C_tensor = sympy.MutableDenseNDimArray(sympy.zeros(dim**4),shape=(dim,dim,dim,dim))
    for i in range(dim):
        for j in range(dim):
            index_1 = conversion[(i,j)]
            for k in range(dim):
                for l in range(dim):
                    index_2 = conversion[(k,l)]
                    C_tensor[i,j,k,l] = C[index_1,index_2]

    return C_tensor

def DoubleContraction(A,B):
    """
    This method performs the double contraction A:B.

    Keyword arguments:
    - A -- Left tensor
    - B -- Right tensor
    """
    rank_A = len(A.shape)
    AB = sympy.tensorproduct(A,B)
    tmp = sympy.tensorcontraction(AB, (rank_A-1,rank_A+1))
    output = sympy.tensorcontraction(tmp, (rank_A-2,rank_A-1))
    return output

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

def SubstituteMatrixValueByVoigtSymbols(where_to_substitute, what_to_substitute, substituted_value, is_strain):
    """
    This method substitutes values into a matrix.
    Note that this is intended to be used to substitute matrix indices by their corresponding Voigt ones

    Keyword arguments:
    - where_to_substitute -- Coordinates where to substitute
    - what_to_substitute -- Components to substitute
    - substituted_value -- Variable to substitute
    - is_strain -- Flag indicating if substituted_value is a strain
    """

    # Check the array to which the substitution is applied
    strain_size = where_to_substitute.shape[0]
    if (where_to_substitute.shape[1] != 1 and where_to_substitute.shape[1] != strain_size):
        raise Exception("Input array is not square (constitutive tensor in Voigt notation) not vector (stress vector in Voigt notation).")

    # Get the Voigt to tensor conversion indices
    conversion = _GetVoigtToTensorConversionIndices(strain_size)

    # Set the factor to be applied to each substitution term
    if is_strain:
        factor_func = lambda i,j : 1.0 if i==j else 2.0
    else:
        factor_func = lambda i,j : 1.0

    # Do the substitution
    dim = 2 if strain_size == 3 else 3
    for lll in range(where_to_substitute.shape[0]):
        for kkk in range(where_to_substitute.shape[1]):
            tmp = where_to_substitute[lll, kkk]
            for i in range(dim):
                for j in range(dim):
                    factor = factor_func(i, j)
                    voigt_index = conversion[(i, j)]
                    tmp = tmp.subs(what_to_substitute[i, j], factor * substituted_value[voigt_index])
            where_to_substitute[lll, kkk] = tmp

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

    For vectors:  `variable_3`   becomes `variable[3]`
    For matrices: `variable_3_7` becomes `variable[3,7]`
    For lists of matrices: `variable_0_3_7` becomes `variable[0][3,7]`

    The accessor for matrices is chosen according to the language (`[]` vs. `()`)
    """
    #Lists of matrices
    pattern = r"_(\d+)_(\d+)_(\d+)"
    replacement = r"[\1][\2,\3]" if language == 'python' else r"[\1](\2,\3)"
    expression = re.sub(pattern, replacement, expression)

    #Matrices
    pattern = r"_(\d+)_(\d+)"
    replacement = r"[\1,\2]" if language == 'python' else r"(\1,\2)"
    expression = re.sub(pattern, replacement, expression)

    # Vectors
    pattern = r"_(\d+)"
    replacement = r"[\1]"
    expression = re.sub(pattern, replacement, expression)

    return expression

def _GetVoigtToTensorConversionIndices(strain_size):
    if strain_size == 3:
        conversion = {
            (0, 0): 0,
            (1, 1): 1,
            (0, 1): 2,
            (1, 0): 2
        }
    elif strain_size == 6:
        conversion = {
            (0, 0) : 0,
            (1, 1) : 1,
            (2, 2) : 2,
            (0, 1) : 3,
            (1, 0) : 3,
            (1, 2) : 4,
            (2, 1) : 4,
            (2, 0) : 5,
            (0, 2) : 5
        }
    else:
        raise ValueError("Wrong strain size {}.".format(strain_size))

    return conversion

def _GetTensorToVoigtConversionIndices(dim):
    if dim == 2:
        conversion = {
            0: (0, 0),
            1: (1, 1),
            2: (0, 1)
        }
    elif dim == 3:
        conversion = {
            0: (0, 0),
            1: (1, 1),
            2: (2, 2),
            3: (0, 1),
            4: (1, 2),
            5: (0, 2)
        }
    else:
        raise ValueError("Wrong dimension {} in input tensor.".format(dim))

    return conversion

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
    - expression -- The expression to generate code from
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
