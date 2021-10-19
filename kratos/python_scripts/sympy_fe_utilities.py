from sympy import *

def DefineMatrix( name, m,n ):
    """ This method defines a symbolic matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    """
    return Matrix( m,n, lambda i,j: var(name+'_%d_%d' % (i,j)) )

def DefineSymmetricMatrix( name, m,n ):
    """ This method defines a symbolic symmetric matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    """
    tmp = DefineMatrix(name,m,n)

    # Impose symm
    for i in range(0,tmp.shape[0]):
        for j in range(i+1,tmp.shape[1]):
            tmp[j,i] = tmp[i,j]

    return tmp

def DefineVector( name, m):
    """ This method defines a symbolic vector

    Keyword arguments:
    name -- Name of variables.
    m -- Number of components.
    """
    return Matrix( m,1, lambda i,j: var(name+'_%d' % (i)) )

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False):
    """ This method defines shape functions and derivatives
    Note that partition of unity is imposed
    the name HAS TO BE --> N and DN

    Keyword arguments:
    nnodes -- Number of nodes
    dim -- Dimension of the space
    impose_partion_of_unity -- Impose the partition of unity

    Note that partition of unity is imposed the name HAS TO BE --> N and DN
    """
    DN = DefineMatrix('DN',nnodes,dim)
    N = DefineVector('N',nnodes)

    #impose partition of unity
    if(impose_partion_of_unity == True):
        N[nnodes-1] = 1
        for i in range(0,nnodes-1):
            N[nnodes-1] -= N[i]

        DN[nnodes-1,:] = -DN[0,:]
        for i in range(1,nnodes-1):
            DN[nnodes-1,:] -= DN[i,:]

    return N,DN

def StrainToVoigt(M):
    """ This method transform the strains matrix to Voigt notation

    Keyword arguments:
    M -- The strain matrix
    """
    if(M.shape[0] == 2):
        vm = Matrix( 3,1, zeros(3,1))
        vm[0,0] = M[0,0]
        vm[1,0] = M[1,1]
        vm[2,0] = 2.0*M[0,1]
    elif(M.shape[0] == 3):
        raise Exception("not implemented yet")
    return vm

def MatrixB(DN):
    """ This method defines the deformation matrix B

    Keyword arguments:
    DN -- The shape function derivatives
    """
    dim = DN.shape[1]
    if(dim == 2):
        strain_size = 3
        nnodes = DN.shape[0]
        B = Matrix( zeros(strain_size, nnodes*dim) )
        for i in range(0,nnodes):
            for k in range(0,dim):
                B[0,i*dim] = DN[i,0]; B[0, i*dim+1] = 0;
                B[1,i*dim] = 0;       B[1, i*dim+1] = DN[i,1];
                B[2,i*dim] = DN[i,1]; B[2, i*dim+1] = DN[i,0];
    elif(dim == 3):
        strain_size = 6
        nnodes = DN.shape[0]
        B = Matrix( zeros(strain_size, nnodes*dim) )
        for i in range(0,nnodes):
            B[ 0, i*3 ] = DN[ i, 0 ];
            B[ 1, i*3 + 1 ] = DN[ i, 1 ];
            B[ 2, i*3 + 2 ] = DN[ i, 2 ];
            B[ 3, i*3 ] = DN[ i, 1 ];
            B[ 3, i*3 + 1 ] = DN[ i, 0 ];
            B[ 4, i*3 + 1 ] = DN[ i, 2 ];
            B[ 4, i*3 + 2 ] = DN[ i, 1 ];
            B[ 5, i*3 ] = DN[ i, 2 ];
            B[ 5, i*3 + 2 ] = DN[ i, 0 ];
    else:
        print("dimension asked in Matrix B is ",dim)
        raise Exception("wrong dimension")
    return B

def grad_sym_voigtform(DN, x):
    """ This method defines a symmetric gradient

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    [nnodes, dim] = x.shape

    B = MatrixB(DN)

    # Put the x components one after the other in a vector
    xvec = Matrix( zeros(B.shape[1], 1 ) );
    for i in range(0,nnodes):
        for k in range(0,dim):
            xvec[i*dim+k] = x[i,k]

    return simplify( B*xvec )

def grad(DN,x):
    """ This method defines a gradient

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    error_msg = """
                    The function grad(DN,x) was removed to avoid misunderstandings.
                    You should use either DfjDxi(DN,f) or DfiDxj(DN,f).

                    DfjDxi(DN,f): returns a matrix D such that D(i,j) = D(fj)/D(xi) - the standard one in fluid dynamics
                    that is:
                        D(f1)/D(x1) D(f2)/D(x1) D(f3)/D(x1)
                        D(f1)/D(x2) D(f2)/D(x2) D(f3)/D(x2)
                        D(f1)/D(x3) D(f2)/D(x3) D(f3)/D(x3)


                    DfiDxj(DN,f): returns a matrix D such that D(i,j) = D(fi)/D(xj) - the standard one in structural mechanics
                    that is:
                        D(f1)/D(x1) D(f1)/D(x2) D(f1)/D(x3)
                        D(f2)/D(x1) D(f2)/D(x2) D(f2)/D(x3)
                        D(f3)/D(x1) D(f3)/D(x2) D(f3)/D(x3)

                    Note that the two gradients are one the transpose of the other.
    """
    raise Exception (error_msg)

def DfjDxi(DN,f):
    """ This method defines a gradient. This returns a matrix D such that D(i,j) = D(fj)/D(xi)

    Keyword arguments:
    DN -- The shape function derivatives
    f-- The variable to compute the gradient
    """
    return simplify(DN.transpose()*f)

def DfiDxj(DN,f):
    """ This method defines a gradient This returns a matrix D such that D(i,j) = D(fi)/D(xj)

    Keyword arguments:
    DN -- The shape function derivatives
    f -- The variable to compute the gradient
    """
    return (DfjDxi(DN,f)).transpose()

def div(DN,x):
    """ This method defines the divergence

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    if(DN.shape != x.shape):
        raise Exception("shapes are not compatible")

    div_x = 0
    for i in range(0,DN.shape[0]):
        for k in range(0,DN.shape[1]):
            div_x += DN[i,k]*x[i,k]

    return Matrix( [ simplify(div_x) ])

def SubstituteMatrixValue( where_to_substitute, what_to_substitute, substituted_value ):
    """ This method substitutes values into a matrix

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    for lll  in range(where_to_substitute.shape[0] ) :
        for kkk  in range(where_to_substitute.shape[1] ) :
            tmp  = where_to_substitute[lll,kkk]
            for i in range(what_to_substitute.shape[0]):
                for j in range(what_to_substitute.shape[1]):
                    #print("what to substitute ",what_to_substitute[i,j])
                    #print("substituted_value ",substituted_value[i,j])
                    tmp = tmp.subs( what_to_substitute[i,j], substituted_value[i,j] )

            where_to_substitute[lll,kkk] = tmp

    return where_to_substitute

def SubstituteScalarValue( where_to_substitute, what_to_substitute, substituted_value ):
    """ This method substitutes values into a scalar

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    for lll  in range(where_to_substitute.shape[0] ) :
        tmp  = where_to_substitute[lll]
        tmp = tmp.subs( what_to_substitute, substituted_value )
        where_to_substitute[lll] = tmp
    return where_to_substitute

def Compute_RHS(functional, testfunc, do_simplifications=False):
    """ This computes the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    do_simplifications -- If apply simplifications
    """
    rhs = Matrix( zeros(testfunc.shape[0],1) )
    for i in range(0,testfunc.shape[0]):
        rhs[i] = diff(functional[0,0], testfunc[i])

        if(do_simplifications):
            rhs[i] = simplify(rhs[i])

    return rhs

def Compute_LHS(rhs, testfunc, dofs, do_simplifications=False):
    """ This computes the LHS matrix

    Keyword arguments:
    rhs -- The RHS vector
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    lhs = Matrix( zeros(testfunc.shape[0],dofs.shape[0]) )
    for i in range(0,lhs.shape[0]):
        for j in range(0,lhs.shape[1]):
            lhs[i,j] = -diff(rhs[i,0], dofs[j,0])

            if(do_simplifications):
                lhs[i,j] = simplify(lhs[i,j])

    return lhs

def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications=False):
    """ This computes the LHS matrix and the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    rhs = Compute_RHS(functional, testfunc, do_simplifications)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    return rhs,lhs

def OutputVector(r, name, mode="python", initial_tabs=3, max_index=30, replace_indices=True, assignment_op="="):
    """ This method converts into text the RHS vector

    Keyword arguments:
    rhs -- The RHS vector
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    initial_spaces = str("")
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")

    outstring = str("")
    for i in range(0,r.shape[0]):
        if(mode == "python"):
            outstring += initial_spaces + name + str("[") + str(i) + str("]") + str(assignment_op) + str(r[i,0]) + str("\n")
        elif(mode=="c"):
            outstring += initial_spaces + name + str("[") + str(i) + str("]") + str(assignment_op) + ccode(r[i,0]) + str(";\n")

    if replace_indices:
        #matrix entries (two indices)
        for i in range(0,max_index):
            for j in range(0,max_index):
                if(mode == "python"):
                    replacement_string = str("[") + str(i) + str(",") + str(j) + str("]")
                elif(mode=="c"):
                    replacement_string = str("(") + str(i) + str(",") + str(j) + str(")")
                to_be_replaced = str("_") + str(i) + str("_") + str(j)
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

        #vector entries(one index)
        for i in range(0,max_index):
            replacement_string = str("[") + str(i) + str("]")
            to_be_replaced = str("_") + str(i)
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring


def OutputMatrix(lhs, name, mode, initial_tabs=3, max_index=30, replace_indices=True, assignment_op="="):
    """ This method converts into text the LHS matrix

    Keyword arguments:
    lhs -- The LHS matrix
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    initial_spaces = str("")
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")

    outstring = str("")
    for i in range(0,lhs.shape[0]):
        for j in range(0,lhs.shape[1]):
            if(mode == "python"):
                outstring += initial_spaces + name + str("[") + str(i) + str(",") + str(j) + str("]") + str(assignment_op) + str(lhs[i,j]) + str("\n")
            elif(mode=="c"):
                outstring += initial_spaces + name + str("(") + str(i) + str(",") + str(j) + str(")") + str(assignment_op) + ccode(lhs[i,j]) + str(";\n")

    if replace_indices:
        #matrix entries (two indices)
        for i in range(0,max_index):
            for j in range(0,max_index):
                if(mode == "python"):
                    replacement_string = str("[") + str(i) + str(",") + str(j) + str("]")
                elif(mode=="c"):
                    replacement_string = str("(") + str(i) + str(",") + str(j) + str(")")
                to_be_replaced = str("_") + str(i) + str("_") + str(j)
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

        #vector entries(one index)
        for i in range(0,max_index):
            replacement_string = str("[") + str(i) + str("]")
            to_be_replaced = str("_") + str(i)
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

def OutputSymbolicVariable(var, mode="python", initial_tabs=3, max_index=30, replace_indices=True):
    """ This method converts into text the LHS matrix (only non-zero terms)

    Keyword arguments:
    var -- The variable to define symbolic
    mode -- The mode of output
    varname -- The name of the variables
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    replace_indices -- If the indixes must be replaced
    """
    initial_spaces = str("")
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")

    outstring = str("")
    if(mode == "python"):
        outstring += initial_spaces+str(var)+str("\n")
    elif(mode=="c"):
        outstring += initial_spaces+ccode(var)+str(";\n")

    if replace_indices:
        #matrix entries (two indices)
        for i in range(0,max_index):
            for j in range(0,max_index):
                if(mode == "python"):
                    replacement_string = str("[") + str(i) + str(",") + str(j) + str("]")
                elif(mode=="c"):
                    replacement_string = str("(") + str(i) + str(",") + str(j) + str(")")
                to_be_replaced = str("_") + str(i) + str("_") + str(j)
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

        #vector entries (one index)
        for i in range(0,max_index):
            replacement_string = str("[") + str(i) + str("]")
            to_be_replaced = str("_") + str(i)
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

def OutputMatrix_CollectingFactors(A, name, mode, initial_tabs=3, max_index=30, optimizations='basic', replace_indices=True, assignment_op="="):
    """ This method collects the constants of the replacement for matrices

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    symbol_name = "c"+name
    A_factors, A_collected = cse(A,numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] #overwrite lhs with the one with the collected components

    Acoefficient_str = str("")
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode,initial_tabs, max_index, replace_indices)
        Acoefficient_str += "const double " + str(varname.__str__()) + " = " + output_value
        #print(output_str)
    A_out = Acoefficient_str + OutputMatrix(A, name, mode, initial_tabs, max_index, replace_indices, assignment_op)
    return A_out

def OutputVector_CollectingFactors(A, name, mode, initial_tabs=3, max_index=30, optimizations='basic', replace_indices=True, assignment_op="="):
    """ This method collects the constants of the replacement for vectors

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    replace_indices -- If the indixes must be replaced
    assignment_op -- The assignment operation
    """
    symbol_name = "c"+name
    A_factors, A_collected = cse(A,numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] #overwrite lhs with the one with the collected components

    Acoefficient_str = str("")
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode, initial_tabs, max_index, replace_indices)
        Acoefficient_str += "const double " + str(varname.__str__()) + " = " + output_value
        #print(output_str)
    A_out = Acoefficient_str + OutputVector(A, name, mode, initial_tabs, max_index, replace_indices, assignment_op)
    return A_out
