# Import base sympy_fe_utilities
import KratosMultiphysics
import KratosMultiphysics.sympy_fe_utilities as sympy_fe_utilities

# Import sympy
import sympy

#TODO: Define base utilities as a class to avoid redefine mehods

def GetSympyVersion():
    """ This method returns the current Sympy version
    """
    return sympy.__version__

def DefineMatrix( name, m, n, mode = "Function" ):
    """ This method defines a symbolic matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    mode -- The type of variable is defined (function or symbol)
    """
    if float(GetSympyVersion()) <= 1.2:
        return sympy.Matrix( m,n, lambda i,j: sympy.var(name+'_%d_%d' % (i,j)) )
    else:
        if mode == "Symbol":
            return sympy.Matrix( m,n, lambda i,j: sympy.var(name+'_%d_%d' % (i,j)) )
        elif mode == "Function":
            return sympy.Matrix( m,n, lambda i,j: sympy.symbols(name+'_%d_%d' % (i,j), cls=sympy.Function))
        else:
            raise Exception("Not implemented yet")

def DefineSymmetricMatrix( name, m, n = -1, mode = "Symbol"):
    """ This method defines a symbolic symmetric matrix

    Keyword arguments:
    name -- Name of variables.
    m -- Number of rows.
    n -- Number of columns.
    """
    # Assuming square matrix
    if n < 0:
        n =  m
    tmp = DefineMatrix(name,m,n)

    # Impose symm
    for i in range(0,tmp.shape[0]):
        for j in range(i+1,tmp.shape[1]):
            tmp[j,i] = tmp[i,j]

    return tmp

def DefineVector( name, m, mode = "Symbol"):
    """ This method defines a symbolic vector

    Keyword arguments:
    name -- Name of variables.
    m -- Number of components.
    mode -- The type of variable is defined (function or symbol)
    """
    if float(GetSympyVersion()) <= 1.2:
        return sympy.Matrix( m,1, lambda i,j: sympy.var(name+'_%d' % (i)) )
    else:
        if mode == "Symbol":
            return sympy.Matrix( m,1, lambda i,j: sympy.var(name+'_%d' % (i)) )
        elif mode == "Function":
            return sympy.Matrix( m,1, lambda i,j: sympy.symbols(name+'_%d' % (i), cls=sympy.Function))
        else:
            raise Exception("Not implemented yet")

def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity = False):
    """ This method defines shape functions and derivatives
    Note that partition of unity is imposed
    the name HAS TO BE --> N and DN

    Keyword arguments:
    nnodes -- Number of nodes
    dim -- Dimension of the space
    impose_partion_of_unity -- Impose the partition of unity
    """
    return sympy_fe_utilities.DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

def DefineCustomShapeFunctions(nnodes, dim, name):
    """ This method defines shape functions and derivatives (custom version)

    Keyword arguments:
    nnodes -- Number of nodes
    dim -- Dimension of the space
    name -- Name of the variable for the shape functions
    """
    DN = DefineMatrix('D'+name,nnodes,dim)
    N = DefineVector(name ,nnodes)

    return N,DN

def StrainToVoigt(M):
    """ This method transform the strains matrix to Voigt notation

    Keyword arguments:
    M -- The strain matrix
    """
    return sympy_fe_utilities.StrainToVoigt(M)

def MatrixB(DN):
    """ This method defines the deformation matrix B

    Keyword arguments:
    DN -- The shape function derivatives
    """
    return sympy_fe_utilities.MatrixB(DN)

def grad_sym_voigtform(DN, x):
    """ This method defines a symmetric gradient

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    return sympy_fe_utilities.grad_sym_voigtform(DN, x)

def grad(DN,x):
    """ This method defines a gradient

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    return sympy_fe_utilities.grad(DN, x)

def DfjDxi(DN,f):
    """ This method defines a gradient. This returns a matrix D such that D(i,j) = D(fj)/D(xi)

    Keyword arguments:
    DN -- The shape function derivatives
    f-- The variable to compute the gradient
    """
    return sympy_fe_utilities.DfjDxi(DN, f)

def DfiDxj(DN,f):
    """ This method defines a gradient This returns a matrix D such that D(i,j) = D(fi)/D(xj)

    Keyword arguments:
    DN -- The shape function derivatives
    f -- The variable to compute the gradient
    """
    return sympy_fe_utilities.DfiDxj(DN, f)

def div(DN,x):
    """ This method defines the divergence

    Keyword arguments:
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """
    return sympy_fe_utilities.div(DN, x)

def DefineJacobian(J, DN, x):
    """ This method defines a Jacobian

    Keyword arguments:
    J -- The Jacobian matrix
    DN -- The shape function derivatives
    x -- The variable to compute the gradient
    """

    [nnodes, dim] = x.shape
    localdim = dim - 1

    if (dim == 2):
        if (nnodes == 2):
            J[0,0] = 0.5 * (x[1,0] - x[0,0])
            J[1,0] = 0.5 * (x[1,1] - x[0,1])
    else:
        if (nnodes == 3):
            J[0,0] = - (x[0,0] + x[1,0])
            J[1,0] = - (x[0,1] + x[1,1])
            J[2,0] = - (x[0,2] + x[1,2])
            J[0,1] = - (x[0,0] + x[2,0])
            J[1,1] = - (x[0,1] + x[2,1])
            J[2,1] = - (x[0,2] + x[2,2])
        else:
            for i in range(dim):
                for j in range(localdim):
                    J[i, j] =  0

            for i in range(nnodes):
                for k in range(dim):
                    for m in range(localdim):
                        J[k,m] += x[i,k] * DN[i,m]

    return J

def DefineCalculateNormals(normal, tangent1, tangent2, J):
    """ This method calculates the normals and tangents

    Keyword arguments:
    normal -- The normal vector
    tangent1 -- The first tangent vector
    tangent2 -- The second tangent vector
    J -- The Jacobian matrix
    """

    [nnodes, dim] = normal.shape

    for node in range(nnodes):
        norm = sympy.sqrt((J.col(0).transpose()*J.col(0))[0,0])
        tangent1[dim*node] = (J.col(0)).transpose()/norm
        if dim == 2:
            normal[node,0] =   tangent1[node,1]
            normal[node,1] = - tangent1[node,0]
        else:
            norm = sympy.sqrt((J.col(1).transpose()*J.col(1))[0,0])
            tangent2[dim*node] = (J.col(1)).transpose()/norm
            normal[dim*node] = (tangent1.row(node)).cross(tangent2.row(node))

    return normal, tangent1, tangent2

def CreateVariableMatrixList(variable_list, variable_matrix):
    """ This method creates a variable list from variable matrix

    Keyword arguments:
    variable_list -- The variable list to be created
    variable_matrix -- The variable matrix
    """
    nnodes = variable_matrix.shape[0]
    dim = variable_matrix.shape[1]
    for i in range(0,nnodes):
        for k in range(0,dim):
            variable_list.append(variable_matrix[i,k])

def CreateVariableVectorList(variable_list, variable_vector):
    """ This method creates a variable list from variable vector

    Keyword arguments:
    variable_list -- The variable list to be created
    variable_vector -- The variable vector
    """
    nnodes = variable_vector.shape[0]
    for i in range(0,nnodes):
        variable_list.append(variable_vector[i])

def DefineDofDependencyScalar(scalar, variable_list):
    """ This method injects a dependency into a scalar from a variable list

    Keyword arguments:
    scalar -- The scalar to be injected
    variable_list -- The variable list injected
    """
    return scalar(*variable_list)

def DefineDofDependencyVector(vector, variable_list):
    """ This method injects a dependency into a vector from a variable list

    Keyword arguments:
    vector -- The vector to be injected
    variable_list -- The variable list injected
    """
    for i in range(0,vector.shape[0]):
        vector[i, 0] = DefineDofDependencyScalar(vector[i, 0], variable_list)
    return vector

def DefineDofDependencyMatrix(matrix, variable_list):
    """ This method injects a dependency into a matrix from a variable list

    Keyword arguments:
    matrix -- The matrix to be injected
    variable_list -- The variable list injected
    """
    for i in range(0,matrix.shape[0]):
        for k in range(0,matrix.shape[1]):
            matrix[i, k] = DefineDofDependencyScalar(matrix[i, k], variable_list)
    return matrix

def SubstituteMatrixValue( where_to_substitute, what_to_substitute, substituted_value ):
    """ This method substitutes values into a matrix

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    return sympy_fe_utilities.SubstituteMatrixValue( where_to_substitute, what_to_substitute, substituted_value )

def SubstituteScalarValue( where_to_substitute, what_to_substitute, substituted_value ):
    """ This method substitutes values into a scalar

    Keyword arguments:
    where_to_substitute -- Coordinates where to substitute
    what_to_substitute -- Components to substitute
    substituted_value -- Variable to substitute
    """
    return sympy_fe_utilities.SubstituteScalarValue( where_to_substitute, what_to_substitute, substituted_value )

def GetShapeFunctionDefinitionLine2D2N(x,xg):
    """ This computes the shape functions on 2D line

    Keyword arguments:
    x -- Definition of line
    xg -- Gauss point
    """
    N = sympy.zeros(2)
    N[0] = (-x[1,0]+ xg[0]-x[1,1]+xg[1]) / (x[0,0] - x[1,0] + x[0,1] - x[1,1])
    N[1] = 1-N[0]

    return N

def GetShapeFunctionDefinitionLine3D3N(x,xg):
    """ This computes the shape functions on 3D line

    Keyword arguments:
    x -- Definition of line
    xg -- Gauss point
    """
    N = sympy.zeros(3)
    N[1] = -(((x[1,2]-x[2,2])*(x[2,0]+x[2,1]-xg[0]-xg[1])-(x[1,0]+x[1,1]-x[2,0]-x[2,1])*(x[2,2]-xg[2]))/(-(x[1,0]+x[1,1]-x[2,0]-x[2,1])*(x[0,2]-x[2,2])+(x[0,0]+x[0,1]-x[2,0]-x[2,1])*(x[1,2]-x[2,2])))
    N[2] = -((x[0,2]*x[2,0]+x[0,2]*x[2,1]-x[0,0]*x[2,2]-x[0,1]*x[2,2]-x[0,2]*xg[0]+x[2,2]*xg[0]-x[0,2]*xg[1]+x[2,2]*xg[1]+x[0,0]*xg[2]+x[0,1]*xg[2]-x[2,0]*xg[2]-x[2,1]*xg[2])/(x[0,2]*x[1,0]+x[0,2]*x[1,1]-x[0,0]*x[1,2]-x[0,1]*x[1,2]-x[0,2]*x[2,0]+x[1,2]*x[2,0]-x[0,2]*x[2,1]+x[1,2]*x[2,1]+x[0,0]*x[2,2]+x[0,1]*x[2,2]-x[1,0]*x[2,2]-x[1,1]*x[2,2]))
    N[0] = 1 - N[1] -N[2]

    return N

def Compute_RHS(functional, testfunc, do_simplifications = False):
    """ This computes the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    do_simplifications -- If apply simplifications
    """
    return sympy_fe_utilities.Compute_RHS(functional, testfunc, do_simplifications)

def Compute_LHS(rhs, testfunc, dofs, do_simplifications = False):
    """ This computes the LHS matrix

    Keyword arguments:
    rhs -- The RHS vector
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    return sympy_fe_utilities.Compute_LHS(rhs, testfunc, dofs, do_simplifications)

def Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications = False):
    """ This computes the LHS matrix and the RHS vector

    Keyword arguments:
    functional -- The functional to derivate
    testfunc -- The test functions
    dofs -- The dofs vectors
    do_simplifications -- If apply simplifications
    """
    return sympy_fe_utilities.Compute_RHS_and_LHS(functional, testfunc, dofs, do_simplifications)

# TODO: Unify with base utilities
def OutputVector(rhs, name, mode="python", initial_tabs = 1,max_index=30, aux_dict={}):
    """ This method converts into text the RHS vector

    Keyword arguments:
    rhs -- The RHS vector
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    aux_dict -- Auxiliar dictionary
    """
    initial_spaces = ""
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")

    outstring = ""
    for i in range(0,rhs.shape[0]):

        if mode == "python":
            outstring += initial_spaces + name + "[" +str(i)+ "]=" +str(rhs[i,0])+"\n"
        elif mode=="c":
            if "// Not supported in C:" in sympy.ccode(rhs[i,0]):
                var = rhs[i,0]
                if  "Derivative" in str(var):
                    for constantname, constantexp in aux_dict.items():
                        var = var.replace(constantname, constantexp)
                    aux_string = str(var)
                    outstring += initial_spaces + name +  "[" +str(i)+ "]=" + "//subsvar_" + aux_string.replace("// Not supported in C:","")+ ";\n"
                else:
                    outstring += initial_spaces + name +  "[" +str(i)+ "]=" + "//subsvar_" + sympy.ccode(rhs[i,0]).split("\n",2)[2]+ ";\n"
            else:
                outstring += initial_spaces + name + "[" +str(i)+ "]=" + sympy.ccode(rhs[i,0])+ ";\n"

    # Matrix entries (two indices)
    for i in range(0,max_index):
        for j in range(0,max_index):
            if mode == "python":
                replacement_string =  "[" +str(i)+ "," +str(j)+ "]"
            elif mode == "c":
                replacement_string =  "(" +str(i)+ "," +str(j)+ ")"
            to_be_replaced  =  "_" +str(i)+ "_" +str(j)
            if not "//subsvar_" in outstring:
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

    # Vector entries(one index(
    for i in range(0,max_index):
        replacement_string =  "[" +str(i)+ "]"
        to_be_replaced  =  "_" +str(i)
        if not "//subsvar_" in outstring:
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

# TODO: Unify with base utilities
def OutputMatrix(lhs, name, mode, initial_tabs = 1, max_index=30, aux_dict={}):
    """ This method converts into text the LHS matrix

    Keyword arguments:
    lhs -- The LHS matrix
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    aux_dict -- Auxiliar dictionary
    """
    initial_spaces = ""
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")
    outstring = ""
    for i in range(0,lhs.shape[0]):
        for j in range(0,lhs.shape[1]):
            if mode == "python":
                outstring += initial_spaces+name +  "[" +str(i)+ "," +str(j)+ "]=" +str(lhs[i,j])+"\n"
            elif mode=="c":
                if "// Not supported in C:" in sympy.ccode(lhs[i,j]):
                    var = lhs[i,j]
                    if  "Derivative" in str(var):
                        for constantname, constantexp in aux_dict.items():
                            var = var.replace(constantname, constantexp)
                        aux_string = str(var)
                        outstring += initial_spaces + name +  "(" +str(i)+ "," + str(j) + ")=" +"//subsvar_" + aux_string.replace("// Not supported in C:","")+ ";\n"
                    else:
                        outstring += initial_spaces + name +  "(" +str(i)+ "," + str(j) + ")=" +"//subsvar_" + sympy.ccode(lhs[i,j]).split("\n",2)[2]+ ";\n"
                else:
                    outstring += initial_spaces + name +  "(" +str(i)+ "," + str(j) + ")=" +sympy.ccode(lhs[i,j]) + ";\n"

    # Matrix entries (two indices)
    for i in range(0,max_index):
        for j in range(0,max_index):
            if mode == "python":
                replacement_string =  "[" +str(i)+ "," +str(j)+ "]"
            elif mode=="c":
                replacement_string =  "(" +str(i)+ "," +str(j)+ ")"
            to_be_replaced  =  "_" +str(i)+ "_" +str(j)
            if not "//subsvar_" in outstring:
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

    # Vector entries(one index(
    for i in range(0,max_index):
        replacement_string =  "[" + str(i) + "]"
        to_be_replaced  =  "_" + str(i)
        if not "//subsvar_" in outstring:
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

# TODO: Unify with base utilities
def OutputVectorNonZero(rhs,name, mode="python", initial_tabs = 1,max_index=30,aux_dict={}):
    """ This method converts into text the RHS vector (only non-zero terms)

    Keyword arguments:
    rhs -- The RHS vector
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    aux_dict -- Auxiliar dictionary
    """
    initial_spaces = ""
    for i in range(0,initial_tabs):
        initial_spaces += "    "

    outstring = ""
    for i in range(0,rhs.shape[0]):

        if mode == "python":
            outstring += initial_spaces + name + "[" + str(i) + "]+=" + str(rhs[i,0]) + "\n"
        elif mode=="c":
            if rhs[i] != 0:
                if "// Not supported in C:" in sympy.ccode(rhs[i,0]):
                    var = rhs[i,0]
                    if  "Derivative" in str(var):
                        for constantname, constantexp in aux_dict.items():
                            var = var.replace(constantname, constantexp)
                        aux_string = str(var)
                        outstring += initial_spaces + name +  "[" + str(i) + "]="  + "//subsvar_"+aux_string.replace("// Not supported in C:","")+ ";\n"
                    else:
                        outstring += initial_spaces + name +  "[" + str(i) + "]+=" + "//subsvar_" + sympy.ccode(rhs[i,0]).split("\n",2)[2]+ ";\n"
                else:
                    outstring += initial_spaces + name +  "[" +str(i)+ "]+=" + sympy.ccode(rhs[i,0])+ ";\n"

    # Matrix entries (two indices)
    for i in range(0, max_index):
        for j in range(0, max_index):
            if mode == "python":
                replacement_string =  "[" + str(i) + "," + str(j) + "]"
            elif mode=="c":
                replacement_string =  "(" + str(i) + "," + str(j) + ")"
            to_be_replaced  =  "_" +str(i)+ "_" +str(j)
            if not "//subsvar_" in outstring:
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

    # Vector entries(one index)
    for i in range(0,max_index):
        replacement_string =  "[" +str(i)+ "]"
        to_be_replaced  =  "_" +str(i)
        if not "//subsvar_" in outstring:
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

# TODO: Unify with base utilities
def OutputMatrixNonZero(lhs, name, mode, initial_tabs = 1, max_index=30,aux_dict={}):
    """ This method converts into text the LHS matrix (only non-zero terms)

    Keyword arguments:
    lhs -- The LHS matrix
    name -- The name of the variables
    mode -- The mode of output
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    aux_dict -- Auxiliar dictionary
    """
    initial_spaces = ""
    for i in range(0,initial_tabs):
        initial_spaces += str("    ")
    outstring = ""
    for i in range(0,lhs.shape[0]):
        for j in range(0,lhs.shape[1]):
            if mode == "python":
                outstring += initial_spaces+name + "[" + str(i) + "," + str(j) + "]=" + str(lhs[i,j]) + "\n"
            elif mode=="c":
                if lhs[i,j] != 0:
                    if "// Not supported in C:" in sympy.ccode(lhs[i,j]):
                        var = lhs[i,j]
                        if  "Derivative" in str(var):
                            for constantname, constantexp in aux_dict.items():
                                var = var.replace(constantname, constantexp)
                            aux_string = str(var)
                            outstring += initial_spaces + name + "(" +str(i)+ "," +str(j) + ")+=" + "//subsvar_" + aux_string.replace("// Not supported in C:","")+ ";\n"
                        else:
                            outstring += initial_spaces + name + "(" +str(i)+ "," +str(j) + ")+=" + "//subsvar_" + sympy.ccode(lhs[i,j]).split("\n",2)[2]+ ";\n"
                    else:
                        outstring += initial_spaces + name + "(" +str(i)+ "," +str(j) + ")+=" + sympy.ccode(lhs[i,j])+ ";\n"

    # Matrix entries (two indices)
    for i in range(0,max_index):
        for j in range(0,max_index):
            if mode == "python":
                replacement_string =  "[" + str(i) + "," + str(j) + "]"
            elif mode=="c":
                replacement_string =  "(" + str(i) + "," + str(j) + ")"
            to_be_replaced  =  "_" +str(i) + "_" +str(j)
            if not "//subsvar_" in outstring:
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring

    # Vector entries(one index(
    for i in range(0,max_index):
        replacement_string =  "[" +str(i)+ "]"
        to_be_replaced  =  "_" +str(i)
        if not "//subsvar_" in outstring:
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

# TODO: Unify with base utilities
def OutputSymbolicVariable(var, mode="python", varname = "",aux_dict={}, initial_tabs = 1,max_index=30):
    """ This method converts into text the LHS matrix (only non-zero terms)

    Keyword arguments:
    var -- The variable to define symbolic
    mode -- The mode of output
    varname -- The name of the variables
    aux_dict -- Auxiliar dictionary
    initial_tabs -- The number of tabulations considered
    max_index -- The maximum index
    """
    initial_spaces = ""
    for i in range(0,initial_tabs):
        initial_spaces += "    "

    outstring = ""

    if mode == "python":
        outstring += initial_spaces + str(var) + "\n"
    elif mode == "c":
        if "// Not supported in C:" in sympy.ccode(var):
            if "Derivative" in str(var):
                for constantname, constantexp in aux_dict.items():
                    var = var.replace(constantname, constantexp)
                aux_string = str(var)
                outstring += initial_spaces + "//subsvar_" + aux_string+ "; // " + aux_string[:].upper() + "\n"
            else:
                aux_string = str(var)
                aux_dict[varname]= aux_string
                outstring += initial_spaces+ " //subsvar_" + aux_string+ "; // " + aux_string[:].upper() + "\n"
        else:
            outstring += initial_spaces + sympy.ccode(var)+ ";\n"

    outstring = SubstituteIndex(outstring, mode, max_index)

    return outstring

def SubstituteIndex(outstring, mode="python",max_index=30):
    """ This method substitutes certain indexes

    Keyword arguments:
    outstring -- The resulting string
    mode -- The mode of output
    max_index -- The maximum index
    """
    # Matrix entries (two indices)
    for i in range(0,max_index):
        for j in range(0,max_index):
            if mode == "python":
                replacement_string =  "[" +str(i)+ "," +str(j)+ "]"
            elif mode=="c":
                replacement_string =  "(" +str(i)+ "," +str(j)+ ")"
            to_be_replaced  =  "_" +str(i)+ "_" +str(j)
            if not "//subsvar_" in outstring:
                newstring = outstring.replace(to_be_replaced, replacement_string)
                outstring = newstring
    # Vector entries(one index(
    for i in range(0,max_index):
        replacement_string =  "[" +str(i)+ "]"
        to_be_replaced  =  "_" +str(i)
        if not "//subsvar_" in outstring:
            newstring = outstring.replace(to_be_replaced, replacement_string)
            outstring = newstring

    return outstring

def DefineVariableLists(variable, name, replacement, dependency, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, typology):
    """ This method define the variables lists

    Keyword arguments:
    variable -- The variable to define
    name -- The name to consider
    replacement --  The replacement
    dependency --  The dependency
    var_strings --  The strings that define the variables
    var_strings --  The strings that substitute the variables
    var_strings_aux_subs --  The strings that auxiliarly substitute the variables
    der_var_strings --  The strings that define the derivatives of the strings
    der_var_list --  The list of derivatives variables
    typology --  The typology considered
    """
    if typology == "scalar":
        var_strings.append(str(variable))
        var_strings_subs.append(replacement)
        var_strings_aux_subs.append(name)
        count = 0
        for dof in dependency:
            der_var_strings.append(str(sympy.diff(variable, dof)))
            der_var_list.append("Delta"+name+str(count))
            count += 1
    elif typology == "vector":
        [nnodes, dim] = variable.shape
        for node in range(nnodes):
            var_strings.append(str(variable[node]))
            var_strings_subs.append(replacement+"["+str(node)+"]")
            var_strings_aux_subs.append(name+"_"+str(node))
            count = 0
            for dof in dependency:
                der_var_strings.append(str(sympy.diff(variable[node], dof)))
                der_var_list.append("Delta"+name+str(count)+"_"+str(node))
                count += 1
    else:
        [nnodes, dim] = variable.shape
        for node in range(nnodes):
            for i in range(dim):
                var_strings.append(str(variable[node,i]))
                var_strings_subs.append(replacement+"("+str(node)+","+str(i)+")")
                var_strings_aux_subs.append(name+"_"+str(node)+"_"+str(i))
                count = 0
                for dof in dependency:
                    der_var_strings.append(str(sympy.diff(variable[node,i], dof)))
                    der_var_list.append("Delta"+name+str(count)+"_"+str(node)+"_"+str(i))
                    count += 1

    return var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list

def Derivatives_CollectingFactors(A,name, mode, initial_tabs = 1, max_index=30, optimizations='basic', postprocess=None, order='canonical'):
    """ This method collects the constants of the replacement for derivatives
    TODO Think about collecting factors also in the derivatives

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    postprocess -- The type of postprocess
    order -- The type of order
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A, sympy.numbered_symbols(symbol_name), optimizations, postprocess, order)
    A = A_collected

    Acoefficient_str = ""
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        Acoefficient_str += "    const double " + varname.__str__() + " = " + value

    return [A, Acoefficient_str]

# TODO: Unify with base utilities
def OutputMatrix_CollectingFactors(A,name, mode, initial_tabs = 1, max_index=30, optimizations='basic'):
    """ This method collects the constants of the replacement for matrices

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A,sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] # Overwrite lhs with the one with the collected components

    aux_dict = {}

    Acoefficient_str = ""
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode, varname, aux_dict)
        Acoefficient_str += "    const double " + varname.__str__() + " = " + output_value
    A_out = Acoefficient_str+"\n" + OutputMatrix(A, name, mode, initial_tabs, max_index, aux_dict)
    return A_out

# TODO: Unify with base utilities
def OutputVector_CollectingFactors(A,name, mode, initial_tabs = 1, max_index=30, optimizations='basic'):
    """ This method collects the constants of the replacement for vectors

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A,sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] # Overwrite lhs with the one with the collected components

    aux_dict = {}

    Acoefficient_str = ""
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode, varname, aux_dict)
        Acoefficient_str += "    const double " + varname.__str__() + " = " + output_value
    A_out = Acoefficient_str + "\n" + OutputVector(A, name, mode, initial_tabs, max_index, aux_dict)
    return A_out

# TODO: Unify with base utilities
def OutputMatrix_CollectingFactorsNonZero(A,name, mode, initial_tabs = 1, max_index=30, optimizations='basic'):
    """ This method collects the constants of the replacement for matrices (only non-zero terms)

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A,sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] # Overwrite lhs with the one with the collected components

    aux_dict = {}

    Acoefficient_str = ""
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode, varname, aux_dict)
        Acoefficient_str += "    const double " + varname.__str__() + " = " + output_value
    A_out = Acoefficient_str + "\n" + OutputMatrixNonZero(A, name, mode, initial_tabs, max_index, aux_dict)
    return A_out

# TODO: Unify with base utilities
def OutputVector_CollectingFactorsNonZero(A,name, mode, initial_tabs = 1, max_index=30, optimizations='basic'):
    """ This method collects the constants of the replacement for vectors (only non-zero terms)

    Keyword arguments:
    A -- The  factors
    name -- The name of the constant
    mode -- The mode of replacement
    initial_tabs -- The number of initial tabulations
    max_index -- The max number of indexes
    optimizations -- The level of optimizations
    """
    symbol_name = "c" + name
    A_factors, A_collected = sympy.cse(A, sympy.numbered_symbols(symbol_name), optimizations)
    A = A_collected[0] # Overwrite lhs with the one with the collected components

    aux_dict = {}

    Acoefficient_str = ""
    for factor in A_factors:
        varname = factor[0]
        value = factor[1]
        output_value = OutputSymbolicVariable(value, mode, varname, aux_dict)
        Acoefficient_str += "    const double " + varname.__str__() + " = " + output_value
    A_out = Acoefficient_str + "\n" + OutputVectorNonZero(A, name, mode, initial_tabs, max_index, aux_dict)
    return A_out
