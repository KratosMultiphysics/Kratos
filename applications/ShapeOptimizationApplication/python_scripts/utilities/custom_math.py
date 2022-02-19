# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics as KM

# Import additional libraries
import math

# ==============================================================================
def SafeConvertVectorToMatrix(_L):
    if IsVector(_L):
        return [_L]
    else:
        return _L

# ------------------------------------------------------------------------------
def IsVector(_X):
    if len(_X)==0 or isinstance(_X[0],(float,int)):
        return True
    else:
        return False

# ------------------------------------------------------------------------------
def IsEmpty(_A):
    _A = SafeConvertVectorToMatrix(_A)

    if len(_A)==0 or len(_A[0])==0:
        return True
    return False

# ------------------------------------------------------------------------------
def ZeroVector(m):
    return [0.0 for i in range(m)]

# ------------------------------------------------------------------------------
def ZeroMatrix(m,n):
    A = [[0.0 for i in range(m)] for j in range(n)]
    return A

# ------------------------------------------------------------------------------
def Ones(n):
    return [1.0 for i in range(n)]

# ------------------------------------------------------------------------------
def RowSize(_A):
    _A = SafeConvertVectorToMatrix(_A)

    if len(_A)==0:
        return 0
    return len(_A[0])

# ------------------------------------------------------------------------------
def CollSize(_A):
    _A = SafeConvertVectorToMatrix(_A)

    return len(_A)

# ------------------------------------------------------------------------------
# horizontal concatenation
def HorzCat(_A,_B):
    _A = SafeConvertVectorToMatrix(_A)
    _B = SafeConvertVectorToMatrix(_B)

    if IsEmpty(_B):
        return _A
    if IsEmpty(_A):
        return _B

    return _A+_B

# ------------------------------------------------------------------------------
# vertical concatenation
def VertCat(_A,_B):
    _A = SafeConvertVectorToMatrix(_A)
    _B = SafeConvertVectorToMatrix(_B)

    if IsEmpty(_B):
        return _A
    if IsEmpty(_A):
        return _B

    coll_size_A = CollSize(_A)
    coll_size_B = CollSize(_B)

    if coll_size_A!=coll_size_B:
        raise ValueError("custom_math::VertCat: Wrong size in vertical concatenation detected!")

    _C = []
    for i in range(coll_size_A):
        _C.append(_A[i]+_B[i])

    if CollSize(_C) == 1:
        return _C[0]
    else:
        return _C

# ------------------------------------------------------------------------------
def Norm2(_X):
    temp_vec = [x**2 for x in _X]
    temp_sum = sum(temp_vec)
    return math.sqrt(temp_sum)

# ------------------------------------------------------------------------------
def NormInf3D(_X):
    temp_vec = [_X[3*i]**2 + _X[3*i+1]**2 + _X[3*i+2]**2 for i in range(int(len(_X)/3))]
    max_squared_value = max(temp_vec)
    return math.sqrt(max_squared_value)

# ------------------------------------------------------------------------------
def Dot(_X, _Y):
    if len(_X) != len(_Y):
        raise RuntimeError("custom_math::Dot: Dot product to be computed but _X and _Y do not have the same dimension!")
    return sum( [_X[i]*_Y[i] for i in range(len(_X))] )

# ------------------------------------------------------------------------------
def Plus(_X, _Y):
    if len(_X)!=len(_Y):
        raise ValueError("custom_math::Plus: Wrong size of input variables!")
    return [ _X[i]+_Y[i] for i in range(len(_X)) ]

# ------------------------------------------------------------------------------
def ScalarVectorProduct(scal, _X):
    return [x*scal for x in _X]

# ------------------------------------------------------------------------------
def ScalarMatrixProduct(scal, _A):
    return [ [ scal*_A[i][j] for j in range(RowSize(_A)) ] for i in range(CollSize(_A)) ]

# ------------------------------------------------------------------------------
def Prod(_A,_B):
    _A = SafeConvertVectorToMatrix(_A)
    _B = SafeConvertVectorToMatrix(_B)

    row_size_A = RowSize(_A)
    row_size_B = RowSize(_B)
    coll_size_A = CollSize(_A)
    coll_size_B = CollSize(_B)

    if IsEmpty(_A) or IsEmpty(_B):
        return []
    if coll_size_A!=row_size_B:
        raise ValueError("cusom_math::Prod: Product detected wrong size in the specified input array!")

    result = [ [ sum(_A[k][i]*_B[j][k] for k in range(coll_size_A)) for i in range(row_size_A) ] for j in range(coll_size_B) ]

    if CollSize(result) == 1:
        return result[0]
    else:
        return result

# ------------------------------------------------------------------------------
def ElemwiseProd(_X,_Y):
    if len(_X)!=len(_Y):
        raise ValueError("custom_math::ElemwiseProd: Wrong size of input vectors!")
    return [_X[i]*_Y[i] for i in range(len(_X))]

# ------------------------------------------------------------------------------
def Trans(_A):
    _A = SafeConvertVectorToMatrix(_A)

    if IsEmpty(_A):
        return []
    return list(map(list, zip(*_A)))
# ------------------------------------------------------------------------------
def Minus(_X,_Y):
    return [ _X[i]-_Y[i] for i in range(len(_X))]

# ------------------------------------------------------------------------------
def TranslateToNewBasis(_A, basis):
    _A = SafeConvertVectorToMatrix(_A)

    if IsEmpty(_A):
        return []

    trans_basis = Trans(basis)
    return Prod(trans_basis,_A)

# ------------------------------------------------------------------------------
def TranslateToOriginalBasis(_A, basis):
    _A = SafeConvertVectorToMatrix(_A)

    if IsEmpty(_A):
        return []

    return Prod(basis,_A)

# ------------------------------------------------------------------------------
def SolveLinearSystem(A,b):
    A = Trans(A)
    n = len(A)
    for i in range(n):
        A[i].append(b[i])
    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x

# ------------------------------------------------------------------------------
# Interior point algorithm
# Nocedal & Wright, Numerical Optimization, 2nd edition, Algorithm 16.4, p.484, Springer
# Solves min x'*x with A*x<=b (x is of size n and A of size m*n)
# Note that original algorithm solves min x'*x with A*x>=b --> negation of input values in the beginning
def QuadProg(A, b, max_itr, tolerance):
    m = RowSize(A)
    n = CollSize(A)

    A = ScalarMatrixProduct(-1.0,A)
    b = ScalarVectorProduct(-1.0,b)

    def rdrp(x,y,l,A,b):
        temp_vec = Prod(Trans(A),l)
        rd = Minus(x,temp_vec)

        temp_vec = Prod(A,x)
        rp = Minus(Minus(temp_vec,y),b)
        return (rd,rp)

    def GradResidu(y,l,A):
        m = RowSize(A)
        n = CollSize(A)
        grad = [ZeroVector(n+2*m) for i in range(n+2*m)]
        for i in range(n):
            grad[i][i] = 1
        for i in range(n):
            for j in range(m):
                grad[n+m+j][i] = -A[i][j]
        for i in range(m):
            for j in range(n):
                grad[j][n+i] = A[j][i]
        for i in range(m):
            grad[n+i][n+i] = -1
        for i in range(m):
            grad[n+i][n+m+i] = l[i]
            grad[n+m+i][n+m+i] = y[i]
        return grad

    #init
    x = ZeroVector(n)
    y = Ones(m) # slack variables
    l = Ones(m) # lagrange multipliers

    rd, rp = rdrp(x,y,l,A,b)
    k = 0
    error = 100000
    while error>tolerance and k <= max_itr:

        # solve affine delta
        gradRes = GradResidu(y,l,A)
        deltaXYLAff = ZeroVector(n+2*m)
        rhs = ScalarVectorProduct(-1, VertCat(rd,VertCat(rp,ElemwiseProd(y,l))))

        deltaXYLAff = SolveLinearSystem(gradRes,rhs)

        # Drop and error if divergence is detected
        if Norm2(deltaXYLAff)>1e20:
            KM.Logger.PrintWarning("ShapeOpt::custom_math::quadprodg", "deltaXYLAff is NAN. The reason is, that feasible domain might be empty. This happens e.g. when the dJdX is parallel to dCdX (like at convergence with a single constraint)")
            exit_code = 2
            return

        dxAff = deltaXYLAff[:n]
        dyAff = deltaXYLAff[n:n+m]
        dlAff = deltaXYLAff[n+m:]

        mu = Dot(y,l)/m
        alphaS = 1
        alphaZ = 1
        for i in range(m):
            if dyAff[i]<0:
                alphaS = min(alphaS, abs(y[i]/dyAff[i]))
            if dlAff[i]<0:
                alphaZ = min(alphaZ, abs(l[i]/dlAff[i]))
        alphaAff = 0.8*min(alphaS,alphaZ)

        muAff = Dot(Plus(y,ScalarVectorProduct(alphaAff,dyAff)),Plus(l,ScalarVectorProduct(alphaAff,dlAff))) / m

        sigma = (muAff/mu)**3

        # solve for delta
        rhs = ScalarVectorProduct(-1, VertCat(rd,VertCat(rp, Plus( Plus(ElemwiseProd(y,l),ElemwiseProd(dyAff,dlAff)) , ScalarVectorProduct(-sigma*mu,Ones(m)) ))))
        deltaXYL = SolveLinearSystem(gradRes,rhs)

        dx = deltaXYL[:n]
        dy = deltaXYL[n:n+m]
        dl = deltaXYL[n+m:]

        alphaS = 1
        alphaZ = 1
        for i in range(m):
            if dy[i]<0:
                alphaS = min(alphaS, abs(y[i]/dy[i]))
            if dl[i]<0:
                alphaZ = min(alphaZ, abs(l[i]/dl[i]))
        alpha = 0.8*min(alphaS,alphaZ)

        # apply delta
        x = Plus(x,ScalarVectorProduct(alpha,dx))
        y = Plus(y,ScalarVectorProduct(alpha,dy))
        l = Plus(l,ScalarVectorProduct(alpha,dl))
        rd,rp = rdrp(x,y,l,A,b)
        k = k+1

        # Determine current error
        error = Norm2(VertCat(rd,rp))

    if k==max_itr:
        raise RuntimeError("custom_math::quadprodg: Suboptimization for projection to halfspaces did not converge in the specified number of iterations")
    else:
        # Exit code of 0 if convergence is reached
        exit_code = 0

    return x, k, error, exit_code

# ------------------------------------------------------------------------------
def PerformBisectioning(func, a, b, target, tolerance, max_itr):

    fa_value_0, fa_is_converged = func(a)
    fb_value_0, fb_is_converged = func(b)

    error = None
    itr = 0

    # Always safe last argument for which function converges and return this if no solution is found
    if fa_is_converged:
        last_allowed_function_argument = a
    elif fb_is_converged:
        last_allowed_function_argument = b

    # Check special cases
    if abs(fa_value_0-target) < tolerance:
        res_function_argument = a
        error = abs(fa_value_0-target)

    elif abs(fb_value_0-target) < tolerance:
        res_function_argument = b
        error = abs(fb_value_0-target)

    elif abs(fa_value_0-fb_value_0) < 1e-13:
        KM.Logger.PrintWarning("ShapeOpt::PerformBisectioning", "Bisectioning intervall yiels to identical function values!")
        res_function_argument = a
        error = abs(fb_value_0-target)

    elif (fa_value_0-target)*(fb_value_0-target)>0:
        KM.Logger.PrintWarning("ShapeOpt::PerformBisectioning", "Bisectioning on function, that has no root in specified intervall!! Returning the argument which yiels a closer value to the target.")

        if abs(fa_value_0-target) < abs(fb_value_0-target):
            res_function_argument = a
            error = fa_value_0-target
        else:
            res_function_argument = b
            error = fb_value_0-target

    # Perform bisectioning if no special case applies
    else:
        p = (a + b)/2

        for itr in range(1,max_itr+1):
            fa_value, fa_is_converged = func(a)
            fp_value, fp_is_converged = func(p)

            if (fa_value-target)*(fp_value-target) < 0:
                b = p
            else:
                a = p

            p = (a + b)/2

            fp_value, fp_is_converged = func(p)
            error = abs(fp_value-target)

            if fp_is_converged:
                last_allowed_function_argument = p

            if itr == max_itr:
                res_function_argument = last_allowed_function_argument
                KM.Logger.PrintWarning("ShapeOpt::PerformBisectioning", "Bisectioning did not converge in the specified maximum number of iterations!")
                break

            elif error < tolerance:
                res_function_argument = p
                break

    return res_function_argument, itr, error

# --------------------------------------------------------------------------
def PerformGramSchmidtOrthogonalization(vector_space):
    V = vector_space
    B = []

    # Orthogonalization
    norm2_V0 = Norm2(V[0])
    B.append( ScalarVectorProduct(1/norm2_V0,V[0]) )
    for v in V[1:]:
        for b in B:
            norm2_b = Norm2(b)
            v = Minus( v , ScalarVectorProduct( Dot(v,b)/norm2_b**2 , b ) )

        # Add only if vector is independent
        norm2_v = Norm2(v)
        if norm2_v>1e-10:
            B.append( ScalarVectorProduct(1/norm2_v,v) )
        else:
            KM.Logger.PrintWarning("ShapeOpt::PerformGramSchmidtOrthogonalization", "Zero basis vector after Gram-Schmidt orthogonalization!")
            B.append(v)

    return B

# ==============================================================================
