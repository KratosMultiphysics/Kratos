"""
 @file solver.py

 PaStiX python object oriented wrapper

 @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.0.1
 @author Pierre Ramet
 @author Mathieu Faverge
 @author Louis Poirel
 @date 2018-07-16

"""
from .pastix import *
from .enum import *
import spm
import scipy.linalg as la

class solver(object):

    def __init__(self, A=None, **kwargs):
        self.iparm, self.dparm = initParam()

        # Verbose level
        self.verbose = kwargs.setdefault("verbose", verbose.No)
        if not self.verbose: # 0 or False
            self.iparm[iparm.verbose] = verbose.Not
        elif self.verbose==2:
            self.iparm[iparm.verbose] = verbose.Yes
        else: # 1 or True or anything else, default value
            self.iparm[iparm.verbose] = verbose.No

        # Set default factotype based on Matrix properties
        self.mtxtype = kwargs.setdefault("mtxtype", mtxtype.General)
        if self.mtxtype == mtxtype.HerPosDef:
            self.factotype = factotype.LLH
        elif self.mtxtype == mtxtype.SymPosDef:
            self.factotype = factotype.LLT
        elif self.mtxtype == mtxtype.Hermitian:
            self.factotype = factotype.LDLH
        elif self.mtxtype == mtxtype.Symmetric:
            self.factotype = factotype.LDLT
        else:
            self.factotype = factotype.LU

        # Factorization type
        self.factotype = kwargs.setdefault("factotype", self.factotype)
        self.iparm[iparm.factorization] = self.factotype

        # Threads
        self.thread_nbr = kwargs.setdefault("thread_nbr", "auto")
        if self.thread_nbr == "auto":
            self.thread_nbr = -1;
        self.iparm[iparm.thread_nbr] = self.thread_nbr

        self.pastix_data = init( self.iparm, self.dparm )
        if A is not None:
            self.setup(A)

    def setup(self, A):
        """
        Setup the solver object for the classic case w/o Schur complement
        """
        self.A = A
        self.spmA = spm.spmatrix(A)
        if self.verbose:
            self.spmA.printInfo()
        task_analyze(self.pastix_data, self.spmA)
        task_numfact(self.pastix_data, self.spmA)

    def solve(self, b, refine=True):
        """
        Solve and refine the full problem Ax = b with refinement when no Schur Complement is involved
        """
        x = b.copy()
        task_solve(self.pastix_data, self.spmA, x)
        if refine:
            task_refine(self.pastix_data, self.spmA, b, x)
        return x

    def schur(self, A, schur_list, full_matrix=True):
        """
        Setup the solver object to work in Schur complement mode

        full_matrix: boolean, default True. If full_matrix is False
        and a symmetric factorization is used, only the lower part of
        the Schur matrix is stored in self.S
        """

        self.A    = A
        self.spmA = spm.spmatrix(A)
        if self.verbose:
            self.spmA.printInfo()

        self.schur_list = schur_list + self.spmA.findBase()
        setSchurUnknownList(self.pastix_data, self.schur_list)

        task_analyze(self.pastix_data, self.spmA)
        task_numfact(self.pastix_data, self.spmA)

        nschur = len(schur_list)
        self.S = np.zeros( (nschur, nschur), order='F', dtype=A.dtype )
        getSchur( self.pastix_data, self.S )
        if full_matrix and (self.factotype != factotype.LU):
            self.S += la.tril(self.S, -1).T

    def schur_forward(self, b, solv_mode=solv_mode.Interface):
        """
        Solves the forward step of the Schur problem:
               A x = b with A = L S L^h

        This step solves  L f = b with f = S L^h x
        """
        self.iparm[iparm.schur_solv_mode] = solv_mode

        x = b.copy()

        # 1- Apply P to b
        subtask_applyorder( self.pastix_data, dir.Forward, x )

        # 2- Forward solve on the non Schur complement part of the system
        if self.factotype in (factotype.LLT, factotype.LLH):
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.NoTrans,
                          diag.NonUnit, x )
        else:
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.NoTrans,
                          diag.Unit, x )

        # 3- Apply the diagonal step if LDL factorization applied
        if self.factotype in (factotype.LDLT, factotype.LDLH):
            subtask_diag( self.pastix_data, x )

        self.x = x
        nschur = len(self.schur_list)
        return x[-nschur:]

    def schur_backward(self, y, b, solv_mode=solv_mode.Interface, refine=True):
        """
        Solves the backward step of the Schur problem:
               A x = b with A = L S L^h

        This step solves L^h x = y with y the solution of L S y = b
        """
        self.iparm[iparm.schur_solv_mode] = solv_mode

        nschur = len(self.schur_list)
        x = self.x
        x[-nschur:] = y
        # 4- Backward solve on the non Schur complement part of the system
        if self.factotype == factotype.LDLH:
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.ConjTrans,
                          diag.Unit, x )
        elif self.factotype == factotype.LDLT:
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.Trans,
                          diag.Unit, x )
        elif self.factotype == factotype.LLH:
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.ConjTrans,
                          diag.NonUnit, x )
        elif self.factotype == factotype.LLT:
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Lower, trans.Trans,
                          diag.NonUnit, x )
        else: # LU
            subtask_trsm( self.pastix_data, side.Left,
                          uplo.Upper, trans.NoTrans,
                          diag.NonUnit, x )

        #  5- Apply P^t to x
        subtask_applyorder( self.pastix_data, dir.Backward, x )

        return x

    def check( self, x, b, x0=None ):
        return self.spmA.checkAxb(x0, b, x)

    def finalize(self):
        finalize(self.pastix_data, self.iparm, self.dparm)
