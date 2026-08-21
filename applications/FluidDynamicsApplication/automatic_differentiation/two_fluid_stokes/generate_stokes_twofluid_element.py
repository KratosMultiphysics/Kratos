import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

do_simplifications = False
dim = 3 #spatial dimensions
mode = "c" #to output to a c++ file

if(dim == 2):
    nnodes = 3
    strain_size = 3
else:
    nnodes = 4
    strain_size = 6

initial_tabs = 0
max_index=30
optimizations='basic'
replace_indices=False

impose_partion_of_unity = False
N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

#defining the unknowns
v = DefineMatrix('v',nnodes,dim) #v(i,j) is velocity of node i component j
vn = DefineMatrix('vn',nnodes,dim) #velocity one step back
vnn = DefineMatrix('vnn',nnodes,dim) #velocity two step back
p = DefineVector('p',nnodes)

#define test functions
w = DefineMatrix('w',nnodes,dim)
q = DefineVector('q',nnodes)

#define other nodally varying data
f = DefineMatrix('f',nnodes,dim)

#constitutive matrix
C = DefineSymmetricMatrix('C',strain_size,strain_size)

#define other symbols
dt = sympy.Symbol('dt', positive=True)
rho = sympy.Symbol('rho', positive=True)
tau1 = sympy.Symbol('tau1', positive=True)
tau2 = sympy.Symbol('tau2')

bdf0 = sympy.Symbol('bdf0')
bdf1 = sympy.Symbol('bdf1')
bdf2 = sympy.Symbol('bdf2')

#interpolaing to the gauss point
fgauss = DefineVector('fgauss',dim)
vgauss = DefineVector('vgauss',dim)
acch = DefineVector('acch',dim)
pgauss = sympy.Symbol('pgauss')

wgauss = w.transpose()*N
qgauss = q.transpose()*N

#computing gradients
grad_v = DN.transpose()*v
grad_q = DN.transpose()*q

grad_p = DefineVector('grad_p',dim)
if(dim == 2):
    stress = DefineVector('stress',3)
elif(dim == 3):
    stress = DefineVector('stress',6)

#v_DV = vgauss.transpose()*grad_v
B = MatrixB(DN)
grad_sym_f = grad_sym_voigtform(DN,f)

div_v = div(DN,v)
div_w = div(DN,w)
grad_sym_w = grad_sym_voigtform(DN,w)

a = bdf0*v+bdf1*vn+bdf2*vnn
grad_sym_res = grad_sym_voigtform(DN,f-rho*a)


#compute galerkin functional
rv_galerkin =  wgauss.transpose()*(fgauss - rho*acch) +  div_w*pgauss - grad_sym_w.transpose()*stress  -qgauss*div_v

rv_stab = rho*div_w*tau2*div_v +  grad_q.transpose()*(rho*tau1*(fgauss - rho*acch - grad_p) )
rv_stab += -tau1*grad_sym_w.transpose()*C*grad_sym_res ##TODO: it might be better to remove this term - it is the symmetric_gradient of the subscale!

rv = rv_galerkin + rv_stab

#define dofs & test function vector
dofs = sympy.zeros(nnodes*(dim+1), 1)
testfunc = sympy.zeros(nnodes*(dim+1), 1)
for i in range(0,nnodes):
    for k in range(0,dim):
        dofs[i*(dim+1)+k] = v[i,k]
        testfunc[i*(dim+1)+k] = w[i,k]
    dofs[i*(dim+1)+dim] = p[i,0]
    testfunc[i*(dim+1)+dim] = q[i,0]
print("dofs = ",dofs)

rhs = Compute_RHS(rv, testfunc, do_simplifications)


##HERE WE MUST SUBSTITUTE EXPRESSIONS
strain = grad_sym_voigtform(DN,v)
rhs_to_derive = rhs.copy()

#vector values
rhs_to_derive = SubstituteMatrixValue( rhs_to_derive, fgauss, f.transpose()*N )
rhs_to_derive = SubstituteMatrixValue( rhs_to_derive, vgauss, v.transpose()*N )
rhs_to_derive = SubstituteMatrixValue( rhs_to_derive, acch, (bdf0*v+bdf1*vn+bdf2*vnn).transpose()*N )
rhs_to_derive = SubstituteMatrixValue( rhs_to_derive, stress, C*strain )
rhs_to_derive = SubstituteMatrixValue( rhs_to_derive, grad_p, DN.transpose()*p )

#scalar values
pgauss_expanded = (p.transpose()*N)[0]
rhs_to_derive = SubstituteScalarValue( rhs_to_derive, pgauss,pgauss_expanded )

##obtain LHS by symbolic derivation
lhs = Compute_LHS(rhs_to_derive, testfunc, dofs, do_simplifications)

#####################################################################
#####################################################################
## define enrichment variables
penr = DefineVector('penr',nnodes)
qenr = DefineVector('qenr',nnodes)
Nenr = DefineVector('Nenr',nnodes)
DNenr = DefineMatrix('DNenr',nnodes, dim)

grad_qenr = DNenr.transpose()*qenr
grad_penr = DNenr.transpose()*penr
penr_gauss = (penr.transpose()*Nenr)[0]

rv_enriched = div_w*penr_gauss + grad_qenr.transpose()*(rho*tau1*(fgauss - rho*acch - grad_p - grad_penr)) + grad_q.transpose()*(rho*tau1*(- grad_penr))
rv_enriched = SubstituteMatrixValue( rv_enriched, fgauss, f.transpose()*N )
rv_enriched = SubstituteMatrixValue( rv_enriched, vgauss, v.transpose()*N )
rv_enriched = SubstituteMatrixValue( rv_enriched, acch, (bdf0*v+bdf1*vn+bdf2*vnn).transpose()*N )
rv_enriched = SubstituteMatrixValue( rv_enriched, stress, C*strain )
rv_enriched = SubstituteMatrixValue( rv_enriched, grad_p, DN.transpose()*p )

dofs_enr = sympy.Matrix( [[penr[0,0]],[penr[1,0]],[penr[2,0]],[penr[3,0]]] )
testfunc_enr = sympy.Matrix( [[qenr[0,0]],[qenr[1,0]],[qenr[2,0]],[qenr[3,0]]] )

##  K V   x    =  b + rhs_eV
##  H Kee penr =  rhs_ee
rhs_eV, V   = Compute_RHS_and_LHS(rv_enriched, testfunc, dofs_enr, do_simplifications)
rhs_ee, H   = Compute_RHS_and_LHS(rv_enriched, testfunc_enr, dofs, do_simplifications)
rhs_ee, Kee = Compute_RHS_and_LHS(rv_enriched, testfunc_enr, dofs_enr, do_simplifications)

#####################################################################
#####################################################################
#simplified way of computing tau
if(dim == 2):
    mu_eq = C[2,2]
elif(dim == 3):
    mu_eq = (C[3,3]+C[4,4]+C[5,5])/3

inv_h2 = 0
for i in range(nnodes):
    for j in range(dim):
        inv_h2 += DN[i,j]**2

tau_denom = 3*mu_eq*inv_h2
tau_denom_out = OutputSymbolicVariable(tau_denom,mode)
#####################################################################
#####################################################################


templatefile = open("stokes_twofluid_cpp_template_3D.cpp")
outstring=templatefile.read()

outstring = outstring.replace("//substitute_lhs",   OutputMatrix_CollectingFactors(lhs,"lhs",mode,initial_tabs, max_index, optimizations,replace_indices) )
outstring = outstring.replace("//substitute_rhs", OutputVector_CollectingFactors(rhs,"rhs",mode,initial_tabs, max_index, optimizations,replace_indices))
outstring = outstring.replace("replace_tau_denom", tau_denom_out)


#####################################################################
#####################################################################
outstring = outstring.replace("//substitute_enrichment_V",   OutputMatrix_CollectingFactors(V,"V",mode,initial_tabs, max_index, optimizations,replace_indices) )
outstring = outstring.replace("//substitute_enrichment_H",   OutputMatrix_CollectingFactors(H,"H",mode,initial_tabs, max_index, optimizations,replace_indices))
outstring = outstring.replace("//substitute_enrichment_Kee", OutputMatrix_CollectingFactors(Kee,"Kee",mode,initial_tabs, max_index, optimizations,replace_indices))
outstring = outstring.replace("//substitute_enrichment_rhs_ee", OutputVector_CollectingFactors(rhs_ee,"rhs_ee",mode,initial_tabs, max_index, optimizations,replace_indices))
#####################################################################
#####################################################################

out = open("stokes_3D_twofluid.cpp",'w')
out.write(outstring)

out.close()

