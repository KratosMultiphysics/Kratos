import write_cnd;
import write_mat;
import write_elem_bas;
import write_cond_bas;
import write_init_bas;




	
def AddScalarBC(name,condtype,meshtype,projectname,domaintype):
	write_cnd.WriteScalarCond(name,condtype,meshtype,projectname,domaintype)	
	write_cond_bas.WriteScalarBC(name,condtype,projectname,domaintype)
	write_init_bas.WriteScalarBC(name,condtype,projectname,domaintype)

def AddVecBC(name,condtype,meshtype,projectname,domaintype):
	write_cnd.WriteVectorCond(name,condtype,meshtype,projectname,domaintype)	
	write_cond_bas.WriteVectorBC(name,condtype,projectname,domaintype)
	write_init_bas.WriteVectorBC(name,condtype,projectname,domaintype)

def AddCondition(name,condtype,projectname,domaintype):
	write_cnd.WriteElemCond(name,condtype,projectname,domaintype)
	write_cond_bas.WriteConditions(name,condtype,projectname,domaintype)
	
def AddElement(name,condtype,projectname,domaintype):
	write_cnd.WriteElemCond(name,condtype,projectname,domaintype)
	write_elem_bas.WriteElems(name,condtype,projectname,domaintype)

def AddFluidMaterial(name,density,viscosity,gravity,projectname):
	write_mat.WriteFluidMat(name,density,viscosity,gravity,projectname)

def AddStructureMaterial(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname):
	write_mat.WriteStructureMat(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname)

####################################################################

######################   ALT  ##########################################
#def AddFluidScalarBC(name,condtype,meshtype,projectname):
#	write_cnd.WriteFluidScalarCond(name,condtype,meshtype,projectname)	
#	write_cond_bas.WriteFluidScalarBC(name,condtype,projectname)
#	write_init_bas.WriteFluidScalarBC(name,condtype,projectname)
#
#def AddFluidVecBC(name,condtype,meshtype,projectname):
#	write_cnd.WriteFluidVectorCond(name,condtype,meshtype,projectname)	
#	write_cond_bas.WriteFluidVectorBC(name,condtype,projectname)
#	write_init_bas.WriteFluidVectorBC(name,condtype,projectname)
#
#def AddStructureScalarBC(name,condtype,meshtype,projectname):
#	write_cnd.WriteStructureScalarCond(name,condtype,meshtype,projectname)	
#	write_cond_bas.WriteStructureScalarBC(name,condtype,projectname)
#	write_init_bas.WriteStructureScalarBC(name,condtype,projectname)
#
#def AddStructureVecBC(name,condtype,meshtype,projectname):
#	write_cnd.WriteStructureVectorCond(name,condtype,meshtype,projectname)
#	write_cond_bas.WriteStructureVectorBC(name,condtype,projectname)
#	write_init_bas.WriteStructureVectorBC(name,condtype,projectname)
#
#
#def AddFluidElement(name,condtype,projectname):
#	write_cnd.WriteFluidElemCond(name,condtype,projectname)
#	write_elem_bas.WriteFluidElems(name,condtype,projectname)
#
#
#def AddStructureElement(name,condtype,projectname):
#	write_cnd.WriteStructureElemCond(name,condtype,projectname)
#	write_elem_bas.WriteStructureElems(name,condtype,projectname)

