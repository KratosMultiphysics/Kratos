import os;

import basicfunctions;

#####################################################################################


def GenerateAllFiles(projectname,projectpath,mode):
	os.mkdir(projectpath+projectname+'.gid')
	if mode==1:
		ImportFiles(projectname,projectpath,1,'archiv/')
		os.chdir(projectpath+projectname+'.gid/')
		GenerateCndFile(projectname+'.cnd',1)
		GenerateMatFile(projectname+'.mat',1)
		GenerateFluidFiles(projectname,projectpath)
	elif mode==2:
		ImportFiles(projectname,projectpath,2,'archiv/')
		os.chdir(projectpath+projectname+'.gid/')
		GenerateCndFile(projectname+'.cnd',2)
		GenerateMatFile(projectname+'.mat',2)
		GenerateStructureFiles(projectname,projectpath)
	elif mode==3:
		ImportFiles(projectname,projectpath,3,'archiv/')
		os.chdir(projectpath+projectname+'.gid/')
		GenerateCndFile(projectname+'.cnd',3)
		GenerateMatFile(projectname+'.mat',3)
		GenerateFluidFiles(projectname,projectpath)
		GenerateStructureFiles(projectname,projectpath)
	else:
		print 'choose mode'
	
#################################################################################################

def GenerateFluidFiles(projectname,projectpath):
	GenerateFluidElemBas('011_'+projectname+'.fluid.elem.bas')
	GenerateFluidCondBas('012_'+projectname+'.fluid.cond.bas')
	GenerateFluidInitBas('013_'+projectname+'.fluid.init.bas')
	
def GenerateStructureFiles(projectname,projectpath):
	GenerateStructureElemBas('021_'+projectname+'.structure.elem.bas')
	GenerateStructureCondBas('022_'+projectname+'.structure.cond.bas')
	GenerateStructureInitBas('023_'+projectname+'.structure.init.bas')


###############################################################################################

################################################################################################

def ImportFiles(projectname,projectpath,mode,archivpath):
	quelle=archivpath+'coupling.prb'
	ziel=projectpath+projectname+'.gid'+'/'+projectname+'.prb'
	basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
	quelle=archivpath+'coupling.bas'
	ziel=projectpath+projectname+'.gid'+'/'+projectname+'.bas'
	basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
	if mode==1:
		quelle=archivpath+'010_coupling.fluid.prop.bas'
		ziel=projectpath+projectname+'.gid'+'/010_'+projectname+'.fluid.prop.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'014_coupling.node.bas'
		ziel=projectpath+projectname+'.gid'+'/014_'+projectname+'.fluid.node.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'fluid.win.bat'
		ziel=projectpath+projectname+'.gid'+'/'+projectname+'.win.bat'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))

		quelle=archivpath+'fluid.unix.bat'
		ziel=projectpath+projectname+'.gid'+'/'+projectname+'.unix.bat'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
	elif mode==2:
		quelle=archivpath+'020_coupling.structure.prop.bas'
		ziel=projectpath+projectname+'.gid'+'/020_'+projectname+'.structure.prop.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'024_coupling.node.bas'
		ziel=projectpath+projectname+'.gid'+'/024_'+projectname+'.structure.node.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'structure.win.bat'
		ziel=projectpath+projectname+'.gid'+'/'+projectname+'.win.bat'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
	elif mode==3:
		quelle=archivpath+'010_coupling.fluid.prop.bas'
		ziel=projectpath+projectname+'.gid'+'/010_'+projectname+'.fluid.prop.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'014_coupling.fluid.node.bas'
		ziel=projectpath+projectname+'.gid'+'/014_'+projectname+'.fluid.node.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'020_coupling.structure.prop.bas'
		ziel=projectpath+projectname+'.gid'+'/020_'+projectname+'.structure.prop.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'024_coupling.structure.node.bas'
		ziel=projectpath+projectname+'.gid'+'/024_'+projectname+'.structure.node.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'fluidstructure.win.bat'
		ziel=projectpath+projectname+'.gid'+'/'+projectname+'.win.bat'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
		
		quelle=archivpath+'030_coupling.coupling.bas'
		ziel=projectpath+projectname+'.gid'+'/030_'+projectname+'.coupling.bas'
		basicfunctions.writefile(ziel,basicfunctions.readfile(quelle))
	
###############################################################################################

################################################################################################

def GenerateCndFile(filename,mode):
	if mode==1:
		filecontent='BOOK: Fluid_Boundary_Conditions\n'
		filecontent=filecontent+'BOOK: Fluid_Element_Type\n'
	elif mode==2:
		filecontent='BOOK: Structure_Boundary_Conditions\n'
		filecontent=filecontent+'BOOK: Structure_Element_Type\n'
	elif mode==3:
		con1='BOOK: Fluid_Boundary_Conditions\n'
		con2='BOOK: Fluid_Element_Type\n'
		con3='BOOK: Fluid_Domain\n'+CndFluidDomainContent()
		con4='BOOK: Structure_Boundary_Conditions\n'
		con5='BOOK: Structure_Element_Type\n'
		con6='BOOK: Structure_Domain\n'+CndStructureDomainContent()
		filecontent=con1+con2+con3+con4+con5+con6
	basicfunctions.writefile(filename,filecontent)
	
def CndFluidDomainContent():
	fluid=DomainCond('Fluid','point')+DomainCond('Fluid','line')+DomainCond('Fluid','surface')+DomainCond('Fluid','volume')
	pfluidb=DomainCond('positive_Fluidboundary','point')+DomainCond('positive_Fluidboundary','line')+DomainCond('positive_Fluidboundary','surface')
	nfluidb=DomainCond('negative_Fluidboundary','point')+DomainCond('negative_Fluidboundary','line')+DomainCond('negative_Fluidboundary','surface')
	return fluid+pfluidb+nfluidb	

def CndStructureDomainContent():
	struc=DomainCond('Structure','point')+DomainCond('Structure','line')+DomainCond('Structure','surface')+DomainCond('Structure','volume')
	strucb=DomainCond('Structureboundary','point')+DomainCond('Structureboundary','line')+DomainCond('Structureboundary','surface')
	return struc+strucb	
	
def DomainCond(name,condtype):
	cond='CONDITION: '+condtype+'_'+name+'\n'
	condt='CONDTYPE: over '+condtype+'s \n'
	condmeshtype='CONDMESHTYPE: over nodes\n'
	question='QUESTION: Domain_Type#CB#('+name+')\nVALUE: '+name+'\n'
	end='END CONDITION\n'
	return cond+condt+condmeshtype+question+end

################################################################################################

def GenerateMatFile(filename,mode):
	if mode==1:
		filecontent='BOOK: Fluid_Material\n'
	elif mode==2:
		filecontent='BOOK: Structure_Material\n'
	elif mode==3:
		filecontent='BOOK: Fluid_Material\n'
		filecontent=filecontent+'BOOK: Structure_Material\n'
	basicfunctions.writefile(filename,filecontent)
	
#################################################################################################

def GenerateFluidElemBas(filename):
	filecontent='// Reading Elements\n\nElementsGroup = fluid_group;\n\n'
	basicfunctions.writefile(filename,filecontent)

#################################################################################################

def GenerateStructureElemBas(filename):
	filecontent='// Reading Elements\n\nElementsGroup = structural_group;\n\n'
	basicfunctions.writefile(filename,filecontent)

#################################################################################################

def GenerateFluidCondBas(filename):
	filecontent='\n\n// Fixing degrees of freedom in nodes\n\n'
	basicfunctions.writefile(filename,filecontent)

#################################################################################################

def GenerateStructureCondBas(filename):
	filecontent='\n\n// Fixing degrees of freedom in nodes\n\n'
	basicfunctions.writefile(filename,filecontent)

#################################################################################################

def GenerateFluidInitBas(filename):
	filecontent='// Fixing degrees of freedom in nodes\n\n'
	basicfunctions.writefile(filename,filecontent)

#################################################################################################

def GenerateStructureInitBas(filename):
	filecontent='// Fixing degrees of freedom in nodes\n\n'
	basicfunctions.writefile(filename,filecontent)

################################################################################################

##############################################################################################

