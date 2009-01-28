import basicfunctions;

################################################################################################

def WriteFluidMat(name,density,viscosity,gravity,projectname):
	newsection=FluidMat(name,density,viscosity,gravity)
	key='BOOK: Fluid_Material\n'
	file=projectname+'.mat'
	filecontent = basicfunctions.splitfile(key,basicfunctions.readfile(file))
	newcontent = filecontent[0] + key + newsection +filecontent[1]
	basicfunctions.writefile(file,newcontent)

def WriteStructureMat(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname):
	newsection=StructureMat(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname)
	key='BOOK: Structure_Material\n'
	file=projectname+'.mat'
	filecontent = basicfunctions.splitfile(key,basicfunctions.readfile(file))
	newcontent = filecontent[0] + key + newsection +filecontent[1]
	basicfunctions.writefile(file,newcontent)


##################################################################################

def FluidMat(name,density,viscosity,gravity):
	section='MATERIAL: '+name+'\n'
	section=section+'QUESTION: ID#CB#(Fluid)\nVALUE: Fluid\nSTATE: HIDDEN\n'
	section=section+'QUESTION: Density\nVALUE: '+density+'\n'
	section=section+'QUESTION: Viscosity\nVALUE: '+viscosity+'\n'
	section=section+'QUESTION: gravity\nVALUE: '+gravity+'\nEND MATERIAL\n'
	return section

def StructureMat(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname):
	section='MATERIAL: '+name+'\n'
	section=section+'QUESTION: ID#CB#(Structure)\nVALUE: Structure\nSTATE: HIDDEN\n'
	section=section+'QUESTION: Density\nVALUE: '+density+'\n'
	section=section+'QUESTION: Young_Module\nVALUE: '+young_module+'\n'
	section=section+'QUESTION: Poisson_Ratio\nVALUE: '+poisson_ratio+'\n'
	section=section+'QUESTION: Thickness\nVALUE: '+thickness+'\n'
	section=section+'QUESTION: Cross_Section_Area\nVALUE: '+cross_section_area+'\nEND MATERIAL\n'
	return section

##################################################################################
#
### Mit Units ###
#
#def FluidMat(name,density,viscosity,gravity):
#	section='MATERIAL: '+name+'\n'
#	section=section+'QUESTION: ID#CB#(Fluid)\nVALUE: Fluid\nSTATE: HIDDEN\n'
#	section=section+'QUESTION: Density#UNITS#\nVALUE: '+density+'kg/m^3\n'
#	section=section+'QUESTION: Viscosity#UNITS#\nVALUE: '+viscosity+'kg/m*s\n'
#	section=section+'QUESTION: gravity#UNITS#\nVALUE: '+gravity+'m/s^2\nEND MATERIAL\n'
#	return section
#
#def StructureMat(name,density,young_module,poisson_ratio,thickness,cross_section_area,projectname):
#	section='MATERIAL: '+name+'\n'
#	section=section+'QUESTION: ID#CB#(Structure)\nVALUE: Structure\nSTATE: HIDDEN\n'
#	section=section+'QUESTION: Density#UNITS#\nVALUE: '+density+'kg/m^3\n'
#	section=section+'QUESTION: Young_Module#UNITS#\nVALUE: '+young_module+'N/mm^2\n'
#	section=section+'QUESTION: Poisson_Ratio\nVALUE: '+poisson_ratio+'\n'
#	section=section+'QUESTION: Thickness#UNITS#\nVALUE: '+thickness+'mm\n'
#	section=section+'QUESTION: Cross_Section_Area#UNITS#\nVALUE: '+cross_section_area+'mm^2\nEND MATERIAL\n'
#	return section
#
###########################################################################################