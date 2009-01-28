import string;

import basicfunctions;

################################################################################################

def WriteScalarCond(name,condtype,meshtype,projectname,domaintype):
	newsection=''
	if (string.count(condtype,'p')>=1):
		newsection=newsection+ScalarCond(name,'point',meshtype,domaintype)
	if (string.count(condtype,'l')>=1):
		newsection=newsection+ScalarCond(name,'line',meshtype,domaintype)
	if (string.count(condtype,'s')>=1):
		newsection=newsection+ScalarCond(name,'surface',meshtype,domaintype)
	if (string.count(condtype,'v')>=1):
		newsection=newsection+ScalarCond(name,'volume',meshtype,domaintype)
	AddBCToCndFile(newsection,projectname,domaintype)
	
	

def WriteVectorCond(name,condtype,meshtype,projectname,domaintype):
	newsection=''
	if (string.count(condtype,'p')>=1):
		newsection=newsection+VectorCond(name,'point',meshtype,domaintype)
	if (string.count(condtype,'l')>=1):
		newsection=newsection+VectorCond(name,'line',meshtype,domaintype)
	if (string.count(condtype,'s')>=1):
		newsection=newsection+VectorCond(name,'surface',meshtype,domaintype)
	if (string.count(condtype,'v')>=1):
		newsection=newsection+VectorCond(name,'volume',meshtype,domaintype)
	AddBCToCndFile(newsection,projectname,domaintype)



def WriteElemCond(name,elemtype,projectname,domaintype):
	newsection=''
	if (string.count(elemtype,'p')>=1):
		newsection=newsection+ElemCond(name,'point')
	if (string.count(elemtype,'l')>=1):
		newsection=newsection+ElemCond(name,'line')
	if (string.count(elemtype,'s')>=1):
		newsection=newsection+ElemCond(name,'surface')
	if (string.count(elemtype,'v')>=1):
		newsection=newsection+ElemCond(name,'volume')
	AddElemToCndFile(newsection,projectname,domaintype)

	
##################################################################################

def ScalarCond(name,condtype,meshtype,domaintype):
	cond='CONDITION: '+condtype+'_'+name+'_('+domaintype+')\n'
	condt='CONDTYPE: over '+condtype+'s \n'
	condmeshtype='CONDMESHTYPE: over '+meshtype+'\n'
	question='QUESTION: Value\nVALUE: 0.0\n'
	end='END CONDITION\n'
	return cond+condt+condmeshtype+question+end

def VectorCond(name,condtype,meshtype,domaintype):
	cond='CONDITION: '+condtype+'_'+name+'_('+domaintype+')\n'
	cond=cond+'CONDTYPE: over '+condtype+'s \n'
	cond=cond+'CONDMESHTYPE: over '+meshtype+'\n'
	cond=cond+'QUESTION:'+name+'_X#CB#(0,1)\nVALUE: 0\nDEPENDENCIES: (0,HIDE,Value_X,#CURRENT#)(1,RESTORE,Value_X,#CURRENT#)\n'
	cond=cond+'QUESTION: Value_X\nVALUE: 0.0\n'
	cond=cond+'QUESTION:'+name+'_Y#CB#(0,1)\nVALUE: 0\nDEPENDENCIES: (0,HIDE,Value_Y,#CURRENT#)(1,RESTORE,Value_Y,#CURRENT#)\n'
	cond=cond+'QUESTION: Value_Y\nVALUE: 0.0\n'
	cond=cond+'QUESTION:'+name+'_Z#CB#(0,1)\nVALUE: 0\nDEPENDENCIES: (0,HIDE,Value_Z,#CURRENT#)(1,RESTORE,Value_Z,#CURRENT#)\n'
	cond=cond+'QUESTION: Value_Z\nVALUE: 0.0\n'
	cond=cond+'END CONDITION\n'
	return cond
	
def ElemCond(name,elemtype):
	cond='CONDITION: '+elemtype+'_'+name+'\n'
	condt='CONDTYPE: over '+elemtype+'s \n'
	condmeshtype='CONDMESHTYPE: over body elements\n'
	question='QUESTION: Element_Type#CB#('+name+')\nVALUE: '+name+'\n'
	end='END CONDITION\n'
	return cond+condt+condmeshtype+question+end


########################################################################################################

def AddBCToCndFile(newsection,projectname,domaintype):
	if domaintype=='Fluid':
		key='BOOK: Fluid_Boundary_Conditions\n'
	if domaintype=='Structure':
		key='BOOK: Structure_Boundary_Conditions\n'
	file=projectname+'.cnd'
	filecontent = basicfunctions.splitfile(key,basicfunctions.readfile(file))
	newcontent = filecontent[0] + key + newsection +filecontent[1]
	basicfunctions.writefile(file,newcontent)

def AddElemToCndFile(newsection,projectname,domaintype):
	if domaintype=='Fluid':
		key='BOOK: Fluid_Element_Type\n'
	if domaintype=='Structure':
		key='BOOK: Structure_Element_Type\n'
	file=projectname+'.cnd'
	filecontent = basicfunctions.splitfile(key,basicfunctions.readfile(file))
	newcontent = filecontent[0] + key + newsection +filecontent[1]
	basicfunctions.writefile(file,newcontent)





























#################################  ALT  ########################################
################################################################################def WriteFluidScalarCond(name,condtype,meshtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[0]
#	restcontent = filecontent[1]+filecontent[2]+filecontent[3]+filecontent[4]+filecontent[5]
#	if (string.count(condtype,'p')>=1):
#		section=section+ScalarCond(name,'point',meshtype,'Fluid')
#	if (string.count(condtype,'l')>=1):
#		section=section+ScalarCond(name,'line',meshtype,'Fluid')
#	if (string.count(condtype,'s')>=1):
#		section=section+ScalarCond(name,'surface',meshtype,'Fluid')
#	if (string.count(condtype,'v')>=1):
#		section=section+ScalarCond(name,'volume',meshtype,'Fluid')
#	filecontent=section+restcontent
#	basicfunctions.writefile(projectname+'.cnd',filecontent)


#def WriteFluidVectorCond(name,condtype,meshtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[0]
#	restcontent = filecontent[1]+filecontent[2]+filecontent[3]+filecontent[4]+filecontent[5]
#	if (string.count(condtype,'p')>=1):
#		section=section+VectorCond(name,'point',meshtype,'Fluid')
#	if (string.count(condtype,'l')>=1):
#		section=section+VectorCond(name,'line',meshtype,'Fluid')
#	if (string.count(condtype,'s')>=1):
#		section=section+VectorCond(name,'surface',meshtype,'Fluid')
#	if (string.count(condtype,'v')>=1):
#		section=section+VectorCond(name,'volume',meshtype,'Fluid')
#	filecontent=section+restcontent
#	basicfunctions.writefile(projectname+'.cnd',filecontent)
	

#def WriteStructureScalarCond(name,condtype,meshtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[3]
#	restcontent1 = filecontent[0]+filecontent[1]+filecontent[2]
#	restcontent2 = filecontent[4]+filecontent[5]
#	if (string.count(condtype,'p')>=1):
#		section=section+ScalarCond(name,'point',meshtype,'Structure')
#	if (string.count(condtype,'l')>=1):
#		section=section+ScalarCond(name,'line',meshtype,'Structure')
#	if (string.count(condtype,'s')>=1):
#		section=section+ScalarCond(name,'surface',meshtype,'Structure')
#	if (string.count(condtype,'v')>=1):
#		section=section+ScalarCond(name,'volume',meshtype,'Structure')
#	filecontent=restcontent1+section+restcontent2
#	basicfunctions.writefile(projectname+'.cnd',filecontent)


#def WriteStructureVectorCond(name,condtype,meshtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[3]
#	restcontent1 = filecontent[0]+filecontent[1]+filecontent[2]
#	restcontent2 = filecontent[4]+filecontent[5]
#	if (string.count(condtype,'p')>=1):
#		section=section+VectorCond(name,'point',meshtype,'Structure')
#	if (string.count(condtype,'l')>=1):
#		section=section+VectorCond(name,'line',meshtype,'Structure')
#	if (string.count(condtype,'s')>=1):
#		section=section+VectorCond(name,'surface',meshtype,'Structure')
#	if (string.count(condtype,'v')>=1):
#		section=section+VectorCond(name,'volume',meshtype,'Structure')
#	filecontent=restcontent1+section+restcontent2
#	basicfunctions.writefile(projectname+'.cnd',filecontent)

#def WriteFluidElemCond(name,elemtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[1]
#	restcontent1 = filecontent[0]
#	restcontent2 = filecontent[2]+filecontent[3]+filecontent[4]+filecontent[5]
#	if (string.count(elemtype,'p')>=1):
#		section=section+ElemCond(name,'point')
#	if (string.count(elemtype,'l')>=1):
#		section=section+ElemCond(name,'line')
#	if (string.count(elemtype,'s')>=1):
#		section=section+ElemCond(name,'surface')
#	if (string.count(elemtype,'v')>=1):
#		section=section+ElemCond(name,'volume')
#	filecontent=restcontent1+section+restcontent2
#	basicfunctions.writefile(projectname+'.cnd',filecontent)
#	
#	
#def WriteStructureElemCond(name,elemtype,projectname):
#	filecontent = basicfunctions.getkeyword('BOOK: ',basicfunctions.readfile(projectname+'.cnd'))
#	section = filecontent[4]
#	restcontent1 = filecontent[0]+filecontent[1]+filecontent[2]+filecontent[3]
#	restcontent2 = filecontent[5]
#	if (string.count(elemtype,'p')>=1):
#		section=section+ElemCond(name,'point')
#	if (string.count(elemtype,'l')>=1):
#		section=section+ElemCond(name,'line')
#	if (string.count(elemtype,'s')>=1):
#		section=section+ElemCond(name,'surface')
#	if (string.count(elemtype,'v')>=1):
#		section=section+ElemCond(name,'volume')
#	filecontent=restcontent1+section+restcontent2
#	basicfunctions.writefile(projectname+'.cnd',filecontent)
#