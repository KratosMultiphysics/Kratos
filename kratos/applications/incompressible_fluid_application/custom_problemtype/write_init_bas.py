import string;

import basicfunctions;

################################################################################################

def WriteScalarBC(name,condtype,projectname,domaintype):
	section=''
	setadd='Set'
	if (string.count(condtype,'v')>=1):
		section = section + InitSetCond(name,'volume',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'s')>=1):
		section = section + InitSetCond(name,'surface',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'l')>=1):
		section = section + InitSetCond(name,'line',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'p')>=1):
		section = section + InitSetCond(name,'point',setadd,domaintype)
		setadd='Add'
	if setadd == 'Add':
		section = section + InitScalarSection(name)
	if domaintype=='Fluid':
		basicfunctions.addtofile('013_'+projectname+'.fluid.init.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofile('023_'+projectname+'.structure.init.bas',section)

def WriteVectorBC(name,condtype,projectname,domaintype):
	section=''
	setadd='Set'
	if (string.count(condtype,'v')>=1):
		section = section + InitSetCond(name,'volume',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'s')>=1):
		section = section + InitSetCond(name,'surface',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'l')>=1):
		section = section + InitSetCond(name,'line',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'p')>=1):
		section = section + InitSetCond(name,'point',setadd,domaintype)
		setadd='Add'
	if setadd == 'Add':
		section = section + InitVectorSection(name)
	if domaintype=='Fluid':
		basicfunctions.addtofile('013_'+projectname+'.fluid.init.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofile('023_'+projectname+'.structure.init.bas',section)


#####################################################################################################

def InitSetCond(name,condtype,setadd,domaintype):
	return '*'+setadd+' cond '+condtype+'_'+name+'_('+domaintype+')'+' *nodes\n'

def InitScalarSection(name):
	section ='*loop nodes  *OnlyInCond\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum]('+name+',0) = *cond(Value);\n'
	section = section + '*end nodes\n'
	return section

def InitVectorSection(name):
	section ='*loop nodes  *OnlyInCond\n'
	section = section + '*if(strcmp(cond('+name+'_X),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum]('+name+'_X,0) = *cond(Value_X);\n'
	section = section + '*endif\n*if(strcmp(cond('+name+'_Y),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum]('+name+'_Y,0) = *cond(Value_Y);\n'
	section = section + '*endif\n*if(strcmp(cond('+name+'_Z),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum]('+name+'_Z,0) = *cond(Value_Z);\n'
	section = section + '*endif\n*end nodes\n'
	return section
