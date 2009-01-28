import string;

import basicfunctions;

################################################################################################

def WriteScalarBC(name,condtype,projectname,domaintype):
	section=''
	setadd='Set'
	if (string.count(condtype,'v')>=1):
		section = section + CondSetCond(name,'volume',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'s')>=1):
		section = section + CondSetCond(name,'surface',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'l')>=1):
		section = section + CondSetCond(name,'line',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'p')>=1):
		section = section + CondSetCond(name,'point',setadd,domaintype)
		setadd='Add'
	if setadd == 'Add':
		section = section + CondScalarSection(name)
	if domaintype=='Fluid':
		basicfunctions.addtofile('012_'+projectname+'.fluid.cond.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofile('022_'+projectname+'.structure.cond.bas',section)

def WriteVectorBC(name,condtype,projectname,domaintype):
	section=''
	setadd='Set'
	if (string.count(condtype,'v')>=1):
		section = section + CondSetCond(name,'volume',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'s')>=1):
		section = section + CondSetCond(name,'surface',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'l')>=1):
		section = section + CondSetCond(name,'line',setadd,domaintype)
		setadd='Add'
	if (string.count(condtype,'p')>=1):
		section = section + CondSetCond(name,'point',setadd,domaintype)
		setadd='Add'
	if setadd == 'Add':
		section = section + CondVectorSection(name)
	if domaintype=='Fluid':
		basicfunctions.addtofile('012_'+projectname+'.fluid.cond.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofile('022_'+projectname+'.structure.cond.bas',section)

def WriteConditions(name,condtype,projectname,domaintype):
	section=''
	if (string.count(condtype,'p')>=1):
		section = section + ConditionSection(name,'point')
	if (string.count(condtype,'l')>=1):
		section = section + ConditionSection(name,'line')
	if (string.count(condtype,'s')>=1):
		section = section + ConditionSection(name,'surface')
	if (string.count(condtype,'v')>=1):
		section = section + ConditionSection(name,'volume')
	if domaintype=='Fluid':
		basicfunctions.addtofilebeginning('012_'+projectname+'.fluid.cond.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofilebeginning('022_'+projectname+'.structure.cond.bas',section)

	
#####################################################################################################

def CondSetCond(name,condtype,setadd,domaintype):
	return '*'+setadd+' cond '+condtype+'_'+name+'_('+domaintype+')'+' *nodes\n'

def CondScalarSection(name):
	section ='*loop nodes  *OnlyInCond\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum].Fix('+name+');\n'
	section = section + '*end nodes\n'
	return section

def CondVectorSection(name):
	section ='*loop nodes  *OnlyInCond\n'
	section = section + '*if(strcmp(cond('+name+'_X),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum].Fix('+name+'_X);\n'
	section = section + '*endif\n*if(strcmp(cond('+name+'_Y),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum].Fix('+name+'_Y);\n'
	section = section + '*endif\n*if(strcmp(cond('+name+'_Z),"1")==0)\n*format "%i%f"\n'
	section = section + 'NODES[*NodesNum].Fix('+name+'_Z);\n'
	section = section + '*endif\n*end nodes\n'
	return section
	
	
def ConditionSection(name,condtype):
	section='*Set cond '+condtype+'_'+name+' *elems\n'
	section = section + '*loop elems *onlyInCond\n*Set var i=0\n*set var j= ElemsNnode\n*format "%i%i%i%i%i%i%i%i"\n'
	section = section + 'CONDITIONS[*ElemsNum] = '+name+'([*\\\n'
	section = section + '*for(i=1;i<j;i=i+1)*\\\n*ElemsConec(*i),*\\\n*end*\\\n*ElemsConec(*ElemsNnode)],*ElemsMat);\n*end elems\n'
	return section