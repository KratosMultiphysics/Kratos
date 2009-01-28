import string;

import basicfunctions;

################################################################################################

def WriteElems(name,condtype,projectname,domaintype):
	section=''
	if (string.count(condtype,'p')>=1):
		section = section + ElemSection(name,'point')
	if (string.count(condtype,'l')>=1):
		section = section + ElemSection(name,'line')
	if (string.count(condtype,'s')>=1):
		section = section + ElemSection(name,'surface')
	if (string.count(condtype,'v')>=1):
		section = section + ElemSection(name,'volume')
	if domaintype=='Fluid':
		basicfunctions.addtofile('011_'+projectname+'.fluid.elem.bas',section)
	if domaintype=='Structure':
		basicfunctions.addtofile('021_'+projectname+'.structure.elem.bas',section)


#####################################################################################################

def ElemSection(name,condtype):
	section='*Set cond '+condtype+'_'+name+' *elems\n'
	section = section + '*loop elems *onlyInCond\n*Set var i=0\n*set var j= ElemsNnode\n*format "%i%i%i%i%i%i%i%i"\n'
	section = section + 'ELEMENTS[*ElemsNum] = '+name+'([*\\\n'
	section = section + '*for(i=1;i<j;i=i+1)*\\\n*ElemsConec(*i),*\\\n*end*\\\n*ElemsConec(*ElemsNnode)],*ElemsMat);\n*end elems\n'
	return section