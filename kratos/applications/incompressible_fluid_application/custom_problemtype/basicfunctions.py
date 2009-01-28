import re;

def getkeyword(keyword,filecontent):
	seperator = keyword+'(?m)'
	m = re.split(seperator,filecontent);
	del m[0]
	for x in range(len(m)):
		m[x]=keyword+m[x]
	return m

def splitfile(keyword,filecontent):
	seperator = keyword+'(?m)'
	m = re.split(seperator,filecontent);
	return m

def addtofilebeginning(filename,content):
	newcontent=content+readfile(filename)
	writefile(filename,newcontent)

def writefile(filename,filecontent):
	datei=open(filename,'w')
	datei.write(filecontent)
	datei.close

def readfile(filename):
	datei = open(filename,'r')
	filecontent = datei.read()
	datei.close()
	return filecontent

def addtofile(filename,filecontent):
	datei=open(filename,'a')
	datei.write(filecontent)
	datei.close