import os
import stat
import re

# Compile regular expresions
book=re.compile(r'(^BOOK:[\s]*(.*?)[\s]*$)',re.MULTILINE)
cnd_book=re.compile(r'^BOOK:[\s]*(.*?)[\s\n]*(?:CONDITION:|TITLE:)',re.MULTILINE)
mat_book=re.compile(r'^BOOK:[\s]*(.*?)[\s\n]*(?:MATERIAL:|TITLE:)',re.MULTILINE)
prb_book=re.compile(r'^BOOK:[\s]*(.*?)[\s\n]*(?:QUESTION:|TITLE:)',re.MULTILINE)
gendata_block=re.compile(r'^PROBLEM DATA[\s\n]*(.*?)[\s\n]*END GENERAL DATA',re.MULTILINE | re.DOTALL)
intdata_block=re.compile(r'^INTERVAL DATA[\s\n]*(.*)\Z',re.MULTILINE | re.DOTALL)

# Create required files

def generate_files(projectname,templates_path):
    """Creates all required files, either copying them from the 'files' folder
or creating empty ones"""
    os.mkdir(projectname+'.gid')
    files_folder = os.path.join(templates_path,'files')
    filelist=os.listdir(files_folder)
    for filename in filelist:
	# print filename
	if (filename.startswith('.')==False):
	    newname=filename.replace('problemtype',projectname)
	    copyfile(os.path.join(files_folder,filename),'./'+projectname+'.gid/'+newname)
    # Check that some relevant files have been created
    if 'problemtype.cnd' not in filelist:
        copyfile('./default/problemtype.cnd','./'+projectname+'.gid/'+projectname+'.cnd')
    if 'problemtype.mat' not in filelist:
        writefile('./'+projectname+'.gid/'+projectname+'.mat','')
    if 'problemtype.prb' not in filelist:
        copyfile('./default/problemtype.prb','./'+projectname+'.gid/'+projectname+'.prb')
    if 'problemtype.tcl' not in filelist:
        copyfile('./default/problemtype.tcl','./'+projectname+'.gid/'+projectname+'.tcl')
    if 'problemtype.bas' not in filelist:
        print 'No problemtype.bas file found in \'files\', an empty file was created\n'
        writefile('./'+projectname+'.gid/'+projectname+'.bas','')
    if 'problemtype.unix.bat' not in filelist:
        print 'WARNING: no .unix.bat file found in \'files\', no .unix.bat file was created\n'
    else:
        pth=os.path.join(os.getcwd(),projectname+'.gid',projectname+'.unix.bat')
        os.chmod(pth,stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
    if 'problemtype.win.bat' not in filelist:
        print 'WARNING: no .win.bat file found in \'files\', no .win.bat file was created\n'
    else:
        pth=os.path.join(os.getcwd(),projectname+'.gid',projectname+'.win.bat')
        os.chmod(pth,stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
    # Note that GiD can work with an empty problemtype.bas file (if there are
    # other .bas files in the same folder), but empty .bat files are useless
    print 'Created problemtype folder'
##    os.chdir(projectname+'.gid/')

# Add default python scripts to the problemtype. Python scripts will be copied
# from the folder 'script_folder', which must be found at the same path as the
# input file

def add_python_scripts(projectname,script_folder):
    problemtype_folder = os.path.join(os.getcwd(),projectname+'.gid')
    script_folder = os.path.join(os.getcwd(),script_folder)
    script_list = os.listdir(script_folder)
    for filename in script_list:
        extension = os.path.splitext(filename)[1]
        if extension == '.py':
            origin = os.path.join(script_folder,filename)
            destination = os.path.join(problemtype_folder,filename)
            copyfile(origin,destination)
            print 'Added Python script',filename

# Delete empty books and generate a custom menu for the problem type

def check_books(projectname,templates_path):
    menuname=projectname.replace('_',' ')
    if menuname[0].islower():
        menuname=menuname[0].capitalize()+menuname[1:]
    menuindex=0
    menucode='\tGiDMenu::Create "'+menuname+'" PRE\n'
    cndbooks=()
    matbooks=()
    prbbooks=()
    # First, find all condition books. Create a menu entry for each non-empty
    # one (except 'Default', which is for conditions that are assigned
    # automatically)
    file_content=readfile(projectname+'.cnd')
    matches=cnd_book.findall(file_content)
    for item in matches:
        if item!='Default':
            menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                      '" [list "'+item+'"] '+str(menuindex)+' PRE "GidOpenConditions \\"'+item+'\\"" "" ""\n'
            menuindex=menuindex+1
        cndbooks=cndbooks+(item,)
    # Rewrite the file, cleaning any unused book
    pieces=book.split(file_content)
    newcode=pieces[0]
    for j in range(0,len(pieces)-3,3):
        if pieces[j+2] in cndbooks:
            newcode=newcode+pieces[j+1]+pieces[j+3]
        else: # The book is unused, forget it
            newcode=newcode+pieces[j+3]
    writefile(projectname+'.cnd',newcode)
    
    # Now, add the material books to the menu
    file_content=readfile(projectname+'.mat')
    matches=mat_book.findall(file_content)
    for item in matches:
        menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                  '" [list "'+item+'"] '+str(menuindex)+' PRE "GidOpenMaterials \\"'+item+'\\"" "" ""\n'
        menuindex=menuindex+1
        matbooks=matbooks+(item,)
    # Rewrite the file, cleaning any unused book
    pieces=book.split(file_content)
    newcode=pieces[0]
    for j in range(0,len(pieces)-3,3):
        if pieces[j+2] in matbooks:
            newcode=newcode+pieces[j+1]+pieces[j+3]
        else: # The book is unused, forget it
            newcode=newcode+pieces[j+3]
    writefile(projectname+'.mat',newcode)
    
    # The property and interval data books
    file_content=readfile(projectname+'.prb')
    properties=gendata_block.findall(file_content)
    if len(properties)>0:
        matches=prb_book.findall(properties[0])
        if len(matches)==0:
            menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                      '" [list "Problem Data"] '+str(menuindex)+' PRE "GidOpenProblemData" "" ""\n'
            menuindex=menuindex+1
        else:
            for item in matches:
                menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                          '" [list "'+item+'"] '+str(menuindex)+' PRE "GidOpenProblemData \\"'+item+'\\"" "" ""\n'
                menuindex=menuindex+1
                prbbooks=prbbooks+(item,)
    intervals=intdata_block.findall(file_content)
    if len(intervals)>0:
        matches=prb_book.findall(intervals[0])
        if len(matches)==0:
            menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                      '" [list "Interval Data"] '+str(menuindex)+' PRE "GidOpenIntervals" "" ""\n'
            menuindex=menuindex+1
        else:
            for item in matches:
                menucode=menucode+'\tGiDMenu::InsertOption "'+menuname+\
                          '" [list "'+item+'"] '+str(menuindex)+' PRE "GidOpenIntervals \\"'+item+'\\"" "" ""\n'
                menuindex=menuindex+1
                prbbooks=prbbooks+(item,)
    # Rewrite the file, cleaning any unused book
    pieces=book.split(file_content)
    newcode=pieces[0]
    for j in range(0,len(pieces)-3,3):
        if pieces[j+2] in prbbooks:
            newcode=newcode+pieces[j+1]+pieces[j+3]
        else: # The book is unused, forget it
            newcode=newcode+pieces[j+3]
    writefile(projectname+'.prb',newcode)
    
    # Add a last function to the menu and finish the code
    menucode=menucode+\
              	'\tGiDMenu::InsertOption "'+menuname+'" [list "---"] '+str(menuindex)+' PRE "" "" ""\n'+\
              	'\tGiDMenu::InsertOption "'+menuname+'" [list "Model Status"] '+str(menuindex+1)+' PRE "cond_report" "" ""\n'+\
              	'\tGiDMenu::UpdateMenus\n'

    return menucode

# Auxiliary functions

def writefile(filename,filecontent):
    datei=open(filename,'w')
    datei.write(filecontent)
    datei.close()

def readfile(filename):
    datei = open(filename,'r')
    filecontent = datei.read()
    datei.close()
    return filecontent

def copyfile(filename1,filename2):
    content = readfile(filename1)
    writefile(filename2,content)
