import re

# Compile required regular expressions
# (see http://www.amk.ca/python/howto/regex/regex.html for an introduction to
# regular expressions in python)

# propertycode is used to detect "*loop properties" loops
propertycode=re.compile(r'^<FOR EACH PROPERTY>.*?<END>',re.DOTALL | re.MULTILINE)

# This identifies property type definitions for materials
prop_def=re.compile(r'<BEGIN\s*([\w]+)\s*>\s*(.*?)\s*<END\s*\1\s*>',re.DOTALL|re.MULTILINE)

# The followning regular expressions identify input for the various kinds of material properties
p_scalar=re.compile(r'(-*[.\d]+)')
p_2vector=re.compile(r'\((-*[.\d]+),(-*[.\d]+)\)')
p_3vector=re.compile(r'\((-*[.\d]+),(-*[.\d]+),(-*[.\d]+)\)')
p_2matrix=re.compile(r'\(\((-*[.\d]+),(-*[.\d]+)\),\((-*[.\d]+),(-*[.\d]+)\)\)')
p_3matrix=re.compile(r'\(\((-*[.\d]+),(-*[.\d]+),(-*[.\d]+)\),\((-*[.\d]+),(-*[.\d]+),(-*[.\d]+)\),\((-*[.\d]+),(-*[.\d]+),(-*[.\d]+)\)\)')
p_text=re.compile(r'([\w]+)')

# marked_text is used to detect leftover code between < > in code that will be cleaned
marked_text=re.compile(r'<.*?>',re.DOTALL | re.MULTILINE)

### Read Entities

# Extract an entity list from an input line
# Expected input:
# - inputlist: A list of words read from an input line
#   it is expected that all entity names will be at the beginning of the list
#   (so all previous words have to be read and removed before calling this function)
# - numline: The number of the line currently being read. It will be used in error messages
# It returns:
# - entities: The new list of entities
# - remlist: any remaining words in inputlist when the first one that is not an entity name is found
# - stop: will be True if an error occurred, telling the main function to stop reading the input file

lastentities='' # A global variable used by read_entities

def read_entities(inputlist,numline,empty_is_error=True):
    global lastentities
    
    entities=list()
    remline=list()
    stop=False
    while len(inputlist)>0:
        word=inputlist[0]
        if word=='point' or word=='line' or word=='surface' or word=='volume':
            entities.append(inputlist.pop(0))
        elif word=='all':
            entities=['point','line','surface','volume']
            inputlist.pop(0)
        elif word=='"':
            if len(lastentities)>0:
                entities=lastentities
                inputlist.pop(0)
            else:
                print 'ERROR: \'"\' used to specify entities in line',numline,'but there aren\'t any previous conditions or model parts to use as a reference'
                stop=True
            remline=inputlist
            break
        else:
            # It is assumed that following words are not related
            remline=inputlist
            break

    if len(entities)==0:
        if empty_is_error==True:
            print 'No entities found in line',numline
            stop=True
    else:
        entities.sort(key=entity_key)
        # Expected order: point,line,surface,volume
        lastentities=entities
        # Update last_entities only if we found new entities            
        
    return entities,remline,stop

### Convert list varibles to tuples
def list_to_tuple(lst):
    tup=()
    for item in lst:
        tup=tup+(item,)
    return tup

def tuple_to_list(tpl):
    lst=list()
    for item in tpl:
        lst.append(item)
    return lst

### Key used to sort a list of entity names by size (point, line,
### surface, volume)
def entity_key(ent_name):
    if ent_name=='point':
        return 0
    elif ent_name=='line':
        return 1
    elif ent_name=='surface':
        return 2
    elif ent_name=='volume':
        return 3

### code class

# The "code" class is used to store code that will be written in any problemtype
# file. four strings are required to create a code instance:
# -the code itself
# -two strings (preffix and suffix) to identify the destination file
# -a string used as a bookmark to know where in the destination file the code
#   will be written
# Code objects have two important methods:
# self.replace returns a new code oject based on the original, in which certain
# words have been replaced. This is useful, for example, to read generic input
# for all vector conditions and later customise this code for displacement or
# velocity
# self.write writes the code in the file (prefix)projecname(suffix)

class code:
    
    def __init__(self,code,dest_file_pref,dest_file_suff,insert_point):
        
        self.code=code
        # GiD or Tcl code (a string)
        
        self.file_prefix=dest_file_pref
        self.file_suffix=dest_file_suff
        # customised code will be written to a file named
        # (prefix)projecname(suffix) Both variables must be strings
        
        self.insertion_point=insert_point
        # A string which matches a certain line of the destination file, used to
        # know where the code must be written

    def filename(self,projectname):
        # Given the name of the project, returns the name of the file where
        # this code piece will be written
        return self.file_prefix+projectname+self.file_suffix

    def replace(self,oldnames,newnames):
        # Create a new code object where all instances of strings in oldnames
        # found in the code are replaced by strings from newnames
        custom_code=self.code
        for oldname,newname in zip(oldnames,newnames):
            custom_code=custom_code.replace(oldname,newname)

        newcode=code(custom_code,self.file_prefix,self.file_suffix,self.insertion_point)
        
        return newcode

    def write(self,projectname):
        """Writes the code to an existing file"""
        filename=self.filename(projectname)
        f=open(filename,'r')
        content=f.read()
        f.close()

        found=content.find(self.insertion_point)
        if found==-1:
            print 'WARNING: String',self.insertion_point,'not found in',filename
            print 'Some code could not be written\n'
            print self.code
        
        f=open(filename,'w')
        if self.insertion_point!='':
            newcontent=content.replace(self.insertion_point,self.code+self.insertion_point)
        else: # If no insertion point was specified, add new code to the end
            newcontent=content+self.code
        f.write(newcontent)
        f.close()

### code_entry class

# Stores all code items related to an object (a condition, an element,...)

class code_entry:

    code_pieces=()
    uses=()
    missing_input=()

    def __init__(self,name,entities): # string,tuple,string 
        self.name=name
        self.entities=entities

    def add_code(self,code,use,additional_input=()): # code,string
        self.code_pieces=self.code_pieces+(code,)
        self.uses=self.uses+(use,)
        self.missing_input=self.missing_input+(additional_input,)

    def find_code(self,target_use): # string
        found=()
        for code,use in zip(self.code_pieces,self.uses):
            if use==target_use:
                found=found+(code,)
        return found

    def find_index(self,target_use): # string
        """returns the index of the first code item with 'use'"""
        i=0
        for use in self.uses:
            if use==target_use:
                return i
            i=i+1

    def find_all(self,target_use): # string
        """returns the index of all code items with 'use'"""
        i=0
        found=()
        for use in self.uses:
            if use==target_use:
                found=found+(i,)
            i=i+1
        return found
        
    def delete_code(self,use): # string
        index=self.find_index(use)
        while index!=None:
            code_pieces=code_pieces[:index]+code_pieces[index+1:]
            uses=uses[:index]+uses[index+1:]
            index=self.find_index(use)

    def replace_code(self,index,newcode,newuse,newextra=()):
        self.code_pieces=self.code_pieces[:index]+(newcode,)+self.code_pieces[index+1:]
        self.uses=self.uses[:index]+(newuse,)+self.uses[index+1:]
        self.missing_input=self.missing_input[:index]+(newextra,)+self.missing_input[index+1:]

### code_container class

# Stores and organizes code entries

class code_container:

    def __init__(self):
        self.entry_list=list()

    def create_entry(self,name,entities):
        """Creates a new entry (will store all code for a specific object)"""
        self.entry_list.append(code_entry(name,entities))

    def find_index(self,name):
        """Returns the index of an entry"""
        i=0
        for entry in self.entry_list:
            if entry.name==name:
                return i
            i=i+1

    def swap_entry(self,index,new_entry):
        """Replaces an entire entry"""
        self.entry_list[index]=new_entry
        
    def add_code(self,entry,code,use,additional_input=()):
        """Adds a new code object to an entry"""
        index=self.find_index(entry)
        self.entry_list[index].add_code(code,use,additional_input)

    def delete_code(self,entry,use):
        """Deletes all code with 'use' from 'entry'"""
        entry_num=self.find_index(name)
        entry_list[entry_num].delete_code(use)

    def find(self,name,use):
        # name!='all', use!='all': returns all code pieces from an entry with a certain use
        # name=='all', use!='all': returns all code pieces from all entries with a certain use
        # name!='all', use=='all': reurns entry named 'name'
        
        if name!='all':
            index=self.find_index(name)
            if use !='all':
                return self.entry_list[index].find_code(use)
            else:
                return self.entry_list[index]
        else:
            # Assumes use!='all' (name=all, use=all -> the whole code_container instance)
            found=()
            for entry in self.entry_list:
                found=found+entry.find_code(use)
            return found

    def clean(self):
        """Removes text between < > in code objects with use 'wip' (work in progress). Changes use to 'write'"""
        for entry in self.entry_list:
            index=entry.find_index('wip')
            while index!=None:
                newcode=entry.code_pieces[index]
                newcode.code=marked_text.sub('',newcode.code)
                entry.replace_code(index,newcode,'write')
                index=entry.find_index('wip')

    def write(self,projectname):
        """Writes all finished code (identified by use=='write') to its file"""
        code_list=self.find('all','write')
        for code in code_list:
            code.write(projectname)

### base class

# This is a reference class to store all conditions, materials, properties, etc
# that will be read from input files. It stores code objects and writes them to
# relevant files.

class base:

    call=None
    definition_file=''
    insert_in=''

    basic_input=()
    additional_input=()

    def parseinput(self,extra_input):
        # The parseinput method is used by the problemtype generator to
        # understand additional input parameters specific to a condition or
        # element type. Obviously, this method should be defined for all derived
        # classes that require extra input parameters, but a simple parseinput
        # method is provided here that returns a tuple containing the input in
        # the same order it is read
        return list_to_tuple(extra_input)

    def valuestring(self,args):
        # Another simple method definition (see comments on parseinput), this
        # one is required to generate Tcl code
        first=True
        result=''
        for item in args:
            if first==True:
                result=str(item)
                first=False
            else:
                result=result+' '+str(item)
        return result

    def bas_entity_code(self,entities,name,mode,modifier=None):
        # Writes *set cond and *add cond for condition loops in .bas files,
        # based on the entities over which a particular condition can be applied
        # mode can be 'node' for properties applied over nodes or 'elem' for
        # properties applied over elements

        # The 'modifier' parameter must be a string, which will be added to all
        # *set and *add lines. This can be used to include instructions like
        # *CanRepeat, *Or and *And to those lines

        setadd='Set'
        bas_header=''

        if mode=='node':
            mesh_type='*nodes'
        elif mode=='elem':
            mesh_type='*elems'

        entity_list=tuple_to_list(entities)
        entity_list.sort(key=entity_key,reverse=True)
        # We want the *set and *add clauses in the following order:
        # volume, surface, line, point
        # Othertwise, conditions assigned to 'bigger' entities will override
        # those assigned to 'smaller' ones
        for entity in entity_list:
            bas_header=bas_header+\
                        '*'+setadd+' cond '+entity+'_'+name+' '+mesh_type
            if modifier!=None:
                bas_header=bas_header+' '+modifier+'\n'
            else:
                bas_header=bas_header+'\n'
            setadd='Add'

        return bas_header

### condition class

# Defines a generic GiD condition definition compatible with this problemtype
# generator's tcl functions

# To create a useful condition, define a new class using this one as a base,
# and define its parseinput() and valuestring() methods. You'll also need
# to define its questions, call, definition_file and additional_input attributes

class condition(base):

    basic_input=('<NAME>','<ENTITY>')
    insert_in='# Add Conditions Here'
    #call=None
    definition_file=''
    questions=''
    apply_over_nodes=True # Set this to False and you have a condition applied over elements instead of over nodes

    def __init__(self):

        self.temp_names=self.basic_input+self.additional_input

        if self.apply_over_nodes:
            condmeshtype='over nodes'
        else:
            condmeshtype='over body elements'

        cndtext='CONDITION: <ENTITY>_<NAME>\n'+\
                 'CONDTYPE: over <ENTITY>s\n'+\
                 'CONDMESHTYPE: '+condmeshtype+'\n'+\
                 'QUESTION: Automatic\n'+\
                 'VALUE: 0\n'+\
                 'STATE: HIDDEN\n'+\
                 self.questions+\
                 'END CONDITION\n\n'

        self.cndcode=code(cndtext,'','.cnd',self.insert_in)

    def add(self,line,numline):
        """Add a new boundary condition to the project. Generates the code for the .cnd file and for any template file"""
        from problemtype import code_db

        # Interpreting the line read from input file
        self.name=line.pop(0)
        self.entities,line,stop=read_entities(line,numline)
        if stop==True:
            return None # Return to the main loop
        values=self.parseinput(line)

        if self.apply_over_nodes:
            mode='node'
        else:
            mode='elem'

        bascode=self.bas_entity_code(self.entities,self.name,mode) # generate *Set and *Add clauses

        # Generate an entry in code_db for the new condition
        code_db.create_entry(self.name,self.entities)
        
        # generate .cnd code
        for entity in self.entities:
            new_names=(self.name,entity)+values
            new_cndcode=self.cndcode.replace(self.temp_names,new_names)
            code_db.add_code(self.name,new_cndcode,'write')

        # retrieve .bas code read from input files
        template_code=code_db.find(self.call,'.bas')
        other_code=code_db.find(self.call,'other')

        # generate .bas code
        for template in template_code:
            temp=template.replace(self.temp_names,new_names)
            temp.code=bascode+temp.code # add *Set and *Add clauses to the template
            code_db.add_code(self.name,temp,'write')

        # generate code for other files (usually tcl code)
        for template in other_code:
            temp=template.replace(self.temp_names,new_names)
            code_db.add_code(self.name,temp,'write')

        # Store a modified self.questions in the code_container (can be used by model parts)
        qtns=self.questions.replace('<NAME>',self.name)
        q_code=code(qtns,'','.cond','NOWHERE')
        code_db.add_code(self.name,q_code,'questions',self.additional_input)

### element class

# A finite element template

# Note that GiD doesn't handle point elements assigned to point entities in the same way as the other
# element types, as they can't be meshed via the mesh menu and materials aren't transferred to points
# Also, point conditions can be looped over 'elems' in bas files. All of this can be solved using tcl,
# but the different syntax in bas files means that they need to be defined as a different custom class,
# with an additional QUESTION field, and reading its template code from a different definition file.
# See the new_kratos definition folder for an example implementation.

class element(base):
    
    basic_input=('<NAME>','<ENTITY>')
    insert_in='# Add Elements Here\n'
    #call=None
    definition_file=''
    questions=''
    apply_over_nodes=False # Set this to true in derived classes for Point Elements

    def __init__(self):

        self.temp_names=self.basic_input+self.additional_input

        if self.apply_over_nodes==True: # When GiD supports point elements, delete this and write 'over body elements' always
            condmeshtype='over nodes'
        else:
            condmeshtype='over body elements'

        cndtext='CONDITION: <ENTITY>_<NAME>\n'+\
                 'CONDTYPE: over <ENTITY>s\n'+\
                 'CONDMESHTYPE: '+condmeshtype+'\n'+\
                 'QUESTION: Automatic\n'+\
                 'VALUE: 0\n'+\
                 'STATE: HIDDEN\n'+\
                 self.questions+\
                 'END CONDITION\n\n'

        self.cndcode=code(cndtext,'','.cnd',self.insert_in)

    def add(self,line,numline):
        """Add a new element type to the project. Generates the code for the .cnd file and for any template file"""
        from problemtype import code_db

        # Interpreting the line read from input file
        self.name=line.pop(0)
        self.entities,line,stop=read_entities(line,numline,empty_is_error=False)
        if stop==True:
            return None # Return to the main loop
        
        if len(self.entities)>1: # Only one entity (point, line, surface, volume) is allowed per element
            print 'WARNING: More than one entity has been specified for',self.name
            print 'Only',self.entities[0],'will be used'
            self.entities=[self.entities[0]]
        elif len(self.entities)==0:
            # Assume that we are defining a special point element.
            # The user will know it because point elements must be built from a different template, so they'll have a different call argument
            self.entities=['point']
            
        values=self.parseinput(line)

        # Generate an entry in code_db for the new condition
        code_db.create_entry(self.name,self.entities)

##        if self.apply_over_nodes:
##            mode='node'
##        else:
##            mode='elem'
        if 'point' in self.entities: ## This if should be removed once GiD supports point elements. Then, always execute the 'else' clause
            # If we are creating a point element, write a call to the point meshing procedure
            tcl_string='\tcreate_point_elems '+self.name+'\n'
            tcl_code=code(tcl_string,'','.tcl','\t# After Mesh Generation')
            code_db.add_code(self.name,tcl_code,'write')
            # For elements to work properly in GID 9.1.1b (and older), point entities are assigned a nodal condition by the user and then recieve an
            # elemental condition when the point elements are created. Because of this, their *set command is different
            bascode='*set cond element_'+self.name+' elem\n'
        else:
            # For other entities, we will want to loop over elements
            bascode=self.bas_entity_code(self.entities,self.name,'elem') # generate *Set and *Add clauses
            

        # retrieve .bas and .tcl code read from input files
        template_code=code_db.find(self.call,'.bas')
        other_code=code_db.find(self.call,'other')

        new_names=(self.name,self.entities[0])+values
        
        # generate .cnd code
        new_cndcode=self.cndcode.replace(self.temp_names,new_names)
        code_db.add_code(self.name,new_cndcode,'write')

        # generate .bas code
        for template in template_code:
            temp=template.replace(self.temp_names,new_names)
            temp.code=bascode+temp.code # add *Set and *Add clauses to the template
            code_db.add_code(self.name,temp,'write')

        # generate code for other files (usually tcl code)
        for template in other_code:
            temp=template.replace(self.temp_names,new_names)
            code_db.add_code(self.name,temp,'write')

        # Store a modified self.questions in the code_container (can be used by model parts)
        qtns=self.questions.replace('<NAME>',self.name)
        q_code=code(qtns,'','.cond','NOWHERE')
        code_db.add_code(self.name,q_code,'questions',self.additional_input)

### material class

# Creates a generic GiD material definition.

class material(base):
    
    basic_input=('<TYPE>','<NAME>')
    basic_types_list=('SCALAR','2DVECTOR','3DVECTOR','2X2MATRIX','3X3MATRIX','TEXT')
    known_prop_types=()
    #call='MATERIAL'

    def __init__(self,line,numline):
        """Add a new type of materials to the project"""
        from problemtype import code_db

        # retrieve .bas code read from input files and organize code related to properties
        bascode=code_db.find(self.call,'.bas')
        templates=self.read_template_code(bascode)

        # Interpreting the line read from input file
        name=line.pop(0)
        self.propnames=()
        self.proptypes=()
        while len(line)>1:
            self.propnames=self.propnames+(line.pop(0),)
            if line[0] in self.basic_types_list:
                if line[0] in self.known_prop_types:
                    self.proptypes=self.proptypes+(line.pop(0),)
                else:
                    print 'ERROR: Property type',line[0],'is not defined in this material\'s definition file (',self.definition_file,')'
            else:
                print 'ERROR: Unknown material type in line',numline,':',line[0]

        if len(line)==1:
            print 'ERROR: Unexpected length of material definitions. Check that all properties have its type specified.'

        self.mattype=name
        code_db.create_entry(name,())

        # Generate .mat code
        mattext='MATERIAL: <NAME>\n'+\
                 'QUESTION: Type\n'+\
                 'VALUE: <TYPE>\n'+\
                 'STATE: HIDDEN\n'

        for property_name,property_type in zip(self.propnames,self.proptypes):
            if property_type=='SCALAR':
                mattext=mattext+\
                         'QUESTION: '+property_name+'\n'+\
                         'VALUE: <VALUE_'+property_name+'>\n'

                self.additional_input=self.additional_input+('<VALUE_'+property_name+'>',)
            elif property_type=='2DVECTOR':
                mattext=mattext+\
                         'QUESTION: '+property_name+'_X\n'+\
                         'VALUE: <VALUE_'+property_name+'_X>\n'+\
                         'QUESTION: '+property_name+'_Y\n'+\
                         'VALUE: <VALUE_'+property_name+'_Y>\n'

                self.additional_input=self.additional_input+('<VALUE_'+property_name+'_X>','<VALUE_'+property_name+'_Y>')
            elif property_type=='3DVECTOR':
                mattext=mattext+\
                         'QUESTION: '+property_name+'_X\n'+\
                         'VALUE: <VALUE_'+property_name+'_X>\n'+\
                         'QUESTION: '+property_name+'_Y\n'+\
                         'VALUE: <VALUE_'+property_name+'_Y>\n'+\
                         'QUESTION: '+property_name+'_Z\n'+\
                         'VALUE: <VALUE_'+property_name+'_Z>\n'

                self.additional_input=self.additional_input+('<VALUE_'+property_name+'_X>','<VALUE_'+property_name+'_Y>','<VALUE_'+property_name+'_Z>')
            elif property_type=='2X2MATRIX':
                mattext=mattext+\
                         'QUESTION: '+property_name+'_XX\n'+\
                         'VALUE: <VALUE_'+property_name+'_XX>\n'+\
                         'QUESTION: '+property_name+'_XY\n'+\
                         'VALUE: <VALUE_'+property_name+'_XY>\n'+\
                         'QUESTION: '+property_name+'_YX\n'+\
                         'VALUE: <VALUE_'+property_name+'_YX>\n'+\
                         'QUESTION: '+property_name+'_YY\n'+\
                         'VALUE: <VALUE_'+property_name+'_YY>\n'

                self.additional_input=self.additional_input+\
                                       ('<VALUE_'+property_name+'_XX>','<VALUE_'+property_name+'_XY>')+\
                                       ('<VALUE_'+property_name+'_YX>','<VALUE_'+property_name+'_YY>')
            elif property_type=='3X3MATRIX':
                mattext=mattext+\
                         'QUESTION: '+property_name+'_XX\n'+\
                         'VALUE: <VALUE_'+property_name+'_XX>\n'+\
                         'QUESTION: '+property_name+'_XY\n'+\
                         'VALUE: <VALUE_'+property_name+'_XY>\n'+\
                         'QUESTION: '+property_name+'_XZ\n'+\
                         'VALUE: <VALUE_'+property_name+'_XZ>\n'+\
                         'QUESTION: '+property_name+'_YX\n'+\
                         'VALUE: <VALUE_'+property_name+'_YX>\n'+\
                         'QUESTION: '+property_name+'_YY\n'+\
                         'VALUE: <VALUE_'+property_name+'_YY>\n'+\
                         'QUESTION: '+property_name+'_YZ\n'+\
                         'VALUE: <VALUE_'+property_name+'_YZ>\n'+\
                         'QUESTION: '+property_name+'_ZX\n'+\
                         'VALUE: <VALUE_'+property_name+'_ZX>\n'+\
                         'QUESTION: '+property_name+'_ZY\n'+\
                         'VALUE: <VALUE_'+property_name+'_ZY>\n'+\
                         'QUESTION: '+property_name+'_ZZ\n'+\
                         'VALUE: <VALUE_'+property_name+'_ZZ>\n'

                self.additional_input=self.additional_input+\
                                       ('<VALUE_'+property_name+'_XX>','<VALUE_'+property_name+'_XY>','<VALUE_'+property_name+'_XZ>')+\
                                       ('<VALUE_'+property_name+'_YX>','<VALUE_'+property_name+'_YY>','<VALUE_'+property_name+'_YZ>')+\
                                       ('<VALUE_'+property_name+'_ZX>','<VALUE_'+property_name+'_ZY>','<VALUE_'+property_name+'_ZZ>')
            elif property_type=='TEXT':
                mattext=mattext+\
                         'QUESTION: '+property_name+'\n'+\
                         'VALUE: <VALUE_'+property_name+'>\n'

                self.additional_input=self.additional_input+('<VALUE_'+property_name+'>',)

        mattext=mattext+'END MATERIAL\n\n'

        self.temp_names=self.basic_input+self.additional_input

        # Create a book for this material family
        matbook='BOOK: '+name+' Materials\n'+\
                 '# Add '+name+' materials here\n'

        code_db.add_code(name,code(matbook,'','.mat',''),'write')

        self.matcode=code(mattext,'','.mat','# Add '+name+' materials here\n')

        bascode=self.write_property_code(templates,bascode,self.propnames,self.proptypes)

        for code_item in bascode:
            code_item=code_item.replace(self.temp_names,(name,'')+self.propnames)
            code_db.add_code(name,code_item,'write')

    def add(self,line,numline):
        """Adds a new material to the project"""

        from problemtype import code_db

        # Interpreting the line read from input file
        line.pop(0) # First word should be 'MATERIAL', can be ignored
        name=line.pop(0)
        prop_values=list_to_tuple(line)

        code_db.create_entry(name,())
        newvalues=(self.mattype,name)

        for value,proptype in zip(prop_values,self.proptypes):
            # Verify that the input belongs to the right property type
            if proptype=='SCALAR':
                m=p_scalar.match(value)
            elif proptype=='2DVECTOR':
                m=p_2vector.match(value)
            elif proptype=='3DVECTOR':
                m=p_3vector.match(value)
            elif proptype=='2X2MATRIX':
                m=p_2matrix.match(value)
            elif proptype=='3X3MATRIX':
                m=p_3matrix.match(value)
            elif proptype=='TEXT':
                m=p_text.match(value)
            if m==None:
                print 'ERROR: Some values for material properties don\'t match the property type:'
                print value,proptype
            else:
                newvalues=newvalues+m.groups()
        
        newcode=self.matcode.replace(self.temp_names,newvalues)

        code_db.add_code(name,newcode,'write')

    def read_template_code(self,bascode):
        """Searches property template definitions (<BEGIN PROPTYPE>code<END PROPTYPE>) in <FOR EACH PROPERTY> loops)"""
        templates=()

        for code_item in bascode:
            prop_blocks=propertycode.findall(code_item.code)

            for entry in prop_blocks:
                
                definitions=dict()
                matches=prop_def.findall(entry)

                for match in matches:
                    if match[0] in self.basic_types_list:
                        definitions[match[0]]=match[1]                        
                        if match[0] not in self.known_prop_types:
                            self.known_prop_types=self.known_prop_types+(match[0],)
                    else:
                        print 'WARNING: Unknown property type found in \'',self.definition_file,'\' definition file'
                        print 'Some code will not be used'
                templates=templates+(definitions,)
            
        return templates

    def write_property_code(self,templates,bascode,property_names,property_types):
        """Replace <FOR EACH PROPERTY> blocks with property code"""

        code_list=()
        finished_code=()
        for definitions in templates:
            newcode=''
            for propname,proptype in zip(property_names,property_types):
                try:
                    propcode=definitions[proptype]
                    newcode=newcode+propcode.replace('<PROPNAME>',propname)+'\n'
                except KeyError:
                    pass
            code_list=code_list+(newcode.rstrip('\n'),) # Remove trailing newline

        i=0
        
        for code_item in bascode:
            newcode=''
            prop_blocks=propertycode.findall(code_item.code)
            for ref in prop_blocks:
                newstring=propertycode.sub(code_list[i],code_item.code,1)
                newcode=code(newstring,code_item.file_prefix,code_item.file_suffix,code_item.insertion_point)
                i=i+1
            finished_code=finished_code+(newcode,)

        return finished_code

### Property class

# A simple single-field propety to use as GiD General Data

class gendata(base):
        
    basic_input=('<NAME>',)#'<VALUE>')
    insert_in='# Properties\n'
    #call='PROPERTY'
    definition_file=''
    questions='QUESTION: <NAME>\n'+\
               'VALUE: <VALUE>\n'

    def __init__(self):

        self.temp_names=self.basic_input+self.additional_input
        
        self.prbcode=code(self.questions,'','.prb',self.insert_in)
        
    def add(self,line,numline):
        """Adds a new item to General Data"""
        from problemtype import code_db

        # Interpreting the line read from input file
        name=line.pop(0)
##        value=line.pop(0)
        if len(line)>0:
            values=self.parseinput(line)
        else:
            values=()

        code_db.create_entry(name,())

##        newvalues=(name,value)+values
        newvalues=(name,)+values
        newcode=self.prbcode.replace(self.temp_names,newvalues)

        code_db.add_code(name,newcode,'write')

        # Read code from in definition files
        template_code=code_db.find(self.call,'.bas')+code_db.find(self.call,'other')

        # Generate code from definition file templates
        for template in template_code:
            temp=template.replace(self.temp_names,newvalues)
            code_db.add_code(name,temp,'write')


### Model Parts

# Describe the class

class part(base):

    basic_input=('<NAME>','<QUESTIONS>')
    insert_in='# Add Parts Here\n'
    call='MODEL PART'

    def basic_init(self,applytolowerentities):

        from problemtype import code_db

        self.conditions_always=[(),(),(),()] # Apply both in 2D and in 3D problems
        self.conditions_2Donly=[(),(),(),()] # Apply in 2D problems only
        self.conditions_3Donly=[(),(),(),()] # Apply in 3D problems only
        # Lists of 4 empty tuples. Each tuple stores the conditions this part applies to volumes, surfaces, lines and points, respectively

        # Generate tcl code
        code_db.create_entry(self.name,self.entities)

        # Tcl: list entities with this part assigned
        for entity in self.entities:
            list_tcl='\tset '+entity+self.name+'list [createlist '+entity+' '+self.name+']\n'
            list_code=code(list_tcl,'','.tcl','\t# End '+entity.capitalize()+' Block')
            code_db.add_code(self.name,list_code,'write')

        # Tcl: assign this condition to the boundary entities of a given entity
        if applytolowerentities==True:
            if 'volume' in self.entities:
                tcl_vol='\tcond_volumetosurface '+self.name+'\n'
                vol_code=code(tcl_vol,'','.tcl','\t# Inherited Surface Parts')
                code_db.add_code(self.name,vol_code,'write')
                
                tcl_surf='\tcond_surfacetoline '+self.name+'\n'
                surf_code=code(tcl_surf,'','.tcl','\t# Inherited Line Parts')
                code_db.add_code(self.name,surf_code,'write')

                tcl_line='\tcond_linetopoint '+self.name+'\n'
                line_code=code(tcl_line,'','.tcl','\t# Inherited Point Parts')
                code_db.add_code(self.name,line_code,'write')
                
            elif 'surface' in self.entities:
                tcl_surf='\tcond_surfacetoline '+self.name+'\n'
                surf_code=code(tcl_surf,'','.tcl','\t# Inherited Line Parts')
                code_db.add_code(self.name,surf_code,'write')
                
                tcl_line='\tcond_linetopoint '+self.name+'\n'
                line_code=code(tcl_line,'','.tcl','\t# Inherited Point Parts')
                code_db.add_code(self.name,line_code,'write')

            elif 'line' in self.entities:
                tcl_line='\tcond_linetopoint '+self.name+'\n'
                line_code=code(tcl_line,'','.tcl','\t# Inherited Point Parts')
                code_db.add_code(self.name,line_code,'write')

        # .cnd code
        for entity in self.entities:
            cndtext='CONDITION: '+entity+'_'+self.name+'\n'+\
                     'CONDTYPE: over '+entity+'s\n'+\
                     'CONDMESHTYPE: over nodes\n'+\
                     'QUESTION: Automatic\n'+\
                     'VALUE: 0\n'+\
                     'STATE: HIDDEN\n'+\
                     '<QUESTIONS>'+\
                     'END CONDITION\n\n'

            cndcode=code(cndtext,'','.cnd',self.insert_in)
            code_db.add_code(self.name,cndcode,'wip')

    def __init__(self,line,numline):

        # Interpreting the line read from input file
        #(name,entities,applytolowerentities)
        self.name=line.pop(0)
        self.entities,line,stop=read_entities(line,numline)
        if stop==True:
            return None # Return to the main loop
        if len(line)==0 or line==['LOWER']:
            applytolowerentities=True
        elif line==['NO','LOWER']:
            applytolowerentities=False
        else:
            print 'WARNING: Unexpected words in line',numline
            print '\'LOWER\' will be assumed for part',self.name
            applytolowerentities=True

        self.basic_init(applytolowerentities)

    def add(self,line,numline):
        """Associate a new condition or element to the part"""
        from problemtype import code_db, condition_dictionary

        # Interpreting the line read from input file (condition,class_obj,*defvalues)
        if line[0] in ('CONDITION','ELEMENT'):
            line.pop(0)
            mode='both'
        elif line[0]=='2D':
            line=line[2:]
            mode='2D'
        elif line[0]=='3D':
            line=line[2:]
            mode='3D'
        condition=line.pop(0)
        # The entities are optinal here, and can be used to specify that a condition will only be applied for some of the entities
        # where it's available
        ents,line,stop=read_entities(line,numline,empty_is_error=False)
        if len(ents)==0:
            ents=('point','line','surface','volume')
        defvalues=condition_dictionary[condition].parseinput(line)

        # Initiate an instance of the condition's class, to be able to access its valuestring method
        cond_class=condition_dictionary[condition]

        # We will work with this part's entry in code_container, replacing it with a new one
        index=code_db.find_index(self.name)
        new_entry=code_db.entry_list[index]

        # Retrieve the code entry related to the condition we want to add
        cond_index=code_db.find_index(condition)
        cond_entry=code_db.entry_list[cond_index]

        # Add the condition to the part's condition array
        for entity in ents:
            i=entity_key(entity)
            if (entity in self.entities) and (entity in cond_entry.entities):
                if mode=='both':
                    self.conditions_always[i]=self.conditions_always[i]+(condition,)
                elif mode=='2D':
                    self.conditions_2Donly[i]=self.conditions_2Donly[i]+(condition,)
                elif mode=='3D':
                    self.conditions_3Donly[i]=self.conditions_3Donly[i]+(condition,)

        k=cond_entry.find_index('questions')
        questions=cond_entry.code_pieces[k].code
        extra_input=cond_entry.missing_input[k]
        
        for old_word,new_word in zip(extra_input,defvalues):
            questions=questions.replace(old_word,new_word)

        numbers=new_entry.find_all('wip')

        if mode=='both':
            modestring=''
        elif mode=='2D':
            modestring='2D'
        else: #elif mode=='3D':
            modestring='3D'
        
        for j in numbers:
            
            stored_code=new_entry.code_pieces[j]

            if stored_code.code.find('<QUESTIONS>')!=-1:
                # This is .cnd code

                newcode='QUESTION: BEGIN'+modestring+condition+'\n'+\
                         'VALUE: BEGIN'+modestring+condition+'\n'+\
                         'STATE: HIDDEN\n'+\
                         questions+\
                         'QUESTION: END'+modestring+condition+'\n'+\
                         'VALUE: END'+modestring+condition+'\n'+\
                         'STATE: HIDDEN\n'+\
                         '<QUESTIONS>'

                cndcode=stored_code.replace(['<QUESTIONS>'],[newcode])
                new_entry.replace_code(j,cndcode,'wip')
                
            elif stored_code.code.find('<CONDITION DATA>')!=-1:
                # This is .tcl code

                tclcode=' BEGIN'+modestring+condition+' '+cond_class.valuestring(line)+' END'+modestring+condition+'<CONDITION DATA>'

                stored_code=stored_code.replace(['<CONDITION DATA>'],[tclcode])
                new_entry.replace_code(j,stored_code,'wip')

        code_db.swap_entry(index,new_entry)

### Default parts

# The following classes are variants of 'part' that can be applied automatically (using Tcl)
# to certain entities of the model.

# all_boundary is applied to the boundary of the model. The normal_dir parameter can be used to
# align the boundary entities' normals: They can be made to point inwards (normal_dir='inwards') or
# outwards (normal_dir='outwards')

class all_boundary(part):

    insert_in='# Add Default Conditions Here\n'
    call='BOUNDARY PART'

    def __init__(self,line,numline):

        from problemtype import code_db

        # Interpreting the line read from input file (name,normal_dir=None)
        self.name=line.pop(0)
        if line == ['INWARDS','NORMALS']:
            normal_dir='inwards'
        elif line == ['OUTWARDS','NORMALS']:
            normal_dir='outwards'
        else:
            normal_dir='NONE'
        
        self.entities=['point','line','surface']

        self.basic_init(True)

        # Align boundary normals
        if normal_dir=='inwards':
            surf_tcl='\t\talignsurfnormals Inwards\n'
            line_tcl='\t\talignlinenormals Inwards\n'
        elif normal_dir=='outwards':
            surf_tcl='\t\talignsurfnormals Outwards\n'
            line_tcl='\t\talignlinenormals Outwards\n'
        else:
            surf_tcl=''
            line_tcl=''

        # Begin Tcl code to assign this condition
        surf_tcl=surf_tcl+'\t\tGiD_AssignData Condition surface_'+self.name+' surfaces "1<CONDITION DATA>" $boundary\n'#+\
                  #'\t\tcond_surfacetoline '+self.name+'\n'' We are giving the same order in surf_code WHERE IS BETTER??
        line_tcl=line_tcl+'\t\tGiD_AssignData Condition line_'+self.name+' lines "1<CONDITION DATA>" $boundary\n'

        # Store tcl code in code_container
        surf_code=code(surf_tcl,'','.tcl','\t\t# 3D Boundary Section')
        code_db.add_code(self.name,surf_code,'wip')

        line_code=code(line_tcl,'','.tcl','\t\t# 2D Boundary Section')
        code_db.add_code(self.name,line_code,'wip')

# default_boundary is a special Model Part automatically assigned to all boundary
# entities without one of the conditions in partlist assigned

class default_boundary(part):

    insert_in='# Add Default Conditions Here\n'
    call='DEFAULT BOUNDARY PART'

    def __init__(self,line,numline):

        from problemtype import code_db

        # Interpreting the line read from input file (name,partlist)
        self.name=line.pop(0)
        line.pop(0) # Should be 'NOT', can be ignored
        partlist=list_to_tuple(line)
        
        self.entities=['point','line','surface']

        self.basic_init(True)

        # Start writing the required Tcl code
        surfacetcl='\tif {[GiD_Info Geometry NumVolumes]>0} {\n'
        linetcl='\tif {[GiD_Info Geometry NumVolumes]==0} {\n'

        for part in partlist:
            surfacetcl=surfacetcl+'\t\tset '+part+'_list [createlist surface '+part+']\n'
            linetcl=linetcl+'\t\tset '+part+'_list [createlist line '+part+']\n'
        
        surfacetcl=surfacetcl+'\t\tforeach entity $boundary {\n'+\
                    '\t\t\tif {\n'
        linetcl=linetcl+'\t\tforeach entity $boundary {\n'+\
                 '\t\t\tif {\n'

        for i in range(len(partlist)):
            surfacetcl=surfacetcl+'\t\t\t\t[lsearch $'+partlist[i]+'_list $entity]==-1'
            linetcl=linetcl+'\t\t\t\t[lsearch $'+partlist[i]+'_list $entity]==-1'
            if i<(len(partlist)-1):
                surfacetcl=surfacetcl+' && \\\n'
                linetcl=linetcl+' && \\\n'

        surfacetcl=surfacetcl+'\n\t\t\t} then {\n'+\
                    '\t\t\t\tlappend surface'+self.name+'list $entity\n'+\
                    '\t\t\t}\n'+\
                    '\t\t}\n'+\
                    '\t\tGiD_AssignData Condition surface_'+self.name+' surfaces "1<CONDITION DATA>" ${surface'+self.name+'list}\n\t}\n'

        linetcl=linetcl+'\n\t\t\t} then {\n'+\
                 '\t\t\t\tlappend line'+self.name+'list $entity\n'+\
                 '\t\t\t}\n'+\
                 '\t\t}\n'+\
                 '\t\tGiD_AssignData Condition line_'+self.name+' lines "1<CONDITION DATA>" ${line'+self.name+'list}\n\t}\n'

        surf_code=code(surfacetcl,'','.tcl','\t# Edit Surface Parts')
        code_db.add_code(self.name,surf_code,'wip')

        line_code=code(linetcl,'','.tcl','\t# Edit Line Parts')
        code_db.add_code(self.name,line_code,'wip')

# all_entities parts are assigned to all entities of a given type

class all_entities(part):
    
    insert_in='# Add Default Conditions Here\n'
    call='ALL ENTITIES PART'

    def __init__(self,line,numline):

        from problemtype import code_db

        # Interpreting the line read from input file (name,entities)
        self.name=line.pop(0)
        self.entities,remline,stop=read_entities(line,numline)
        if stop==True:
            return None # Return to the main loop

        self.basic_init(False)

        # Tcl code to assign the part to all entities
        if 'volume' in self.entities:
            tcode='\tGiD_AssignData Condition volume_'+self.name+' volumes "1<CONDITION DATA>" ${volumelist}\n'
            insert='\t# Assign Volume Parts\n'
            tclcode=code(tcode,'','.tcl',insert)
            code_db.add_code(self.name,tclcode,'wip')
        if 'surface' in self.entities:
            tcode='\tGiD_AssignData Condition surface_'+self.name+' surfaces "1<CONDITION DATA>" ${surfacelist}\n'
            insert='\t# Edit Surface Parts\n'
            tclcode=code(tcode,'','.tcl',insert)
            code_db.add_code(self.name,tclcode,'wip')
        if 'line' in self.entities:
            tcode='\tGiD_AssignData Condition line_'+self.name+' lines "1<CONDITION DATA>" ${linelist}\n'
            insert='\t# Edit Line Parts\n'
            tclcode=code(tcode,'','.tcl',insert)
            code_db.add_code(self.name,tclcode,'wip')
        if 'point' in self.entities:
            tcode='\tGiD_AssignData Condition point_'+self.name+' points "1<CONDITION DATA>" ${pointlist}\n'
            insert='\t# Edit Point Parts\n'
            tclcode=code(tcode,'','.tcl',insert)
            code_db.add_code(self.name,tclcode,'wip')

# default_entities parts are assigned to all entities of selected types except those which belong to some other parts
# warning: other automatically assigned parts involved in the definitions must have been defined before to avoid problems in Tcl code

class default_entities(part):

    insert_in='# Add Default Conditions Here\n'
    call='DEFAULT ENTITIES PART'

    def __init__(self,line,numline):

        from problemtype import code_db

        # Interpreting the line read from input file (name,partlist)
        self.name=line.pop(0)
        self.entities,remline,stop=read_entities(line,numline)
        if stop==True:
            return None # Return to the main loop
        
        line.pop(0) # Should be 'NOT', can be ignored
        partlist=list_to_tuple(line)

        self.basic_init(False)

        for entity in self.entities:
            for part in partlist:
                tclcode='\tset '+part+'_list [createlist '+entity+' '+part+']\n'
            tclcode=tclcode+'\t\tforeach entity $'+entity+'list {\n'+\
                    '\t\t\tif {\n'

            for i in range(len(partlist)):
                tclcode=tclcode+'\t\t\t[lsearch $'+partlist[i]+'_list $entity]==-1'
                if i<(len(partlist)-1):
                    tclcode=tclcode+' && \\\n'
                    
                tclcode=tclcode+'\n\t\t\t} then {\n'+\
                         '\t\t\t\tlappend '+entity+self.name+'list $entity\n'+\
                         '\t\t\t}\n'+\
                         '\t\t}\n'+\
                         '\t\tGiD_AssignData Condition '+entity+'_'+self.name+' '+entity+'s "1<CONDITION DATA>" ${'+entity+self.name+'list}\n\t}\n'

        tcl_code=code(tclcode,'','.tcl','\t# Edit Surface Parts')
        code_db.add_code(self.name,tcl_code,'wip')

# An additional class: Option
# This class allows the user to choose an element while assigning or in Problem Data.
# While meshing, the element will be used in all entitites with the condition assigned
# Its intended use is as a complement to the model parts
# If an element is choosen while assigning the model part, the entity will be meshed with that element
# If an element is not coosen (or if the part is assigned automatically using Tcl),
# the entity will be meshed with the element specyfied in Problem Data

class option(condition):

    insert_in='# Add Default Conditions Here\n'
    call='OPTION'
    definition_file='' # This condition doesn't have a definition file, it is an internal device and it shouldn't appear in any template file

    questions='QUESTION: <NAME>#CB#(<OPTIONS>)\n'+\
               'VALUE: <DEFVALUE>\n'

    additional_input=('<OPTIONS>','<DEFVALUE>')

    def add(self,line,numline):
        """Add a new boundary condition to the project. Generates the code for the .cnd file and for any template file"""
        from problemtype import code_db,auto_conditions

        # Interpreting the line read from input file
        self.name=line.pop(0)
        self.entities,line,stop=read_entities(line,numline)
        if stop==True:
            return None # Return to the main loop

        if len(self.entities)>1: # Only one entity (point, line, surface, volume) is allowed per element
            print 'WARNING: More than one entity has been specified for',self.name
            print 'Only',self.entities[0],'will be used'
            entities=(self.entities[0],)
        
        values=self.internal_parseinput(line)
        # Note that the input is not exacty the same for parts that use it or for the condition itself
        # The menu in Problem Data doesn't have a 'Use default' option, as it sets the default value

        # Add the available elements to the auto_conditions dictionary, as they can be assigned automatically now
        # All GiD conditions in auto_conditions will be cleaned before meshing
        auto_elements=values[0].split(',')
        for elem in auto_elements:
            try:
                auto_conditions[elem].append(self.entities[0])
            except KeyError:
                auto_conditions[elem]=[self.entities[0]]

        # Generate an entry in code_db for the new condition
        code_db.create_entry(self.name,self.entities)

        bascode=self.bas_entity_code(self.entities,self.name,'node') # generate *Set and *Add clauses

        newvalues=(self.name,self.entities[0])+values

        # generate .cnd code
        new_cndcode=self.cndcode.replace(self.temp_names,newvalues)
        code_db.add_code(self.name,new_cndcode,'write')

        # Generate an entry in code_db for the new condition
        code_db.create_entry(self.name,self.entities)

        prb_code='QUESTION:<NAME>#CB#(<OPTIONS>)\n'+\
                  'VALUE: <DEFVALUE>\n'

        for old,new in zip(self.temp_names,newvalues):
            prb_code=prb_code.replace(old,new)

        prb_entry=code(prb_code,'','.prb','# Options')
        code_db.add_code(self.name,prb_entry,'write')

        # Generate Tcl Code to assign the proper element
        tcl_code='\tassign_element_choice '+self.name+' '+self.entities[0]+' '+values[0].replace(',',' ')+'\n'
        
        tcl_entry=code(tcl_code,'','.tcl','\t# Select Elements from Options')
        code_db.add_code(self.name,tcl_entry,'write')

        # Store a modified self.questions in the code_container (can be used by model parts)
        qtns=self.questions.replace('<NAME>',self.name)
        qtns=qtns.replace('<OPTIONS>',values[0])
        q_code=code(qtns,'','.cond','NOWHERE')
        code_db.add_code(self.name,q_code,'questions',self.additional_input)

    def internal_parseinput(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No available elements found for this option'
            return
        defvalue=str(extra_input[0])
        first=True
        for val in extra_input:
            if first==False:
                values=values+','+str(val)
            else:
                values=str(val)
                first=False
                
        return values,defvalue
    
    def parseinput(self,extra_input):
        if len(extra_input)==0:
            return '','Use_Default'
        else:
            if len(extra_input)>1:
                print 'More values than expected found, only the first one will be used'
                print '\t'+str(extra_input)

            return '',extra_input[0]

    def valuestring(self,extra_input):
        if len(extra_input)==0:
            extra_input='Use_Default'
        elif len(extra_input)>1:
            print 'More values than expected found, only the first one will be used'
            print '\t'+str(extra_input)
            extra_input=extra_input[0]
            
        return extra_input
