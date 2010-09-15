import read_tools
import file_functions
import core_definitions
import tcl_functions
import os

def generate(input_file,base_folder = '.'):

    global code_db
    global condition_dictionary
    global auto_conditions

    # Define some variables
    condition_dictionary=dict()
    property_dictionary=dict()
    material_dictionary=dict()
    part_dictionary=dict()
    current_matrial=None
    current_group=None
    currentconds=list()
    auto_conditions=dict()
    script_folder = None
    # tcl variables
    point_code=''
    point_end=''
    line_code=''
    line_end=''
    surf_code=''
    surf_end=''
    vol_code=''
    vol_end=''
    
    source=open(input_file,'r')

    # Read the first line of the input file
    line=source.readline()
    numline=1

    stop=False # If something goes wrong, stop will be set to True and this
    # function will end, diplaying an error message

    # Loop over the lines of the input file
    while line!='':
        # readline() only returns '' after the last line has been read
        
        # Throw away everything written after a '#' (It's a comment)
        line=line.split('#')
        line=line[0]

        # Delete the newline at the end of the line and divide the line in words
        line=line.strip('\n')
        line=line.split()

        # Check the first word(s) to assign the instruction to a category
        if len(line)>0:

            while line[-1]=='\\':
                # \ means that the instruction continues in the next line
                line.pop() # Remove final '\'
                
                # Read next line from the input file
                newline=source.readline()
                numline=numline+1
                newline=newline.split('#')
                newline=newline[0]
                newline=newline.strip('\n')
                newline=newline.split()
                line.extend(newline)

            # The input file is understood here
            
            if line[0]=='PROBLEMTYPE':
                # The name of the problemtype that will be generated
                projectname=line[1]
                print 'Problem type name: ',projectname
                
            elif line[0]=='DEFINITION' and line[1]=='FOLDER':
                # Read input files from the 'template_name' folder
                template_name=line[2]
                template_path = os.path.join(os.getcwd(),template_name)
                print 'Definition folder: ',template_name
                code_db,condition_classes,property_classes,material_classes,part_classes=read_tools.read_definitions(template_path)

            elif line == ['USE','KRATOS','DEFINITIONS']:
                # Use default templates for Kratos
                template_path = os.path.join(base_folder,'kratos_definitions')
                print 'Using default Kratos definitions'
                code_db,condition_classes,property_classes,material_classes,part_classes=read_tools.read_definitions(template_path)

            elif line[0]=='USE' and line[1]=='PYTHON' and line[2]=='SCRIPTS':
                # Include the python scripts found in folder as default Kratos
                # scripts fot the problemtype
                # Folder must be in the same path as the input file
                script_folder = line[3]
                
            elif line[0]=='DEFINE':
                # Create a material or part data entry
                line.pop(0)
                name,class_obj,remline=read_tools.identify_command(line,numline,material_classes,part_classes)
                if name==None:
                    stop=True
                else:
                    print 'Reading',name
                    current_material=name
                    if class_obj.call in material_classes:
                        material_dictionary[name]=class_obj(remline,numline)
                    else: # class_obj in part_classes
                        part_dictionary[name]=class_obj(remline,numline)
                
            elif line[0]=='ADD':
                # Data related to a material or part entry
                #print 'ADDING ',line
                if current_material!=None:
                    line.pop(0)
                    if current_material in material_dictionary.keys():
                        material_dictionary[current_material].add(line,numline)
                    else: # current_material in part_dictionary.keys()
                        part_dictionary[current_material].add(line,numline)
                else:
                    print 'ERROR: Found material or part data in line',numline
                    print 'but no materials or parts have been defined yet'
                    stop=True
                #print part_dictionary[current_material].conditions_always
                #print part_dictionary[current_material].conditions_2Donly
                #print part_dictionary[current_material].conditions_3Donly
                    
            elif line[0]=='GROUP':
                # Group model parts and define their interaction

                # Check if other groups have been defined before
                if current_group==None:
                    code_db.create_entry('tcl',())
                if current_group!=None:
                    # If a group was defined before, strore its code and reset the code strings
                    point_code=point_code+point_end
                    pc=core_definitions.code(point_code,'','.tcl','\t# End Point Block')
                    code_db.add_code('tcl',pc,'write')
                    point_code=''
                    point_end=''
                    
                    line_code=line_code+line_end
                    lc=core_definitions.code(line_code,'','.tcl','\t# End Line Block')
                    code_db.add_code('tcl',lc,'write')
                    line_code=''
                    line_end=''
                    
                    surf_code=surf_code+surf_end
                    sc=core_definitions.code(surf_code,'','.tcl','\t# End Surface Block')
                    code_db.add_code('tcl',sc,'write')
                    surf_code=''
                    surf_end=''
                    
                    vol_code=vol_code+vol_end
                    vc=core_definitions.code(vol_code,'','.tcl','\t# End Volume Block')
                    code_db.add_code('tcl',vc,'write')
                    vol_code=''
                    vol_end=''

                    # Delete some data stored in each part, to avoid generating duplicate code later
                    i=0
                    for entity in ('point','line','surface','volume'):
                        if entity in currententities:
                            for part in current_group:
                                store_auto_data(part_dictionary[part],auto_conditions)
                                part_dictionary[part].conditions_always[i]=()
                                part_dictionary[part].conditions_2Donly[i]=()
                                part_dictionary[part].conditions_3Donly[i]=()
                        i=i+1                            

                separator=line.index('IN')
                current_group=line[1:separator]
                entities=line[(separator+1):]

                print 'Reading Model Part group:',
                for part in current_group:
                    print part,
                print ''

                # Check that entities are valid
                currententities,remline,stop=core_definitions.read_entities(entities,numline)
                if len(remline)>0:
                    print 'Unexpected word "',remline[0],'" in line',numline
                    stop=True

                # Check that all parts in the group are valid
                group_parts=()
                for part in current_group:
                    if part not in part_dictionary:
                        print 'Unknown Model Part',part,'in line',numline
                        stop=True
                    for entity in currententities:
                        if (entity in part_dictionary[part].entities)==False:
                            print 'Error in line',numlinePart,'while creating group:\n',\
                                  part,'is not defined for',entity,'but group is'
                            stop=True
                            
                    if stop==False:
                        group_parts=group_parts+(part,)

                if ('point' in currententities)==True:
                    point_code,point_end=tcl_functions.basic_group_code('point',group_parts,part_dictionary)
                if ('line' in currententities)==True:
                    line_code,line_end=tcl_functions.basic_group_code('line',group_parts,part_dictionary)
                if ('surface' in currententities)==True:
                    surf_code,surf_end=tcl_functions.basic_group_code('surface',group_parts,part_dictionary)
                if ('volume' in currententities)==True:
                    vol_code,vol_end=tcl_functions.basic_group_code('volume',group_parts,part_dictionary)
            
            elif line.count('ASSIGN')>0:
                # Data related to part interaction
                if current_group==None:
                    print 'ERROR: Found group related data in line',numline
                    print 'but no groups have been defined yet'
                    stop=True

                separator=line.index('ASSIGN')

                # Part1 ... PartN "ASSIGN" [Cond1 ... CondN "FROM" Parti] and/or [Cond Values]

                parts=line[0:separator] # Part1 ... PartN
                actions=line[(separator+1):] # Everithing after "ASSIGN"

                for part in parts:
                    if current_group!=None:
                        if part not in current_group:
                            print '\nError  in line',numline,': Part',part,'does not belong to current group\n'
                            stop=True

                if 'point' in currententities:
                    point_code=point_code+tcl_functions.if_block('point',current_group,*parts)
                if 'line' in currententities:
                    line_code=line_code+tcl_functions.if_block('line',current_group,*parts)
                if 'surface' in currententities:
                    surf_code=surf_code+tcl_functions.if_block('surface',current_group,*parts)
                if 'volume' in currententities:
                    vol_code=vol_code+tcl_functions.if_block('volume',current_group,*parts)

                print '\tReading combination',
                for part in parts:
                    print part,
                print ''

                while len(actions)>0:#for i in range(len(actions)):
                    word=actions.pop(0)
                    if word in condition_dictionary.keys():
                        currentconds.append(word)
                    elif word=='FROM':
                        # extract next word from the list, (must be a model part name used before in the line) and write tcl code
                        part=actions.pop(0)
                        if part in parts:
                            print '\t\tAssign',
                            for c in currentconds:
                                print c,
                            print 'from',part
                            if 'point' in currententities:
                                point_code=point_code+tcl_functions.assign_from_part('point',part,currentconds,part_dictionary[part])
                            if 'line' in currententities:
                                line_code=line_code+tcl_functions.assign_from_part('line',part,currentconds,part_dictionary[part])
                            if 'surface' in currententities:
                                surf_code=surf_code+tcl_functions.assign_from_part('surface',part,currentconds,part_dictionary[part])
                            if 'volume' in currententities:
                                vol_code=vol_code+tcl_functions.assign_from_part('volume',part,currentconds,part_dictionary[part])
                            currentconds=list()

                    else:                                                        
                        if len(currentconds)==1:
                            # Assign condition, treat all words before next condition name as default values
                            remline=list()
                            remline.extend(actions)
                            condvalues=(word,)
                            for wrd in remline:
                                if wrd in condition_dictionary:
                                    break
                                else:
                                    condvalues=condvalues+(actions.pop(0),)
                                    # We remove the words from the "actions" list as we read them, so the position of the current word in the list remains constant

                            values='1 '+condition_dictionary[currentconds[0]].valuestring(condvalues)

                            if 'point' in currententities:
                                point_code=point_code+'\t\t\tGiD_AssignData Condition point_'+currentconds[0]+' points "'+values+'" $point\n'
                            if 'line' in currententities:
                                line_code=line_code+'\t\t\tGiD_AssignData Condition line_'+currentconds[0]+' lines "'+values+'" $line\n'
                            if 'surface' in currententities:
                                surf_code=surf_code+'\t\t\tGiD_AssignData Condition surface_'+currentconds[0]+' surfaces "'+values+'" $surface\n'
                            if 'volume' in currententities:
                                vol_code=col_code+'\t\t\tGiD_AssignData Condition volume_'+currentconds[0]+' volumes "'+values+'" $volume\n'

                            print '\t\tAssign',currentconds[0]
                            currentconds=list()
                            condvalues=()
                        else:
                            print 'Wrong syntax in line',numline
                            stop=True

            else:
                # Define a condition-like entry (or error in line)
                name,class_obj,remline=read_tools.identify_command(line,numline,condition_classes,property_classes)
                if name==None:
                    stop=True
                else:
                    print 'Reading',name
                    if class_obj.call in condition_classes:
                        condition_dictionary[name]=class_obj()
                        condition_dictionary[name].add(remline,numline)
                    else: # class_obj in property_classes
                        property_dictionary[name]=class_obj()
                        property_dictionary[name].add(remline,numline)                    

        # If an error is found, stop reading
        if stop==True:
            print 'An error occurred, file generation stopped'
            break
        
        # Read next line from the input file
        line=source.readline()
        numline=numline+1

    # Once the entire input file has been read, close it and start writing files
    source.close()

    # If the input file was correct, finish code and write it
    if stop==False:
        # If groups where defined, finish code for the last group
        if current_group!=None:
            point_code=point_code+point_end
            pc=core_definitions.code(point_code,'','.tcl','\t# End Point Block')
            code_db.add_code('tcl',pc,'write')

            line_code=line_code+line_end
            lc=core_definitions.code(line_code,'','.tcl','\t# End Line Block')
            code_db.add_code('tcl',lc,'write')

            surf_code=surf_code+surf_end
            sc=core_definitions.code(surf_code,'','.tcl','\t# End Surface Block')
            code_db.add_code('tcl',sc,'write')

            vol_code=vol_code+vol_end
            vc=core_definitions.code(vol_code,'','.tcl','\t# End Volume Block')
            code_db.add_code('tcl',vc,'write')

            # Delete some data stored in each part, to avoid generating duplicate code later
            i=0
            for entity in ('point','line','surface','volume'):
                if entity in currententities:
                    for part in current_group:
                        store_auto_data(part_dictionary[part],auto_conditions)
                        part_dictionary[part].conditions_always[i]=()
                        part_dictionary[part].conditions_2Donly[i]=()
                        part_dictionary[part].conditions_3Donly[i]=()
                i=i+1

        # Write tcl code for every entity+part combination that hasn't had its conditions assigned using groups
        if current_group==None:
            code_db.create_entry('tcl',())

        point_code='\tforeach point $pointlist {\n'
        line_code='\tforeach line $linelist {\n'
        surf_code='\tforeach surface $surfacelist {\n'
        vol_code='\tforeach volume $volumelist {\n'
        used_point=False
        used_line=False
        used_surf=False
        used_vol=False
        
        for part,part_obj in part_dictionary.iteritems():
            store_auto_data(part_obj,auto_conditions)
            if len(part_obj.conditions_always[0])>0 or \
               len(part_obj.conditions_2Donly[0])>0 or \
               len(part_obj.conditions_3Donly[0])>0:
                point_code=point_code+tcl_functions.basic_single_code('point',part,part_obj)
                used_point=True
            if len(part_obj.conditions_always[1])>0 or \
               len(part_obj.conditions_2Donly[1])>0 or \
               len(part_obj.conditions_3Donly[1])>0:
                line_code=line_code+tcl_functions.basic_single_code('line',part,part_obj)
                used_line=True
            if len(part_obj.conditions_always[2])>0 or \
               len(part_obj.conditions_2Donly[2])>0 or \
               len(part_obj.conditions_3Donly[2])>0:
                surf_code=surf_code+tcl_functions.basic_single_code('surface',part,part_obj)
                used_surf=True
            if len(part_obj.conditions_always[3])>0 or \
               len(part_obj.conditions_2Donly[3])>0 or \
               len(part_obj.conditions_3Donly[3])>0:
                vol_code=vol_code+tcl_functions.basic_single_code('volume',part,part_obj)
                used_vol=True

        if used_point==True:
            point_code=point_code+'\t}\n'
            pc=core_definitions.code(point_code,'','.tcl','\t# End Point Block')
            code_db.add_code('tcl',pc,'write')
            
        if used_line==True:
            line_code=line_code+'\t}\n'
            lc=core_definitions.code(line_code,'','.tcl','\t# End Line Block')
            code_db.add_code('tcl',lc,'write')

        if used_surf==True:
            surf_code=surf_code+'\t}\n'
            sc=core_definitions.code(surf_code,'','.tcl','\t# End Surface Block')
            code_db.add_code('tcl',sc,'write')

        if used_vol==True:
            vol_code=vol_code+'\t}\n'
            vc=core_definitions.code(vol_code,'','.tcl','\t# End Volume Block')
            code_db.add_code('tcl',vc,'write')

        # To avoid problems when remeshing, all data assigned via tcl must be erased before generating a new mesh
        # Generate code to clean model parts:
        clean_tcl=''
        for name,part_obj in part_dictionary.iteritems():
            clean_tcl=clean_tcl+'\tcleanautomatic '+name
            for entity in part_obj.entities:
                clean_tcl=clean_tcl+' '+entity
            clean_tcl=clean_tcl+'\n'
        # Generate code to clean conditions and elements
        for condition,entities in auto_conditions.iteritems():
            clean_tcl=clean_tcl+'\tcleanautomatic '+condition
            for entity in set(entities):
                clean_tcl=clean_tcl+' '+entity
            clean_tcl=clean_tcl+'\n'
        # Write code
        clean_code=core_definitions.code(clean_tcl,'','.tcl','\t# End Reset Block')
        code_db.add_code('tcl',clean_code,'write')

        # Write files
        print '\nGenerating Files...\n'
        file_functions.generate_files(projectname,template_path)
        # if a folder containing default python scripts for GiD was provided, include
        # them in the problemtype
        if script_folder != None:
            file_functions.add_python_scripts(projectname,script_folder)

        # Once the problemtype folder has been created and filled with default
        # files, move there (set it as current working directory)
        os.chdir(projectname+'.gid/')
        
        code_db.clean()
        code_db.write(projectname)

        # Purge unused books (GiD won't work properly with them) and create a custom Tcl Menu
        menustring=file_functions.check_books(projectname,template_path)
        menucode=core_definitions.code(menustring,'','.tcl','\t# Custom Menu\n')
        menucode.write(projectname)
        
        endstring='\nPROBLEMTYPE GENERATED\n\nCopy the '+projectname+'.gid folder to GiD\'s problemtypes folder.'
        print endstring

def store_auto_data(part_obj,auto_dict):
    for entity,i in zip(('point','line','surface','volume'),range(4)):
        for condition in part_obj.conditions_always[i]:
            try:
                auto_dict[condition].append(entity)
            except KeyError:
                auto_dict[condition]=[entity]
        for condition in part_obj.conditions_2Donly[i]:
            try:
                auto_dict[condition].append(entity)
            except KeyError:
                auto_dict[condition]=[entity]
        for condition in part_obj.conditions_3Donly[i]:
            try:
                auto_dict[condition].append(entity)
            except KeyError:
                auto_dict[condition]=[entity]

if __name__=='__main__':
    import sys
    import os

    # Path to the folder containing this script
    basefolder = os.path.split(sys.argv[0])[0]
    try:
        inputfile = sys.argv[1]
    except: inputfile = None
    if inputfile:
        import problemtype
        problemtype.generate(inputfile,basefolder)
    else:
        print 'Use python problemtype.py <input file>'            
