import re
import sys
import os
import inspect
import core_definitions

# Compile required regular expressions
# (see http://www.amk.ca/python/howto/regex/regex.html for an introduction to
# regular expressions in python)

# The following re are used to parse input files
code_blocks=re.compile("""
^file[\s]+  # block starts with "file" and one or more whitespace characters
(.*?)problemtype(.*?)$  # the filename is any remaining text before the end
[\s]*(.*?)$             # the code itself (with or without an insert point)
(?=[\s]*^file[\s]+|[\s]*\Z)
# block ends before the next block starts or at the end of the file
""",re.VERBOSE | re.DOTALL | re.MULTILINE)

code_pieces=re.compile("""
^where[\s]+(.*?)$               # starts with "where" and the insertion point
[\n]+(.*?)$                     # code
(?=[\s]*^where[\s]+|[\s]*\Z)    # ends with a new "where" or at string's end
""",re.VERBOSE | re.DOTALL | re.MULTILINE)


### Read the definition of the problemtype from def_file_name.py

def read_definitions(definitions_path):

    # Create a dictionary for each basic data type    
    condition_dict=dict() # GiD conditional data (conditions and elements)
    property_dict=dict() # GiD general data
    material_dict=dict() # GiD materials
    part_dict=dict() # Model parts
    
    template_code=core_definitions.code_container()
    # A core_definitions.code_container object, stores code read from definition files

    sys.path.append(definitions_path)
    import new_classes
    class_list=inspect.getmembers(new_classes,inspect.isclass)

    for class_name,class_item in class_list:
        template_code.create_entry(class_item.call,list())
        definitions_folder = os.path.join(definitions_path,'definitions')

        if class_item.definition_file!='': # If the class uses a definition file
            def_file = os.path.join(definitions_folder,class_item.definition_file)
            f=open(def_file)
            templates=f.read()
            f.close()

            # Parse the file using regular expressions
            matches=code_blocks.findall(templates)
            # Store read input as code templates
            for item in matches:
                pref=item[0]
                suff=item[1]
                if suff[-4:]=='.bas': # This is a GiD template file
                    filetype='.bas'
                else: #elif suff[-4:]!='.bas':
                    filetype='other'
                pieces=code_pieces.findall(item[2])
                if len(pieces)==0:
                    # A single code block without "where" statements
                    bascode=item[2]+'\n'
                    newcode=core_definitions.code(bascode,pref,suff,'')
                    template_code.add_code(class_item.call,newcode,filetype)
                else:
                    # multiple code blocks in the same file
                    for code_part in pieces:
                        bascode=code_part[1]+'\n'
                        newcode=core_definitions.code(bascode,pref,suff,code_part[0])
                        template_code.add_code(class_item.call,newcode,filetype)

        # Assign all defined data to one of the four dictionaries: conditions/elements, materials
        # properties and model parts. This classification is done based on the base classes of each
        # defined class. Note that conditions, elements and properties require a single line in the
        # input file to be defined, while materials and model parts require several lines.

        base_classes=inspect.getmro(class_item)
        if (core_definitions.condition in base_classes) or (core_definitions.element in base_classes):
            condition_dict[class_item.call]=class_item
        elif core_definitions.gendata in base_classes:
            property_dict[class_item.call]=class_item
        elif core_definitions.material in base_classes:
            material_dict[class_item.call]=class_item
        elif core_definitions.part in base_classes:
            part_dict[class_item.call]=class_item

    # Add some extra model part classes defined in core_definitions.py to their dictionary
    # Only classes with a defined 'call' attribute will be used (Default call attribute is None)
    default_classes=inspect.getmembers(core_definitions,inspect.isclass)
    for class_name,class_item in default_classes:
        
        base_classes=inspect.getmro(class_item)
        try:
            call_attr=class_item.call
        except AttributeError:
            call_attr=None
            
        if call_attr!=None:
            if (core_definitions.condition in base_classes) or (core_definitions.element in base_classes):
                condition_dict[class_item.call]=class_item
            elif core_definitions.gendata in base_classes:
                property_dict[class_item.call]=class_item
            elif core_definitions.material in base_classes:
                material_dict[class_item.call]=class_item
            elif core_definitions.part in base_classes:
                part_dict[class_item.call]=class_item

    # Add the 'option' class from core_definitions
    condition_dict[core_definitions.option.call]=core_definitions.option

    return template_code,condition_dict,property_dict,material_dict,part_dict

### Basic command identifyinig: command and entities

# Identify, based on the first words of a line from input file, which kind of
# information is being given

def identify_command(line,numline,*class_dict):
    for dictionary in class_dict:
        for command,class_obj in dictionary.iteritems():
            expected_words=command.split()
            com_len=len(expected_words)
            head=line[0:com_len]
            if expected_words==head:
                name=line[com_len]
                remline=line[com_len:]
                return name,class_obj,remline
    # The following lines are only executed if the command is unknown
    print 'Unknown or unexpected instruction in line',numline
    return None,None,None
