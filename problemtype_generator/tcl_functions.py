def basic_group_code(entity,group,parts_dict): #(string,list,dict)
    
    head='\tforeach '+entity+' $'+entity+'list {\n'+\
          '\t\tif {\n'

    i=0
    # First Check: the entity hasn't any model part from this group assigned (do nothing)
    for part in group:
        i=i+1
        head=head+'\t\t\t[lsearch ${'+entity+part+'list} $'+entity+'] == -1'

        if i==len(group):
            head=head+'\n'
        elif i>0:
            head=head+' && \\\n'

    head=head+'\t\t} then {\n' # Do nothing

    # Now check for entities belonging to a single part
    for part in group:
        
        i=0

        head=head+'\t\t} elseif {\n'

        for currentpart in group:
            i=i+1

            if part==currentpart:
                head=head+'\t\t\t[lsearch ${'+entity+currentpart+'list} $'+entity+'] != -1'
            else:
                head=head+'\t\t\t[lsearch ${'+entity+currentpart+'list} $'+entity+'] == -1'

            if i==len(group):
                head=head+'\n'
            elif i>0:
                head=head+' && \\\n'

        head=head+'\t\t} then {\n'+\
              '\t\t\tset '+part+'pos [lsearch ${'+entity+part+'list} $'+entity+']\n'+\
              '\t\t\tlreplace ${'+entity+part+'list} ${'+part+'pos} ${'+part+'pos}\n'

        entities=['point','line','surface','volume']
        index=entities.index(entity)

        for condition in parts_dict[part].conditions_always[index]:
            head=head+'\t\t\tcondfrompart '+part+' '+condition+' always '+entity+' $'+entity+'\n'
        if len(parts_dict[part].conditions_2Donly[index])>0:
            head=head+'\t\t\tif {[GiD_Info Geometry NumVolumes]==0} {\n'
            for condition in parts_dict[part].conditions_2Donly[index]:
                head=head+'\t\t\t\tcondfrompart '+part+' '+condition+' only2D '+entity+' $'+entity+'\n'
            head=head+'\t\t\t}\n'
        if len(parts_dict[part].conditions_3Donly[index])>0:
            head=head+'\t\t\tif {[GiD_Info Geometry NumVolumes]>0} {\n'
            for condition in parts_dict[part].conditions_2Donly[index]:
                head=head+'\t\t\t\tcondfrompart '+part+' '+condition+' only3D '+entity+' $'+entity+'\n'
            head=head+'\t\t\t}\n'

    # Specific part combinations will be added between the head and end parts

    # Final case: unexpected combination of parts
    end='\t\t} else {\n'+\
         '\t\t\tputs "Unexpected combination of Model Parts in '+entity+' $'+entity+'"\n\t\t}\n\t}\n'

    return head,end

def basic_single_code(entity,part,parts_obj): #(string,string,class instance)

    # Check for entities belonging to the part
    head='\t\tif {[lsearch ${'+entity+part+'list} $'+entity+'] != -1} then {\n'+\
          '\t\t\tset '+part+'pos [lsearch ${'+entity+part+'list} $'+entity+']\n'+\
          '\t\t\tlreplace ${'+entity+part+'list} ${'+part+'pos} ${'+part+'pos}\n'

    entities=['point','line','surface','volume']
    index=entities.index(entity)

    for condition in parts_obj.conditions_always[index]:
        head=head+'\t\t\tcondfrompart '+part+' '+condition+' always '+entity+' $'+entity+'\n'
    if len(parts_obj.conditions_2Donly[index])>0:
        head=head+'\t\t\tif {[GiD_Info Geometry NumVolumes]==0} then {\n'
        for condition in parts_obj.conditions_2Donly[index]:
            head=head+'\t\t\t\tcondfrompart '+part+' '+condition+' only2D '+entity+' $'+entity+'\n'
        head=head+'\t\t\t}\n'
    if len(parts_obj.conditions_3Donly[index])>0:
        head=head+'\t\t\tif {[GiD_Info Geometry NumVolumes]>0} then {\n'
        for condition in parts_obj.conditions_3Donly[index]:
            head=head+'\t\t\t\tcondfrompart '+part+' '+condition+' only3D '+entity+' $'+entity+'\n'
        head=head+'\t\t\t}\n'
    head=head+'\t\t}\n'

    return head
            
def if_block(entity,group,*parts):# (str,list,*str)

    code='\t\t} elseif {\n'

    i=0

    for currentpart in group:
        i=i+1
        
        if (currentpart in parts)==True:
            code=code+'\t\t\t[lsearch ${'+entity+currentpart+'list} $'+entity+'] != -1'
        else:
            code=code+'\t\t\t[lsearch ${'+entity+currentpart+'list} $'+entity+'] == -1'

        if i==len(group):
            code=code+'\n'
        elif i>0:
            code=code+' && \\\n'

    code=code+'\t\t} then {\n'

    for part in parts:
        code=code+'\t\t\tset '+part+'pos [lsearch ${'+entity+part+'list} $'+entity+']\n'+\
              '\t\t\tlreplace ${'+entity+part+'list} ${'+part+'pos} ${'+part+'pos}\n'
    
    return code

def assign_from_part(entity,part,conds,part_obj):# (string,string,list,class instance)

    entities=['point','line','surface','volume']
    index=entities.index(entity)

    code=''

    # Build code for conditions to be applied both in 2D and 3D problems
    for condition in conds:
        if condition in part_obj.conditions_always[index]:
            code=code+'\t\t\tcondfrompart '+part+' '+condition+' always '+entity+' $'+entity+'\n'
    # Build code for conditions to be applied in 2D problems only
    if len(part_obj.conditions_2Donly[index])>0:
        code='\t\t\tif {[GiD_Info Geometry NumVolumes]==0} {\n'
        for condition in conds:
            if condition in part_obj.conditions_2Donly[index]:
                code=code+'\t\t\t\tcondfrompart '+part+' '+condition+' only2D '+entity+' $'+entity+'\n'
        code=code+'\t\t\t}\n'
    # Build code for conditions to be applied in 3D problems only
    if len(part_obj.conditions_3Donly[index])>0:
        code='\t\t\tif {[GiD_Info Geometry NumVolumes]>0} {\n'
        for condition in conds:
            if condition in part_obj.conditions_3Donly[index]:
                code=code+'\t\t\t\tcondfrompart '+part+' '+condition+' only3D '+entity+' $'+entity+'\n'
        code=code+'\t\t\t}\n'

    return code
