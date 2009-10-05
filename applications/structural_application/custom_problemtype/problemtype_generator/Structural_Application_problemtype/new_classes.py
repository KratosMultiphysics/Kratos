import core_definitions
# This file contains the basic data classes

class condition_scalar(core_definitions.condition):
    # A scalar condition
    call='SCALAR CONDITION'
    definition_file='scalar_condition'
    insert_in='# Nodal Values'

    questions='QUESTION: Fixed#CB#(1,0)\n'+\
               'VALUE: <TYPE>\n'+\
               'QUESTION: <NAME>\n'+\
               'VALUE: <DEFVALUE>\n'

    additional_input=('<TYPE>','<DEFVALUE>')

    def parseinput(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No default values found for this condition'
            return
        elif extra_input[0]=='fixed':
            fixed='1'
            value=extra_input[1]
            if len(extra_input)>2:
                print 'More default values than expected found, only the first two will be used'
                print '\t'+str(extra_input)
        elif extra_input[0]=='free':
            fixed='0'
            value=extra_input[1]
            if len(extra_input)>2:
                print 'More default values than expected found, only the first two will be used'
                print '\t'+str(extra_input)
        else:
            fixed='1'
            value=extra_input[0]
            if len(extra_input)>1:
                print 'More default values than expected found, only the first one will be used'
                print '\t'+str(extra_input)
            extra_input=extra_input[0:2]
        return fixed,value
    
    def valuestring(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No values found for this condition'
            return
        elif extra_input[0]=='fixed':
            fixed='1'
            value=extra_input[1]
            if len(extra_input)>2:
                print 'More values than expected found, only the first two will be used'
                print '\t'+str(extra_input)
        elif extra_input[0]=='free':
            fixed='0'
            value=extra_input[1]
            if len(extra_input)>2:
                print 'More values than expected found, only the first two will be used'
                print '\t'+str(extra_input)
        else:
            fixed='1'
            value=extra_input[0]
            if len(extra_input)>1:
                print 'More values than expected found, only the first one will be used'
                print '\t'+str(extra_input)
            extra_input=extra_input[0:2]
        return fixed+' '+value

class condition_flag(condition_scalar):
    # A flag condition
    call='FLAG CONDITION'
    definition_file='scalar_condition'
    insert_in='# Nodal Values'

    questions='QUESTION: Fixed#CB#(1,0)\n'+\
               'VALUE: <TYPE>\n'+\
               'QUESTION: <NAME>#CB#(<VALUES>)\n'+\
               'VALUE: <DEFVALUE>\n'

    additional_input=('<TYPE>','<DEFVALUE>','<VALUES>')

    def parseinput(self,extra_input):
        # This method processes input from flag definitions in the input file,
        # and returns a tuple containing any additional parameters required
        # (ordered as in additional_input)
        if extra_input[0]=='free':
            cond_type='0'
            extra_input=extra_input[1:]
        elif extra_input[0]=='fixed':
            cond_type='1'
            extra_input=extra_input[1:]
        else:
            cond_type='1'
        if len(extra_input)==0:
            print 'ERROR: No default values found for this condition'
            return
        defvalue=str(extra_input[0])
        first=True
        for val in extra_input:
            if first==False:
                values=values+','+str(val)
            else:
                values=str(val)
                first=False
        return cond_type,defvalue,values

class condition_vector(core_definitions.condition):
    # A vector (initial or boundary) condition. It has 7 question fields:
    # fixed: 1 (boundary condition) or 0 (initial condition)
    # For each direction (X,Y,Z):
    # - <NAME>_X: A flag. 1 if the variable is prescribed in this direction
    # - X_Value: Prescribed value along direction X
    call='VECTOR CONDITION'
    definition_file='vector_condition'
    insert_in='# Nodal Values'

    questions='QUESTION: <NAME>_X#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,X_Value,0.0,SET,Fix_X,#CURRENT#)(#DEFAULT#,RESTORE,X_Value,#CURRENT#,RESTORE,Fix_X,#CURRENT#)\n'+\
               'QUESTION: Fix_X#CB#(1,0)\n'+\
               'VALUE: <FIXX>\n'+\
               'QUESTION: X_Value\n'+\
               'VALUE: <VALX>\n'+\
               'QUESTION: <NAME>_Y#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,Y_Value,0.0,SET,Fix_Y,#CURRENT#)(#DEFAULT#,RESTORE,Y_Value,#CURRENT#,RESTORE,Fix_Y,#CURRENT#)\n'+\
               'QUESTION: Fix_Y#CB#(1,0)\n'+\
               'VALUE: <FIXY>\n'+\
               'QUESTION: Y_Value\n'+\
               'VALUE: <VALY>\n'+\
               'QUESTION: <NAME>_Z#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,Z_Value,0.0,SET,Fix_Z,#CURRENT#)(#DEFAULT#,RESTORE,Z_Value,#CURRENT#,RESTORE,Fix_Z,#CURRENT#)\n'+\
               'QUESTION: Fix_Z#CB#(1,0)\n'+\
               'VALUE: <FIXZ>\n'+\
               'QUESTION: Z_Value\n'+\
               'VALUE: <VALZ>\n'

    additional_input=('<FIXX>','<VALX>','<FIXY>','<VALY>','<FIXZ>','<VALZ>')

    def parseinput(self,extra_input):
        if len(extra_input)<=4:
            if extra_input[0]=='fixed':
                fixed='1'
                extra_input=extra_input[1:]
            elif extra_input[0]=='free':
                fixed='0'
                extra_input=extra_input[1:]
            else: # no fixed/free argument
                fixed='1'

            if len(extra_input)<3:
                print 'ERROR: Not enough default values found for this condition'
                print '\t'+str(extra_input)
                return
            elif len(extra_input)>3:
                'More default values than expected found, some will be ignored'
                print '\t'+str(extra_input)
                extra_input=extra_input[0:3]

            ValX=extra_input[0]
            ValY=extra_input[1]
            ValZ=extra_input[2]

            return fixed,ValX,fixed,ValY,fixed,ValZ

        elif len(extra_input)==6:
            if extra_input[0]=='fixed':
                fixX='1'
            elif extra_input[0]=='free':
                fixX='0'
            ValX=extra_input[1]
            
            if extra_input[2]=='fixed':
                fixY='1'
            elif extra_input[2]=='free':
                fixY='0'
            ValY=extra_input[3]

            if extra_input[4]=='fixed':
                fixZ='1'
            elif extra_input[4]=='free':
                fixZ='0'
            ValZ=extra_input[5]

            return fixX,ValX,fixY,ValY,fixZ,ValZ

    def valuestring(self,extra_input):
        if len(extra_input)<=4:
            if extra_input[0]=='fixed':
                fixed='1'
                extra_input=extra_input[1:]
            elif extra_input[0]=='free':
                fixed='0'
                extra_input=extra_input[1:]
            else: # no fixed/free argument
                fixed='1'

            if len(extra_input)<3:
                print 'ERROR: Not enough default values found for this condition'
                print '\t'+str(extra_input)
                return # This should crash the problem type generator
            elif len(extra_input)>3:
                'More default values than expected found, some will be ignored'
                print '\t'+str(extra_input)
                extra_input=extra_input[0:3]

            ValX=extra_input[0]
            ValY=extra_input[1]
            ValZ=extra_input[2]

            return '1 '+fixed+' '+ValX+' 1 '+fixed+' '+ValY+' 1 '+fixed+' '+ValZ
        
        elif len(extra_input)==6:
            if extra_input[0]=='fixed':
                fixX='1 '
            elif extra_input[0]=='free':
                fixX='0 '
            ValX=extra_input[1]
            
            if extra_input[2]=='fixed':
                fixY='1 '
            elif extra_input[2]=='free':
                fixY='0 '
            ValY=extra_input[3]

            if extra_input[4]=='fixed':
                fixZ='1 '
            elif extra_input[4]=='free':
                fixZ='0 '
            ValZ=extra_input[5]

            return '1 '+fixX+ValX+' 1 '+fixY+ValY+' 1 '+fixZ+ValZ

class condition_default(core_definitions.condition):
    # A single-field condition. If it's not assigned by the user, it will be set
    # to it's default value (which can be modified from the properies menu)

    call='DEFAULT CONDITION'
    definition_file='default_condition'
    insert_in='# Nodal Values'
    questions='QUESTION: <NAME>\n'+\
               'VALUE: <VALUE>\n'

    additional_input=('<VALUE>','<ENTITIES>')

    def valuestring(self,value):
        return str(value)

    def parseinput(self,extra_input):
        entitylist=''
        first=True
        for entity in self.entities:
            if first==False:
                entitylist=entitylist+' '
            else:
                first=False
            entitylist=entitylist+entity
        return extra_input[0],entitylist

class element(core_definitions.element):
    # A finite element type
    call='ELEMENT'
    definition_file='element'
    insert_in='# Elements'

    additional_input=('<ELEMTYPE>','<TCL_COND>')

    def parseinput(self,extra_input):
        ent_name=self.entities[0]
        if len(extra_input)>1 and extra_input!=('Only','Points'):
            print 'Additional input must be the element type, but the following was found instead:'
            print '\t'+str(extra_input)
            print 'Only '+extra_input[0]+' will be used'
            extra_input=(extra_input[0],)
        elif extra_input==('Only','Points'):
            extra_input=('OnlyPoints',)
        if len(extra_input)==1 and (ent_name in ('surface','volume')):
            if ent_name=='surface' and (extra_input[0] not in ('Triangle','Quadrilateral','Linear','Circle')):
                print 'ERROR: Only \'Triangle\', \'Quadrilateral\', \'Linear\' or \'Circle\' can be used as surface element types'
                return
            if ent_name=='volume' and (extra_input[0] not in ('Tetrahedra','Linear','Hexahedra','Prism','OnlyPoints','Sphere')):
                print 'ERROR: Only \'Tetrahedra\', \'Linear\', \'Hexahedra\', \'Prism\', \'Only Points\' or \'Sphere\' can be used as volume element types'
                return
            elemtype=extra_input[0]
            tcl_cond='[llength [GiD_Info list_entities conditions '+ent_name+'_'+self.name+' ]]>0'
        else: #elif len(extra_input)==0: or this is a line element
            tcl_cond='0'
            elemtype='None'

        return elemtype,tcl_cond

class point_element(core_definitions.element):
    # Due to limitations in GiD, elements applied over point entities must be build from their own templates
    # Note that this doesn't apply to volumes meshed with the 'OnlyPoints' element type, those can and should use the standard element template
    call='POINT ELEMENT'
    definition_file='point_element'
    insert_in='# Elements'
    apply_over_nodes=True # As of GiD 9.1.1b, point conditions can't be applied over elements

    questions='QUESTION: ID#FUNC#(NumEntity)\n'+\
               'VALUE: 0\n'+\
               'STATE: HIDDEN\n'
##    questions='QUESTION: Material\n'+\
##               'VALUE: None\n'+\
##               'TKWIDGET: TkwidgetGetMaterial\n'

    def valuestring(self,extra_input):
        return '0'
##        return 'None'
    
class face_condition(element):
    # A condition to identify model faces (surfaces in models with volumes, lines in those without volumes)
    call='FACE CONDITION'
    definition_file='face_condition'
    insert_in='# Conditions'

class point_condition(point_element):
    # See comment on point_element
    call='POINT CONDITION'
    definition_file='point_condition'
    insert_in='# Conditions'

class material(core_definitions.material):
    # Template for material properties
    call='MATERIAL'
    definition_file='material'

class python_property(core_definitions.gendata):
    call='PROPERTY'
    definition_file='py_property'
    additional_input=('<VALUE>','<VAR_NAME>')

class python_string(core_definitions.gendata):
    call='TEXT PROPERTY'
    definition_file='py_string'
    additional_input=('<VALUE>','<VAR_NAME>')

class python_flag(core_definitions.gendata):
    call='FLAG PROPERTY'
    definition_file='py_property'

    questions='QUESTION: <NAME>#CB#(<VALUES>)\n'+\
               'VALUE: <DEFVALUE>\n'
    
    additional_input=('<DEFVALUE>','<VALUES>','<VAR_NAME>')

    def parseinput(self,extra_input):
        if len(extra_input)<2:
            print 'ERROR: Not enough values found for this flag property'
            return

        varname=str(extra_input[-1])

        defvalue=str(extra_input[0])
        first=True
        for val in extra_input[0:-1]:
            if first==False:
                values=values+','+str(val)
            else:
                values=str(val)
                first=False
        return defvalue,values,varname

class python_text_flag(python_flag):
    call='TEXT FLAG PROPERTY'
    definition_file='py_string'

class file_path(core_definitions.gendata):
    call='FILE SELECTION'
    definition_file='py_filepath'
    additional_input=('<VALUE>','<VAR_NAME>')
    
    questions='QUESTION: <NAME>\n'+\
               'VALUE: <VALUE>\n'+\
               'TKWIDGET: TkwidgetFilePath\n'

class scalar_face_condition(core_definitions.condition):
    # This is a scalar condition that can be applied over elements
    call='SCALAR FACE VALUE'
    definition_file='scalar_conditional_var'
    insert_in='# Conditional Data'
    apply_over_nodes=False

    questions='QUESTION: <NAME>\n'+\
               'VALUE: <DEFVALUE>\n'

    additional_input=('<DEFVALUE>',)

    def parseinput(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No default values found for this condition'
            return
        elif len(extra_input)>1:
            print 'More default values than expected found, only the first one will be used'
            print '\t'+str(extra_input)
        return (extra_input[0],)

    def valuestring(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No values found for this condition'
            return
        elif len(extra_input)>1:
            print 'More values than expected found, only the first one will be used'
            print '\t'+str(extra_input)
        return extra_input[0]

class flag_face_condition(scalar_face_condition):
    # This is a flag condition that can be applied over elements
    call='FLAG FACE VALUE'
    definition_file='scalar_conditional_var'
    insert_in='# Conditional Data'
    apply_over_nodes=False

    questions='QUESTION: <NAME>#CB#(<VALUES>)\n'+\
               'VALUE: <DEFVALUE>\n'

    additional_input=('<DEFVALUE>','<VALUES>')

    def parseinput(self,extra_input):
        if len(extra_input)==0:
            print 'ERROR: No default values found for this condition'
            return
        defvalue=str(extra_input[0])
        first=True
        for val in extra_input:
            if first==False:
                values=values+','+str(val)
            else:
                values=str(val)
                first=False
        return defvalue,values

class vector_face_condition(core_definitions.condition):
    # This is a copy of the flag condition, but applied over the element
    call='VECTOR FACE VALUE'
    definition_file='vector_conditional_var'
    insert_in='# Conditional Data'
    apply_over_nodes=False

    questions='QUESTION: <NAME>_X#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,X_Value,0.0)(#DEFAULT#,RESTORE,X_Value,#CURRENT#)\n'+\
               'QUESTION: X_Value\n'+\
               'VALUE: <VALX>\n'+\
               'QUESTION: <NAME>_Y#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,Y_Value,0.0)(#DEFAULT#,RESTORE,Y_Value,#CURRENT#)\n'+\
               'QUESTION: Y_Value\n'+\
               'VALUE: <VALY>\n'+\
               'QUESTION: <NAME>_Z#CB#(1,0)\n'+\
               'VALUE: 1\n'+\
               'DEPENDENCIES: (0,SET,Z_Value,0.0)(#DEFAULT#,RESTORE,Z_Value,#CURRENT#)\n'+\
               'QUESTION: Z_Value\n'+\
               'VALUE: <VALZ>\n'

    additional_input=('<VALX>','<VALY>','<VALZ>')

    def parseinput(self,extra_input):
        
        if len(extra_input)<3:
            print 'ERROR: Not enough default values found for this condition'
            print '\t'+str(extra_input)
            return
        elif len(extra_input)>3:
            'More default values than expected found, some will be ignored'
            print '\t'+str(extra_input)
            extra_input=extra_input[0:3]

        ValX=extra_input[0]
        ValY=extra_input[1]
        ValZ=extra_input[2]

        return fixed,ValX,ValY,ValZ

    def valuestring(self,extra_input):
            
        if len(extra_input)<3:
            print 'ERROR: Not enough default values found for this condition'
            return # This should crash the problem type generator
        elif len(extra_input)>3:
            'More default values than expected found, some will be ignored'
            extra_input=extra_input[0:3]

        ValX=extra_input[0]
        ValY=extra_input[1]
        ValZ=extra_input[2]

        return fixed+' 1 '+ValX+' 1 '+ValY+' 1 '+ValZ
    
class scalar_elemental_value(scalar_face_condition):
    # This is a copy of the scalar condition applied over condition faces
    call='SCALAR ELEMENTAL VALUE'
    definition_file='scalar_elemental_var'
    insert_in='# Elemental Data'

class flag_elemental_condition(flag_face_condition):
    # This is a copy of the flag condition applied over condition faces
    call='FLAG ELEMENTAL VALUE'
    definition_file='scalar_elemental_var'
    insert_in='# Elemental Data'

class vector_elemental_condition(vector_face_condition):
    # This is a copy of the flag condition applied over condition faces
    call='VECTOR ELEMENTAL VALUE'
    definition_file='vector_elemental_var'
    insert_in='# Elemental Data'

class kratos_solver(core_definitions.gendata):
    call='SOLVER'
    definition_file='property_text'
    additional_input=('<VAR_NAME>',)
    questions='QUESTION: Solver#CB#(StaticSolver, DynamicSolver ,ParallelSolver)\n'+\
               'VALUE: <VAR_NAME>\n'

