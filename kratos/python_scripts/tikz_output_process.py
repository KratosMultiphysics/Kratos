# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager

import os

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TikZOutputProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class TikZOutputProcess(KratosMultiphysics.Process):
    """This process generates TikZ source files of 2D meshes

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"                    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_control_type"                : "step",
            "output_interval"                    : 1.0,
            "folder_name"                        : "tikZ_Output",
            "save_output_files_in_folder"        : true,
            "cmyk_colors_fill"                   : [0, 1.0, 1.0, 1.0],
            "include_nodes_label"                : true,
            "include_elements_label"             : true,
            "include_axis"                       : true,
            "include_grid"                       : true,
            "landscape_mode"                     : true,
            "include_caption"                    : false,
            "sans_fonts"                         : false,
            "prefix_nodes"                       : "",
            "prefix_elements"                    : ""
        }
        """
        )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[self.settings["model_part_name"].GetString()]

        # Checking dimension of the problem
        if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] != 2:
            raise Exception("Expected 2 dimensional problem")

        # Warning: we may be changing the parameters object here:
        self.TranslateLegacyVariablesAccordingToCurrentStandard(settings)

        if self.settings["save_output_files_in_folder"].GetBool():
            if self.model_part.GetCommunicator().MyPID() == 0:
                folder_name = self.settings["folder_name"].GetString()
                if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
                    kratos_utils.DeleteDirectoryIfExisting(folder_name)
                if not os.path.isdir(folder_name):
                    os.mkdir(folder_name)
            self.model_part.GetCommunicator().Barrier()

        self.output_interval = self.settings["output_interval"].GetDouble()
        self.output_control = self.settings["output_control_type"].GetString()
        self.next_output = 0.0
        self.step_count = 0

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        old_name = 'output_frequency'
        new_name = 'output_interval'

        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

    def ExecuteInitialize(self):
        """ This method is executed in order to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.output_control == "time":
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # First we retrieve the template argument
        include_nodes_label = self.settings["include_nodes_label"].GetBool()
        include_elements_label = self.settings["include_elements_label"].GetBool()
        include_axis = self.settings["include_axis"].GetBool()
        include_grid = self.settings["include_grid"].GetBool()
        landscape_mode = self.settings["landscape_mode"].GetBool()
        include_caption = self.settings["include_caption"].GetBool()
        sans_fonts = self.settings["sans_fonts"].GetBool()
        cmyk_colors_fill = self.settings["cmyk_colors_fill"].GetVector()
        template_string = GetTemplateString(self.model_part, cmyk_colors_fill, include_nodes_label, include_elements_label, include_axis, include_grid, landscape_mode, include_caption, sans_fonts)

        # Generate nodal data
        prefix_nodes = self.settings["prefix_nodes"].GetString()
        nodal_data = ""
        counter = 0
        for node in self.model_part.Nodes:
            if counter > 0:
                nodal_data += "\n"
            if prefix_nodes != "":
                nodal_data += prefix_nodes + "_"
            nodal_data += "{" + str(node.Id) + "}"
            nodal_data += "    " + "{0:.12g}".format(node.X)
            nodal_data += "    " + "{0:.12g}".format(node.Y)
            counter += 1

        # Generate elemental data
        prefix_elements = self.settings["prefix_elements"].GetString()
        elemental_data = ""
        counter = 0
        for elem in self.model_part.Elements:
            if counter > 0:
                elemental_data += "\n"
            if prefix_elements != "":
                elemental_data += prefix_elements + "_"
            elemental_data += "{" + str(elem.Id) + "}"
            for node in elem.GetNodes():
                if prefix_nodes != "":
                    elemental_data += prefix_nodes + "_"
                elemental_data += "    {" + str(node.Id) + "}"
            counter += 1

        outputstring = template_string.replace("%%REPLACE_NODE_DATA", nodal_data)
        outputstring = outputstring.replace("%%REPLACE_ELEMENT_DATA", elemental_data)
        folder_name = self.settings["folder_name"].GetString()
        model_part_name = self.settings["model_part_name"].GetString()
        name_file = os.path.join(folder_name, model_part_name)
        if self.output_control == "time":
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            with open(name_file + "_TIME_" + str(time) + ".tex", 'w') as output:
                output.write(outputstring)
        else:
            step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            with open(name_file + "_STEP_" + str(step) + ".tex", 'w') as output:
                output.write(outputstring)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.step_count += 1

    def PrintOutput(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # TODO: Add post process with values

        # Schedule next output
        time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        if self.output_interval > 0.0: # Note: if == 0, we'll just always print
            if self.output_control == "time":
                while GetPrettyTime(self.next_output) <= time:
                    self.next_output += self.output_interval
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_interval

    def IsOutputStep(self):
        """ This method determines if corresponds to output

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.output_control == "time":
            time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= GetPrettyTime(self.next_output))
        else:
            return ( self.step_count >= self.next_output )

def GetTemplateString(model_part, cmyk_colors_fill, include_nodes_label = True, include_elements_label = True, include_axis = True, include_grid = True, landscape_mode = True, include_caption = False, sans_fonts = False):
    sans_font_string = ""
    if sans_fonts:
        sans_font_string = "\\usepackage{verbatim}\n\\usepackage{avant}\n\\usepackage[scaled]{helvet}\n\\usepackage[helvet]{sfmath}\n\\usepackage[EULERGREEK]{sansmath}\n\\usepackage{floatrow}\n\DeclareFloatFont{henry}{\\sffamily\\sansmath}\n\\floatsetup[figure]{font=henry}\n\\renewcommand{\\familydefault}{\\sfdefault}\n"
    template_string = "\documentclass{article}\n\\usepackage{tikz,listofitems,readarray,filecontents}\n\\usepackage{pdflscape}\n\\pagenumbering{gobble}\n\\usepackage{xcolor}\n\\usetikzlibrary{calc}\n" + sans_font_string + "\\begin{filecontents*}{nodedata.dat}\n%%REPLACE_NODE_DATA\n\end{filecontents*}\n\\begin{filecontents*}{elementdata.dat}\n%%REPLACE_ELEMENT_DATA\n\end{filecontents*}\n\n\definecolor{CustomColor}{cmyk}{" + str(cmyk_colors_fill[0]) + "," + str(cmyk_colors_fill[1]) + "," + str(cmyk_colors_fill[2]) + "," + str(cmyk_colors_fill[3]) + "}\n\n\\newcommand\coord[2][]{%\n  \edef\comparenode{#2}%\n  \\foreachitem\zzz\in\\noddat[]{%\n    \edef\\testnode{\\noddat[\zzzcnt,1]}%\n    \ifx\\testnode\comparenode\n      \\xaddtomacro\\tmp{(\\noddat[\zzzcnt,2]#1,\\noddat[\zzzcnt,3]#1)}\\fi\n  }%\n}\n\n\makeatletter\let\\addtomacro\g@addto@macro\makeatother\n\\newcommand\\xaddtomacro[2]{%\n  \edef\\xtmp{#2}%\n  \expandafter\\addtomacro\expandafter#1\expandafter{\\xtmp}%\n}\n\n\\newcommand\drawmesh[1][{\\filldraw[fill=CustomColor,fill opacity=0.1]}]{%\n  \def\\tmp{}%\n  \\foreachitem\z\in\eledat[]{%\n    \\addtomacro\\tmp{#1}%\n    \\foreachitem\zz\in\eledat[\zcnt]{%\n      \ifnum\zzcnt=1\\relax\else\n        \ifnum\zzcnt<\listlen\eledat[\zcnt]\\relax\n          \ifnum\zzcnt=2\\relax\coord{\zz}\\fi\n          \\addtomacro\\tmp{--}%\n          \coord{\eledat[\zcnt,\\the\\numexpr\zzcnt+1\\relax]}%\n        \else\n          \\addtomacro\\tmp{--}%\n          \coord{\eledat[\zcnt,2]}%\n        \\fi\n      \\fi\n    }%\n    \\addtomacro\\tmp{;}%\n  }%\n  \\tmp%\n}\n\n\\newcommand\labelnodes[1][\\node at]{%\n  \\foreachitem\z\in\\noddat[]{%\n    #1 (\\noddat[\zcnt,2],\\noddat[\zcnt,3]){%\n%% ALTERNATIVE 1\n%      \\textcolor{red}{$\\noddat[\zcnt,1]$}};\n%% ALTERNATIVE 2\n      \\fboxsep=0pt\\relax\n      \colorbox{white}{\color{red}$\\noddat[\zcnt,1]$}};\n%%\n  }%\n}\n\n\\newcommand\labelelements[1][\\node at]{%\n  \\foreachitem\z\in\eledat[]{%\n    \def\\tmp{#1 }%\n    \\addtomacro\\tmp{($}\n    \\foreachitem\zz\in\eledat[\zcnt]{%\n      \ifnum\zzcnt=1\\relax\else\n        \ifnum\zzcnt=2\\relax\else\\addtomacro\\tmp{ + }\\fi%\n        \coord[{/\\the\\numexpr\listlen\eledat[\zcnt]-1\\relax}]{%\n          \eledat[\zcnt,\zzcnt]}%\n      \\fi\n    }%\n    \\addtomacro\\tmp{$)}%\n    \\xaddtomacro\\tmp{{\\noexpand\\textcolor{blue!70!green}{$\eledat[\zcnt,1]$}};}%\n   \\tmp\n  }%\n}\n\n\\newcommand\\readmesh[2]{%\n  \ignoreemptyitems%\n  \\readarraysepchar{,}%\n  \ifx\\relax#1\\relax\else\\readdef{#1}\\nodedata\\fi\n  \ifx\\relax#2\\relax\else\\readdef{#2}\elementdata\\fi\n  \setsepchar{,/ }%\n  \\readlist*\\noddat{\\nodedata}%\n  \\readlist*\eledat{\elementdata}%\n}\n\\begin{document}\n%% FILE INPUT\n\\readmesh{nodedata.dat}{elementdata.dat}"

    if landscape_mode:
        template_string += "\n\\begin{landscape}"
    template_string += "\n\n\\begin{figure}[ht]\n\centering\n\\begin{tikzpicture}[scale=1.5]\n  \drawmesh"
    if include_nodes_label:
        template_string += "\n  \labelnodes"
    if include_elements_label:
        template_string += "\n  \labelelements"

    # Calculate limits of the graphic with the model part
    min_x = 0.0
    max_x = 0.0
    min_y = 0.0
    max_y = 0.0
    for node in model_part.Nodes:
        node_x = node.X
        node_y = node.Y
        if min_x > node_x:
            min_x = node_x
        if max_x < node_x:
            max_x = node_x
        if min_y > node_y:
            min_y = node_y
        if max_y < node_y:
            max_y = node_y

    separation_x = max_x/4
    separation_y = max_y/4
    minor_speration = min(separation_x, separation_y)

    # Increasin for a 30% the size of the mesh
    min_x *= 1.3
    max_x *= 1.3
    min_y *= 1.3
    max_y *= 1.3

    if include_axis:
        template_string += "\n  \draw[thick,->] (0,0) -- (0," + str(max_y) + ");\n  \draw[thick,->] (0,0) -- (" + str(max_x) + ",0) node[anchor=north west] {x axis};\n  \draw[thick,->] (0,0) -- (0," + str(max_y) + ") node[anchor=south east] {y axis};\n  \\foreach \\x in {0," + str(separation_x) + "," + str(2 * separation_x) + "," + str(3 * separation_x) + "," + str(4 * separation_x) + "}\n  \draw (\\x cm,1pt) -- (\\x cm,-1pt) node[anchor=north] {$\\x$};\n  \\foreach \y in {0," + str(separation_y) + "," + str(2 * separation_y) + "," + str(3 * separation_y) + "," + str(4 * separation_y) + "}\n  \draw (1pt,\y cm) -- (-1pt,\y cm) node[anchor=east] {$\y$};"
    if include_grid:
        template_string += "\n  \draw[step=" + str(minor_speration) + "cm,gray,very thin] (" + str(min_x) + "," + str(min_y) + ") grid (" + str(max_x) + "," + str(max_y) + ");"
    template_string += "\n\end{tikzpicture}"
    if include_caption:
        template_string += "\n\caption{Finite element mesh}"
    template_string += "\n\end{figure}"
    if landscape_mode:
        template_string += "\n\\end{landscape}"
    template_string += "\n\end{document}"
    return template_string

def GetPrettyTime(time):
    pretty_time = "{0:.12g}".format(time)
    pretty_time = float(pretty_time)
    return pretty_time
