

proc WriteMaterials { basename dir problemtypedir TableDict} {

    ## Start ProjectParameters.json file
    set filename [file join $dir materials.py]
    set FileVar [open $filename w]

    puts $FileVar ""
    puts $FileVar "from __future__ import print_function, absolute_import, division"
    puts $FileVar "from KratosMultiphysics import *"
    puts $FileVar "from KratosMultiphysics.SolidMechanicsApplication import *"
    puts $FileVar "from KratosMultiphysics.FemToDemApplication import *"
    puts $FileVar ""
    puts $FileVar "def AssignMaterial(Properties):"
    puts $FileVar ""

    set Groups [GiD_Info conditions Body_Part groups]
    #W "Grupos  [lindex [lindex $Groups 1] 3]"
    #puts $FileVar "testeo  [lindex [lindex $Groups 1] 3]"

    for {set i 0} {$i < [llength $Groups]} {incr i} {

        incr PropertyId
        puts $FileVar "    prop_id = $PropertyId"
        puts $FileVar "    prop = Properties\[prop_id\]"
        puts $FileVar "    mat = LinearElastic3DLaw()"
        puts $FileVar "    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())"
        puts $FileVar ""

    }





    close $FileVar

}