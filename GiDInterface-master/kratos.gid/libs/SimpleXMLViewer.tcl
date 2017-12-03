#########################################################
################### Javi Garate #########################
###################### CIMNE ############################
#########################################################


proc Kratos::ViewXML {root} {

    package require BWidget
    package require tdom
    
    set w .gid.win_example
    InitWindow $w [= "A simple XML Viewer"] ExampleCMAS "" "" 1
    
    Tree $w.trjg -yscrollcommand "$w.yrjg set" -height 100
    scrollbar $w.yrjg -ori vert -command "$w.trjg yview"
    grid $w.yrjg -column 1 -row 0 -sticky nse
    grid $w.trjg -column 0 -row 0 -sticky nsew
    
    
    grid rowconfigure $w 1 -weight 1 -minsize 100
    grid columnconfigure $w 0 -weight 1 -minsize 500
    
    after 5 Kratos::recurseInsert $w.trjg $root root
 
}
 proc Kratos::recurseInsert {w node parent} {
    set name [$node nodeName]
    
        set text <$name
        foreach att [$node attributes] {
            catch {append text " $att=\"[$node getAttribute $att]\""}
        }
        append text >
        set fill blue
    
    $w insert end $parent $node -text $text -fill $fill -open 1
    foreach child [$node childNodes] {recurseInsert $w $child $node}
 }