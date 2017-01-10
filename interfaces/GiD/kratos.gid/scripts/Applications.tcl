##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval apps {
    variable activeApp
    variable appList
}

proc apps::Init { } {
    variable activeApp
    variable appList
    set activeApp ""
    set appList [list ]
}

proc apps::setActiveApp {appid} {
    variable activeApp
    variable appList
    
    foreach app $appList {
        if {[$app getName] eq $appid} {
            set activeApp $app
            $app activate
            break
        }
    }
    spdAux::activeApp $appid
}

proc apps::getActiveApp { } {
    variable activeApp;
    return $activeApp
}
proc apps::setActiveAppSoft { appid } {
    variable activeApp
    variable appList
    #W "set active app $appid in $appList"
    foreach app $appList {
        #W [$app getName]
        if {[$app getName] eq $appid} {
            set activeApp $app
            break
        }
    }
}

proc apps::getActiveAppId { } {
    variable activeApp;
    set id ""
    catch {
        set id [$activeApp getName]
    }
    return $id
}

proc apps::getAppById { id } {
    variable appList
    set appR ""
    foreach app $appList {
        if {[$app getName] eq $id} {set appR $app; break}
    }
    return $appR
}

proc apps::NewApp {appid publicname prefix} {
    variable appList
    set ap [App new $appid]
    $ap setPublicName $publicname
    $ap setPrefix $prefix
    lappend appList $ap
    return $ap
}

proc apps::getAllApplicationsName {} {
    variable appList
    
    set appnames [list ]
    foreach app $appList {
        lappend appnames [$app getPublicName]
    }
    return $appnames
}

proc apps::getAllApplicationsID {} {
    variable appList
    
    set appnames [list ]
    foreach app $appList {
        lappend appnames [$app getName]
    }
    return $appnames
}

proc apps::getImgFrom { appName {img "logo" } } {
    return [gid_themes::GetImageModule [getImgPathFrom $appName $img] ""]
}
proc apps::getImgPathFrom { appName {img "logo" } } {
    variable appList
    
    set imagespath ""
    foreach app $appList {
        if {[$app getName] eq $appName} {set imagespath [expr {$img == "logo" ? [$app getIcon] : [$app getImagePath $img] }]; break}
    }
    return $imagespath
    #return [Bitmap::get [file native $imagespath]]
}

proc apps::getMyDir {appName} {
    return [file join $::Kratos::kratos_private(Path) apps $appName]
}

proc apps::getCurrentUniqueName {un} {
    return [ExecuteOnCurrentXML getUniqueName $un]
}
proc apps::getAppUniqueName {appName un} {
    variable appList
    foreach app $appList {
        if {[$app getName] eq $appName} {return [$app executexml getUniqueName $un]}
    }
}

proc apps::ExecuteOnCurrentXML { func args} {
    variable activeApp
    if {$activeApp ne ""} {
        return [$activeApp executexml $func {*}$args]
    }
}
proc apps::ExecuteOnApp {appid func args} {
    set response ""
    catch {
        set app [getAppById $appid]
        set response [$app execute $func {*}$args]   
    }
    return $response
}
proc apps::ExecuteOnCurrentApp {func args} {
    variable activeApp
    set response [ExecuteOnApp [$activeApp getName] $func {*}$args]
    return $response
}
proc apps::LoadAppById {appid} {
    variable appList
    foreach app $appList {
        if {[$app getName] eq $appid} {
            $app activate
            break
        }
    }
}

proc apps::isPublic {appId} {
    set app [getAppById $appId]
    if {$app eq ""} {return 0}
    return [$app isPublic]
}

proc apps::CheckElemState {elem inputid {arg ""} } {
    variable activeApp
    
    return [$activeApp executexml CheckElemState $elem $inputid $arg]
}


# Clase App
catch {App destroy}
oo::class create App {
    variable publicname
    variable name
    variable imagepath
    variable writeModelPartEvent
    variable writeParametersEvent
    variable writeCustomEvent
    variable prefix
    variable release
    
    constructor {n} {
        variable name
        variable publicname
        variable imagepath
        variable writeModelPartEvent
        variable writeParametersEvent
        variable writeCustomEvent
        variable prefix
        variable public
        
        set name $n
        set publicname $n
        set imagepath [file nativename [file join $::Kratos::kratos_private(Path) apps $n images] ]
        set writeModelPartEvent $n
        append writeModelPartEvent "::write"
        append writeModelPartEvent "::writeModelPartEvent"
        set writeParametersEvent $n
        append writeParametersEvent "::write"
        append writeParametersEvent "::writeParametersEvent"
        set writeCustomEvent $n
        append writeCustomEvent "::write"
        append writeCustomEvent "::writeCustomFilesEvent"
        set prefix ""
        set public 0
    }
    
    method activate { } {
        variable name
        set dir [file join $::Kratos::kratos_private(Path) apps $name]
        set fileName [file join $dir start.tcl]
        apps::loadAppFile $fileName
        gid_groups_conds::add_images_dir [file join $dir images]
    }
    
    method getPrefix { } {variable prefix; return $prefix}
    method setPrefix { p } {variable prefix; set prefix $p}
    
    method getPublicName { } {variable publicname; return $publicname}
    method setPublicName { pn } {variable publicname; set publicname $pn}
    
    method getName { } {variable name; return $name}
    
    method getIcon { } {return [my getImagePath logo.png]}
    method getImagePath { imgName } {variable imagepath; return [file nativename [file join $imagepath $imgName] ]}
    
    method getWriteModelPartEvent { } {variable writeModelPartEvent; return $writeModelPartEvent}
    
    method getWriteParametersEvent { } {variable writeParametersEvent; return $writeParametersEvent}
    
    method getWriteCustomEvent { } {variable writeCustomEvent; return $writeCustomEvent}
    
    method executexml { func args } {
        #W "func $func "
        variable name
        set f ${name}::xml::${func}
        $f {*}$args
	}
    method execute { func args } {
        #W "func $func "
        variable name
        set f ${name}::${func}
        $f {*}$args
	}
    
    method setPublic {v} {
        variable public
        set public $v
    }
    method isPublic { } {
        variable public
        return $public
    }
}

proc apps::loadAppFile {fileName} {
	uplevel 2 [list source $fileName]
}

apps::Init
