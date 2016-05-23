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
    W "set active app $appid in $appList"
    foreach app $appList {
        W [$app getName]
        if {[$app getName] eq $appid} {
            set activeApp $app
            $app activate
            break
        }
    }
    W "apps activada $appid"
    spdAux::activeApp $appid
}

proc apps::getActiveApp { } {
    variable activeApp;
    return $activeApp
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
    variable appList;
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

proc apps::getImgFrom { appName } {
    variable appList
    
    set imagespath ""
    foreach app $appList {
        if {[$app getName] eq $appName} {set imagespath [$app getIcon]; break}
    }
    return [gid_themes::GetImageModule $imagespath ""]
    #return [Bitmap::get [file native $imagespath]]
}

proc apps::getMyDir {appName} {
    return [file join $::Kratos::kratos_private(Path) apps $appName]
}

proc apps::getCurrentUniqueName {un} {
    return [ExecuteOnCurrent $un getUniqueName]
}
proc apps::getAppUniqueName {appName un} {
    variable appList
    foreach app $appList {
        if {[$app getName] eq $appName} {return [$app executexml getUniqueName $un]}
    }
}

proc apps::ExecuteOnCurrent {elem func} {
    variable activeApp
    if {$activeApp ne ""} {
        return [$activeApp executexml $func $elem]
    }
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

proc apps::isRelease {appId} {
    set app [getAppById $appId]
    if {$app eq ""} {return 0}
    return [$app isRelease]
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
        variable release
        
        set name $n
        set publicname $n
        set imagepath [file nativename [file join $::Kratos::kratos_private(Path) apps $n images logo.gif] ]
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
        set release 0
    }
    
    method activate { } {
        variable name
        set dir [file join $::Kratos::kratos_private(Path) apps $name]
        set fileName [file join $dir start.tcl]
        apps::loadAppFile $fileName
        set func $name
        append func "::LoadMyFiles"
        eval $func
    }
    
    method getPrefix { } {variable prefix; return $prefix}
    method setPrefix { p } {variable prefix; set prefix $p}
    
    method getPublicName { } {variable publicname; return $publicname}
    method setPublicName { pn } {variable publicname; set publicname $pn}
    
    method getName { } {variable name; return $name}
    
    method getIcon { } {variable imagepath; return $imagepath}
    
    method getWriteModelPartEvent { } {variable writeModelPartEvent; return $writeModelPartEvent}
    
    method getWriteParametersEvent { } {variable writeParametersEvent; return $writeParametersEvent}
    
    method getWriteCustomEvent { } {variable writeCustomEvent; return $writeCustomEvent}
    
    method executexml { func args } {
        #W "func $func "
        variable name
        set f ${name}::xml::${func}
        $f {*}$args
	}
    
    method setRelease {v} {
        variable release
        set release $v
    }
    method isRelease { } {
        variable release
        return $release
    }
}

proc apps::loadAppFile {fileName} {

	uplevel 2 [list source $fileName]
}

apps::Init
