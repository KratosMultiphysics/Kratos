namespace eval ::DEM {

}

proc ::DEM::Init { } {
    set Kratos::must_quit 1
    if {[GidUtils::GiveProblemTypeFullname G-DEMPack/kratos.gid] ne ""} {
        GiD_Process Mescape Data Defaults ProblemType G-DEMPack/kratos escape 
    } else {
        W "Dem is not installed at G-DEMPack/kratos.gid"
        W "Please, go to Data -> Problemtype -> Internet retrieve and download it"
        spdAux::deactiveApp DEM
        apps::ClearActiveApp
    }

}
::DEM::Init
