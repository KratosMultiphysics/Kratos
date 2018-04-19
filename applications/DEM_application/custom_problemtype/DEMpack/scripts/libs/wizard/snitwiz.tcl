# snitwiz.tcl --
 #
 #      This file implements package snitwiz, which is a wizard
 #       building package built with snit
 #
 # Copyright (c) 2004 Steve Cassidy
 # 
 #

 package require snit
 #package require msgcat

 #if {[package present msgcat] != "package msgcat not present"} {namespace import msgcat::mc}

 # load message catalogs
 #::msgcat::mcload [file join [file dirname [info script]] .. msgs]

 package provide snitwiz 0.1 


 namespace eval ::snitwiz {
 #namespace export -clear *

     # create a default image
     # image create photo snitwiz::feather 


 } # end of namespace snitwiz


 snit::widget ::snitwiz::wizard {

     hulltype frame

     variable buttonFrame
     variable layoutFrame
     variable sep
     variable helpButton
     variable backButton
     variable nextButton
     variable cancelButton
     variable finishButton

     variable currentstep
     variable disabledSteps {}; # a list of steps that won't be shown


     option -title {Wizard}
     option -showhelp 1

     option -steps {}  ;# the order of steps

     # when we get new steps we make sure each knows it's parent
     # and set the current step to the first    
     onconfigure -steps {value} {
         set options(-steps) $value
             foreach s $value {
                 $s configure -parent $self
             }
             set currentstep 0
     }

     constructor {args} {

         # The dialog is composed of two areas: the row of buttons and the
         # area with the dynamic content. To make it look the way we want it to
         # we'll use another frame for a visual separator
         install buttonFrame using frame $win.buttonFrame -bd 0
         install layoutFrame using frame $win.layoutFrame -bd 0\
             -height 100 \
             -width 300 
         install sep using frame $win.sep1 -class WizSeparator \
             -height 2 -bd 2 -relief groove

         pack $layoutFrame -side top -expand 1 -fill both
         pack $sep -side top -expand 0 -fill x
         pack $buttonFrame -side top -expand 0 -fill x

         # make all of the buttons
         install helpButton using ttk::button $buttonFrame.helpButton \
             -text [= "Help"] \
             -default normal \
             -command [mymethod help]
        tooltip::tooltip $buttonFrame.helpButton [= "Show the help window"]
 
         install backButton  using ttk::button $buttonFrame.backButton \
             -text "< [= Back]" \
             -default normal \
             -width 8 \
             -command [mymethod back]
        tooltip::tooltip $buttonFrame.backButton [= "Go to the previous step"]

         install nextButton  using ttk::button $buttonFrame.nextButton \
             -text "[= Next] >" \
             -default normal \
             -width 8 \
             -command [mymethod next]
        tooltip::tooltip $buttonFrame.nextButton [= "Go to the next step"]

         install finishButton using ttk::button $buttonFrame.finishButton \
             -text [= Exit] \
             -default normal \
             -width 8 \
             -command [mymethod finish]
        
        tooltip::tooltip $buttonFrame.finishButton [= "Exit the wizard"]

         install cancelButton using ttk::button $buttonFrame.cancelButton \
             -text [= Cancel]   \
             -default normal \
             -width 8 \
             -command [mymethod cancel]
        tooltip::tooltip $buttonFrame.finishButton [= "Cancel all defined properties"]

         # pack the buttons
         if {[$self cget -showhelp]} {
            # pack $helpButton -side left -padx 4 -pady 8
         }
         #pack $cancelButton -side right -padx 4 -pady 8
         pack $finishButton -side right -pady 6 -padx 4
         pack $nextButton -side right -pady 6
         pack $backButton -side right -pady 6

        pack $self -side top -expand 1 -fill both

        # $self configure minsize 550 430
        # $self configure maxsize 860 780
     
        $self configurelist $args      
     }

    destructor {
      if {[winfo exists [winfo parent $self]]} {
        #::WarnWin "si exists [winfo parent $self]"
        destroy [winfo parent $self]
      }
     }
    method delete {} {
     
     # ::WarnWin "self:$self"
     if {[winfo exists $self]} {
        destroy $self
     }
     
    }
     method layoutFrame {} {
         return $layoutFrame
     }

     ## methods to handle the button commands
     method help {} {
         # event generate $self <<WizHelp>>
        # ::WinUtils::GiDToolbarState Disable
        # ::WinUtils::GiDWindowState Close
     }

     method next {} {
         event generate $self <<WizNext>>
         update idletasks
        ::kwiz::NextEvent $currentstep
         # go to the next non-disabled step
         foreach step [lrange [$self cget -steps] [expr $currentstep+1] end] {
             incr currentstep
             if {[lsearch $disabledSteps $step] < 0} {
                 $self buildStep $step
                 return
             }
         }
     }

     method back {} {
         event generate $self <<WizBack>>
         update idletasks
        ::kwiz::BackEvent $currentstep
         ## go to the previous non-disabled step
         for {incr currentstep -1} {$currentstep >= 0} {incr currentstep -1} {
             set step [lindex [$self cget -steps] $currentstep]
             if {[lsearch $disabledSteps $step] < 0} {
                 $self buildStep $step
                 return
             }
         }
     }

     method finish {} {
         event generate $self <<WizFinish>>
         if {[winfo exists $self]} {
            destroy $self
         }
     }

     method cancel {} {
         exit
     }

     method buttonstate {name state} {

         switch $name {
             next {$nextButton configure -state $state}
             back {$backButton configure -state $state}
             finish {$finishButton configure -state $state}
             cancel {$cancelButton configure -state $state}
             help {$helpButton configure -state $state}
         }
     }

     ## defining and rendering steps 
     ## a step is an object of type wizardstep
     ## it is rendered by calling it's 'render' method
     ## 

     method buildStep {step} {
         # destroy the existing layout
         eval destroy [winfo children $layoutFrame]

         # run the step rendering method
         $step render $layoutFrame $self

         # enable/disable buttons as required
         if {$currentstep >= [llength [$self cget -steps]]-1} {
             # disable next
             $self buttonstate next disabled
         } else {
             $self buttonstate next normal
         }
         if {$currentstep == 0} {
             $self buttonstate back disabled
         } else {
             $self buttonstate back normal
         }
     }
    method gotostep { step } {
	
	set i 0
	foreach s [$self cget -steps] {
		if {$s == $step } {
			set currentstep $i
		}
		incr i
	}
	
	$self buildStep $step
    }
            
     ##  start --
     #
     # set the wizard going with the first step
     #
     method start {} {
        set currentstep 0
        $self buildStep [lindex $options(-steps) 0]
     }

     ## enable -- 
     #  
     # enable a step, if previously disabled
     method enable {step} {

         set pos [lsearch $disabledSteps $step]
         if {$pos >= 0} {
             set disabledSteps [lreplace $disabledSteps $pos $pos]
         }
     }

     ## disable -- 
     #  
     # disable a step so that it won't appear in the sequence
     #   
     method disable {step} {

         set pos [lsearch $disabledSteps $step]
         if {$pos < 0} {
             lappend disabledSteps $step
         }
     }
 }

 snit::widget ::snitwiz::wizardLayout-basic {
     option -icon {}
     option -title ""
     option -subtitle ""
     option -pretext ""
     option -posttext ""

     variable titleframe
     variable title
     variable subtitle
     variable icon
     variable sep
     variable pretext
     variable posttext
     variable layoutFrame

     onconfigure -icon {value} {
         $icon configure -image $value
         set options(-icon) $value
     }

     onconfigure -title {value} {
         $title configure -text $value
         set options(-title) $value
     }

     onconfigure -subtitle {value} {
         $subtitle configure -text $value
         set options(-subtitle) $value
     }

     onconfigure -pretext {value} {
         $pretext configure -text $value
         set options(-pretext) $value
     }

     onconfigure -posttext {value} {
         $posttext configure -text $value
         set options(-posttext) $value
     }

     constructor {args} {

         install titleframe using frame $win.tf -bd 4 -relief flat \
             -background white
         install title using label $titleframe.t\
             -background white -width 40\
             -text [$self cget -title]\
             -anchor w -justify left
         install subtitle using label $titleframe.st\
             -height 1\
             -background white\
             -padx 15\
             -width 40 \
             -text [$self cget -subtitle]\
             -anchor w -justify left
         install icon using label $titleframe.icon\
             -borderwidth 0 \
             -image $options(-icon) \
             -background white \
             -anchor c
         install sep using frame $win.sep -class WizSeparator\
             -height 2 -bd 2 -relief groove

         set labelfont [font actual [$title cget -font]]
         $title configure -font [concat $labelfont -weight bold]

         install pretext using label $win.pretext -width 40\
             -text [$self cget -pretext]\
             -anchor w -justify left

         install layoutFrame using frame $win.layout -bd 0 -height 200 
         install posttext using label $win.posttext -width 40\
             -text [$self cget -posttext]\
             -anchor w -justify left

         grid $title $icon -sticky ew
         grid $subtitle ^ -sticky ew

         grid $titleframe -sticky new
         grid $sep -sticky new
         grid $pretext -sticky ew
         grid $layoutFrame -sticky nsew
         grid $posttext -sticky ew

         grid columnconfigure $win 0 -weight 1

         grid rowconfigure $win 0 -weight 0   ;# titleframe
         grid rowconfigure $win 1 -weight 0   ;# sep
         grid rowconfigure $win 2 -weight 0   ;# pretext
         grid rowconfigure $win 3 -weight 1   ;# layoutFrame
         grid rowconfigure $win 4 -weight 0   ;# posttext

         $self configurelist $args
     }

     method layoutFrame {} {
         return $layoutFrame
     }
 }


 snit::type ::snitwiz::wizardstep {

     option -title {}
     option -subtitle {}
     option -icon {}
     option -pretext {}
     option -posttext {}
     option -body {}
     option -parent {}
     option -layout basic

     ## render method calls the body code to generate the frame
     method render {frame wizard} {
         # first make the layout
         set layout [::snitwiz::wizardLayout-$options(-layout) $frame.wiz\
                         -title $options(-title)\
                         -subtitle $options(-subtitle)\
                         -icon $options(-icon)\
                         -pretext $options(-pretext)\
                         -posttext $options(-posttext)]

         pack $layout -side top -fill both -expand 1

         ## setup the special variable win as the container frame for
         ## the content note that we also have $wizard as the name of
         ## the containing wizard
         set win [$layout layoutFrame]

         # now render the step itself
         eval [$self cget -body]
     }
 }



 # require_all_keys --
 #
 #   Disable the Next button in a wizard unless all
 #   keys in the given array have values
 #
 # Arguments:
 #   wiz    -- a wizard widget
 #   arrayvar  -- array name
 #   keys   -- list of keys to check 
 # Results:
 #   Puts a variable trace on the given 
 #   array variables which will enable the next button
 #   only if all keys have non-null values.
 #
 proc ::snitwiz::require_all_keys {wiz arrayvar keys} {
     upvar $arrayvar array

     foreach k $keys {
         trace add variable array($k) write \
             [list ::snitwiz::require_all_keys_tracefn $wiz $keys]
     }
 }


 proc ::snitwiz::require_all_keys_cleanup {wiz arrayvar keys} {
     upvar $arrayvar array

     foreach k $keys {
         trace remove variable array($k) write \
             [list ::snitwiz::require_all_keys_tracefn $wiz $keys]
     }
 }

proc ::snitwiz::disablenext { {what "disabled"} } {
   buttonstate next $what
}

proc ::snitwiz::gotonext { } {
   #next
}


 # helper proc used by require_all_keys
 proc ::snitwiz::require_all_keys_tracefn {wiz keys arrayvar key op} {
     upvar $arrayvar array

     foreach k $keys {
         if {$array($k) eq ""}  {
             $wiz buttonstate next disable
             return
         }
     }
     # reset to normal only if all ok
     $wiz buttonstate next normal
 }