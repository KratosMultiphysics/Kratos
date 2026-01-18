#!/usr/bin/env perl

###############################################################################
#
#   Script: genfort.pl
#
###############################################################################
#
#   A script that converts <libmmg*.h> file into a Fortran include file.
#
#   Usage:
#    > ./genfort.pl -f libmmg*.h
#    >    converts the libmmg*.h into a fortran header.
#    >     -h          Shows this help
#    >     -r <size>   Defines REAL kind (4 or 8 usually)
#    >     -s <size>   Defines MMG_INT  kind (4 or 8 usually)
#
#   Authors:
#     Xavier Lacoste - lacoste@labri.fr
#     Algiane Froehly - algiane@gmail.com
###############################################################################
use POSIX;
use strict;
use Getopt::Std;

###############################################################################
# Group: Variables

#
#   integer: real
#     Kind of the reals
#
#   integer: int
#     Kind of the MMG_INT
#
#   string: fichier
#     Path to the <libmmg*.h> file
#
#   string: format
#     Format to print PARAMETERs.
#
#   hash: opts
#     Options given to the script
my $real     = 0;
my $int      = 0;
my $fichier;
#my $format = "MMG_INTEGER, PARAMETER :: %-30s = %d";
my $format = "#define %-30s %d";
my $formatbyval = "#define %-30s \%val(%d)";
my $definebyval = "#define MMG5_ARG_%-30s \%val(%d)\n";
my %opts;

###############################################################################
# Group: Functions

#
# Function: Usage
#
# Prints usage.
#
sub Usage {

    print "./genfort.pl -f libmmg*.h\n";
    print "  converts the file libmmg*.h into a fortran header.\n";
    print "   -h          Shows this help\n";
    print "   -r <size>   Defines COEF/REAL kind (4 or 8 usually)\n";
    print "   -s <size>   Defines INT kind (4 or 8 usually)\n";
}
#
# Function: printTab
#
# Print *chaine* with *tabcount* indentations.
#
# If *comm* is 0, it will also replace INT and REAL by the
# Correct value.
#
# Parameters:
#   chaine   - String to print
#   tabcount - Number of indentations to add.
#   comm     - Indicate if we are in a comments section.
#
sub printTab # ($chaine, $tabcount, $comm)
{
    my ($chaine, $tabcount, $comm) = @_;
    for (my $i = 0; $i < $tabcount; $i++)
    {
        $chaine = sprintf("%s",$chaine);
    }
    if ($comm == 0)
    {
        if ($int != 0)
        {
            $chaine =~ s/MMG5_INTEGER,/INTEGER(KIND=$int),/g;
        }
        else
        {
            $chaine =~ s/MMG5_INTEGER,/INTEGER,/g;
        }
        if ($real != 0)
        {
            $chaine =~ s/REAL,/REAL(KIND=$real),   /g;
        }
        else
        {
            $chaine =~ s/REAL,/REAL,   /g;
        }
        if ($chaine =~ /\/\*(.*)/)
        {
            # Fortran comment
            $chaine =~ s/\/\*/\! \/\*/g;
        }
    }
    print $chaine;
}

#
# Function: Convert
#
# Main function.
#
# Converts the header <libmmg*.h> file into a Fortran include.
#
sub Convert {

    my $startcom  = 0;
    my $opendef   = 0;
    my $startenum = 0;
    my $countenum = 0;
    my $byvalue   = 0;
    my $chaine;
    my $tabcount = 0;
    my $interfaceprinted = 0;
    my $modulename;

    open (APIc, $fichier);

    foreach my $line ( <APIc> )
    {
        if ($startcom == 0)
        {
            if ($startenum == 0)
            {
                if ($line =~ /^\/\*/)
                {
                    # We are in a comment area
                    if ($line =~ /^[^#]/)
                    {
                        $startcom = 1;
                        $chaine = sprintf("! %s", $line);
                        if ($line =~ /\*\//) {
                            #    remove the "*/" pattern
                            #    $chaine =~ s/\*\///;
                            $startcom = 0;
                        }
                        printTab( $chaine, $tabcount, 1);
                    }
                }
                elsif($line =~ /^\s*enum MMG5_arg/)
                {
                    # We start an enum and we want fortran to pass directly its
                    # value and not its reference.
                    $startenum = 1;
                    $countenum = 0;
                    $byvalue = 1;
                    #$countenum = $countenum + 1 if (/PARAM/);
                }
                elsif($line =~ /^\s*enum/)
                {
                    # We start an enum
                    $startenum = 1;
                    $countenum = 0;
                    $byvalue = 0;
                    #$countenum = $countenum + 1 if (/PARAM/);
                }
                elsif ( $line =~ /^\s*$/ )
                {
                    # Discard line and replace it by a white line
                    print "\n";
                }
                elsif($line =~ /\#ifndef/ )
                {
                    printTab($line,0,0 );
                }
                elsif ($line =~ /\#ifdef USE_POINTMAP/ ) {
                    $opendef++;
                    $chaine = sprintf("! %s",$line);
                    printTab($chaine,0,0 );
                }
                elsif ($line =~ /\#ifdef __cplusplus/ ) {
                    $opendef++;
                    $chaine = sprintf("! %s",$line);
                    printTab($chaine,0,0 );
                }
                elsif($line =~ /\#endif/ )
                {
                    if ( $opendef ) {
                        $opendef--;
                        $chaine = sprintf("! %s",$line);
                        printTab($chaine,0,0 );
                    }
                    else {
                        printTab($line,0,0 );
                    }
                }
                elsif ($line =~ /\#define MMG5_ARG_(\w*)\s+(.*)/)
                {
                    $chaine = sprintf($definebyval,$1,$2);
                    printTab($chaine,1,0 );
                }
                elsif ($line =~ /\#define/)
                {
                    printTab($line,1,0 );
                }
                elsif($line =~ /typedef/)
                {
                    if ( $line =~ /{/ ) {
                        # We start a typedef area, we want to comment it
                        while (<APIc>)
                        {
                            if (/[^}]/)
                            {
                                $chaine = sprintf("! %s", <APIc>);
                                printTab( $chaine, 1, 1);
                                redo unless eof();
                            }
                        }
                        $chaine = sprintf("! %s", $line);
                        printTab( $chaine, 1, 1);
                    }
                    else
                    {
                        $chaine = sprintf("! %s", $line);
                        printTab( $chaine, 1, 1);
                    }
                }
                else
                {
                    $chaine = sprintf("! %s", $line);
                    printTab( $chaine, 1, 1);
                }
            }
            else
            {
                if ($line =~ /}/)
                {
                    $startenum = 0;
                    $countenum = 0;
                }
                elsif($line =~ /[ ]*{$/)
                {
                    # bracket line, do nothing
                }
                elsif($line =~ /[ ]*([^ |^,]*)[ ]*,?/)
                {
                    my $key = $line;
                    chomp $key;
                    $key =~ s/,//g;

                    if ( $key =~ /(.*)(\/\*.*\*\/)/ )
                    {
                        if ( $byvalue ) {
                            # Fortran code must pass a value and not a
                            # reference over a value
                            $chaine = sprintf($formatbyval, $1,$countenum);
                        }
                        else {
                            $chaine = sprintf($format, $1, $countenum);
                        }
                        $chaine = "$2\n$chaine\n";
                    }
                    else
                    {
                        if ( $byvalue ) {
                            # Fortran code must pass a value and not a
                            # reference over a value
                            $chaine = sprintf($formatbyval, $key, $countenum,);
                        }
                        else {
                            $chaine = sprintf($format, $key, $countenum);
                        }
                        $chaine = "$chaine\n";
                    }
                    printTab($chaine,$tabcount, 0);
                    $countenum++;
                }


            }
        }
        else
        {
            if ($line =~ /^[ \*]*> (.*)\\n/ )
            {
                if ($interfaceprinted == 0)
                {
                    if ( $startcom == 1 ) {
                        $chaine = "!  */\nINTERFACE\n";
                    }
                    else {
                        $chaine = "INTERFACE\n";
                    }
                    printTab($chaine, $tabcount);
                    $tabcount = $tabcount+1;
                    $interfaceprinted = 1;
                }

                $chaine = sprintf("%s\n", $1);
                if ( $interfaceprinted == 1 && $1=~/.*END SUBROUTINE/  )
                {
                    $chaine = sprintf("%s%s\n",$chaine,"END INTERFACE");
                    $interfaceprinted = 0;
                }

                printTab($chaine, $tabcount, 0);
            }
            elsif ($line =~ /(.*)\*\//)
            {
                $startcom = 0;
                $chaine = sprintf("! %s\n", $line);
                printTab($chaine, $tabcount, 1);
            }
            elsif($line =~ /^\ \*/ )
            {
                $chaine = sprintf("! %s", $line);
                printTab($chaine, $tabcount, 1);
            }
            elsif($line =~ /^\*\*/ )
            {
                $chaine = sprintf("! %s", $line);
                printTab($chaine, $tabcount, 1);
            }
            elsif($line =~ /(.*)/)
            {
                $chaine = sprintf("! %s\n", $line);
                if ($line =~ /Mmg's constants/)
                {
                    my $chaine2 .= "!   Enum: KINDS\n";
                    $chaine2 .= "!\n";
                    $chaine2 .= "!   Type kinds\n";
                    $chaine2 .= "!\n";
                    $chaine2 .= "!   Contains:\n";
                    $chaine2 .= "!     MMG5_DATA_PTR_T - Kind to use for POINTERS\n";
                   # $chaine2 .= "!     MMG5_REAL_KIND  - Kind to use for REAL\n";
                   # $chaine2 .= "!     MMG5_INT_KIND   - Kind to use for INT\n";
                    $chaine2 .= "!\n";
                    if ($real != 0)
                    {
                        $chaine2 .= "INTEGER, PARAMETER :: MMG5_REAL_KIND                = $real\n";
                    }
                    if ($int != 0)
                    {
                        $chaine2 .= "INTEGER, PARAMETER :: MMG5_INT_KIND                = $int\n";
                    }
                    $tabcount --;
                    printTab($chaine2, $tabcount, 0);
                }

                printTab($chaine, $tabcount, 1);
            }
        }


    }

    close APIc;
}


getopts("hf:r:i:",\%opts);

if ( defined $opts{r} ){
    $real = $opts{r};
}
if ( defined $opts{i} ){
    $int = $opts{i};
}

if ( defined $opts{f} ){
    $fichier = $opts{f};
}
else {
    Usage();
    exit;
}

if ( defined $opts{h} ){
    Usage();
    exit;
}

Convert();
