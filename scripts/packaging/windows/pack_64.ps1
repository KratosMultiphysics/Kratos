# This script packs and upload a compiled version of kratos.
# Instructions: Set the proper paths under "Defines" and run.

# Defines
$KRATOS_MAJOR = 5
$KRATOS_MINOR = 3
$KRATOS_PATCH = 0
$KRATOS_ARCHY = 64
$KRATOS_NEWSD = "Beta Release"
$KRATOS_TVERS = "$KRATOS_MAJOR.$KRATOS_MINOR"
$KRATOS_TDATE = "$((Get-Date).AddDays(-1).ToString('dd-MM-yyy'))"

$KRATOS_RLS = "C:\KratosStage\Pack64"			                      # Local folder where the release will be packed
$KRATOS_SRC = "C:\KratosSource"		                        			# Kratos source
$KRATOS_CMP = "C:\KratosStage\Release64"			                  # Kratos install directory
$KRATOS_AUX = "C:\KratosStage\LibsForPack64"                    # Kratos additional libraries.
$GIDINT_SRC = "C:\GiDInterface\kratos.gid"											# GID Interface
$GIDINT_KTS = "kratos.gid\exec"			        										# GiD Kratos Interface install dir
$GIDINT_DEM = "UNUSED"			            												# GiD DEM Interface install dir
$DEPLOY_DIR = "C:\KratosDeploy"               									# Deploy directory (masterdisc).
$env:Path  += ";C:\Program Files\7-Zip" 	                    	# Winrar/winzip binaries
$env:Path  += ";C:\Users\Pooyan\AppData\Local\Atlassian\SourceTree\git_local\bin\" # Git bin

$GIDINT_ARR = @("$GIDINT_KTS")					# List of all the places where kratos has to be copied

# Obtain the Git hash of the current version
$GIT_HASH = git.exe --git-dir=C:\KratosSource\.git rev-parse --short HEAD
$GIT_BRAN = git.exe --git-dir=C:\KratosSource\.git rev-parse --abbrev-ref HEAD
$GIT_NUMB = "$GIT_HASH($GIT_BRAN)"

# Clean old releases in the stage folder to prevent problems
Remove-Item -Force -Recurse "$KRATOS_RLS\*"

# Copy Kratos, Kratos.Gid and additional libraries to the local release directory
xcopy /s/e/h/y "$GIDINT_SRC\*"   "$KRATOS_RLS\kratos.gid\"

for ($i=0; $i -lt $GIDINT_ARR.length; $i++) {
	$GITINT_DIR = $GIDINT_ARR[$i]

	xcopy /s/e/h/y "$KRATOS_CMP"   "$KRATOS_RLS\$GITINT_DIR\Kratos\"
	xcopy /s/e/h/y "$KRATOS_AUX\*" "$KRATOS_RLS\$GITINT_DIR\Kratos\"

	# Remove trash (.lib, .svn)
	get-childitem "$KRATOS_RLS\$GITINT_DIR\" -Force -Include *.svn -Recurse | remove-item -Force -Recurse
	get-childitem "$KRATOS_RLS\$GITINT_DIR\" -Force -Include *.lib -Recurse | remove-item -Force -Recurse
}

# Update the xml file with the correct svn release number
[xml]$kratosxml = Get-content "$KRATOS_RLS\kratos.gid\kratos.xml"

$kratosxml.Infoproblemtype.Program.Version               = "$KRATOS_MAJOR.$KRATOS_MINOR.$KRATOS_PATCH"
# $kratosxml.Infoproblemtype.Program.NewsInVersion.version = "$KRATOS_TVERS"
# $kratosxml.Infoproblemtype.Program.NewsInVersion.date    = "$KRATOS_TDATE"
# $kratosxml.Infoproblemtype.Program.NewsInVersion.'#text' = "$KRATOS_NEWSD"

# The ones in the subproblemtypes
# for ($i=0; $i -lt $GIDINT_ARR.length; $i++) {
#	  $GITINT_DIR = $GIDINT_ARR[$i]
#	  $kratosxml.save("$KRATOS_RLS\$GITINT_DIR\kratos.xml")
# }

# The one in the general problemtype, that has the image route in a different place...
# $kratosxml.Infoproblemtype.Program.ImageFileBrowser = "Kratos/kratos.gid/images/Classic/ImageFileBrowser.gif"
$kratosxml.save("$KRATOS_RLS\kratos.gid\kratos.xml")

# Create the the zip file
7z a -r "$KRATOS_RLS\kratos-$KRATOS_MAJOR.$KRATOS_MINOR.$KRATOS_PATCH-$GIT_NUMB-win-$KRATOS_ARCHY.zip" "$KRATOS_RLS\kratos.gid"

# Rm the old release in the deploy directory (Optional)
# Remove-Item "$DEPLOY_DIR\*.zip"

# Move the current release to the deploy directory
# TODO: remove "echo F" and find the proper flag.
echo f Yes | xcopy "$KRATOS_RLS\kratos-$KRATOS_MAJOR.$KRATOS_MINOR.$KRATOS_PATCH-$GIT_NUMB-win-$KRATOS_ARCHY.zip" "$DEPLOY_DIR\kratos-$KRATOS_MAJOR.$KRATOS_MINOR.$KRATOS_PATCH-$GIT_NUMB-win-$KRATOS_ARCHY.zip"
