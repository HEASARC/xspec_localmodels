# Support routines common to more than one cluster mass model
#
############################################################
#
# Number of beta model parameters (same for all models)
#
proc fixedPars {} {
    return 4
}
#
############################################################
#
# Mixing model name
#
proc mixName {} {
    tclout compinfo 1
    return [lindex $xspec_tclout 0]
}
#
############################################################
#
# Number of mixing model parameters
#
proc mixNpars {} {
    tclout compinfo 1
    return [lindex $xspec_tclout 2]
}
#
############################################################
#
# Maximum number of shells
#
proc maxShells {} {
    set mixname [mixName]
    return [${mixname}MaxShells [mixNpars]]
}
#
############################################################
#
# Front end sanity checks and numbers for accessing mass mixing
# models.
# Checks for:
# a known mixing component
# parameter names beta, switch, kT*
# Sets values:
# mxName to the maximum number of shells
# shName to the number of active shells
# mxmName to the name of the mixing model
#
proc frontend {mxName shName mxmName} {
    upvar $mxName maxshells $shName shells $mxmName mixname
    # Mixing component must be a known one
    set mixname [mixName]
    if {![${mixname}Known]} {
	puts stderr "frontend: expecting a known mass model, got $mixname"
	return 1
    }

    # Total number of mixing model parameters
    set mixpars [mixNpars]
    # Number of beta model parameters
    set fixedpars [fixedPars]

    # Check names of last two fixed params and the first temperature
    tclout pinfo [expr $fixedpars - 1]
    set betname [lindex $xspec_tclout 0]
    tclout pinfo $fixedpars
    set swname [lindex $xspec_tclout 0]
    tclout pinfo [expr $fixedpars + 1]
    set ktfirst [lindex $xspec_tclout 0]
    if {[string compare $betname beta] != 0 \
	    || [string compare $swname switch] != 0 \
	    || ![string match {kT*} $ktfirst]} {
	puts stderr "frontend: unexpected parameter names for $mixname"
	return 1
    }

    # Maximum number of shells
    set maxshells [maxShells]
    # Get number of active shells
    tclout datagrp
    set shells $xspec_tclout
    if {$shells > $maxshells} {
	puts stderr "frontend: too many shells - edit model.dat to increase"
	return 1
    }
    return 0
}
#
############################################################
#
# Find a parameter by name and the name of its units.
# Search can be limited to a single model component.
# Returns component and parameter number of first found.
#
proc findParName {parName unitName {compNum 0}} {
    # Find the component with the named parameter
    if {$compNum} {
	set compStart $compNum
	set compEnd $compNum
    } else {
	tclout modcomp
	set compStart 1
	set compEnd $xspec_tclout
    }
    for {set icomp $compStart} {$icomp <= $compEnd} {incr icomp} {
	tclout compinfo $icomp
	set pstart [lindex $xspec_tclout 1]
	set pend [expr $pstart + [lindex $xspec_tclout 2] - 1]
	for {set ipar $pstart} {$ipar <= $pend} {incr ipar} {
	    tclout pinfo $ipar
	    if {[string compare [lindex $xspec_tclout 0] $parName] == 0
	        && [string compare [lindex $xspec_tclout 1] $unitName] == 0} {
		return [list $icomp [expr $ipar - $pstart]]
	    }
	}
    }
    return {}
}
#
############################################################
#
# Parameter linking and freezing required by clmass, monomass and
# nfwmass.
# Arguments give settings for the beta model (defaults should
# match model.dat).
#
proc massmod_start {{rinner 0} {a 0} {beta 0.7} {sw 0}} {
    # Sanity checks, etc
    if {[frontend maxshells shells mixname]} {
	return
    }

    # Identify thermal component
    set tmod [findParName kT keV]
    if {[llength $tmod] != 2} {
	puts stderr "massmod_start: no thermal component"
	return
    }
    set thermcomp [lindex $tmod 0]
    set thermoff [lindex $tmod 1]
    tclout compinfo $thermcomp
    set thermname [lindex $xspec_tclout 0]

    # Link thermal model temperatures to clmass temperatures
    set fixedpars [fixedPars]
    for {set igrp 1} {$igrp <= $shells} {incr igrp} {
	tclout compinfo $thermcomp $igrp
	if {[string compare [lindex $xspec_tclout 0] $thermname] != 0} {
	    puts stderr "massmod_start failed redundant check"
	    return
	}
	set thpar [expr [lindex $xspec_tclout 1] + $thermoff]
	set mixtpar [expr $fixedpars + $igrp]
	newpar $thpar =$mixtpar
    }

    # Freeze unused shell parameters
    if {$shells < $maxshells} {
	set ftstart [expr $fixedpars + $shells + 1]
	set ftend [expr $fixedpars + $maxshells]
	if {[${mixname}HasDens]} {
	    set fdstart [expr $ftstart + $maxshells]
	    set fdend [expr $ftend + $maxshells]
	    freeze $ftstart-$ftend $fdstart-$fdend
	} else {
	    freeze $ftstart-$ftend
	}
    }

    # Set and freeze beta model params 
    newpar 1 $rinner
    newpar 2 $a
    newpar 3 $beta
    newpar 4 $sw
    freeze 2 3
    # Freeze density for the outermost shell when using the beta model
    if {$sw != 0 && [${mixname}HasDens]} {
	set extraden [expr $fixedpars + $maxshells + $shells]
	freeze $extraden
    }
}
#
############################################################
#
# Untie model parameters between data groups.
# $pars is a list of parameter name, unit name pairs.
# The first model component found to contain the first named parameter,
# unit pair is the only component affected.  All other params must belong
# to the same model component.
#
proc breakLinks {{pars {{kT keV} {norm ""}}}} {
    # Find the component and offset of the first
    set p [lindex $pars 0]
    set primePar [findParName [lindex $p 0] [lindex $p 1]]
    if {[llength $primePar] != 2} {
	puts stderr "breakLinks: failed to find [lindex $p 0] in model"
	return
    }
    # Model component containing the primary parameter
    set comp [lindex $primePar 0]
    set unlink [lindex $primePar 1]
    # Offsets of the remaining parameters
    for {set i 1} {$i < [llength $pars]} {incr i} {
	set p [lindex $pars $i]
	set pName [lindex $p 0]
	set pUnit [lindex $p 1]
	set pLoc [findParName $pName $pUnit $comp]
	if {[llength $pLoc] != 2} {
	    puts stderr "breakLinks: no $pName in model component $comp"
	    return
	}
	lappend unlink [lindex $pLoc 1]
    }
    # Untie all from the first
    tclout datagrp
    set ndg $xspec_tclout
    set cmd untie
    for {set igrp 2} {$igrp <= $ndg} {incr igrp} {
	tclout compinfo $comp $igrp
	set pbase [lindex $xspec_tclout 1]
	for {set j 0} {$j < [llength $unlink]} {incr j} {
	    lappend cmd [expr $pbase + [lindex $unlink $j]]
	}
    }
    eval $cmd
    return
}
#
############################################################
#
# If a model parameter is variable, append it to the list
#
proc psave {saveName pnum} {
    tclout plink $pnum
    # Ignore linked parameters
    if {[string compare [lindex $xspec_tclout 0] F] == 0} {
	tclout param $pnum
	# Ignore frozen parameters
	if {[lindex $xspec_tclout 1] > 0} {
	    upvar $saveName saves
	    lappend saves [list $pnum [lindex $xspec_tclout 0]]
	}
    }
    return
}
#
############################################################
#
# Save all parameters that can vary - ie those that are untied
# and thawed.
#
proc parsave {} {
# There should only be one copy of the mixing model, so deal
# with that first.
    tclout compinfo 1
    set mixpars [lindex $xspec_tclout 2]
    for {set ipar 1} {$ipar <= $mixpars} {incr ipar} {
	# Save the parameter if it is variable
 	psave saves $ipar
    }
    # Get number of model componnents
    tclout modcomp
    set ncomp $xspec_tclout
# Save the rest of the variable parameters.
    tclout datagrp
    set ndg $xspec_tclout
    for {set igrp 1} {$igrp <= $ndg} {incr igrp} {
	for {set icomp 2} {$icomp <= $ncomp} {incr icomp} {
	    tclout compinfo $icomp $igrp
	    set pbase [lindex $xspec_tclout 1]
	    set pnum [lindex $xspec_tclout 2]
	    set plast [expr $pbase + $pnum - 1]
	    for {set ipar $pbase} {$ipar <= $plast} {incr ipar} {
		psave saves $ipar
	    }
	}
    }
    return $saves
}
#
############################################################
#
# Restore params from list made by parsave
#
proc restore {savName} {
    upvar $savName saves
    set n [llength $saves]
    tclout chatter
    set oldchatter $xspec_tclout
    chatter 5 5
    for {set i 0} {$i < $n} {incr i} {
	set pair [lindex $saves $i]
	newpar [lindex $pair 0] [lindex $pair 1]
    }
    chatter [lindex $oldchatter 0] [lindex $oldchatter 1]
    return
}
#
############################################################
#
# Get the shell boundaries
#
proc shellBounds {} {
    # Go through spectra to get rout for each shell
    tclout datasets
    set ndata $xspec_tclout
    for {set j 1} {$j <= $ndata} {incr j} {
	tclout datagrp $j
	set dgno $xspec_tclout
	tclout xflt $j
	# Make sure we have XFLT0001 at least
	if {[lindex $xspec_tclout 0] == 0} {
	    puts stderr "shellBounds: no XFLT keywords for dataset $j"
	    return
	}
	set rout [lindex $xspec_tclout 1]
	# Get new entries and check those already seen
	if {[info exists dg($rout)]} {
	    if {$dg($rout) != $dgno} {
		puts stderr "shellBounds: mismatched shells"
		return
	    }
	} else {
	    set dg($rout) $dgno
	}
    }
    # Inner radius is first parameter
    tclout param 1
    set rprev [lindex $xspec_tclout 0]
    # List shell boundaries and group number in radial order
    foreach r [lsort -real [array names dg]] {
	lappend rbounds [list $rprev $r $dg($r)]
	set rprev $r
    }
    return $rbounds
}
#
############################################################
#
# Data group number for the i^{th} shell.
# irad is shell number (counting from 0)
# rbName is the name of the list returned by shellBounds
#
proc dgOfRadial {irad rbName} {
    upvar $rbName rbounds
    return [lindex [lindex $rbounds $irad] 2]
}
#
############################################################
#
# Inner radius of the i^{th} shell.
# irad is shell number (counting from 0)
# rbName is the name of the list returned by shellBounds
#
proc rinOfRadial {irad rbName} {
    upvar $rbName rbounds
    return [lindex [lindex $rbounds $irad] 0]
}
#
############################################################
#
# Outer radius of the i^{th} shell.
# irad is shell number (counting from 0)
# rbName is the name of the list returned by shellBounds
#
proc routOfRadial {irad rbName} {
    upvar $rbName rbounds
    return [lindex [lindex $rbounds $irad] 1]
}
#
############################################################
#
# Get redshift from the first thermal model component
#
proc getz {} {
    # Identify thermal model component
    set tmod [findParName kT keV]
    if {[llength $tmod] != 2} {
	puts stderr "getz: no thermal component"
	return
    }
    set thermcomp [lindex $tmod 0]
    # Parameter number for the redshift
    set zpar [findParName Redshift "" $thermcomp]
    if {[llength $zpar] != 2} {
	puts stderr "getz: no Redshift in component $thermcomp"
	return
    }
    tclout compinfo $thermcomp
    set zpnum [expr [lindex $xspec_tclout 1] + [lindex $zpar 1]]
    tclout param $zpnum
    return [lindex $xspec_tclout 0]
}
#
############################################################
#
# Lists names and numbers of the active density related params
# in radial order (for clmass and monomass).
# rbName is the name at the calling level of the rbounds list as
# returned by shellBounds
#
proc denPars {rbName} {
    upvar $rbName rbounds
    set fixedpars [fixedPars]
    set dpbase [expr $fixedpars + [maxShells]]
    set dpno [llength $rbounds]
    # Density parameter for outermost shell is ignored when
    # the beta model is in use
    tclout param $fixedpars
    if {[lindex $xspec_tclout 0]} {
	incr dpno -1
    }
    for {set i 0} {$i < $dpno} {incr i} {
	set pnum [expr $dpbase + [dgOfRadial $i rbounds]]
	tclout pinfo $pnum
	# Append a list giving parameter name and number
	lappend dpnames [list [lindex $xspec_tclout 0] $pnum]
    }
    return $dpnames
}
#
############################################################
#
# Angular diameter distance and critical density
#
proc cospars {z {cosfile stdcosmo.pars}} {
    set cno [exec ./cosinfo $cosfile $z]
    regexp {^Angular diameter distance: (\S+) m\n} $cno junk angdd
    regexp {\nCritical density: (\S+) kg m\^{-3}} $cno junk rhocrit
    return [list $angdd $rhocrit]
}
#
############################################################
#
# Length, density and mass units for the model in SI
# dang = angular diameter distance (m)
# radialScale = arcsec per distance unit (for shell radii)
#
proc phyScale {dang radialScale} {
    set arcsec [expr $phyconst::pi / (180.0 * 3600.0)]
    set runit [expr $radialScale * $arcsec * $dang]
    # $$ M = M' {r \over r'} {kT \over (kT)'} {1 \over \mu m_H G} $$
    set Munit [expr $runit * $phyconst::keV \
		   / ($phyconst::mu * $phyconst::m_H * $phyconst::G_newton)]
    set denunit [expr $Munit / ($runit * $runit * $runit)]
    return [list $runit $denunit $Munit]
}
#
############################################################
#
# Set physical scales used for fitting.
# radiascale - unit of shell radii in arcsec
# rbName - list of shell boundaries made by shellBounds
# rmName - name for model length unit
# mmName - name for model mass unit
# rcmName - name for critical density given in model units
#
proc setscales {radialScale rbName rmName mmName rcmName {redshift -1}} {
    # Shell boundaries
    upvar $rbName rbounds 
    set rbounds [shellBounds]
    # Get redshift from model parameters if not provided by caller
    if {$redshift < 0.0} {
	set redshift [getz]
    }
    # Angular diameter distance and critical density
    set cpl [cospars $redshift]
    set dang [lindex $cpl 0]
    set dMpc [expr $dang / $phyconst::megaparsec]
    puts "Angular diameter distance: $dMpc Mpc"
    set rhocrit [lindex $cpl 1]
    puts "Critical density: $rhocrit kg m^{-3}"
    # Model scales in physical units
    set phscale [phyScale $dang $radialScale]
    upvar $rmName runit
    set runit [lindex $phscale 0]
    puts "Distance unit: [expr $runit / $phyconst::kiloparsec] kpc"
    set denunit [lindex $phscale 1]
    puts "Model density unit: $denunit kg m^{-3}"
    upvar $mmName Munit
    set Munit [lindex $phscale 2]
    puts "Model mass unit: [expr $Munit / $phyconst::Msun] Msun"
    upvar $rcmName rcrit
    set rcrit [expr $rhocrit / $denunit]
    puts "Critical density in model units: $rcrit"
}
#
############################################################
#
# Set mass model temperatures to a fixed value.
# kTval = gas temperature in keV
# nshell = number of active shells
#
proc setMMkT {kTval nshell} {
    # Convenient to use 1 based count here
    set fixed [fixedPars]
    for {set ishell 1} {$ishell <= $nshell} {incr ishell} {
	newpar [expr $fixed + $ishell] $kTval
    }
}
#
############################################################
#
# Estimate virial radius (Mpc) from temperature using scaling
# of Evrard et al (1996).
# kTval = cluster gas temperature in keV
#
proc rVirEstimate {kTval} {
    return [expr 2.74 * sqrt (0.1 * $kTval)]
}
#
############################################################
#
# Overdensity at the virial radius
#
proc virialOverDensity {} {
    return 200
}
#
############################################################
#
# Use the perfect isothermal sphere to estimate mean mass densities
# for shells in model units.
# rbName = name of list giving radial boundaries
# kTval = gas temperature in keV
# rUnit = model distance unit (m)
# rhoCritModel = critical density in model units
#
proc denIsothermal {rbName kTval rUnit rhoCritModel} {
    upvar $rbName rb
    set nshell [llength $rb]

    # Estimate the virial radius in Mpc
    set rvMpc [rVirEstimate $kTval]
    # Virial radius in model units
    set rVirMod [expr $rvMpc * $phyconst::megaparsec / $rUnit]
    # Density norm in model units, assuming mean density within
    # the virial radius is 200 times the critical density and
    # that the cluster is a perfect isothermal sphere.
    set overden [virialOverDensity]
    set rhoNormMod [expr $overden * $rhoCritModel * $rVirMod * $rVirMod]

    set denList {}
    for {set ishell 0} {$ishell < $nshell} {incr ishell} {
	set rin [rinOfRadial $ishell rb]
	set rout [routOfRadial $ishell rb]
	set rhoMod [expr $rhoNormMod / ($rout * ($rout + $rin) + $rin * $rin)]
	lappend denList $rhoMod
    }
    return $denList
}
#
############################################################
#
# Form of mass profile for NFW potential
#
proc NFWmassForm {x} {
    if {$x < 0.001} {
        return [expr $x * $x * (0.5 - $x * ((2.0 / 3.0) - $x * (0.75 \
		     - $x * (0.8 * $x))))]
    }
    return [expr log (1 + $x) - $x / (1 + $x)]
}
#
############################################################
#
# Use NFW parameters (in model units) to estimate shell mass
# densities.
# rbName = name of list of radial boundaries
# nfwa = NFW length scale in model units
# conc = concentration parameter
# rhoCritModel = critical density in model units
#
proc denNFW {rbName nfwa conc rhoCritModel} {
    upvar $rbName rb
    set nshell [llength $rb]

    set overden [virialOverDensity]
    set rhoNormMod [expr $overden * $rhoCritModel * $conc * $conc * $conc \
			/ [NFWmassForm $conc]]

    set denList {}
    for {set ishell 0} {$ishell < $nshell} {incr ishell} {
	set xin [expr [rinOfRadial $ishell rb] / $nfwa]
	set xout [expr [routOfRadial $ishell rb] / $nfwa]
	set rhoMod [expr $rhoNormMod * \
		    ([NFWmassForm $xout] - [NFWmassForm $xin]) \
		    / (($xout - $xin) * ($xout * ($xout + $xin) + $xin * $xin))]
	lappend denList $rhoMod
    }
    return $denList
}
#
############################################################
#
# Use kT and c to estimate NFW parameters (model units).
# kTval = temperature guess (keV)
# conc = concentratio parameter guess
# rUnit = physical size of model scale unit
# rhoCritModel = critical density in model units
#
proc NFWguess {kTval conc rUnit rhoCritModel} {
    # Virial radius in Mpc
    set rvMpc [rVirEstimate $kTval]
    # Virial radius in model units
    set rVirMod [expr $rvMpc * $phyconst::megaparsec / $rUnit]

    # NFW scale length in model units
    set NFWa [expr $rVirMod / $conc]
    # NFW norm in model units
    set overden [virialOverDensity]
    set NFWnorm [expr (4 * $phyconst::pi / 3) * $overden * $rhoCritModel \
		     * $rVirMod * $rVirMod * $conc / [NFWmassForm $conc]]
    return [list $NFWa $NFWnorm]
}
#
############################################################
#
# Set initial parameters for the mass models, using estimates
# of temperature and concentration parameter for an NFW model.
# rbName = name of list giving radial boundaries
# kTval = gas temperature in keV
# conc = concentration parameter
# rUnit = model distance unit (m)
# rhoCritModel = critical density in model units
#
proc mixNFWstartPars {rbName kTval conc rUnit rhoCritModel} {
    upvar $rbName rb
    set mixmod [mixName]
    return [${mixmod}NFWstartPars rb $kTval $conc $rUnit $rhoCritModel]
}
#
############################################################
#
# Apply fixed mass constraint to the shell contributing
# most to the mass for clmass and monomass.
# mass - fixed mass in model units
# mcName - list of density params and associated volumes
#
proc applyConstraint {mass mcName} {
    upvar $mcName mcdat
    set dmmax 0.0
    set pivot -1
    # Find shell that contributes most to the mass
    foreach pv $mcdat {
	set pnum [lindex $pv 0]
	tclout param $pnum
	set den [lindex $xspec_tclout 0]
	set vol [lindex $pv 1]
	set dm [expr $vol * $den]
	if {$dm > $dmmax} {
	    set dmmax $dm
	    set pivot $pnum
	    set pvol $vol
	}
    }
    # Build and apply constraint string
    set onpvol [expr 1.0 / $pvol]
    set straint [format " = %#.7g" [expr $mass * $onpvol]]
    foreach pv $mcdat {
	set pnum [lindex $pv 0]
	if {$pnum != $pivot} {
	    set vol [lindex $pv 1]
	    set dfac [expr $vol * $onpvol]
	    set straint [format "%s - %#.7g * %d" $straint $dfac $pnum]
	}
    }
    newpar $pivot $straint
    return $pivot
}
#
############################################################
#
# Total mass inside radius rcut for clmass and monomass.
#
proc massOfR {rcut rbName} {
    upvar $rbName rbounds
    tclout compinfo 1
    set mixname [lindex $xspec_tclout 0]
    if {[${mixname}HasDens]} {
	# Density parameters and their related volumes
	set mcdat [${mixname}ConstraintData $rcut rbounds]
	set mass 0.0
	foreach pv $mcdat {
	    set pnum [lindex $pv 0]
	    tclout param $pnum
	    set dden [lindex $xspec_tclout 0]
	    set vol [lindex $pv 1]
	    set mass [expr $mass + $vol * $dden]
	}
    } else {
	set mass [${mixname}MofR $rcut rbounds]
    }
    return $mass
}    
#
############################################################
#
# Fit model with mass inside rfix equal to mfix for clmass
# or monomass.
#
proc fitConstrained {mfix rfix rbName} {
    upvar $rbName rbounds
    tclout compinfo 1
    set mixname [lindex $xspec_tclout 0]
    # Density parameters and their related volumes
    set mcdat [${mixname}ConstraintData $rfix rbounds]
    set oldpivot -1
    tclout chatter
    set oldchatter $xspec_tclout
    chatter 5 5
    while {[set pivot [applyConstraint $mfix mcdat]] != $oldpivot} {
	fit 100
	untie $pivot
	set oldpivot $pivot
    }
    tclout stat
    untie $pivot
    chatter [lindex $oldchatter 0] [lindex $oldchatter 1]
    return $xspec_tclout
}
#
############################################################
#
# Find all solutions where average density is rhobar.
#
proc avdsol {rhobar rbName} {
    upvar $rbName rbounds
    # Get list of shell densities
    tclout compinfo 1
    set mixname [lindex $xspec_tclout 0]
    set dens [${mixname}ShellDensities rbounds]

    # Find all solutions
    set vnorm [expr 4.0 * $phyconst::pi / 3.0]
    set mtot 0.0
    foreach sh $rbounds d $dens {
	if {[string length $d] > 0} {
	    set ri [lindex $sh 0]
	    set numer [expr $mtot - $d * $vnorm * $ri * $ri * $ri]
	    set denom [expr $rhobar - $d]
	    if {$numer * $denom > 0.0} {
		set r [expr pow ($numer / ($vnorm * $denom), 1.0 / 3.0)]
		if {$r >= $ri && $r < [lindex $sh 1]} {
		    if {$rhobar > $d} {
			lappend sols $r
		    } else {
			lappend sols -$r
		    }
		}
	    }
	    set rout [lindex $sh 1]
	    set mtot [expr $mtot + $d * $vnorm * ($rout - $ri) \
			  * ($rout * ($rout + $ri) + $ri * $ri)]
	}
    }
    if {[info exists sols]} {
	return $sols
    }
    return
}
#
############################################################
#
# For finding mass confidence range at fixed radius, gives
# radius from fixed radius and mass.
#
proc fixedRad {rfix mass} {
    return $rfix
}
#
############################################################
#
# For mass confidence range at fixed mean density, gives
# radius from density and mass.
#
proc fixedDen {rhofix mass} {
    return [expr pow (0.75 * $mass / ($phyconst::pi * $rhofix), 1.0 / 3.0)]
}
#
############################################################
#
# Radius and mass from fixed radius at best fit.
#
proc rmSolFixedRad {rfix rbName} {
    upvar $rbName rbounds
    set bestMass [massOfR $rfix rbounds]
    return [list $rfix $bestMass]
}
#
############################################################
#
# Compute radius and mass from the mean density at best fit
#
proc rmSolFixedDen {rhofix rbName} {
    upvar $rbName rbounds
    set sols [avdsol $rhofix rbounds]
    # Deal with possible multiple solutions
    set np 0
    set nm 0
    foreach s $sols {
	if {$s > 0.0} {
	    incr np
	    if {![info exists rbest]} {
		set rbest $s
	    }
	} else {
	    incr nm
	}
    }
    if {![info exists rbest]} {
	puts stderr "rmSolFixedDen: no radii with mean density $rhofix (-$nm)"
	return 
    } elseif {$np > 1 || $nm > 0} {
	puts stderr "rmSolFixedDen: multiple solutions: -$nm +$np"
    }
    set vnorm [expr 4.0 * $phyconst::pi / 3.0]
    set bestMass [expr $vnorm * $rhofix * $rbest * $rbest * $rbest]
    return [list $rbest $bestMass]
}
#
############################################################
#
# Search for overdensity radius at which the best fit statistic
# exceeds a target value.
#
# targetStat - value of fit statistic to bracket
# fixed - fixed radius or mean density
# epsStat - target precision
# min - starting mass
# sin - fit stat corresponding to starting radius
# step - factor by which to step the mass
# rbName - name of list of shell boundaries
# rprocName - routine to compute radius
# bestStat - best fit statistic, if set
# tinyDs - change in fit statistic for slow progress
# maxSlow - maximum number of slow steps
#
proc bracket {targetStat fixed epsStat min sin step rbName rprocName 
	      {bestStat 0} {tinyDs 0} {maxSlow 5}} {
    upvar $rbName rbounds
    if {$tinyDs <= 0.0} {
	# May need to tune this
	set tinyDs [expr 2.0 * $epsStat]
    }
    # Outer loop reduces step size until results seem reliable
    while {1} {
	set slowCount 0
	# Inner loop steps outward
	while {1} {
	    set pin [parsave]
	    set mtry [expr $min * $step]
	    set rtry [$rprocName $fixed $mtry]
	    # Apply constraint and fit
	    set stry [fitConstrained $mtry $rtry rbounds]
	    puts "radius, mass, stat: $rtry $mtry $stry"
	    if {$stry < $bestStat - $epsStat} {
		# New best fit - return with params at the new fit
		return -1
	    }
	    if {$stry > $targetStat + $epsStat} {
		# Apparent bracket found
		break;
	    }

	    if {$stry - $sin < $tinyDs 
		&& $targetStat - $stry > 2.0 * ($stry - $sin)} {
		# Repeated tiny steps when far from the target means 
		# step is too small
		incr slowCount
		if {$slowCount >= $maxSlow} {
		    puts "Increase step"
		    set step [expr $step * $step]
		}
	    } else {
		set slowCount 0
	    }
	    set min $mtry
	    set sin $stry
	}

	# Retry fit
	set schk [fitConstrained $mtry $rtry rbounds]
	if {$schk < $bestStat - $epsStat} {
	    return -1
	}
	puts "Refit outer limit: $rtry $mtry $schk"
	if {$schk < $targetStat} {
	    # Not there yet
	    set min $mtry
	    set sin $schk
	} elseif {$stry - $schk < $epsStat} {
	    # Assume outer limit is OK
	    return [list $min $sin $pin $mtry $schk [parsave]]
	} else {
	    # Return to inner point and rerun fit
	    restore pin
	    set rin [$rprocName $fixed $min]
	    set sin [fitConstrained $min $rin rbounds]
	    if {$sin < $bestStat - $epsStat} {
		return -1
	    }
	    puts "revisit inner, radius, mass, stat: $rin $min $sin"
	    set step [expr sqrt ($step)]
	}
    }
}
#
############################################################
#
# Find one confidence limit for the mass.
#
# targetStat - value of fit statistic at confidence limit
# min - starting mass (best fit on entry)
# sin - fit stat corresponding to min
# step - factor by which to step the mass
# fixed - fixed radius or mean density
# rbName - name of list of shell boundaries
# epsStat - target precision
# rprocName - routine to compute radius
#
proc MassLimit {targetStat min sin step fixed rbName epsStat rprocName} {
    upvar $rbName rbounds
    set asymmetry 30.0
    set statBest $sin
    # Bracket the mass limit
    set bkret [bracket $targetStat $fixed $epsStat $min $sin $step \
		   rbounds $rprocName $statBest]
    if {$bkret == -1} {
	return -1
    }
    set min [lindex $bkret 0]
    set sin [lindex $bkret 1]
    set pin [lindex $bkret 2]
    set mout [lindex $bkret 3]
    set sout [lindex $bkret 4]
    set pout [lindex $bkret 5]
    puts "Bounds: $min $sin $mout $sout"
    # Restore parameters for inner limit, since it does better searching out
    restore pin
    puts "Refining"
    # Use binary search, since analytical methods interact badly with
    # fickle fit results
    while {abs ($sout - $sin) > $epsStat} {
	# Refine limits on the mass at the confidence limit
	set mtry [expr 0.5 * ($min + $mout)]
	set rtry [$rprocName $fixed $mtry]
	set stry [fitConstrained $mtry $rtry rbounds]
	if {$stry < $statBest - $epsStat} {
	    return -1
	} elseif {$stry > $targetStat} {
	    set mout $mtry
	    set sout $stry
	    set pout [parsave]
	} elseif {$asymmetry * ($stry - $sin) < $sout - $stry} {
	    # Fit stat at the outer limit is suspect (also gets stry < sin)
	    puts "Got $rtry $mtry $stry, searching out again"
	    set bkret [bracket $targetStat $fixed $epsStat $mtry $stry \
			   [expr $mout / $min] rbounds $rprocName $statBest]
	    if {$bkret == -1} {
		return -1
	    }
	    set min [lindex $bkret 0]
	    set sin [lindex $bkret 1]
	    set pin [lindex $bkret 2]
	    set mout [lindex $bkret 3]
	    set sout [lindex $bkret 4]
	    set pout [lindex $bkret 5]
	    restore pin
	} else {
	    set min $mtry
	    set sin $stry
	}
	puts "Bounds: $min $sin $mout $sout"
    }
    return [list $mout $pout]
}
#
############################################################
#
# Confidence range for the mass at fixed radius or density.
# Must be at best fit on entry.
#
# quantFixed - fixed quantity (radius or density)
# fixed - either the fixed radius or the fixed mean density
# rbName - shell boundaries
# delStat - increase in fit statistic at confidence limits
# grow - initial mass step when bracketing the limit
# epsStat - target accuracy
#
proc MassConf {quantFixed fixed rbName {delStat 1.0} {grow 1.03}
		 {epsStat 0.01}} {
    # Save best fit info
    set bestFitPars [parsave]
    # Compute constrained mass and radius at best fit
    upvar $rbName rbounds
    array set rmsProc {density rmSolFixedDen radius rmSolFixedRad}
    set t [$rmsProc($quantFixed) $fixed rbounds]
    set rbest [lindex $t 0]
    set bestMass [lindex $t 1]
    # Best fit statistic
    tclout stat
    set bestStat $xspec_tclout
    puts "Constrained best fit radius, mass and statistic:"
    puts "  $rbest $bestMass $bestStat"

    # Lower confidence limit
    set targetStat [expr $bestStat + $delStat]
    array set rProc {density fixedDen radius fixedRad}
    set radProc $rProc($quantFixed)
    set lowlim [MassLimit $targetStat $bestMass $bestStat [expr 1.0 / $grow] \
		    $fixed rbounds $epsStat $radProc]
    if {$lowlim == -1} {
 	puts stderr "Better fit found.  Rerun fit now."
 	return
    }

    # Output density parameters at the lower limit.
    # Break parameter values out into an array.
    foreach prs [lindex $lowlim 1] {
 	set idx [lindex $prs 0]
 	set pvals($idx) "[lindex $prs 1]"
    }
    puts "Density parameters:"
    set ct 0
    set diffpars [denPars rbounds]
    foreach dp $diffpars {
 	set pnum [lindex $dp 1]
 	set val $pvals($pnum)
 	if {$ct} {
 	    puts -nonewline ", [lindex $dp 0] $val"
 	} else {
 	    puts -nonewline "[lindex $dp 0] $val"
 	}
 	incr ct
    }
    set mlow [lindex $lowlim 0]
    set rlow [$radProc $fixed $mlow]
    puts "\nLower bound radius, mass: $rlow $mlow"

    # Search up
    puts "Searching up"
    restore bestFitPars
    set uplim [MassLimit $targetStat $bestMass $bestStat $grow $fixed \
		   rbounds $epsStat $radProc]
    if {$uplim == -1} {
	puts stderr "Better fit found.  Rerun fit now."
	return
    }

    # Output density parameters at upper limit
    foreach prs [lindex $uplim 1] {
 	set idx [lindex $prs 0]
 	set pvals($idx) "[lindex $prs 1]"
    }
    puts "Density parameters:"
    set ct 0
    foreach dp $diffpars {
	set pnum [lindex $dp 1]
	set val $pvals($pnum)
	if {$ct} {
	    puts -nonewline ", [lindex $dp 0] $val"
	} else {
	    puts -nonewline "[lindex $dp 0] $val"
	}
	incr ct
    }
    set mhi [lindex $uplim 0]
    set rhi [$radProc $fixed $mhi]
    puts "\nUpper bound radius, mass: $rhi $mhi"
    
    restore bestFitPars
    tclout chatter
    set oldchat $xspec_tclout
    chatter 5 5
    fit
    chatter [lindex $oldchat 0] [lindex $oldchat 1]

    return [list $bestMass $mlow $mhi]
}
