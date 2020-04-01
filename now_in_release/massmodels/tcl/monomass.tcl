# Procedures specific to the monomass model.
#
############################################################
#
# Allows generic check that the model is known
#
proc monomassKnown {} {
    return 1
}
#
############################################################
#
# Deduce maximum possible shells from total number of parameters
#
proc monomassMaxShells {mixpars} {
    return [expr ($mixpars - 4) / 2]
}
#
############################################################
#
# Model has density parameters
#
proc monomassHasDens {} {
    return 1
}
#
############################################################
#
# Info for applying a fixed mass constraint at rcut.  Makes
# a list in radial order, with each element listing a
# parameter number and the associated volume.
# rcut - radius at which mass is fixed
# rbName refers to an ordered list of rbounds like that 
#        returned by shellBounds
#
proc monomassConstraintData {rcut rbName} {
    upvar $rbName rbounds
    set diffpars [denPars rbounds]
    set ri [rinOfRadial 0 rbounds]
    set vnorm [expr 4.0 * $phyconst::pi / 3.0]
    foreach shell $rbounds dp $diffpars {
	if {[string length $dp] > 0} {
	    # Outer radius for current shell
	    set ro [lindex $shell 1]
	    if {$ro > $rcut} {
		set ro $rcut
	    }
	    # Volume between rcut and rinner to which the density
	    # difference for the current shell contributes
	    set vin [expr $vnorm * ($ro - $ri) \
			 * ($ro * ($ro + $ri) + $ri * $ri)]
	    # Record parameter number and effective volume
	    lappend mcdat [list [lindex $dp 1] $vin]
	}
    }
    return $mcdat
}
#
############################################################
#
# List shell densities in radial order
#
proc monomassShellDensities {rbName} {
    upvar $rbName rbounds
    set diffpars [denPars rbounds]
    set dprev 0.0
    # Sum density differences to make a list of shell densities
    for {set i [expr [llength $diffpars] - 1]} {$i >= 0} {incr i -1} {
	set pnum [lindex [lindex $diffpars $i] 1]
	tclout param $pnum
	set d [expr [lindex $xspec_tclout 0] + $dprev]
	if {[info exists dens]} {
	    set dens [linsert $dens 0 $d]
	} else {
	    set dens $d
	}
	set dprev $d
    }
    return $dens
}
#
############################################################
#
# For monomass, set initial model parameters to values estimated
# for a perfect isothermal sphere.
# rbName = name of list giving radial boundaries
# kTval = gas temperature in keV
# rUnit = model distance unit (m)
# rhoCritModel = critical density in model units
#
proc monomassISstartPars {rbName kTval rUnit rhoCritModel} {
    upvar $rbName rb

    # Set temperatures
    set nshell [llength $rb]
    setMMkT $kTval $nshell

    # Active density parameters
    set denpars [denPars rb]
    set nactive [llength $denpars]

    # Mass density estimates
    set dens [denIsothermal rb $kTval $rUnit $rhoCritModel]

    set ishell 0
    foreach dp $denpars {
	set delden [lindex $dens $ishell]
	incr ishell
	if {$ishell < $nactive} {
	    set delden [expr $delden - [lindex $dens $ishell]]
	}
	set lowdel [expr 0.001 * $delden]
	set hidel [expr 1000 * $delden]
	# newpar arguments: <par number> <par value> <delta> <hard low>
	# <soft low> <soft high> <hard high>
	newpar [lindex $dp 1] $delden $lowdel $lowdel $lowdel $hidel $hidel
    }
}
#
############################################################
#
# For monomass, set initial model parameters to values estimated
# for an NFW potential.
# rbName = name of list giving radial boundaries
# kTval = gas temperature in keV
# conc = concentration parameter
# rUnit = model distance unit (m)
# rhoCritModel = critical density in model units
#
proc monomassNFWstartPars {rbName kTval conc rUnit rhoCritModel} {
    upvar $rbName rb

    # Set temperatures
    set nshell [llength $rb]
    setMMkT $kTval $nshell

    # Active density parameters
    set denpars [denPars rb]
    set nactive [llength $denpars]

    # Mass density estimates
    set nfwa [lindex [NFWguess $kTval $conc $rUnit $rhoCritModel] 0]
    set dens [denNFW rb $nfwa $conc $rhoCritModel]

    set ishell 0
    foreach dp $denpars {
	set delden [lindex $dens $ishell]
	incr ishell
	if {$ishell < $nactive} {
	    set delden [expr $delden - [lindex $dens $ishell]]
	}
	set lowdel [expr 0.001 * $delden]
	set hidel [expr 1000 * $delden]
	# newpar arguments: <par number> <par value> <delta> <hard low>
	# <soft low> <soft high> <hard high>
	newpar [lindex $dp 1] $delden $lowdel $lowdel $lowdel $hidel $hidel
    }
}
