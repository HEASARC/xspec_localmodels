# Procedures specific to nfwmass model.
#
############################################################
#
# Allows generic check that the model is known
#
proc nfwmassKnown {} {
    return 1
}
#
############################################################
#
# Deduce maximum possible shells from total number of parameters
#
proc nfwmassMaxShells {mixpars} {
    return [expr $mixpars - 6]
}
#
############################################################
#
# Model has density parameters
#
proc nfwmassHasDens {} {
    return 0
}
#
############################################################
#
# Mass in model units
#
proc nfwmassMofR {r rbName} {
    tclout compinfo 1
    set mixpars [lindex $xspec_tclout 2]
    set pnuma [expr $mixpars - 1]
    set pnumpot $mixpars
    tclout pinfo $pnuma
    set pa [lindex $xspec_tclout 0]
    tclout pinfo $pnumpot
    set ppot [lindex $xspec_tclout 0]
    if {![string match nfwa $pa] || ![string match nfwpot $ppot]} {
	puts stderr "nfwmassMofR: wrong model parameters"
	return
    }
    tclout param $pnuma
    set nfwa [lindex $xspec_tclout 0]
    tclout param $pnumpot
    set nfwpot [lindex $xspec_tclout 0]
    return [expr $nfwpot * $nfwa * [NFWmassForm [expr $r / $nfwa]]]
}
#
############################################################
#
# For nfwmass, set initial model parameters estimated from a
# temperature and concentration parameter.
# rbName = name of list giving radial boundaries
# kTval = gas temperature in keV
# conc = concentration parameter
# rUnit = model distance unit (m)
# rhoCritModel = critical density in model units
#
proc nfwmassNFWstartPars {rbName kTval conc rUnit rhoCritModel} {
    upvar $rbName rb

    # Set temperatures
    set nshell [llength $rb]
    setMMkT $kTval $nshell

    # Estimate NFW parameters
    set nfwpar [NFWguess $kTval $conc $rUnit $rhoCritModel]

    # Set nfwa
    set pno [expr [fixedPars] + [maxShells] + 1]
    set nfwa [lindex $nfwpar 0]
    set alow [expr 0.01 * $nfwa]
    set ahi [expr 100 * $nfwa]
    # newpar arguments: <par number> <par value> <delta> <hard low>
    # <soft low> <soft high> <hard high>
    newpar $pno $nfwa $alow $alow $alow $ahi $ahi

    # Set nfwnorm
    incr pno
    tclout param $pno
    set normmax [lindex $xspec_tclout 5]
    set nfwnorm [lindex $nfwpar 1]
    if {$nfwnorm > $normmax} {
	puts stderr "Estimated nfwnorm exceeds hard limit - unusable"
	return
    }
    set normhi [lindex $xspec_tclout 4]
    set normlo [expr 0.01 * $nfwnorm]
    set normup [expr 100 * $nfwnorm]
    if {$normup < $normhi} {
	set normhi $normup
	set normmax $normup
    } elseif {$normup < $normmax} {
	set normmax $normup
    }
    newpar $pno $nfwnorm $normlo $normlo $normlo $normhi $normmax
}
