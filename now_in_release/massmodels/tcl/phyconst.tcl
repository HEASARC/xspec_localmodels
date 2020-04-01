# Provide physical constants
#
namespace eval phyconst {

    # Values defined in the file of constants (in order)
    set bcdict {
	{{Speed of light} {c_light}} 
	{{Planck constant} {hbar}}
	{{Newton constant} {G_newton}}
	{{Boltzmann constant} {k_boltzmann}}
	{{Fine structure constant} {alpha_fs}}
	{{Electron charge} {q_electron}}
	{{Electron mass} {m_electron}}
	{{Proton mass} {m_proton}}
	{{Alpha particle mass} {m_alpha}}
	{{Tropical year \(days\)} {year_days}}
	{{Astronomical unit} {Astr_unit}}
	{{GM for Sun} {GMsun}}
    }

    # Slurp constants from file
    set constFile Basic_Constants
    set fp [open $constFile r]
    set slurp [read $fp]
    close $fp

    # Break into lines and use to set values for export.
    # NB: file order must match bcdict.
    set lines [split $slurp "\n"]
    foreach cl $lines cf $bcdict {
	# Ignore unmatched lines
	if {[string length $cf] > 0} {
	    # Construct regexp to match entry
	    set label [lindex $cf 0]
	    set name [lindex $cf 1]
	    set rgs "$label = (\\S+)"
	    # Extract value
	    regexp "$rgs" $cl junk x
	    # Export under name from bcdict
	    variable $name $x
	}
    }

    # Provide pi
    variable pi 3.1415926535897932384

    # Derived constants
    # Helium mass fraction
    variable f_He 0.25
    # Hydrogen mass fraction
    variable f_H [expr 1.0 - $f_He]
    # Permeability of free space
    variable mu_0 [expr 4e-7 * $pi]
    # Permittivity of free space
    variable epsilon_0 [expr 1.0 / ($mu_0 * $c_light * $c_light)]
    variable year_secs [expr $year_days * 24.0 * 3600.0]
    variable parsec [expr $Astr_unit * 180.0 * 3600.0 / $pi]
    variable kiloparsec [expr $parsec * 1000.0]
    variable megaparsec [expr $parsec * 1e6]
    variable Msun [expr $GMsun / $G_newton]
    variable m_H [expr $m_proton + $m_electron]
    variable m_He [expr $m_alpha + 2.0 * $m_electron]
    variable keV [expr 1000.0 * $q_electron]
    set tfac [expr 2.0 * $f_H + 3.0 * $f_He * $m_H / $m_He]
    # Mean mass per particle
    variable mu [expr 1.0 / $tfac]
    set xfac [expr $f_H + 2.0 * $f_He * $m_H / $m_He]
    # n_tot / n_e
    variable ntot_ne [expr $tfac / $xfac]
    # n_H / n_e
    variable nh_ne [expr $f_H / $xfac]
    # rho / n_e
    variable rho_ne [expr $m_H / $xfac]
    # Radiation density constant
    proc radiation_density {} {
	variable k_boltzmann
	variable c_light
	variable hbar
	set t [expr $k_boltzmann / ($c_light * $hbar)]
	set pi $phyconst::pi
       	return [expr ($pi * $pi / 15.0) * $k_boltzmann * $t * $t * $t]
    }
    variable SB_density [radiation_density]
    # Stefan-Boltzmann constant
    variable SB_sigma [expr $SB_density * $c_light / 4.0]
    # Classical electron radius
    proc electron_radius {} {
	variable q_electron
	variable c_light
	variable m_electron
	variable epsilon_0
	set pi $phyconst::pi
	return [expr $q_electron * $q_electron / (4.0 * $pi \
		        * $epsilon_0 * $m_electron * $c_light * $c_light)]
    }
    variable r_electron [electron_radius]
    variable sigma_thomson [expr (8.0 * $pi / 3.0) \
				* $r_electron * $r_electron]
}
