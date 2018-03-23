#!/usr/bin/perl
#
# Runs initpackage and adds GSL to the build by editing extra flags
# into the Makefile.
#
# If you include GSL in the main build for xspec, this script is
# pointless.

# Porting checks:
# Your path to perl (1st line above).
# Compile flags for including the GSL headers, in $GSLinclude.
# Linker flags for GSL, in $GSLlink.

use strict;

if (-f "Makefile") {
    # Overkill?
    my $cleancmd = "hmake clean";
    system $cleancmd;
}

my $initcmd = "initpackage clmass model.dat " . $ENV {"PWD"};
system $initcmd || die "Did you heainit?  Failed $initcmd\n";

# Compile flags giving path to GSL header files.
# Something like "-I/usr/local/include" if needed.
my $GSLinclude = "";

# Link flags to incorporate GSL libraries for your system.
# Something like "-L/usr/local/lib -lgsl -lgslcblas" if the
# default below fails.
# If the string is long, you may want to split it over several
# lines to match the format of the Makefile, ie
# "lots of stuff \\\n\t\t\t  and more stuff"
my $GSLlink = "-lgsl -lgslcblas";

my $mf = "Makefile";
my $mfo = $mf . ".$$";
my $maxlen = 72;
open MAK, "<$mf" || die "Failed to open $mf\n";
open EMK, ">$mfo" || die "Failed to open $mfo\n";
while (<MAK>) {
    if (/^HD_CXXFLAGS/) {
	&addtodef ($GSLinclude, $maxlen);
    } elsif (/^HD_SHLIB_LIBS/) {
	&addtodef ($GSLlink, $maxlen);
    } else {
	print EMK;
    }
}

unlink $mf;
rename $mfo, $mf;


############################################################
#
sub addtodef {
    my ($mess, $lmax) = @_;

    # Emit continuation lines un altered
    while (/^.*\\$/) {
	    print EMK;
	    # Just in case we hit EOF
	    last if (!defined ($_ = <MAK>));
    }

    # Edit and emit final line
    chomp;
    # Expand tabs to determine visible lengths
    my $t = $_;
    $t =~ s/\t/        /g;
    my $lenin = length ($t);
    $t = $mess;
    $t =~ s/\t/        /g;
    my $lenadd = length ($t);
    if (length ($_) == 0) {
	print EMK "\t\t\t  $mess\n\n";
    } elsif ($lenin + $lenadd < $lmax) {
	print EMK "$_ $mess\n";
    } else {
	print EMK "$_ \\\n\t\t\t  $mess\n";
    }
}
