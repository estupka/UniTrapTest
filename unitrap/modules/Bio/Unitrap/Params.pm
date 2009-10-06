package Bio::Unitrap::Params;

use strict;

sub params_check () {
	my ($param, $value, $USAGE) = @_;

	if (!$value) {
	    die($USAGE . "\n\tMust specify a valid parameter for $param\n");
	}
}

1;
