#! /usr/bin/perl -w

# Guglielmo Roma
# guglielmo.roma@gmail.com

$|=1;
use strict;
BEGIN {
    require "/home/roma/src/scripts/unitrap2/unitrap_conf.pl";
};

use strict;
use DBI;
use Getopt::Long;

my $USAGE = "insert_go_.pl [-debug] [-h help]";

my (
     	$traphost, $trapuser, $trappass, $trapdbname, 
		$debug, $help);

&GetOptions (	'debug'	=> \$debug,
				'help|h' => \$help);

#### Configuration options
my %conf =  %::conf;

#connect to the trapdb
$traphost = $conf{'traphost'};
$trapuser = $conf{'trapuser'};
$trapdbname = $conf{'trapdbname'};
$trappass = $conf{'trappass'};
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect to $trapdbname";


#my $go_class = "molecular function";
#my $go_class = "biological process";
my $go_class = "cell location";

my $file = shift @ARGV;
open (FH, $file);
while (<FH>){
	my $row = $_;
	chomp $row;
	my ($ens_id,$accession,$desc,$source) = split /,/;
	
	$source =~ s/\n//;

	if ($debug) {print $desc."\t".$source."\t".$accession ."\t".$go_class."\n";}	
	if ($accession) {
		my $st =  "insert into trapens_go set ensembl_id = \"$ens_id\", accession = \"$accession\", source = \"$source\", description = \"$desc\", go_class = \"$go_class\"";
		if ($debug) {print "SQL CODE: $st\n";}
		my $sth = $trapdb->prepare($st);
		$sth->execute();
	}
}
