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

my (	$traphost, $trapuser, $trappass, $trapdbname, 
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

my $file = shift @ARGV;
open (FH, $file);
while (<FH>){
	my $row = $_;
	chomp $row;
	my ($ens_hs_id, $omim_id, $desc) = split /\t/;
	
	$desc =~ s/\n//;

	if ($debug) {print $ens_hs_id."\t".$omim_id."\t".$desc."\n";}	

	if ($omim_id) {
		my $st =  "select ortholog_id from ortholog where hs_ensembl_id = \"$ens_hs_id\"";
		if ($debug) {print "SQL CODE: $st\n";}
		my $sth = $trapdb->prepare($st);
		$sth->execute();
		
		while (my $rhref = $sth->fetchrow_hashref) {
			my $orth_id = $rhref->{'ortholog_id'};
			
			my $st_ins =  "insert into omim set ortholog_id = \"$orth_id\", accession = \"$omim_id\", description = \"$desc\"";
			if ($debug) {print "SQL CODE: $st_ins\n";}
			my $sth_ins = $trapdb->prepare($st_ins);
			$sth_ins->execute();
		}	
	}
}
