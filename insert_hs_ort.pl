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
use Bio::EnsEMBL::Registry;

my $USAGE = "insert_go_.pl [-debug] [-h help]";

my (	 $ens_hs_host, $ens_hs_user, $ens_hs_pass, $ens_hs_dbname, 
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
$ens_hs_host = $conf{'ens_hs_dbhost'};
$ens_hs_user = $conf{'ens_hs_dbuser'};
$ens_hs_dbname = $conf{'ens_hs_dbname'};
$ens_hs_pass = $conf{'ens_hs_pass'};

my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect to $trapdbname";
#my $ens_hs_db = DBI->connect("DBI:mysql:database=$ens_hs_dbname;host=$ens_hs_host;port=3306", $ens_hs_user, $ens_hs_pass) || die "Can't connect;";

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    	-host => $ens_hs_host,
        -user => $ens_hs_user,
		-pass => $ens_hs_pass,
		-database => $ens_hs_dbname
);

my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "gene" );

my $file = shift @ARGV;
open (FH, $file);
while (<FH>){
	my $row = $_;
	chomp $row;
	my ($ens_id,$ens_hs_id) = split /\t/;
	
	$ens_hs_id =~ s/\n//;

	if ($debug) {print $ens_id."\t".$ens_hs_id."\n";}	

	if ($ens_hs_id) {
		my $gene_hs = $gene_adaptor->fetch_by_stable_id($ens_hs_id);
		my $ort_name = $gene_hs->external_name;
		
		my $st =  "insert into ortholog set mm_ensembl_id = \"$ens_id\", hs_ensembl_id = \"$ens_hs_id\", ort_name = \"$ort_name\"";
		if ($debug) {print "SQL CODE: $st\n";}
		my $sth = $trapdb->prepare($st);
		$sth->execute();
	}
}
