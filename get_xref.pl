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

my $USAGE = "get_xref.pl [-debug] [-h help]";

my (	$trap_id, 
     	$traphost, $trapuser, $trappass, $trapdbname, 
		$enshost, $ensuser, $enspass, $ensdbname,
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
$enshost = $conf{'enshost'};
$ensuser = $conf{'ensuser'};
$ensdbname = $conf{'ensdbname'};
$enspass = $conf{'enspass'};

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    	-host => $enshost,
        -user => $ensuser,
		-pass => $enspass,
		-database => $ensdbname
);

#### Connecting to unitrap_db
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect: ", $DBI::errst;
my $ensdb = DBI->connect("DBI:mysql:database=$ensdbname;host=$enshost;port=3306", $ensuser, $enspass) || die "Can't connect: ", $DBI::errst;

$trapdb->{RaiseError}=1;

my $st1 = "select distinct ensembl_id from trapens where (source = 'ENSMUSG' || source = 'REFSEQ')";
if ($debug) {print "SQL CODE: $st1\n";}
my $sth = $trapdb->prepare($st1);
$sth->execute();

while (my @rows = $sth->fetchrow_array) {
	my $ensid = $rows[0];
	$debug && print "---> $ensid\n";

	my $gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene');
	my $gene = $gene_adaptor->fetch_by_stable_id ($ensid);

	my @db_names = qw(RefSeq_dna RefSeq_peptide UniGene MGI miRBase MIM_MORBID GO Uniprot/SWISSPROT Uniprot/SPTREMBL);
	
	foreach my $db (@db_names) {
		print "DB $db\n";

		my @dblinks = @{ $gene->get_all_DBLinks($db) };
		
		foreach my $accession (@dblinks){
			print "Accession $accession\n";
			
			$accession =~ s/$db://;
			
			print "new Accession $accession\n";

			my %toinsert;
			$toinsert{'accession'} = $accession;
			$toinsert{'ensembl_id'} = $ensid;
			$toinsert{'dbname'} = $db;
			
			my $st = &prepare_stmt($trapdb, \%toinsert);
			my $returned = &insert_set($trapdb, $st, 'xref');
		}
	}
}

############################ sub

sub prepare_stmt {
    my ($dbh,  $par) = @_;
    my %params = %{$par};
    #construct statment nd quote
    my $stmt;
    
    foreach my $f (keys %params) {
        $stmt .= " , " if $stmt;
        $stmt .= $f . " = \"" . $params{$f} ."\" ";
    }
    return $stmt;
}

sub insert_set {
    my ($dbh, $stmt, $table_name) = @_;
    #print "here" . $stmt . "\n";
    my $s = "INSERT INTO $table_name SET  $stmt";
    print STDERR "### doing insert  $s ###\n";
    my $sth = $dbh->prepare($s);
    $sth->execute() || warn "insert failed : $DBI::errstr";
    my $dbi = $sth->{'mysql_insertid'};
    #print STDERR "the table $table_name last inserted id is $dbi \n";
    return $dbi;
}
