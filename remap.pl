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
use Bio::Unitrap::Db;
use Bio::Tools::Run::AnalysisFactory::Pise;
use Bio::SeqIO;
use Bio::EnsEMBL::Registry;

my $USAGE = "remap_cluster.pl id lastid";

# configuration options
my %conf =  %::conf;
my $tmp_dir = $conf{'tmp_dir'};
my $hit_db = $conf{'hit_db'};
my $debug = $conf{'debug'};
my ($idcount, $id);

#connect to the trapdb
my $traphost = $conf{'traphost'};
my $trapuser = $conf{'trapuser'};
my $trapdbname = $conf{'trapdbname'};
my $trappass = $conf{'trappass'};

#### Connecting to unitrap_db
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect!";

#### Connecting to Ensembl core database
my $enshost = $conf{'enshost'};
my $ensuser = $conf{'ensuser'};
my $ensdbname = $conf{'ensdbname'};
my $enspass = $conf{'enspass'};

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    	-host => $enshost,
        -user => $ensuser,
		-pass => $enspass,
		-database => $ensdbname
);

########################
#### Get BATCHARG params
#########################
die print "Must provide id and idcount as BATCHARG" if scalar(@ARGV) == 0;
chomp $ARGV[0];
chomp $ARGV[1];

$id = $ARGV[0];
$idcount = $ARGV[1];

unless ($id && $idcount) {
    print STDERR " no trap id \n";
    exit;
}

print STDOUT "ID $id\tID count $idcount\n";

my $tid;
my $lastid = $id + $idcount - 1;

foreach ($tid = $id; $tid <= $lastid; $tid++) {
	my $st1 = "select distinct ensembl_id, gene_chr, gene_start, gene_end from mutated_gene where mutated_gene_id = $tid;";
	$debug && print "SQL CODE: $st1\n";
	my $sth = $trapdb->prepare($st1);
	$sth->execute();
	
	while (my @rows = $sth->fetchrow_array){
		my $ensid = $rows[0];
		my $chr = $rows[1];
		my $start = $rows[2];
		my $end = $rows[3];
		
		$debug && print "Ensembl_id: $ensid\n";
		
		my $region;
		if ($chr =~ /NT/) {
			$region = "supercontig";
		} else {
			$region = "chromosome";
		}
		
		my $sliceadaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $sliceadaptor->fetch_by_region ($region,$chr,$start,$end);
		
		my $fastafile = $tmp_dir."fasta/".$ensid.".genomic.fa";
		
		####create seqio object
		my $out = Bio::SeqIO->new(      -file => ">$fastafile" , '-format' => 'Fasta');
		
		my $seqobj = Bio::PrimarySeq->new ( -seq => $slice->seq,
						   -id  => $ensid,
						   -desc => "genomic region (Chr $chr: $start..$end +)",
						   -alphabet => 'dna'
						   );
		
		### Write GENOMIC sequence in the fasta dir
		$out->write_seq ($seqobj);
		
		### Launch remap in order to get the HTML page to display in the UNITRAP website
		$debug && print STDOUT "remap -sequence $fastafile -enzymes 'EcoRI,BamHI,XbaI,HindIII' -HTML Y -sitelen 4 -outfile ".$tmp_dir."remap/".$ensid.".html";
		system ("remap -sequence $fastafile -enzymes 'EcoRI,BamHI,XbaI,HindIII' -HTML Y -sitelen 4 -outfile ".$tmp_dir."remap/".$ensid.".html");
		
		### Working to retrieve fragments
		my @enzymes = qw(EcoRI BamHI XbaI HindIII);
		#my @enzymes = qw(EcoRI);
		
		foreach my $enzyme (@enzymes) {
			my $this_start = 1;
			
			### Launch restrict in order to get fragment information
			my $outfile = $tmp_dir."restrict/".$ensid.".$enzyme.restrict";
			my $outpng = $tmp_dir."restrict/png/".$ensid.".$enzyme.restrict";
			$debug && print STDOUT "restrict -sequence $fastafile -enzymes '$enzyme' -sitelen 4 -outfile $outfile";
			system ("restrict -sequence $fastafile -enzymes '$enzyme' -sitelen 4 -outfile $outfile");
			
			#### Parse the restrict output
			open (FILE, $outfile);
			open (OUT, ">$outfile.parsed");
			
			my $last_check;
			while (defined (my $line=<FILE>)) {
				chomp($line);
				#regular expression to match start and the end of the sequence
				if ($line =~ /^#\sSequence.*from:\s+(\d+).*to:\s+(\d+)/) {
					my ($seqstart, $seqend) = ($1,$2);
					
					#print the DNA start and the beginning in the input file
					print OUT "Start\t$seqstart\nEnd\t$seqend\n\ngroup\n";
				}
				
				#regular expression to match start, end and name of the enzyme
				if ($line =~ /^\s+([0-9]+)\s+([0-9]+)\s+(\w+)\s+/) {
					my ($enz_start,$enz_end,$this_enzyme) = ($1,$2,$3);
					print OUT "label\nTick\t$enz_start\t8\n$this_enzyme\nendlabel\n";
					
					print STDOUT "+>>>>> THIS START $this_start\n";
					print STDOUT "+>>>>> $this_enzyme: $enz_start\t$enz_end\n";
					
					######################
					#to convert between coordinates
					print STDOUT "SUBSLICE: $this_start..".($enz_start)."\n";
					my $subslice =  $slice->sub_Slice($this_start,$enz_start-1);
					my $sliceproj = $subslice->project($region, $hit_db);
					
					foreach my $seg (@$sliceproj) {
						my $sl = $seg->to_Slice();
						print $subslice->seq_region_name(), ':', $seg->from_start(), '-',$seg->from_end(), ' -> ', $sl->seq_region_name(), ':', $sl->start(), '-',$sl->end(), $sl->strand(), "\n",$sl->seq(),"\n\n";
						
						### insert fragment into db
						my $fragment_seq;
						if ($this_start == 1) {
						    $fragment_seq = '...'.$sl->seq();
						} else {
						    $fragment_seq = $sl->seq();
						}
						
						my %mutated_gene_fragment;
						$mutated_gene_fragment{'ensembl_id'} = $ensid;
						$mutated_gene_fragment{'chr'} = $chr;
						$mutated_gene_fragment{'start'} = $sl->start();
						$mutated_gene_fragment{'end'} = $sl->end();
						$mutated_gene_fragment{'sequence'} = $fragment_seq;
						$mutated_gene_fragment{'enzyme'} = $enzyme;
						$debug && print Dumper %mutated_gene_fragment;
						
						my $st_mutated_gene_fragment = &prepare_stmt($trapdb, \%mutated_gene_fragment);
						my $mutated_gene_fragment_id = &insert_set($trapdb, $st_mutated_gene_fragment, 'mutated_gene_fragment');
						
						if ($mutated_gene_fragment_id == 0) {
							print "ID=0\tExiting\n";
							exit;
						}
						
						### internal probe design
						my $fragfile = $tmp_dir."/fragments/probedesign.".$mutated_gene_fragment_id;
						my $probefile = $tmp_dir."/probedesign/probedesign.".$mutated_gene_fragment_id.".out";

						####create seqio object
						my $out2 = Bio::SeqIO->new(      -file => ">$fragfile", '-format' => 'Fasta');
						
						my $seqobj2 = Bio::PrimarySeq->new ( -seq => $sl->seq(),
										   -id  => $mutated_gene_fragment_id,
										   -alphabet => 'dna'
										   );
										   
						### Write fragment sequence
						$out2->write_seq ($seqobj2);
						
						$debug && print STDOUT "eprimer3 -sequence $fragfile -task 1 -productsizerange 300-500 -outfile $probefile";
						system ("eprimer3 -sequence $fragfile -task 1 -productsizerange 300-500 -outfile $probefile");
					}
					
					## change next start to the enzyme end
					$this_start = $enz_start;
					$last_check = 1;
				} elsif ($last_check == 1 && $line eq '') {
					### retrieve sequence from last restriction site to the end of the slice
					#to convert between coordinates
					print STDOUT "SUBSLICE: $this_start..".length($slice->seq())."\n";
					my $subslice =  $slice->sub_Slice($this_start,length($slice->seq()));
					my $sliceproj = $subslice->project($region, $hit_db);
					
					foreach my $seg (@$sliceproj){
						my $sl = $seg->to_Slice();
						print $subslice->seq_region_name(), ':', $seg->from_start(), '-',$seg->from_end(), ' -> ',
						$sl->seq_region_name(), ':', $sl->start(), '-',$sl->end(), $sl->strand(), "\n",$sl->seq(),"\n\n";
						
						my $fragment_seq = $sl->seq().'...';
						
						my %mutated_gene_fragment;
						$mutated_gene_fragment{'ensembl_id'} = $ensid;
						$mutated_gene_fragment{'chr'} = $chr;
						$mutated_gene_fragment{'start'} = $sl->start();
						$mutated_gene_fragment{'end'} = $sl->end();
						$mutated_gene_fragment{'sequence'} = $fragment_seq;
						$mutated_gene_fragment{'enzyme'} = $enzyme;
						$debug && print Dumper %mutated_gene_fragment;
						
						my $st_mutated_gene_fragment = &prepare_stmt($trapdb, \%mutated_gene_fragment);
						my $mutated_gene_fragment_id = &insert_set($trapdb, $st_mutated_gene_fragment, 'mutated_gene_fragment');
						
						if ($mutated_gene_fragment_id == 0) {
							print "ID=0\tExiting\n";
							exit;
						}
						
						### internal probe design
						my $fragfile = $tmp_dir."/fragments/probedesign.".$mutated_gene_fragment_id;
						
						####create seqio object
						my $out2 = Bio::SeqIO->new(      -file => ">$fragfile", '-format' => 'Fasta');
						
						my $seqobj2 = Bio::PrimarySeq->new ( -seq => $sl->seq(),
										   -id  => $mutated_gene_fragment_id,
										   -alphabet => 'dna'
										   );
										   
						### Write fragment sequence
						$out2->write_seq ($seqobj2);
						
						$debug && print STDOUT "eprimer3 -sequence $fragfile -task 1 -productsizerange 300-500 -outfile $fragfile.out";
						system ("eprimer3 -sequence $fragfile -task 1 -productsizerange 300-500 -outfile $fragfile.out");
					}
				}
			}
			
			print OUT "endgroup";
			close (OUT);
			
			### launch lindna
			print STDOUT "lindna $outfile.parsed -goutfile $outpng -graph png -auto";
			system ("lindna $outfile.parsed -goutfile $outpng -graph png -auto");
		}
	}
}

#######################
#
#     SUBS 
#
#######################

sub prepare_stmt {
    my ($dbh,  $par) = @_;
    my %params = %{$par};
    #construct statment nd quote
    my $stmt;
    
    foreach my $f (keys %params) {
        $stmt .= " , " if $stmt;
        $stmt .= $f . " = " . $dbh->quote ($params{$f});
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
