#! /usr/bin/perl -w 
#
# Guglielmo Roma
# guglielmo.roma@gmail.com
#
# Run this script to assign a trapped gene to each trapping event, predict vector insertion, and create the UniTraps
# 

$|=1;
use strict;
BEGIN {
    require "/home/roma/src/scripts/unitrap2/unitrap_conf.pl";
};

use DBI;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::Seq;
use Bio::Unitrap::File;
use Bio::Tools::Primer3;
use Bio::EnsEMBL::Registry;

my ( $enshost, $ensuser, $enspass, $ensdbname, $ensestdbname,
     $traphost, $trapuser, $trappass, $trapdbname,
     $blastnexec, $blastdb, $hit_db, $idcount, $id, $tmp_dir, $mysql_path, $repeatMasker_dir, $debug, $debugSQL);

GetOptions( 
	    "enshost=s",    	\$enshost,
	    "ensuser=s",    	\$ensuser,
	    "enspass=s",    	\$enspass,
	    "ensdbname=s",  	\$ensdbname,
	    "ensestdbname=s",  	\$ensestdbname,
	    "traphost=s",   	\$traphost,
	    "trapuser=s",   	\$trapuser,
	    "trapdbname=s", 	\$trapdbname,
	    "trappass=s",   	\$trappass,
	    "blastdb=s",    	\$blastdb,
	    "hit_db=s",    		\$hit_db,
	    "id=s", 	    	\$id,
	    "idcount=s",    	\$idcount,
	    "debug",        	\$debug,
	    "debugSQL",        	\$debugSQL
	    );

#### Configuration options
my %conf =  %::conf;
$enshost = $conf{'enshost'};
$ensuser = $conf{'ensuser'};
$ensdbname = $conf{'ensdbname'};
$ensestdbname = $conf{'ensestdbname'};
$enspass = $conf{'enspass'};
$traphost = $conf{'traphost'};
$trapuser = $conf{'trapuser'};
$trapdbname = $conf{'trapdbname'};
$trappass = $conf{'trappass'};
$blastdb = $conf{'blastdb'};
$hit_db = $conf{'hit_db'};
$blastnexec ||=  $::conf{'blastnexec'};
$tmp_dir = $conf{'tmp_dir'};
$mysql_path = $conf{'mysql_path'};
$repeatMasker_dir = $conf{'repeatMasker_dir'};
$debug = $conf{'debug'};
$debugSQL = $conf{'debugSQL'};

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

#### Connecting to unitrap_db
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect: ", $DBI::errst;

#### Connecting to Ensembl core database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(	-host => $enshost,	-user => $ensuser,	-pass => $enspass,	-database => $ensdbname	);

my $ensdb = DBI->connect("DBI:mysql:database=$ensdbname;host=$enshost;port=3306", $ensuser, $enspass) || die "Can't connect: ", $DBI::errst;
									
$trapdb->{RaiseError}=1;
$trapdb->{RaiseError}=1;

my $error;
my $status;

my $tid;
my $lastid = $id + $idcount - 1;

foreach ($tid = $id; $tid <= $lastid; $tid++) {
	my $sth_map = $trapdb->prepare("select tm.trapmap_id, tm.start, tm.end, tm.hit_id, tm.chosen from trapmap tm where tm.trap_id = $tid");
	$debug && print STDOUT "SQL CODE: select tm.trapmap_id, tm.start, tm.end, tm.hit_id, tm.chosen from trapmap tm where tm.trap_id = $tid\n";
	my $nres = $sth_map->execute;
	
	### If this trap id is in the database continue
	if ($nres >= 1) {
		while (my $rhref = $sth_map->fetchrow_hashref) {
			my $trapmap_id = $rhref->{'trapmap_id'};
			my $start = $rhref->{'start'};
			my $end = $rhref->{'end'};
			my $hit_id = $rhref->{'hit_id'};
			my $chosen = $rhref->{'chosen'};
			
			my $strand;
			my @features;
			
			print STDOUT "### TRAP: $tid\tTRAPMAP $trapmap_id\n";
			
			my $sth_block = $trapdb->prepare("select tb.trapblock_id, tb.start, tb.end, tb.strand from trapblock tb where tb.trapmap_id=$trapmap_id");
			$debug && print STDOUT "SQL CODE: select tb.trapblock_id, tb.start, tb.end, tb.strand from trapblock tb where tb.trapmap_id=$trapmap_id\n";
			my $nres_block = $sth_block->execute;
			
			if ($nres_block >= 1) {
				while (my $rhref_block = $sth_block->fetchrow_hashref) {
					my $trapblock_id = $rhref_block->{'trapblock_id'};
					my $tbstart = $rhref_block->{'start'};
					my $tbend = $rhref_block->{'end'};
					$strand = $rhref_block->{'strand'};
					
					### Create a feature for each trap_block
					my $feature = Bio::SeqFeature::Generic->new(	-display_name => $trapblock_id,
																	-start => $tbstart,
																	-end => $tbend,
																	-strand => $strand	);
					push @features,$feature;
					
					$debug && print STDOUT "Trapblock $trapblock_id \n";
				}
			} else {
				print STDERR "This trapmap ($trapmap_id) has no blocks!\n";
			}
			
			if ($chosen == 1) {
				&assing_trapped_gene ($tid, $trapmap_id, $strand);
			}
		}
		
		$status = 1;

	} else {
		print STDERR "this trap ($tid) was not mapped\n";
		my $trap_id = $tid;

		$error = "not mapped trap";
	}
}

#######################
#
#     SUBS 
#
#######################

sub assing_trapped_gene {
    my ($trap_id, $tm_chosen, $strand) = @_;
    my $sth_map = $trapdb->prepare("select distinct m.trapmap_id, t.mol_type, t.race, t.project_id, t.vector_name from trapmap m, trap t where t.trap_id = m.trap_id and m.trap_id = $trap_id and m.trapmap_id = $tm_chosen and m.chosen = 1;");
    my $res = $sth_map->execute;
    
    while (my $maps = $sth_map->fetchrow_hashref) {
		my $trapmap_id = $maps->{'trapmap_id'};
		my $mol_type = $maps->{'mol_type'};
		my $race = $maps->{'race'};
		my $project_id = $maps->{'project_id'};
		my $vector_name = $maps->{'vector_name'};
		my $returned;
		
		### Looking for the right strand
		my $sql_strand;
		
		### Take into account that 5'RACE from tigem needed to be reversed
		if ($project_id == 1 && $race eq "5'") {
			if ($strand == 1) {
				$sql_strand = " c.strand = '-1' ";
			} elsif ($strand == -1) {
				$sql_strand = " c.strand = '1' ";
			} else {
				print STDERR "Strand $strand is not 1 or -1\n";
				print STDERR "TRAP $trap_id\tTRAPMAP $trapmap_id\nExiting\n";
				exit;
			}
		} else {
			$sql_strand = " c.strand = '$strand' ";
		}
		
		my ($refseq, $ensg, $fantom, $unigene, $ensestg, $enscdna, $ensest, $abinitio, $tclg);
		
		$debug && print STDOUT "\n######### TRAP $trap_id - TRAPMAP $trapmap_id\n";
		
		if ($res) {
			$debug && print STDOUT "1) Check for RefSeq ---";
			my $sql1 = "select c.label, c.ensmusg, c.ensmuse, c.ensmust, c.dS_s, c.dE_e , m.hit_id, t.start, t.end, m.trapmap_id, m.trap_id, c.refseq_id from trapblock_ensg_annotation c, trapblock t, trapmap m where t.trapmap_id = m.trapmap_id and c.trap_id = $trap_id and m.trapmap_id = $trapmap_id and c.trapblock_id = t.trapblock_id and m.chosen = '1' and c.refseq_id != '' and $sql_strand order by t.start";
			$debug && print STDOUT "SQL CODE: $sql1\n";
			$refseq = &check_ens_annotation ($sql1, 'REFSEQ', $mol_type, $ensdb, 'trapblock_ensg_annotation', $race, $project_id, $vector_name, 1);
			
			if ($refseq) {
				my $st = &prepare_stmt($trapdb, $refseq);
				$returned = &insert_set($trapdb, $st, 'trapens');
				next;
			}
			
			unless ($refseq) {
				$debug && print STDOUT "\n2) Check for ENSG ---";
				my $sql = "select c.label, c.ensmusg, c.ensmuse, c.ensmust, c.dS_s, c.dE_e , m.hit_id, t.start, t.end, m.trapmap_id, m.trap_id from trapblock_ensg_annotation c, trapblock t, trapmap m where t.trapmap_id = m.trapmap_id and c.trap_id = $trap_id and m.trapmap_id = $trapmap_id and c.trapblock_id = t.trapblock_id and m.chosen = '1' and $sql_strand order by t.start";
				$debug && print STDOUT "SQL CODE: $sql\n";
				$ensg = &check_ens_annotation ($sql, 'ENSMUSG', $mol_type, $ensdb, 'trapblock_ensg_annotation', $race, $project_id, $vector_name, 0);
			}
			
			if ($ensg) {
				my $st = &prepare_stmt($trapdb, $ensg);
				$returned = &insert_set($trapdb, $st, 'trapens');
				next;
			}
			
			unless ($ensg) {
				$debug && print STDOUT "\n3) Check for UNIGENE --- ";
				my $sql = "select c.trap_id, c.trapblock_id, c.display_name, m.trapmap_id from trapblock_unigene_annotation c, trapblock t, trapmap m where m.trapmap_id = t.trapmap_id and c.trap_id = $trap_id and m.trapmap_id = $trapmap_id and t.trapblock_id = c.trapblock_id and $sql_strand ";
				$debug && print STDOUT "SQL CODE: $sql\n";
				$unigene = &check_other_annotation($sql, 'UNIGENE');
			}
			
			if ($unigene) {
				my $st = &prepare_stmt($trapdb, $unigene);
				$returned = &insert_set($trapdb, $st, 'trapens');
				next;
			}
			
			unless ($unigene) {
				$debug && print STDOUT "\n6) Check for CDNA --- ";
				my $sql = "select c.trap_id, c.trapblock_id, c.display_name, m.trapmap_id from trapblock_cdna_annotation c, trapblock t, trapmap m where m.trapmap_id = t.trapmap_id and c.trap_id = $trap_id and m.trapmap_id = $trapmap_id and $sql_strand ;";
				$debug && print STDOUT "SQL CODE: $sql\n";
				$enscdna = &check_other_annotation($sql, 'CDNA');
			}
			
			if ($enscdna) {
				my $st = &prepare_stmt($trapdb, $enscdna);
				$returned = &insert_set($trapdb, $st, 'trapens');
				next;
			}
			
			unless ($enscdna) {
				$debug && print STDOUT "\n7) Check for EST --- ";
				my $sql  = "select c.trap_id, c.trapblock_id, c.display_name, m.trapmap_id from trapblock_ensest_annotation c, trapblock t, trapmap m where m.trapmap_id = t.trapmap_id and c.trap_id = $trap_id and m.trapmap_id = $trapmap_id and t.trapblock_id = c.trapblock_id and $sql_strand ";
				$debug && print STDOUT "SQL CODE: $sql\n";
				$ensest = &check_other_annotation($sql, 'EST');
			}
			
			if ($ensest) {
				my $st = &prepare_stmt($trapdb, $ensest);
				$returned = &insert_set($trapdb, $st, 'trapens');
				next;
			} else {
				$debug && print STDOUT "NO ANNOTATION FOUND!/n";
			}			
		}
		
		$debug && print STDOUT "\tdone!\n";
	}
}

sub check_ens_annotation {
	my ($sql, $source, $mol_type, $dbh, $table, $race, $project_id, $vector_name, $refseq_check) = @_;
	my $sth1 = $trapdb->prepare($sql);
	my $retr = $sth1->execute;
	my $amb = 0;
	my %toinsert;
	my $c;
	my $g;
	
	### FIND the ortholog, OMIM and xref, too.
	### first time gene is inserted, desing primer for pcr, too.
	
	while (my $href = $sth1->fetchrow_hashref) {
		$toinsert{'trap_id'} = $href->{'trap_id'};
		$toinsert{'trapmap_id'} = $href->{'trapmap_id'};
		$toinsert{'hit_id'} = $href->{'hit_id'};
		
		if ($href->{'ensmusg'}) {
			$toinsert{'ensembl_id'} = $href->{'ensmusg'};
			$g = 1;
			$c++;
		}
		
		if ($refseq_check == 1) {
			$toinsert{'refseq_id'} = $href->{'refseq_id'};
		} else {
			$toinsert{'refseq_id'} = '';
		}
		
		if ($mol_type eq "mRNA" && !$href->{'ensmuse'}) {
			$amb = 1;
		}
	}
	
	$toinsert{'annotation_ambiguous'} = $amb;
	$toinsert{'number_trapblocks'} = $retr;
	$toinsert{'number_annotated_trapblocks'} = $c;
	$toinsert{'source'} = $source;
	
	if ($mol_type) {
		$toinsert {'mol_type'} = $mol_type;
	}
	
	if ($g) {
		$debug && print STDOUT "-->>Finding the vector insertion site\n";
		
		if ($mol_type eq 'mRNA') {
			my $insert = &find_insertion_site (\%toinsert, $dbh, $table, $race, $project_id, $vector_name, $refseq_check);
			
			if ($insert) {
				%toinsert = %$insert;
				
				### Build Unitrap only if a RefSeq or an Ensembl Gene has been trapped
				### only if trap is mRNA, and there are both flanking and trapped exons
				if ($toinsert {'insertion_ambiguous'} == 0 && $toinsert {'flanking_exon_id'} && $toinsert {'exon_id'} && ($source eq "REFSEQ" || $source eq "ENSMUSG")) {
					### create a new UNITRAP
					my $unitrap_id = &build_unitrap (\%toinsert, $source, $debug);
					
					my %toinsert_trap_unitrap;
					$toinsert_trap_unitrap{'unitrap_id'} = $unitrap_id;
					$toinsert_trap_unitrap{'trap_id'} = $toinsert{'trap_id'};
					$toinsert_trap_unitrap{'trapmap_id'} = $toinsert{'trapmap_id'};
					
					$debug && print Dumper %toinsert_trap_unitrap;
					
					my $st_trap_unitrap = &prepare_stmt($trapdb, \%toinsert_trap_unitrap);
					my $trap_unitrap_id = &insert_set($trapdb, $st_trap_unitrap, 'trap_unitrap');
					
					&calculate_mutated_protein_foreach_transcript ($unitrap_id, $toinsert{'ensembl_id'}, $toinsert{'putative_insertion_start'}, $toinsert{'putative_insertion_end'}, $toinsert{'hit_id'}, $debug);
					&calculate_primers_for_vector_validation ($unitrap_id, $toinsert{'ensembl_id'}, $toinsert{'putative_insertion_start'}, $toinsert{'putative_insertion_end'}, $toinsert{'hit_id'}, $debug);
				}
			}
		} else {
			### Unitraps from GENOMIC DNA traps and SPLK!
			my $insert = &find_insertion_site_gdna (\%toinsert, $dbh, $table, $race, $project_id, $vector_name, $refseq_check);
			
			if ($insert) {
				%toinsert = %$insert;
				
				### Build Unitrap only if a RefSeq or an Ensembl Gene has been trapped
				### only if trap is mRNA, and there are both flanking and trapped exons
				if ($toinsert {'insertion_ambiguous'} == 0 && $toinsert {'flanking_exon_id'} && $toinsert {'exon_id'} && ($source eq "REFSEQ" || $source eq "ENSMUSG")) {
					### create a new UNITRAP
					my $unitrap_id = &build_unitrap (\%toinsert, $source, $debug);
					
					my %toinsert_trap_unitrap;
					$toinsert_trap_unitrap{'unitrap_id'} = $unitrap_id;
					$toinsert_trap_unitrap{'trap_id'} = $toinsert{'trap_id'};
					$toinsert_trap_unitrap{'trapmap_id'} = $toinsert{'trapmap_id'};
					
					$debug && print Dumper %toinsert_trap_unitrap;
					
					my $st_trap_unitrap = &prepare_stmt($trapdb, \%toinsert_trap_unitrap);
					my $trap_unitrap_id = &insert_set($trapdb, $st_trap_unitrap, 'trap_unitrap');
					
					&calculate_mutated_protein_foreach_transcript ($unitrap_id, $toinsert{'ensembl_id'}, $toinsert{'putative_insertion_start'}, $toinsert{'putative_insertion_end'}, $toinsert{'hit_id'}, $debug);
					&calculate_primers_for_vector_validation ($unitrap_id, $toinsert{'ensembl_id'}, $toinsert{'putative_insertion_start'}, $toinsert{'putative_insertion_end'}, $toinsert{'hit_id'}, $debug);
				}
			}
		}
		
		if (($source eq "ENSMUSESTG" || $source eq "GENSCAN") && !$toinsert{'exon_id'}) {
			return undef;
		} else {
			return (\%toinsert);
		}
	} else {
		return undef;
	}
}

sub build_unitrap () {
	my ($toinsert, $source, $debug) = @_;
	my %toinsert = %{$toinsert};
	
	### check if unitrap exists
	my $sql = "select unitrap_id, accession, chr, start, end from unitrap where ensembl_id = '".$toinsert{'ensembl_id'}."' and ((exon_id = '".$toinsert{'exon_id'}."' and flanking_exon_id = '".$toinsert{'flanking_exon_id'}."') or (flanking_exon_id = '".$toinsert{'exon_id'}."' and exon_id = '".$toinsert{'flanking_exon_id'}."'));";
	$debug && print STDOUT "$sql\n";
	my $sth = $trapdb->prepare($sql);
	my $res = $sth->execute;
	
	if ($res == 1) {
		my $rhref = $sth->fetchrow_hashref;
		my $unitrap_id = $rhref->{'unitrap_id'};
		my $accession = $rhref->{'accession'};
		my $chr = $rhref->{'chr'};
		my $start = $rhref->{'start'};
		my $end = $rhref->{'end'};
		
		$debug && print STDOUT "###### UNITRAP_ID already exists: $unitrap_id\t";
		
		if ($chr ne $toinsert{'hit_id'} || $start != $toinsert{'putative_insertion_start'} || $end != $toinsert{'putative_insertion_end'}) {
			print STDOUT "CHR $chr vs. ".$toinsert{'hit_id'}."\t";
			print STDOUT "START $start vs. ".$toinsert{'putative_insertion_start'}."\t";
			print STDOUT "END $end vs. ".$toinsert{'putative_insertion_end'}."\t";
			print STDOUT "Should EXIT..."; #exit;
		}
		
		return $unitrap_id;
	} else {
		my %toinsert_unitrap;
		
		my $uni_accession = &get_next_accession ($trapdb, $debug);
		
		$toinsert_unitrap {'accession'} = $uni_accession;
		$toinsert_unitrap {'chr'} = $toinsert{'hit_id'};
		$toinsert_unitrap {'start'} = $toinsert{'putative_insertion_start'};
		$toinsert_unitrap {'end'} = $toinsert{'putative_insertion_end'};
		$toinsert_unitrap {'hit_db'} = $hit_db;
		
		#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($toinsert{'ensembl_id'});

		my ($gene_adaptor, $gene);	
		if ($toinsert->{'ensembl_id'} =~ 'ENSMUSESTG') {
			$gene_adaptor  = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Gene');
			$gene = $gene_adaptor->fetch_by_stable_id ($toinsert->{'ensembl_id'});	
		} elsif ($toinsert->{'ensembl_id'} =~ 'ENSMUSG') {
			$gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene');
			$gene = $gene_adaptor->fetch_by_stable_id ($toinsert->{'ensembl_id'});		
		} else {
			print "No ENSG, no ENSESTG. exiting...\n";
			exit; 
		}

		$toinsert_unitrap {'ensembl_id'} = $toinsert{'ensembl_id'};
		
		if ($toinsert {'name'} ne '') {
			$toinsert_unitrap {'gene_name'} = $toinsert {'name'};
		}
		
		$toinsert_unitrap {'gene_strand'} = $toinsert {'ens_strand'};
		$toinsert_unitrap {'gene_description'} = $toinsert {'description'};		
		
		$toinsert_unitrap {'gene_start'} = $gene->start;
		$toinsert_unitrap {'gene_end'} = $gene->end;
		
		$toinsert_unitrap {'gene_type'} = $source;
		
		if ($toinsert{'refseq_id'} ne '') {
			$toinsert_unitrap {'refseq_id'} = $toinsert{'refseq_id'};
		}
		
		$toinsert_unitrap {'exon_id'} = $toinsert{'exon_id'};
		$toinsert_unitrap {'flanking_exon_id'} = $toinsert{'flanking_exon_id'};
		
		my $st_unitrap = &prepare_stmt($trapdb, \%toinsert_unitrap);
		my $unitrap_id = &insert_set($trapdb, $st_unitrap, 'unitrap');
		return $unitrap_id;
	}
}

sub calculate_primers_for_vector_validation () {
	my ($unitrap_id, $ensembl_id, $vec_start, $vec_end, $chr, $debug) = @_;
	my $path = "/tmp/";
	my $region;
	
	if ($chr =~ /NT/) {
		$region="supercontig";
	} else {
		$region = "chromosome";
	}
	
	#my $sliceadaptor = $ensdb->get_SliceAdaptor;

	my ($sliceadaptor, $geneadaptor, $gene);
	if ($ensembl_id =~ 'ENSMUSESTG') {
		$sliceadaptor = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Slice' );
		$geneadaptor  = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Gene' );
		$gene = $geneadaptor->fetch_by_stable_id ($ensembl_id);
	} elsif ($ensembl_id =~ 'ENSMUSG') {
		$sliceadaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		$geneadaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
		$gene = $geneadaptor->fetch_by_stable_id ($ensembl_id);
	} else {
		print "No ENSG, no ENSESTG. exiting...\n";
		exit; 
	}
	
	################################## pcr on the DNA with these alternative REVERSE primers
	################################################## manipulating primers sequences
	# reverse primers for the rt-pcr on the DNA
	############################################### reverse primer 1
	my $L232 = 'GATGTGCTGCAAGGCGATTA';
	
	my $L232_obj = Bio::Seq->new( 	-display_id => 'L232',
									-seq => $L232);
	
	my $L232_revcomp_obj = $L232_obj->revcom();
	
	# reverse complement of primer 1
	my $L232_revcomp = $L232_revcomp_obj->seq();
	############################################### reverse primer 2
	my $LacZ = 'CCAGGGTTTTCCCAGTCACG';
	
	my $LacZ_obj = Bio::Seq->new( 	-display_id => 'LacZ',
									-seq => $LacZ);
	
	my $LacZ_revcomp_obj = $LacZ_obj->revcom();
	
	# reverse complement of primer 2
	my $LacZ_revcomp = $LacZ_revcomp_obj->seq();
	
	$debug && print STDERR "############################## DNA PCR\nREVERSE PRIMER 1(L232): $L232\nREVERSE PRIMER 2(LacZ): $LacZ\n";
	
	#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($ensembl_id);
		
	my $gene_strand = $gene->strand;
	
	$debug && print STDOUT "### Gene strand $gene_strand\n";
	
	# aggiungo alla sequenza intronica 50 basi dell'esone che si trova di fronte al vettore
	my ($slice_start, $slice_end);
	
	if ($gene_strand == 1) {
		$slice_start = $vec_start - 50;
		$slice_end = $vec_end;
	} elsif ($gene_strand == -1) {
		$slice_start = $vec_start; 
		$slice_end = $vec_end + 50;
	}
	
	$debug && print STDOUT "###SLICE COORDS: CHR $chr, START $slice_start, END $slice_end\n\n";
	
	if ($slice_start >= $slice_end) {
		print "Problem... START > END\n"; 
		exit;
	}
	
	# prendo dalla slice la sequenza "repeat masked"
	my $slice = $sliceadaptor->fetch_by_region($region, $chr, $slice_start, $slice_end);
	my $rm_slice = $slice->get_repeatmasked_seq();
	my $rm_seq = lc($rm_slice->seq);
	my $total_seq_length = length ($rm_seq);
	
	# se il gene si trova sullo strand negativo, calcolo la sequenza reverse complement per girarla in direzione 5'->3'.
	unless ($gene_strand == 1) {
		my $rm_seq_obj = Bio::Seq->new( 	-display_id => 'seq',
											-seq => $rm_seq);

		my $rm_revcomp_seq_obj = $rm_seq_obj->revcom();
		my $rm_revcomp_seq = $rm_revcomp_seq_obj->seq();
		$rm_seq = $rm_revcomp_seq;
	}
	
	$debug && print STDERR "Sequence for PCR\tLength: $total_seq_length\n";
		
	my $j;
	my ($pcr1,$pcr2);
			
	# Partendo dalla posizione 0 (che coincide SEMPRE con l'inizio delle 50 basi aggiunte di esone)
	# e per intervalli di 3000 basi, calcolo delle sottosequenze di una lunghezza pari a 100 basi
	# su cui trovare i primers.
	#													 |- 3000 -|
	#													__  basi  _
	#													| |		 | |						
	#	 ________________                ________________ |      | |          _________________               ______________
	#	|  		|--------------|trapped exon   |||------|-|----------|                |--------------|             \  
	#	|________________|	 	|________________||      | |          |________________|              |_____________/
	#													| |		 | |	  XXXX	 
	#													| |		 | |    <-LacZ reverse primer
	#													|_|		 |_|
	#												    100      100
	#												   basi     basi
	#					forward primers da cercare	    ->		  ->
	#
	# successivamente aggiungo a queste sottosequenze una sequenza reverse complement dei primer (ex. LacZ, L232) che dovranno funzionare come reverse.  
	
	for ($j=0; $j<$total_seq_length; $j+=3000) {
		my $sub_seq = substr ($rm_seq, $j, 100);
		my $seq_for_pcr_1 = ">$unitrap_id-$j\n".$sub_seq.$L232_revcomp."\n";
		my $seq_for_pcr_2 = ">$unitrap_id-$j\n".$sub_seq.$LacZ_revcomp."\n";
		
		$pcr1 .= $seq_for_pcr_1;
		$pcr2 .= $seq_for_pcr_2;
		
		$debug && print STDERR "100-bases sequence starting from position $j\n$seq_for_pcr_1\n$seq_for_pcr_2\n";
	}
	
	my $file1 = "$path"."$unitrap_id"."_"."L232";
	my $outfile1 = "$path"."$unitrap_id"."_"."L232.out";
	my $file2 = "$path"."$unitrap_id"."_"."LacZ";
	my $outfile2 = "$path"."$unitrap_id"."_"."LacZ.out";
	
	Bio::Unitrap::File->writeDataToFile ($file1, $pcr1);
	Bio::Unitrap::File->writeDataToFile ($file2, $pcr2);
	
	$debug && print STDERR "eprimer3 -sequence '$file1' -outfile '$outfile1' -task 2 -reverseinput '$L232'\n";
	$debug && print STDERR "eprimer3 -sequence '$file2' -outfile '$outfile2' -task 2 -reverseinput '$LacZ'\n";
		
	system ("eprimer3 -sequence '$file1' -outfile '$outfile1' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 2 -reverseinput '$L232'");
	system ("eprimer3 -sequence '$file2' -outfile '$outfile2' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 2 -reverseinput '$LacZ'");
	
	&parse_and_insert ($outfile1, $slice_start, $slice_end, $gene_strand, $L232, "L232", $chr, $debug);
	&parse_and_insert ($outfile2, $slice_start, $slice_end, $gene_strand, $LacZ, "LacZ", $chr, $debug);
	
	system ("rm $file1");
	system ("rm $outfile1");
	system ("rm $file2");
	system ("rm $outfile2");
}

sub calculate_mutated_protein_foreach_transcript () {
	my ($unitrap_id, $ensembl_id, $vec_start, $vec_end, $chr, $debug) = @_;

	#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($ensembl_id);

	my ($geneadaptor, $gene);
	if ($ensembl_id =~ 'ENSMUSESTG') {
		$geneadaptor  = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Gene' );
		$gene = $geneadaptor->fetch_by_stable_id ($ensembl_id);
	} elsif ($ensembl_id =~ 'ENSMUSG') {
		$geneadaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
		$gene = $geneadaptor->fetch_by_stable_id ($ensembl_id);
	} else {
		print "No ENSG, no ENSESTG. exiting...\n";
		exit; 
	}

	my $gene_strand = $gene->strand;
	my @transcripts = @{$gene->get_all_Transcripts()};
	my $vec_ins_point = $vec_start + (($vec_end - $vec_start)/2);
	
	$debug && print STDOUT "VEC POINT: $vec_ins_point\n";
	
	foreach my $transcript (@transcripts) {
		my $transcript_id = $transcript->stable_id();
		my @exons = @{$transcript->get_all_Exons()};
		my ($protein_seq, $mutagenicity);
		my $not_mutated_exons_num = 0;
		my $mutated_exons_num = 0;
		my $not_mutated_domains_num = 0;
		my $mutated_domains_num = 0;
		
		if ($gene_strand == 1) {
			@exons = sort {$a->end <=> $b->end} @exons;
		} else {
			@exons = sort {$b->end <=> $a->end} @exons;
		}
		
		foreach my $exon (@exons) {
			my $exon_start = $exon->start();
			my $exon_end = $exon->end();
			
			# Retrieves the portion of the transcripts peptide encoded by this exon			
		 	my $peptide = $exon->peptide($transcript)->seq;
			
			if (($gene_strand == 1 && $vec_ins_point >= $exon_end) || ($gene_strand == -1 && $vec_ins_point <= $exon_start)) {
				$protein_seq .= uc ($peptide);  # Non mutated amminoacides - uppercase
				$not_mutated_exons_num++;
			} else {
				$protein_seq .= lc ($peptide);	# Mutated amminoacides - lowercase
				$mutated_exons_num++;
			}
			
			$debug && print STDOUT "PEPTIDE $peptide\n";
		}
		
		my @mutated;
		my @not_mutated;
		
		if ($protein_seq) {
			my @protDB = qw(PFam);
			
			foreach my $db (@protDB) {
				my @domains_f = @{$transcript->translation->get_all_ProteinFeatures($db)};
				
				foreach my $domain (@domains_f) {
					my $term = $domain->display_id;				
					my $domain_start = $domain->start;
					my $domain_end = $domain->end;
					my $region_acc_no = $domain->interpro_ac;
					my $idesc = $domain->idesc;
					
					$debug && print STDOUT "DOMAIN $region_acc_no >> TAG $idesc\n";
					
					my @gencoords = $transcript->pep2genomic($domain_start, $domain_end);
					
					foreach my $obj (@gencoords) {
						if ($obj->isa("Bio::EnsEMBL::Mapper::Gap")){next;}
						my $genomic_start = $obj->start;
						my $genomic_end = $obj->end;
						my $genomic_strand = $obj->strand;
						
						if (($gene_strand == 1 && $genomic_end <= $vec_ins_point) || ($gene_strand == -1 && $genomic_start >= $vec_ins_point)) {
							if ($idesc." ($term)") {
								push @not_mutated, $idesc." ($term)";							
							}

							$debug && print STDOUT "THIS DOMAIN IS UPSTREAM THE INSERTION SITE - NOT MUTATED\n";
						} else {
							if ($idesc." ($term)") {
								push @mutated, $idesc." ($term)";							
							}
							
							$debug && print STDOUT "THIS DOMAIN IS DOWNSTREAM THE INSERTION SITE - MUTATED\n";
						}
					}
				}
			}
			$debug && print STDOUT Dumper @mutated;
			$debug && print STDOUT Dumper @not_mutated;
			
 			my $not_mutated = &unique_array (\@not_mutated, $debug);
			my $mutated = &unique_array (\@mutated, $debug);
			
			if ($not_mutated) {
				@not_mutated = @{$not_mutated};
			} else {
				@not_mutated = '';
			}
			
			if ($mutated) {
				@mutated = @{$mutated};
			} else { 
				@mutated = '';
			}
			
                        my $not_mutated_domains = join (",", @not_mutated);
                        my $mutated_domains = join (",", @mutated);
			
			my ($not_mutated_domains_num,$mutated_domains_num);	
			if ($not_mutated_domains ne '') {
				$not_mutated_domains_num = scalar(@not_mutated); 
			} else {
				$not_mutated_domains_num = 0;
			}
			
			if ($mutated_domains ne '') {
				$mutated_domains_num = scalar(@mutated);
			} else {
				$mutated_domains_num = 0;
			}
			
			$debug && print "MUTATED DOMAINS: ";
			$debug && print Dumper @mutated;
			$debug && print "\nNOT MUTATED DOMAINS: ";
			$debug && print Dumper @not_mutated;
			$debug && print "\n";
			$debug && print "not_mutated = $not_mutated_domains_num\t";
			$debug && print "mutated = $mutated_domains_num";
			$debug && print "\n";
			
			my $rank = 0;
			
			if ($not_mutated_domains_num == 0 && $mutated_domains_num > 0) {
				$mutagenicity = "high"; # The vector insertion occurred in the first introns and they knock-out most/all protein domains of the protein";
				$rank = 1;
			} elsif ($not_mutated_domains_num >= 1 && $mutated_domains_num >= 1) {		
				$mutagenicity = "medium"; # At least one domain is still intact so cannot be clear that the function is completely ablished";
				$rank = 2;
			} elsif ($not_mutated_domains_num >= 1 && $mutated_domains_num == 0) {		
				$mutagenicity = "low"; # The vector insertion is after all functional protein domains";
				$rank = 3;
			}
			
			my %toinsert;
			$toinsert{'unitrap_id'} = $unitrap_id;
			$toinsert{'protein_seq'} = $protein_seq;
			$toinsert{'mutagenicity'} = $mutagenicity;
			$toinsert{'rank'} = $rank;
			$toinsert{'transcript'} = $transcript_id;
			$toinsert{'ensembl_id'} = $ensembl_id;
			$toinsert{'not_mutated_exons_num'} = $not_mutated_exons_num;
			$toinsert{'mutated_exons_num'} = $mutated_exons_num;
			
			$toinsert{'not_mutated_domains'} = $not_mutated_domains;
			$toinsert{'mutated_domains'} = $mutated_domains;
			
			$toinsert{'not_mutated_domains_num'} = $not_mutated_domains_num;
			$toinsert{'mutated_domains_num'} = $mutated_domains_num;
			
			my $st_protein_mutagenicity = &prepare_stmt($trapdb, \%toinsert);
			my $protein_mutagenicity_id = &insert_set($trapdb, $st_protein_mutagenicity, 'protein_mutagenicity');		
		}
	}
}

sub parse_and_insert () {
	my ($file, $slice_start, $slice_end, $gene_strand, $reverse_seq, $reverse_name, $chr, $debug) = @_;
	
	# parse primer3 output to get some data
	my $p3 = Bio::Tools::Primer3->new(-file=>$file);
	
	# get all the results
	my $all_results = $p3->all_results;
	
	my %primers;
	my $unitrap_id;
	my $sub_slice_start;
	
	# NB. bioperl parser works only on a canonic primer3 output (in this case I'm looking only for forward primers and somenthing goes wrong!!!)
	# make some variations
	
	foreach my $key (keys %{$all_results}) {				
		my $line = $key;
		
		if ($line =~ /RESULTS/ && $line =~ /#/) {
			my @words = split (" ", $line);
			($unitrap_id, $sub_slice_start) = split ("-", $words[4]);
			$debug && print STDERR "UNITRAP ID: $unitrap_id\tINTERVAL $sub_slice_start\tREVERSE IS $reverse_name\n";
		}
		
		if ($unitrap_id && $line !~ /#/ && $line =~ /PRIMER/) {
			$line =~ s/\s+/&/g;
			my @splitted = split (/([&])/, $line);
			
			my $start = $splitted[6];
			my $length = $splitted[8];
			my $tm = $splitted[10];
			my $gc = $splitted[12];
			my $seq = $splitted[14];

			$debug && print STDERR "PARSED LINE: start $start, length $length, tm $tm, gc $gc, seq $seq\n\n";

			my $forward_genomic_start;
			
			if ($gene_strand == 1) { 
				$forward_genomic_start = $slice_start + $start + $sub_slice_start; 
			} elsif ($gene_strand == -1) { 
				$forward_genomic_start = $slice_end - $start - $sub_slice_start - $length;
			}
			
			my $forward_genomic_end = $forward_genomic_start + $length;

			$debug && print STDERR ">>> $forward_genomic_start - $forward_genomic_end ||| Slice start $slice_start, Start $start, SUB SLICE START $sub_slice_start, length $length\n\n";
			
			if ($gene_strand == 1) {
				$gene_strand = '+';
			} elsif ($gene_strand == -1) {
				$gene_strand = '-';
			}
			
			if ($tm >= 58) {
				my $sql = "insert into pcr set unitrap_id = '$unitrap_id', hit_id = '$chr', forward = '$seq', forward_start = '$forward_genomic_start', forward_end = '$forward_genomic_end', forward_strand = '$gene_strand', reverse = '$reverse_seq', reverse_name = '$reverse_name', gc_perc = '$gc', tm = '$tm'";
				$debug && print STDERR "SQL CODE: $sql\n";
				my $sth = $trapdb->prepare($sql);
				$sth->execute();
				last;
			}
		}
	}
}

sub find_insertion_site {
	my ($toinsert, $dbh, $table, $race, $project_id, $vector_name, $refseq_check) = @_;
	my %insert = %{$toinsert};
	
	my ($gene_adaptor, $gene);	
	if ($toinsert->{'ensembl_id'} =~ 'ENSMUSG') {
		$gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene');
		$gene = $gene_adaptor->fetch_by_stable_id ($toinsert->{'ensembl_id'});		
	} else {
		print "No ENSG. exiting...\n";
		exit; 
	}
	
	my $gene_strand = $gene->strand;
	my $sql1 = "select t.trapblock_id, e.ensmusg, e.ensmuse, e.ensmust, m.hit_id, t.start, t.end, t.strand, e.exon_coverage from $table e, trapblock t, trapmap m where t.trapmap_id = m.trapmap_id and e.trap_id = m.trap_id and e.trap_id = ".$toinsert->{'trap_id'}." and m.trapmap_id = ".$toinsert->{'trapmap_id'}." and e.trapblock_id = t.trapblock_id and m.chosen = 1 ";
	
	if ($refseq_check == 1) {
		$sql1 .= " and e.refseq_id != '' ";
	}
	
	$sql1 .= "order by t.start;";
	print STDOUT "SQL CODE: $sql1\n";
	my $sth1 = $trapdb->prepare($sql1);
	my $retr1 = $sth1->execute;
	
	if ($retr1 > 0) {
		my @blocks;
		
		while (my $h = $sth1->fetchrow_hashref()) {
			push @blocks, $h;
		}
		
		my $case;							# case 1: trapped_exon and upstream_exon (i.e., polyA trap - 3' race)
											# case 2: trapped_exon and downstream_exon (i.e., 5' race with reversed trap)
		my $blocks = scalar (@blocks);		# case 3: trapped_exon and downstream_exon (i.e., 5' race with not-reversed trap)
		my $last_block_num = $blocks - 1;
		my $trap_strand = $blocks[0]->{'strand'};
		
		if ($trap_strand == 1) {
			$insert {'trap_strand'} = '+';
		} elsif ($trap_strand == -1) {
			$insert {'trap_strand'} = '-';
		}
		
		my @exons = @{$gene->get_all_Exons()};
		my $exons = scalar (@exons);
		
		if ($gene_strand == 1) {
			@exons = sort {$a->end <=> $b->end} @exons;
		} else {
			@exons = sort {$b->end <=> $a->end} @exons;
		}
		
		$insert {'total_number_exons'} = $exons;
		my ($exon_id, $flanking_exon_id, $tag_point);
		$debug && print STDOUT "RACE = $race\tTRAP STRAND: $trap_strand\n GENE STRAND: $gene_strand\tGene contains $exons exons.\n";
		
		### Default race is 3'
		if ($race eq 'na') {
			$race = "3'";
		}
		
		### Checking for the right trapped exon... according to race type
		### If a trapped exon is not identified, a tag point is defined from the trapblock
		if ($race eq "5'" && $gene_strand eq "1") {
			$exon_id = $blocks[$last_block_num]->{'ensmuse'};
			
			unless ($exon_id) {
				$insert {'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
				$tag_point = $blocks[$last_block_num]->{'end'};
			}
			
			if ($gene_strand eq $trap_strand) {
				$case = 2;
			} elsif ($gene_strand ne $trap_strand) {
				### This is more or less the case in which the sequence-tag is not reversed
				### Note this should not happen for TIGEM at this point, since I have chosen 
				### all the trapped genes located on the opposite strand where race=5' and project is TIGEM
				$case = 3;
			}
		} elsif ($race eq "5'" && $gene_strand eq "-1") {
			$exon_id = $blocks[0]->{'ensmuse'};
			
			unless ($exon_id) {
				$insert {'putative_insertion_end'} = $blocks[0]->{'start'};
				$tag_point = $blocks[0]->{'start'};
			}
			
			if ($gene_strand eq $trap_strand) {
				$case = 2;
			} elsif ($gene_strand ne $trap_strand) {
				### This is more or less the case in which the sequence-tag is not reversed
				### Note this should not happen for TIGEM at this point, since I have chosen 
				### all the trapped genes located on the opposite strand where race=5' and project is TIGEM
				$case = 3;
			}
		}  elsif ($race eq "3'" && $gene_strand eq "-1") {
			$exon_id = $blocks[$last_block_num]->{'ensmuse'};
			
			unless ($exon_id) {
				$insert {'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
				$tag_point = $blocks[$last_block_num]->{'end'};
			}
			$case = 1;
		}  elsif ($race eq "3'" && $gene_strand eq "1") {
			$exon_id = $blocks[0]->{'ensmuse'};
			
			unless ($exon_id) {
				$insert {'putative_insertion_end'} = $blocks[0]->{'start'};
				$tag_point = $blocks[0]->{'start'};
			}
			$case = 1;
		}
		
		$debug && print STDOUT "CASE = $case...\n";
		
		### Once the trapped exon or the tag_point is identified, we identify the flaking exon according to the case
		if ($case) {
			$insert {'insertion_case'} = $case;
			if ($exon_id) {
				$insert {'exon_id'} = $exon_id;
				my $i = 1;
				my (%exon_hash, $trapped_exon_num, $trapped_exon_obj);
				
				foreach my $exon (@exons) {
					my $this_exon_id = $exon->stable_id();
					my $exon_start = $exon->start();
					my $exon_end = $exon->end();
					
					if ($exon_start > $exon_end) {$debug && print STDERR "EXON START > END"; exit;}
					
					if ($this_exon_id eq $exon_id) {
						print "OK.FOUND EXON TRAPPED: $exon_id\n";
						$trapped_exon_num = $i;
						$trapped_exon_obj = $exon;
						$insert {'trapped_exon_rank'} = $trapped_exon_num;
					}
					
					$exon_hash {$i} {'exon'} = $this_exon_id;
					$exon_hash {$i} {'obj'} = $exon;
					
					$debug && print STDERR "I: $i - TRAPPED EXON NUM $trapped_exon_num\tCASE:$case\tTHIS EXON:$this_exon_id\tTRAPPED EXON: $exon_id\n";
					$i++;
				}
				
				$debug && print STDERR "Gene has $exons EXONS; TRAPPED EXON NUM $trapped_exon_num\n";
				
				if ($trapped_exon_num == 1 && $case == 1) {
					$debug && print STDOUT "TRAPPED EXON = FIRST EXON of the gene, but CASE is 1... really STRANGE!!! Could be a missannotated gene\n";
					
					if ($exons == 1) {
						$debug && print STDOUT "1 exon-gene\n";
					}
					
					if ($gene_strand == '1') {
						$insert {'putative_insertion_end'} = $blocks[0]->{'start'};
					} elsif ($gene_strand == '-1') {
						$insert {'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
					}

					$insert {'insertion_ambiguous'} = 1;
				} elsif ($trapped_exon_num == $exons && ($case == 3 || $case == 2)) {
					$debug && print STDOUT "TRAPPED EXON = LAST EXON of the gene, but CASE is 2 or 3... really STRANGE!!! Could be a missannotated gene\n";
					
					if ($gene_strand == '1') {
						 $insert {'putative_insertion_start'} = $blocks[$last_block_num]->{'end'};
					} elsif ($gene_strand == '-1') {
						 $insert {'putative_insertion_end'} = $blocks[0]->{'start'};
					}
					
					$insert {'insertion_ambiguous'} = 1;
				} else {
					$debug && print STDERR "TRAPPED EXON = MIDDLE EXON of the gene... Checking for the intron\n";
					my ($flanking_exon_id, $flanking_exon_obj);
					
					# NB:
					# case 1: intron between trapped_exon and upstream_exon
					# case 2: intron between trapped_exon and downstream_exon
					# case 3: intron between trapped_exon and downstream_exon
					
					### Calculate the flanking_exon
					if ($case == 1) {
						 my $flanking_exon_num = $trapped_exon_num - 1;
						 
						 $debug && print STDOUT "Flanking_exon_num $flanking_exon_num\n";
						 $flanking_exon_id = $exon_hash {$flanking_exon_num} {'exon'};
						 $flanking_exon_obj = $exon_hash {$flanking_exon_num} {'obj'};
						 
						 $debug && print STDOUT "Vector is between $exon_id (the trapped exon) and $flanking_exon_id (upstream exon)\n";	 
					} elsif ($case == 2 || $case == 3) {
						 my $flanking_exon_num = $trapped_exon_num + 1;
						 
						 print "Flanking_exon_num $flanking_exon_num\n";
						 $flanking_exon_id = $exon_hash {$flanking_exon_num} {'exon'};
						 $flanking_exon_obj = $exon_hash {$flanking_exon_num} {'obj'};
						 
						 $debug && print STDOUT "Vector is between $exon_id (the trapped exon) and $flanking_exon_id (downstream exon)\n";		
					}
					
					### Handling the intron object
					if ($exon_id && $flanking_exon_id) {
						$insert {'flanking_exon_id'} = $flanking_exon_id;
						my $intron;

						if ($case == 1) {
							 $intron = new Bio::EnsEMBL::Intron ($flanking_exon_obj,$trapped_exon_obj);
						} elsif ($case == 2 || $case == 3) {
							 $intron = new Bio::EnsEMBL::Intron ($trapped_exon_obj,$flanking_exon_obj);
						}
						my $intron_start = $intron->start();
						my $intron_end = $intron->end();
						
						if ($intron_start < $intron_end) {
							 $insert {'putative_insertion_start'} = $intron_start;
							 $insert {'putative_insertion_end'} = $intron_end;
						} else {
							 $insert {'putative_insertion_start'} = $intron_end;
							 $insert {'putative_insertion_end'} = $intron_start;
						}
						
						$insert {'insertion_ambiguous'} = 0;
						
						$debug && print STDOUT "The insertion is in the INTRON between $exon_id & $flanking_exon_id - $intron_start..$intron_end\n";
					}
				}
			} elsif ($tag_point) {
				$debug && print STDOUT "No exon_id. But, tag_point\n";
				my ($pre_exon_id, $pre_exon_obj, $next_exon_id, $next_exon_obj);
				
				foreach my $exon (@exons) {
					my $this_exon_id = $exon->stable_id();
					my $exon_start = $exon->start();
					my $exon_end = $exon->end();
					
					if ($exon_start > $exon_end) {print STDERR "EXON START > END"; exit;}
					$debug && print STDOUT "TAG $tag_point - CURRENT EXON: $this_exon_id\n";
					
					if ($gene_strand == '1') {
						if ($tag_point >= $exon_end) {
							$pre_exon_id = $this_exon_id;
							$pre_exon_obj = $exon;
							next;
						} else {
							$next_exon_id = $this_exon_id;
							$next_exon_obj = $exon;
							last;
						}
					} elsif ($gene_strand == '-1') {
						if ($tag_point <= $exon_start) {
							$pre_exon_id = $this_exon_id;
							$pre_exon_obj = $exon;
							next;
						} else {
							$next_exon_id = $this_exon_id;
							$next_exon_obj = $exon;
							last;
						}
					}
				}
				
				$debug && print STDOUT "TAG-POINT: $tag_point\tPREV-EXON:$pre_exon_id\tNEXT-EXON:$next_exon_id\n";
				$insert {'insertion_ambiguous'} = 1;
				
				if ($pre_exon_id && $next_exon_id) {
					if ($gene_strand == '1') {
						if ($case == '1') {
							$insert {'putative_insertion_start'} = $pre_exon_obj->end;
						} elsif ($case == '2' || $case == '3') {
							$insert {'putative_insertion_end'} = $next_exon_obj->start;
						}
					} elsif ($gene_strand == '-1') {
						if ($case == '1') {
							$insert {'putative_insertion_end'} = $pre_exon_obj->start;
						} elsif ($case == '2' || $case == '3') {
							$insert {'putative_insertion_start'} = $next_exon_obj->end;
						}
					}
				}
				
				$insert {'insertion_ambiguous'} = 1;
				$insert {'label'} = 'novel_exon';
			} else {
				print STDERR "No exon_id, no tag_point.\n";
				exit;
			}
		} else {
			print STDERR "No case!!!\n";
			exit;
		}
	} else {
		print STDERR "No block!?!\n";
		exit;
	}
	
	# if ($insert {'putative_insertion_end'} && !$insert {'putative_insertion_start'}) {
		# $insert {'putative_insertion_start'} = $insert {'putative_insertion_end'} - 50;
	# } elsif (!$insert {'putative_insertion_end'} && $insert {'putative_insertion_start'}) {
		# $insert {'putative_insertion_end'} = $insert {'putative_insertion_start'} + 50;
	# }
	
	if ($gene_strand == 1) {
		$insert {'ens_strand'} = '+';
	} elsif ($gene_strand == -1) {
		$insert {'ens_strand'} = '-';
	}
	
	if ($gene->external_name) {
		$insert {'name'} = $gene->external_name;
	}
	
	if ($gene->description) {
		$insert {'description'} = $gene->description;
	}
	
	return (\%insert);
}

### find insertion site for genomic DNA traps (SPLK and genomic dna) - which corresponds to the coordinates of the trapped intron!
### Althoguh there is no race for genomic DNA traps, we use race 3' and 5' to distinguish the side of the sequencing the genomic DNA.
### Race is set for SPLK traps, but not for the genomic taps. Default is 5' race!
sub find_insertion_site_gdna {
	my ($toinsert, $dbh, $table, $race, $project_id, $vector_name, $refseq_check) = @_;
	my %insert = %{$toinsert};

	#my $gene = $dbh->get_GeneAdaptor->fetch_by_stable_id($toinsert->{'ensembl_id'});

	my ($gene_adaptor, $gene);	
	if ($toinsert->{'ensembl_id'} =~ 'ENSMUSESTG') {
		$gene_adaptor  = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Gene');
		$gene = $gene_adaptor->fetch_by_stable_id ($toinsert->{'ensembl_id'});	
	} elsif ($toinsert->{'ensembl_id'} =~ 'ENSMUSG') {
		$gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene');
		$gene = $gene_adaptor->fetch_by_stable_id ($toinsert->{'ensembl_id'});		
	} else {
		print "No ENSG, no ENSESTG. exiting...\n";
		exit; 
	}

	my $gene_strand = $gene->strand;

	my $sql1 = "select t.trapblock_id, e.ensmusg, e.ensmuse, e.ensmust, m.hit_id, t.start, t.end, t.strand, e.exon_coverage from $table e, trapblock t, trapmap m where t.trapmap_id = m.trapmap_id and e.trap_id = m.trap_id and e.trap_id = ".$toinsert->{'trap_id'}." and m.trapmap_id = ".$toinsert->{'trapmap_id'}." and e.trapblock_id = t.trapblock_id and m.chosen = 1 ";
	$sql1 .= "order by t.start;";
	print STDOUT "SQL CODE: $sql1\n";
	my $sth1 = $trapdb->prepare($sql1);
	my $retr1 = $sth1->execute;
	
	if ($retr1 > 0) {
		my @blocks;
		
		while (my $h = $sth1->fetchrow_hashref()) {
			push @blocks, $h;
		}
		
		my $blocks = scalar (@blocks);
		my $last_block_num = $blocks - 1;
		my $trap_strand = $blocks[0]->{'strand'};
		
		if ($trap_strand == 1) {
			$insert {'trap_strand'} = '+';
		} elsif ($trap_strand == -1) {
			$insert {'trap_strand'} = '-';
		}
		
		my @exons = @{$gene->get_all_Exons()};
		my $exons = scalar (@exons);
		
		if ($gene_strand == 1) {
			@exons = sort {$a->end <=> $b->end} @exons;
		} else {
			@exons = sort {$b->end <=> $a->end} @exons;
		}
		
		$insert {'total_number_exons'} = $exons;
		
		### 1) Identification of the tag_point based on the trap strand!
		$debug && print STDOUT "Genomic DNA or SPLK => identification of the tag_point\n";
		
		### Default race is 5'. I know race is not the most appropiate wrod, but I idnetify as race the side of the sequencing.
		if ($race eq 'na') {
			$race = "5'";
		}
		
		$debug && print STDOUT "RACE = $race\tTRAP STRAND: $trap_strand\n GENE STRAND: $gene_strand\tGene contains $exons exons.\n";
		
		my ($tag_point);
		if ($trap_strand eq "1") {
			if ($race eq "3'") {
				$tag_point = $blocks[0]->{'start'};
			} elsif ($race eq "5'") {
				$tag_point = $blocks[$last_block_num]->{'end'};
			}
		} elsif ($trap_strand eq "-1") {
			if ($race eq "3'") {
				$tag_point = $blocks[$last_block_num]->{'end'};
			} elsif ($race eq "5'") {
				$tag_point = $blocks[0]->{'start'};
			}
		}
		
		### 2) Identification of the trapped intron based on the tag_point!
		$debug && print STDOUT "Genomic DNA or SPLK => identification of the intron in which the trap insertion occurred. Based on the tag_point\n";
		
		if ($tag_point) {
			my ($pre_exon_id, $pre_exon_obj, $next_exon_id, $next_exon_obj);
			foreach my $exon (@exons) {
				my $this_exon_id = $exon->stable_id();
				my $exon_start = $exon->start();
				my $exon_end = $exon->end();
				
				if ($exon_start > $exon_end) {print STDERR "EXON START > END"; exit;}
				$debug && print STDOUT "TAG $tag_point - CURRENT EXON: $this_exon_id\n";
				
				if ($gene_strand == '1') {
					if ($tag_point >= $exon_end) {
						$pre_exon_id = $this_exon_id;
						$pre_exon_obj = $exon;
						next;
					} else {
						$next_exon_id = $this_exon_id;
						$next_exon_obj = $exon;
						last;
					}
				} elsif ($gene_strand == '-1') {
					if ($tag_point <= $exon_start) {
						$pre_exon_id = $this_exon_id;
						$pre_exon_obj = $exon;
						next;
					} else {
						$next_exon_id = $this_exon_id;
						$next_exon_obj = $exon;
						last;
					}
				}
			}
			
			$debug && print STDOUT "TAG-POINT: $tag_point\tPREV-EXON:$pre_exon_id\tNEXT-EXON:$next_exon_id\n";
			
			if ($pre_exon_id) {
				$insert {'exon_id'} = $pre_exon_id;
			}
			
			if ($next_exon_id) {
				$insert {'flanking_exon_id'} = $next_exon_id;
			}
			
			if ($pre_exon_id && $next_exon_id) {
				my $intron;
				if ($gene_strand == '1') {
					$intron = new Bio::EnsEMBL::Intron ($pre_exon_obj,$next_exon_obj);
				} elsif ($gene_strand == '-1') {
					$intron = new Bio::EnsEMBL::Intron ($next_exon_obj,$pre_exon_obj);
				}
				
				my $intron_start = $intron->start();
				my $intron_end = $intron->end();
				
				if ($intron_start < $intron_end) {
					 $insert {'putative_insertion_start'} = $intron_start;
					 $insert {'putative_insertion_end'} = $intron_end;
				} else {
					 $insert {'putative_insertion_start'} = $intron_end;
					 $insert {'putative_insertion_end'} = $intron_start;
				}
				
				$insert {'insertion_ambiguous'} = 0;
				$debug && print STDOUT "The insertion is in the INTRON between $pre_exon_id & $next_exon_id - $intron_start..$intron_end\n";
			} else {
				### We are in the case in which we don't hvae both pre- and next- exons.
				$insert {'insertion_ambiguous'} = 1;
			}
		} else {
			print STDERR "No tag_point.\n";
			exit;
		}
	} else {
		print STDERR "No block!?!\n";
		exit;
	}
	
	if ($gene_strand == 1) {
		$insert {'ens_strand'} = '+';
	} elsif ($gene_strand == -1) {
		$insert {'ens_strand'} = '-';
	}
	
	if ($gene->external_name) {
		$insert {'name'} = $gene->external_name;
	}
	
	if ($gene->description) {
		$insert {'description'} = $gene->description;
	}
	
	return (\%insert);
}

sub unique_array () {
	my ($list, $debug) = @_;
	$debug && print STDOUT "Creating a UNIQUE array!\n";
	
	my @list = @{$list};
	my %seen = ();
	my @uniq = ();
	
	foreach my $item (@list) {
		if ($item ne '') {
			unless ($seen{$item}) {
				# if we get here, we have not seen it before
				$seen{$item} = 1;
				push(@uniq, $item);
			}
		}
	}

	if (@uniq) {
		return \@uniq;
	} else {
		return undef;
	}
}

sub check_other_annotation {
    my ($sql, $source) = @_;
    my $sth = $trapdb->prepare($sql);
    my $retr = $sth->execute;
    my @array;
    my %toinsert;
    my $b =0;
    
    while (my $rhref = $sth->fetchrow_hashref) {
		my $one = 0;
		
		unless ($one == 1) {
			if ($rhref->{'display_name'}) {
				push @array,  $rhref->{'display_name'};
				$one =1;
				$b =1;
			}
		}
		
		$toinsert{'trap_id'} = $rhref->{'trap_id'};
		$toinsert{'trapmap_id'} = $rhref->{'trapmap_id'};
    }
    
    $toinsert{'annotation_ambiguous'} = 1;
    if (@array) {
	    my $unique_array = &unique_array (\@array, $debug);
	    my @unique_array = @{$unique_array};
	    print Dumper @unique_array;
	    
	    if (@unique_array) {
		    @array = @unique_array;
	    }
    }
    
    $toinsert{'ensembl_id'} = join(' - ', @array);
    $toinsert{'source'} = $source;
    
    if ($b == 1) {
		return (\%toinsert);
    } else {
		return undef;
    }
}

sub get_next_accession () {
	my ($trapdb, $debug) = @_;
	my $uni_accession;
	
	my $sql_acc = "select unitrap.accession from unitrap order by unitrap_id desc limit 1;";
	$debug && print STDOUT "SQL CODE: $sql_acc\n";
	my $sth_acc = $trapdb->prepare($sql_acc);
	$sth_acc->execute() || die "insert failed : $DBI::errstr";
	
	my @last_uni_acc = @{$sth_acc->fetchall_arrayref};
	print STDOUT Dumper @last_uni_acc;
	
	my $last = $last_uni_acc[0][0];
		
	if (!$last) {
		$uni_accession = "UNI1";
	} else {
	    my @numbers = split (/[A-Za-z]+/, $last);
		
		my $number = $numbers[1];
		my $next_number = $number + 1;
		
		$uni_accession = "UNI".$next_number;
	}
	
	$debug && print STDOUT "LAST unitrap accession: $last - NEW unitrap accession: $uni_accession\n";
	
	return $uni_accession;
}

### A refseq_dna ID is preferred to a refseq_peptide, refseq predictions are discarded!
sub get_refseq_by_gene_stable_id () {
        my ($stable_id) = @_;
		
		#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($stable_id);
		my $geneadaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
		my $gene = $geneadaptor->fetch_by_stable_id ($stable_id);
        
        my @dblinks = @{$gene->get_all_DBLinks()};
	
        my $refseq_id;
	
        foreach my $link (@dblinks) {
                $debug && print STDOUT "DB ".$link->database."\n";
		
                if ($link->database eq 'RefSeq_dna') {
                        $refseq_id = $link->primary_id;
                        if ($refseq_id =~ /NM/) {
				last;
			}
                } elsif ($link->database eq 'RefSeq_peptide') {
                        $refseq_id = $link->primary_id;
                }
        }
	
        $debug && print STDOUT "GENE $stable_id - REFSEQ $refseq_id\n";
	
        if ($refseq_id) {
                return $refseq_id;
        } else {
                return undef;
        }
}

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

sleep 5;

END {
    if ($@) {
		$error .= $@;
    }
    
    if ($DBI::errstr) {
		$error .= $DBI::errstr;
    }
    
    if ($!) {
		$error .= $!;
    }

    unless ($status) {
		print STDERR "FAILED!\n";
		print STDERR " trap id  $tid my contain fatal errors!!!!\n";

		my %error;
		$error{'trap_id'} = $tid ;
		$error{'error_note'} = $error;
		my $stmt = &prepare_stmt($trapdb, \%error);
		my $id = &insert_set($trapdb,$stmt, 'error'); 
    }
    $trapdb->disconnect;
}
