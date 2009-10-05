#! /usr/bin/perl -w 
#
# Guglielmo Roma
#
# guglielmo.roma@gmail.com
#
# Calculate mutated transcripts and translations by vector insertion

$|=1;
use strict;
BEGIN {
    require "/home/roma/src/scripts/unitrap2/unitrap_conf.pl";
};

use Bio::Range;
use Bio::RangeI;

use DBI;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Root::IO;
use Bio::SearchIO;
#use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::DBSQL::GeneAdaptor;
#use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::Tools::Primer3;
use Bio::EnsEMBL::Registry;

my ( $enshost, $ensuser, $enspass, $ensdbname, 
     $idcount, $id, $traphost, $trapuser, $trappass, $trapdbname,
     $tmp_dir, $blastdb, $blastnexec, $debug);

GetOptions( 
		"enshost=s",    \$enshost,
	    "ensuser=s",    \$ensuser,
	    "enspass=s",    \$enspass,
	    "ensdbname=s",  \$ensdbname,
	    "traphost=s",   \$traphost,
	    "trapuser=s",   \$trapuser,
	    "trapdbname=s", \$trapdbname,
	    "trappass=s",   \$trappass,
	    "id=s", 		\$id,
	    "idcount=s",     \$idcount,	    
	    "debug",        \$debug,
	    );

#### Configuration options
my %conf =  %::conf;
$enshost = $conf{'enshost'};
$ensuser = $conf{'ensuser'};
$ensdbname = $conf{'ensdbname'};
$enspass = $conf{'enspass'};
$traphost = $conf{'traphost'};
$trapuser = $conf{'trapuser'};
$trapdbname = $conf{'trapdbname'};
$trappass = $conf{'trappass'};
$debug = $conf{'debug'};
$blastdb = $conf{'blastdb'};
$blastnexec = $conf{'blastnexec'};
$tmp_dir = $conf{'tmp_dir'};

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

my $error;
my $status;

########################
##### Get BATCHARG params
##########################
die print "Must provide id and idcount as BATCHARG" if scalar(@ARGV) == 0;
chomp $ARGV[0];
chomp $ARGV[1];

$id = $ARGV[0];
$idcount = $ARGV[1];

unless ($id && $idcount) {
    print STDERR " no ids \n";
    exit;
}

### utilizzo un reverse arbitrario che mi servira' per disegnare un forward sulla sequenza dell'esone
my $LacZ = 'CCAGGGTTTTCCCAGTCACGAA';
my $LacZ_obj = Bio::Seq->new(   -display_id => 'LacZ',
				-seq => $LacZ);
my $LacZ_revcomp_obj = $LacZ_obj->revcom();

# reverse complement of primer 2
my $LacZ_revcomp = $LacZ_revcomp_obj->seq();

my $sth_e = $trapdb->prepare("select distinct u.ensembl_id from unitrap u left outer join mutated_gene m on m.ensembl_id=u.ensembl_id where m.mutated_gene_id is null limit $id,$idcount");
$debug && print STDOUT "SQL CODE: select distinct u.ensembl_id from unitrap u left outer join mutated_gene m on m.ensembl_id=u.ensembl_id where m.mutated_gene_id is null limit $id,$idcount;\n";
$sth_e->execute() || die "insert failed : $DBI::errstr";

while (my $ens = $sth_e->fetchrow_hashref) {
	my $ensembl_id = $ens->{'ensembl_id'};

	my $sth = $trapdb->prepare("select distinct u.ensembl_id, u.gene_name, u.chr, u.gene_start, u.gene_end, u.gene_strand, u.gene_description, u.gene_type, u.refseq_id, count(unitrap_id) as c from unitrap u left outer join mutated_gene m on m.ensembl_id=u.ensembl_id where m.mutated_gene_id is null and u.ensembl_id = \"$ensembl_id\" group by ensembl_id");
	$debug && print STDOUT "SQL CODE: select distinct u.ensembl_id, u.gene_name, u.chr, u.gene_start, u.gene_end, u.gene_strand, u.gene_description, u.gene_type, u.refseq_id, count(unitrap_id) as c from unitrap u left outer join mutated_gene m on m.ensembl_id=u.ensembl_id where m.mutated_gene_id is null and u.ensembl_id = \"$ensembl_id\" group by ensembl_id;\n";
	$sth->execute() || die "insert failed : $DBI::errstr";
	
	while (my $hit_ids = $sth->fetchrow_hashref) {
		my $chr = $hit_ids->{'chr'};
		
		#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($ensembl_id);
		my $gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
		my $gene = $gene_adaptor->fetch_by_stable_id ($ensembl_id);
		my @transcripts = @{$gene->get_all_Transcripts()};
		
		my $transcript_num = scalar(@transcripts);
		my @exons = @{$gene->get_all_Exons()};
	
		my $gene_strand = $hit_ids->{'gene_strand'};
	
		my %toinsert;
		$toinsert{'ensembl_id'} = $ensembl_id;
	
		$toinsert{'gene_name'} = $hit_ids->{'gene_name'};
		$toinsert{'gene_description'} = $hit_ids->{'gene_description'};
		$toinsert{'gene_chr'} = $chr;
		$toinsert{'gene_start'} = $hit_ids->{'gene_start'};
		$toinsert{'gene_end'} = $hit_ids->{'gene_end'};
		$toinsert{'gene_strand'} = $gene_strand;
		$toinsert{'gene_type'} = $hit_ids->{'gene_type'};
		$toinsert{'refseq_id'} = $hit_ids->{'refseq_id'};
		$toinsert{'unitrap_num'} = $hit_ids->{'c'};
		$toinsert{'transcript_num'} = $transcript_num;
		
		my $st_gene = &prepare_stmt($trapdb, \%toinsert);
		my $mutated_gene_id = &insert_set($trapdb, $st_gene, 'mutated_gene');
			
		my $j = 0;
		
		foreach my $exon (@exons) {
			my $exon_start = $exon->start();
			my $exon_end = $exon->end();
			
			my %toinsert;
			$toinsert{'exon_id'} = $exon->stable_id;
			$toinsert{'ensembl_id'} = $ensembl_id;
			$toinsert{'exon_chr'} = $chr;
			$toinsert{'exon_start'} = $exon_start;
			$toinsert{'exon_end'} = $exon_end;
			$toinsert{'exon_strand'} = $gene_strand;
			$toinsert{'exon_num'} = $j+1;
			$toinsert{'exon_seq'} = $exon->seq->seq;
			
			my $st_exon = &prepare_stmt($trapdb, \%toinsert);
			my $mutated_gene_exon_id = &insert_set($trapdb, $st_exon, 'mutated_gene_exon');		
			
			my $exon_strand;
			
			if ($gene_strand eq '+') {
				$exon_strand = 1;
			} elsif ($gene_strand eq '-') {
				$exon_strand = -1;		
			}
			
			#### cerco il forward
			my $seq1 = $exon->seq->seq.lc($LacZ_revcomp);
			&eprimer3_run_parse_test_and_insert ($chr, $ensembl_id, $exon->stable_id, $mutated_gene_exon_id, $seq1, $tmp_dir, $blastnexec, $blastdb, $j+1, $exon_start, $exon_end, $exon_strand, $LacZ, '', 'forward', length($exon->seq->seq), $debug);
			
			#### cerco il reverse
			my $seq2 = lc($LacZ).$exon->seq->seq;
			&eprimer3_run_parse_test_and_insert ($chr, $ensembl_id, $exon->stable_id, $mutated_gene_exon_id, $seq2, $tmp_dir, $blastnexec, $blastdb, $j+1, $exon_start, $exon_end, $exon_strand, '', $LacZ, 'reverse', length($LacZ), $debug);
			
			&set_best_primers ($mutated_gene_exon_id, $debug);
			$j++;
		}
		
		&calculate_mutated_protein_foreach_transcript ($trapdb, $ensembl_id, $chr, $debug);
	
		my $sth_update = $trapdb->prepare("update mutated_gene set finished = 1 where ensembl_id = '$ensembl_id'");
		$sth_update->execute;
	}
}

#################
#  subroutines  #
#################

sub calculate_mutated_protein_foreach_transcript () {
	my ($trapdb, $ensembl_id, $chr, $debug) = @_;
	
	#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($ensembl_id);
	my $gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
	my $gene = $gene_adaptor->fetch_by_stable_id ($ensembl_id);
	my $gene_strand = $gene->strand;
	
	my @transcripts = @{$gene->get_all_Transcripts()};
	
	foreach my $transcript (@transcripts) {
		my $transcript_id = $transcript->stable_id();
		
		my @exons = @{$transcript->get_all_Exons()};
		
		if ($gene_strand == 1) {
			@exons = sort {$a->end <=> $b->end} @exons;
		} else {
			@exons = sort {$b->end <=> $a->end} @exons;
		}

		my %exon_hash;
		my $k = 0;
		
		foreach my $exon (@exons) {
			my $exon_start = $exon->start();
			my $exon_end = $exon->end();
			
			# Retrieves the portion of the transcripts peptide encoded by this exon			
		 	my $peptide = $exon->peptide($transcript)->seq;
						
			$exon_hash {$k} {'exon'} = $exon;
			$exon_hash {$k} {'pep'} = $peptide;
			$exon_hash {$k} {'transc'} = $exon->seq->seq;
			
			$debug && print STDOUT "K: $k\n";
			
			$k++;
		}
		
		my ($transcript_seq,$peptide_seq);
		my ($uni_transcript_seq,$uni_peptide_seq);
		
		my $unitraps = 0;
		my $j;
		
		for ($j = 0; $j < $k; $j++) {
			if ($exon_hash {$j+1}{'exon'}) {
				my $intron = new Bio::EnsEMBL::Intron ($exon_hash {$j}{'exon'},$exon_hash {$j+1}{'exon'});
				
				my $intron_start = $intron->start();
				my $intron_end = $intron->end();
				
				$transcript_seq .= $exon_hash {$j}{'transc'};
				$peptide_seq .= $exon_hash {$j}{'pep'};
				$uni_transcript_seq .= $exon_hash {$j}{'transc'};
				$uni_peptide_seq .= $exon_hash {$j}{'pep'};
				
				my $sql = "select distinct unitrap_id, accession from unitrap where ensembl_id = '$ensembl_id' and start <= '$intron_end' and end >= '$intron_start';";
				my $sth = $trapdb->prepare($sql);
				$debug && print STDOUT "SQL CODE: $sql\n";
				$sth->execute;
				
				while (my $uni = $sth->fetchrow_hashref) {
					my $unitrap_id = $uni->{'unitrap_id'};
					my $accession = $uni->{'accession'};
				
					$uni_transcript_seq .= " - $accession -";
					$uni_peptide_seq .= " - $accession -";
				}
				
				$uni_transcript_seq .= " ";
				$uni_peptide_seq .= " ";
				
				my %toinsert;
				$toinsert{'exon_id'} = $exon_hash {$j}{'exon'}->stable_id;
				$toinsert{'ensembl_id'} = $ensembl_id;
				$toinsert{'transcript_id'} = $transcript_id;
				$toinsert{'exon_num'} = $j+1;
				$toinsert{'peptide_seq'} = $exon_hash {$j}{'pep'};
				$toinsert{'exon_seq'} = $exon_hash {$j}{'transc'};
				
				my $st_exon = &prepare_stmt($trapdb, \%toinsert);
				my $mutated_transcript_exon_id = &insert_set($trapdb, $st_exon, 'mutated_transcript_exon');
			} else {
				$transcript_seq .= $exon_hash {$j}{'transc'};
				$peptide_seq .= $exon_hash {$j}{'pep'};
				
				$uni_transcript_seq .= $exon_hash {$j}{'transc'};
				$uni_peptide_seq .= $exon_hash {$j}{'pep'};
				
				my %toinsert;
				$toinsert{'exon_id'} = $exon_hash {$j}{'exon'}->stable_id;
				$toinsert{'ensembl_id'} = $ensembl_id;
				$toinsert{'transcript_id'} = $transcript_id;
				$toinsert{'exon_num'} = $j+1;
				$toinsert{'peptide_seq'} = $exon_hash {$j}{'pep'};
				$toinsert{'exon_seq'} = $exon_hash {$j}{'transc'};
				
				my $st_exon = &prepare_stmt($trapdb, \%toinsert);
				my $mutated_transcript_exon_id = &insert_set($trapdb, $st_exon, 'mutated_transcript_exon');		
			}
		}
		
		my %toinsert;
		$toinsert{'ensembl_id'} = $ensembl_id;
		$toinsert{'transcript_id'} = $transcript_id;
		$toinsert{'peptide_seq'} = $peptide_seq;
		$toinsert{'transcript_seq'} = $transcript_seq;
		$toinsert{'uni_peptide_seq'} = $uni_peptide_seq;
		$toinsert{'uni_transcript_seq'} = $uni_transcript_seq;
		$toinsert{'exon_num'} = $k;
		
		my $st_transcript = &prepare_stmt($trapdb, \%toinsert);
		my $mutated_transcript_id = &insert_set($trapdb, $st_transcript, 'mutated_transcript');
	}
}

sub set_best_primers () {
	my ($mutated_gene_exon_id, $debug) = @_;

	#### selecting best forward
	my $sth_f = $trapdb->prepare("select mutated_gene_exon_primer_id from mutated_gene_exon_primer where mutated_gene_exon_id = '$mutated_gene_exon_id' and type = 'forward' order by blast_hits limit 1;");
	$debug && print STDOUT "SQL CODE: select mutated_gene_exon_primer_id from mutated_gene_exon_primer where mutated_gene_exon_id = '$mutated_gene_exon_id' and type = 'forward' order by blast_hits limit 1;\n";
	$sth_f->execute;
	
	my $rhref_f = $sth_f->fetchrow_hashref;
	my $mutated_gene_exon_primer_id = $rhref_f->{'mutated_gene_exon_primer_id'};
	
	my $sql2 = "update mutated_gene_exon_primer set best_primer = '1' where mutated_gene_exon_primer_id = '$mutated_gene_exon_primer_id';";
	$debug && print STDOUT "SQL CODE: $sql2\n";
	my $sth2 = $trapdb->prepare($sql2);
	$sth2->execute() || die "insert failed : $DBI::errstr";	

	#### selecting best reverse
	my $sth_r = $trapdb->prepare("select mutated_gene_exon_primer_id from mutated_gene_exon_primer where mutated_gene_exon_id = '$mutated_gene_exon_id' and type = 'reverse' order by blast_hits limit 1;");
	$debug && print STDOUT "SQL CODE: select mutated_gene_exon_primer_id from mutated_gene_exon_primer where mutated_gene_exon_id = '$mutated_gene_exon_id' and type = 'reverse' order by blast_hits limit 1;\n";
	$sth_r->execute;
	
	my $rhref_r = $sth_r->fetchrow_hashref;
	my $mutated_gene_exon_primer_id_rev = $rhref_r->{'mutated_gene_exon_primer_id'};
	
	my $sql3 = "update mutated_gene_exon_primer set best_primer = '1' where mutated_gene_exon_primer_id = '$mutated_gene_exon_primer_id_rev';";
	$debug && print STDOUT "SQL CODE: $sql3\n";
	my $sth3 = $trapdb->prepare($sql3);
	$sth3->execute() || die "insert failed : $DBI::errstr";	
}

sub eprimer3_run_parse_test_and_insert () {
	my ($chr, $ensembl_id, $exon_id, $block_id, $sequence, $tmp_dir, $blastnexec, $blastdb, $exon_num, $block_start, $block_end, $strand, $reverse_seq, $forward_seq, $type, $ref_length, $debug) = @_;
	
	#create tmp file
	my $io = Bio::Root::IO->new(-flush=>1);
	my ($fh, $infile) = $io->tempfile();
	
	$debug && print "\n###############\n";
	print STDOUT "Checking for a $type primer on exon $exon_num: $block_start, $block_end, $strand\n";
	$debug && print "###############\nFile $infile\nSequence $sequence\n";
	
	#create seqio object
	my $out = Bio::SeqIO->new(      -file => ">$infile" , '-format' => 'Fasta');
	my $seqobj = Bio::Seq->new( 	-display_id => $ensembl_id."-".$block_id,
									-seq  => $sequence);
	
	$out->write_seq ($seqobj);
	$out->close();
	
	my $seqio = Bio::SeqIO->new(-file => $infile);
	my $seq = $seqio -> next_seq;
	
	#create tmp file
	my $io2 = Bio::Root::IO->new(-flush=>1);
	my ($fh2, $outfile) = $io2->tempfile();
	
	if ($type eq 'forward') {
		print STDOUT "eprimer3 -sequence '$infile' -outfile '$outfile' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 2 -reverseinput '$reverse_seq'";
		system ("eprimer3 -sequence '$infile' -outfile '$outfile' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 2 -reverseinput '$reverse_seq'");
		
	} elsif ($type eq 'reverse') {
		print STDOUT "eprimer3 -sequence '$infile' -outfile '$outfile' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 3 -forwardinput '$forward_seq'";
		system ("eprimer3 -sequence '$infile' -outfile '$outfile' -mintm 58 -maxtm 64 -otm 60 -minsize 20 -maxsize 24 -osize 20 -task 3 -forwardinput '$forward_seq'");
	}
	
	my %primer_info;
	
	open (FILE, $outfile);

	my $f_num = 1;
	my $r_num = 1;
	
	while (defined (my $line=<FILE>)) {
		print "Line $line\n";
		
		if ($line !~ /#/ && $line =~ /PRIMER/) {
			my $primer_info = &get_primer_info ($ensembl_id."-".$block_id, $type, $line, $blastnexec, $blastdb, $debug);
			%primer_info = %{$primer_info};
			
			if (($type eq 'forward' && $primer_info {'end'} <= $ref_length) || ($type eq 'reverse' && $primer_info {'start'} > $ref_length)) {
				my %exon_primer;
				
				$exon_primer{'ensembl_id'} = $ensembl_id;
				$exon_primer{'exon_id'} = $exon_id;
				$exon_primer{'mutated_gene_exon_id'} = $block_id;
				
				########## Generic INFO
				$exon_primer{'chr'} = $chr;
				$exon_primer{'type'} = $type;
				$exon_primer{'exon_num'} = $exon_num;
				$exon_primer{'sequence'} = $primer_info {'sequence'};
				$exon_primer{'start'} = $primer_info {'start'};
				$exon_primer{'end'} = $primer_info {'end'};
				$exon_primer{'length'} = $primer_info {'length'};
				
				$exon_primer{'gc_perc'} = $primer_info {'gc_perc'};
				$exon_primer{'tm'} = $primer_info {'tm'};
				$exon_primer{'blast_hits'} = $primer_info {'blast_hits'};
				
				if ($type eq 'forward') {
					########### FORWARD
					$exon_primer{'name'} = "FP".$ensembl_id.".".$exon_num.".".$f_num;
					
					if ($strand eq 1) {
						$exon_primer{'genomic_start'} = $block_start + $primer_info {'start'} - 1;
						$exon_primer{'genomic_end'} = $block_start + $primer_info {'end'} - 1;
					} elsif ($strand eq -1) {
						$exon_primer{'genomic_start'} = $block_end - $primer_info {'end'} + 1;
						$exon_primer{'genomic_end'} = $block_end - $primer_info {'start'} + 1;
					}
					
					$exon_primer{'genomic_strand'} = $strand;
					
					$f_num++;
				} elsif ($type eq 'reverse') {
					########### REVERSE
					$exon_primer{'name'} = "RP".$ensembl_id.".".$exon_num.".".$r_num;
					
					if ($strand eq 1) {
						$exon_primer{'genomic_start'} = $block_start + $primer_info {'start'} - 1;
						$exon_primer{'genomic_end'} = $block_start + $primer_info {'end'} - 1;
					} elsif ($strand eq -1) {
						$exon_primer{'genomic_start'} = $block_end - $primer_info {'end'} + 1;
						$exon_primer{'genomic_end'} = $block_end - $primer_info {'start'} + 1;
					}
					
					if ($strand eq 1) {
						$exon_primer{'genomic_strand'} = -1;
					} elsif ($strand eq -1) {
						$exon_primer{'genomic_strand'} = 1;
					}
					
					$r_num++;
				}
				
				my $st_exon_primer = &prepare_stmt($trapdb, \%exon_primer);
				my $exon_primer_id = &insert_set($trapdb, $st_exon_primer, 'mutated_gene_exon_primer');
			}
		}
	}
	
	close (FILE);
	
	system ("rm $infile");
	system ("rm $outfile");
}

sub readFileData () {
    my ($file, $debug) = @_;	
	my $content;
	
	if ($debug) {print "Reading the file $file\n";}
	
	open (FILE, $file);
	while (defined (my $line=<FILE>)) {$content .= $line;}
	close (FILE);
	return $content;
}

sub writeDataToFile () {
    my ($file, $data, $debug) = @_;
    
	$debug && print STDOUT "Writing in the file $file\n";
	
	open (FILE, ">$file");
	print FILE "$data";
	close (FILE);
}

sub get_primer_info () {
	my ($name, $type, $line, $blastnexec, $blastdb, $debug) = @_;
	
	print ">>> TYPE $type\n";
	
	$line =~ s/\s+/&/g;
	my @splitted = split (/([&])/, $line);
	
	#print Dumper @splitted;
	
#	my $type = $splitted[2];
	my $start = $splitted[6];
	my $length = $splitted[8];
	my $tm = $splitted[10];
	my $gc = $splitted[12];
	my $seq = $splitted[14];
	my $end = $start + $length - 1;
	
	$debug && print STDOUT "PARSED LINE: type $type, start $start, end $end, length $length, tm $tm, gc $gc, seq $seq\n\n";
	
	############# blasting primers for checking...			
	#create tmp file
	my $io = Bio::Root::IO->new(-flush=>1);
	my ($fh, $primer_file) = $io->tempfile();
	
	#create seqio object
	my $out = Bio::SeqIO->new(      -file => ">$primer_file" , '-format' => 'Fasta');
	my $seqobj = Bio::Seq->new( 	-display_id => $name."-".$type,
									-seq  => $seq);
	my $seqlength = $seqobj->length;
	
	#if the length of the sequence is smaller than the word size W stop
	if ($length < 11) {
		$debug && print STDERR "length of sequence $name = $seqlength is to small to handle, exiting\n";
		$error = "sequence length $seqlength";
		exit;
	}
	
	$out->write_seq ($seqobj);
	$out->close();	
	
	#run wublast
	my $fh2;
	
	print STDOUT "blastn $blastdb $primer_file\n";
	open ($fh2, "$blastnexec $blastdb $primer_file E=10 B=100 filter=dust W=15 wink=15 M=+1 N=-3 Q=10 -nogaps X=5 topcomboN=1 warnings cpus 2 sort_by_totalscore span1 -noseqs |") || die $!;
	my $searchio;
	
	eval{
		$searchio = new Bio::SearchIO->new (-fh => $fh2, -format => 'blast') || die $!;
	};
	
	if ($@) {
		print STDERR "could not create searchio: $@\n";
		$error = $@;
	}
	
	my $num_hits;
	
	while (my $result = $searchio->next_result) {
		$num_hits = $result->num_hits;
		
		$debug && print STDOUT "### NAME: ".$result->query_name."\tlength: ".$result->query_length."\tnum_hits: ".$result->num_hits."\n";
	}
	
	close ($fh);
	close ($fh2);

	my %primer;
	
	print STDOUT "$type => $num_hits\nSEQ: $seq\n";
	
	$primer {'start'} = $start;
	$primer {'end'} = $end;
	$primer {'length'} = $length;
	$primer {'sequence'} = $seq;
	$primer {'tm'} = $tm;
	$primer {'gc_perc'} = $gc;
	$primer {'blast_hits'} = $num_hits;
	
	return \%primer;
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

sub insert_set {
    my ($dbh, $stmt, $table_name) = @_;
    #print "here" . $stmt . "\n";
    my $s = "INSERT INTO $table_name SET  $stmt";
    print STDOUT "### doing insert  $s ###\n";
    my $sth = $dbh->prepare($s);
    $sth->execute() || warn "insert failed : $DBI::errstr";
    my $dbi = $sth->{'mysql_insertid'};
    #print STDERR "the table $table_name last inserted id is $dbi \n";
    return $dbi;
}
