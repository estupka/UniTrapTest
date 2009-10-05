#! /usr/bin/perl -w 
#
# Guglielmo Roma
# roma@tigem.it
#
# 2) This script aligns the sequence-tags to the mouse genome,
# retrieves the relevant annotation of the trapped region from
# the Ensembl database
# 

$|=1;
use strict;
BEGIN {
    require "/home/roma/src/scripts/unitrap/unitrap_conf.pl";
};

use DBI;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Root::IO;
use Bio::SearchIO;
use Bio::Tools::RepeatMasker;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::GeneAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::SeqFeature::Collection;

use Bio::Unitrap::File;

use Bio::Tools::Primer3;

my ( $enshost, $ensuser, $enspass, $ensdbname, $ensestdbname,
     $traphost, $trapuser, $trappass, $trapdbname,
     $fantomhost, $fantomuser, $fantompass, $fantomdbname,
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
	    "fantomhost=s",   \$fantomhost,
	    "fantomuser=s",   \$fantomuser,
	    "fantomdbname=s", \$fantomdbname,
	    "fantompass=s",   \$fantompass,
	    "blastdb=s",    	\$blastdb,
	    "hit_db=s",    	\$hit_db,
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
$fantomhost = $conf{'fantom3host'};
$fantomuser = $conf{'fantom3user'};
$fantomdbname = $conf{'fantom3dbname'};
$fantompass = $conf{'fantom3pass'};
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
#exit;

#### Connecting to unitrap_db
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect: ", $DBI::errst;

#### Connecting to Ensembl core database
my $ensdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $enshost,                                                                                     
                                            -user => $ensuser,                                                                                        
                                            -pass => $enspass,                                                                                        
                                            -dbname => $ensdbname) || die "Can't connect: ", $DBI::errst;

#### Connecting to Ensembl otherfeatures database
my $ensestdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $enshost,                                                                                     
                                            -user => $ensuser,                                                                                        
                                            -pass => $enspass,                                                                                        
                                            -dbname => $ensestdbname) || die "Can't connect: ", $DBI::errst;  

#### Connecting to Fantom3 database
my $fantom3db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $fantomhost,                                                                                     
                                            -user => $fantomuser,                                                                                        
                                            -pass => $fantompass,                                                                                        
                                            -dbname => $fantomdbname) || die "Can't connect to $fantomdbname";


$trapdb->{RaiseError}=1;
$trapdb->{RaiseError}=1;
$trapdb->{RaiseError}=1;
$trapdb->{RaiseError}=1;

my $error;
my $status;

my $tid;
my $lastid = $id + $idcount - 1;

my $do_masking = 0; # Mask the sequence tag, if do_masking != 0
my $do_mapping = 1; # Map the sequence-tag and get the annotation, if do_mapping != 0

my @tids;

foreach ($tid = $id; $tid <= $lastid; $tid++) {
	my $sth = $trapdb->prepare("select trap_id,sequence,trap_name,mol_type from trap where trap_id = $tid and checked = 0");
	$debug && print STDOUT "SQL CODE: select trap_id,sequence,trap_name,mol_type from trap where trap_id = $tid and checked = 0\n";
	my $nres = $sth->execute;
	
	push (@tids, $tid);
	
	### If this trap id is in the database continue
	if ($nres == 1) {
		my $rhref = $sth->fetchrow_hashref;
		my $trap_id = $rhref->{'trap_id'};
		my $trapseq = $rhref->{'sequence'};
		my $trap_name = $rhref->{'trap_name'};
		my $mol_type = $rhref->{'mol_type'};
		my $trapseq_masked;
		
		### remove non-word charachters
		$trapseq =~ s/\W//g;
		
		$debug && print STDOUT "###### $trap_id\t";
		
		# if ($do_masking) {
			# ##########################
			# ##### Repeat Masking #####
			# ### Create tmp file
			# my $io1 = Bio::Root::IO->new(-flush=>1);
			# my ($fh1, $file1) = $io1->tempfile();
			# 
			# ### Create seqio object
			# my $out1 = Bio::SeqIO->new(-file => ">$file1", '-format' => 'Fasta');
			# my $seqobj1 = Bio::Seq->new( -display_id => $trap_id, -seq => $trapseq);
			# $out1->write_seq ($seqobj1);
			# 
			# ### Run RepeatMasker
			# system ($repeatMasker_dir."RepeatMasker -x $file1");
			# my $rpt_parser = new Bio::Tools::RepeatMasker(-file => "$file1.out");
			# 
			# while (my $rpt = $rpt_parser->next_result) {
				# my $rpt_trap_name = $rpt->seq_id;
				# my $rpt_score = $rpt->score;
				# my $rpt_start = $rpt->start;
				# my $rpt_end = $rpt->end;
				# my $rpt_tag = $rpt->primary_tag;
				# 
				# #############################################################
				# # $trapseq_masked = ???? how to retrieve the masked seq????
 				# 
				# $debug && print STDOUT "SEQ ID $rpt_trap_name\nSCORE: $rpt_score\nSTART: $rpt_start\nEND $rpt_end\nTAG $rpt_tag\n\n";
				# 
				# ### Put repeat-info into a hash
				# my %traprptinsert;
				# $traprptinsert{'trap_id'}= $trap_id;
				# $traprptinsert{'start'}  = $rpt_start;
				# $traprptinsert{'end'} = $rpt_end;
				# $traprptinsert{'score'} = $rpt_score;
				# $traprptinsert{'tag'} = $rpt_tag;
				# $debug && print STDOUT Dumper \%traprptinsert;
				# 
				# ### Insert repeat
				# my $stmt_rpt = Bio::Genetrap::Db->prepare_stmt($trapdb, \%traprptinsert);
				# my $traprpt_id = Bio::Genetrap::Db->insert_set($trapdb, $stmt_rpt, "traprpt", $debug);
			# }
		# }
		# 
		# if ($trapseq_masked ne '') {
			# $trapseq = $trapseq_masked;
		# }
		
		if ($do_mapping) {
			##########################
			##### Genome Mapping #####
			### Create tmp file
			my $io = Bio::Root::IO->new(-flush=>1);
			my ($fh, $file) = $io->tempfile();
			
			### Create seqio object
			my $out = Bio::SeqIO->new(-file => ">$file", '-format' => 'Fasta');
			my $seqobj = Bio::Seq->new( -display_id => $trap_id,-seq => $trapseq);
			my $length = $seqobj->length;
			
			### If the length of the sequence is smaller than the word size W stop
			if ($length < 11) {
				$debug && print STDERR "length of sequence $trap_id = $length is to small to handle, exiting\n";
				$error = "sequence length $length";
				exit;
			}
			
			$out->write_seq ($seqobj);
			
			### Run wublast
			my $fh2;
			$debug && print STDOUT "blastn $blastdb $file\n";
			
			#if ($mol_type eq 'genomic DNA') {
			#	open ($fh2, "$blastnexec $blastdb $file M=1 N=-1 Q=10 E=1e-05 -nogap topcomboN=1 warnings cpus 2 sort_by_totalscore span1 |")|| die $!;
			#} else {
				open ($fh2, "$blastnexec $blastdb $file M=1 N=-1 Q=3 R=3 E=1e-05 topcomboN=1 hspsepsmax=1000000 warnings cpus 2 sort_by_totalscore span1 |")|| die $!;
			#}
			
			#####################################
			##### Insert best genomic locus #####
			### Parse wublast output
			my $searchio;
			eval{$searchio = new Bio::SearchIO->new (-fh => $fh2, -format => 'blast') || die $!;};
			
			if ($@) {
				print STDERR "Could not create searchio: $@\n";
				$error = $@;
			}
			
			my $res;
			my $fi = sprintf ("%.2f",0);
			my $full;
			my @chosen_array;
			
			while (my $result = $searchio->next_result) {
				$debug && print STDOUT "### TRAP ID: ".$result->query_name."\tlength: ".$result->query_length."\tnum_hits: ".$result->num_hits."\n";
				my ($chosen,$multiple,$score,$hsps);
				
				while (my $hit = $result->next_hit) {
					$debug && print STDOUT "\tHIT==> ".$hit->name." score ".$hit->raw_score." significance: ".$hit->significance." num_hsps: ".$hit->num_hsps." frac_identical ".$hit->frac_aligned_query." ".$hit->num_unaligned_query."\n";
					$res = 1;
					
					### Cut-off of 96% identity on the hit
					#($hit->frac_identical < 0.96) && next;
					
					($hit->frac_aligned_query < $fi) && next;
					
					if ($hit->frac_aligned_query >= 1) {
						if ($chosen) {
							if ($hit->frac_identical < $score) { ### Case 1: hit having a score lower than prev hit! NOT CHOSEN!
								last;
							} elsif ($hsps == 1 && $hit->num_hsps > 1) { ### Case 2: hit having a score equal or higher than the prev hit and more hsps (i.e, exons => gene vs pseudogene). CHOSEN!
								$chosen = $hit;
								$score = $chosen->frac_identical;
								$hsps = $chosen->num_hsps;
								
								undef (@chosen_array);
								push @chosen_array, $chosen;
								last;
							} elsif ($hsps > 1 && ($hit->num_hsps == 1)) { ### Case 3: hit having a score equal or higher than the prev hit, but less hsps (i.e, exons => pseudogene vs gene). NOT CHOSEN, but inserting
								$multiple = 1;
								
								&insert_hit_into_mysql($hit, $hit_db, $trap_id, undef, $multiple);
								last;
							} elsif ($hit->significance > $chosen->significance) { ### Case 4: hit having a score equal or higher than the prev hit and a significance (pvalue) higher than the prev hit. NOT CHOSEN
								last;
							} elsif ($hit->significance < $chosen->significance) { ### Case 5: hit having a score equal or higher than the prev hit and a significance (pvalue) lower than the prev hit. CHOSEN
								$chosen = $hit;
								$score = $chosen->frac_identical;
								$hsps = $chosen->num_hsps;
								
								undef (@chosen_array);
								push @chosen_array, $chosen;
								last;
							} elsif ($hit->significance == $chosen->significance) { ### Case 6: hit having a significance (pvalue) equal than the prev hit. Both are CHOSEN!!!
								$debug && print STDOUT "We have more than one hit with frac_identical, multiple mapping";
								print STDOUT "\tHIT==> ".$chosen->name." score ".$chosen->raw_score." significance: ".$chosen->significance." num_hsps: ".$chosen->num_hsps." frac_identical ".$chosen->frac_aligned_query." ".$chosen->frac_identical." ".$chosen->num_unaligned_query."\n";
								$multiple=1;
								
								### In case of multiple mapping, more hits are chosen
								push @chosen_array, $hit;
								next;
							}
						} else {
							$chosen = $hit;
							$score = $chosen->frac_identical;
							$hsps = $chosen->num_hsps;
							$full=1;
							
							push @chosen_array, $chosen;
							next;
						}
					}
					
					$full && last;
					
					if ($hit->frac_aligned_query >= $fi) {
						if ($chosen) {
							if ($hit->frac_identical < $score) {
								last;
							} elsif (($hsps == 1) && ($hit->num_hsps > 1)) {
								$chosen = $hit;
								$score = $chosen->frac_identical;
								$hsps = $chosen->num_hsps;
								
								undef (@chosen_array);
								push @chosen_array, $chosen;
								last;
							} elsif (($hsps > 1) && ($hit->num_hsps == 1)) {
								$multiple = 1;
								
								&insert_hit_into_mysql($hit, $hit_db, $trap_id, undef, $multiple);
								last;
							} elsif ($hit->significance > $chosen->significance) {
								last;
							} elsif ($hit->significance < $chosen->significance) {
								$chosen = $hit;
								$score = $chosen->frac_identical;
								$hsps = $chosen->num_hsps;
								
								undef (@chosen_array);
								push @chosen_array, $chosen;
								last;
							}  elsif ($hit->significance == $chosen->significance) {
								$debug && print STDOUT "We have more than one hit with frac_identical, multiple mapping!";
								$debug && print STDOUT "\tHIT==> ".$chosen->name." score ".$chosen->raw_score." significance: ".$chosen->significance." num_hsps: ".$chosen->num_hsps ." frac_identical ".$chosen->frac_aligned_query." ".$chosen->frac_identical." ".$chosen->num_unaligned_query."\n";
								$multiple = 1;
								
								### In case of multiple mapping, more hits are chosen
								push @chosen_array, $hit;
								next;
							}
						} elsif ($hit->frac_aligned_query > $fi) {
							$chosen = $hit;
							$fi = $hit->frac_aligned_query;
							$score = $chosen->frac_identical;
							$hsps = $hit->num_hsps;
							
							push @chosen_array, $chosen;
							next;
						}
					}
					
					$chosen && last;
				}
				
				if ($chosen) {
					$multiple && print STDOUT "Multiple mappings: ".scalar(@chosen_array)."\n";
					
					foreach my $chosen_hit (@chosen_array) {
						my ($tm_chosen, $feats)  = &insert_hit_into_mysql ($chosen_hit, $hit_db, $trap_id,  1, $multiple);
						$debug && print STDOUT "TRAPMAP CHOSEN: $tm_chosen\nChosen is: HIT==> ".$chosen_hit->name.":".$chosen_hit->start("hit")."..".$chosen_hit->end("hit")."\n Chosen info: Score ".$chosen_hit->raw_score." Significance: ".$chosen_hit->significance." Num_hsps: ".$chosen_hit->num_hsps." Frac_identical ".$chosen_hit->frac_aligned_query." ".$chosen_hit->frac_identical." ".$chosen_hit->num_unaligned_query."\n";
					}
				} elsif ($res) {
					print STDOUT "No hits above $fi\n";
				}
			}
			
			close ($fh2);
			
			$res || print STDOUT "No hits on the genome\n";
			
			$status = 1;
			my $sth4 = $trapdb->prepare("update trap set checked = 1 where trap_id = $trap_id");
			my $cs = $sth4->execute;
			
			print STDOUT "Updated trap $trap_id out $cs\n";
		}
	} else {
		print STDERR "No results, $tid does not exist\n";
		
		my $trap_id = $tid;
		$error = "not a valid trap_id";
	}
}

############# Correction
### Discard trapmap when the percentage identity of all the trapblocks is lower than 96%
### Set chosen = 0 where the trapmap has no trapblocks
my $sql_update = "update trapmap set chosen=0 where trapmap_id not in (select distinct trapmap_id from trapblock) and trap_id in (".join(",", @tids).");";
$debug && print STDOUT "SQL CODE: $sql_update\n";
my $sth_update = $trapdb->prepare($sql_update);
$sth_update->execute;

#######################
#
#     SUBS 
#
#######################

sub insert_hit_into_mysql {
	my ($hit, $hit_db, $trap_id, $chosen, $multiple) = @_;
	
	### Store trapmap information for chosen
	my %insert;
	$insert{'trap_id'} = $trap_id;
	$insert{'hit_id'} = $hit->name;
	$insert{'hit_db'} = $hit_db;
	$insert{'start'} = $hit->start("hit");
	$insert{'end'}  = $hit->end("hit");
	$insert{'frac_aligned_query'} = $hit->frac_aligned_query;
	$insert{'num_hsps'} = $hit->num_hsps;
	$insert{'frac_identical'} = $hit->frac_identical;
	$insert{'score'} = $hit->score;
	$insert{'significance'} = $hit->significance;
	
	if ($chosen) {
		$insert{'chosen'} = 1;
	}
	
	if ($multiple) {
		$insert{'multiple'} = 1;
	}
	
	my $stmt = &prepare_stmt($trapdb, \%insert);
	my $trapmap_id = &insert_set ($trapdb, $stmt, 'trapmap');
	
	### Insert hsps
	my @features;
	
	while (my $hsp = $hit->next_hsp) {
		### Use a cut-off of 96% identity for hsps
		### Simply do not store if percent_identity >= 96
		if ($hsp->percent_identity >= 96) {
			my %store;
			$store{'start'} = $hsp->start('hit');
			$store{'end'} = $hsp->end('hit');
			$store{'trapmap_id'}= $trapmap_id; 
			$store{'perc_ident'} = $hsp->percent_identity;
			$store{'strand'}  = $hsp->strand('query');
			$store{'trap_id'} = $trap_id;
			my $stmt = &prepare_stmt($trapdb, \%store);
			my $trapblock_id = &insert_set ($trapdb, $stmt, 'trapblock');
			
			### Create a feature for each trap_block
			my $feature = Bio::SeqFeature::Generic->new(	-display_name => $trapblock_id,
									-start => $hsp->start('hit'),
									-end => $hsp->end('hit'),
									-strand => $hsp->strand('hit')	);
			push @features,$feature;
			
			$debug && print STDOUT "Trapblock inserted $trapblock_id \n";
		}
	}
	
	return ($trapmap_id, \@features);
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
