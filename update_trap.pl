#!/usr/bin/perl -w
#
# Guglielmo Roma
# guglielmoroma@libero.it
#
#   modified by Vincenza Maselli
#   v.maselli@cancer.ucl.ac.uk
#
# 1) This script executes a query in GenBank and retrieves the trap data in order to update the "trap" table
#

use strict;

BEGIN {
    require "unitrap_conf.pl";
};

use DBI;
use Getopt::Long;
use Data::Dumper;
use Net::SMTP;

use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::Tools::RepeatMasker;

use Bio::Unitrap::Db;

my $USAGE = "update_trap.pl [-c check] [-mindate] [-h help]";
my ($check, $mindate, $help);

&GetOptions(    	'check|c'			        => \$check,
			                    'mindate|m=s'    => \$mindate,
			                    'output|o'            => \$output,
			                    'help|h'           		=> \$help);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime (time); 
$year = $year+1900;
$mon = $mon+1;

if (!$mindate) {
	print OUT "You must specify a minimum date to retrieve from.\n";
	exit;
}

if ($help) {
	print OUT "HELP: 
	      -c = check y if ... n if ...
	      -mindate = minimum date to retrieve from
	      -o = the complete path of the output file, optional, only if debug is required (see configuration file to set up the debug)
	      example: nohup perl update_trap.pl -c yes -mindate 2000/01/01 -o ~/unitrap/update_trap.out &\n";
	exit;
}

#### Configuration options
my %conf =  %::conf;
my $debug = $conf{'debug'};
my $debugSQL = $conf{'debugSQL'};
my $traphost = $conf{'traphost'};
my $trapuser = $conf{'trapuser'};
my $trapdbname = $conf{'trapdbname'};
my $trappass = $conf{'trappass'};

if (!$output) {
	print OUT "The output will be write into the default directory $ENV{'HOME'}/tmp/unitrap/.\n" id $debug;
}

my $fileout = "$ENV{'HOME'}/tmp/unitrap/";
if ($output){$fileout = $output}
open (OUT, ">>$fileout") || die "$! $fileout\n"; 



$debug && print OUT "Trying to download new traps from gss db\n";


#### Dumping UniTrap structure (with data)

my $unitrapdb_obj = Bio::Unitrap::Db;

#$unitrapdb_obj->exec_dump($conf{'mysql_path'}, $conf{'trapuser'}, $conf{'trappass'}, $conf{'trapdbname'}, $conf{'traphost'}, $conf{'tmp_dir'}.$conf{'trapdbname'}."_"."$mday"."_"."$mon"."_"."$year".".sql", "0", $debug);

#### Connecting to unitrap_db
my $trap_db = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect.";




#### Reminder: Tigem traps available in Gss are not updated, so we can retrieve them from the last UNITRAP database

#### Retrieving all the gene-trap projects from the trapdb, except tigem

#my $sql_projects qq{select project_id, project_name, wordkey from project where project_name != 'tigem' and project_id=12 order by project_id};
my $sql_projects qq{select project_id, project_name, wordkey from project where project_name != 'tigem' order by project_id};

$debug && print OUT "SQL CODE: $sql_projects\n";
my $sth_projects = $trap_db->prepare($sql_projects);
my $n_projects = $sth_projects->execute;
$debug && print OUT "Number of projects in ".$conf{'trapdbname'}.": $n_projects\n";

my $base_keyword = 'gene trap and mouse AND GSS and ';
#my $base_keyword = 'gene trap and mouse AND gbdiv_gss[PROP] and ';

while (my $rhref_projects = $sth_projects->fetchrow_hashref) {
	my $p_id = $rhref_projects->{'project_id'};
	my $p_name = $rhref_projects->{'project_name'};
	my $wordkey = $rhref_projects->{'wordkey'};
	print OUT "wordkey: $wordkey\n";
	
	#### Querying GenBank
	my $seqio = &genbank_exec($base_keyword.$wordkey, $mindate);
	while (my $seq = $seqio->next_seq) {
		############################ 1 - Retrieving trap info from NCBI GenBank #############
		#### Trap_name is supposed to be the 1st word in the description
		my @desc = split (" ", $seq->desc());
		my $trap_name = $desc[0];
		
		#### Building a hash containing trap info
		my %trapinsert;
		$trapinsert{'trap_name'}= $trap_name;
		
		$debug && print Dumper @desc;
		$debug && print OUT "Project: $p_name\tTrap_name: $trap_name\n";
		
		##### Some traps from egtc cannot be parsed correctly from GSS. Need some corrections
		if ($trap_name eq 'Mus' && $p_name eq 'egtc') {
			$trap_name = $desc[8];
			if ($trap_name eq 'genomic') {
				$trap_name = $desc[7];
			}
			
			$debug && print OUT "EGTC exception => TRAP NAME $trap_name\n";
			$trap_name =~ s/,//;
			$debug && print "EGTC exception => NEW NAME $trap_name\n";
			
			if ($trap_name =~ ":") {
				my @words = split(":", $trap_name);
				$trap_name = $words[1];
			}
			
			$trapinsert{'trap_name'}= $trap_name;
		}
		
    ##### Some traps from tigm cannot be parsed correctly from GSS. Need some corrections
    if ($trap_name eq 'Mus' && $p_name eq 'tigm') {
			$trap_name = $desc[3];
		
			$debug && print OUT "TIGM exception => TRAP NAME $trap_name\n";
			$trap_name =~ s/,//;
			$debug && print "TIGM exception => NEW NAME $trap_name\n";
			
			$trapinsert{'trap_name'}= $trap_name;
		}

		my $trap_id;
		my $exists = 0;
		my $sth_check;
		
		if ($check == 1) {
			my $sql_check qq{select distinct trap_id from trap where trap_name = ?};
			$debug && print OUT "SQL CODE: $sql_check\n";
			$sth_check = $trap_db->prepare($sql_check);
			$exists = $sth_check->execute($trap_name) || die "select failed : $DBI::errstr";
		}
		
		my $splk=0;

		if ($exists == 0) {
			#### Extracting the race-type (3' or 5')
			my $race;
			#$seq is Bio::Seq object; $seq->annotation method returns a Bio::AnnotationCollectionI
			# get_Annotation returns a list of Bio::AnnotationI
			foreach my $value ($seq->annotation->get_Annotations("comment")) {
				my $comment = $value->as_text;
				# NOTE THAT $tag_name is not used in this loop
				my $tag_name = $value->tagname;
				
				if ($comment =~ /5'/) {
					$race = "5'";
				} elsif ($comment =~ /3'/) {
					$race = "3'";
				} else {
					$race = "na";
				}
				
				if ($comment =~ /Splinkerette/) {
					$splk = 1;
				}
				
				$debug && print OUT "Comment: $comment\nRACE $race\n";
			}
			
			my $gb_id = $seq->primary_id;
			my $gb_locus = $seq->accession_number;
			my $sequence =  $seq->seq;
			my $seq_length = $seq->length;
			my $nrepeat = $seq->seq =~ tr/N/N/;
			my $percent_masked = ($nrepeat/$seq_length)*100;
			
			$debug && print OUT "SEQUENCE: $sequence\nGB-ID $gb_id\nGB-LOCUS $gb_locus\n\nSEQ LENGTH $seq_length\nN-REPEAT $nrepeat\n%MASKED $percent_masked\n";
			
			#### checking for sequence quality
			my $seq_length_not_N = $seq_length - $nrepeat;
			
			my $max_frag_length_N_splitted;
			 
			# check if there are repeated regions
			# if yes, check the longer non repeated fragment
			if ($seq_length_not_N < $seq_length) {
				my @fragments = split (/N+/, $sequence);
				foreach my $fragment (@fragments) {
					$debug && print OUT "Fragment is $fragment\n";
					my $frag_length = length($fragment);
					# check for the maximum length of non repated fragments
					if ($frag_length > $max_frag_length_N_splitted) {
						$max_frag_length_N_splitted = $frag_length;
					}
				}
			} 
			else {
				$max_frag_length_N_splitted = $seq_length;
			}
			
			$debug && print OUT "Length not-N: $seq_length_not_N\nMax Fragment Length (splitted by N): $max_frag_length_N_splitted\n";
			
			$trapinsert{'project_id'} = $p_id;
			$trapinsert{'sequence'} = $sequence;
			$trapinsert{'seq_length'} = $seq_length;
			$trapinsert{'nrepeat'} = $nrepeat;
			$trapinsert{'seq_length_not_N'}= $seq_length_not_N;
			$trapinsert{'max_frag_length_N_splitted'}  = $max_frag_length_N_splitted;
			$trapinsert{'gb_id'} = $gb_id;
			$trapinsert{'gb_locus'} = $gb_locus;
			$trapinsert{'race'} = $race;
			$trapinsert{'n_percent_masked'} = $percent_masked;
			$trapinsert{'splk'} = $splk;
			
			my @features = $seq->get_all_SeqFeatures;
			
			foreach my $feat (@features) {
				$debug && print OUT "TAG ".$feat->primary_tag."\n";
				
				### Collecting other trap info, such as vector_name and mol_type
				if ($feat->primary_tag eq 'source') {
					if ($feat->has_tag('note')) {
						#### Extracting the mol_type
						my @mol = $feat->get_tag_values('mol_type');
						my $mol_type = join (",", @mol);
						$debug && print OUT "mol_type: $mol_type\n";
						
						#### Extracting the vector_name
						my @notes = $feat->get_tag_values('note');
						my $str = join (",", @notes). "\n";
						$debug && print OUT "NOTE: $str\n";
						
						my ($vec, $vector_desc) = split (":", $str);
						my ($vector_name, $other) = split (";", $vector_desc);
						$vector_name =~ s/ //g;
						$vector_name =~ s/\n//g;
						$debug && print OUT "vector name: $vector_name\n";
						
						### The way to parse the GenBank format is different for the egtc project, i.e. trap_name needs to be modified
						if ($p_name eq 'egtc') {
							my @clone = $feat->get_tag_values('clone');
							$trap_name = $clone[0];
							$trapinsert{'trap_name'}= $trap_name;
						}
						
						$trapinsert{'vector_name'} = $vector_name;
						$trapinsert{'mol_type'} = $mol_type;
						$debug && print OUT Dumper \%trapinsert;
					}
				}
			}
			
			# Inserting trap info into the trap_db
			my $stmt = $unitrapdb_obj->prepare_stmt($trap_db, \%trapinsert);
			$trap_id = $unitrapdb_obj->insert_set($trap_db, $stmt, "trap", $debugSQL);
		} 
		# trap already exists
		else {
			my $rhref_trap_id = $sth_check->fetchrow_hashref;
			$trap_id = $rhref_trap_id->{'trap_id'};
			$debug && print OUT "### Message => This trap already exists: $trap_id!!!\n";
		}
	}
}

###################
#   subroutines   #
###################

# Execute a query to 