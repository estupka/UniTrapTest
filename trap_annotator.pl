#! /usr/bin/perl -w 
#
# Guglielmo Roma
# guglielmo.roma@gmail.com
#
# This script retrieves the relevant annotation of the trapped region from the Ensembl database
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
use Bio::SeqIO;
use Bio::Seq;
use Bio::Root::IO;
use Bio::SearchIO;
use Bio::Tools::RepeatMasker;
use Bio::SeqFeature::Collection;

use Bio::EnsEMBL::Registry;
use Bio::Unitrap::File;

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

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    	-host => $enshost,
        -user => $ensuser,
		-pass => $enspass,
		-database => $ensdbname
);

my $registry_est = 'Bio::EnsEMBL::Registry';
$registry_est->load_registry_from_db(
        -host => $enshost,
		-user => $ensuser,
		-pass => $enspass,
		-database => $ensestdbname
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

#### Connecting to unitrap_db
my $trapdb = DBI->connect("DBI:mysql:database=$trapdbname;host=$traphost;port=3306", $trapuser, $trappass) || die "Can't connect: ", $DBI::errst;
my $ensdb = DBI->connect("DBI:mysql:database=$ensdbname;host=$enshost;port=3306", $ensuser, $enspass) || die "Can't connect: ", $DBI::errst;
my $ensestdb = DBI->connect("DBI:mysql:database=$ensestdbname;host=$enshost;port=3306", $ensuser, $enspass) || die "Can't connect: ", $DBI::errst;
#my $fantom3db = DBI->connect("DBI:mysql:database=$fantomdbname;host=$fantomhost;port=3306", $fantomuser, $fantompass) || die "Can't connect: ", $DBI::errst;
my $fantom3db;

$trapdb->{RaiseError}=1;
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
			
			###############################
			##### Retrieve Annotation #####
			### Querying the Ensembl database
			my $t2e_status = &annotate_with_ensembl (\@features, $hit_id, $start, $end, $tid);
			
			my $sql_update = "update trapmap set checked=1 where trapmap_id = $trapmap_id;";
			$debug && print STDOUT "SQL CODE: $sql_update\n";
			my $sth_update = $trapdb->prepare($sql_update);
			$sth_update->execute;
		}
		
		$status = 1;
		

	} else {
		print STDERR "not mapped trap ($tid) \n";
		my $trap_id = $tid;
		$error = "not mapped trap";
	}
}

#######################
#
#     SUBS 
#
#######################

sub annotate_with_ensembl {
    my ($fs, $chr, $start, $end, $trap_id) = @_;
    
    ### Ensembl genes parameters (included RefSeq genes)
    my $do_ensgene = 1;
    
    ### Estgenes parameters
    my $do_ensestgene = 0;
    
    ### Genescan parameters
    my $do_genescan = 0;
    
    ### FANTOM3 parameters
    my $do_fantom = 0;
    my $fantom_coverage = 0;
    my $fantom_perc_id = 96;
    
    ### UniGene parameters
    my $do_unigene = 1;
    my $unigene_coverage = 0;
    my $unigene_perc_id = 96;
    
    ### cDNA parameters
    my $do_mouse_cdna = 1;
    my $cdna_coverage = 66;
    my $cdna_perc_id = 96;
    
    ### Ests parameters
    my $do_ensest = 1;
    my $est_coverage = 66;
    my $est_perc_id = 96;
    
    ### TCLG parameters
    my $do_tclg = 0;
    my $tclg_coverage = 0;
    
    ### Repeats parameters
    my $do_ensrepeats = 0;
    
    my $dodumper;
    my @features = @{$fs};
    my $region;
    
    if ($chr =~ /NT/) {
		$region = "supercontig";
   	} else {
		$region = "chromosome";
    }
    
    $debug && print STDOUT "Got ".scalar @features." blocks to analyze\n";
    
    ### Reset coordinates to be local
    foreach my $f (@features) {
		my $s = $f->start();
		$f->start($s-$start+1);
		my $e = $f->end();
		$f->end($e-$start+1);
	}
	
	### Get overlap with ensembl genes (both RefSeq and ENSEMBL)
	if ($do_ensgene) {
		my $hashref;
		
		#my $slice = $ensdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
		
		$debug && print STDOUT "Annotating with $ensdbname \n";
		$debug && print STDOUT "Got slice for chr ".$slice->seq_region_name." start ".$slice->start." to ".$slice->end."\n";
		$debug && print STDOUT "Trap from 1 to ".$slice->length."\n";
		
		foreach my $f (@features) {
			$debug && print STDOUT "Trapblock ".$f->start." end ".$f->end."\n";
			$hashref->{$f->display_name} = &annotate_with_ensembl_gene ($f, $slice);
		}
		
		$dodumper && print Dumper $hashref;
		$hashref = &resolve_ensembl_conflicts($hashref);
		
		### Add strand to the hash
		foreach my $f (keys %{$hashref}) {
			foreach my $g (keys %{$hashref->{$f}}){
				if ($g =~ "ENSMUS") {
					
					#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id ($g);
					my $gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
					my $gene = $gene_adaptor->fetch_by_stable_id ($g);
					
					$hashref->{$f}{$g}{'strand'} = $gene->strand;
				}
			}
		}
		
		$dodumper && print Dumper $hashref;
		my $stmt =  &prepare_sql($trap_id, $hashref,'trapblock_ensg_annotation');
	}
		
	### Get overlap with ensembl est-genes
	if ($do_ensestgene) {
		my $hashref;
		$debug && print STDOUT "Annotating with $ensestdbname \n";
		
		#my $slestgene= $ensestdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'otherfeatures', 'Slice' );
		my $slestgene = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
		
		foreach my $f (@features) {
			$debug && print STDOUT "Trapblock ".$f->start." end ".$f->end."\n";
			$hashref->{$f->display_name} = &annotate_with_ensembl_gene ($f, $slestgene);
		}
		
		$dodumper && print Dumper $hashref;
		$hashref = &resolve_ensembl_conflicts($hashref);
		
		### Add strand to the hash
		#	foreach my $f (keys %{$hashref}) {
		#		foreach my $g (keys %{$hashref->{$f}}){
		#			if ($g =~ "ENSMUS") {
		#				#my $gene = $ensestdb->get_GeneAdaptor->fetch_by_stable_id ($g);
		#				my $gene_adaptor  = $registry_est->get_adaptor( 'Mouse', 'Core', 'Gene' );
		#				my $gene = $gene_adaptor->fetch_by_stable_id ($g);
		#
		#				$hashref->{$f}{$g}{'strand'} = $gene->strand;
		#			}
		#		}
		#	}
		
		$dodumper && print Dumper $hashref;
		my $stmt =  &prepare_sql ($trap_id, $hashref, 'trapblock_ensestg_annotation');
	}
	
	### Get overlap with genescan gene predictions
	if ($do_genescan) {
		#my $slice = $ensdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
	
		$debug && print STDOUT "Annotating with $ensdbname using prediction_transcripts (i.e., genescan)\n";
		$debug && print STDOUT "Got slice for chr ".$slice->seq_region_name." start ".$slice->start." to ".$slice->end."\n";
		$debug && print STDOUT "Trap from 1 to ".$slice->length."\n";
		
		foreach my $f (@features) {
			$debug && print STDOUT "Trapblock ".$f->start." end ".$f->end."\n";
			&annotate_with_ensembl_genescan ($f, $slice, 'trapblock_abinitio_annotation', $trap_id);
		}
	}
	
	### Get overlap with Fantom3 (previoulsy mapped with our in-house pipeline with a cut-off of 96% identity)
	if ($do_fantom) {
		my $sql_fantom = "select f.name, fb.start, fb.end, fb.strand, fb.perc_ident from feature f, featuremap fm, featureblock fb where f.feature_id=fm.feature_id and fm.featuremap_id=fb.featuremap_id and fb.start <= $end and fb.end >= $start and fm.hit_id='$chr'";
		my $sth_fantom = $fantom3db->prepare($sql_fantom);
		$debug && print STDOUT "SQL CODE: $sql_fantom\n";
		$sth_fantom->execute;
		
		$debug && print STDOUT "Annotating with features from $fantomdbname using Fantom\n";
		
		my @dnadnafeat;
		while (my $fantom = $sth_fantom->fetchrow_hashref) {
			if ($fantom->{'perc_ident'} >= $fantom_perc_id) {
				### Handle coordinates to be local
				my $feature = Bio::SeqFeature::Generic->new (	-display_name => $fantom->{'name'},
										-start => $fantom->{'start'}-$start+1,
										-end => $fantom->{'end'}-$start+1,
										-strand => $fantom->{'strand'});
				push @dnadnafeat, $feature;
			}
		}
		
		my $col = new Bio::SeqFeature::Collection();
		my $totaladded = $col->add_features(\@dnadnafeat);
		
		foreach my $f (@features) {
			$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
			
			my @subset = $col->features_in_range (-range=>$f, -contain=> 0);
			
			foreach my $s (@subset) {
				$debug && print STDOUT "In range: ".$s->display_name."\t".$s->start."\t".$s->end."\n";
				
				my $posrelative = &relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
				my $inter = $f->intersection($s);
				
				$debug && print STDOUT "Intersection ".$inter->length."\n";
				my $coverage = ($inter->length/$s->length)*100;
				$debug && print STDOUT "Coverage ".$coverage."\n";
				
				if ($coverage > $fantom_coverage) {
					my %fantom_ins;
					$fantom_ins{'display_name'} = $s->display_name;
					$fantom_ins{'dS_s'} = $posrelative->{'dS_s'};
					$fantom_ins{'dE_e'} = $posrelative->{'dE_e'};
					$fantom_ins{'label'} = $posrelative->{'label'};
					$fantom_ins{'comment'} = $posrelative->{'comment'};
					$fantom_ins{'coverage'} = ($inter->length/$s->length)*100;
					$fantom_ins{'trapblock_id'} = $f->display_name;
					$fantom_ins{'trap_id'} = $trap_id;
					$fantom_ins{'strand'} = $s->strand;
					
					my $st = &prepare_stmt ($trapdb, \%fantom_ins);
					my $insid = &insert_set ($trapdb, $st, 'trapblock_fantom_annotation');
					
					$debug && print STDOUT "Inserted trapblock_fantom_annotation with id $insid\n";
				}
			}
		}
	}
	
	### Get overlap with ensembl UniGene
	if ($do_unigene) {
		#my $slice = $ensdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
	
		$debug && print STDOUT "Annotating with dnadna align features from $ensdbname using UniGene\n";
		
		my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('Unigene')};
		my @dnadnafeat;
		
		#ensembl dnaalign features is not compatible with Bio::RangeI
		foreach my $dnadna (@dna_dna_align_feats) {
			if ($dnadna->percent_id >= $unigene_perc_id && $dnadna->hseqname =~ /Mm/) {
		#		print $dnadna->hseqname . "\t".$dnadna->hstart ."\t". $dnadna->hend . "\t". $dnadna->hstrand . "\t". $dnadna->percent_id. "\n";
		#		print $dnadna->seqname . "\t". $dnadna->start ."\t". $dnadna->end . "\t". $dnadna->strand . "\t". $dnadna->score. "\n";
				
				my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
										-start => $dnadna->hstart,
										-end => $dnadna->hend,
										-strand => $dnadna->hstrand);
				push @dnadnafeat ,$feature;
			}
		}
		
		my $col = new Bio::SeqFeature::Collection();
		my $totaladded = $col->add_features(\@dnadnafeat);
		
		foreach my $f (@features) {
			$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
			
			my @subset = $col->features_in_range (-range=>$f, -contain=> 0);
			
			foreach my $s (@subset) {
				$debug && print STDOUT "In range: ".$s->display_name."\t".$s->start."\t".$s->end."\n";
				
				my $posrelative = &relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
				my $inter = $f->intersection($s);
				
				$debug && print STDOUT "Intersection ".$inter->length."\n";
				my $coverage = ($inter->length/$s->length)*100;
				$debug && print STDOUT "Coverage ".$coverage."\n";
				
				if ($coverage > $unigene_coverage) {
					my %unigene_ins;
					$unigene_ins{'display_name'} = $s->display_name;
					$unigene_ins{'dS_s'} = $posrelative->{'dS_s'};
					$unigene_ins{'dE_e'} = $posrelative->{'dE_e'};
					$unigene_ins{'label'} = $posrelative->{'label'};
					$unigene_ins{'comment'} = $posrelative->{'comment'};
					$unigene_ins{'coverage'} = ($inter->length/$s->length)*100;
					$unigene_ins{'trapblock_id'} = $f->display_name;
					$unigene_ins{'trap_id'} = $trap_id;
					$unigene_ins{'strand'} = $s->strand;
					
					my $st = &prepare_stmt ($trapdb, \%unigene_ins);
					my $insid = &insert_set ($trapdb, $st, 'trapblock_unigene_annotation');
					
					$debug && print STDOUT "Inserted trapblock_unigene_annotation with id $insid\n";
				}
			}
		}
	}
	
	### Get overlap with ensembl mouse-cdnas
	if ($do_mouse_cdna) {
		#my $slice = $ensdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
		
		$debug && print STDOUT "Annotating with dnadna align features from $ensdbname using cDNA\n";
		
		my @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures ('mouse_cdna')};
		my @dnadnafeat;
		
		#ensembl dnaalign features is not compatible with Bio::RangeI
		foreach my $dnadna (@dna_dna_align_feats) {
			if ($dnadna->percent_id >= $cdna_perc_id) {
		#		print $dnadna->hseqname . "\t".$dnadna->hstart ."\t". $dnadna->hend . "\t". $dnadna->hstrand . "\t". $dnadna->percent_id. "\n";
		#		print $dnadna->seqname . "\t". $dnadna->start ."\t". $dnadna->end . "\t". $dnadna->strand . "\t". $dnadna->score. "\n";
				
				my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
										-start => $dnadna->hstart,
										-end => $dnadna->hend,
										-strand => $dnadna->hstrand);
				push @dnadnafeat ,$feature;
			}
		}
		
		my $col = new Bio::SeqFeature::Collection();
		my $totaladded = $col->add_features(\@dnadnafeat);
		
		foreach my $f (@features) {
			$debug && print STDERR "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
			
			my @subset = $col->features_in_range (-range=>$f, -contain=> 0);
			
			foreach my $s (@subset) {
				$debug && print STDOUT "In range: ".$s->display_name."\t".$s->start."\t".$s->end."\n";
				
				my $posrelative = &relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
				my $inter = $f->intersection($s);
				
				$debug && print STDOUT "Intersection ".$inter->length."\n";
				my $coverage = ($inter->length/$s->length)*100;
				$debug && print STDOUT "Coverage ".$coverage."\n";
				
				if ($coverage > $cdna_coverage) {
					my %cdna_ins;
					$cdna_ins{'display_name'} = $s->display_name;
					$cdna_ins{'dS_s'} = $posrelative->{'dS_s'};
					$cdna_ins{'dE_e'} = $posrelative->{'dE_e'};
					$cdna_ins{'label'} = $posrelative->{'label'};
					$cdna_ins{'comment'} = $posrelative->{'comment'};
					$cdna_ins{'coverage'} = ($inter->length/$s->length)*100;
					$cdna_ins{'trapblock_id'} = $f->display_name;
					$cdna_ins{'trap_id'} = $trap_id;
					$cdna_ins{'strand'} = $s->strand;
					
					my $st = &prepare_stmt ($trapdb, \%cdna_ins);
					my $insid = &insert_set ($trapdb, $st, 'trapblock_cdna_annotation');
					
					$debug && print STDOUT "Inserted trapblock_cdna_annotation with id $insid\n";
				}
			}
		}
	}
	
	### Get overlap with ensembl mouse-ests
	if ($do_ensest) {
		$debug && print STDOUT "Annotating with dnadna align features from $ensestdbname using Ests\n";
		
		#my $slest= $ensestdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry_est->get_adaptor( 'Mouse', 'otherfeatures', 'Slice' );
		my $slest = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
	
		my @dna_dna_align_feats = @{$slest->get_all_DnaAlignFeatures ('mouse_est')};
		my @dnadnafeat;
		
		### Ensembl dnaalign features is not compatible with Bio::RangeI
		foreach my $dnadna (@dna_dna_align_feats) {
			if ($dnadna->percent_id >= $est_perc_id) {
				my $feature = Bio::SeqFeature::Generic->new (	-display_name => $dnadna->hseqname,
										-start => $dnadna->hstart,
										-end => $dnadna->hend,
										-strand => $dnadna->hstrand);
				push @dnadnafeat ,$feature;
			}
		}
		
		my $col = new Bio::SeqFeature::Collection();
		my $totaladded = $col->add_features(\@dnadnafeat);
		
		foreach my $f (@features) {
			$debug && print STDOUT "FEATURE: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
			
			my @subset = $col->features_in_range(-range=>$f, -contain=> 0);
			
			foreach my $s (@subset) {
				$debug && print STDOUT "in range: ".$s->display_name."\t".$s->start."\t".$s->end."\n";
				
				my $posrelative = &relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
				my $inter = $f->intersection($s);
				
				$debug && print STDOUT "intersection ".$inter->length."\n";
				my $coverage = ($inter->length/$s->length)*100;
				$debug && print STDOUT "Coverage ".$coverage."\n";
				
				if ($coverage > $est_coverage) {
					my %est_ins;
					$est_ins{'display_name'} = $s->display_name;
					$est_ins{'dS_s'} = $posrelative->{'dS_s'};
					$est_ins{'dE_e'} = $posrelative->{'dE_e'};
					$est_ins{'label'} = $posrelative->{'label'};
					$est_ins{'comment'} = $posrelative->{'comment'};
					$est_ins{'coverage'} = ($inter->length/$s->length)*100;;
					$est_ins{'trapblock_id'} = $f->display_name;
					$est_ins{'trap_id'} = $trap_id;
					$est_ins{'strand'} = $s->strand;
					
					my $st = &prepare_stmt ($trapdb, \%est_ins);
					my $insid= &insert_set ($trapdb, $st, 'trapblock_ensest_annotation');
					
					$debug && print STDOUT "inserted trapblock_ensest_annotation with id $insid\n";
				}
			}
		}
	}
	
	### Get overlap with TrapclusterGene (build from our in-house pipeline)
	if ($do_tclg) {
		my $sql_tclg = "select n.accession, nb.start, nb.end, nb.strand from novel_gene n, novel_genemap nm, novel_geneblock nb where n.novel_gene_id=nm.novel_gene_id and nm.novel_genemap_id=nb.novel_genemap_id and nb.start <= $end and nb.end >= $start and nm.hit_id='$chr';";
		my $sth_tclg = $trapdb->prepare($sql_tclg);
		$debug && print STDOUT "SQL CODE: $sql_tclg\n";
		$sth_tclg->execute;
		
		$debug && print STDOUT "Annotating with features from $trapdbname using TrapCLusterGene\n";
		
		my @dnadnafeat;
		while (my $tclg = $sth_tclg->fetchrow_hashref) {
			### Handle coordinates to be local
			my $feature = Bio::SeqFeature::Generic->new (	-display_name => $tclg->{'accession'},
									-start => $tclg->{'start'}-$start+1,
									-end => $tclg->{'end'}-$start+1,
									-strand => $tclg->{'strand'});
			push @dnadnafeat, $feature;
		}
		
		my $col = new Bio::SeqFeature::Collection();
		my $totaladded = $col->add_features(\@dnadnafeat);
		
		foreach my $f (@features) {
			$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
			
			my @subset = $col->features_in_range (-range=>$f, -contain=> 0);
			
			foreach my $s (@subset) {
				$debug && print STDOUT "In range: ".$s->display_name."\t".$s->start."\t".$s->end."\n";
				
				my $posrelative = &relative_positioning ($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
				my $inter = $f->intersection($s);
				
				$debug && print STDOUT "Intersection ".$inter->length."\n";
				my $coverage = ($inter->length/$s->length)*100;
				$debug && print STDOUT "Coverage ".$coverage."\n";
				
				if ($coverage > $tclg_coverage) {
					my %tclg_ins;
					$tclg_ins{'display_name'} = $s->display_name;
					$tclg_ins{'dS_s'} = $posrelative->{'dS_s'};
					$tclg_ins{'dE_e'} = $posrelative->{'dE_e'};
					$tclg_ins{'label'} = $posrelative->{'label'};
					$tclg_ins{'comment'} = $posrelative->{'comment'};
					$tclg_ins{'coverage'} = ($inter->length/$s->length)*100;
					$tclg_ins{'trapblock_id'} = $f->display_name;
					$tclg_ins{'trap_id'} = $trap_id;
					$tclg_ins{'strand'} = $s->strand;
					
					my $st = &prepare_stmt ($trapdb, \%tclg_ins);
					my $insid = &insert_set ($trapdb, $st, 'trapblock_tclg_annotation');
					
					$debug && print STDOUT "Inserted trapblock_tclg_annotation with id $insid\n";
				}
			}
		}
	}
	
	### Get overlap with repeats
	if ($do_ensrepeats) {
		$debug && print STDOUT "Checking for reapet features from $ensdbname\n";
		
		my @repeatfeat;
		#my $slice = $ensdb->get_SliceAdaptor->fetch_by_region ($region,$chr,$start,$end);
		my $slice_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Slice' );
		my $slice = $slice_adaptor->fetch_by_region ($region,$chr,$start,$end);
	
		foreach my $repeat (@{$slice->get_all_RepeatFeatures()}) {
			$debug && print STDOUT "Repeat Name: ".$repeat->seqname." - Start  : ".$repeat->start." - End:".$repeat->end."Strand:".$repeat->strand." - Score:".$repeat->score." - Hit start:".$repeat->hstart." - Hit end:".$repeat->hend." - Consensus:".$repeat->display_id."\n";
			
			### Checking only for more_complexity repeats
			if ($repeat->display_id ne 'Low_complexity' && $repeat->display_id ne 'Simple_repeat') {
				my $feature = Bio::SeqFeature::Generic->new(	-display_name => $repeat->display_id,
										-start => $repeat->hstart,
										-end => $repeat->hend,
										-strand => $repeat->hstrand);
				push @repeatfeat ,$feature;
			}
		}
		
		if (@repeatfeat) {
			my $col = new Bio::SeqFeature::Collection();
			my $totaladded = $col->add_features(\@repeatfeat);
			
			foreach my $f (@features) {
				$debug && print STDOUT "Feature: ".$f->display_name."\t".$f->start."\t".$f->end."\n";
				
				my @subset = $col->features_in_range(-range=>$f, -contain=> 0);
				
				foreach my $s (@subset) {
					$debug && print STDOUT "In range: " . $s->display_name . "\t". $s->start . "\t". $s->end . "\n";
					my $posrelative = &relative_positioning($s->display_name, $s->start, $s->end, $f->display_name, $f->start, $f->end);
					
					my $inter = $f->intersection($s);
					$debug && print STDOUT "Intersection ". $inter->length . "\n";
					
					my $coverage = ($inter->length/$s->length)*100;
					
					my %rep_ins;
					$rep_ins{'display_name'} = $s->display_name;
					$rep_ins{'dS_s'} = $posrelative->{'dS_s'};
					$rep_ins{'dE_e'} = $posrelative->{'dE_e'};
					$rep_ins{'label'} = $posrelative->{'label'};
					$rep_ins{'comment'} = $posrelative->{'comment'};
					$rep_ins{'coverage'} = ($inter->length/$s->length)*100;;
					$rep_ins{'trapblock_id'}= $f->display_name;
					$rep_ins{'trap_id'} = $trap_id;
					
					my $st = &prepare_stmt($trapdb, \%rep_ins);
					my $insid= &insert_set($trapdb, $st, 'trapblock_repeat_annotation');
					$debug && print STDOUT "Inserted trapblock_repeat_annotation with id $insid\n";
				}
			}
		} else {
			$debug && print STDOUT "No repeats\n";
		}
	}
}

sub annotate_with_ensembl_gene {
	my ($f, $slice) = @_;
	my $genecount = 0;
	my $transcriptcount = 0;
	my $hashref;
	my $yes_exon;
	my $yes_gene;
	
	foreach my $gene (@{$slice->get_all_Genes}) {
		my $gene_feature = Bio::SeqFeature::Generic->new (	-display_name => $gene->stable_id,
															-start => $gene->start,
															-end => $gene->end,
															-strand => $gene->strand);
		my $ig = $f->intersection($gene_feature);
		
		if ($ig) {
			my $iglength = $ig->length;
			$yes_gene = 1;
			my $rank = 0;
			$genecount++;
			
			$debug && print STDOUT "\tGot a intersection between gene ". $gene->stable_id . " and Trapblock ". $f->display_name . "\n";
			
			foreach my $transcript (@{$gene->get_all_Transcripts}) {
				my $transcript_feature = Bio::SeqFeature::Generic->new (-display_name => $transcript->stable_id,
																		-start => $transcript->start,
																		-end => $transcript->end,
																		-strand => $transcript->strand);
				my $it = $f->intersection($transcript_feature);
				
				if ($it) {
					$transcriptcount++;
					my $rank = 0;
					my $exoncount = 0;
					
					foreach my $exon (@{$transcript->get_all_Exons}) {
						my $exon_feature = Bio::SeqFeature::Generic->new (	-display_name => $exon->stable_id,
																			-start => $exon->start,
																			-end => $exon->end,
																			-strand => $exon->strand);
						my $i = $f->intersection($exon_feature);
						
						if ($i) {
							$yes_exon = 1;
							$exoncount++;
							
							#we should only have one exon, after all this is trapblocks, small things...
							$debug && print STDOUT "\t\tGot a intersection between exon ". $exon->stable_id . " and Trapblock ". $f->display_name . "\n";
							$debug && print STDOUT "\t\tExon start: ".$exon->start." end: ".$exon->end."\n";
							$debug && print STDOUT "\t\tTrapblock  start: ".$f->start." end: ".$f->end."\n";
							$debug && print STDOUT "\t\tIntersection start: ".$i->start." end: ".$i->end." length " .  $i->length . "\n";
							$debug && print STDOUT "\t\tExon length: ".$exon->length." and TrapBlock length: ". $f->length. "\n";
							
							my $ilength = $i->length;
							my $ecoverage = ($ilength/$exon->length)*100;
							$debug && print STDOUT "\t \tthe exon " . $exon->stable_id. " coverage by trapblock ". $f->display_name ." is $ecoverage\n"; 
							
							$hashref->{$gene->stable_id}{$transcript->stable_id}{$exon->stable_id}{'coverage'} = $ecoverage;
							
							my $posrelative = &relative_positioning ($exon->stable_id, $exon->start, $exon->end, $f->display_name, $f->start, $f->end);
							$hashref->{$gene->stable_id}{$transcript->stable_id}{$exon->stable_id}{'relative'} = $posrelative;
							
							#rank the exons, give points for coverage and sharing of splice sites
							$rank  += $ecoverage;
							
							if ($hashref->{$gene->stable_id}{$transcript->stable_id}{$exon->stable_id}{'relative'}{'dE_e'} == 0) {
								$rank += 100;
							}
							
							if ($hashref->{$gene->stable_id}{$transcript->stable_id}{$exon->stable_id}{'relative'}{'dS_s'} == 0) {
								$rank += 100;
							}
							
							$hashref->{$gene->stable_id}{$transcript->stable_id}{'rank'} = $rank;
							
							my $transcoverage = ($ilength/$transcript->length)*100;
							$hashref->{$gene->stable_id}{$transcript->stable_id}{'coverage'} = $transcoverage * 100;
						}
					}
					
					if ($exoncount > 1) {
						$debug && print STDOUT "\t\tgot more then one exon from the same transcript... exit\n";
						#should exit...
					}
				}
			}
			
			unless ($yes_exon) {
				$iglength += 1000; #for the perl issue to solve!!!
				$debug && print STDOUT "\t\tgot zero exons, can't say which transcript is correct, not storing\n";
				
				my $genecoverage = ($iglength/$gene->length)*100;
				$debug && print STDOUT "\t\tNOTE: the intersection with the gene is $iglength and the coverage is $genecoverage\n";
				
				$hashref->{$gene->stable_id}{'no_transcript'} = $genecoverage;
			}
		}
	}
	
	unless ($yes_gene) {
		$debug && print STDOUT "\tno genes in this slice!!\n";
		
		$hashref->{'no_gene'}++;
		return $hashref;
	}
	
	return $hashref; 
}

sub resolve_ensembl_conflicts {
	my $hashref = shift;
	
	#resolve transcript conflicts by the exon rank
	my (%seengene, %addedranks, %addedtrans_coverage, %genecoverage);
	
	foreach my $f (keys %{$hashref}) {
		foreach my $g (keys %{$hashref->{$f}}){
			last if ($g eq 'no_gene');
			
			$seengene{$g}++;
			
			foreach my $t (keys %{$hashref->{$f}{$g}}) {
				if ($t eq 'no_transcript') {
					my $v = $hashref->{$f}{$g}{'no_transcript'};
					$genecoverage{$g} = $v;
					
					last;
				} else {
					my $r = $hashref->{$f}{$g}{$t}{'rank'};
					
					$addedranks{$g} += $r;
					$addedtrans_coverage{$t} += $hashref->{$f}{$g}{$t}{'coverage'};
				}
			}
		}
    }
    
    if (scalar (keys %seengene) > 1) {
		my $mostseen = &check_if_topexists (\%seengene);
		my $highestrank = &check_if_topexists (\%addedranks);
		
		if ($mostseen && $highestrank) {
			if ($highestrank  eq $mostseen) {	
				$debug && print STDERR "\t\tThe most seen IS the highest ranking\n";
				
				foreach my $f (keys %{$hashref}) {
					foreach my $g (keys %{$hashref->{$f}}) {
						last if ($g eq 'no_gene');
						
						if ($g ne $mostseen) {
							$debug && print STDERR "\t\t deleting $g\n";
							delete $hashref->{$f}{$g};
						} else {
							my $ct;
							
							foreach my $t (keys %{$hashref->{$f}{$g}}) {
								$ct++;
							}
							
							if ($ct > 1) {
								my $toptrans_coverage = &check_if_topexists(\%addedtrans_coverage);
								if ($toptrans_coverage) {
									foreach my $t (keys %{$hashref->{$f}{$g}}) {
										if ($t ne $toptrans_coverage) {
											next if ($t eq 'no_transcript');
											
											$debug && print STDOUT "\t\t deleting $t\n";
											delete  $hashref->{$f}{$g}{$t};
										}
									}
								}
							}
						}
					}
				}
			} else {
				$debug && print STDOUT "\t\tThe most seen ($mostseen)  is not the  the highestranking ($highestrank) using the rank only\n";
				
				foreach my $f (keys %{$hashref}) {
					foreach my $g (keys %{$hashref->{$f}}) {
						last if ($g eq 'no_gene');
						
						if ($g ne $highestrank) {
							$debug && print STDERR "\t\t deleting $g\n";
							delete $hashref->{$f}{$g};
						} else {
							my $ct;
							
							foreach my $t (keys %{$hashref->{$f}{$g}}) {
								$ct++;
							}
							
							if ($ct > 1) {
								my $toptrans_coverage = &check_if_topexists (\%addedtrans_coverage);
								
								if ($toptrans_coverage) {
									foreach my $t (keys %{$hashref->{$f}{$g}}) {
										if ($t ne $toptrans_coverage) {
											next if ($t eq 'no_transcript');
											
											$debug && print STDOUT "\t\t deleting $t\n";
											delete  $hashref->{$f}{$g}{$t};
										}
									}
								}	
							}
						}
					}
				}
			}
		} elsif ($highestrank) {
			$debug && print STDOUT "\t\tThere is no mostseen but the highest ranking is $highestrank\n";
			
			foreach my $f (keys %{$hashref}) {
				foreach my $g (keys %{$hashref->{$f}}) {
					last if ($g eq 'no_gene');
					
					if ($g ne $highestrank) {
						$debug && print STDOUT "\t\t deleting $g\n";
						delete $hashref->{$f}{$g};
					} else {
						my $ct;
						
						foreach my $t (keys %{$hashref->{$f}{$g}}) {
							$ct++;
						}
						
						if ($ct > 1) {
							my $toptrans_coverage = &check_if_topexists (\%addedtrans_coverage);
							
							if ($toptrans_coverage) {
								foreach my $t (keys %{$hashref->{$f}{$g}}) {
									if ($t ne $toptrans_coverage) {
										next if ($t eq 'no_transcript');
										
										$debug && print STDOUT "\t\t deleting $t\n";
										delete  $hashref->{$f}{$g}{$t};
									}
								}
							}
						}
					}
				}
			}
		} elsif ($mostseen) {
			$debug && print STDOUT "\t\tI dont think this a good parameter... $mostseen... but using\n";
			
			foreach my $f (keys %{$hashref}) {
				foreach my $g (keys %{$hashref->{$f}}) {
					last if ($g eq 'no_gene');
					
					if ($g ne $mostseen) {
						$debug && print STDOUT "\t\t deleting $g\n";
						delete $hashref->{$f}{$g};
					} else {
						my $ct;
						
						foreach my $t (keys %{$hashref->{$f}{$g}}) {
							$ct++;
						}
						
						if ($ct > 1) {
							my $toptrans_coverage = &check_if_topexists (\%addedtrans_coverage);
							
							if ($toptrans_coverage) {
								foreach my $t (keys %{$hashref->{$f}{$g}}) {
									if ($t ne $toptrans_coverage) {
										next if ($t eq 'no_transcript');
										
										$debug && print STDOUT "\t\t deleting $t\n";
										delete  $hashref->{$f}{$g}{$t};
									}
								}
							}
						}
					}
				}
			}
		} else {
			$debug && print STDOUT "\t\treached here because got more than one gene with no transcripts\n using gene coverage to choose\n";
			
			my $topgene = &check_if_topexists (\%genecoverage);
			
			if ($topgene) {
				foreach my $f (keys %{$hashref}) {
					foreach my $gt (keys %{$hashref->{$f}}) {
						last if ($gt eq 'no_gene');
						
						if ($gt ne $topgene) {
							$debug && print STDOUT "\t\t deleting $gt\n";
							delete  $hashref->{$f}{$gt};
						}
					}
				}
			} else {
				$debug && print STDOUT "\t\tnothing to do, i think this is just not the right annotation, delete remaining\n";
				
				foreach my $f (keys %{$hashref}) {
					foreach my $gt (keys %{$hashref->{$f}}) {
						last if ($gt eq 'no_gene');
						
						$debug && print STDOUT "\t\t deleting $gt\n";
						delete  $hashref->{$f}{$gt};
					}
					
					$hashref->{$f}{'no_gene'}++;
				}
			}
		}
    } elsif (scalar (keys %seengene) == 1) {
		$debug && print STDOUT "\t\t only one gene, normal \n";
		
		foreach my $f (keys %{$hashref}) {
			foreach my $g (keys %{$hashref->{$f}}) {
				last if ($g eq 'no_gene');
				
				my $c;
				
				foreach my $t (keys %{$hashref->{$f}{$g}}) {
					$c++;
				}
				
				if ($c > 1) {
					my $toptrans_coverage = &check_if_topexists (\%addedtrans_coverage);
					
					if ($toptrans_coverage) {
						foreach my $t (keys %{$hashref->{$f}{$g}}) {
							if ($t ne $toptrans_coverage) {
								next if ($t eq 'no_transcript');
								
								$debug && print STDOUT "\t\t deleting $t\n";
								delete  $hashref->{$f}{$g}{$t};
							}
						}
					}
				}
			}
		}
    } else {
		$debug && print STDOUT "\t\t no genes, leave it, just to confirm \n";
		return $hashref;
    }
    
    return $hashref;
}

sub annotate_with_ensembl_genescan {
	my ($f, $slice, $table, $trap_id) = @_;
	
	foreach my $ptrans (@{$slice->get_all_PredictionTranscripts}) {
		my $ptrans_feature = Bio::SeqFeature::Generic->new (	-display_name => $ptrans->stable_id,
									-start => $ptrans->start,
									-end => $ptrans->end,
									-strand => $ptrans->strand);
		my $ig = $f->intersection($ptrans_feature);
		
		my $exon_check;
		
		foreach my $exon (@{$ptrans->get_all_Exons}) {
			my $exon_feature = Bio::SeqFeature::Generic->new (	-display_name => $exon->display_id,
										-start => $exon->start,
										-end => $exon->end,
										-strand => $exon->strand);
			my $i = $f->intersection($exon_feature);
			
			if ($i) {
				### We should have only one exon, after all this is trapblock, small thing...
				$debug && print STDOUT "\t\tGot a intersection between exon ".$exon->display_id." and Trapclusterblock ".$f->display_name."\n";
				$debug && print STDOUT "\t\tExon start: ".$exon->start." end: ".$exon->end."\n";
				$debug && print STDOUT "\t\tTrapclusterblock  start: ".$f->start." end: ".$f->end."\n";
				$debug && print STDOUT "\t\tIntersection start: ".$i->start." end: ".$i->end." length ".$i->length."\n";
				$debug && print STDOUT "\t\tExon length: ".$exon->length." and TrapclusterBlock length: ".$f->length."\n";
				
				print STDOUT "Genescan Transcript ".$ptrans->display_id."\tExon ".$exon->display_id."EXON strand: ".$exon->strand."\n";
				
				### Calculate the exact overlap between the features (i.e., the trapblock and the overlapped exon) to choose the right trapped gene
				my $posrelative = &relative_positioning ($exon->display_id, $exon->start, $exon->end, $f->display_name, $f->start, $f->end);
				
				my $ilength = $i->length;
				my $ecoverage = ($ilength/$exon->length)*100;
				$debug && print STDOUT "\t \tthe exon " . $exon->display_id. " coverage by trapblock ". $f->display_name ." is $ecoverage\n"; 
				
				my %ins;
				$ins{'display_name'} = $ptrans->display_id;
				$ins{'exon_name'} = $exon->display_id;
				$ins{'trapblock_id'} = $f->display_name;
				$ins{'trap_id'} = $trap_id;
				$ins{'exon_coverage'} = $ecoverage;
				$ins{'dS_s'} = $posrelative->{'dS_s'};
				$ins{'dE_e'} = $posrelative->{'dE_e'};
				$ins{'label'} = $posrelative->{'label'};
				$ins{'comment'} = $posrelative->{'comment'};
				$ins{'type'} = $ptrans->analysis->logic_name();
				$ins{'strand'} = $ptrans->strand;
				
				my $st = &prepare_stmt ($trapdb, \%ins);
				my $insid = &insert_set ($trapdb, $st, $table);
				$exon_check = 1;
			}
		}
		
		unless ($exon_check) {
			my %ins;
			$ins{'display_name'} = $ptrans->display_id;
			$ins{'trapblock_id'} = $f->display_name;
			$ins{'trap_id'} = $trap_id;
			$ins{'label'} = 'intronic';
			$ins{'strand'} = $ptrans->strand;
			$ins{'type'} = $ptrans->analysis->logic_name();
			
			my $st = &prepare_stmt ($trapdb, \%ins);
			my $insid = &insert_set ($trapdb, $st, $table);
			$debug && print STDOUT "Inserted $table with id $insid\n";
		}
	}
}

sub relative_positioning {
    my ($f1_name, $f1_start, $f1_end, $f2_name, $f2_start, $f2_end) =@_;
    #f1_start may be negative, start of coord system is relative to the slice,  not the genome or this trapblock
    my $dE_e = $f1_end - $f2_end; 
    my $dS_s = $f1_start - $f2_start;
    my %hash;
    $hash{'dE_e'} = $dE_e;
    $hash{'dS_s'} = $dS_s;
    $hash{'tb_start'} = $f2_start;
    $hash{'tb_end'} = $f2_end;
    $hash{'feat_start'}=  $f1_start;
    $hash{'feat_end'} = $f1_end;
    
    #the impossible cases
    if ($f1_start > $f1_end || $f2_start > $f2_end || $f1_end < $f2_start){
	$debug && print STDERR "your setting are not valid \tEXIT\n";
	$debug && print STDERR "exon start: $f1_start ; exon end : $f1_end ; tb start: $f2_start ; tb end: $f2_end\n";
	exit;
    }
	
    if ($f1_start == $f2_start && $f1_end == $f2_end && $dE_e == 0 && $dS_s ==0){
	#  trapblock s|-----------|e
	#  feature   S|-----------|E
	$debug && print STDOUT "\t\tPOSITION: the trapblock matchs the feature perfectly\n";
	$hash{'label'} = 'perfect';
	$hash{'comment'} = 'the trapblock matchs the feature perfectly';
    } elsif ($f1_start >  $f2_start && $f1_end < $f2_end && $dE_e < 0 && $dS_s > 0){
	#  trapblock s|-----------|e
	#  feature      S|----|E
	$debug && print STDOUT "\t\tPOSITION: trapblock contains the feature\n";
	$hash{'label'} = 'trapblock_contains_feature';
	$hash{'comment'} = 'the trapblock contains the feature';
    } elsif ($f1_start <  $f2_start && $f1_end > $f2_end && $dE_e > 0 && $dS_s < 0){
	#  f2 trapblock     s|------|e
	#  f1 feature S|-----------------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock is contained in  the feature\n";
	$hash{'label'} = 'trapblock_contained_in_feature';
	$hash{'comment'} = 'the trapblock is contained in the feature';
    } elsif ($f1_start ==  $f2_start && $f1_end > $f2_end && $dE_e > 0 && $dS_s == 0){
	#  f2 trapblock s|------------|e
	#  f1 feature   S|-----------------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock and feature share start but trapblock end is whitin feature\n";
	$hash{'label'} = 'share_start_tb_end_whitin';
	$hash{'comment'} = 'the trapblock and the feature share start, but the trapblock end is whitin feature';
    } elsif ($f1_start ==  $f2_start && $f1_end < $f2_end && $dE_e < 0 && $dS_s == 0){
	#  f2 trapblock s|------------------|e
	#  f1 feature   S|-----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock and feature share start but end of trapblock extends outside of the  feature\n";
	$hash{'label'} = 'share_start_tb_end_extends';
	$hash{'comment'} = 'the trapblock and the feature share start, but the trapblock end extends outside of the feature';
    } elsif ($f1_start >  $f2_start && $f1_end == $f2_end && $dE_e == 0 && $dS_s >  0){
	#  f2 trapblock s|------------------|e
	#  f1 feature          S|-----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock and feature share end but start of feature is whitin the  trapblock\n";
	$hash{'label'} = 'share_end_tb_start_whitin';
	$hash{'comment'} = 'the trapblock and the feature share end, but the feature start is whitin the trapblock';
    } elsif ($f1_start <  $f2_start && $f1_end == $f2_end && $dE_e == 0 && $dS_s <  0){
	#  f2 trapblock       s|------------|e
	#  f1 feature  S|-------------------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock and feature share end but start of feature extends outside the  trapblock\n";
	$hash{'label'} = 'share_end_tb_start_extends';
	$hash{'comment'} = 'the trapblock and the feature share end, but the feature start extends outside the trapblock';
    } elsif ($f1_start <  $f2_start && $f1_end < $f2_end && $dE_e < 0 && $dS_s <  0){
	#  f2 trapblock       s|------------|e
	#  f1 feature  S|----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock overlaps feature at start of trapblock\n";
	$hash{'label'} = 'overlap in tb start';
	$hash{'comment'} = 'the trapblock overlaps the feature at start of trapblock';
    } elsif ($f1_start >  $f2_start && $f1_end > $f2_end && $dE_e > 0 && $dS_s >  0){
	#  f2 trapblock       s|------------|e
	#  f1 feature               S|----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock overlaps feature at end of trapblock\n";
	$hash{'label'} = 'overlap at tb end';
	$hash{'comment'} = 'the trapblock overlaps the feature at end of trapblock';
     } elsif ($f1_end == $f2_start){
	#  f2 trapblock            s|------------|e
	#  f1 feature   S|----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock start is the end of the feature \n";
	$hash{'label'} = 'tb_start_matches_feat_end';
	$hash{'comment'} = 'the trapblock start is the end of the feature';
    } elsif ($f1_start == $f2_end){
	#  f2 trapblock s|------------|e
	#  f1 feature                S|----------|E
	$debug && print STDOUT "\t\tPOSITION: trapblock end is the start of the feature\n";
	$hash{'label'} = 'tb_end_matches_feat_start';
	$hash{'comment'} = 'the trapblock end is the start of the feature';
    } else {
	$debug && print STDOUT "\t\tPOSITION: something is fishy  exit\n";
	exit;
    }
    
    return \%hash;
}

sub check_if_topexists {
	#badly written :(
	my $hashref = shift;
	my %hash = %{$hashref};
	
	if (scalar (keys %hash) == 1) {
		my $id;
		
		foreach my $r (keys %hash) {
			$id = $r;
		}
		
		$debug && print STDOUT "\t\t only one element, return it: $id\n";
		
		return $id;
	} elsif (scalar (keys %hash) > 1) {
		my @array;
		
		foreach my $i (keys %hash) {
			push @array, $hash{$i};
		}
		
		my @sorted  = sort @array;
		
		if ($sorted[0] == $sorted[1]) {
			$debug && print STDOUT "\t\t no top, first is the same as the second\n";
			
			return undef;
		} else {
			my @sorted_keys = sort {$hash{$b} <=> $hash{$a}} keys %hash;
			my $top = $sorted_keys[0];
			
			return $top;
		}
	} else {
		$debug && print STDOUT "\t\t hash is not defined ?\n";
		
		return undef;
	}
}

sub prepare_sql {
    my $c = 1;
    my ($trap_id, $hashref, $table) = @_;
    
    foreach my $f (keys %{$hashref}) {
		my %hash;
		
		#print STDERR "%%%%% $f \n";
		$hash{'trapblock_id'} = $f;
		$hash{'trap_id'} = $trap_id;
		$c++;
		
		foreach my $g (keys %{$hashref->{$f}}) {
			if ($g =~ /ENSMUS/) {
				$hash{'ensmusg'} = $g;
				$hash{'strand'} = $hashref->{$f}{$g}{'strand'};
				
				if ($g =~ /ENSMUSG/) {
					my $refseq_id = &get_refseq_by_gene_stable_id ($g);
					$debug && print STDOUT "Gene: $g - RefSeq: $refseq_id\n";
					if ($refseq_id) {
						$hash{'refseq_id'} = $refseq_id;
					}
				}
			} else {
				last;
			}
			
			foreach my $t (keys %{$hashref->{$f}{$g}}) {
				if ($t =~ /ENSMUS/) {
					$hash {'ensmust'}= $t;
				} else {
					next;
				}
				
				foreach my $e (keys %{$hashref->{$f}{$g}{$t}}) {
					if ($e =~ /ENSMUS/) {
						$hash{'ensmuse'} = $e;
						$hash{'exon_coverage'} = $hashref->{$f}{$g}{$t}{$e}{'coverage'} ;
						$hash{'dS_s'} = $hashref->{$f}{$g}{$t}{$e}{'relative'}{'dS_s'};
						$hash{'dE_e'} =  $hashref->{$f}{$g}{$t}{$e}{'relative'}{'dE_e'};
						$hash{'label'} =  $hashref->{$f}{$g}{$t}{$e}{'relative'}{'label'};
						$hash{'comment'} =  $hashref->{$f}{$g}{$t}{$e}{'relative'}{'comment'};
					} elsif($e eq 'rank') {
						$hash {'rank'}= $hashref->{$f}{$g}{$t}{'rank'};
					} else {
						next;
					}
				}
			}
			
			unless ($hash{'label'}) {
				$hash{'label'} = 'intronic';
			}
		}
		
		$debug && print Dumper \%hash;
		my $stmt = &prepare_stmt ($trapdb, \%hash);
		$debug && print STDERR $stmt;
		my $id = &insert_set ($trapdb,$stmt, $table); 
		$debug && print STDERR "inserted id : $id \n";
    }
    
    return 1;
}

### A refseq_dna ID is preferred to a refseq_peptide, refseq predictions are discarded!
sub get_refseq_by_gene_stable_id () {
        my ($stable_id) = @_;
	#my $gene = $ensdb->get_GeneAdaptor->fetch_by_stable_id($stable_id);
        my $gene_adaptor  = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
	my $gene = $gene_adaptor->fetch_by_stable_id ($stable_id);
	
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
		print STDERR " trap id  $tid may contain fatal errors!!!!\n";

		my %error;
		$error{'trap_id'} = $tid ;
		$error{'error_note'} = $error;
		my $stmt = &prepare_stmt($trapdb, \%error);
		my $id = &insert_set($trapdb,$stmt, 'error'); 
    }
    $trapdb->disconnect;
}
