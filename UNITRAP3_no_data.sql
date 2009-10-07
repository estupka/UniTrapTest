-- MySQL dump 10.11
--
-- Host: hal9000    Database: UNITRAP3
-- ------------------------------------------------------
-- Server version	5.0.72

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `error`
--

DROP TABLE IF EXISTS `error`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `error` (
  `error_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `error_note` text NOT NULL,
  PRIMARY KEY  (`error_id`),
  KEY `start_id` (`trap_id`)
) ENGINE=MyISAM AUTO_INCREMENT=9737 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `maxicluster`
--

DROP TABLE IF EXISTS `maxicluster`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `maxicluster` (
  `maxicluster_id` int(10) unsigned NOT NULL auto_increment,
  `sequence` mediumtext NOT NULL,
  `accession` varchar(40) NOT NULL default '',
  `mol_type` varchar(40) default NULL,
  PRIMARY KEY  (`maxicluster_id`),
  UNIQUE KEY `acc` (`accession`),
  KEY `mol` (`mol_type`)
) ENGINE=MyISAM AUTO_INCREMENT=21797 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `maxiclusterblock`
--

DROP TABLE IF EXISTS `maxiclusterblock`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `maxiclusterblock` (
  `maxiclusterblock_id` int(10) unsigned NOT NULL auto_increment,
  `maxiclustermap_id` int(10) NOT NULL default '0',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `sequence` mediumtext NOT NULL,
  `checked` tinyint(4) NOT NULL default '0',
  PRIMARY KEY  (`maxiclusterblock_id`),
  KEY `t` (`maxiclustermap_id`),
  KEY `c` (`start`,`end`),
  KEY `start` (`start`),
  KEY `end` (`end`)
) ENGINE=MyISAM AUTO_INCREMENT=72393 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `maxiclustermap`
--

DROP TABLE IF EXISTS `maxiclustermap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `maxiclustermap` (
  `maxiclustermap_id` int(10) unsigned NOT NULL auto_increment,
  `maxicluster_id` int(10) unsigned NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '0',
  `hit_db` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  PRIMARY KEY  (`maxiclustermap_id`),
  UNIQUE KEY `clusthit` (`maxicluster_id`,`start`,`end`,`strand`),
  KEY `hit_id` (`hit_id`),
  KEY `strand` (`strand`),
  KEY `end` (`end`),
  KEY `start` (`start`),
  KEY `h` (`hit_id`)
) ENGINE=MyISAM AUTO_INCREMENT=21797 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_gene`
--

DROP TABLE IF EXISTS `mutated_gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_gene` (
  `mutated_gene_id` int(10) unsigned NOT NULL auto_increment,
  `ensembl_id` varchar(40) NOT NULL default '',
  `gene_name` varchar(40) NOT NULL default '',
  `gene_description` varchar(255) default NULL,
  `gene_strand` enum('+','-') default NULL,
  `gene_type` varchar(40) default NULL,
  `gene_start` int(10) NOT NULL default '0',
  `gene_end` int(10) NOT NULL default '0',
  `refseq_id` varchar(40) NOT NULL default '',
  `unitrap_num` int(4) NOT NULL default '0',
  `transcript_num` int(4) NOT NULL default '0',
  `gene_chr` varchar(40) NOT NULL default '',
  `finished` tinyint(1) NOT NULL default '0',
  `word` text,
  `disease` mediumtext,
  `human_ortholog` tinyint(1) NOT NULL default '0',
  `omim` tinyint(1) NOT NULL default '0',
  `mutagenicity` varchar(20) NOT NULL default '',
  `public` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`mutated_gene_id`),
  KEY `ensembl_id` (`ensembl_id`),
  KEY `unitrap_num` (`unitrap_num`),
  KEY `trascript_num` (`transcript_num`),
  KEY `gene_name` (`gene_name`),
  KEY `gene_description` (`gene_description`),
  KEY `gene_start` (`gene_start`),
  KEY `gene_end` (`gene_end`),
  KEY `gene_strand` (`gene_strand`),
  KEY `gene_type` (`gene_type`),
  FULLTEXT KEY `word` (`word`)
) ENGINE=MyISAM AUTO_INCREMENT=12126 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_gene_exon`
--

DROP TABLE IF EXISTS `mutated_gene_exon`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_gene_exon` (
  `mutated_gene_exon_id` int(10) unsigned NOT NULL auto_increment,
  `exon_id` varchar(40) default NULL,
  `ensembl_id` varchar(40) default NULL,
  `exon_seq` mediumtext,
  `exon_num` int(4) NOT NULL default '0',
  `exon_start` int(10) NOT NULL default '0',
  `exon_end` int(10) NOT NULL default '0',
  `exon_chr` varchar(40) NOT NULL default '',
  `exon_strand` enum('+','-') default NULL,
  PRIMARY KEY  (`mutated_gene_exon_id`),
  KEY `ex` (`exon_id`),
  KEY `e` (`ensembl_id`),
  KEY `n` (`exon_num`)
) ENGINE=MyISAM AUTO_INCREMENT=178826 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_gene_exon_primer`
--

DROP TABLE IF EXISTS `mutated_gene_exon_primer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_gene_exon_primer` (
  `mutated_gene_exon_primer_id` int(10) unsigned NOT NULL auto_increment,
  `mutated_gene_exon_id` int(10) NOT NULL default '0',
  `ensembl_id` varchar(40) NOT NULL default '0',
  `exon_id` varchar(40) NOT NULL default '0',
  `exon_num` int(10) NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `sequence` mediumtext NOT NULL,
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `length` int(10) NOT NULL default '0',
  `genomic_start` int(10) NOT NULL default '0',
  `genomic_end` int(10) NOT NULL default '0',
  `genomic_strand` int(10) NOT NULL default '0',
  `gc_perc` float NOT NULL default '0',
  `tm` float NOT NULL default '0',
  `type` varchar(255) NOT NULL default '',
  `blast_hits` int(10) NOT NULL default '0',
  `chr` varchar(40) NOT NULL default '',
  `best_primer` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`mutated_gene_exon_primer_id`),
  KEY `mutated_gene_exon_id` (`mutated_gene_exon_id`),
  KEY `exon_num` (`exon_num`),
  KEY `name` (`name`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `length` (`length`),
  KEY `genomic_start` (`genomic_start`),
  KEY `genomic_end` (`genomic_end`),
  KEY `genomic_strand` (`genomic_strand`),
  KEY `gc_perc` (`gc_perc`),
  KEY `tm` (`tm`),
  KEY `type` (`type`),
  KEY `blast_hits` (`blast_hits`),
  KEY `exon_id` (`exon_id`)
) ENGINE=MyISAM AUTO_INCREMENT=993083 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_gene_fragment`
--

DROP TABLE IF EXISTS `mutated_gene_fragment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_gene_fragment` (
  `mutated_gene_fragment_id` int(10) unsigned NOT NULL auto_increment,
  `ensembl_id` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `chr` varchar(40) NOT NULL default '',
  `enzyme` varchar(20) NOT NULL default '',
  `restriction_site` varchar(20) NOT NULL default '',
  `sequence` mediumtext NOT NULL,
  PRIMARY KEY  (`mutated_gene_fragment_id`),
  KEY `ensembl_id` (`ensembl_id`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `chr` (`chr`),
  KEY `enzyme` (`enzyme`),
  KEY `rs` (`restriction_site`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_transcript`
--

DROP TABLE IF EXISTS `mutated_transcript`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_transcript` (
  `mutated_transcript_id` int(10) unsigned NOT NULL auto_increment,
  `ensembl_id` varchar(40) default NULL,
  `transcript_id` varchar(40) default NULL,
  `transcript_seq` mediumtext,
  `peptide_seq` mediumtext,
  `exon_num` int(4) NOT NULL default '0',
  `uni_transcript_seq` mediumtext,
  `uni_peptide_seq` mediumtext,
  PRIMARY KEY  (`mutated_transcript_id`),
  KEY `e` (`ensembl_id`),
  KEY `t` (`transcript_id`),
  KEY `n` (`exon_num`)
) ENGINE=MyISAM AUTO_INCREMENT=24985 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutated_transcript_exon`
--

DROP TABLE IF EXISTS `mutated_transcript_exon`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutated_transcript_exon` (
  `mutated_transcript_exon_id` int(10) unsigned NOT NULL auto_increment,
  `exon_id` varchar(40) default NULL,
  `ensembl_id` varchar(40) default NULL,
  `transcript_id` varchar(40) default NULL,
  `exon_seq` mediumtext,
  `peptide_seq` mediumtext,
  `exon_num` int(4) NOT NULL default '0',
  PRIMARY KEY  (`mutated_transcript_exon_id`),
  KEY `ex` (`exon_id`),
  KEY `e` (`ensembl_id`),
  KEY `t` (`transcript_id`),
  KEY `n` (`exon_num`)
) ENGINE=MyISAM AUTO_INCREMENT=302874 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel_gene`
--

DROP TABLE IF EXISTS `novel_gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel_gene` (
  `novel_gene_id` int(10) unsigned NOT NULL auto_increment,
  `accession` varchar(40) default NULL,
  `sequence` mediumtext,
  `type` varchar(40) NOT NULL default '',
  `link_to_ensembl` mediumtext,
  `link_to_ucsc` mediumtext,
  `overlapping` varchar(40) default NULL,
  `all_exons_more_than_40_bases` tinyint(1) NOT NULL default '0',
  `novel_geneblocks` int(10) NOT NULL default '0',
  `trapclusters` int(10) NOT NULL default '0',
  `traps` int(10) unsigned NOT NULL default '0',
  `hypertrapped` tinyint(1) NOT NULL default '0',
  `overlapping_ensembl` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`novel_gene_id`),
  KEY `a` (`accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel_geneblock`
--

DROP TABLE IF EXISTS `novel_geneblock`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel_geneblock` (
  `novel_geneblock_id` int(10) unsigned NOT NULL auto_increment,
  `novel_gene_id` int(10) NOT NULL default '0',
  `novel_genemap_id` int(10) NOT NULL default '0',
  `trapclusterblock_id` int(10) NOT NULL default '0',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  `sequence` mediumtext NOT NULL,
  `overlapping` varchar(40) default NULL,
  `intronic` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`novel_geneblock_id`),
  KEY `t` (`novel_genemap_id`),
  KEY `c` (`start`,`end`),
  KEY `overlapping` (`overlapping`),
  KEY `intronic` (`intronic`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel_genemap`
--

DROP TABLE IF EXISTS `novel_genemap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel_genemap` (
  `novel_genemap_id` int(10) unsigned NOT NULL auto_increment,
  `novel_gene_id` int(10) unsigned NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '0',
  `hit_db` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  PRIMARY KEY  (`novel_genemap_id`),
  UNIQUE KEY `clusthit` (`novel_gene_id`,`start`,`end`,`strand`),
  KEY `hit_id` (`hit_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `omim`
--

DROP TABLE IF EXISTS `omim`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `omim` (
  `omim_id` int(10) unsigned NOT NULL auto_increment,
  `ortholog_id` varchar(40) NOT NULL default '',
  `accession` varchar(40) default NULL,
  `description` varchar(255) default NULL,
  PRIMARY KEY  (`omim_id`),
  KEY `t` (`accession`,`ortholog_id`,`description`)
) ENGINE=MyISAM AUTO_INCREMENT=1872 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ortholog`
--

DROP TABLE IF EXISTS `ortholog`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ortholog` (
  `ortholog_id` int(10) unsigned NOT NULL auto_increment,
  `hs_ensembl_id` varchar(40) default NULL,
  `mm_ensembl_id` varchar(40) NOT NULL default '',
  `ort_name` varchar(40) default NULL,
  PRIMARY KEY  (`ortholog_id`),
  KEY `ensembl_id` (`hs_ensembl_id`),
  KEY `mm_ensembl_id` (`mm_ensembl_id`)
) ENGINE=MyISAM AUTO_INCREMENT=10939 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pcr`
--

DROP TABLE IF EXISTS `pcr`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pcr` (
  `pcr_id` int(10) unsigned NOT NULL auto_increment,
  `unitrap_id` int(10) NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '0',
  `forward` varchar(40) NOT NULL default '0',
  `forward_start` int(10) NOT NULL default '0',
  `forward_strand` enum('+','-') default NULL,
  `reverse` varchar(40) NOT NULL default '0',
  `gc_perc` float default NULL,
  `tm` float default NULL,
  `reverse_name` varchar(10) default NULL,
  `forward_end` int(10) NOT NULL default '0',
  PRIMARY KEY  (`pcr_id`),
  KEY `u` (`unitrap_id`),
  KEY `c` (`hit_id`),
  KEY `reverse_name` (`reverse_name`)
) ENGINE=MyISAM AUTO_INCREMENT=386390 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `project`
--

DROP TABLE IF EXISTS `project`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `project` (
  `project_id` int(10) unsigned NOT NULL auto_increment,
  `project_name` varchar(40) NOT NULL default '',
  `base_url` varchar(250) default NULL,
  `url` varchar(250) NOT NULL default '',
  `wordkey` varchar(255) NOT NULL default '',
  `private` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`project_id`),
  UNIQUE KEY `project_name` (`project_name`)
) ENGINE=MyISAM AUTO_INCREMENT=13 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_mutagenicity`
--

DROP TABLE IF EXISTS `protein_mutagenicity`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_mutagenicity` (
  `protein_mutagenicity_id` int(10) unsigned NOT NULL auto_increment,
  `unitrap_id` int(10) NOT NULL default '0',
  `protein_seq` mediumtext,
  `mutagenicity` varchar(50) default NULL,
  `ensembl_id` varchar(50) default NULL,
  `transcript` varchar(50) default NULL,
  `not_mutated_exons_num` int(10) NOT NULL default '0',
  `mutated_exons_num` int(10) NOT NULL default '0',
  `not_mutated_domains_num` int(10) NOT NULL default '0',
  `mutated_domains_num` int(10) NOT NULL default '0',
  `not_mutated_domains` mediumtext,
  `mutated_domains` mediumtext,
  `rank` tinyint(1) NOT NULL default '0',
  `longest_transcript` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`protein_mutagenicity_id`),
  KEY `u` (`unitrap_id`),
  KEY `t` (`transcript`),
  KEY `e` (`ensembl_id`),
  KEY `m` (`mutagenicity`),
  KEY `nme` (`not_mutated_exons_num`),
  KEY `me` (`mutated_exons_num`),
  KEY `nmd` (`not_mutated_domains_num`),
  KEY `md` (`mutated_domains_num`)
) ENGINE=MyISAM AUTO_INCREMENT=495654 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tcl_downstream`
--

DROP TABLE IF EXISTS `tcl_downstream`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tcl_downstream` (
  `tcl_downstream_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `type` varchar(40) NOT NULL default '',
  `chr` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  `slice_start` int(10) NOT NULL default '0',
  `slice_end` int(10) NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `ensembl_id` varchar(255) NOT NULL default '',
  `close_trapcluster_id` int(10) unsigned NOT NULL default '0',
  `refseq_id` varchar(40) NOT NULL default '',
  PRIMARY KEY  (`tcl_downstream_id`),
  KEY `tcl` (`trapcluster_id`),
  KEY `type` (`type`),
  KEY `chr` (`chr`),
  KEY `start` (`start`),
  KEY `end` (`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tcl_upstream`
--

DROP TABLE IF EXISTS `tcl_upstream`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tcl_upstream` (
  `tcl_upstream_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `type` varchar(40) NOT NULL default '',
  `chr` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  `slice_start` int(10) NOT NULL default '0',
  `slice_end` int(10) NOT NULL default '0',
  `name` varchar(255) NOT NULL default '',
  `ensembl_id` varchar(255) NOT NULL default '',
  `close_trapcluster_id` int(10) unsigned NOT NULL default '0',
  `refseq_id` varchar(40) NOT NULL default '',
  PRIMARY KEY  (`tcl_upstream_id`),
  KEY `tcl` (`trapcluster_id`),
  KEY `type` (`type`),
  KEY `chr` (`chr`),
  KEY `start` (`start`),
  KEY `end` (`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trap`
--

DROP TABLE IF EXISTS `trap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trap` (
  `trap_id` int(10) unsigned NOT NULL auto_increment,
  `trap_name` varchar(40) NOT NULL default '',
  `project_id` int(10) NOT NULL default '0',
  `sequence` mediumtext NOT NULL,
  `processing` varchar(10) default NULL,
  `clone_id` varchar(40) default NULL,
  `user` varchar(20) default NULL,
  `verified` varchar(5) default NULL,
  `quality` varchar(16) default NULL,
  `timestp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
  `vector_name` varchar(40) default NULL,
  `cell_line` varchar(40) default NULL,
  `vector_type` varchar(10) default NULL,
  `mapping_type` tinyint(1) default NULL,
  `note` varchar(255) default NULL,
  `gb_id` int(10) default NULL,
  `gb_locus` varchar(30) default NULL,
  `checked` tinyint(1) NOT NULL default '0',
  `race` enum('5''','3''','na') default 'na',
  `mol_type` varchar(40) default NULL,
  `seq_length` int(2) unsigned default '0',
  `seq_length_not_N` int(2) unsigned default '0',
  `max_frag_length_N_splitted` int(2) unsigned default '0',
  `x_masked_seq` mediumtext,
  `nrepeat` int(2) default '0',
  `xrepeat` int(2) default '0',
  `x_percent_masked` int(2) default '0',
  `n_percent_masked` int(2) default '0',
  `mapped` tinyint(1) NOT NULL default '0',
  `splk` tinyint(1) NOT NULL default '0',
  `old_tigem_trap_name` varchar(40) NOT NULL default '',
  PRIMARY KEY  (`trap_id`),
  UNIQUE KEY `trap` (`trap_name`,`project_id`),
  KEY `project_id` (`project_id`),
  KEY `clone_id` (`clone_id`),
  KEY `checked` (`checked`),
  KEY `mol` (`mol_type`),
  KEY `race` (`race`),
  KEY `mapped` (`mapped`),
  KEY `splk` (`splk`)
) ENGINE=MyISAM AUTO_INCREMENT=901289 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trap_maxicluster`
--

DROP TABLE IF EXISTS `trap_maxicluster`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trap_maxicluster` (
  `trap_maxicluster_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) default NULL,
  `maxicluster_id` int(10) default NULL,
  `trapmap_id` int(10) NOT NULL default '0',
  PRIMARY KEY  (`trap_maxicluster_id`),
  UNIQUE KEY `tum` (`trap_id`,`maxicluster_id`,`trapmap_id`),
  KEY `u` (`trap_id`),
  KEY `t` (`maxicluster_id`)
) ENGINE=MyISAM AUTO_INCREMENT=189465 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trap_trapcluster`
--

DROP TABLE IF EXISTS `trap_trapcluster`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trap_trapcluster` (
  `trap_trapcluster_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) default NULL,
  `trapcluster_id` int(10) default NULL,
  `trapmap_id` int(10) NOT NULL default '0',
  PRIMARY KEY  (`trap_trapcluster_id`),
  UNIQUE KEY `tum` (`trap_id`,`trapcluster_id`,`trapmap_id`),
  KEY `u` (`trap_id`),
  KEY `t` (`trapcluster_id`)
) ENGINE=MyISAM AUTO_INCREMENT=189465 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trap_unitrap`
--

DROP TABLE IF EXISTS `trap_unitrap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trap_unitrap` (
  `trap_unitrap_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) default NULL,
  `trapmap_id` int(10) NOT NULL default '0',
  `unitrap_id` int(10) default NULL,
  PRIMARY KEY  (`trap_unitrap_id`),
  UNIQUE KEY `tum` (`trap_id`,`unitrap_id`,`trapmap_id`),
  KEY `t` (`trap_id`),
  KEY `u` (`unitrap_id`)
) ENGINE=MyISAM AUTO_INCREMENT=202941 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock`
--

DROP TABLE IF EXISTS `trapblock`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock` (
  `trapblock_id` int(10) unsigned NOT NULL auto_increment,
  `trapmap_id` int(10) NOT NULL default '0',
  `trap_id` int(10) NOT NULL default '0',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` enum('-1','1') default NULL,
  `perc_ident` float NOT NULL default '0',
  PRIMARY KEY  (`trapblock_id`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `trapmap_id` (`trapmap_id`),
  KEY `trap_id` (`trap_id`)
) ENGINE=MyISAM AUTO_INCREMENT=777067 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_abinitio_annotation`
--

DROP TABLE IF EXISTS `trapblock_abinitio_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_abinitio_annotation` (
  `trapblock_abinitio_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) default NULL,
  `exon_name` varchar(40) default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `exon_coverage` float default NULL,
  `type` varchar(40) default NULL,
  PRIMARY KEY  (`trapblock_abinitio_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `display_name` (`display_name`),
  KEY `label` (`label`),
  KEY `exon_coverage` (`exon_coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_cdna_annotation`
--

DROP TABLE IF EXISTS `trapblock_cdna_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_cdna_annotation` (
  `trapblock_cdna_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `coverage` float default NULL,
  `comment` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`trapblock_cdna_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `display_name` (`display_name`),
  KEY `strand` (`strand`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM AUTO_INCREMENT=275169 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_ensest_annotation`
--

DROP TABLE IF EXISTS `trapblock_ensest_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_ensest_annotation` (
  `trapblock_ensest_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `strand` enum('-1','1') default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_ensest_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `strand` (`strand`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`),
  KEY `display_name` (`display_name`),
  KEY `trapblock_id` (`trapblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=17796333 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_ensestg_annotation`
--

DROP TABLE IF EXISTS `trapblock_ensestg_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_ensestg_annotation` (
  `trapblock_ensestg_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `rank` float NOT NULL default '0',
  `ensmusg` varchar(40) default NULL,
  `ensmust` varchar(40) default NULL,
  `ensmuse` varchar(40) default NULL,
  `strand` enum('-1','1') default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `exon_coverage` float default NULL,
  PRIMARY KEY  (`trapblock_ensestg_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `ensmusg` (`ensmusg`),
  KEY `ensmust` (`ensmust`),
  KEY `ensmuse` (`ensmuse`),
  KEY `label` (`label`),
  KEY `exon_coverage` (`exon_coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_ensg_annotation`
--

DROP TABLE IF EXISTS `trapblock_ensg_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_ensg_annotation` (
  `trapblock_ensg_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `rank` float NOT NULL default '0',
  `ensmusg` varchar(40) default NULL,
  `ensmust` varchar(40) default NULL,
  `ensmuse` varchar(40) default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `exon_coverage` float default NULL,
  `strand` enum('-1','1') default NULL,
  `refseq_id` varchar(50) NOT NULL default '',
  PRIMARY KEY  (`trapblock_ensg_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `strand` (`strand`),
  KEY `ensmusg` (`ensmusg`),
  KEY `ensmust` (`ensmust`),
  KEY `ensmuse` (`ensmuse`),
  KEY `label` (`label`),
  KEY `exon_coverage` (`exon_coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`),
  KEY `refseq_id` (`refseq_id`),
  KEY `trapblock_id` (`trapblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=777067 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_fantom_annotation`
--

DROP TABLE IF EXISTS `trapblock_fantom_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_fantom_annotation` (
  `trapblock_fantom_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_fantom_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `display_name` (`display_name`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_repeat_annotation`
--

DROP TABLE IF EXISTS `trapblock_repeat_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_repeat_annotation` (
  `trapblock_repeat_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_repeat_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_tcl_annotation`
--

DROP TABLE IF EXISTS `trapblock_tcl_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_tcl_annotation` (
  `trapblock_tcl_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_tcl_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `display_name` (`display_name`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_tclg_annotation`
--

DROP TABLE IF EXISTS `trapblock_tclg_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_tclg_annotation` (
  `trapblock_tclg_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_tclg_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `display_name` (`display_name`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapblock_unigene_annotation`
--

DROP TABLE IF EXISTS `trapblock_unigene_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapblock_unigene_annotation` (
  `trapblock_unigene_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) unsigned NOT NULL default '0',
  `trapblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `strand` enum('-1','1') default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapblock_unigene_annotation_id`),
  KEY `trap_id` (`trap_id`),
  KEY `trapblock_id` (`trapblock_id`),
  KEY `strand` (`strand`),
  KEY `display_name` (`display_name`),
  KEY `label` (`label`),
  KEY `coverage` (`coverage`),
  KEY `dS_s` (`dS_s`),
  KEY `dE_e` (`dE_e`)
) ENGINE=MyISAM AUTO_INCREMENT=69396 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapcluster`
--

DROP TABLE IF EXISTS `trapcluster`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapcluster` (
  `trapcluster_id` int(10) unsigned NOT NULL auto_increment,
  `accession` varchar(40) NOT NULL default '',
  `overlapping` varchar(40) default NULL,
  `maxicluster_id` int(10) default NULL,
  `sequence` mediumtext NOT NULL,
  `trapclusterblocks` int(10) unsigned NOT NULL default '0',
  `traps` int(10) unsigned NOT NULL default '0',
  `checked` tinyint(1) NOT NULL default '0',
  `hypertrapped` tinyint(1) NOT NULL default '0',
  `intronic` tinyint(1) NOT NULL default '0',
  `link_to_ensembl` mediumtext,
  `link_to_ucsc` mediumtext,
  `flanking` tinyint(1) NOT NULL default '0',
  `isolated` tinyint(1) NOT NULL default '0',
  `checked_flanking` tinyint(1) NOT NULL default '0',
  `antisense` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`trapcluster_id`),
  UNIQUE KEY `acc` (`accession`),
  KEY `overlapping` (`overlapping`),
  KEY `isolated` (`isolated`),
  KEY `antisense` (`antisense`),
  KEY `checked_flanking` (`checked_flanking`),
  KEY `checked` (`checked`),
  KEY `flanking` (`flanking`),
  KEY `hypertrapped` (`hypertrapped`),
  KEY `traps` (`traps`),
  KEY `trapclusterblocks` (`trapclusterblocks`),
  KEY `intronic` (`intronic`)
) ENGINE=MyISAM AUTO_INCREMENT=38315 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapcluster_novel_gene`
--

DROP TABLE IF EXISTS `trapcluster_novel_gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapcluster_novel_gene` (
  `trapcluster_novel_gene_id` int(10) unsigned NOT NULL auto_increment,
  `novel_gene_id` int(10) default NULL,
  `trapcluster_id` int(10) default NULL,
  PRIMARY KEY  (`trapcluster_novel_gene_id`),
  UNIQUE KEY `nt` (`novel_gene_id`,`trapcluster_id`),
  KEY `n` (`novel_gene_id`),
  KEY `t` (`trapcluster_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock`
--

DROP TABLE IF EXISTS `trapclusterblock`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock` (
  `trapclusterblock_id` int(10) unsigned NOT NULL auto_increment,
  `trapclustermap_id` int(10) NOT NULL default '0',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `sequence` mediumtext NOT NULL,
  `overlapping` varchar(40) default NULL,
  `intronic` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`trapclusterblock_id`),
  KEY `t` (`trapclustermap_id`),
  KEY `c` (`start`,`end`),
  KEY `overlapping` (`overlapping`),
  KEY `intronic` (`intronic`)
) ENGINE=MyISAM AUTO_INCREMENT=72393 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_abinitio_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_abinitio_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_abinitio_annotation` (
  `trapclusterblock_abinitio_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) default NULL,
  `exon_name` varchar(40) default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `strand` enum('-1','1') default NULL,
  `exon_coverage` float default NULL,
  `type` varchar(40) default NULL,
  PRIMARY KEY  (`trapclusterblock_abinitio_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=89368 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_cdna_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_cdna_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_cdna_annotation` (
  `trapclusterblock_cdna_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `coverage` float default NULL,
  `strand` enum('-1','1') default NULL,
  `comment` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`trapclusterblock_cdna_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=33924 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_ensest_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_ensest_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_ensest_annotation` (
  `trapclusterblock_ensest_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `strand` enum('-1','1') default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapclusterblock_ensest_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1679745 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_ensestg_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_ensestg_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_ensestg_annotation` (
  `trapclusterblock_ensestg_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `rank` float NOT NULL default '0',
  `ensmusg` varchar(40) default NULL,
  `ensmust` varchar(40) default NULL,
  `ensmuse` varchar(40) default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `strand` enum('-1','1') default NULL,
  `comment` varchar(255) default NULL,
  `exon_coverage` float default NULL,
  PRIMARY KEY  (`trapclusterblock_ensestg_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_ensg_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_ensg_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_ensg_annotation` (
  `trapclusterblock_ensg_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `rank` float NOT NULL default '0',
  `ensmusg` varchar(40) default NULL,
  `ensmust` varchar(40) default NULL,
  `ensmuse` varchar(40) default NULL,
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `exon_coverage` float default NULL,
  `strand` enum('-1','1') default NULL,
  `refseq_id` varchar(50) NOT NULL default '',
  PRIMARY KEY  (`trapclusterblock_ensg_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`),
  KEY `refseq_id` (`refseq_id`)
) ENGINE=MyISAM AUTO_INCREMENT=72393 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_fantom_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_fantom_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_fantom_annotation` (
  `trapclusterblock_fantom_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `strand` enum('-1','1') default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapclusterblock_fantom_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclusterblock_unigene_annotation`
--

DROP TABLE IF EXISTS `trapclusterblock_unigene_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclusterblock_unigene_annotation` (
  `trapclusterblock_unigene_annotation_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `trapclusterblock_id` int(10) unsigned NOT NULL default '0',
  `display_name` varchar(40) NOT NULL default '',
  `dS_s` int(10) default NULL,
  `dE_e` int(10) default NULL,
  `label` varchar(60) default NULL,
  `comment` varchar(255) default NULL,
  `strand` enum('-1','1') default NULL,
  `coverage` float default NULL,
  PRIMARY KEY  (`trapclusterblock_unigene_annotation_id`),
  KEY `trapcluster_id` (`trapcluster_id`),
  KEY `strand` (`strand`),
  KEY `trapclusterblock_id` (`trapclusterblock_id`)
) ENGINE=MyISAM AUTO_INCREMENT=9550 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapclustermap`
--

DROP TABLE IF EXISTS `trapclustermap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapclustermap` (
  `trapclustermap_id` int(10) unsigned NOT NULL auto_increment,
  `trapcluster_id` int(10) unsigned NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '0',
  `hit_db` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `strand` varchar(10) NOT NULL default '',
  PRIMARY KEY  (`trapclustermap_id`),
  UNIQUE KEY `clusthit` (`trapcluster_id`,`start`,`end`,`strand`),
  KEY `hit_id` (`hit_id`)
) ENGINE=MyISAM AUTO_INCREMENT=38315 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapens`
--

DROP TABLE IF EXISTS `trapens`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapens` (
  `trapens_id` int(10) unsigned NOT NULL auto_increment,
  `trapmap_id` int(10) unsigned NOT NULL default '0',
  `trap_id` int(10) unsigned NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '',
  `ensembl_id` varchar(255) default NULL,
  `insertion_ambiguous` tinyint(1) NOT NULL default '0',
  `number_trapblocks` int(10) unsigned default NULL,
  `number_annotated_trapblocks` int(10) unsigned default NULL,
  `source` varchar(40) NOT NULL default '',
  `trap_strand` enum('-','+') default NULL,
  `total_number_exons` int(10) default NULL,
  `putative_insertion_start` int(10) default NULL,
  `putative_insertion_end` int(10) default NULL,
  `insertion_case` int(1) default NULL,
  `trapped_exon_rank` int(10) default NULL,
  `exon_id` varchar(40) default NULL,
  `flanking_exon_id` varchar(40) default NULL,
  `name` varchar(40) NOT NULL default '',
  `description` varchar(255) NOT NULL default '',
  `ens_strand` enum('-','+') default NULL,
  `mol_type` varchar(40) NOT NULL default '',
  `annotation_ambiguous` tinyint(1) NOT NULL default '1',
  `label` varchar(40) NOT NULL default '',
  `refseq_id` varchar(40) NOT NULL default '',
  PRIMARY KEY  (`trapens_id`),
  UNIQUE KEY `un` (`trap_id`,`trapmap_id`),
  KEY `trapmap_id` (`trapmap_id`),
  KEY `ensembl_id` (`ensembl_id`),
  KEY `hit_id` (`hit_id`),
  KEY `exon_id` (`exon_id`),
  KEY `flanking_exon_id` (`flanking_exon_id`),
  KEY `putative_insertion_start` (`putative_insertion_start`),
  KEY `putative_insertion_end` (`putative_insertion_end`),
  KEY `trapped_exon_rank` (`trapped_exon_rank`),
  KEY `total_number_exons` (`total_number_exons`),
  KEY `number_trapblocks` (`number_trapblocks`),
  KEY `number_annotated_trapblocks` (`number_annotated_trapblocks`),
  KEY `label` (`label`),
  KEY `mol_type` (`mol_type`),
  KEY `annotation_ambiguous` (`annotation_ambiguous`)
) ENGINE=MyISAM AUTO_INCREMENT=310429 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapens_go`
--

DROP TABLE IF EXISTS `trapens_go`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapens_go` (
  `trapens_go_id` int(10) unsigned NOT NULL auto_increment,
  `accession` varchar(20) NOT NULL default '',
  `description` varchar(255) default NULL,
  `source` varchar(20) default NULL,
  `ensembl_id` varchar(40) NOT NULL default '',
  `go_class` varchar(255) default NULL,
  PRIMARY KEY  (`trapens_go_id`),
  KEY `a` (`accession`),
  KEY `d` (`description`),
  KEY `s` (`source`),
  KEY `c` (`go_class`),
  KEY `e` (`ensembl_id`)
) ENGINE=MyISAM AUTO_INCREMENT=73755 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trapmap`
--

DROP TABLE IF EXISTS `trapmap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `trapmap` (
  `trapmap_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) NOT NULL default '0',
  `hit_id` varchar(40) NOT NULL default '',
  `hit_db` varchar(40) NOT NULL default '',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `frac_aligned_query` float NOT NULL default '0',
  `num_hsps` int(10) NOT NULL default '0',
  `frac_identical` float NOT NULL default '0',
  `score` float NOT NULL default '0',
  `significance` double default NULL,
  `chosen` tinyint(1) NOT NULL default '0',
  `multiple` tinyint(1) default NULL,
  `chosen_filter` tinyint(1) NOT NULL default '0',
  `checked` tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (`trapmap_id`),
  UNIQUE KEY `hit` (`trap_id`,`hit_id`,`start`,`end`,`frac_aligned_query`,`num_hsps`,`frac_identical`,`score`,`chosen`),
  KEY `hit_id` (`hit_id`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `trap_id` (`trap_id`),
  KEY `c` (`chosen`),
  KEY `cf` (`chosen_filter`)
) ENGINE=MyISAM AUTO_INCREMENT=994641 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `traprpt`
--

DROP TABLE IF EXISTS `traprpt`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `traprpt` (
  `traprpt_id` int(10) unsigned NOT NULL auto_increment,
  `trap_id` int(10) default NULL,
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `score` float default NULL,
  `tag` varchar(25) default NULL,
  PRIMARY KEY  (`traprpt_id`),
  KEY `t` (`trap_id`),
  KEY `a` (`tag`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `unitrap`
--

DROP TABLE IF EXISTS `unitrap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `unitrap` (
  `unitrap_id` int(10) unsigned NOT NULL auto_increment,
  `accession` varchar(40) NOT NULL default '',
  `hit_db` varchar(40) NOT NULL default '',
  `chr` varchar(40) NOT NULL default '0',
  `start` int(10) NOT NULL default '0',
  `end` int(10) NOT NULL default '0',
  `ensembl_id` varchar(40) NOT NULL default '',
  `gene_name` varchar(40) NOT NULL default '',
  `gene_description` varchar(255) default NULL,
  `gene_strand` enum('+','-') default NULL,
  `gene_type` varchar(40) default NULL,
  `gene_start` int(10) NOT NULL default '0',
  `gene_end` int(10) NOT NULL default '0',
  `refseq_id` varchar(40) NOT NULL default '',
  `ambiguous` tinyint(1) NOT NULL default '0',
  `exon_id` varchar(255) NOT NULL default '',
  `flanking_exon_id` varchar(255) NOT NULL default '',
  PRIMARY KEY  (`unitrap_id`),
  UNIQUE KEY `uni` (`accession`),
  UNIQUE KEY `unihit` (`unitrap_id`,`chr`,`start`,`end`,`gene_strand`),
  KEY `chr` (`chr`),
  KEY `start` (`start`),
  KEY `end` (`end`),
  KEY `ensembl_id` (`ensembl_id`),
  KEY `gene_name` (`gene_name`),
  KEY `gene_description` (`gene_description`),
  KEY `gene_strand` (`gene_strand`),
  KEY `gene_type` (`gene_type`)
) ENGINE=MyISAM AUTO_INCREMENT=38729 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `xref`
--

DROP TABLE IF EXISTS `xref`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL auto_increment,
  `ensembl_id` varchar(40) NOT NULL default '',
  `dbname` varchar(40) default NULL,
  `accession` varchar(40) default NULL,
  PRIMARY KEY  (`xref_id`),
  KEY `ensembl_id` (`ensembl_id`),
  KEY `dbname` (`dbname`),
  KEY `accession` (`accession`)
) ENGINE=MyISAM AUTO_INCREMENT=271385 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2009-09-09 16:37:13
