BEGIN {
    package main;
    use lib "/var/perl/bioperl-run";
    use lib "/var/perl/unitrap/modules";
    #use lib "/var/perl/bioperl-1.5.2";
    use lib "/var/perl/roma/bioperl";
    use lib "/var/perl/ensembl-52/modules";
    use lib "/var/perl/bioperl-live/core/";
    use lib "/var/perl/biomart-perl/lib";
    use lib "/var/perl/biomart-perl";

    $ENV{'BLASTFILTER'} = '/usr/local/bio/wublast/filter';
    $ENV{'BLASTMAT'} = '/usr/local/bio/wublast/matrix/aa';
    $ENV{'BLASTDATADIR'} = 'data/blat/db/mus_musculus_52';
    
    %conf = (
             'enshost' => "dbmaster",
             'ensuser' => "mysql_dev",
             'enspass' => "dEvEl0pEr",
             'ensdbname' => "mus_musculus_core_52_37e",
             'ensestdbname' => "mus_musculus_otherfeatures_52_37e",
             'traphost' => "hal9000",
             'trapuser' => "mysql_dev",
             'trapdbname' => "UNITRAP3",
             'trappass' => "dEvEl0pEr",
             'fantom3host' => "hal9000",
	     'fantom3user' => "mysql_dev",
	     'fantom3dbname' => "fantom3_ann_43",
	     'fantom3pass' => "dEvEl0pEr",
             'ens_hs_dbhost' => "dbmaster",
	     'ens_hs_dbuser' => "mysql_dev",
	     'ens_hs_dbname' => "homo_sapiens_core_52_36n",
	     'ens_hs_pass' => "dEvEl0pEr",
             'marthost' => "dbmaster",
	     'martuser' => "mysql_dev",
	     'martdbname' => "ensembl_mart_52",
	     'martpass' => "dEvEl0pEr",
	     'hit_db' => "NCBIM37", 
	     'ensembl_base_url' => "http://may2009.archive.ensembl.org/Mus_musculus/",
	     'ucsc_db_name' => "mm9",
	     'blastdb' => '/data/blat/db/mus_musculus_48/mus_musculus_e48',
	    # 'blastdb' => '/data/blat/db/mus_musculus_e43/repeat_masked/mus_musculus_e43_rm_gr',
            #'blastprotdb' => '/data/blastdbs/guglielmo/uniref90.fasta',
             'blastnexec' =>  '/usr/local/bio/wublast/blastn',
            #'dashost' => "elia.tigem.it",
            #'dasuser' => "mysql_dev",
            #'dasdbname' => "current_trap_das",
            #'daspass' => 'dEvEl0pEr',
             'mysql_path' => "/usr/bin/",
             'tmp_dir' => "/home/roma/src/scripts/data/tmp/",
             'repeatMasker_dir' => "",
	     'debug' => 1,
	     'debugSQL' => 1
        );
}
1;
