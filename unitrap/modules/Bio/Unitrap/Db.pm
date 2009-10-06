package Bio::Unitrap::Db;

use strict;
use DBI;
use Carp;

sub db_connection {
    my ($self, $host, $db, $user, $password, $debug) = @_;
	
    my $dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(	-host=> $host,
							-user=> $user,
							-pass=> $password,
					       		-dbname=> $db	);
					     							
    $debug && print "SQL CONNECTION:\nHOST=$host, DB=$db, USER=$user, PASS=$password\n\n";
    
    return $dbh;
}

sub prepare_stmt {
    my ($self,$dbh,$par,$debug) = @_;
    
    my %params = %{$par};
    #construct statment and quote
    my $stmt;
    foreach my $f (keys %params){
	$stmt .= " , " if $stmt;
	$stmt .= $f . " = \"$params{$f}\"";
    }
    
    $debug && print "DB: prepare_stmt=>STMT $stmt\n";
    
    return $stmt;
}


sub insert_set {
    my ($self, $dbh, $stmt, $table_name, $debug)=@_;
    
    my $s = "INSERT INTO $table_name SET  $stmt";
    $debug && print STDOUT "SQL CODE: $s\n";
    my $sth = $dbh->prepare($s);
    $sth->execute();# || die "insert failed : $DBI::errstr";
    my $dbi = $sth->{'mysql_insertid'};
    $debug && print STDOUT "DB:insert_set => Table $table_name last inserted id is $dbi\n";
    
    return $dbi;
}


sub check_return {
    my ($self,$value, $table_name, $table_column, $sth, $debug) = @_;
    
    my $r = $sth->rows;
    
    if ($r == 0) {
	$debug && print STDERR "DB:check_return => CANNOT FIND THIS $value IN THE $table_name TABLE AND $table_column COLUMN\n";
	return undef;
    } else {
	return 1;
    }
}


sub get_table_id {
    my ($self, $dbh, $name, $table_name, $column_name, $debug) = @_;
    
    my $stmt = "select ". $table_name ."_id from $table_name where $column_name = ?";
    my $sth = $dbh->prepare($stmt);
    $sth->execute($name) || die "select failed : $DBI::errstr";
    my @array = $sth->fetchrow_array();
    my $check = $self->check_return ($name, $table_name ,$column_name, $sth, $debug);

    if ($check) {
	return $array[0];
    } else {
	return undef;
    }
}

sub update_set {
    my ($self, $dbh, $stmt, $where, $table_name, $debug) = @_;
    
    my $s = "UPDATE $table_name SET  $stmt WHERE $where";
    $debug && print STODUT "SQL CODE: $s\n";
    my $sth = $dbh->prepare($s);
    $sth->execute() || die "update failed : $DBI::errstr";
    my $dbi = $sth->rows;
    
    return $dbi;
}

sub delete_from_table {
    my ($self, $dbh, $table_name, $where) = @_;
    
    my $s = "DELETE FROM $table_name WHERE $where";
    my $sth = $dbh->prepare($s);
    $sth->execute() || die "delete failed : $DBI::errstr";
    my $dbi = $sth->rows;
    
    return $dbi;
}

sub select_from_table {
    my ($self, $dbh, $fields, $table_name, $condition, $debug) = @_;
    my $s = "SELECT $fields FROM $table_name ";
    
    if ($condition) {
    	$s .= "$condition";
    }
   
    $debug && print STDOUT "SQL CODE: $s\n";
    
    my $sth = $dbh->prepare($s);
    $sth->execute() || die "delete failed : $DBI::errstr";	

    my @array = $sth->fetchall_arrayref();
    if (@array) {return \@array;}
}

sub create_db () {
    my ($self, $path, $user, $pass, $db, $host, $debug) = @_;

    $debug && print STDOUT $path."mysqladmin -u $user -p$pass -h $host CREATE $db\n";
    system $path."mysqladmin -u $user -p$pass -h $host CREATE $db";	
}

sub exec_dump () {
    my ($self, $path, $user, $pass, $db, $host, $file, $no_data, $debug) = @_;
    my $attr;

    if ($no_data) {
        $attr = "--no_data";
    }
	
    $debug && print STDOUT $path."mysqldump -u $user -p$pass -h $host $attr $db > $file\n";
    system $path."mysqldump -u $user -p$pass -h $host $attr $db > $file";	
}

sub exec_command_sql () {
    my ($self, $path, $user, $pass, $db, $host, $file, $sql, $debug) = @_;

    $debug && print STDOUT $path."mysql -u $user -h $host -p$pass $db -e \"$sql\" > $file\n";
    system $path."mysql -u $user -h $host -p$pass $db -e \"$sql\" > $file";
}

sub exec_import () {
    my ($self, $path, $user, $pass, $db, $host, $file, $debug) = @_;

    $debug && print STDOUT $path."mysql -u $user -p$pass -h $host $db < $file\n";
    system $path."mysql -u $user -p$pass -h $host $db < $file";
}

sub import_into_table () {
    my ($self, $path, $user, $pass, $db, $host, $fields, $file, $debug) = @_;

    $debug && print "$path"."mysqlimport -u $user -p$pass -h $host $db -c $fields $file\n";
    system $path."mysqlimport -u $user -p$pass -h $host $db -c $fields $file";
}

1;
