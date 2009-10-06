package Bio::Unitrap::File;

use strict;
use Data::Dumper;

sub readFileData () {
    my ($self, $file, $debug) = @_;	
	my $content;

	if ($debug) {print "Reading the file $file\n";}
	
	open (FILE, $file);
	while (defined (my $line=<FILE>)) {$content .= $line;}
	close (FILE);
	return $content;
}

sub writeDataToFile () {
    my ($self, $file, $data, $debug) = @_;

	if ($debug) {print "Writing in the file $file\n";}
	
	open (FILE, ">$file");
	print FILE "$data";
	close (FILE);	
}

sub appendDataToFile () {
    my ($self, $file, $data, $debug) = @_;

	if ($debug) {print "Writing in the file $file\n";}
	
	open (FILE, ">>$file");
	print FILE "$data";
	close (FILE);	
}

sub fromFastaFileToHash () {
	my ($self, $file, $debug) = @_;
		
	my $cnt = $self->readFileData ($file, $debug);	
	my %hash = $self->fromFastaFormatToHash ($cnt, $debug);

	return %hash;
}

# Convert a multifasta text in a hash of sequences. the Key is the sequence name and the value is the sequence
sub fromFastaFormatToHash () {
    my ($self, $text, $debug) = @_;	
	
	if ($debug) {print "Moving the content of the fasta file into a hash\n";}
	
	my @array = split (/>/, $text);
	my ($sequence, %hash);
	
	foreach my $seq (@array) {
		if ($seq) {
			my @lines = split ("\n", $seq);
			my $description = shift @lines;
			
			my @words = split (" ", $description);
			my $name = $words[0];
			my $desc = $words[1];
			
			foreach my $line (@lines) {$sequence = $sequence.$line;}
			
			$name =~ s/\,//g;
			
			$hash {$name}{'sequence'} = $sequence;
			$hash {$name}{'desc'} = $desc;
			
			undef $name;
			undef $sequence;
			undef $desc;
		}
	}
	return %hash; 
}

1;