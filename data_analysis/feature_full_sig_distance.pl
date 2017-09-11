#! /bash/usr/env perl

# Will create a combined matrice and combined distance to signal graph. Assumes the data is from a human

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use Path::Class;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
use Data::Dumper;

#My variables
my $help = 0;
my $man = 0;
my $base_name;
my $out;
my $loop_file;
my $bin_length;
my $rbin;

# Create Chrom Array
my @chrom;
$chrom[23] = "Y";
$chrom[22] = "X";
for ( my $i = 0; $i < scalar(@chrom); $i++ ) {
	$chrom[$i] = $i + 1;
}


#Read in the variables from the command line
GetOptions( 'man'		=>	\$man,
			'help|h'	=>	\$help,
			'base_file_name|bf=s'	=>	\$base_name,
			'out_dir|o=s'	=>	\$out,
			'loop_file|l=s'	=>	\$loop_file,
			'bin_length|bl=i'	=> \$bin_length,
			'rcode_bin|r=s'	=>	\$rbin,
		   ) || die("There was an error in the command line arguements");

# Pod Usage for the manual and the help pages
if ($help) { pod2usage(0) }
if ($man) { pod2usage(-verbose => 3) }

# Setup Logging Environment
my $logger = get_logger();

## MAIN ##

my $count_files_href = get_count_file_names_for_each_chrom(\@chrom, $base_name);
create_matrice_files($count_files_href, $out);
merge_matrice_and_create_dist_vs_sig_graph(\@chrom,$out);

## Subroutines ##
sub get_count_file_names_for_each_chrom {
	my ( $chrom_aref, $bn ) = @_;
	$logger->info("Geting the count files for each chromosome");
	my @dir_full = split( /\//, $bn);
	my $just_name = pop(@dir_full);
	my $dir = join("/", @dir_full);
	
	my @full_dir_files;
	opendir(my $DIR, $dir) || die("$dir could not be found\n");
	while (readdir($DIR)) {
		next if ( $_ =~ /^\./ );
		push @full_dir_files, $_;
	}
	closedir($DIR);
	print(Dumper(@full_dir_files))
	#fill a hash with file names
	my %file;
	foreach my $chromo ( @$chrom_aref ) {
		foreach my $fname ( @full_dir_files ) {
			next if ( $fname =~ /^\./ );
			if ( $chromo =~ qr/$bin_length\Kb/ && $chromo =~ qr/chr$chromo/ ) {
				$file{$chromo} = "$dir/$fname";
			}
		}
	}
	print Dumper(%file);
	#return hash with count file names
	return \%file;
}

sub create_matrice_files {
	my ( $count_href, $o ) = @_;
	my @jobs;
	foreach my $key ( keys %$count_href ) {
		my $file = $count_href->{$key};
		my $cmd = "Rscript $rbin/get_feature_specific_data.R --loop_file $loop_file --count_file $file --chromosome $key --out_dir $o";
		push @jobs, $cmd;
	}
	submit_and_stall(\@jobs, "matrice", $o, "20");
	
	return;
}

sub merge_matrice_and_create_dist_vs_sig_graph {
	my ( $chrom_aref, $o ) = @_;
	$logger->info("Combining all the signal to distance matrices and creating a signal graph");
	
	my $cmd = "bsub -M 20 -J final -o $o/final.out Rscript $rbin/graph_tot_dist_sig.R $o"; #Command that will be ran
	my @array = ("background_", "TADs_", "loop_flare");
	foreach my $chromo ( @$chrom_aref ) {
		foreach my $name ( @array ) {
			my $full_name = "$out/" . $name . $chromo . ".txt";
			$cmd .= " $full_name";
		}
	}
	
	#cmd should be complete so submit the command here
	system($cmd);
	return;
}

sub submit_and_stall {
	my ( $job_aref, $job_name, $out_directory, $mem) = @_;
	$logger->info("Submitting and Stalling until all $job_name jobs finish\n");
	print Dumper($job_aref);
	foreach my $cmd (@$job_aref) {
		my $bsub_cmd = "bsub -M $mem -J $job_name -o $out_directory/$job_name.out $cmd";
		system($bsub_cmd);
	}
	
	#Create Active Jobs File
	`echo start > $out_directory/ACTIVE_JOBS`;
	
	my $minutes = 0;
	while ( -s "$out_directory/ACTIVE_JOBS" ) {
		sleep 60;
		my $cmd = "bjobs -J " . $job_name . " > $out_directory/ACTIVE_JOBS";
		`$cmd`;
		$minutes++;
	}
	
	`rm $out_directory/ACTIVE_JOBS`;
	$logger->info("All $job_name jobs have finished in ~$minutes minutes");
	
	return 1;
}

__END__
=head1	TITLE

=head1	VERSION

=head1	INCLUDED MODULES

Getopt::Long;
Pod::Usage;
Carp;
Path::Class;
Log::Log4Perl qw(easy);
Log::Log4perl::CommandLine qw(:all);
Data::Dumper;

=head1	SYNOPSIS

=head1	PARAMETERS

=head1	CONFIGURATION AND ENVIRONMENT

=head1	DEPENDENCIES

=head1	INCOMPATIBILITIES

	None reported.

=head1	BUGS AND LIMITATIONS

	Please report any bugs or feature requests
	
=head1	AUTHOR

Nicholas Colaianni
Contact via email C<< <ncolaian@live.unc.edu> >>

=head1	LICENCE AND COPYRIGHT

=head1	DISCLAIMER OF WARRANTY