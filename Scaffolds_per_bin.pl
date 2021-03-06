#!/usr/bin/env perl
###############################################################################
#
#    Scaffolds_per_bin.pl
#
#    Reads in scaffold data and summarizes into table
#
#    Copyright (C) Inka Vanwonterghem
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
use Data::Dumper;

#CPAN modules
use Bio::SeqIO;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################

#open the sspace output file:
#   scaffold1    15      521     contig15       contig521    BIN_1   BIN_5      5
#   scaffold50   256     21445   contig256      contig21445  BIN_240 BIN_2144   10
my $scaffolds = openRead($global_options->{'sspace'});

# write to a file
my $perBIN = openWrite($global_options->{'output_perBIN'});

#global variables
my %linestore;
my %bin2scaff;
my $scaff;
my $bin1;
my $bin2;

#read the sspace output file and create:
#a hash of arrays containing all the scaffold lines
#a hash of hashes containing all the scaffolds associated with one bin

foreach my $line (<$scaffolds>) {
    chomp $line;
    my @array = split(/\t/,$line);
    #create the hash of arrays of scaffold ids
    $scaff = $array[0];
    if (exists $linestore{$scaff}) {
        push (@{$linestore{$scaff}},$line);
    }
    else {
        my @tmp_array = ();
        push(@tmp_array, $line);
        $linestore{$scaff} = \@tmp_array;
    }
    #create the hash of hashes of bin ids
    $bin1 = $array[5];
    $bin2 = $array[6];
    if (exists $bin2scaff{$bin1}) {
        $bin2scaff{$bin1}->{$scaff} = 0;
    }
    else {
        my %tmp_hash;
        $tmp_hash{$scaff} = 0;
        $bin2scaff{$bin1} = \%tmp_hash;
    }
    if (exists $bin2scaff{$bin2}) {
        $bin2scaff{$bin2}->{$scaff} = 0;
    }
    else {
        my %tmp_hash;
        $tmp_hash{$scaff} = 0;
        $bin2scaff{$bin2} = \%tmp_hash;
    }
}

foreach my $bin (keys %bin2scaff) {
    print $perBIN "\n$bin\n";
    my @scaffolds = keys(%{$bin2scaff{$bin}});
    foreach my $scaffold (@scaffolds) {
        my $lines = $linestore{$scaffold};
        foreach my $line (@{$lines}){
            #print OUT Data::Dumper->Dump($line)."\n";
            print $perBIN "$line\n";
        }
    }
}

close ($scaffolds);
close ($perBIN);


######################################################################
# CUSTOM SUBS
######################################################################

sub plus2 {
  my ($in) = @_;
  return $in + 2;

}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    # "long|single<type>"
    #
    # :i number (integer)
    # :f number (decimal)
    # :s string
    # +  flag
    #
    my @standard_options = ( "help|h+", "sspace|s:s", "output_perBIN|o:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    #if(!exists $options{''} ) { printParamError (""); }
    if(!exists $options{'sspace'} ) { printParamError ("You need to supply the sspace summary file"); }
    if(!exists $options{'output_perBIN'} ) { printParamError ("You need to supply an output file"); }
    
    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub createDefault
{
    #-----
    # Set default values for parameters
    #
    my ($option_name, $default_value) = @_;
    if(not exists $global_options->{$option_name})
    {
      $global_options->{$option_name} = $default_value;
    }
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          },
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;

    isCommandInPath($cmd, $failure_type);

    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};

    my $cmd_str = $cmd . " " . $param_str;

    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to
    # checkAndRunCommand
    #
    my $ref = shift;

    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
 $0
 Copyright (C) Inka Vanwonterghem

 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

    Scaffolds_per_bin.pl

=head1 COPYRIGHT

   copyright (C) Inka Vanwonterghem

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    This script takes the sspace output (after parsing with script Jason -
    Scaffolding_binned_contigs_Jason) and creates two new files: one summary
    of all the inter-intra-unbinned links for each bin (total and percentage),
    one summary of all the links for each bin

=head1 SYNOPSIS

    Scaffolds_per_bin.pl -s <SSPACE.txt> <output_perBIN.txt> [-help|h]

    -sspace -s                  sspace CSV from Jason's scripts
    -output_perBIN -o          output file to write the summary to
    [-help -h]                  Displays basic usage information

=cut

