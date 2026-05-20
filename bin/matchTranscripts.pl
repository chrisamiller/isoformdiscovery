#!/usr/bin/perl
#
# matchTranscripts.pl — Filter a TSV by matching transcript IDs in the first column.
#
# Usage:
#   matchTranscripts.pl <id_list.txt> <data.tsv> > filtered.tsv
#
# Arguments:
#   [0] id_list.txt  — one transcript ID per line (first column used if tab-delimited)
#   [1] data.tsv     — TSV file to filter; first column may be a comma-separated list
#                      of transcript IDs (ESPRESSO abundance/classification format)
#
# Output: Lines from data.tsv where at least one transcript ID in column 0
#         matches an ID in the id_list. Output includes the full original line.
#
# Used as a general-purpose matching utility when filtering ESPRESSO output tables
# to a specific transcript set (e.g., after SQANTI3 filtering).

use warnings;
use strict;
use IO::File;

my %transcripts;

# Load transcript IDs from first column of id_list
my $inFh = IO::File->new( $ARGV[0] ) || die "can't open $ARGV[0]: $!\n";
while ( my $line = $inFh->getline ) {
    chomp $line;
    my @F = split "\t", $line;
    $transcripts{ $F[0] } = 1;
}
close $inFh;

# Filter data file: keep rows where any transcript ID in col 0 matches
$inFh = IO::File->new( $ARGV[1] ) || die "can't open $ARGV[1]: $!\n";
while ( my $line = $inFh->getline ) {
    chomp $line;
    my @F = split "\t", $line;
    my @ids = split ",", $F[0];
    for my $id (@ids) {
        if ( $transcripts{$id} ) {
            print "$line\n";
            last;
        }
    }
}
close $inFh;
