#!/usr/bin/perl
#
# filter_gtf.pl — Extract GTF records whose transcript_id matches an ID list.
#
# Usage:
#   filter_gtf.pl <id_list.txt> <annotation.gtf> > filtered.gtf
#
# Arguments:
#   [0] id_list.txt  — one transcript ID per line (first column used if tab-delimited)
#   [1] annotation.gtf — GTF file to filter
#
# Output: GTF lines where transcript_id attribute matches an ID in the list.
#         Comment lines (beginning with #) are passed through unchanged.
#
# Used in RESTORE_ENSEMBL to re-add Ensembl transcripts that the ML filter
# classified as Artifact. Only transcript_id matching is performed; gene lines
# are included if any transcript in that gene is kept.

use warnings;
use strict;
use IO::File;

my %transcripts;

# Read transcript IDs from first column of id_list
my $inFh = IO::File->new( $ARGV[0] ) || die "can't open $ARGV[0]: $!\n";
while ( my $line = $inFh->getline ) {
    chomp $line;
    my @F = split "\t", $line;
    $transcripts{ $F[0] } = 1;
}
close $inFh;

# Stream GTF; output lines matching the transcript list
my $removed = 0;
$inFh = IO::File->new( $ARGV[1] ) || die "can't open $ARGV[1]: $!\n";
while ( my $line = $inFh->getline ) {
    chomp $line;
    if ( $line =~ /^#/ ) {
        print "$line\n";
        next;
    }
    if ( $line =~ /transcript_id "([^"]+)"/ ) {
        if ( defined $transcripts{$1} ) {
            print "$line\n";
        } else {
            my @F = split "\t", $line;
            $removed++ if $F[2] eq "transcript";
        }
    }
}
close $inFh;
print STDERR "removed $removed transcripts\n";
