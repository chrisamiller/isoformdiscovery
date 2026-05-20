#!/usr/bin/perl
#
# remove_nonfsm_data.pl — Filter ESPRESSO intermediate files to retain only
# Full Splice Match (FSM) reads, then filter the per-sample BAMs accordingly.
#
# Usage:
#   remove_nonfsm_data.pl <compatible_isoform.tsv> <espresso_chr_dir> <chr.tsv.updated>
#
# Arguments:
#   [0] compatible_isoform.tsv  — per-read transcript assignments from ESPRESSO_Q
#                                  (filtered/round-2 run); col 0 = read_id, col 2 = class
#   [1] espresso_chr_dir        — ESPRESSO working directory for this chromosome
#                                  (contains per-sample subdirs 0/, 1/, ..., sj.list, etc.)
#   [2] chr.tsv.updated         — ESPRESSO input TSV (bam_path <TAB> sample_id <TAB> index)
#
# Output: creates <espresso_chr_dir>/fsm/ with FSM-filtered copies of all intermediate
#         files and per-sample BAMs. Also writes a new chr.tsv.updated pointing to the
#         FSM-filtered BAMs.
#
# Bug fix applied (vs original scripts/remove_nonfsm_data.pl):
#   In sj.list processing, the column-8 block previously iterated over $F[7] instead of
#   $F[8]. This caused FSM read IDs from col 7 to be used when filtering col 8, meaning
#   col 8 read IDs were never filtered. Fixed below.

use warnings;
use strict;
use IO::File;
use File::Basename;
use File::Spec;

# ── Build FSM read set from compatible_isoform.tsv ─────────────────────────
print STDERR "reading fsm read ids\n";
my %reads;
my $inFh = IO::File->new( $ARGV[0] ) || die "can't open $ARGV[0]: $!\n";
while ( my $line = $inFh->getline ) {
    chomp $line;
    my @F = split "\t", $line;
    if ( $F[2] eq "FSM" ) {
        $reads{ $F[0] } = 1;
    }
}
close $inFh;
print STDERR "done reading fsm read ids\n";

my $dir = $ARGV[1];
unless ( -e "$dir/fsm" ) {
    mkdir "$dir/fsm" or die "cannot mkdir $dir/fsm: $!";
}

# ── Filter *_read_final.txt files ──────────────────────────────────────────
print STDERR "filtering read_final files in $dir\n";
my @read_final_files = glob "$dir/*/*_read_final.txt";
foreach my $file (@read_final_files) {
    print STDERR "  $file\n";
    my $subdir = basename( dirname($file) );
    mkdir "$dir/fsm/$subdir" unless -e "$dir/fsm/$subdir";
    open( my $outFh, '>', "$dir/fsm/$subdir/" . basename($file) )
        or die "can't open output: $!";
    my $inFh2 = IO::File->new($file) || die "can't open $file: $!";
    while ( my $line = $inFh2->getline ) {
        chomp $line;
        my @F = split "\t", $line;
        print $outFh "$line\n" if $reads{ $F[0] };
    }
    close $inFh2;
    close $outFh;
}

# ── Filter sam.list3 files ──────────────────────────────────────────────────
print STDERR "filtering sam.list3 files in $dir\n";
my @sam_files = glob "$dir/*/sam.list3";
foreach my $file (@sam_files) {
    print STDERR "  $file\n";
    my $subdir = basename( dirname($file) );
    mkdir "$dir/fsm/$subdir" unless -e "$dir/fsm/$subdir";
    open( my $outFh, '>', "$dir/fsm/$subdir/" . basename($file) )
        or die "can't open output: $!";
    my $inFh2 = IO::File->new($file) || die "can't open $file: $!";
    while ( my $line = $inFh2->getline ) {
        chomp $line;
        my @F = split "\t", $line;
        print $outFh "$line\n" if $reads{ $F[2] };
    }
    close $inFh2;
    close $outFh;
}

# ── Filter sj.list files ────────────────────────────────────────────────────
# sj.list columns (0-indexed):
#   5 = read count for col 7 junction strand
#   6 = read count for col 8 junction strand
#   7 = comma-separated read IDs for strand 1  (trailing comma)
#   8 = comma-separated read IDs for strand 2  (trailing comma)
print STDERR "filtering sj.list files in $dir\n";
my @sj_files = glob "$dir/*/sj.list";
foreach my $file (@sj_files) {
    print STDERR "  $file\n";
    my $subdir = basename( dirname($file) );
    mkdir "$dir/fsm/$subdir" unless -e "$dir/fsm/$subdir";
    open( my $outFh, '>', "$dir/fsm/$subdir/" . basename($file) )
        or die "can't open output: $!";
    my $inFh2 = IO::File->new($file) || die "can't open $file: $!";
    while ( my $line = $inFh2->getline ) {
        chomp $line;
        my @F = split "\t", $line;

        # Filter column 7 (strand-1 read IDs)
        unless ( $F[7] eq "NA" ) {
            my @names = split ",", $F[7];
            my @keep  = grep { $reads{$_} } @names;
            $F[5] = scalar @keep;
            $F[7] = @keep ? join( ",", @keep ) . "," : "NA";
        }

        # Filter column 8 (strand-2 read IDs)
        # BUG FIX: original code split $F[7] here instead of $F[8]
        unless ( $F[8] eq "NA" ) {
            my @names = split ",", $F[8];    # was $F[7] — now correctly $F[8]
            my @keep  = grep { $reads{$_} } @names;
            $F[6] = scalar @keep;
            $F[8] = @keep ? join( ",", @keep ) . "," : "NA";
        }

        print $outFh join( "\t", @F ) . "\n";
    }
    close $inFh2;
    close $outFh;
}

# ── Filter per-sample BAMs ──────────────────────────────────────────────────
mkdir "$dir/fsm/bams" unless -e "$dir/fsm/bams";

my ( @bams, @samps, @nums );
my $inFh2 = IO::File->new( $ARGV[2] ) || die "can't open $ARGV[2]: $!";
while ( my $line = $inFh2->getline ) {
    chomp $line;
    my @F = split "\t", $line;
    push @bams,  $F[0];
    push @samps, $F[1];
    push @nums,  $F[2];
}
close $inFh2;

open( my $configFh, '>', "$dir/fsm/" . basename( $ARGV[2] ) )
    or die "can't open output config: $!";

for ( my $i = 0; $i < @bams; $i++ ) {
    my $samname = "$dir/fsm/bams/$samps[$i].sam";
    my $bamname = "$dir/fsm/bams/$samps[$i].bam";
    print STDERR "writing FSM BAM for $samps[$i]\n";

    open( my $outFh, '>', $samname ) or die "can't open $samname: $!";
    open( my $bamFh, "samtools view -h $bams[$i] |" )
        or die "can't run samtools: $!";
    while ( my $line = <$bamFh> ) {
        chomp $line;
        my @F = split "\t", $line;
        print $outFh "$line\n" if $reads{ $F[0] } || $line =~ /^\@/;
    }
    close $bamFh;
    close $outFh;

    system("samtools view -Sb -o $bamname $samname") == 0
        or die "samtools view failed: $!";
    unlink $samname;

    print $configFh join( "\t",
        File::Spec->rel2abs($bamname),
        $samps[$i],
        $nums[$i]
    ) . "\n";
}
close $configFh;
