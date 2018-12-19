#!/usr/bin/perl

use strict;
use warnings;

my $OUT_PREFIX = $ARGV[1];

# prepare BED files for directed trans interactions
my $filename = $OUT_PREFIX . "_directed_trans.bed";
open(my $directed_trans, '>', $filename) or die "Could not open file '$filename' $!";

$filename = $OUT_PREFIX . "_undirected_trans.bed";
open(my $undirected_trans, '>', $filename) or die "Could not open file '$filename' $!";

# prepare BED files for directed trans interactions
$filename = $OUT_PREFIX . "_directed_strong_long.bed";
open(my $directed_strong_long, '>', $filename) or die "Could not open file '$filename' $!";

# prepare BED files for directed trans interactions
$filename = $OUT_PREFIX . "_undirected_strong_long.bed";
open(my $undirected_strong_long, '>', $filename) or die "Could not open file '$filename' $!";

# command for extracting FASTA from BED
# bedtools getfasta -name -fi ../VPV_data/mm9/mm9.fa -bed xxx_undirected_trans.bed > xxx_undirected_trans.fasta

# command for replacing all small letters in FASTA wit N

# sed -e '/^>/! s/[[:lower:]]/N/g' in.fasta > out.fasta



my $LONG_RANGE_CUTOFF = 10000;

########## sub definitions first

sub getInteractionDist {

    my ($frag1_chr, $frag1_sta, $frag1_end, $frag2_chr, $frag2_sta, $frag2_end) = @_;

    print "Distance between fragments on different chromosomes is not defined. I'm going to exit now!" and exit if $frag1_chr ne $frag2_chr;

    my $return_dist;

    if($frag1_end < $frag2_sta) {
        $return_dist = $frag2_sta - $frag1_end;
    } else {
        $return_dist = $frag1_sta - $frag2_end;
    }

    return $return_dist;

}


open(my $fh, "<", $ARGV[0])
	or die "Can't open < input.txt: $!";

my $N_CIS_SIMPLE_READ_PAIRS = 0;
my $N_CIS_TWISTED_READ_PAIRS = 0;
my $N_TRANS_SIMPLE_READ_PAIRS = 0;
my $N_TRANS_TWISTED_READ_PAIRS = 0;

my $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS = 0;
my $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS = 0;
my $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS = 0;
my $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS = 0;

my $N_CIS_DIRECTED_INTERACTION = 0;
my $N_CIS_UNDIRECTED_INTERACTION = 0;
my $N_TRANS_DIRECTED_INTERACTION = 0;
my $N_TRANS_UNDIRECTED_INTERACTION = 0;

my $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION = 0;
my $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION = 0;
my $N_CIS_LONG_RANGE_DIRECTED_INTERACTION = 0;
my $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION = 0;

my $N_AT_LEASTS_TWO_READ_PAIR_INTERACTION = 0;
my $SUM_AT_LEASTS_TWO_READ_PAIR_INTERACTION = 0;
my $AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION;

my $N_WEAK_DIRECTED_INTERACTION = 0;
my $N_WEAK_UNDIRECTED_INTERACTION = 0;
my $N_STRONG_DIRECTED_INTERACTION = 0;
my $N_STRONG_UNDIRECTED_INTERACTION = 0;
my $N_STRONG_DIRECTED_LONG_RANGE_INTERACTION = 0;

my $N_SIMPLE_WEAK_READ_PAIRS = 0;
my $N_TWISTED_WEAK_READ_PAIRS = 0;
my $N_SIMPLE_STRONG_READ_PAIRS = 0;
my $N_TWISTED_STRONG_READ_PAIRS = 0;

my $N_EXCLUDED = 0;
my $N_INCLUDED = 0;

while (my $row = <$fh>) {

    chomp $row;

    my @TMP = split(/\t/,$row);

    my $frag1_chr = $TMP[0];
    my $frag1_sta = $TMP[1];
    my $frag1_end = $TMP[2];
    my $frag1_act = $TMP[3];

    my $frag2_chr = $TMP[4];
    my $frag2_sta = $TMP[5];
    my $frag2_end = $TMP[6];
    my $frag2_act = $TMP[7];

    my @TMP2 = split(/:/,$TMP[8]);

    my $n_simple = $TMP2[0];
    my $n_twisted = $TMP2[1];

    if($n_simple + $n_twisted < 2) {$N_EXCLUDED++; next;} # restrict all analyses to interactions that might be directed
    if($frag1_chr eq "chrM" || $frag2_chr eq "chrM") {$N_EXCLUDED++; next;} # exclude interactions involving chrM
    else { $N_INCLUDED++; }

    my $directed = 0;
    if((1 < $n_simple + $n_twisted) && ((2*$n_simple <= $n_twisted) || (2*$n_twisted <= $n_simple))) {
        $directed = 1;
    }

    if(1 < $n_simple + $n_twisted) {
        $N_AT_LEASTS_TWO_READ_PAIR_INTERACTION++;
        $SUM_AT_LEASTS_TWO_READ_PAIR_INTERACTION = $SUM_AT_LEASTS_TWO_READ_PAIR_INTERACTION + $n_simple + $n_twisted;
    }

    if($frag1_chr eq $frag2_chr) {
        $N_CIS_SIMPLE_READ_PAIRS = $N_CIS_SIMPLE_READ_PAIRS + $n_simple;
        $N_CIS_TWISTED_READ_PAIRS = $N_CIS_TWISTED_READ_PAIRS + $n_twisted;

        my $dist = getInteractionDist($frag1_chr, $frag1_sta, $frag1_end, $frag2_chr, $frag2_sta, $frag2_end);
        if($dist <= $LONG_RANGE_CUTOFF) {
            $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS = $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS + $n_simple;
            $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS = $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS + $n_twisted;
            if($directed) {
                $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION++;
            } else {
                $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION++;
            }
        } else {
            $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS = $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS + $n_simple;
            $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS = $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS + $n_twisted;
            if($directed) {
                $N_CIS_LONG_RANGE_DIRECTED_INTERACTION++;
            } else {
                $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION++;
            }
        }
        if($directed) {
            $N_CIS_DIRECTED_INTERACTION++;
        } else {
            $N_CIS_UNDIRECTED_INTERACTION++;
        }
    } else {
        $N_TRANS_SIMPLE_READ_PAIRS = $N_TRANS_SIMPLE_READ_PAIRS + $n_simple;
        $N_TRANS_TWISTED_READ_PAIRS = $N_TRANS_TWISTED_READ_PAIRS + $n_twisted;

        if($directed) {
            $N_TRANS_DIRECTED_INTERACTION++;
            print $directed_trans $frag1_chr . "\t" . $frag1_sta . "\t" . $frag1_end . "\tA:" . $n_simple . ":" . $n_twisted . ":" . $N_TRANS_DIRECTED_INTERACTION . "\n";
            print $directed_trans $frag2_chr . "\t" . $frag2_sta . "\t" . $frag2_end . "\tB:" . $n_simple . ":" . $n_twisted . ":" . $N_TRANS_DIRECTED_INTERACTION . "\n";
        } else {
            if($N_TRANS_UNDIRECTED_INTERACTION <= $N_TRANS_DIRECTED_INTERACTION) { # also write out undirected interactions as control
                print $undirected_trans $frag1_chr . "\t" . $frag1_sta . "\t" . $frag1_end . "\tA:" . $n_simple . ":" . $n_twisted . ":" . $N_TRANS_UNDIRECTED_INTERACTION . "\n";
                print $undirected_trans $frag2_chr . "\t" . $frag2_sta . "\t" . $frag2_end . "\tB:" . $n_simple . ":" . $n_twisted . ":" . $N_TRANS_UNDIRECTED_INTERACTION . "\n";
            }
            $N_TRANS_UNDIRECTED_INTERACTION++;
        }
    }
}

$AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION = $SUM_AT_LEASTS_TWO_READ_PAIR_INTERACTION/$N_AT_LEASTS_TWO_READ_PAIR_INTERACTION;

close $fh;
open($fh, "<", $ARGV[0])
	or die "Can't open < input.txt: $!";

while (my $row = <$fh>) {

    chomp $row;

    my @TMP = split(/\t/,$row);

    my $frag1_chr = $TMP[0];
    my $frag1_sta = $TMP[1];
    my $frag1_end = $TMP[2];
    my $frag1_act = $TMP[3];

    my $frag2_chr = $TMP[4];
    my $frag2_sta = $TMP[5];
    my $frag2_end = $TMP[6];
    my $frag2_act = $TMP[7];

    my @TMP2 = split(/:/,$TMP[8]);

    my $n_simple = $TMP2[0];
    my $n_twisted = $TMP2[1];

    my $dist=-1;
    if($frag1_chr eq $frag2_chr) {
        $dist = getInteractionDist($frag1_chr, $frag1_sta, $frag1_end, $frag2_chr, $frag2_sta, $frag2_end);
    }

    my $directed = 0;
    if((1 < $n_simple + $n_twisted) && ((2*$n_simple <= $n_twisted) || (2*$n_twisted <= $n_simple))) {
        $directed = 1;
    }

    if($frag1_chr eq "chrM" || $frag2_chr eq "chrM") {next;} # exclude interactions involving chrM

    if(2 <= $n_simple + $n_twisted) { # restrict analysis to interaction that might be directed

        if($AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION < $n_simple + $n_twisted) {

            if($directed) {
                $N_STRONG_DIRECTED_INTERACTION++;
                if($dist!=-1 && $LONG_RANGE_CUTOFF<$dist) {
                    $N_STRONG_DIRECTED_LONG_RANGE_INTERACTION++;
                    if(7<$n_simple + $n_twisted){
                        print $directed_strong_long $frag1_chr . "\t" . $frag1_sta . "\t" . $frag1_end . "\tA:" . $n_simple . ":" . $n_twisted . ":" . $N_STRONG_DIRECTED_LONG_RANGE_INTERACTION . "\n";
                        print $directed_strong_long $frag2_chr . "\t" . $frag2_sta . "\t" . $frag2_end . "\tB:" . $n_simple . ":" . $n_twisted . ":" . $N_STRONG_DIRECTED_LONG_RANGE_INTERACTION . "\n";
                    }
                }
            } else {
                $N_STRONG_UNDIRECTED_INTERACTION++;
                if($N_STRONG_UNDIRECTED_INTERACTION<$N_STRONG_DIRECTED_LONG_RANGE_INTERACTION) {
                    if(7<$n_simple + $n_twisted){
                        print $undirected_strong_long $frag1_chr . "\t" . $frag1_sta . "\t" . $frag1_end . "\tA:" . $n_simple . ":" . $n_twisted . ":" . $N_STRONG_UNDIRECTED_INTERACTION . "\n";
                        print $undirected_strong_long $frag2_chr . "\t" . $frag2_sta . "\t" . $frag2_end . "\tB:" . $n_simple . ":" . $n_twisted . ":" . $N_STRONG_UNDIRECTED_INTERACTION . "\n";
                        }
                }
            }

            $N_SIMPLE_STRONG_READ_PAIRS = $N_SIMPLE_STRONG_READ_PAIRS + $n_simple;
            $N_TWISTED_STRONG_READ_PAIRS = $N_TWISTED_STRONG_READ_PAIRS + $n_twisted;

        } else {

            if($directed) {
                $N_WEAK_DIRECTED_INTERACTION++;
            } else {
                $N_WEAK_UNDIRECTED_INTERACTION++;
            }

            $N_SIMPLE_WEAK_READ_PAIRS = $N_SIMPLE_WEAK_READ_PAIRS + $n_simple;
            $N_TWISTED_WEAK_READ_PAIRS = $N_TWISTED_WEAK_READ_PAIRS + $n_twisted;
        }
    }
}

# Finished counting

print("==============================================================================\n\n");

print("An interaction is defined to be directed, if it consists of more than one read pair,\n");
print("and either the number of simple is twice as much as the number of twisted read pairs,\n");
print("or the otherway around. I.e. 2:0 and 0:2 are the smallest directed interactions.\n\n");

print("This analysis is restricted to interactions with at least two read pairs,\n");
print("because a directed interaction requires at least two read pairs.\n\n");

print("N_EXCLUDED: $N_EXCLUDED (Interactions with only one read pair cannot be directed.)\n");
print("N_INCLUDED: $N_INCLUDED\n\n");

print("==============================================================================\n");
print("1A. Simple and twisted read pairs within cis and trans interactions:\n\n");

print("N_CIS_SIMPLE_READ_PAIRS: $N_CIS_SIMPLE_READ_PAIRS\n");
print("N_CIS_TWISTED_READ_PAIRS: $N_CIS_TWISTED_READ_PAIRS\n");
print("N_TRANS_SIMPLE_READ_PAIRS: $N_TRANS_SIMPLE_READ_PAIRS\n");
print("N_TRANS_TWISTED_READ_PAIRS: $N_TRANS_TWISTED_READ_PAIRS\n\n");

my $N_READ_PAIRS_TOTAL = $N_CIS_SIMPLE_READ_PAIRS + $N_CIS_TWISTED_READ_PAIRS + $N_TRANS_SIMPLE_READ_PAIRS + $N_TRANS_TWISTED_READ_PAIRS;

my $N_CIS_READ_PAIRS_TOTAL = $N_CIS_SIMPLE_READ_PAIRS + $N_CIS_TWISTED_READ_PAIRS;
my $N_TRANS_READ_PAIRS_TOTAL = $N_TRANS_SIMPLE_READ_PAIRS + $N_TRANS_TWISTED_READ_PAIRS;

print("N_READ_PAIRS_TOTAL: $N_READ_PAIRS_TOTAL\n");
print("N_CIS_READ_PAIRS_TOTAL: $N_CIS_READ_PAIRS_TOTAL\n");
print("N_TRANS_READ_PAIRS_TOTAL: $N_TRANS_READ_PAIRS_TOTAL\n\n");

print("The distribution of simple and twisted read pairs seems to be independent of cis and trans.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("ReadPairs1A <- as.table(rbind(c($N_CIS_SIMPLE_READ_PAIRS, $N_CIS_TWISTED_READ_PAIRS), c($N_TRANS_SIMPLE_READ_PAIRS, $N_TRANS_TWISTED_READ_PAIRS)))\n");
print("dimnames(ReadPairs1A) <- list(location = c(\"Cis\", \"Trans\"),orientation = c(\"Simple\",\"Twisted\"))\n");
print("fisher.test(ReadPairs1A, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("Xsq <- chisq.test(ReadPairs1A,correct=FALSE)\n\n");

print("All tests fail to discard the null hypothesis probably because of large sample size. However, the independence is obvious:\n");

my $sumOfCis = $N_CIS_SIMPLE_READ_PAIRS + $N_CIS_TWISTED_READ_PAIRS;
my $sumOfTrans = $N_TRANS_SIMPLE_READ_PAIRS + $N_TRANS_TWISTED_READ_PAIRS;

my $F_CIS_SIMPLE_READ_PAIRS = $N_CIS_SIMPLE_READ_PAIRS/$sumOfCis;
my $F_CIS_TWISTED_READ_PAIRS = $N_CIS_TWISTED_READ_PAIRS/$sumOfCis;
my $F_TRANS_SIMPLE_READ_PAIRS = $N_TRANS_SIMPLE_READ_PAIRS/$sumOfTrans;
my $F_TRANS_TWISTED_READ_PAIRS = $N_TRANS_TWISTED_READ_PAIRS/$sumOfTrans;

print("F_CIS_SIMPLE_READ_PAIRS: $F_CIS_SIMPLE_READ_PAIRS\n");
print("F_CIS_TWISTED_READ_PAIRS: $F_CIS_TWISTED_READ_PAIRS\n\n");
print("F_TRANS_SIMPLE_READ_PAIRS: $F_TRANS_SIMPLE_READ_PAIRS\n");
print("F_TRANS_TWISTED_READ_PAIRS: $F_TRANS_TWISTED_READ_PAIRS\n\n");

print("------------------------------------------------------------------------------\n");
print("1B. Directed and undirected interactions within cis and trans:\n\n");

print("N_CIS_DIRECTED_INTERACTION: $N_CIS_DIRECTED_INTERACTION\n");
print("N_TRANS_DIRECTED_INTERACTION: $N_TRANS_DIRECTED_INTERACTION\n");
print("N_CIS_UNDIRECTED_INTERACTION: $N_CIS_UNDIRECTED_INTERACTION\n");
print("N_TRANS_UNDIRECTED_INTERACTION: $N_TRANS_UNDIRECTED_INTERACTION\n\n");

my $N_CIS = $N_CIS_DIRECTED_INTERACTION + $N_CIS_UNDIRECTED_INTERACTION;
my $N_TRANS = $N_TRANS_DIRECTED_INTERACTION + $N_TRANS_UNDIRECTED_INTERACTION;
my $N_DIRECTED = $N_CIS_DIRECTED_INTERACTION + $N_TRANS_DIRECTED_INTERACTION;
my $N_UNDIRECTED = $N_CIS_UNDIRECTED_INTERACTION + $N_TRANS_UNDIRECTED_INTERACTION;
my $N_TOTAL = $N_CIS + $N_TRANS;

print("N_CIS: $N_CIS\n");
print("N_TRANS: $N_TRANS\n");
print("N_DIRECTED: $N_DIRECTED\n");
print("N_UNDIRECTED: $N_UNDIRECTED\n");
print("N_TOTAL: $N_TOTAL\n\n");

my $fraction_directed_cis = $N_CIS_DIRECTED_INTERACTION / ($N_CIS_DIRECTED_INTERACTION + $N_CIS_UNDIRECTED_INTERACTION);
my $fraction_directed_trans = $N_TRANS_DIRECTED_INTERACTION / ($N_TRANS_DIRECTED_INTERACTION + $N_TRANS_UNDIRECTED_INTERACTION);

print("Fraction of directed interactions within cis: $fraction_directed_cis\n");
print("Fraction of directed interactions within trans: $fraction_directed_trans\n\n");

my $factor = $fraction_directed_cis/$fraction_directed_trans;

print("Within cis interactions, directed interactions occur $factor times more often than in trans interactions,\n");
print("which, is on the suspicion of containing more artificial random interactions introduced through cross-ligation than cis.\n\n");

my $fraction_undirected_cis = $N_CIS_UNDIRECTED_INTERACTION / ($N_CIS_UNDIRECTED_INTERACTION + $N_CIS_DIRECTED_INTERACTION);
my $fraction_undirected_trans = $N_TRANS_UNDIRECTED_INTERACTION / ($N_TRANS_UNDIRECTED_INTERACTION + $N_TRANS_DIRECTED_INTERACTION);

print("Fraction of _undirected_ interactions within cis: $fraction_undirected_cis\n");
print("Fraction of _undirected_ interactions within trans: $fraction_undirected_trans\n\n");

$factor = $fraction_undirected_cis/$fraction_undirected_trans;

print("Within cis, _undirected_ interactions occur $factor times more often than in trans.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("Interactions1B <- as.table(rbind(c($N_CIS_DIRECTED_INTERACTION, $N_CIS_UNDIRECTED_INTERACTION), c($N_TRANS_DIRECTED_INTERACTION, $N_TRANS_UNDIRECTED_INTERACTION)))\n");
print("dimnames(Interactions1B) <- list(location = c(\"Cis\", \"Trans\"),orientation = c(\"Directed\",\"Undirected\"))\n");
print("fisher.test(Interactions1B, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("chisq.test(Interactions1B,correct=FALSE)\n\n");

print("==============================================================================\n");
print("2A. Simple and twisted read pairs within short and long range interactions:\n\n");

print("A read pair is defined to be in short range,\n");
print("if the distance between the corresponding restriction fragments (shortest) is less or equal than $LONG_RANGE_CUTOFF.\n\n");

print("N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS: $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS\n");
print("N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS: $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS\n");
print("N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS: $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS\n");
print("N_CIS_TWISTED_LONG_RANGE_READ_PAIRS: $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS\n\n");

#$N_CIS_READ_PAIRS_TOTAL = $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS + $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS + $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS + $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS;
#print("N_CIS_READ_PAIRS_TOTAL: $N_CIS_READ_PAIRS_TOTAL\n\n");

print("The distribution of simple and twisted read pairs seems to be independent of long and short range.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("ReadPairs2A <- as.table(rbind(c($N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS, $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS), c($N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS, $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS)))\n");
print("dimnames(ReadPairs2A) <- list(range = c(\"Short\", \"Long\"),orientation = c(\"Simple\",\"Twisted\"))\n");
print("fisher.test(ReadPairs2A, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("chisq.test(ReadPairs2A,correct=FALSE)\n\n");

print("All tests fail to discard the null hypothesis probably because of large sample size. There seems to be a slight preference for twisted read pairs within short range:\n");

my $sumOfShort = $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS + $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS;
my $sumOfLong = $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS + $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS;

my $F_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS = $N_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS/$sumOfShort;
my $F_CIS_TWISTED_SHORT_RANGE_READ_PAIRS = $N_CIS_TWISTED_SHORT_RANGE_READ_PAIRS/$sumOfShort;
my $F_CIS_SIMPLE_LONG_RANGE_READ_PAIRS = $N_CIS_SIMPLE_LONG_RANGE_READ_PAIRS/$sumOfLong;
my $F_CIS_TWISTED_LONG_RANGE_READ_PAIRS = $N_CIS_TWISTED_LONG_RANGE_READ_PAIRS/$sumOfLong;

print("F_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS: $F_CIS_SIMPLE_SHORT_RANGE_READ_PAIRS\n");
print("F_CIS_TWISTED_SHORT_RANGE_READ_PAIRS: $F_CIS_TWISTED_SHORT_RANGE_READ_PAIRS\n\n");
print("F_CIS_SIMPLE_LONG_RANGE_READ_PAIRS: $F_CIS_SIMPLE_LONG_RANGE_READ_PAIRS\n");
print("F_CIS_TWISTED_LONG_RANGE_READ_PAIRS: $F_CIS_TWISTED_LONG_RANGE_READ_PAIRS\n\n");

print("------------------------------------------------------------------------------\n");
print("2B. Directed and undirected cis interactions within short and long range interactions:\n\n");

print("N_CIS_SHORT_RANGE_DIRECTED_INTERACTION: $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION\n");
print("N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION: $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION\n");
print("N_CIS_LONG_RANGE_DIRECTED_INTERACTION: $N_CIS_LONG_RANGE_DIRECTED_INTERACTION\n");
print("N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION: $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION\n\n");

my $N_SHORT = $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION + $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION;
my $N_LONG = $N_CIS_LONG_RANGE_DIRECTED_INTERACTION + $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION;
$N_DIRECTED = $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION + $N_CIS_LONG_RANGE_DIRECTED_INTERACTION;
$N_UNDIRECTED = $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION + $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION;
$N_TOTAL = $N_SHORT + $N_LONG;

print("N_SHORT: $N_CIS\n");
print("N_LONG: $N_TRANS\n");
print("N_DIRECTED: $N_DIRECTED\n");
print("N_UNDIRECTED: $N_UNDIRECTED\n");
print("N_TOTAL: $N_TOTAL\n\n");

my $fraction_directed_short_range = $N_CIS_SHORT_RANGE_DIRECTED_INTERACTION / ($N_CIS_SHORT_RANGE_DIRECTED_INTERACTION + $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION);
my $fraction_directed_long_range = $N_CIS_LONG_RANGE_DIRECTED_INTERACTION / ($N_CIS_LONG_RANGE_DIRECTED_INTERACTION + $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION);

print("Fraction of directed interactions within long range: $fraction_directed_long_range\n");
print("Fraction of directed interactions within short range: $fraction_directed_short_range\n\n");

$factor = $fraction_directed_long_range/$fraction_directed_short_range;

print("Within long range interactions, directed interactions occur $factor times more often than in short range interactions,\n");
print("which is on the suspicion of containing more random interactions as compared to long range.\n\n");

print("Unexpected result. This might be due to fewer degrees of freedom in short range.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("Interactions2B <- as.table(rbind(c($N_CIS_SHORT_RANGE_DIRECTED_INTERACTION, $N_CIS_SHORT_RANGE_UNDIRECTED_INTERACTION), c($N_CIS_LONG_RANGE_DIRECTED_INTERACTION, $N_CIS_LONG_RANGE_UNDIRECTED_INTERACTION)))\n");
print("dimnames(Interactions2B) <- list(range = c(\"Short\", \"Long\"),orientation = c(\"Directed\",\"Undirected\"))\n");
print("fisher.test(Interactions2B, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("chisq.test(Interactions2B,correct=FALSE)\n\n");



print("==============================================================================\n");
print("3A. Simple and twisted read pairs within weak and strong interactions:\n\n");

print("An interaction is defined to be weak,\n");
print("if it contains a smaller or equal number of read pairs than the average number of read pairs in interactions with at least two read pairs,\n");
print("otherwise it is defined to be strong.\n\n");

print("AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION: $AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION\n\n");

print("N_SIMPLE_WEAK_READ_PAIRS: $N_SIMPLE_WEAK_READ_PAIRS\n");
print("N_TWISTED_WEAK_READ_PAIRS: $N_TWISTED_WEAK_READ_PAIRS\n");
print("N_SIMPLE_STRONG_READ_PAIRS: $N_SIMPLE_STRONG_READ_PAIRS\n");
print("N_TWISTED_STRONG_READ_PAIRS: $N_TWISTED_STRONG_READ_PAIRS\n\n");

print("The distribution of simple and twisted read pairs seems to be independent of weak and strong.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("ReadPairs3A <- as.table(rbind(c($N_SIMPLE_WEAK_READ_PAIRS, $N_TWISTED_WEAK_READ_PAIRS), c($N_SIMPLE_STRONG_READ_PAIRS, $N_TWISTED_STRONG_READ_PAIRS)))\n");
print("dimnames(ReadPairs3A) <- list(strength = c(\"Weak\", \"Strong\"),orientation = c(\"Simple\",\"Twisted\"))\n");
print("fisher.test(ReadPairs3A, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("chisq.test(ReadPairs3A,correct=FALSE)\n\n");

print("All tests fail to discard the null hypothesis probably because of large sample size.\n");
print("There seems to be a slight preference for twisted read pairs within short range:\n\n");

my $sumOfWeak = $N_SIMPLE_WEAK_READ_PAIRS + $N_TWISTED_WEAK_READ_PAIRS;
my $sumOfStrong = $N_SIMPLE_STRONG_READ_PAIRS + $N_TWISTED_STRONG_READ_PAIRS;

my $F_SIMPLE_WEAK_READ_PAIRS = $N_SIMPLE_WEAK_READ_PAIRS/$sumOfWeak;
my $F_TWISTED_WEAK_READ_PAIRS = $N_TWISTED_WEAK_READ_PAIRS/$sumOfWeak;
my $F_SIMPLE_STRONG_READ_PAIRS = $N_SIMPLE_STRONG_READ_PAIRS/$sumOfStrong;
my $F_TWISTED_STRONG_READ_PAIRS = $N_TWISTED_STRONG_READ_PAIRS/$sumOfStrong;

print("F_SIMPLE_WEAK_READ_PAIRS: $F_SIMPLE_WEAK_READ_PAIRS\n");
print("F_TWISTED_WEAK_READ_PAIRS: $F_TWISTED_WEAK_READ_PAIRS\n\n");
print("F_SIMPLE_STRONG_READ_PAIRS: $F_SIMPLE_STRONG_READ_PAIRS\n");
print("F_TWISTED_STRONG_READ_PAIRS: $F_TWISTED_STRONG_READ_PAIRS\n\n");


print("------------------------------------------------------------------------------\n");
print("3B. Directed and undirected interactions within weak and strong interactions:\n\n");

print("N_WEAK_DIRECTED_INTERACTION: $N_WEAK_DIRECTED_INTERACTION\n");
print("N_WEAK_UNDIRECTED_INTERACTION: $N_WEAK_UNDIRECTED_INTERACTION\n");
print("N_STRONG_DIRECTED_INTERACTION: $N_STRONG_DIRECTED_INTERACTION\n");
print("N_STRONG_UNDIRECTED_INTERACTION: $N_STRONG_UNDIRECTED_INTERACTION\n\n");

my $N_WEAK = $N_WEAK_DIRECTED_INTERACTION + $N_WEAK_UNDIRECTED_INTERACTION;
my $N_STRONG = $N_STRONG_DIRECTED_INTERACTION + $N_STRONG_UNDIRECTED_INTERACTION;
$N_DIRECTED = $N_WEAK_DIRECTED_INTERACTION + $N_STRONG_DIRECTED_INTERACTION;
$N_UNDIRECTED = $N_WEAK_UNDIRECTED_INTERACTION + $N_STRONG_UNDIRECTED_INTERACTION;
$N_TOTAL = $N_WEAK + $N_STRONG;

print("N_WEAK: $N_WEAK\n");
print("N_STRONG: $N_STRONG\n");
print("N_DIRECTED: $N_DIRECTED\n");
print("N_UNDIRECTED: $N_UNDIRECTED\n");
print("N_TOTAL: $N_TOTAL\n\n");

my $fraction_directed_weak = $N_WEAK_DIRECTED_INTERACTION / ($N_WEAK_DIRECTED_INTERACTION + $N_WEAK_UNDIRECTED_INTERACTION);
my $fraction_directed_strong = $N_STRONG_DIRECTED_INTERACTION / ($N_STRONG_DIRECTED_INTERACTION + $N_STRONG_UNDIRECTED_INTERACTION);

print("Fraction of directed interactions within weak: $fraction_directed_weak\n");
print("Fraction of directed interactions within strong: $fraction_directed_strong\n\n");

$factor = $fraction_directed_strong/$fraction_directed_weak;

print("Within the strong interactions, directed interactions occur $factor times more often than in weak interactions,\n");
print("which are on the suspicion of containing more random interactions as compared to strong interactions.\n\n");

print("Fisher's exact (use these lines in R):\n");
print("Interactions3B <- as.table(rbind(c($N_WEAK_DIRECTED_INTERACTION, $N_WEAK_UNDIRECTED_INTERACTION), c($N_STRONG_DIRECTED_INTERACTION, $N_STRONG_UNDIRECTED_INTERACTION)))\n");
print("dimnames(Interactions3B) <- list(strength = c(\"Weak\", \"Strong\"),orientation = c(\"Directed\",\"Undirected\"))\n");
print("fisher.test(Interactions3B, alternative = \"two.sided\")\n\n");

print("Chi-Squared (use this line in R):\n");
print("chisq.test(Interactions3B,correct=FALSE)\n\n");

print("------------------------------------------------------------------------------\n");
print("3D. Directed and strong long range interactions:\n\n");

print("N_STRONG_DIRECTED_LONG_RANGE_INTERACTION: $N_STRONG_DIRECTED_LONG_RANGE_INTERACTION\n\n");


close $directed_trans;
close $undirected_trans;
close $directed_strong_long;
close $undirected_strong_long;

my $bed_filename = $OUT_PREFIX . "_directed_trans.bed";
my $fasta_filename_dir = $OUT_PREFIX . "_directed_trans.fasta";

# extract FASTA
system("bedtools getfasta -name -fi ../VPV_data/mm9/mm9.fa -bed $bed_filename > $fasta_filename_dir");

# replace small letters with N
my $fasta_filename_dir_no_rep = $OUT_PREFIX . "_directed_trans_no_rep.fasta";
system("sed -e '/^>/! s/[[:lower:]]/N/g' $fasta_filename_dir > $fasta_filename_dir_no_rep");

# the same for the control
$bed_filename = $OUT_PREFIX . "_undirected_trans.bed";
my $fasta_filename_undir = $OUT_PREFIX . "_undirected_trans.fasta";
my $fasta_filename_undir_no_rep = $OUT_PREFIX . "_undirected_trans_no_rep.fasta";
system("bedtools getfasta -name -fi ../VPV_data/mm9/mm9.fa -bed $bed_filename > $fasta_filename_undir");
system("sed -e '/^>/! s/[[:lower:]]/N/g' $fasta_filename_undir > $fasta_filename_undir_no_rep");

# run DREME for sequences with and without repeats
my $dreme_directory = $OUT_PREFIX . "_DREME_NOREP";
#system("dreme -oc $dreme_directory -p $fasta_filename_dir_no_rep -n $fasta_filename_undir_no_rep -eps");

$dreme_directory = $OUT_PREFIX . "_DREME";
#system("dreme -oc $dreme_directory -p $fasta_filename_dir -n $fasta_filename_undir -eps");


$bed_filename = $OUT_PREFIX . "_directed_strong_long.bed";
$fasta_filename_dir = $OUT_PREFIX . "_directed_strong_long.fasta";

# extract FASTA
system("bedtools getfasta -name -fi ../VPV_data/mm9/mm9.fa -bed $bed_filename > $fasta_filename_dir");
$fasta_filename_dir_no_rep = $OUT_PREFIX . "_directed_strong_long_no_rep.fasta";
system("sed -e '/^>/! s/[[:lower:]]/N/g' $fasta_filename_dir > $fasta_filename_dir_no_rep");

$bed_filename = $OUT_PREFIX . "_undirected_strong_long.bed";
$fasta_filename_undir = $OUT_PREFIX . "_undirected_strong_long.fasta";
$fasta_filename_undir_no_rep = $OUT_PREFIX . "_undirected_strong_long_no_rep.fasta";

system("bedtools getfasta -name -fi ../VPV_data/mm9/mm9.fa -bed $bed_filename > $fasta_filename_undir");
system("sed -e '/^>/! s/[[:lower:]]/N/g' $fasta_filename_undir > $fasta_filename_undir_no_rep");

system("dreme -oc $dreme_directory -p $fasta_filename_dir_no_rep -n $fasta_filename_undir_no_rep -eps");







