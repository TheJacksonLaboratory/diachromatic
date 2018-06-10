#!/usr/bin/perl

use strict;
use warnings;

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

my $N_SIMPLE_WEAK_READ_PAIRS = 0;
my $N_TWISTED_WEAK_READ_PAIRS = 0;
my $N_SIMPLE_STRONG_READ_PAIRS = 0;
my $N_TWISTED_STRONG_READ_PAIRS = 0;

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

    if($n_simple + $n_twisted < 2) {next;} # restrict all analyses to interactions that might be directed

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
        } else {
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

    my $directed = 0;
    if((1 < $n_simple + $n_twisted) && ((2*$n_simple <= $n_twisted) || (2*$n_twisted <= $n_simple))) {
        $directed = 1;
    }

    if(2 <= $n_simple + $n_twisted) { # restrict analysis to interaction that might be directed

        if($AVG_AT_LEASTS_TWO_READ_PAIR_INTERACTION < $n_simple + $n_twisted) {

            if($directed) {
                $N_STRONG_DIRECTED_INTERACTION++;
            } else {
                $N_STRONG_UNDIRECTED_INTERACTION++;
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

print("A interaction is defined to be directed, if it consists of more than one read pair,\n");
print("and either the number of simple is twice as much as the number of twisted read pairs,\n");
print("or the othernway around.\n\n");

print("N_CIS_DIRECTED_INTERACTION: $N_CIS_DIRECTED_INTERACTION\n");
print("N_TRANS_DIRECTED_INTERACTION: $N_TRANS_DIRECTED_INTERACTION\n");
print("N_CIS_UNDIRECTED_INTERACTION: $N_CIS_UNDIRECTED_INTERACTION\n");
print("N_TRANS_UNDIRECTED_INTERACTION: $N_TRANS_UNDIRECTED_INTERACTION\n\n");

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

print("This analysis is restricted to interactions with at least two read pairs,\n");
print("because a directed interaction requires at least two read pairs (should this also be done for cis/trans and short/long?).\n\n");

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