#!/usr/bin/perl

my $TRUNC_FILE = $ARGV[0];
my $ALIGN_FILE = $ARGV[1];
my $COUNT_FILE = $ARGV[2];

my @HEADER;
my @NUMBERS;

# parse trunc file
open(my $fh, "<", $TRUNC_FILE)
	or die "Can't open < input.txt: $!";

while (my $row = <$fh>) {
    chomp $row;
    if($row =~ m/:/) {
        my @TMP = split(/:/,$row);
        my @TMP2 = split(' ',$TMP[1]);
        chomp($TMP[0]);
        chomp($TMP2[0]);
        $TMP[0] =~ s/\t//g;
        $TMP2[0] =~ s/\t//g;
        push(@HEADER,$TMP[0]);
        push(@NUMBERS,$TMP2[0]);

    }
}
close($fh);

# parse align file
open(my $fh, "<", $ALIGN_FILE)
	or die "Can't open < input.txt: $!";

while (my $row = <$fh>) {
    chomp $row;
    if($row =~ m/:/ && $row !~ m/array/ && $row !~ m/Note/ && $row !~ m/Fractions/) {
        my @TMP = split(/:/,$row);
        my @TMP2 = split(' ',$TMP[1]);
        chomp($TMP[0]);
        chomp($TMP2[0]);
        $TMP[0] =~ s/\t//g;
        $TMP2[0] =~ s/\t//g;
        push(@HEADER,$TMP[0]);
        push(@NUMBERS,$TMP2[0]);
    }
}
close($fh);

# parse align file
open(my $fh, "<", $COUNT_FILE)
	or die "Can't open < input.txt: $!";

while (my $row = <$fh>) {
    chomp $row;
    if($row =~ m/:/ && $row !~ m/array/ && $row !~ m/Note/) {
        my @TMP = split(/:/,$row);
        my @TMP2 = split(' ',$TMP[1]);
        chomp($TMP[0]);
        chomp($TMP2[0]);
                $TMP[0] =~ s/\t//g;
                $TMP2[0] =~ s/\t//g;
        push(@HEADER,$TMP[0]);
        push(@NUMBERS,$TMP2[0]);
    }
}
close($fh);



# print header line
foreach my $n (@HEADER) {
    print("$n\t");
}
print("\n");
# print number line
foreach my $n (@NUMBERS) {
    print("$n\t");
}

print("\n");