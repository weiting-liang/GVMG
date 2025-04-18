use warnings;
use strict;

die "Usage: perl $0 [fasta] [out]\n" unless @ARGV == 2;

my ($in_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

if ($in_f =~ /\.gz$/) { open IN, "gzip -dc $in_f |" or die $!;
} else { open IN, $in_f or die $!;
}
open OT,">$out_f" or die $!;

my ($head, $len, $non) = ('', 0, 0);

while (<IN>) {
        chomp;

        if (/^>(\S+)/) {
                print OT "$head\t$len\t$non\n" if $len > 0;
                ($head, $len, $non) = ($1, 0, 0);
                next;
        }
        $len += length $_;
        s/N|n//g; $non += length $_;
}
print OT "$head\t$len\t$non\n" if $len > 0;
close IN;
close OT;
