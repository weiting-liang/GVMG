($taxo,$clus)=@ARGV;
open I,"$taxo";
%hash;while(<I>){chomp;@xx=split /\t/,$_; $hash{$xx[0]}{$xx[2]}=$xx[1];}close I;
open I,"$clus";
while(<I>){
        chomp;
        @xj=split /\s+/,$_;
        if(@xj==1){print "$_\n";next;}
        my @array=split/,/,$xj[-1];
        my %count=();my %iden=();
        my  $jie=0;
        foreach my $ij (@array){
                if(exists $hash{$ij}){
                        $jie++; foreach my $k (keys %{$hash{$ij}}){$count{$k}++;$iden{$k}+=$hash{$ij}{$k};}
                }
        }
        my @temp=sort {$count{$b} <=> $count{$a} or $iden{$b} <=> $iden{$a}} keys %count;
        print "$xj[0]\t$xj[-2]";  print "\t$jie";
        if($count{$temp[0]}/$xj[-2]<0.4){print "\tUnclassified\t$count{$temp[0]}\n";next;}
        print "\t$temp[0]\t$count{$temp[0]}\t",$count{$temp[0]}/$xj[-2],"\t",$iden{$temp[0]}/$count{$temp[0]};
        for my $k(1..$#temp){
                if($count{$temp[$k-1]}/$xj[-2]-$count{$temp[$k]}/$xj[-2]>0.2){last;}
                print "\t$temp[$k]\t$count{$temp[$k]}\t",$count{$temp[$k]}/$xj[-2],"\t",$iden{$temp[$k]}/$count{$temp[$k]};
        }
        my $jie=0; for my $k(1..$#temp){$jie=$jie+ $count{$temp[$k]}} print "$jie\n";
}
close I;
