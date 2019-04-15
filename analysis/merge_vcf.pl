#!/usr/bin/perl
##by Li Lei, 20170809, in St.Paul;
#this is to merge vcf files for Connor. Since vcf-merge can not handle the the position out of 500 MB. 
# THhe prerequest is to make sure all of the SNPs in the 1-9 columns are identical. You can run this commandline to make sure if they are the same:  diff -y <(grep -v "#" sorted_filtered_round2.vcf|cut -f 1,2,3,4,5,9) <(grep -v "#" sorted_filtered_NAM_9k.vcf|cut -f 1,2,3,4,5,9)|less -S
#usage: perl /panfs/roc/groups/9/morrellp/llei/Introgressed_line/script/merge_vcf.pl sorted_filtered_round2.vcf sorted_filtered_NAM_9k.vcf >merged_NAM_9k_round2.vcf

use strict;
use warnings;
use Data::Dumper;
my ($file1, $file2) = @ARGV;

my %SNPposHash;


open(SNPID,  "$file1") or die "Could not open $file1";

foreach my $row (<SNPID>){
       if ($row =~ /^#/) {
           next;
       }
       else{
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $SNP_id = $rtemp[2];
           $SNPposHash{$SNP_id} = $rtemp[9]."\t".$rtemp[10]."\t".$rtemp[11]."\t".$rtemp[12]."\t".$rtemp[13]."\t".$rtemp[14]."\t".$rtemp[15]."\t".$rtemp[16]."\t".$rtemp[17]."\t".$rtemp[18]."\t".$rtemp[19]."\t".$rtemp[20]."\t".$rtemp[21]."\t".$rtemp[22]."\t".$rtemp[23]."\t".$rtemp[24]."\t".$rtemp[25]."\t".$rtemp[26]."\t".$rtemp[27]."\t".$rtemp[28]."\t".$rtemp[29]."\t".$rtemp[30]."\t".$rtemp[31]."\t".$rtemp[32]."\t".$rtemp[33]."\t".$rtemp[34];
      }
        #print "$rtemp[0]\n";
}
close (SNPID);
#print Dumper(\%gidhash);

open(SIGSNP,  "$file2") or die "Could not open $file2";
foreach my $row (<SIGSNP>){
          chomp $row;
  if ($row =~ /^##/){
           print "$row\n";
  }
  elsif($row =~ /^#CHROM/){
    print "$row\t0_M109\t0_WBDC016\t0_WBDC020\t0_WBDC028\t0_WBDC032\t0_WBDC035\t0_WBDC042\t0_WBDC061\t0_WBDC082\t0_WBDC092\t0_WBDC103\t0_WBDC115\t0_WBDC142\t0_WBDC150\t0_WBDC172\t0_WBDC173\t0_WBDC182\t0_WBDC227\t0_WBDC234\t0_WBDC255\t0_WBDC292\t0_WBDC302\t0_WBDC336\t0_WBDC340A\t0_WBDC348A\t0_WBDC350\n";
  }
  else{
        my @rtemp1 = split(/\t/,$row);
           #print "$rtemp1[0]\n";
           my $SNP = $rtemp1[2];
        if(exists $SNPposHash{$SNP}){
          print "$row\t$SNPposHash{$SNP}\n"
        }
        else{
           print "$SNP\tNOT_EXISTS\n";
        }
  }
}
close (SNPID);
