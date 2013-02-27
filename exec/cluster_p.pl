#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min/;

my ($file_merge,$file_cluster)=@ARGV;
my ($len_f,$len_r,$start,$end,$count,@items,%tags);
my ($temp,$chr,@pos_start,@cluster,$i,$line,$strand);

open(FILE_IN,$file_merge) or die "read alignment file ".$file_merge." failed!\n";

while (<FILE_IN>) # read and compute length of each tag
{
  @items=split("\t",$_);

  $len_f=0;
  while ($items[4]=~/([0-9]+)([MD])/g) {$len_f+=$1;}    
  $len_r=0;
  while ($items[8]=~/([0-9]+)([MD])/g) {$len_r+=$1;}
  
  if ($items[6]=~/-1/) {$items[6]=$len_f;}
  if ($items[10]=~/-1/) {$items[10]=$len_r."\n";}

  $start=min($items[3],$items[7]);
  $end=max($items[3]+$len_f-1,$items[7]+$len_r-1);  

  $count=()=$items[4]=~/(I)/g;
  $count+=()=$items[6]=~/(\^[A-Z]+|[ATCG])/g;
  $count+=()=$items[8]=~/(I)/g;
  $count+=()=$items[10]=~/(\^[A-Z]+|[ATCG])/g;

  $temp=$tags{$items[1]."_".$items[2]}->{$start}->{$end-$start+1};
 
  # information is stored as mutant count, tag name, pos_f, CIGAR_f, seq_f, MD_f, pos_r, CIGAR_r, seq_r, MD_r 
  if (exists ${$temp}[0])
  {
    if (${$temp}[0]<$count) {@$temp=($count,@items[0,3..10]);} 
  }else
  {
    $tags{$items[1]."_".$items[2]}->{$start}->{$end-$start+1}=[$count,@items[0,3..10]];
  }
}

close(FILE_IN);
unlink($file_merge);

open(FILE_OUT,">".$file_cluster) or die "write to temporary cluster file failed!\n";

foreach $temp (keys %tags) # form clusters
{
  ($chr,$strand)=split("_",$temp);

  @pos_start=keys %{$tags{$temp}};
  @pos_start=sort {$a <=> $b} @pos_start; # @pos_start contains sorted start sites
  map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[0]}}; # @cluster contains the name of the tags that form one single cluster
  $end=$pos_start[0]-1+max(keys %{$tags{$temp}->{$pos_start[0]}});

  for $i (1..$#pos_start)
  {
    if ($end>=$pos_start[$i])
    {
      $end=max($end,$pos_start[$i]-1+max(keys %{$tags{$temp}->{$pos_start[$i]}}));
      map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[$i]}};
    }else
    {
      if ($#cluster>0) 
      {
        print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
        map {print FILE_OUT join("\t",@{$_}[1..9]);} @cluster;
      }
      @cluster=();
      map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[$i]}};
      $end=$pos_start[$i]-1+max(map {abs($_)} keys %{$tags{$temp}->{$pos_start[$i]}});
    }
  }

  if ($#cluster>0) 
  {
    print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
    map {print FILE_OUT join("\t",@{$_}[1..9]);} @cluster;
  }
  @cluster=();
}

close(FILE_OUT);

exit;




