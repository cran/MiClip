#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max/;

my ($file_tag,$file_cluster)=@ARGV;
my ($temp,$count,@items,%tags,$len,$chr,@pos_start,$end,@cluster,$i,$line,$strand);
my $md_field=-1;

open(FILE_IN,$file_tag) or die "read alignment file ".$file_tag." failed!\n";

while ($line=<FILE_IN>) # search for MD field
{
  if ($line=~/^@/) {next;}

  @items=split("\t",$line);
  unless ($items[1]==0 || $items[1]==16) {next;}
  $i=0;

  foreach (@items)
  {
    if ($_=~/MD\:Z\:/) {$md_field=$i;}
    $i++;
  }

  if ($md_field>0) {last;}
}

if ($md_field<0) {die "no md field found in sam file!\nMaybe your file is not single-end?\n";}
seek FILE_IN,0,0;

while (<FILE_IN>) # read and compute length of each tag
{
  if ($_=~/^@/) {next;}
  @items=split("\t",$_);
  if ($items[5]=~/N/) {next;} # delete gapped mapping

  if ($items[1]==0 || $items[1]==16)
  {
    $len=0;
    while ($items[5]=~/([0-9]+)([MD])/g) {$len+=$1;}    

    $count=()=$items[5]=~/(I)/g;
    $count+=()=$items[$md_field]=~/(\^[A-Z]+|[ATCG])/g;
   
    $temp=$tags{$items[1]."_".$items[2]}->{$items[3]}->{$len};

    if ($items[$md_field]!~/MD/) 
    {
      $items[$md_field]=$len;
    }else
    {
      $items[$md_field]=~/MD\:Z\:(.*)$/;
      $items[$md_field]=$1;
    }

    if (exists ${$temp}[0])
    {
      if (${$temp}[0]<$count) {@$temp=($count,$items[0],$items[3],$items[5],$items[9],$items[$md_field]);} 
      # information of each tag is stored as mutant count, tag name, pos, CIGAR, sequence, MD
      # strand, chr, pos and length are stored as key values
    }else
    {
      $tags{$items[1]."_".$items[2]}->{$items[3]}->{$len}=[$count,$items[0],$items[3],$items[5],$items[9],$items[$md_field]];
    }
  }
}

close(FILE_IN);

open(FILE_OUT,">".$file_cluster) or die "write to temporary cluster file failed!\n";

foreach $temp (keys %tags) # form clusters
{
  ($strand,$chr)=split("_",$temp);
  if ($strand==0) 
  {
    $strand="+";
  }else
  {
    $strand="-";
  }

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
        map {print FILE_OUT ${$_}[1]."\t".${$_}[2]."\t".${$_}[3]."\t".${$_}[4]."\t".${$_}[5]."\n";} @cluster;
      }
      @cluster=();
      map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[$i]}};
      $end=$pos_start[$i]-1+max(map {abs($_)} keys %{$tags{$temp}->{$pos_start[$i]}});
    }
  }

  if ($#cluster>0) 
  {
    print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
    map {print FILE_OUT ${$_}[1]."\t".${$_}[2]."\t".${$_}[3]."\t".${$_}[4]."\t".${$_}[5]."\n";} @cluster;
  }
  @cluster=();
}

close(FILE_OUT);

exit;




