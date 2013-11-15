#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/floor/;
use List::Util qw/min max/;
$SIG{__WARN__}=sub{die "SAM file format error!\n"};

# $step is the step resolution desired
# $file_cluster is the combined cluster file
# $mut_type is the type of mutation wanted, bracketed by "" and separated by ",", e.g. "T2C,A2G"
# $file_count is the output file

my ($step,$file_cluster,$mut_type,$file_count)=@ARGV;
my ($chr,$strand,$tag_number,%mut_type);
my (@items,@starts,@ends,@muts,@bins,@bins_mut);
my ($start,$end,$i,$j);
my $region_id=1;

open(FILE_IN,$file_cluster) or die "read cluster file ".$file_cluster." failed!\n";
open(FILE_OUT,"> ".$file_count) or die "can't write to file ".$file_count."!\n";

mut_hash($mut_type,\%mut_type); # generate mutation type hash

while (<FILE_IN>)
{
  ($chr,$strand,$tag_number)=split("\t",$_);
  
  @muts=();
  @starts=();
  @ends=();
  @bins=();

  foreach $i (1..$tag_number)
  {
    @items=split("\t",<FILE_IN>);
    read_tag($mut_type,$strand,$items[1],$items[2],$items[3],$items[4],\@starts,\@ends,\@muts);
  }
    
  $start=$step*floor(min(@starts)/$step);  # define cluster region
  $end=$step*(floor(max(@ends)/$step)+1)-1;

  foreach ($start..$end) {$bins[$_-$start]=0;} # calculate base coverage on each base

  foreach $i (0..($tag_number-1)) 
  {
    foreach (($starts[$i]-$start)..($ends[$i]-$start)) {$bins[$_]++;}
  }

  foreach ($start..$end) {$bins_mut[$_-$start]=0;} # calcualte mutant base coverage on each base

  foreach (@muts)
  {
    $bins_mut[$_-$start]++;
  }

  foreach $i (0..(($end-$start+1)/$step-1)) 
  {
    print FILE_OUT $region_id."\t".$chr."\t".($i*$step+$start)."\t".$strand."\t";
      
    @items=();
    foreach $j (($i*$step)..($i*$step+$step-1)) {push @items,$bins[$j];}
    print FILE_OUT join("\t",@items)."\t";

    @items=();
    foreach $j (($i*$step)..($i*$step+$step-1)) {push @items,$bins_mut[$j];}
    print FILE_OUT join("\t",@items)."\n";      
  }  # print

  $region_id++;
}

close(FILE_IN);
close(FILE_OUT);

sub read_tag # read tag file, only maintain tags in the %tags hash
{
  my ($mut_type,$strand,$pos,$CIGAR,$seq,$MD,$starts_ref,$ends_ref,$muts_ref)=@_;  
  my (@items,$len,@mut_pos);

  $len=0;
  while ($CIGAR=~/([0-9]+)([MD])/g) {$len+=$1;}
  push @$starts_ref,$pos;
  push @$ends_ref,($pos+$len-1); # record position information       
 
  @mut_pos=();
  read_mut($mut_type,$CIGAR,$seq,$MD,\@mut_pos,$strand);
  map {$_+=($pos-1);} @mut_pos;
  push @$muts_ref,@mut_pos; # record mutation information
}

sub read_mut
{
  my ($mut_type,$CIGAR,$seq,$MD,$mut_pos_ref,$strand)=@_;
  my @mut_pos=();
  my $ref_pos=0;
  my $tag_pos=0;
  my ($match,$temp);

  if ($CIGAR=~/([0-9]+)S.*[0-9]+M/) {$seq=substr($seq,$1);} # offset soft-clipping
  $CIGAR=~s/[0-9]+H//g; # offset hard-clipping

  while ($CIGAR=~/([0-9]+)([MDI])/g)
  {
    if ($2 eq "M") 
    {
      $ref_pos+=$1;
      $tag_pos+=$1;
    }elsif ($2 eq "I")
    {
      {if (exists $mut_type{"Ins"}) {push @$mut_pos_ref,($ref_pos+1);}}
      substr($seq,$tag_pos,$1)="";
    }else
    { 
      {if (exists $mut_type{"Del"}) {push @$mut_pos_ref,($ref_pos+1);}}
      $ref_pos+=$1;
    }
  }  

  $ref_pos=0;
  $tag_pos=0;

  while ($MD=~/([0-9]+|[ACGTN]|\^[ACGTN]+)/g) 
  {
    $match=$1;
    if ($match=~/[0-9]+/) 
    {
      $ref_pos+=$match;
      $tag_pos+=$match;      
    }elsif ($match=~/^[ACGTN]$/)
    {
      $ref_pos+=1;
      $temp=substr($seq,$tag_pos,1);

      if ($strand eq "-")
      {
        $match=transform($match);
        $temp=transform($temp);
      }

      $tag_pos+=1;
      if (exists $mut_type{$match."2".$temp}) {push @$mut_pos_ref,$ref_pos;}
    }else
    {
      $ref_pos+=length($match)-1;
    }
  }
}

sub transform # negative strand to positive strand
{
  my $base=$_[0];
  
  if ($base eq "A")
  {
    $base="T";
  }elsif ($base eq "T")
  {
    $base="A";
  }elsif ($base eq "C")
  {
    $base="G";
  }elsif ($base eq "G")
  {
    $base="C";
  }

  return($base);
}

sub mut_hash
{
  my ($mut_type,$mut_type_ref)=@_;

  $mut_type=~s/ //g;
  $mut_type=~s/\n//;

  if ($mut_type=~/all/)
  {
    %$mut_type_ref=("A2C"=>1,"A2T"=>1,"A2G"=>1,"T2A"=>1,"T2G"=>1,"T2C"=>1,"C2A"=>1,"C2T"=>1,"C2G"=>1,"G2A"=>1,"G2T"=>1,"G2C"=>1,"Del"=>1,"Ins"=>1);
  }else
  {
    map {$mut_type_ref->{$_}=1;} split(",",$mut_type);
  }
}









