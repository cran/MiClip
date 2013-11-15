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
my (@items,@starts,@ends,@muts,@bins,@bins_mut,@f,@r);
my ($start,$end,$i,$j);
my $region_id=1;

open(FILE_IN,$file_cluster) or die "read cluster file ".$file_cluster." failed!\n";
open(FILE_OUT,"> ".$file_count) or die "can't write to file ".$file_count."!\n";

mut_hash($mut_type,\%mut_type);

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
    @f=@items[1..4];
    @r=@items[5..8];
    read_tag($mut_type,$strand,\@f,\@r,\@starts,\@ends,\@muts);
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
  my ($mut_type,$strand,$f_ref,$r_ref,$starts_ref,$ends_ref,$muts_ref)=@_;  
  my ($len_f,$len_r,@mut_f,@mut_r,$start,$end);

  $len_f=0;
  while ($$f_ref[1]=~/([0-9]+)([MD])/g) {$len_f+=$1;}
  $len_r=0;
  while ($$r_ref[1]=~/([0-9]+)([MD])/g) {$len_r+=$1;}

  $start=min($$f_ref[0],$$r_ref[0]);
  $end=max($$f_ref[0]+$len_f-1,$$r_ref[0]+$len_r-1);

  push @$starts_ref,$start;
  push @$ends_ref,$end; # record position information       
 
  $strand=($strand eq "+")?1:0;

  @mut_f=();
  read_mut($mut_type,$f_ref,\@mut_f,$strand);
  map {$_+=($$f_ref[0]-1);} @mut_f;

  @mut_r=();
  read_mut($mut_type,$r_ref,\@mut_r,$strand);
  map {$_+=($$r_ref[0]-1);} @mut_r;  

  push @$muts_ref,(uniq (@mut_f,@mut_r)); # record mutation information
}

sub read_mut
{
  my ($mut_type,$tag_ref,$mut_ref,$strand)=@_;
  my @mut_pos=();
  my $ref_pos=0;
  my $tag_pos=0;
  my ($match,$temp);

  my $CIGAR=$$tag_ref[1];
  my $seq=$$tag_ref[2];
  my $MD=$$tag_ref[3];

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
      {if (exists $mut_type{"Ins"}) {push @$mut_ref,($ref_pos+1);}}
      substr($seq,$tag_pos,$1)="";
      $tag_pos+=$1;
    }else
    { 
      {if (exists $mut_type{"Del"}) {push @$mut_ref,($ref_pos+1);}}
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

      if ($strand==0)
      {
        $match=transform($match);
        $temp=transform($temp);
      }

      $tag_pos+=1;
      if (exists $mut_type{$match."2".$temp}) {push @$mut_ref,$ref_pos;}
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

sub uniq
{
  my %temp=();
  map {$temp{$_}=1} @_;
  keys %temp;
}


