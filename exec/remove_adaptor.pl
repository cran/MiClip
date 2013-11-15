#!/usr/bin/perl
use strict;
use warnings;
$SIG{__WARN__}=sub{die "Sequencing file format error!\n"};

my ($format,$file,$adaptor,$len_limit,$mismatch)=@ARGV;
my $ad_len=0;
my ($i,$j,$k,$dist,@char,@read,%regex,$new_len,$mis_j);

if (length($adaptor)<=10) {die "The length of the adaptor is too short!\n";}

#####  construct adaptor hash  #####

$adaptor=uc $adaptor;
$adaptor=~s/ //g;

while ($adaptor=~/([ATCGN])/g)
{
  if ($1 eq "N") 
  {
    $regex{$ad_len}={"N"=>1,"A"=>1,"T"=>1,"C"=>1,"G"=>1};
  }else
  {
    $regex{$ad_len}={"N"=>1,$1=>1};
  }
  $ad_len++;
}

#####  transform file  #####

if ($format eq "fastq")
{
  open(FILE_IN,$file) or die "Cannot open the sequencing file ".$file."!\n";
  open(FILE_OUT,">".$file.".removed") or die "Cannot write to the folder of the sequencing file!\n";

  $i=0;
  @read=();

  while (<FILE_IN>)
  {
    $i++;
    push @read,$_;
  
    if ($i % 4==0)
    {
      @char=split(//,substr($read[1],length($read[1])-1-$ad_len),$ad_len);
      $char[$#char]=~s/\n//;
      $new_len=length($read[1]);

LINE1: foreach $j (0..$#char)
      {
        $dist=0;
        $mis_j=$mismatch*($#char-$j+1);

        #print "j: ".$j." end:".$#char." mis:".$mis_j."\n";

        foreach $k ($j..$#char)
        {
          unless (exists $regex{$k-$j}->{$char[$k]}) 
          {
            $dist++;
            if ($dist>$mis_j) {next LINE1;}
          }
          #print " k:".$k." char:".$char[$k]." dist:".$dist."\n";
        }

        $new_len=length($read[1])-1-$ad_len+$j;

        if ($new_len>=$len_limit)
        {
          $read[0]=~s/=[0-9]+\n/=$new_len\n/;
          $read[2]=~s/=[0-9]+\n/=$new_len\n/;
          $read[1]=substr($read[1],0,$new_len)."\n";
          $read[3]=substr($read[3],0,$new_len)."\n";
          last;
        }
      }
   
      if ($new_len>=$len_limit) {print FILE_OUT @read;}
      @read=();
    }
  } 

  close(FILE_IN);
  close(FILE_OUT);
} elsif ($format eq "fasta")
{
  open(FILE_IN,$file) or die "Cannot open the sequencing file ".$file."!\n";
  open(FILE_OUT,">".$file.".removed") or die "Cannot write to the folder of the sequencing file!\n";

  $i=0;
  @read=();

  while (<FILE_IN>)
  {
    $i++;
    push @read,$_;

    if ($i % 2==0)
    {
      @char=split(//,substr($read[1],length($read[1])-1-$ad_len),$ad_len);
      $char[$#char]=~s/\n//;
      $new_len=length($read[1]);

LINE2: foreach $j (0..$#char)
      {
        $dist=0;
        $mis_j=$mismatch*($#char-$j+1);

        foreach $k ($j..$#char)
        {
          unless (exists $regex{$k-$j}->{$char[$k]})
          {
            $dist++;
            if ($dist>$mis_j) {next LINE2;}
          }
        }

        $new_len=length($read[1])-1-$ad_len+$j;

        if ($new_len>=$len_limit)
        {
          $read[1]=substr($read[1],0,$new_len)."\n";
          last;
        }
      }

      if ($new_len>=$len_limit) {print FILE_OUT @read;}
      @read=();
    }
  }

  close(FILE_IN);
  close(FILE_OUT);
}else
{
  print "Format of sequencing file is not recognized, or this format can not be processed by this helper function!\n";
}








