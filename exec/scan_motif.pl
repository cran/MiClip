#!/usr/bin/perl
#use strict;
use warnings;

# The first parameter is the sequence file. Each line has a name and a sequence, separated by a tab
# The second parameter is the motif. All in capital letters.
# The third parameter is the length of the motif.
# The fourth parameter is the name of the output file. 

my ($seq_file,$motif,$length,$fn)=@ARGV; # the sequences to be searched, the regex for the motif, length of the motif and the output filename
my ($flag,$name,$seq,$line,$num,$matches);
my $total_motif=0;
my $total_seq=0;
my $seqs=0;
my $regex=qr/$motif/;
my @motifs=();

open(FILE_IN,$seq_file);
open(FILE_OUT,"> ".$fn);

while ($line=<FILE_IN>)
{
#  if ($line=~/peak/) {next;}

  ($name,$seq)=split("\t",$line);
  $seqs++;

  $num=0;
  $flag=0;
  $matches="";

  print FILE_OUT $name."\t";

  while ($seq=~/($regex)/g)
  {
    $flag=1;
    $total_motif++;
    $num++;
    $matches.="_".(pos($seq)-$length);  
    push @motifs,$1;
  }
  
  $total_seq+=$flag;
  print FILE_OUT $matches."\n";
}

close(FILE_IN);
close(FILE_OUT);

print "Input file: ".$seq_file."\n";
print "Total number of input regions is ".$seqs."\n";
print "Total number of motifs found is ".$total_motif." ".(100*$total_motif/$seqs)."%\n";
print "Total number of sequences with motifs found is ".$total_seq." ".(100*$total_seq/$seqs)."%\n";
print "Motif:\n";

no strict "refs";
no warnings;

my (@letters,$i);
my @bases=qw(A T G C);

foreach $i (0..$length-1)
{
  foreach (@bases)
  {
    ${$_}[$i]=0;
  }
}

foreach $motif (@motifs)
{
  @letters=split("",$motif);
  foreach $i (0..$length-1)
  {
    ${$letters[$i]}[$i]++;
  }
}

foreach $i (0..$length-1)
{
  print "base ".($i+1)." ";

  foreach (@bases)
  {
    print $_."(";
    printf("%.3f",${$_}[$i]/($A[$i]+$T[$i]+$G[$i]+$C[$i]));
    print ") ";
  }  
  print "\n";
}

exit;
