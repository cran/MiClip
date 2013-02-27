#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/floor/;

my ($suffix_f,$suffix_r,$file_tag,$file_merge)=@ARGV;
my ($i,%tags,@items,$tag,$name);
my $md_field=-1;
$suffix_f='^(.*)'.$suffix_f.'$';
$suffix_r='^(.*)'.$suffix_r.'$';
$suffix_f=qr/$suffix_f/;
$suffix_r=qr/$suffix_r/;

open(FILE_IN,$file_tag) or die "read alignment file ".$file_tag." failed!\n";

while (<FILE_IN>) # search for MD field
{
  if ($_=~/^@/) {next;}

  @items=split("\t",$_);
  $i=0;

  foreach (@items)
  {
    if ($_=~/MD\:Z\:/) {$md_field=$i;}
    $i++;
  }

  if ($md_field>0) {last;}
}

if ($md_field==-1) {die "no MD field found in SAM file!\nMaybe your file is not paired-end?\n";}
seek FILE_IN,0,0;

while (<FILE_IN>) # find mates
{
  if ($_=~/^@/) {next;}
  @items=split("\t",$_);
  if ($items[5]=~/N/) {next;} # delete gapped mapping
  unless ($items[1]==147 || $items[1]==163 || $items[1]==83 || $items[1]==99) {next;}
  
  if (defined $items[$md_field] && $items[$md_field]=~/MD/)
  {
    $items[$md_field]=~s/MD\:Z\://;
  }else
  {
    $items[$md_field]=-1;
  }

  # information is stored as pos, CIGAR, seq and MD

  if ($items[1]==83)
  {
    $items[0]=~/$suffix_f/;
    $name=$1."_".($items[3]+$items[7]);
    $tags{$name}->{"strand"}="-";
    $tags{$name}->{"chr"}=$items[2];
    $tags{$name}->{"f"}=[$items[3],$items[5],$items[9],$items[$md_field]];
  }elsif ($items[1]==99)
  {
    $items[0]=~/$suffix_f/;
    $name=$1."_".($items[3]+$items[7]);
    $tags{$name}->{"strand"}="+";
    $tags{$name}->{"chr"}=$items[2];
    $tags{$name}->{"f"}=[$items[3],$items[5],$items[9],$items[$md_field]];
  }elsif ($items[1]==147)
  {
    $items[0]=~/$suffix_r/;
    $name=$1."_".($items[3]+$items[7]);
    $tags{$name}->{"r"}=[$items[3],$items[5],$items[9],$items[$md_field]];
  }else
  {
    $items[0]=~/$suffix_r/;
    $name=$1."_".($items[3]+$items[7]);
    $tags{$name}->{"r"}=[$items[3],$items[5],$items[9],$items[$md_field]];
  }
}

close(FILE_IN);

open(FILE_OUT,">".$file_merge) or die "can't write to file ".$file_merge."!\n";

foreach $tag (keys %tags) # write mates into file
{
  if (defined $tags{$tag}->{"f"} && $tags{$tag}->{"r"} && abs($tags{$tag}->{"f"}->[0]-$tags{$tag}->{"r"}->[0])<500)
  {
    # tag name, chr, strand, f_pos, f_CIGAR, f_seq, f_MD, r_pos, r_CIGAR, r_seq and r_MD
    print FILE_OUT $tag."\t".$tags{$tag}->{"chr"}."\t".$tags{$tag}->{"strand"}."\t".join("\t",@{$tags{$tag}->{"f"}})."\t".join("\t",@{$tags{$tag}->{"r"}})."\n";
  }
}

close(FILE_OUT);

exit;

#10010011 R->n 147
#10100011 R->p 163
#01010011 F->n 83
#01100011 F->p 99
