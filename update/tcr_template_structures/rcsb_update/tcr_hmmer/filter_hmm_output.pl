#!/usr/bin/perl

use strict;

my $hmmout = $ARGV[0];
my $cutoff = $ARGV[1];
my $skip_file = $ARGV[2];  # file of TCR delta chains

if (($hmmout eq "") || ($cutoff eq "")) { die("usage: filter_hmm_output.pl hmmout cutoff [skip_file]\n"); }

open(HMM, $hmmout) || die("unable to open file: $hmmout\n");
my @hmm_lines = <HMM>;
close(HMM);
#chomp(@hmm_lines);

my %skip_pdb = ();
if ($skip_file ne "")
{
    open(SKP, $skip_file) || die("unable to open file: $skip_file\n");
    my @skip_lines = <SKP>;
    close(SKP);
    foreach my $line (@skip_lines)
    {
	my @fields = split(" ", $line);
	my $pdb = $fields[8];
	if ($pdb eq "") { next; } # this should probably not happen
	$skip_pdb{$pdb} = 1;
    }
}

my $line_num = 0;
my $line = $hmm_lines[$line_num];
while ((!($line =~ /E-value  score  bias    E-value  score  bias    exp  N  Sequence Description/)) && ($line_num < @hmm_lines)) { $line = $hmm_lines[++$line_num]; }
if ($line_num >= @hmm_lines) { die("error: score start not found!!\n"); }
$line_num += 2;
$line = $hmm_lines[$line_num];
while ((substr($line, 0, 35) ne "  ------ inclusion threshold ------") && ($line_num < @hmm_lines))
{
    my @fields = split(" ", $line);
    if (($fields[0] <= $cutoff) && ($fields[3] <= $cutoff) && ($skip_pdb{$fields[8]} != 1)) { print $line; }
    #if (($fields[0] <= $cutoff) && ($fields[3] <= $cutoff)) { print $line; }
    $line = $hmm_lines[++$line_num];
}
