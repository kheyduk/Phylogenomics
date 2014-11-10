#!/usr/bin/perl
use strict;
use Bio::SeqIO;
#script to extract both informative and non-informative SNPs and their position from alignments
#usage: perl getSNPs.pl .ending
#run in directory where fasta alignements are. Note that this script assumes fasta format - if different, change line X accordingly
#output: SNPS.positions.txt lists every position in every alignment that is either informative (at least two different SNPs each in at least 
#two taxa) or simply has a SNP in it, SNPS.counts.txt lists informative and polymorphic sites per file.

my $ending = $ARGV[0]; #please give an ending that all your alignment files have in commong (ie, .fasta).
my @files = glob("*$ending"); 

my %PI;
my %poly;
my %lengths;

foreach my $file (@files) {
	my $length;
	my @headers;
	my %hash = ();
	my $fasta = Bio::SeqIO->new(-format => 'fasta', -file => $file);
		while (my $io_obj = $fasta -> next_seq() ) {
			my $header = $io_obj->id();
			push (@headers, $header);
			my $seq = $io_obj->seq();
			$length = $io_obj->length();
			my @bases = split("", $seq);
			$hash{$header} = \@bases;
			$lengths{$file} = $length;
			}
		
		my @pos = (0..$length);
		foreach my $count (@pos) { #iterate over each base until you reach the end of alignment
			my @position = ();
			#print "$count\t$length\n";
			#print "$count\n";
			foreach my $header (@headers) {
				my $base = @{$hash{$header}}[$count];
				if ($base =~ "-") { #ignore bases which have missing data
					next;
					}
				else {
					#print "$base\t";
					push (@position, $base); #put all bases from all taxa at that position into array
					}
				}
			#print "@position\n";
			
			my %counts = ();
			my @unique = grep !$counts{$_}++, @position; #how many unique elements?
			my $bases = @unique;
			if ($bases == 1) { 
				next; #monomorphic
				}
			else {
				my $counter = 0;
				my $countSNP = 0;
				foreach my $unique (@unique) {
					my @greps = grep (/$unique/, @position);
					#print "$count\t@greps\n";
					my $countSNP = @greps;
					#print "$count\t$countSNP\n";
					if ($countSNP >= 2) {
						$counter++;
						}
					else {
						next;
						}
					}
				#	print "$count\t$counter\n";
					if ($counter >= 2) {
						$PI{$file}{($count+1)} = 0;
						}
					else {
						$poly{$file}{($count+1)} = 0;
						}
					
				}
			}
		}
		
my %sizePI;
my %sizepoly;
open POS, ">>SNPS.positions.txt";
#print out polyinformative sites:
print POS "parsimony informative sites:\n";
my $informcounter = 0;
foreach my $file (sort keys %PI) {
	foreach my $position (sort {$a<=>$b} keys %{$PI{$file}}) {
		print POS "$file\t$position\n";
		my $values = keys %{$PI{$file}};
		$sizePI{$file} = $values;
		}
	}
	
#print poly but non informative:
print POS "polymorphic, noninformative sites:\n";
foreach my $file (sort keys %poly) {
	foreach my $position (sort {$a<=>$b} keys %{$poly{$file}}) {
		print POS "$file\t$position\n";
		my $values = keys %{$poly{$file}};
		$sizepoly{$file} = $values;
		}
	}
close POS;

open OUT, ">>SNPS.counts.txt";
print OUT "fileID\tinformative\tpoly noninformative\ttotal length\n";
foreach my $file (@files) {
	print OUT "$file\t$sizePI{$file}\t$sizepoly{$file}\t$lengths{$file}\n";
	}
close OUT;

	
		