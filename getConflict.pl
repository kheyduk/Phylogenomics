#!/usr/bin/perl
#written by Karolina Heyduk, UGA kheyduk@plantbio.uga.edu. 2014

use strict;
#use 5.010;
use Bio::TreeIO;
use Bio::Tree::Tree;
use IO::String;
use Array::Utils qw(:all);
use Data::Dumper;
my $genefile = $ARGV[0]; #file containing all gene trees, newick format
my $speciestree = $ARGV[1]; #Species tree

my @trees;
my %species;
my %genes;
my %Gtree;
my %sptips;

#put all gene trees into array, also count number of tips per tree
my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $genefile);
while (my $tree = $treeio->next_tree) {
	my @nodes;
	my @tipids;
	push(@trees, $tree);
	push(@nodes, $tree->get_nodes);
	foreach my $node (@nodes) {
		if ($node->is_Leaf) {
			my $id = $node->id;
			push (@tipids, $id);
			}
		else {
			next;
			}
		}
	#my $length = @tipids;
	#print "$length\t@tipids\t";
	my @sortedtaxa = sort @tipids;
	#my $length2 = @sortedtaxa;
	$Gtree{$tree} = \@sortedtaxa;
	@nodes = ();
	}

#print Dumper %SPtree;
my @SPnodes;
my $SPtree;
#get node and tip information from species tree
my $speciesio = Bio::TreeIO->new(-format => 'newick', -file => $speciestree); 
while (my $species = $speciesio->next_tree) {
		$SPtree = $species;
		push (@SPnodes, ($species->get_nodes));
		foreach my $SPnode (@SPnodes) {
	
			if ($SPnode->is_Leaf) {
				next;
				}
			else {
				my @SPtaxas;
				my @SPtips;
				my @rawids;
				push(@SPtips, ($SPnode->get_all_Descendents));
				foreach my $SPtip (@SPtips) {
					if ($SPtip->is_Leaf) {
						push (@rawids, $SPtip);
						my $SPtaxa = $SPtip->id;
						push(@SPtaxas, $SPtaxa);
						}
					else {
						next;
						}
					}
				my @rawsort = sort @rawids;
				$sptips{$SPnode} = \@rawsort;
				my @sorted = sort @SPtaxas;
				$species{$SPnode} = \@sorted;
				}
			}
		}

		
foreach my $tree (@trees) { #for every gene tree, get descendants of each node
	$tree->move_id_to_bootstrap;
	my @nodes;
	push (@nodes, ($tree->get_nodes)); 
	foreach my $node (@nodes) {
		if ($node->is_Leaf) {
			next;
			}
		else {
			my @IDs;
			push (@IDs, $node->get_all_Descendents);
			my @ids2;
			my $BS;
			my $BS = $node->bootstrap;
			foreach my $ID (@IDs) {
				if ($ID->is_Leaf) {
					my $taxa = $ID->id;
					push(@ids2, $taxa);
					}
				else {
					next;
					}
				}
			my @sorted = sort @ids2;
			#print "@sorted\n";
			if (length($BS)<1) {
				next;
				}
			else {
				$genes{$tree}{$node}{$BS} = \@sorted;
				}
			}
			
		}
	}
my %mono;
my %poly;
my %polymiss;
my %agreemiss;
my %missing;

foreach my $SPnode (sort keys %species) {
my $speciestips = join('',@{$species{$SPnode}});
	foreach my $tree (@trees) {
	#foreach my $tree (sort keys %genes) {
		#foreach my $node (sort keys %{$genes{$tree}}) {
			#foreach my $BS (sort {$a<=>$b} keys %{$genes{$tree}{$node}}) {	
				#my $genetips = join('', @{$genes{$tree}{$node}{$BS}});
					my @int = intersect(@{$Gtree{$tree}}, @{$species{$SPnode}});
					my $intlength = @int;
					my $SPnodelength = @{$species{$SPnode}};
					if ($SPnodelength == @int) { #all taxa present
						#print "$speciestips\t$SPnodelength\t@int\n";
						my @LCAdes;
						my @LCAleaf;
						my @LCAids;
						my @raw;
						foreach my $int (@int) {
							my $rawid = $tree->find_node(-id=> $int);
							push (@raw, $rawid);
							}
						my $LCA = $tree->get_lca(-nodes=>\@raw);
						my $lcaBS = $LCA->bootstrap;
						push (@LCAdes, $LCA->get_all_Descendents);
						foreach my $LCAdes (@LCAdes) {
							if ($LCAdes->is_Leaf) {
								push (@LCAleaf, $LCAdes);
								my $id = $LCAdes->id;
								push (@LCAids, $id);
								}
							else {
								next;
								}
							}
						if (@LCAleaf == $SPnodelength) { #check monophyly 
							if ($lcaBS >= 80) {
									my $level = 80;
									$mono{$speciestips}{$level}{$tree} = $lcaBS;
									}
								elsif (($lcaBS < 80)&&($lcaBS >= 50)) {
									my $level = 50;
									$mono{$speciestips}{$level}{$tree}=$lcaBS
									}
								elsif (($lcaBS < 50)&&($lcaBS >=20)) {
									my $level = 49;									
									$mono{$speciestips}{$level}{$tree}=$lcaBS;
									}
								elsif ($lcaBS < 20) {
									my $level = 20;
									$mono{$speciestips}{$level}{$tree} = $lcaBS;
									}
							}
						elsif (@LCAleaf > $SPnodelength) { #polyphyletic
							my $level = 0;
							#if ($lcaBS > 50) { #**COMMENT BACK IN FOR SUPPORTED CONFLICT. Change to whatever Bs value you want
							$poly{$speciestips}{$level}{$tree} = $lcaBS;
							#	} #**
							#else { #**
							#	next; #**
							#	} #**
							}
						}
					else { #missing taxa
						my @hashint;
						my $workingtree;
						my @LCAnodes;
						my @LCAtips;
						if (@int > 1) {
						foreach my $int (@int) {
							#print "@int\n";
							my $id = $tree->find_node(-id=> $int);
							push (@hashint, $id);
								}
							my $LCA = $tree->get_lca(-nodes => \@hashint); 
							my $lcaBS = $LCA->bootstrap; 
							push (@LCAnodes, $LCA->get_all_Descendents);
							#print "@LCAnodes\n";
							foreach my $LCAnode (@LCAnodes) {
								if ($LCAnode->is_Leaf) {
									push(@LCAtips, $LCAnode);
										}
									}
							if (@LCAtips == @int) { #monophyly with missing
								if ($lcaBS >= 80) {
									my $level = 80;
									$agreemiss{$speciestips}{$level}{$tree} = $lcaBS;
									}
								elsif (($lcaBS < 80)&&($lcaBS >= 50)) {
									my $level = 50;
									$agreemiss{$speciestips}{$level}{$tree}=$lcaBS
									}
								elsif (($lcaBS < 50)&&($lcaBS > 20)) {
									my $level = 49;									
									$agreemiss{$speciestips}{$level}{$tree}=$lcaBS
									}
								elsif ($lcaBS < 20) {
									my $level = 20;									
									$agreemiss{$speciestips}{$level}{$tree}=$lcaBS
									}
								}
							elsif (@LCAtips > @int) {
								my $level = 0;
								#if ($lcaBS > 50) { #**COMMENT BACK IN FOR SUPPORTED CONFLICT. Change to whatever Bs value you want
								$polymiss{$speciestips}{$level}{$tree} = $lcaBS;
								#	} #**
							#else { #**
							#	next; #**
							#	} #**
								}
							}
						
					elsif (@int = 1) {
						#print "@int\n";
						$missing{$speciestips}{$tree} = 0;
					}	
					}
				}
			}

print "Agreements:\n";	
foreach my $speciestips (sort keys %mono) {
	foreach my $level (sort keys %{$mono{$speciestips}}) {
		my $count = 0;
		my $BSsum;
		foreach my $tree (sort keys %{$mono{$speciestips}{$level}}) {
			$BSsum+=$mono{$speciestips}{$level}{$tree};
			$count++;
			}
	my $avgBS = ($BSsum/$count);
	print "$speciestips\t$level\t$count\t$avgBS\n";
	}
}

#print output - conflict (no missing)
print "Conflict - no missing:\n";
foreach my $speciestips (sort keys %poly) {
	foreach my $level (sort keys %{$poly{$speciestips}}) {
		my $count = 0;
		my $BSsum;
		foreach my $tree (sort keys %{$poly{$speciestips}{$level}}) {
				$BSsum+=$poly{$speciestips}{$level}{$tree};
				$count++;
			}
		my $avgBS = ($BSsum/$count);
		print "$speciestips\t$level\t$count\t$avgBS\n";
	}
}

print "Agreement - missing:\n";
foreach my $speciestips (sort keys %agreemiss) {
	foreach my $level (sort keys %{$agreemiss{$speciestips}}) {
			my $count = 0;
			my $BSsum = 0;
			foreach my $tree (sort keys %{$agreemiss{$speciestips}{$level}}) {
				$count++;
				$BSsum+=$agreemiss{$speciestips}{$level}{$tree}
				}
			my $avgBS = ($BSsum/$count);
			print "$speciestips\t$level\t$count\t$avgBS\n";
	}
}

print "Conflict - missing:\n";
#print Dumper %conmissing;

foreach my $speciestips (sort keys %polymiss) {
	foreach my $level (sort keys %{$polymiss{$speciestips}}) {
		my $count = 0;
		my $BSsum = 0;
		foreach my $tree (sort keys %{$polymiss{$speciestips}{$level}}) {
			$count++;
			$BSsum+=$polymiss{$speciestips}{$level}{$tree};
			}
		my $avgBS = ($BSsum/$count);
		print "$speciestips\t$level\t$count\t$avgBS\n";
		}
	}
	
print "Missing data:\n";
foreach my $speciestips (sort keys %missing) {
	my $count = 0;
	foreach my $tree (sort keys %{$missing{$speciestips}}) {
		$count++;
		}
	print "$speciestips\t$count\n";
	}						
			