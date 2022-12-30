
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/Mathomics
# Copyright (c)
# year: 2017
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
	print "$0 usage : -a  -b  -c\n";
	print "Error in use\n";
	exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  ) {
	usage;
}



open(FILE,$opts{a}) or die "cannot open file $opts{a}\n";
my @svs=();
while(my $line=<FILE>){
	chomp $line;
	if($line=~/^#/){
	  print $line."\n";
	next;
	}
	my @d=split /\t/,$line;
	

	my ($h)=parse_tags($d[7]);

	my $hash=();
	if($h->{SVTYPE} eq "DEL"){
	$hash->{$_}++ foreach (split("",$d[3]));
	#skip homopolymer deletions
	next if(scalar(keys %{$hash})==1);
	my $s=0;
	foreach my $l(keys %{$hash}){
		if($hash->{$l}/length($d[3]) > 0.8){
			$s=1;
		}
	}
	next if($s==1);
	}elsif($h->{SVTYPE} eq "INS"){
	$hash->{$_}++ foreach (split("",$d[4]));
	#skip homopolymer deletions
	next if(scalar(keys %{$hash})==1);
	my $s=0;
	foreach my $l(keys %{$hash}){
		if($hash->{$l}/length($d[4]) > 0.8){
			$s=1;
		}
	}
	next if($s==1);
	}else{
		next;
	}
	print $line."\n";

}

sub parse_tags{
	my ($d)=@_;
	my @data=split(";",$d);
	my $h=();
	foreach my $l (@data){
		my ($a,$b)=split("=",$l);
		$h->{$a}=$b;
	}
	return $h;
}	
