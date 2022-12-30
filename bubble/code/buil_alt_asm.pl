
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

my ($int)=load_intervals($opts{a});
#print Dumper($int);
open(FASTA,$opts{b}) or die "cannot open file $opts{b}\n"; 
my @aux = undef;
my ($name, $seq, $qual);
while (($name, $seq, $qual) = readfq(\*FASTA, \@aux)) {

	if(defined $int->{$name}){
		my @s=split //,$seq;
		foreach my $b (@{$int->{$name}}){
			if($b->{t} eq "DEL"){
				for(my $i=$b->{s}; $i < $b->{e}; $i++){
					$s[$i]="N";
				}
			}else{
				#we add and insertion
				$s[$b->{s}]=$b->{alt};
			}
		
		}
		$seq=join("",@s);
		$seq=~s/N//g;
	}
	print ">".$name."\n".$seq."\n";
}



sub readfq {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if ($aux == undef);
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}





sub load_intervals{
	my ($f)=@_;
	open(FILE,$opts{a}) or die "cannot open file $opts{f}\n";
	my $h=();
	while(my $line=<FILE>){
	chomp $line;
	my @data=split /\t/,$line;
	my $tmp=();
	#WSC100042       169878  170530  16111   DEL     32      -652    2       AATTTCCCCTGAGTAACTCAGGGTGAATACCGGGAATTGAGAATGCTGTGTTCAGAGGGAT

	$tmp->{ctg}=$data[0];
	$tmp->{s}=$data[1];
	$tmp->{e}=$data[2];
	$tmp->{t}=$data[4];
	$tmp->{ref}=$data[8];
	$tmp->{alt}=$data[9];
	
	push(@{$h->{$tmp->{ctg}}},$tmp);
	}
	return $h;
}



