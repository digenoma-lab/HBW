
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
	next if($line=~/^#/);
	my @d=split /\t/,$line;
	

	my ($h)=parse_tags($d[7]);

	my $hash=();
	if($h->{SVTYPE} eq "DEL"){
	$hash->{$_}++ foreach (split("",$d[3]));
	#skip homopolymer deletions
	next if(scalar(keys %{$hash})==1);
	my $s=0;
	foreach my $l(keys %{$hash}){
		if($hash->{$l}/length($d[3]) > 0.9){
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
		if($hash->{$l}/length($d[4]) > 0.9){
			$s=1;
		}
	}
	next if($s==1);
	}else{
		next;
	}
	#          'Kurtosis_quant_start' => '-1.170500',
         # 'STRANDS' => '+-',
         # 'END' => '10158',
         # 'RE' => '15',
         # 'SVLEN' => '-39',
          #'SVMETHOD' => 'Snifflesv1.0.11',
          #'PRECISE' => undef,
          #'Kurtosis_quant_stop' => '-1.275763',
          #'SUPTYPE' => 'AL',
          #'SVTYPE' => 'DEL',
         # 'CHR2' => 'WUC105328',
          #'STD_quant_stop' => '3.619392',
          #'STD_quant_start' => '2.449490'

	#print join(" ",abs($h->{SVLEN}),$h->{SVTYPE},$h->{RE})."\n";
	#WUC104949	18514	93	AAGGCCAGGAGGATGTTATTGCATGCAGAGAAAAGAGCATGAGTACATAGCAGGGTCCAGC	N	.	PASS	PRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=WUC104949;END=18575;STD_quant_start=0.000000;STD_quant_stop=0.000000;Kurtosis_quant_start=16.318228;Kurtosis_quant_stop=7.197710;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-61;STRANDS=+-;RE=47	GT:DR:DV	./.:.:47
	next if($h->{RE} < 15 || $h->{RE} > 60);
	if(abs($h->{SVLEN}) >= 30){
	#print $line."\n";
	#print join("\t",$d[0],$d[1],$h->{END},$d[2],$h->{SVLEN},$h->{SVTYPE},$h->{RE},$d[3],$d[4])."\n"
	my $tmp=();
	$tmp->{ctg}=$d[0];
	$tmp->{start}=$d[1];
	$tmp->{stop}=$h->{END};
	$tmp->{ref}=$d[3];
	$tmp->{alt}=$d[4];
	$tmp->{id}=$d[2];
	$tmp->{tags}=$h;
	push(@svs,$tmp);
	}
}

my @s_svs=sort {$a->{ctg} cmp $b->{ctg} || $a->{start} <=> $b->{start} || $a->{stop} <=> $b->{stop}} @svs;
#we build the cluster of SVs
my $cluster=0;
$s_svs[0]->{cluster}=$cluster;	
my $d=5000;
my $c_count=();
$c_count->{$cluster}->{count}++;
push(@{$c_count->{$cluster}->{index}},0);
for(my $i=1 ; $i <= $#s_svs; $i++){
	if($s_svs[$i]->{start}-$s_svs[$i-1]->{start} < $d and $s_svs[$i]->{ctg} eq $s_svs[$i-1]->{ctg}){
	}else{
	$cluster++;
	}
	$s_svs[$i]->{cluster}=$cluster;
	$c_count->{$cluster}->{count}++;
	push(@{$c_count->{$cluster}->{index}},$i);
}
# we select the best SV of each cluster, the one with higher read support
my $selected_sv=();
foreach my $c(keys %{$c_count}){
	if($c_count->{$c}->{count} == 1){
		$selected_sv->{$s_svs[$_]->{id}}=1 foreach(@{$c_count->{$c}->{index}});
	}elsif($c_count->{$c}->{count} < 5){
		my $mx=0;
		my $id=0;
		foreach(@{$c_count->{$c}->{index}}){
			if($s_svs[$_]->{tags}->{RE} > $mx){
				$mx=$s_svs[$_]->{tags}->{RE};
				$id=$s_svs[$_]->{id};
			}
		}
		$selected_sv->{$id}=1;
	}
}

#we print only sv with a maximum of 5 clusters
foreach my $sv(@s_svs){
	#if($c_count->{$sv->{cluster}}->{count} < 4){
	#	print join("\t",$sv->{ctg},$sv->{start},$sv->{stop},$sv->{id},$sv->{tags}->{SVTYPE},$sv->{tags}->{RE},$sv->{tags}->{SVLEN}, $sv->{cluster}, defined $selected_sv->{$sv->{id}} ? 1:0)."\n";
	#}
	if(defined $selected_sv->{$sv->{id}}){
		print join("\t",$sv->{ctg},$sv->{start},$sv->{stop},$sv->{id},$sv->{tags}->{SVTYPE},$sv->{tags}->{RE},$sv->{tags}->{SVLEN}, $sv->{cluster},$sv->{ref}, $sv->{alt},  defined $selected_sv->{$sv->{id}} ? 1:0)."\n";
	}
	
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
