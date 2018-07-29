#!/usr/bin/perl

use strict;
use File::Basename;
my $verbose=0;
my $dir = dirname(__FILE__);
my $aachangefile=$ARGV[0];
my $snapdir=$dir."/currentyeast";

#EXONIC:YAL069W:YAL069W:exon1:c.G3A:p.M1I,nonsynonymous SNV
#
# input is file with list of ANNOVAR annotations like the one above
#

print "#GENE\tAACHANGE\tSNAPSCORE\tSNAPSCORE_znorm\tSNAPSCORE_scaled\n";

open(AA,$aachangefile) or die $!;
while(my $str=<AA>){
    chomp($str);
    next if $str!~/\S/;
    my ($exonic,$yad,$yad2,$exon,$cdna,$protein)=split(/\:/,$str);
    if($protein =~ /nonsyn/){
	$protein=~s/,nonsynonymous SNV//;
	my $yadfile=$snapdir."/$yad.snap.tsv";
	print STDERR "# FILE $yadfile\n# PROTEIN: $protein" if $verbose;
	if(-e $yadfile){
	    my ($score,$scorenorm)=getSNAPscore($protein,$yadfile);
	    my $scorescaled=sprintf("%.2f",($score+100)/200);
	    print "$yad\t$protein\t$score\t$scorenorm\t$scorescaled\n";
	    print STDERR "# FILE exists\n" if $verbose;
	}else{
	    print STDERR "# FILE $yadfile does not exist\n" if $verbose;
	    print "$yad\t$protein\tNA\tNA\tNA\n";
	}
    }else{
	print STDERR "# NOT SYNONYMOUS: $str\n";
    }
}
close(AA);

sub getSNAPscore{
    my $aachange= $_[0];
    my $file    = $_[1];
    my $score=0;
    my $scorenorm=0;

    $aachange=~/p.([A-Z])([0-9]+)([A-Z])/;
    my $ref=$1;
    my $pos=$2;
    my $alt=$3;

    my %aas=("A"=>1,"C"=>2,"D"=>3,"E"=>4,"F"=>5,"G"=>6,"H"=>7,"I"=>8,"K"=>9,"L"=>10,"M"=>11,"N"=>12,"P"=>13,"Q"=>14,"R"=>15,"S"=>16,"T"=>17,"V"=>18,"W"=>19,"X"=>20,"Y"=>21);

    my $ok=0;
    open(FILE,$file) or die "# ERROR: cannot read $file\n";
    while(my $str=<FILE>){
	chomp($str);
	next if $str!~/\S/;
	next if $str=~/^#/;

	my @split=split(/\t/,$str);
	if ($split[0] eq "$ref$pos:"){
	    $ok=1;
	    die "# ERROR: $alt not known !!\n" if !exists($aas{$alt});
	    print STDERR "# REF=$ref ALT= $alt $aas{$alt} $split[$aas{$alt}]\n" if $verbose; 
	    $score=$split[$aas{$alt}];
	    #$scorenorm=sprintf("%.2f",$score/$split[22]);
	    my @array-0;
	    for(my $i=1;$i<=21;$i++){
		next if $i eq 20;
		push @array,$split[$i];
	    }
	    my ($mean,$std)=getMeanStd(\@array);
	    $scorenorm=sprintf("%.2f",($score-$mean)/$std);
	}

	last if $ok
    }
    close(FILE);

    print STDERR "# AACHANGE: $aachange ref=$ref alt=$alt pos=$pos snapscore=$score snapscorenorm=$scorenorm\n" if $verbose;

    return ($score,$scorenorm);
}

# getMeanStd calculates the mean and the standard deviation 
# # for given array
#
#
sub getMeanStd{

    my @values = @{$_[0]};
    my $n      = scalar(@values);

    my $mean = 0;
    my $std  = 0;

    return (0,0) if $n==0;

    foreach(@values){
	$mean+=$_;
    }
    $mean/=$n;
    
    foreach(@values){
	$std+=($_-$mean)*($_-$mean);
    }
    if ($n>1) {
	$std /= ($n-1);
	$std  = sqrt($std);
    }else{
	$std=0;
    }

    return ($mean,$std);
}

