#!/usr/bin/perl 
#===============================================================================
#
#         FILE: hmm.pl
#
#        USAGE: ./hmm.pl  
#
#  DESCRIPTION: Hidden Markov Model
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zhang Xiao'ou (kepbod)
#      COMPANY: PICB
#      VERSION: 1.0
#      CREATED: 12/27/2011 15:17:03
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

print "\nPlease input all the possible internal states:";
chomp(my $line=<STDIN>);
my @u = split /\s/,$line;
print "There are ".($#u+1)." internal states: @u\n\n";
my %u = map {$u[$_],$_} 0..$#u;

print "Please input all the possible emissions:";
chomp($line=<STDIN>);
my @x = split /\s/,$line;
print "There are ".($#x+1)." emissions: @x\n\n";
my %x = map {$x[$_],$_} 0..$#x;

RE1: print "Please input pi vector:";
chomp($line=<STDIN>);
my @pro=split /\s/,$line;
my $psum=0;
$psum+=$_ for @pro;
if ($#u!=$#pro || $psum ne "1"){
    print "Your input wrong!!!\n";
    goto RE1;
}
my %pi=map {$u[$_],$pro[$_]} 0..$#u;
print "The initial state probabilities: ";
print "$pi{$_}($_) " for @u;
print "\n\n";

my $a=[];
R1: print "Please input the A matrix(end with \"exit\"):\n";
while (<STDIN>) {
    chomp;
    if ($_ eq "exit"){
        last;
    }
    my @line = split;
    my $asum=0;
    $asum+=$_ for @line;
    if ($asum ne "1" || $#line!=$#u){
        print "Your input wrong!!!\n";
        $a=[];
        goto R1;
    }
    push @$a,\@line;
}
if ($#{$a} != $#u) {
    print "Your input wrong!!!\n";
    $a=[];
    goto R1;
}
print "The A matrix is:\n";
print "@{$a->[$_]}\n" for 0..$#{$a};
print "\n";

my $b=[];
R2: print "Please input the B matrix file(end with \"exit\"):\n";
while (<STDIN>) {
    chomp;
    if ($_ eq "exit") {
        last;
    }
    my @line = split;
    my $bsum=0;
    $bsum+=$_ for @line;
    if ($bsum ne "1" || $#line!=$#x){
        print "Your input wrong!!!\n";
        $b=[];
        goto R2;
    }
    push @$b,\@line;
}
if ($#{$b} != $#u) {
    print "Your input wrong!!!\n";
    $b=[];
    goto R2;
}
print "The B matrix is:\n";
print "@{$b->[$_]}\n" for 0..$#{$b};
print "\n";

RE2: print "\nPlease input the sequence:";
chomp($line=<STDIN>);
my @s=split /\s/,$line;
for (@s) {
    if (!defined($x{$_})){
        print "Your input wrong!!!";
        goto RE2;
    }
}

RE3: print "\nPlease choose 1(Forward Algorithm) 2(Backward Algorithm) 3(Viterbi Algorithm) 4(Posterior Decoding Algorithm) 5(EXIT):";
chomp(my $choose=<STDIN>);
if ($choose==1) {
    &fa(1);
}
elsif ($choose==2) {
    &ba(1);
}
elsif ($choose==3) {
    &va;
}
elsif ($choose==4) {
    &pda;
}
elsif ($choose==5) {
    exit;
}
else{
    print "Your input wrong!!!\n";
    goto RE3;
}

sub fa{
    my $fm=[];
    my @ini;
    push @ini,$pi{$_}*$b->[$u{$_}][$x{$s[0]}] for (@u);
    push @$fm,\@ini;
    for my $e (1..$#s) {
        my @step;
        for my $i (0..$#u){
            my $sum=0;
            for (0..$#u){
                $sum+=$fm->[$e-1][$_]*$a->[$_][$i];
            }
            $sum*=$b->[$i][$x{$s[$e]}];
            push @step,$sum;
        }
        push @$fm,\@step;
    }
    if ($_[0]==1){
        print "\nForward Algorithm:\n";
        print "@{$fm->[$_]}\n" for 0..$#{$fm};
        my $pro;
        $pro+=$fm->[-1][$_] for (0..$#u);
        print "The probability: $pro\n";
        print "\n";
        goto RE3;
    }
    else {
        return $fm;
    }
}

sub ba{
    my $bm=[];
    my @ini;
    push @ini,1 for (@u);
    push @$bm,\@ini;
    for my $e (reverse(0..$#s-1)) {
        my @step;
        for my $i (0..$#u) {
            my $sum=0;
            for (0..$#u) {
                $sum+=$bm->[0][$_]*$a->[$i][$_]*$b->[$_][$x{$s[$e+1]}];
            }
            push @step,$sum;
        }
        unshift @$bm,\@step;
    }
    if ($_[0]==1) {
        print "\nBackward Algorithm:\n";
        print "@{$bm->[$_]}\n" for 0..$#{$bm};
        my $pro;
        $pro+=$bm->[0][$_]*$pi{$u[$_]}*$b->[$_][$x{$s[0]}] for (0..$#u);
        print "The probability: $pro\n";
        print "\n";
        goto RE3;
    }
    else {
        return $bm;
    }
}
sub va{
    my $vm1=[];
    my $vm2=[];
    my @ini;
    push @ini,$pi{$_}*$b->[$u{$_}][$x{$s[0]}] for (@u);
    push @$vm1,\@ini;
    for my $e (1..$#s) {
        my (@step1,@step2);
        for my $i (0..$#u){
            my $max=0;
            my $sign;
            for (0..$#u) {
                my $tem=$vm1->[$e-1][$_]*$a->[$_][$i];
                if ($tem > $max) {
                    $max=$tem;
                    $sign=$_;
                }
            }
            push @step1,$max*$b->[$i][$x{$s[$e]}];
            push @step2,$sign;
        }
        push @$vm1,\@step1;
        push @$vm2,\@step2;
    }
    my $lmax=0;
    my $lsign;
    for (0..$#{$vm1->[-1]}) {
        if ($vm1->[-1][$_] > $lmax) {
            $lmax=$vm1->[-1][$_];
            $lsign=$_;
        }
    }
    print "\nViterbi Algorith:\n";
    print "@{$vm1->[$_]}\n" for 0..$#{$vm1};
    my @ls;
    push @ls,$lsign;
    unshift @ls,$vm2->[$_][$ls[0]] for reverse(0..$#{$vm2});
    print "The internal sequence:";
    print $u[$_]." " for @ls;
    print "\n";
    goto RE3;
}

sub pda{
    my $fm=&fa(0);
    my $bm=&ba(0);
    my @ls;
    my $pdm=[];
    for my $e (0..$#s) {
        my @step;
        my $max=0;
        my $sign;
        for my $i (0..$#u){
            my $tem=$fm->[$e][$i]*$bm->[$e][$i];
            push @step,$tem;
            if ($tem > $max) {
                $max=$tem;
                $sign=$i
            }
        }
        push @$pdm,\@step;
        push @ls,$sign;
    }
    print "\nPosterior Decoding Algorithm:\n";
    print "@{$pdm->[$_]}\n" for 0..$#{$pdm};
    print "The internal sequence:";
    print $u[$_]." " for @ls;
    print "\n";
    goto RE3;
}
