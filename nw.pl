#!/usr/bin/perl 
#===============================================================================
#
#         FILE: nw.pl
#
#        USAGE: ./nw.pl  
#
#  DESCRIPTION: Needleman-Wunsch algorithm
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Zhang Xiao'ou
#      COMPANY: 
#      VERSION: 1.0
#      CREATED: 12/22/2011 21:56:01
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use Data::Dumper;

my %n=("A"=>0,"G"=>1,"C"=>2,"T"=>3);

print "Please input the similarity matrix file (SEQ:A G C T):";
chomp(my $sf = <STDIN>);
open my $f,"<","$sf" or die "Can't open $sf: $!";
my $sm = [];
while (<$f>) {
    chomp;
    my @line=split;
    push @$sm,\@line;
}
close $f;

RE1: print "Please input the gap penalty:";
chomp(my $d=<STDIN>);
if ($d !~ /^(-)?\d\d?$/) {
    print "The gap penalty must be an integer!!!\n";
    goto RE1;
}

RE2: print "Please input sequence1:";
chomp(my $s1 = <STDIN>);
if ($s1 !~ /^[A|G|C|T][A|G|C|T]*[A|G|C|T]$/){
    print "Your sequence is wrong!!!\n";
    goto RE2;
}
my $lenA=length($s1);

RE3: print "Please input sequence2:";
chomp(my $s2 = <STDIN>);
if ($s2 !~ /^[A|G|C|T][A|G|C|T]*[A|G|C|T]$/){
    print "Your sequence is wrong!!!\n";
    goto RE3;
}
my $lenB=length($s2);

my $fm=[];

$fm->[$_][0] = $_*$d for 0..$lenA;
$fm->[0][$_] = $_*$d for 0..$lenB;
for my $i (1..$lenA) {
    for my $j (1..$lenB) {
        my $a=substr($s1,$i-1,1);
        my $b=substr($s2,$j-1,1);
        my $match=$fm->[$i-1][$j-1] + $sm->[$n{$a}][$n{$b}];
        my $delete=$fm->[$i-1][$j] + $d;
        my $insert=$fm->[$i][$j-1] + $d;
        $fm->[$i][$j]=$match>$insert ? $match : $insert>$delete ? $insert : $delete;
    }
}

my (@a1,@a2);
my ($i,$j)=($lenA,$lenB);
while ($i>0 && $j>0) {
    my $s = $fm->[$i][$j];
    my $sd = $fm->[$i-1][$j-1];
    my $su = $fm->[$i][$j-1];
    my $sl = $fm->[$i-1][$j];
    if ($s==$sd+$sm->[$n{substr($s1,$i-1,1)}][$n{substr($s2,$j-1,1)}]){
        unshift @a1,substr($s1,$i-1,1);
        unshift @a2,substr($s2,$j-1,1);
        $i--;
        $j--;
    }
    elsif ($s==$sl + $d) {
        unshift @a1,substr($s1,$i-1,1);
        unshift @a2,"-";
        $i--;
    }
    else {
        unshift @a1,"-";
        unshift @a2,substr($s2,$j-1,1);
        $j--;
    }
}
while ($i>0) {
    unshift @a1,substr($s1,$i-1,1);
    unshift @a2,"-";
    $i--;
}
while ($j>0) {
    unshift @a1,"-";
    unshift @a2,substr($s2,$j-1,1);
    $j--;
}

print "\t\t".join("\t",split(//,$s2))."\n";
print "\t".join("\t",@{$fm->[0]})."\n";
print substr($s1,$_-1,1)."\t".join("\t",@{$fm->[$_]})."\n" for 1..$lenA;
print "\n";
print "@a1\n";
print "@a2\n";
