use strict;
use warnings;
use Data::Dumper;

my($dir,$head,$out)=@ARGV;

open IN,"$head" or die "$!";
my $h;
while(<IN>){
	chomp;
	$h = "Sample\t" . $_;
}
close IN;

my @files=glob "$dir/*.txt";

open OUT,">$out" or die "$!";
print OUT $h,"\n";
my $empty = "NA\t"x170;
$empty =~ s/\t$//;
for my $i(@files){
	my $mark=`basename $i`;
	$mark=~s/.TNscope.*txt//;
	chomp $mark;
	open IN,"$i" or die "$!";
	my $n = 0;
	<IN>;
	while(<IN>){
		chomp;
		$n++;
		print OUT $mark,"\t",$_,"\n";
	}
	close IN;
	if ($n == 0){
		print OUT $mark,"\t",$empty,"\n";
	}
	$n = 0;
}
close OUT;
