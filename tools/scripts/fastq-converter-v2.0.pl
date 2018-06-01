#!/usr/bin/perl
#
#  File: fastq-converter.pl
#  Created: Dec 1st, 2009
#
#  Author: Xiao Yang <isuxyang@gmail.com>
#
#  Copyright 2009 Xiao Yang
#
#  fastq-converter.pl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  fastq-converter.pl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#  You should have received a copy of the GNU Lesser General Public License
#  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
#
 

# converting each .fastq file in directory [InDIR] into 
# .fa fasta file as well as corresponding raw .q quality file (added @ in front of each entry)
# and converted .q quality file with (dec)
# When setting [Flag] = 0, all "N"s will be converted to "A" if for any window with size w that
# contains this "N", there exists in total <= d number of "N"s

use List::Util qw(min max);
use Data::Dumper;

my $numArgs = $#ARGV + 1;

if (($numArgs != 3) && ($numArgs != 6)) {
	print "\n--------------------------------------------------------------------\n";
	print "usage: ./fastaq-converter.pl [InDIR] [OutDIR] [Flag] [w] [d] [char]\n\n";
	print "Flag --- 1) Flag = 1: keep reads containing ACGT(acgt) only\n";
	print "\t 2) Flag = 2: keep reads intact\n";
	print "\t 3) Flag = 0: convert ambiguous characters to [char] (ACGTacgt);\n\t need to specify optional paramters (w,d,char)\n";
	print "w : window size\n";
	print "d : maximum number of ambiguous characters allowed in length w window\n";
	print "All files in [InDIR] with extension .fastq will be processed\n";
	print "--------------------------------------------------------------------\n\n";
	exit;
}

my $dir = $ARGV[0];
my $odir = $ARGV[1];
my $flag = $ARGV[2];
if ($flag != 1 && $flag != 0 && $flag != 2) {
	print "please enter correct flag value (0 or 1 or 2)\n"; exit;
}
my $w = 0; #window size
my $d = 0;
if ($flag == 0) {
	$w = $ARGV[3];
	$d = $ARGV[4];
	$ch = $ARGV[5];
}

opendir(ODIR, $odir) or die "couldn't open $odir: $!\n";
opendir(DIR, $dir) or die "couldn't open $dir: $!\n";
my @files = readdir DIR;
closedir DIR;

my $readRank = 0;
my $discarded = 0;
my $fileNum = 0;
foreach $file (@files) {

	if ($file =~ /fastq$/) {	

		print "\nprocess file: $fileNum: $file\n\n";
		$fileNum ++;		
		my $inpath = $dir.$file;
		my @fileName = split(/\./, $file);
		my $outFa = $odir.$fileName[0].".fa";
		my $outRawQual = $odir.$fileName[0].".q";
        #my $outQual = $odir.$fileName[0].".q";


		open (FILE, "<$inpath") or die "unable to open file $inpath to read\n";
		open (OFa, ">$outFa") or die "can't open file $outFa to write\n";
		open (ORQual, ">$outRawQual") or die "can't open file $outRawQual to write\n";
		#open (OQual, ">$outQual") or die "can't open file $outQual to write\n";

		my $ambig_changed_cnt = 0;
		while (<FILE>){
			my @line = split(/\s+/);
			my $readName = $line[0];
			$readName =~ s/\@//;
			#fa file
			my $read = <FILE>;
			my $skipline = <FILE>;			
			#qual file
			my $qual = <FILE>;			

			chomp($read);

			if ($flag == 1){ # discarding reads containing ambiguous characters				
				 if ($read =~ /[^ACGTacgt]/){
					++ $discarded;
					next;
				}
			}
			
			if ($flag == 0) { # convert ambig chars to A satisfying (w,d) constraint
				
				my @indices;
				my @tochange;
				#print $read."\n";
				while ($read =~ /[^ACGTacgt]/g){ push @indices, (pos($read) - 1);}

				for (my $i = 0; $i < scalar @indices; ++ $i){
					my $found = 0;
					for (my $j = max(0, $i - $d); $j <= $i; ++ $j){
						if ($j + $d >= scalar @indices){
							last;
					    } 
						if ($indices[$j + $d] - $indices[$j] <= $w){
							$found = 1;
							last;
						}
					}
					if ($found == 0){
						push @tochange, $indices[$i];
					}
				}

				$ambig_changed_cnt += scalar @tochange;

				for (my $i = 0; $i < scalar @tochange; ++ $i){
					#substr($read, $tochange[$i], 1) = 'a';		
					substr($read, $tochange[$i], 1) = $ch ;
				}											
			}

			print OFa ">$readRank"."\n".$read."\n";
			print ORQual ">$readRank"."\n@".$qual;

			$readRank ++;					
		}

		print "Total number of ambiguous nucleotides(e.g., N) changed to $ch: $ambig_changed_cnt \n";	

		close(FILE);
		close(OFa);
		close(ORQual);
	}
}
if ($fileNum == 0) {
	print "No .fastq files have been found in the directory... make sure the file extension is .fastq\n\n";
	exit;
}
print "\nDone!! $readRank reads remaining; $discarded discarded; $fileNum files\n\n";

close(DIR);
close(ODIR);

#--------------------------------------------
 

