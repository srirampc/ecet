#!/usr/bin/perl

# This program takes as input Shrec output (corrected) file (.fa),
# original short reads file (.fa), 
# Then output to [errReport]
#
#  File : SHREC-Analy.pl
#  Created on December 1, 2011
#  Author: Xiao Yang <isuxyang@gmail.com>
#
#  This file is part of Error Correction Review Toolkit.
#  Error Correction Review Toolkit is free software: you can 
#  redistribute it and/or modify  it under the terms of the GNU 
#  Lesser General Public License as published by 
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Error Correction Review Toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
#

 use Scalar::Util qw(looks_like_number);
 use Data::Dumper;
 use Switch;

 my $numArgs = $#ARGV + 1;
 
 if($numArgs != 3) {	 
  print "\nusage: shrec-analy.pl [shrec-corrected.fa] [sr-ori.fa] [O/errReport] \n\n";
  exit;
 }

 print "\nStart... \n\n";
 
 open(ShrecFILE, "<$ARGV[0]") || die "can't open shrec output file: $ARGV[0]\n";
 open(SrFILE, "<$ARGV[1]") || die "can't open original short reads file: $ARGV[1]\n";
 open(ErrFILE, ">$ARGV[2]") || die "can't open O/errReport file: $ARGV[2]\n";

 my $srCnt = 0; 
 my $stopFlag = 0;  # whether one of the file is ending 
 my $correctedCnt = 0;
 my $progress = 0; 

 my $batch  = 100000;
 print "progress: ($batch reads/unit) \t";
 while (<ShrecFILE>){
	$progress ++;
	if ($progress % $batch == 0) {
		print ($progress/$batch);
		print "\t";
	}
	if (/>/ && /corrected/){

		$correctedCnt ++; 

		# this is shrec v1.0 format
		#my @shrecline = split(/\s+/);
		#my $shrecReadID = $shrecline[3];
		# to accomodate shrec v2.0 output format changed to the following lines:
		$_ =~ />(\d+)\s+/;
		my $shrecReadID = $1;	
		#print $shrecReadID."\t";	
	 
		my $shrecRead;
		if (!eof(ShrecFILE)) { #check eof of ShrecFILE; deal with corrupted Shrecfile
			$shrecRead = <ShrecFILE>;
		}
		else { last; }

		if (looks_like_number($shrecReadID)){		
			#search in SrFile this 
			while (<SrFILE>){				
				my $header = $_;
				my $srRead;
				if (!eof(SrFILE)){ #check eof of SrFILE
					$srRead = <SrFILE>;
				}
				else { $stopFlag = 1; last; }				
				
				if ($srCnt == $shrecReadID) { # compare and output to [errFile]

					$srRead = lc($srRead); 
					$shrecRead = lc($shrecRead); #convert to lower cases
					chomp($srRead);
					chomp($shrecRead);
					#print "--------------------------------\n";
					#print $srRead."\n";
					#print $shrecRead."\n";
					#if ($correctedCnt == 4) {
					#	exit;
					#}
					my @charsFrom = split(//, $shrecRead); #reference genome, here, corrected Reads
					my @charsTo = split(//, $srRead);
					
					print ErrFILE $shrecReadID."\t".hd($srRead, $shrecRead);					

					for (my $i = 0; $i < length($srRead); ++ $i){

						if ($charsFrom[$i] ne $charsTo[$i]){ # found mismatch, write to output
							print ErrFILE "\t", $i, "\t", char_to_num($charsFrom[$i]), "\t", 
									char_to_num($charsTo[$i]), "\t", 0;							
						}
					}							
					print ErrFILE "\n";	
					++ $srCnt;					
				    last;
				}
				elsif ($srCnt > $shrecReadID) {
					print "Err: can't find shrecReadID $shrecReadID \n"; exit;
				}
				++ $srCnt;
			}
		}	
		if ($stopFlag == 1) { last; }
	}
 }

 print "\n\ntotal corrected reads by Shrec = ", $correctedCnt."\n\n";

 print "program finished successfully!\n\n";
 #---------------------------------------
 # hamming distance
 #---------------------------------------
 sub hd
 {
     #String length is assumed to be equal
     my ($k,$l) = @_;
     my $len = length ($k);
     my $num_mismatch = 0;

     for (my $i=0; $i<$len; $i++)
     {
      ++$num_mismatch if substr($k, $i, 1) ne substr($l, $i, 1);
     }
     return $num_mismatch;
 }

 sub char_to_num {
	switch($_[0]) {
		case 'a' { return 0; }
		case 'c' { return 1; }
		case 'g' { return 2; }
		case 't' { return 3; }
	}
 }
 
 close(ShrecFILE);
 close(SrFILE);
 close(ErrFile);
