#!/usr/bin/perl

# This program takes as input HiTec output (corrected) file (.fa),
# original short reads file (.fa), 
# Then output to [errReport]
#
#  File :hitec-analy.pl
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
  print "\nusage: perl hitec-analy.pl [hitec-corrected.fa] [sr-ori.fa] [O/errReport] \n\n";
  exit;
 }

 print "\nStart... \n\n";
 
 open(HitecFILE, "<$ARGV[0]") || die "can't open shrec output file: $ARGV[0]\n";
 open(SrFILE, "<$ARGV[1]") || die "can't open original short reads file: $ARGV[1]\n";
 open(ErrFILE, ">$ARGV[2]") || die "can't open O/errReport file: $ARGV[2]\n";

 my $srCnt = 0; 
 my $stopFlag = 0;  # whether one of the file is ending 
 my $correctedCnt = 0;
 my $progress = 0; 

 my $batch  = 100000;
 print "progress: ($batch reads/unit) \t";
 while (<HitecFILE>){
	$progress ++;
	if ($progress % $batch == 0) {
		print ($progress/$batch);
		print "\t";
	}
	if (/>/){

		$_ =~ />(\d+)/;
		my $hitecReadID = $1;	
	 
		my $hitecRead;
		if (!eof(HitecFILE)) { #check eof of HitecFILE;
			$hitecRead = <HitecFILE>;
		}
		else { last; }

		if (looks_like_number($hitecReadID)){		
			#search in SrFile this 
			while (<SrFILE>){				
				my $header = $_;
				my $srRead;
				if (!eof(SrFILE)){ #check eof of SrFILE
					$srRead = <SrFILE>;
				}
				else { $stopFlag = 1; last; }				
				
				if ($srCnt == $hitecReadID) { # compare and output to [errFile]

					$srRead = lc($srRead); 
					$hitecRead = lc($hitecRead); #convert to lower cases
					chomp($srRead);
					chomp($hitecRead);

					if ($srRead ne $hitecRead){

						my @charsFrom = split(//, $hitecRead); #reference genome, here, corrected Reads
						my @charsTo = split(//, $srRead);
					
						my $hammingdist = hd($srRead, $hitecRead);	
				
						if ($hammingdist > 0) {
				
							$correctedCnt ++;

							print ErrFILE $hitecReadID."\t".$hammingdist;					

							for (my $i = 0; $i < length($srRead); ++ $i){

								if ($charsFrom[$i] ne $charsTo[$i]){ # found mismatch, write to output
									print ErrFILE "\t", $i, "\t", char_to_num($charsFrom[$i]), "\t", 
											char_to_num($charsTo[$i]), "\t", 0;							
								}
							}							
							print ErrFILE "\n";																		
						}					
					}
				}
				elsif ($srCnt > $hitecReadID) {
					print "Err: can't find hitecReadID $hitecReadID \n"; exit;
				}

				++ $srCnt;
				last;
			} # while
		}# if (looks_like_number)
 	
		if ($stopFlag == 1) { last; }
	}
 }

 print "\n\ntotal corrected reads by Hitec = ", $correctedCnt."\n\n";

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
 
 close(HitecFILE);
 close(SrFILE);
 close(ErrFile);
