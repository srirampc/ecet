#!/usr/bin/perl

# This program takes as input HiTec output (corrected) file (.fa),
# original short reads file (.fa), 
# Then output to [errReport]
#
#  File : coral-analy.pl
#  Created on December 1, 2011
#  Author : Xiao Yang <xyang@gmail.com>
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
  print "\nusage: perl coral-analy.pl [coral-corrected.fa] [sr-ori.fa] [O/errReport] \n\n";
  exit;
 }

 print "\nStart... \n\n";
 
 open(CoralFILE, "<$ARGV[0]") || die "can't open coral output file: $ARGV[0]\n";
 open(SrFILE, "<$ARGV[1]") || die "can't open original short reads file: $ARGV[1]\n";
 open(ErrFILE, ">$ARGV[2]") || die "can't open O/errReport file: $ARGV[2]\n";

 my $srCnt = 0; 
 my $stopFlag = 0;  # whether one of the file is ending 
 my $correctedCnt = 0;
 my $progress = 0; 

 my $batch  = 100000;
 print "progress: ($batch reads/unit) \t";
 while (<CoralFILE>){
	$progress ++;
	if ($progress % $batch == 0) {
		print ($progress/$batch);
		print "\t";
	}
	if (/>/){
		my $header = $_;

 		my $coralRead;
		my $srRead;
		if (!eof(CoralFILE)) { #check eof of CoralFILE;
			$coralRead = <CoralFILE>;
			if (!eof(SrFILE)){
				my $tmp = <SrFILE>; # skip the header
				$srRead = <SrFILE>;
			}
		}
		else { last; }

		if ($header =~ />\s\w+\s(\d+)\s\(\s(\d+)/ || 
                    $header =~ />(\d+)\s\(\s(\d+)/ ){
			my $coralReadID = $1;
			my $rhd = $2;  #reported hd by coral

			$srRead = lc($srRead);
			$coralRead = lc($coralRead); #convert to lower cases
            chomp($srRead);
            chomp($coralRead);

			if ($srRead ne $coralRead){ #correction found
				
				if (length($srRead) != length($coralRead)){
					print "\n ins/del found in read $coralReadID\n";
				}
				my $thd = hd($srRead, $coralRead); # true calculated hd
				if ($thd >= 5){
					print "thd=$thd  "."id=$coralReadID\n";
				}
				if ($thd != $rhd) {	# recording by coral is wrong 					
					print $coralReadID."   ".$thd."   ".$rhd."\n";					
				}

				$correctedCnt ++;

				my @charsFrom = split(//, $coralRead); #reference genome, here, corrected Reads
				my @charsTo = split(//, $srRead);
				
                print ErrFILE $coralReadID."\t".$thd;

		        for (my $i = 0; $i < length($srRead); ++ $i){

					if ($charsFrom[$i] ne $charsTo[$i]){ # found mismatch, write to output
						print ErrFILE "\t", $i, "\t", char_to_num($charsFrom[$i]), "\t", 
						char_to_num($charsTo[$i]), "\t", 0;
                    }
				} 
				print ErrFILE "\n";
			}
		}
        else {
			print "\nHeader=".$header."-- err occurred parsing header\n";
			exit;
		}
	} #if (/>/)
 } # while ()

 print "\n\ntotal corrected reads by Coral = ", $correctedCnt."\n\n";

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
