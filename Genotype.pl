#! /usr/bin/perl -w
use strict;
use IPC::Open3;
use IO::Select;
use JSON qw( decode_json );

my ($out, $err) = ( '', '' );

my $loci_f = shift;
my %loci = ();
my @loci_arr = ();

&read_locilist ($loci_f, \@loci_arr, \%loci);

my $curr_loci = shift(@loci_arr); 
print "$curr_loci curr_loci\n";

my $var_f = shift;

my $pid_var = open3(\*VARCHLD_IN, \*VARCHLD_OUT, \*VARCHLD_ERR, "bzcat $var_f")
#my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, 'cat')
    or die "open3() failed $!";

my $reader_var = IO::Select->new(\*VARCHLD_OUT);

my @readyvar = $reader_var->can_read();
my $fhvar = $readyvar[0];

my $cnt = 0;
open (OUT, ">Genotype") || die;


#4858615	2	1	chr5	77396835	77396838	ref	TTC	TTC	564	563	PASS
#4858615	2	2	chr5	77396835	77396838	del	TTC		564	563	PASS
my $loci = "";
while (<$fhvar>) {
	#print $_;
	if (/^\d+\t+(\d+)\t+(\S+)\t+(chr\S+)\t+(\d+)\t+(\d+)\t(\S*?)\t(\S*?)\t(.*?)\t\d*\t\d*\t(\S*)/) {
		$cnt ++;
		my ($ploidy, $allele, $chr, $begin, $end, $type, $ref, $var, $qual) = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
		
		if ($cnt%10000 == 0) {
			print "$cnt finished ($chr, $begin, $end, $ref, $var)";
		}

		#if ($end != 35683240) {
		#	next;
		#}
		my $match_flag = &isMatch($chr, $begin, $end, $curr_loci);
		if ($match_flag > 0) {
			next;
		} elsif ($match_flag < 0) {
			if (@loci_arr == 0) {
				last;
			}

			while ($match_flag < 0) {
				print "($chr, $begin, $end, $curr_loci) $match_flag bf\n";
				if (@loci_arr > 0) {
					$curr_loci = shift(@loci_arr);
				} else {
					last;
				}
				$match_flag = &isMatch($chr, $begin, $end, $curr_loci);
				print "($chr, $begin, $end, $curr_loci) $match_flag af\n";
			}
		}

		if ($match_flag == 0) {
			print "match ($chr, $begin, $end, $ref, $var)\n";

			my ($rsID, $rs_ref) = split (/\s+/, $loci{"$curr_loci"});

			if ($allele eq 'all') {
				if (($ref eq '=') && ($var eq '=')) {
					print OUT "$rsID\t$rs_ref$rs_ref\n";
				} elsif ($var eq '?') {
					print OUT "$rsID\tNR\n";
				} else {
					print "Impossible allele all\n";
				}
			} elsif (($allele == 1) && ($ploidy == 1)) {
				print "Skip ($allele == 1) && ($ploidy == 1)\n";
			} elsif (($allele == 1) && ($ploidy == 2)) {
				if (($type ne 'ref') && ($type ne 'snp')) {
					$loci = "not a snp $type $var";
				} else {
					$loci = $var;
				}
			} elsif (($allele == 2) && ($ploidy == 2)) {
				if ((($type ne 'ref') && ($type ne 'snp')) || 
					($loci =~ /not a snp/)) {
					#print OUT "$rsID\tNA\t$type $var\t$loci\n";
					print OUT "$rsID\tNA\n";
				} else {
					if ($loci cmp $var) {
						$loci = $var.$loci;
					} else {
						$loci .= $var;
					}

					print OUT "$rsID\t$loci\n";
				}

				
			} else {
				print "Impossible\n";
			}
		}
	}
#	exit;
}
close (OUT);
exit;

sub isMatch () {
	my ($chr, $begin, $end, $curr_loci) = @_;
	my $rtn = 0;

	my ($loci_chr, $loci_begin, $loci_end) = split (/\s+/, $curr_loci);

	$chr =~ s/chrX/chr23/g;
	$loci_chr =~ s/chrX/chr23/g;
	$chr =~ s/chr//g;
	$loci_chr =~ s/chr//g;

	if ($chr < $loci_chr) {
		$rtn = 1;
	} elsif ($chr > $loci_chr) {
		$rtn = -1;
	} else {
		if ($end < $loci_end) {
			$rtn = 1;
		} elsif ($end == $loci_end) {
			$rtn = 0;
		} else {
			if ($begin >= $loci_end) {
				$rtn = -1;
			} else {
				$rtn = 0;
			}
		}
	}

	return $rtn;
}

exit;


sub parse_varfileLine () {
	my ($line) = @_;
	my ($chr, $begin, $end, $ref, $var);

	if ($line =~ /^\d+\s+\d+\s+\S+\s+(chr\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
		($chr, $begin, $end, $ref, $var) = ($1, $2, $3, $4, $5);
	}

	return ($chr, $begin, $end, $ref, $var);
}


sub read_locilist () {
	my ($loci_f, $loci_arr_r, $loci_r) = @_;

	open (LOCI, "<$loci_f") || die;
	while (<LOCI>) {
		chomp;

		if (/^(rs\S+)\t(chr\S+)\t(\d+)\t(\d+)\t(\S+)/) {
			my ($rsID, $chr, $begin, $end, $ref) = ($1, $2, $3, $4, $5);

			$loci_r->{$chr." ".$begin." ".$end} = $rsID." ".$ref;
			push (@$loci_arr_r, $chr." ".$begin." ".$end);
		}
	}
	close (LOCI);
}
