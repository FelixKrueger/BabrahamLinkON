#!/usr/bin/perl
use warnings;
use strict;

die "Please provide one or more J-end FastQ files for the mouse IgK mispriming correction!\n\n" unless (@ARGV);
warn "Performing mispriming correction for the mouse IgK locus\n\n";
# %fr contains sequences to find and replace
my %fr = read_find_and_replace_sequences();

# Setting count to 0
foreach my $seq (keys %fr){
    $fr{$seq} -> {count} = 0;
    #warn "$seq\t$fr{$seq}->{seq}\t$fr{$seq}->{count}\n";
}

foreach my $in (@ARGV){
    pre_process_bait_file_for_mispriming($in);
}

### Printing a detailed summary
warn "\nSequence detected\tReplaced with\tCount\n";
print REPORT "\nSequence detected\tReplaced with\tCount\n";

# sorting alphabetically
foreach my $seq (sort {$a cmp $b}  keys %fr){
    warn "$seq\t$fr{$seq}->{seq}\t$fr{$seq}->{count}\n";
    print REPORT "$seq\t$fr{$seq}->{seq}\t$fr{$seq}->{count}\n";
}

sub pre_process_bait_file_for_mispriming {

    ### this subroutine will look at the start of every read 2 reads for any of the sequences strored in %fr. If we find a mispriming event we simply replace 
    ### the bait sequence with the sequence stored in $fr{seq} and write this out to a mispriming corrected version of the Read 2 file. This is then processed normally.

    my $infile = shift;
    my $mispriming_replaced;
    
    my $outfile_read2 = $infile;
    $outfile_read2 =~ s/\.gz$//;
    $outfile_read2 =~ s/\.fq$//;
    $outfile_read2 =~ s/\.fastq$//;
    my $report_file = $outfile_read2;
    $outfile_read2 .= '.mispriming_corrected.fq.gz';
    $report_file   .= '.mispriming_correction.txt';

  #  die "Outfile $outfile_read2 already existed\n\n" if (-e $outfile_read2);

    warn "Reading from bait file '$infile'\n";
    if ($infile =~ /gz$/){
	open (MISPRIMING,"zcat $infile |") or die "Couldn't read from gzipped file '$infile'\n";
    }
    else{
	open (MISPRIMING,$infile) or die "Couldn't read from file '$infile'\n";
    }
    warn "Writing mispriming corrected version of the bait file to '$outfile_read2'\n";
    warn "Writing mispriming correction report file to '$report_file'\n\n";
    open (OUT,"| gzip -c - > $outfile_read2" ) or die "Can't write to file '$outfile_read2': $!\n";
    open (REPORT,'>',$report_file) or die "Can't write report file '$report_file': $!\n";

    my $total = 0;
    my $this_is_j2 = 0;
    my $j_too_short= 0;
    
    while (1){
	
	my $id    = <MISPRIMING>;
	my $seq   = <MISPRIMING>;
	my $three = <MISPRIMING>;
	my $qual  = <MISPRIMING>;

	last unless ($qual);
	$total++;

	if ($total%250000 == 0 ){
	    warn "Processed $total sequences so far...\n";
	}
	chomp ($seq);
	my $replacement_seq = $seq;
	chomp $qual;
	my $replacement_qual = $qual;

	foreach my $fr (keys %fr){
	    # warn "$fr\t$fr{$fr}->{seq}\t$fr{$fr}->{count}\n"; sleep(1);
	    
	    if ($seq =~ /^$fr/){
		# replacing the J2 with the J4 bait sequence. Note that TCCATA is part of the J4 bait we are looking for
		$replacement_seq =~ s/^$fr/$fr{$fr}->{seq}/;

		#warn "TBR     $fr{$fr}->{seq}\n";
		# warn "Seq   : $fr\n";
		#warn "before: $replacement_qual\n";
		substr ($replacement_qual,0,length($fr),''); # deleting the old sequence
		#warn "middle: $replacement_qual\n"; sleep(1);
		$replacement_qual = 'I' x length(($fr{$fr}->{seq})).$replacement_qual;
		#warn "after:  $replacement_qual\n";
	
		++$mispriming_replaced;
		++$fr{$fr} -> {count}; ## increasing this find/replace counter
		
		#warn "original sequence:\t$seq\nreplaced sequence:\t$replacement_seq\nreplaced quality:\t$replacement_qual\n\n"; sleep(1);
		last; # exiting the loop after the first replacement
	    }
	}

	### printing out the replaced sequence (which is the same if there was no replacement)
	if (length $replacement_qual != length($replacement_seq)){
	    die "Length is not the same:\n$replacement_seq\n$replacement_qual\n\n";sleep(1);
	}
		
	print OUT $id;
	print OUT "$replacement_seq\n";
	print OUT $three;
	print OUT "$replacement_qual\n";
	
    }

    my $percent_mispriming = sprintf ("%.1f",$mispriming_replaced/$total*100);

    warn "\nPre-processing report for bait file $infile (Read 2)\n";
    warn "========================================================\n";
    warn "Total number of sequences analysed:\t$total\n";
    warn "Cases detected as mispriming in total: $mispriming_replaced ($percent_mispriming%)\n";

    print REPORT "\nPre-processing report for bait file $infile (Read 2)\n";
    print REPORT "========================================================\n";
    print REPORT "Total number of sequences analysed:\t$total\n";
    print REPORT "Cases detected as mispriming in total: $mispriming_replaced ($percent_mispriming%)\n";

    close OUT or warn "Failed to close read2 mispriming corrected file\n";
}

sub read_find_and_replace_sequences{

    # warn "Now storing sequences to find/replace\n\n";
 
    ### sequences from file 'mouse_IgK_mispriming_correction.txt'
    my %fr;
    ## FOr the run in question we had to trim the first 4bp off because of poor basepair qualities. Thus we also need to shorten the find and replace sequences
    ## I'll leave the full length version in as well in case we need them again at a later point

    #Jk1-Jk2 full length
    # $fr{'TTTGATTTCCAGCTTGGTGCCTCCT'}          -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCT';
    # $fr{'TTTGATTTCCAGCTTGGTCCCCCCT'}          -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCT';
    # $fr{'TTTGATTTCCAGCTTGGTGCCCCCT'}          -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCT';
    #Jk1-Jk2
    $fr{'ATTTCCAGCTTGGTGCCTCCT'}          -> {seq} = 'TTTCCAGCTTGGTCCCCCCT';
    $fr{'ATTTCCAGCTTGGTCCCCCCT'}          -> {seq} = 'TTTCCAGCTTGGTCCCCCCT';
    $fr{'ATTTCCAGCTTGGTGCCCCCT'}          -> {seq} = 'TTTCCAGCTTGGTCCCCCCT';

    #Jk1 mm full length
    # $fr{'TTTGATTTCCACCTTGGTGCCTCCACCGAACGTC'} -> {seq} = 'TTTGATTTCCAGCTTGGTGCCTCCACCGAACGTC';
    #Jk1 mm 
    $fr{'ATTTCCACCTTGGTGCCTCCACCGAACGTC'} -> {seq} = 'ATTTCCAGCTTGGTGCCTCCACCGAACGTC';
 

    #Jk1-Jk5 full length
    # $fr{'TTTGATTTCCAGCTTGGTCCCAGC'}           -> {seq} = 'CAGCTCCAGCTTGGTCCCAGC';
    # $fr{'TTTGATTTCCAGCTTGGTCACAGC'}           -> {seq} = 'CAGCTCCAGCTTGGTCCCAGC';
    # $fr{'TTTGATTTCCAGCTTGGTGCCAGC'}           -> {seq} = 'CAGCTCCAGCTTGGTCCCAGC';
    #Jk1-Jk5
    $fr{'ATTTCCAGCTTGGTCCCAGC'}           -> {seq} = 'TCCAGCTTGGTCCCAGC';
    $fr{'ATTTCCAGCTTGGTCACAGC'}           -> {seq} = 'TCCAGCTTGGTCCCAGC';
    $fr{'ATTTCCAGCTTGGTGCCAGC'}           -> {seq} = 'TCCAGCTTGGTCCCAGC';
  
    #Jk2-Jk1 full length
    # $fr{'TTTATTTCCAGCTTGGTGCCTCCA'}           -> {seq} = 'TTTGATTTCCAGCTTGGTGCCTCCA';
    #Jk2-Jk1
    $fr{'TTTCCAGCTTGGTGCCTCCA'}           -> {seq} = 'ATTTCCAGCTTGGTGCCTCCA';
  
    #Jk2-Jk4 full length
    # $fr{'TTTATTTCCAGCTTGGTCCCCGAGCC'}         -> {seq} = 'CGTTTTATTTCCAACTTTGTCCCCGAGCC';
    # $fr{'TTTATTTCCAGCTTGGTCCCAGCC'}           -> {seq} = 'CGTTTTATTTCCAACTTTGTCCCCGAGCC';
    #Jk2-Jk4
    $fr{'TTTCCAGCTTGGTCCCCGAGCC'}         -> {seq} = 'TTATTTCCAACTTTGTCCCCGAGCC';
    $fr{'TTTCCAGCTTGGTCCCAGCC'}           -> {seq} = 'TTATTTCCAACTTTGTCCCCGAGCC';
    
    #Jk2-Jk5 full length
    # $fr{'TTTATTTCCAGCTTGGTCCCAGCACC'}         -> {seq} = 'CAGCTCCAGCTTGGTCCCAGCACC';
    # $fr{'TTTATTTCCAGCTTGGTGCCAGCACC'}         -> {seq} = 'CAGCTCCAGCTTGGTCCCAGCACC';
    # $fr{'TTTATTTCCAGCTTGGTCCAAGCACC'}         -> {seq} = 'CAGCTCCAGCTTGGTCCCAGCACC';
    #Jk2-Jk5
    $fr{'TTTCCAGCTTGGTCCCAGCACC'}         -> {seq} = 'TCCAGCTTGGTCCCAGCACC';
    $fr{'TTTCCAGCTTGGTGCCAGCACC'}         -> {seq} = 'TCCAGCTTGGTCCCAGCACC';
    $fr{'TTTCCAGCTTGGTCCAAGCACC'}         -> {seq} = 'TCCAGCTTGGTCCCAGCACC';
    
    #Jk4-Jk1 full length
    # $fr{'CGTTTTATTTCCAACTTTGTCCCCGAACGTC'}    -> {seq} = 'TTTGATTTCCAGCTTGGTGCCTCCACCGAACGTC';
    #Jk4-Jk1
    $fr{'TTATTTCCAACTTTGTCCCCGAACGTC'}    -> {seq} = 'ATTTCCAGCTTGGTGCCTCCACCGAACGTC';
   
    #Jk4-Jk2 full length
    # $fr{'CGTTTTATTTCCAACTTTGTCCCCCCT'}        -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCT';
    #Jk4-Jk2
    $fr{'TTATTTCCAACTTTGTCCCCCCT'}        -> {seq} = 'TTTCCAGCTTGGTCCCCCCT';
   
    #Jk4-Jk5 full length
    # $fr{'CGTTTTATTTCCAACTTTGTCCCCGAACGTG'}    -> {seq} = 'CAGCTCCAGCTTGGTCCCAGCACCGAACGTG';
    # $fr{'CGTTTTATTTCCAACTTTGTCCCCGAACCGAACGTG'}->{seq} = 'CAGCTCCAGCTTGGTCCCAGCACCGAACGTG';
    # $fr{'CGTTTTATTTCCAACTTTGTCCCAGCACCGAACGTG'}->{seq} = 'CAGCTCCAGCTTGGTCCCAGCACCGAACGTG';
    # $fr{'CGTTTTATTTCCAACTTTGTCCCCGCACCGAACGTG'}->{seq} = 'CAGCTCCAGCTTGGTCCCAGCACCGAACGTG';
    #Jk4-Jk5
    $fr{'TTATTTCCAACTTTGTCCCCGAACGTG'}    -> {seq} = 'TCCAGCTTGGTCCCAGCACCGAACGTG';
    $fr{'TTATTTCCAACTTTGTCCCCGAACCGAACGTG'}->{seq} = 'TCCAGCTTGGTCCCAGCACCGAACGTG';
    $fr{'TTATTTCCAACTTTGTCCCAGCACCGAACGTG'}->{seq} = 'TCCAGCTTGGTCCCAGCACCGAACGTG';
    $fr{'TTATTTCCAACTTTGTCCCCGCACCGAACGTG'}->{seq} = 'TCCAGCTTGGTCCCAGCACCGAACGTG';
  
    #Jk5-Jk1 full length
    # $fr{'CAGCTCCAGCTTGGTGCCTCCA'}             -> {seq} = 'TTTGATTTCCAGCTTGGTGCCTCCA';
    #Jk5-Jk1
    $fr{'TCCAGCTTGGTGCCTCCA'}             -> {seq} = 'ATTTCCAGCTTGGTGCCTCCA';
   
    #Jk5-Jk2 full length
    # $fr{'CAGCTCCAGCTTGGTCCCAGCTCC'}           -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCTCC';
    # $fr{'CAGCTCCAGCTTGGTCCCCCCTCC'}           -> {seq} = 'TTTATTTCCAGCTTGGTCCCCCCTCC';
    #Jk5-Jk2
    $fr{'TCCAGCTTGGTCCCAGCTCC'}           -> {seq} = 'TTTCCAGCTTGGTCCCCCCTCC';
    $fr{'TCCAGCTTGGTCCCCCCTCC'}           -> {seq} = 'TTTCCAGCTTGGTCCCCCCTCC';
   
    #Jk5-Jk4 full length
    # $fr{'CAGCTCCAGCTTGGTCCCAGCCG'}            -> {seq} = 'CGTTTTATTTCCAACTTTGTCCCCGAGCCG';
    # $fr{'CAGCTCCAGCTTGGTCCCCGAGCCG'}          -> {seq} = 'CGTTTTATTTCCAACTTTGTCCCCGAGCCG';
    #Jk5-Jk4
    $fr{'TCCAGCTTGGTCCCAGCCG'}            -> {seq} = 'TTATTTCCAACTTTGTCCCCGAGCCG';
    $fr{'TCCAGCTTGGTCCCCGAGCCG'}          -> {seq} = 'TTATTTCCAACTTTGTCCCCGAGCCG';

    return %fr;
}
