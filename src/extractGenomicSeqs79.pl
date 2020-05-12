#!/usr/bin/perl
my ($chrDir, $inF) = @ARGV;

$chrDir = $chrDir;
@chrFiles = ("chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa",
             "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", "chr11.fa", "chr12.fa",
             "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", "chr18.fa",
             "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chr23.fa", "chrX.fa",
             "chrY.fa");

%chrSeqs = ();

foreach $chr (@chrFiles) {

   open (CHR,$chrDir."/$chr");
   my $line = <CHR>;  # Throw away defline

   my $seq = "";
   while ($line = <CHR>) {

      chomp( $line );
      $seq = $seq.$line;

   }
   close (CHR);

   $chr =~ s/.fa$//g;
   $chrSeqs{ $chr } = $seq;
}
print "Read in data for ".scalar(keys %chrSeqs)." chromosomes\n";

# Read in guides and print out gene, strand, guide, sequence context, + pam
# Format: [23N]20mer[NGG][33N]
open (IN, $inF);
open (OUT, ">gRNA_seqs79_pam.txt");

my $line = <IN>;
while ($line = <IN>) {

   chomp($line);
   my ($CHROMOSOME, $START, $STOP, $STRAND, $GUIDE, $GENE) = split("\t", $line);

   if ($STRAND eq "-") {

      $START = $START - 36;
      $SEQUENCE = uc(substr($chrSeqs{$CHROMOSOME},$START,79));
      $revcomp = reverse($SEQUENCE);
      $revcomp =~ tr/ACGTNacgt/TGCANtgca/;
      $pam = substr($revcomp,43,3);
      print OUT $GENE."\t".$STRAND."\t".$GUIDE."\t".$revcomp."\t".$pam."\n";

   } else {

      $START = $START - 23;
      $SEQUENCE = uc(substr($chrSeqs{$CHROMOSOME},$START,79));
      $pam = substr($SEQUENCE,43,3);
      print OUT $GENE."\t".$STRAND."\t".$GUIDE."\t".$SEQUENCE."\t".$pam."\n";

   }
}
close (IN);
close (OUT);
