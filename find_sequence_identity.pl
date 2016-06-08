#!/project/michal/apps/perl/bin/perl
use strict;
use warnings ;
use Benchmark;
use File::Slurp;
use Cwd;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Algorithm::NeedlemanWunsch;
use List::Util qw( min max );
use AI::NaiveBayes1;
use YAML;

# This scripts finds global identity of any two proteins.
# Input is the pdb of both the proteins.
# Output is the global seq identity.

# $ARGV[0] : pdb of the receptor
# $ARGV[1] : pdb of the ligand



##section1###############################################################################################
## Translate PDB files to sequence string file
## usage > pdb2seq(@input_pdb_array);

 sub pdb2seq
  {
  my %aa=qw(
	ALA 	A
	CYS 	C
	ASP 	D
	GLU 	E 
	PHE 	F
	GLY	G
	HIS	H
	ILE	I
	LYS	K
	LEU	L
	MET	M
	ASN	N
	PRO	P
	GLN	Q
	ARG	R
	SER	S
	THR	T
	VAL	V
	TRP	W
	TYR	Y);

  my @pdb_in = @_;my $residue ; my $resno;
  my $oldresno=-1;my $seq=(); #print "@pdb_in";
  foreach (@pdb_in)
   {
    if (/^ATOM/)
    {
     my $type = substr($_,13,2); 
     if ($type eq "CA")	
      {	
       my $res = substr($_, 17, 3); 
       chomp($residue=$aa{$res});$residue=~ s/^\s+//;$residue=~s/\s+$//;
       my $resno=substr($_, 22, 4);
       if ($resno>$oldresno)
         {
         $seq=$seq.$residue ; 
         $oldresno=$resno;
         }
       }
     }
    }
   return $seq;
   }

##section2###############################################################################################
# Gives the global sequence identity between two sequence strings.
# usage > get_identity(seq_string1, seq_string2);

 use vars qw( %nwmat4 @nwseq3 @nwseq4 $nwseq3o $nwseq4o );
 my %nwmat1 = (); 
 my @nwmat3 = qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
 
 while ( my $wdat1 = <DATA> )
 {
  chomp $wdat1;
  if ( length($wdat1) == 70 and $wdat1 ne '   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X' )
  {
   my $nwr1 = substr($wdat1, 0, 1);   
   for ( my $xg = 0; $xg < 23; $xg++ )
   {
    $nwmat1{$nwr1.$nwmat3[$xg]} = substr($wdat1, 1 + $xg * 3, 3) * 1;
   }
  }
 }
 
 my @nwseq1 = ();my @nwseq2 = ();
 my $nwseq1o = '';my $nwseq2o = '';
 
 sub blosum 
 { 
  my ($anw, $bnw) = @_;  
  my $snw = 0;  
  $snw = $nwmat1{$anw.$bnw} if ( exists $nwmat1{$anw.$bnw} );  
  return ($snw);
 }
 
 my $matcher1 = Algorithm::NeedlemanWunsch->new(\&blosum);
 $matcher1->gap_open_penalty(-10);
 $matcher1->gap_extend_penalty(-2);
 
 sub prepend_align1 
 { 
  my ($i, $j) = @_;
  $nwseq1o = $nwseq1[$i] . $nwseq1o;
  $nwseq2o = $nwseq2[$j] . $nwseq2o;
 }
 
 sub prepend_first_only1 
 { 
  my $i = shift;  
  $nwseq1o = $nwseq1[$i] . $nwseq1o;
  $nwseq2o = "-$nwseq2o";
 }
 
 sub prepend_second_only1 
 {
  my $j = shift;
  $nwseq1o = "-$nwseq1o";
  $nwseq2o = $nwseq2[$j] . $nwseq2o;
 }
 
 sub get_identity 
 {
  my ($iseq1, $iseq2) = @_;  
  my $iss = 0.0; 
  @nwseq1 = split(//, $iseq1);
  @nwseq2 = split(//, $iseq2);
  
  $nwseq1o = '';$nwseq2o = '';
 
  my $score = $matcher1->align(\@nwseq1, \@nwseq2, { align => \&prepend_align1, shift_a => \&prepend_first_only1, shift_b => \&prepend_second_only1, });
  
  my @nwseq1a = split(//, $nwseq1o);
  my @nwseq2a = split(//, $nwseq2o);
  
  my $niseq1 = @nwseq1a;my $niseq2 = @nwseq2a;
  my $isid1 = 0; my $isid2 = 0;
  
  if ( $niseq1 == $niseq2 )
  {
   for ( my $ixa = 0; $ixa < $niseq1; $ixa++ )
   {
    $isid1++ if ( $nwseq1a[$ixa] ne '-' and $nwseq2a[$ixa] ne '-' and $nwseq1a[$ixa] eq $nwseq2a[$ixa] );
    $isid2++ if ( $nwseq1a[$ixa] ne '-' and $nwseq2a[$ixa] ne '-' );
   }
  }
  
  $iss = $isid1 / $isid2 if ( $isid2 );  
  return $iss;
 }

##section3###############################################################################################
 my $seq_global = 0;
 my @receptor = read_file($ARGV[0]);
 my @ligand = read_file($ARGV[1]);
 my $rec_seq = pdb2seq(@receptor);
 my $lig_seq = pdb2seq(@ligand);
 $seq_global = get_identity($rec_seq,$lig_seq);

 print "$seq_global\n";
 
#########################################################################################################


__DATA__
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1
B -2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  2 -1
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  2  5 -1
X -1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1
