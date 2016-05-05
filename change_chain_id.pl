#!/usr/bin/perl -w

## $ARGV[0] : Input PDB
## $ARGV[1] : Desired Chain ID

 use strict;
 use warnings;
 use File::Slurp;

 my @file = read_file($ARGV[0]);
 foreach (@file)
  {
   if (/^ATOM/)
    {
     my $x1 = substr($_,0,21);
     my $x2 = substr($_,22,55);
     my $y= sprintf("$x1\%s$x2\%5.2f \%5.2f \n","$ARGV[1]",1,0);
     print $y;
   }
  }
