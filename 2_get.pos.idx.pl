#!/usr/bin/perl -w

use strict;

my $fr   = shift;

my %chr_pos = ();
my @chr_key=();
my @tm=();
my $pos_marker= 0; 
my $pos_value = 0;
my $chr_no = "chr0";
open(R,$fr);
while(<R>)
     { if( $_ =~ /\#/) {$pos_value = tell(R); next;}
       @tm = split(/\t/,$_);
       $pos_marker = int($tm[1]/10000);
       if($chr_no ne $tm[0]) { $chr_no = $tm[0];
                               if($pos_marker>0)
                                 { $chr_pos{$tm[0].'-0'} = $pos_value;
                                   push(@chr_key, $tm[0]."-0");
                                  }
                             }          
       if( !exists($chr_pos{$tm[0]."-".$pos_marker}) )
         { $chr_pos{$tm[0]."-".$pos_marker} = $pos_value;
           push(@chr_key, $tm[0]."-".$pos_marker);
          }
       $pos_value = tell(R);  
      }
open(W,"> your_file_name.pos.idx");
foreach (@chr_key) {printf W "%s\t%s\n", $_, $chr_pos{$_};}
close(W);   
