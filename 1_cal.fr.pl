#!/usr/bin/perl -w
use strict;
my $vcf  = shift;   
my $fr   = shift;

my ($i,$j);
my @tm=();

my @sp_name  = qw (Niv1 Ruf1 Niv2 Ruf2);
my @pop_name = qw (Niv1b Niv2a Niv2b Niv2c Ruf1b Ruf1a Ruf2a Ruf2b);
my %pop_code=();

$pop_code{$pop_name[0]}='Niv1b_';
$pop_code{$pop_name[1]}='Niv2a_';  
$pop_code{$pop_name[2]}='Niv2b_';
$pop_code{$pop_name[3]}='Niv2c_';
$pop_code{$pop_name[4]}='Ruf1b_';
$pop_code{$pop_name[5]}='Ruf1a_';
$pop_code{$pop_name[6]}='Ruf2a_';
$pop_code{$pop_name[7]}='Ruf2b_';

my %pop_sam=();
   @{$pop_sam{'ALL'}}=();
   @{$pop_sam{$sp_name[0]}}=();
   @{$pop_sam{$sp_name[1]}}=();
   @{$pop_sam{$sp_name[2]}}=();
   @{$pop_sam{$sp_name[3]}}=();
foreach(@pop_name){@{$pop_sam{$_}}=();}   
my %sam_pos=();
open(R,$vcf) || die "Reading $vcf failed ..\n";
open(W,">".$fr) || die "Writing $vcf failed ..\n";
while(<R>)
     { chomp $_;
       if($_ =~ /^\#\#/){next;}
       if($_ =~ /^\#CHROM/) 
         { @tm = split(/\t/, $_);
           printf W "%s\t%s\t%s\t%s\t%s", $tm[0],$tm[1],$tm[3],$tm[4],"ALL";
           foreach(@sp_name) {printf W "\t%s", $_;}
           print W "\n";
           for $i(9..scalar(@tm)-1) 
                 { $sam_pos{$tm[$i]}=$i;
                   push(@{$pop_sam{'ALL'}},$tm[$i]);
                   for $j(0..scalar(@pop_name)-1)                 
                         {if($tm[$i] =~ /^$pop_code{$pop_name[$j]}/)
                            {push(@{$pop_sam{$pop_name[$j]}}, $tm[$i]);}
                         }
                  }        
            push(@{$pop_sam{$sp_name[0]}},@{$pop_sam{$pop_name[0]}});
            push(@{$pop_sam{$sp_name[1]}},@{$pop_sam{$pop_name[4]}});
            push(@{$pop_sam{$sp_name[1]}},@{$pop_sam{$pop_name[5]}});
            push(@{$pop_sam{$sp_name[2]}},@{$pop_sam{$pop_name[1]}});
            push(@{$pop_sam{$sp_name[2]}},@{$pop_sam{$pop_name[2]}});
            push(@{$pop_sam{$sp_name[2]}},@{$pop_sam{$pop_name[3]}});
            push(@{$pop_sam{$sp_name[3]}},@{$pop_sam{$pop_name[6]}});
            push(@{$pop_sam{$sp_name[3]}},@{$pop_sam{$pop_name[7]}});          
           next;
          }
       @tm = split(/\t/, $_);
       if($tm[4] =~ /\*/) {next;}
       my @allele=split(/,/,$tm[4]);
       my $aln=scalar(@allele);
       printf W "%s\t%s\t%s\t%s", $tm[0], $tm[1], $tm[3], $tm[4];
       my @seg=();
       foreach(@{$pop_sam{'ALL'}}) { push(@seg, $tm[$sam_pos{$_}]);} 
       my @fr=cal_fr($aln, \@seg);
       my $frn=$fr[0];
       for $i(1..scalar(@fr)-1){$frn .= ':'.$fr[$i];}
       printf W "\t%s", $frn;
       

       
       foreach(@sp_name)
              { @seg=();
                foreach(@{$pop_sam{$_}}) {push(@seg, $tm[$sam_pos{$_}]);}
                @fr=cal_fr($aln, \@seg);
                $frn=$fr[0];
                for $i(1..scalar(@fr)-1){$frn .= ':'.$fr[$i];}
                printf W "\t%s", $frn;
              }  
       print W "\n"; 
       }                                              
close(R);
close(W);
         
sub   cal_fr  
       { my ($an, $seg)=@_;
         my %ale=();
         my $i=0;
         for $i(0..$an) {$ale{$i}=0;}
         foreach (@$seg)
                 { if( $_ =~ /^\./) {next;} 
                   my @tk=split(/:/,$_);
                   my @tn=();
                   if($tk[0] =~ /\//){@tn=split(/\//, $tk[0]);}
                   else              {@tn=split(/\|/, $tk[0]);}                           
                   $ale{$tn[0]}++; 
                   $ale{$tn[1]}++;
                 }
         my @sfr=();
         for  $i(0..$an) { push (@sfr, $ale{$i}); }
         return @sfr;
        }
                 
                         
             


