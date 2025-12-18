#!/usr/bin/perl -w

my $seg_name  = shift;   # 基因段名称
my $seg_chr   = shift;   # 染色体
my $seg_start = shift;   # 起始位置
my $seg_end   = shift;   # 结束位置
my $seg_wnd   = shift;   # 窗口大小
my $seg_step  = shift;   # 滑窗步长
my $out = shift;        # 输出文件

my $snp_fre   = "XX.fr";
my $snp_pos   = "XX.pos.idx";

my @tm = ();
my @tn = ();

my ($i, $j);

my @gene = ();
my $nw = int(($seg_end - $seg_start + 1 - $seg_wnd) / $seg_step);  # 根据步长调整窗口数
   
for $i (0..$nw) {
    $gene[$i]->[0] = $seg_name . '-' . $i;
    $gene[$i]->[1] = $seg_chr;
    $gene[$i]->[2] = $seg_start + $i * $seg_step;
    $gene[$i]->[3] = $seg_start + $i * $seg_step + $seg_wnd - 1;
    if ($gene[$i]->[3] > $seg_end) { $gene[$i]->[3] = $seg_end; }
}

my %pos_idx = ();
open(R, $snp_pos) || die "reading $gene_info failed ..\n";
while (<R>) {
    chomp $_;
    @tm = split(/\t/, $_);
    $pos_idx{$tm[0]} = $tm[1];
}
close(R);

open(R, $snp_fre) || die "reading $gene_info failed ..\n";
my $headline = <R>;
chomp $headline;
@tm = split(/\t/, $headline);
my @pop_name = qw (Niv1 Ruf1 Niv2 Ruf2);
my %pop_fr = ();
my %pop_pos = ();
for $i (5..scalar(@tm)-1) { $pop_pos{$tm[$i]} = $i; }
$i = 0;

open(W ,">".$out) || die " writing failed ..\n";
printf W "%s\t%s\t%s\t%s\t%s", 'segment','chr','start_pos','end_pos','snp_num';
foreach (@pop_name) {
    printf W "\t%s-pi(1kb)", $_;
    printf W "\t%s-theta-W(1kb)", $_;
    printf W "\t%s-Tajima-D/sample", $_;
}

# 所有两两组合
my @pair = (
    ['Niv1','Ruf1'],
    ['Niv1','Niv2'],
    ['Niv1','Ruf2'],
    ['Ruf1','Niv2'],
    ['Ruf1','Ruf2'],
    ['Niv2','Ruf2']
);

foreach my $p (@pair) {
    printf W "\t%s-%s-dxy(1kb)", $p->[0], $p->[1];
    printf W "\t%s-%s-fst", $p->[0], $p->[1];
}


print W "\n";

foreach (@gene) {
    foreach (@pop_name) { @{$pop_fr{$_}} = (); }
    my $cr = $_->[1];
    my $start_pos = $_->[2];
    if ($start_pos < 1) { $start_pos = 1; }
    my $end_pos = $_->[3];
    my $len = $end_pos - $start_pos + 1;
    printf W "%s\t%s\t%s\t%s", $_->[0], $cr, $start_pos, $end_pos;

    my $pos_marker = int($start_pos / 10000);
    while (!exists($pos_idx{$cr . '-' . $pos_marker})) { $pos_marker--; }
    seek(R, $pos_idx{$cr . '-' . $pos_marker}, 0);
    my $snm = 0;
    while (<R>) {
        chomp $_;
        @tm = split(/\t/, $_);
        if ($tm[0] eq $cr && $start_pos <= $tm[1] && $tm[1] <= $end_pos) {
            $snm++;
            foreach (@pop_name) { push(@{$pop_fr{$_}}, $tm[$pop_pos{$_}]); }
        }
        if (($tm[0] eq $cr && $tm[1] > $end_pos) || $tm[0] ne $cr) { last; }
    }
    printf W "\t%s", $snm;
    my %k = ();
    my %ss = ();
    my %sn = ();
    my %pi = ();
    my %th = ();
    my %tj = ();
    foreach (@pop_name) {
        $k{$_} = cal_k(\@{$pop_fr{$_}});
        $ss{$_} = cal_ss(\@{$pop_fr{$_}});
        $sn{$_} = cal_sn(\@{$pop_fr{$_}});
        $pi{$_} = 1000 * $k{$_} / $len;
        $th{$_} = 1000 * cal_theta($ss{$_}, $sn{$_}, $len);
        $tj{$_} = cal_tajimaD($k{$_}, $ss{$_}, $sn{$_});
    }
		my %dxy = ();
		my %fst = ();
		foreach my $p (@pair) {
			my ($p1, $p2) = @$p;
			my $dxy_val = cal_dxy(\@{$pop_fr{$p1}}, \@{$pop_fr{$p2}}, $len);
			$dxy_val = 1000 * $dxy_val;
			my $fst_val = 0;
			if ($dxy_val > 0) {
				$fst_val = (2 * $dxy_val - $pi{$p1} - $pi{$p2}) / (2 * $dxy_val);
				if ($fst_val < 0) { $fst_val = 0; }
			}
			$dxy{"$p1-$p2"} = $dxy_val;
			$fst{"$p1-$p2"} = $fst_val;
		}

    foreach (@pop_name) {
        printf W "\t%10.5f", $pi{$_};
        printf W "\t%10.5f", $th{$_};
        printf W "\t%10.5f/%d", $tj{$_}, $sn{$_};
    }
	foreach my $p (@pair) {
		my ($p1, $p2) = @$p;
		printf W "\t%10.5f", $dxy{"$p1-$p2"};
		printf W "\t%10.5f", $fst{"$p1-$p2"};
	}

    print W "\n";
}

# 计算k值
sub cal_k {
    my ($fr) = @_;
    my $k = 0;
    my @sfr = ();
    foreach (@$fr) {
        @sfr = split(/:/, $_);
        my $an = scalar(@sfr);
        my ($i, $j);
        my $mch = 0;
        for $i (0..$an-2) { for $j ($i+1..$an-1) { $mch += $sfr[$i] * $sfr[$j]; } }
        my $sun = 0;
        foreach (@sfr) { $sun += $_; }
        my $ek = 0;
        if ($sun > 1) { $ek = 2 * $mch / ($sun * ($sun - 1)); }
        $k += $ek;
    }
    return $k;
}

# 计算ss值
sub cal_ss {
    my ($fr) = @_;
    my $ss = 0;
    my @sfr = ();
    foreach (@$fr) {
        @sfr = split(/:/, $_);
        my $i = 0;
        foreach (@sfr) { if ($_ > 0) { $i++; } }
        $ss += $i - 1;
    }
    return $ss;
}

# 计算sn值
sub cal_sn {
    my ($fr) = @_;
    my $n = 0;
    my $fn = 0;
    my @sfr = ();
    foreach (@$fr) {
        @sfr = split(/:/, $_);
        my $i = 0;
        my $j = 0;
        foreach (@sfr) { if ($_ > 0) { $i++; } }
        foreach (@sfr) { $j += $_; }
        if ($i - 1 > 0) { $n++; $fn += $j; }
    }
    my $sn = 0;
    if ($n > 0) { $sn = int($fn / $n + 0.5); }
    return $sn;
}

sub cal_dxy 
      { my ($p1, $p2, $ln) = @_;
        my $dxy = 0;
        my @sfr1 = ();
	my @sfr2 = ();
	my ($i, $j);
	my $sn=scalar(@$p1);
	for $i(0..$sn-1)
	      { @sfr1=split(/:/,$p1->[$i]);               
                @sfr2=split(/:/,$p2->[$i]);
                my $mch=0;
                my $an=scalar(@sfr1);
                my $am1=0;
                my $am2=0;
                foreach(@sfr1){$am1 += $_;}
                foreach(@sfr2){$am2 += $_;}
                my $am=$am1*$am2;
                for $i(0..$an-1)
                      {for $j(0..$an-1)
                             { if($i==$j){next;}
                               else      {$mch += $sfr1[$i]*$sfr2[$j];}
                             }
                      }
                if($am>0){$dxy += $mch/$am;}
               }
        $dxy = $dxy/$ln;
        return $dxy;
       }             
                             


# calcute theta-W 
sub cal_theta 
      { my ($ss, $sn, $ln) = @_;   # $ss: number of segregating sites  $sn: sample size  $ln: sequence length 
        my $a1=0;
        my $i=0;
        my $theta=0;
        if($ss>0 && $sn>0)
          {    for $i(1..$sn-1) {$a1 += 1/$i;}
               $theta = $ss/$a1/$ln;
          }     
        return $theta;
       }

# calcute tajima-D 
sub cal_tajimaD
     {  my ($k, $ss, $sn) = @_;
        my ($a1, $b1, $c1, $e1) = (0, 0, 0, 0);
        my ($a2, $b2, $c2, $e2) = (0, 0, 0, 0);
        my ($i, $j);
        my $td=0;
        if($ss>1 && $sn>0)
          { for $i(1..$sn-1) {$a1 += 1/$i;}
            $b1 = ($sn+1)/(3*($sn-1));
            $c1 = $b1-1/$a1;
            $e1 = $c1/$a1;
            for $i(1..$sn-1) {$a2 += 1/($i*$i);}
            $b2 = 2*($sn*$sn+$sn+3)/(9*$sn*($sn-1));
            $c2 = $b2-($sn+2)/($a1*$sn)+$a2/($a1*$a1);
            $e2 = $c2/($a1*$a1+$a2);
            $td=($k-$ss/$a1)/sqrt($e1*$ss+$e2*$ss*($ss-1));
           }
        return $td;
      }         
        
        
       
