#!/usr/bin/perl -w
use strict;

my ($in, $gff, $out) = @ARGV;

# 1. 读取区间文件并构建按染色体分组的排序区间列表
open my $in_fh, '<', $in or die "Can't open $in: $!";
my %intervals;
while (<$in_fh>) {
    chomp;
    my @tn = split /\t/;
    die "Invalid SSW line: $_" if @tn < 3;
    push @{$intervals{$tn[0]}}, [$tn[1], $tn[2]];
}
close $in_fh;

# 2. 对每个染色体的区间排序并合并重叠
foreach my $chr (keys %intervals) {
    # 按起始位置排序
    @{$intervals{$chr}} = sort { $a->[0] <=> $b->[0] } @{$intervals{$chr}};
    
    # 合并重叠区间
    my @merged;
    my $current = $intervals{$chr}[0];
    for my $i (1..$#{$intervals{$chr}}) {
        my $next = $intervals{$chr}[$i];
        if ($current->[1] >= $next->[0]) {
            $current->[1] = $next->[1] if $next->[1] > $current->[1]; # 扩展区间
        } else {
            push @merged, $current;
            $current = $next;
        }
    }
    push @merged, $current;
    $intervals{$chr} = \@merged;
}

# 3. 处理GFF文件并查询重叠
open my $gff_fh, '<', $gff or die "Can't open $gff: $!";
open my $out_fh, '>', $out or die "Can't open $out: $!";
my %results;

while (<$gff_fh>) {
    chomp;
    next if /^#/;  # 跳过注释行
    my @tm = split /\t/;
    die "Invalid GFF line: $_" if @tm < 4;
    
    my $chr = $tm[0];
    my $feat_start = $tm[1] - 1;  # GFF转BED坐标（1-based → 0-based）
    my $feat_end   = $tm[2] - 1;
    my $attr = $tm[3];
    
    next unless exists $intervals{$chr};
    
    # 4. 二分查找重叠区间
    my $arr_ref = $intervals{$chr};
    my $low = 0;
    my $high = $#$arr_ref;
    while ($low <= $high) {
        my $mid = int(($low + $high) / 2);
        my ($int_start, $int_end) = @{$arr_ref->[$mid]};
        
        # 判断重叠：特征区间与目标区间相交
        if ($feat_start <= $int_end && $feat_end >= $int_start) {
            # 找到重叠后，向两侧继续查找所有重叠区间
            my $i = $mid;

            # 向左扫描
            while ($i >= 0) {
                my ($left_start, $left_end) = @{$arr_ref->[$i]};
                last if $feat_end < $left_start;
                if ($feat_start <= $left_end && $feat_end >= $left_start) {
                    my $key = join "\t", $chr, $left_start, $left_end, $attr;
                    $results{$key} = 1;
                }
                $i--;
            }

            # 向右扫描
            $i = $mid + 1;
            while ($i <= $#$arr_ref) {
                my ($right_start, $right_end) = @{$arr_ref->[$i]};
                last if $feat_start > $right_end;
                if ($feat_start <= $right_end && $feat_end >= $right_start) {
                    my $key = join "\t", $chr, $right_start, $right_end, $attr;
                    $results{$key} = 1;
                }
                $i++;
            }

            last;  # 已处理完所有重叠，退出二分查找循环
        }

        # 调整搜索范围
        elsif ($feat_end < $int_start) {
            $high = $mid - 1;
        }
        else {
            $low = $mid + 1;
        }
    }
}

# 5. 按染色体+坐标排序输出
for my $key (sort {
    my @a = split /\t/, $a;
    my @b = split /\t/, $b;
    $a[0] cmp $b[0] || $a[1] <=> $b[1] || $a[2] <=> $b[2]
} keys %results) {
    print $out_fh "$key\n";
}

close $_ for ($gff_fh, $out_fh);
