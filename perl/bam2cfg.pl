#!/usr/bin/env perl
#Create a BreakDancer configuration file from a set of bam files

use strict;
use warnings;
use Getopt::Std;
use Statistics::Descriptive;
use GD::Graph::histogram;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
#use lib '/gscuser/kchen/1000genomes/analysis/scripts/';
use AlnParser;

$| = 1; # enable AUTOFLUSH modE

my %opts = (q=>35, n=>10000, v=>1, c=>4, b=>50, s=>50);
getopts('q:n:c:b:p:s:hmf:gCv:', \%opts);
die("
Usage:   bam2cfg.pl <bam files>
Options:
         -q INT    Minimum mapping quality [$opts{q}]
         -m        Using mapping quality instead of alternative mapping quality
         -s        Minimal mean insert size [$opts{s}]
         -C        Change default system from Illumina to SOLiD
         -c FLOAT  Cutoff in unit of standard deviation [$opts{c}]
         -n INT    Number of observation required to estimate mean and s.d. insert size [$opts{n}]
         -v FLOAT  Cutoff on coefficients of variation [$opts{v}]
         -f STRING A two column tab-delimited text file (RG, LIB) specify the RG=>LIB mapping, useful when BAM header is incomplete
	 -b INT	   Number of bins in the histogram [$opts{b}] 
         -g        Output mapping flag distribution
         -h        Plot insert size histogram for each BAM library
\n
") unless (@ARGV);

my $AP=new AlnParser(platform=>$opts{C});

my %cRGlib; my %clibs;
if($opts{f}){
  open(RGLIB,"<$opts{f}") || die "unable to open $opts{f}\n";
  while(<RGLIB>){
    chomp;
    my ($rg,$lib)=split;
    $clibs{$lib}=1;
    $cRGlib{$rg}=$lib;
  }
}

foreach my $fbam(@ARGV){
  my %RGlib;
  my %RGplatform;
  my %libs;
  my %insert_stat;
  my %readlen_stat;
  my %libpos;
  my %flagHgram;
  my $recordcounter=0;
  my $expected_max;
  if(defined $opts{f}){
    %RGlib=%cRGlib;
    %libs=%clibs;
  }

  my $ppos = 0;
  my $lastchr = "";
  stderr_log('Processing bam: ', $fbam);
  my $samtools_pid = open(my $bam, "samtools view -h $fbam |")
    || die "unable to open $fbam\n";

  while(<$bam>){
    chomp;
    if(/^\@RG/){  #getting RG=>LIB mapping from the bam header
      my ($id)=($_=~/ID\:(\S+)/);
      my ($lib)=($_=~/LB\:(\S+)/);
      my ($platform)=($_=~/PL\:(\S+)/);
      my ($sample)=($_=~/SM\:(\S+)/);
      my ($insertsize)=($_=~/PI\:(\d+)/);
      #if(defined $insertsize && $insertsize>0){
	#$lib=$sample . '_'. $lib;
	$libs{$lib}=1;
	$RGlib{$id}=$lib;
        $RGplatform{$id}=$platform;
      #}
    }
    else{
      next if(/^\@/);
      my @libas=keys %libs;
      my @selected_libs=keys %insert_stat;
      if($#libas<0){ 
	if($#selected_libs>=0){
          stderr_log('selected_libs is : ', scalar @selected_libs);
	  last;
	}
	else{
	  $libs{'NA'}=1;
	  $RGlib{'NA'}='NA';
	  $RGplatform{'NA'}=($opts{C})?'solid':'illumina';
	}
      }
      if(!defined $expected_max || $expected_max<=0){
	$expected_max=3*($#libas+1)*$opts{n};
      }

      if ($recordcounter > $expected_max) {
         stderr_log('$recordcounter > $expected_max: ', $recordcounter, " > ", $expected_max);
         last;
      }

      my $t=$AP->in($_,'sam',\%RGplatform,$opts{m});
      $lastchr = $t->{chr} unless defined $lastchr;
      $ppos = 0 if $t->{chr} ne $lastchr;
      $lastchr = $t->{chr};
      die "Please sort bam by position\n" if($t->{pos}<$ppos);
      $ppos=$t->{pos};
      my $lib=defined($t->{readgroup})?$RGlib{$t->{readgroup}}:'NA';  #when multiple libraries are in a BAM file
      next unless(defined $lib && $libs{$lib});
      $readlen_stat{$lib}=Statistics::Descriptive::Full->new() if(!defined $readlen_stat{$lib});
      $readlen_stat{$lib}->add_data($t->{readlen});
      next if ($t->{qual}<=$opts{q});  #skip low quality mapped reads
      $recordcounter++;
      $libpos{$lib}++;
      if(defined $t->{readgroup}){
	$flagHgram{$t->{readgroup}}{$t->{flag}}++;
	$flagHgram{$t->{readgroup}}{all}++;
      }

      my $nreads=(defined $insert_stat{$lib})?$insert_stat{$lib}->count():1;
      if($nreads/$libpos{$lib}<1e-4){  #single-end lane
	delete $libs{$lib};
	delete $insert_stat{$lib};
      }
      next unless(($t->{flag}==18 || $t->{flag}==20) && $t->{dist}>=0);

      $insert_stat{$lib}=Statistics::Descriptive::Full->new() if(!defined $insert_stat{$lib});
      $insert_stat{$lib}->add_data($t->{dist});
      if($insert_stat{$lib}->count()>$opts{n}){
	delete $libs{$lib};
      }
    }
  }

  stderr_log('Closing BAM file');
  close_samtools($bam, $samtools_pid);

  my %stdms;
  my %stdps;

  foreach my $lib(keys %insert_stat){
    my $readlen=$readlen_stat{$lib}->mean();
    my @isize=$insert_stat{$lib}->get_data();
    my $mean=$insert_stat{$lib}->mean();
    my $std=$insert_stat{$lib}->standard_deviation();

    delete $insert_stat{$lib};
    my $insertsize=Statistics::Descriptive::Full->new();
    foreach my $x(@isize){
      next if($x>$mean+5*$std);
      $insertsize->add_data($x);
    }

    $mean=$insertsize->mean();
    $std=$insertsize->standard_deviation();
    next if($mean<$opts{s});
    my $cv=$std/$mean;
    if($cv>=$opts{v}){
      print STDERR "Coefficient of variation $cv in library $lib is larger than the cutoff $opts{v}, poor quality data, excluding from further analysis.\n";
      next;
    }

    my $num=$insertsize->count();
    next if($num<100);

    my ($stdm,$stdp)=(0,0);
    my ($nstdm,$nstdp)=(0,0);
    foreach my $x($insertsize->get_data()){
      if($x>$mean){
	$stdp+=($x-$mean)**2;
	$nstdp++;
      }
      else{
	$stdm+=($x-$mean)**2;
	$nstdm++;
      }
    }
    $stdm=sqrt($stdm/($nstdm-1));
    $stdp=sqrt($stdp/($nstdp-1));

    $stdms{$lib}=$stdm;
    $stdps{$lib}=$stdp;
    $insert_stat{$lib}=$insertsize;
  }

  foreach my $rg(keys %RGlib){
    my $lib=$RGlib{$rg};
    my $platform=$RGplatform{$rg} || 'illumina';  #default illumina
    next unless($insert_stat{$lib});
    my $readlen=$readlen_stat{$lib}->mean();
    my $mean=$insert_stat{$lib}->mean();
    my $std=$insert_stat{$lib}->standard_deviation();
    my $num=$insert_stat{$lib}->count();

    my $upper=$mean+$opts{c}*$stdps{$lib} if(defined $opts{c});
    my $lower=$mean-$opts{c}*$stdms{$lib} if(defined $opts{c});
    $lower=0 if(defined $lower && $lower<0);

    printf "readgroup\:%s\tplatform:%s\tmap\:%s\treadlen\:%.2f\tlib\:%s\tnum:%d",$rg,$platform,$fbam,$readlen,$lib,$num;
    printf "\tlower\:%.2f\tupper\:%.2f",$lower,$upper if(defined $upper && defined $lower);
    printf "\tmean\:%.2f\tstd\:%.2f",$mean,$std;


    # compute the normality
    my @data=$insert_stat{$lib}->get_data();
    my $n_data = @data;
    @data = sort {$a <=> $b} @data;
    my $n_data_ = @data;
    my $p_value = ShapiroWilk(\@data);
    if($p_value>0) {
      my $log_p_value = log($p_value)/log(10);
      printf "\tSWnormality\:%.2f",$log_p_value;
    }
    elsif($p_value == -1) { 
      printf "\tSWnormality\:data not qualified -1"; 
    } 
    elsif($p_value == -2.1){
      printf "\tSWnormality\:data not qualified -2.1";
    }
    elsif($p_value == -2.2){
      printf "\tSWnormality\:data not qualified -2.2";
    }
    elsif($p_value == -2.3){
      printf "\tSWnormality\:data not qualified -2.3";
    }
    elsif($p_value == 0) { 
      printf "\tSWnormality\:minus infinity";
    }



    if($opts{g}){
      printf "\tflag:";
      foreach my $f(sort keys %{$flagHgram{$rg}}){
	next if($f eq 'all');
	printf "%d(%.2f%%)",$f,($flagHgram{$rg}{$f} || 0)*100/$flagHgram{$rg}{all};
      }
      printf "%d",$flagHgram{$rg}{all};
    }
    printf "\texe:samtools view\n";
  }

  if($opts{h}){  # plot insert size histogram for each library
    foreach my $lib(keys %insert_stat){
      my $graph = new GD::Graph::histogram(1000,600);
      my $library="$fbam.$lib";
      $graph->set(
		  x_label         => 'Insert Size (bp)',
		  y_label         => 'Count',
		  title           => $library,
		  x_labels_vertical => 1,
		  bar_spacing     => 0,
		  shadow_depth    => 1,
		  shadowclr       => 'dred',
		  transparent     => 0,
		  histogram_bins   => $opts{b},
		 ) or warn $graph->error;
      my @data=$insert_stat{$lib}->get_data();
      my $gd = $graph->plot(\@data) or die $graph->error;

      $library=~s/.*\///g;
      my $imagefile="$library.insertsize_histogram.png";
      open(IMG, ">$imagefile") or die $!;
      binmode IMG;
      print IMG $gd->png;

      my $datafile="$library.insertsize_histogram";
      open(OUT,">$datafile");
      foreach my $x(@data){
	print OUT "$x\n";
      }
      close(OUT);
    }
  }
}

###################################################### Shapiro Wilk Normality test (including functions ppnd, alnorm, poly, min, sign and asin) #########################################################
sub ShapiroWilk {

# read input
my ($x) = @_;

my $n=@$x;
my $w = 0;
my $ifault = 2;
my $init = 0;
my $n2 = $n/2;
my @a;
my $n1 = $n;

my $upper = 1;
my @c1 = (0.0, 0.221157, -0.147981, -2.07119, 4.434685, -2.706056);
my @c2 = (0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633);
my @c3 = (0.5440, -0.39978, 0.025054, -0.6714*(10**(-3)));
my @c4 = (1.3822, -0.77857, 0.062767, -0.0020322);
my @c5 = (-1.5861, -0.31082, -0.083751, 0.0038915);
my @c6 = (-0.4803, -0.082676, 0.0030302);
my @c7 = (0.164, 0.533);
my @c8 = (0.1736, 0.315);
my @c9 = (0.256, -0.00635);
my @g  = (-2.273, 0.459);
my $z90 = 1.2816;
my $z95 = 1.6449;
my $z99 = 2.3263;
my $zm = 1.7509;
my $zss = 0.56268;
my $bf1 = 0.8378; 
my $xx90 = 0.556;
my $xx95 = 0.622;
my $zero = 0.0;
my $one = 1.0;
my $two = 2.0;
my $three = 3.0;
my $sqrth = 0.70711;
my $qtr = 0.25;
my $th = 0.375;
my $small = 1*(10**(-19));
my $pi6 = 1.909859;
my $stqr = 1.047198;

my $summ2;
my $ssumm2;
my $fac;
my $rsn;
my $an;
my $an25;
my $a1;
my $a2;
my $delta;
my $range;
my $sa;
my $sx;
my $ssx;
my $ssa;
my $sax;
my $asa;
my $xsx;
my $ssassx;
my $w1;
my $y;
my $xx;
my $xi;
my $gamma;
my $m;
my $s;
my $ld;
my $bf;
my $z90f;
my $z95f;
my $z99f;
my $zfm;
my $zsd;
my $zbar;
my $ncens;
my $nn2;
my $i;
my $i1;
my $j;

my $pw  =  $one;
if($w >= $zero) {$w = $one;}
$an = $n;
$ifault = 3;
$nn2 = $n/2;
if($n2 < $nn2) {return (-1);}
$ifault = 1;
if($n < 3) {return (-1);}

# If INIT is false, calculates coefficients for the test

if(!$init){
  if($n == 3) {
    $a[1-1] = $sqrth;}	
  else {
    $an25 = $an + $qtr;
    $summ2 = $zero;
    for(my $i = 1; $i <= $n2; $i++){
      ($a[$i-1], $ifault) = ppnd(($i - $th)/$an25);
      my $a_val = $a[$i-1];
      $summ2 = $summ2 + $a[$i-1] ** 2;
    }
    $summ2 = $summ2 * $two;
    $ssumm2 = $summ2 ** 0.5;
    $rsn = $one / ($an ** 0.5);
    my $poly_result = poly_(\@c1, 6, $rsn);
    $a1 = $poly_result - $a[1-1] / $ssumm2;
    
# Normalize coefficients
    
    if($n > 5) {
      $i1 = 3;
      $a2 = -$a[2-1]/$ssumm2 + poly_(\@c2, 6, $rsn);
      my $a_val1_new = $a[1-1];
      my $a_val2_new = $a[2-1];
      my $temp1 = $summ2 - $two * $a[1-1] ** 2 - $two * $a[2-1] ** 2;
      my $temp2 = $one - $two * $a1 ** 2 - $two * $a2 ** 2;
      $fac = (($summ2 - $two * $a[1-1] ** 2 - $two * $a[2-1] ** 2)/($one - $two * $a1 ** 2 - $two * $a2 ** 2)) ** 0.5;
      $a[1-1] = $a1;
      $a[2-1] = $a2;
    }
    else {
      $i1 = 2;
      $fac = (($summ2 - $two * $a[1-1] ** 2)/ ($one - $two * $a1 ** 2)) ** 0.5;
      $a[1-1] = $a1;
    }
    for($i = $i1; $i <= $nn2; $i++){
      $a[$i-1] = -$a[$i-1]/$fac;
      my $a_val = $a[$i-1];
    }
    for($i = 1; $i <= $nn2; $i++){
      my $a_value = $a[$i-1];
    }
  }
  $init = 1;
}
if($n1 < 3) {return (-1);}
$ncens = $n - $n1;
$ifault = 4;
if($ncens < 0 || ($ncens > 0 && $n < 20)) {return (-1);}
$ifault = 5;
$delta = $ncens/$an;
if($delta > 0.8) {return (-2.1);}

# If W input as negative, calculate significance level of -W

if($w < $zero) {
  $w1 = $one + $w;
  $ifault = 0;
}
else{
# Check for zero range
  $ifault = 6;
  $range = @$x[$n1-1] - @$x[1-1];
  if($range < $small) {return (-2.2);}

# Check for correct sort order on range - scaled X

  $ifault = 7;
  $xx = @$x[1-1]/$range;
  $sx = $xx;
  $sa = -$a[1-1];
  $j = $n - 1;
  for($i = 2; $i <= $n1; $i++){
    $xi = @$x[$i-1]/$range;
    if($xx-$xi > $small){
      return (-2.3);
    }
    $sx = $sx + $xi;    
    if($i != $j) {
      my $sign_value = sign($i - $j);
      my $min_value = $a[min($i, $j)-1];
      my $min_ = min($i, $j);
      $sa = $sa + sign($i - $j) * $a[min($i, $j)-1];
    } 
    $xx = $xi;
    $j = $j - 1;
  }
  $ifault = 0;
  if($n > 5000) {$ifault = 2};

# Calculate W statistic as squared correlation between data and coefficients
  $sa = $sa/$n1;
  $sx = $sx/$n1;
  $ssa = $zero;
  $ssx = $zero;
  $sax = $zero;
  $j = $n;
  for($i = 1; $i <= $n1; $i ++){
    if($i != $j) {
      my $sign_result = sign($i - $j);
      my $min_result = min($i, $j);
      $asa = sign($i - $j) * $a[min($i, $j)-1] - $sa;
    } 
    else {
      $asa = -$sa;
    }
    $xsx = @$x[$i-1]/$range - $sx;
    $ssa = $ssa + $asa * $asa;
    $ssx = $ssx + $xsx * $xsx;
    $sax = $sax + $asa * $xsx;
    $j = $j - 1;
  }

# W1 equals (1-W) claculated to avoid excessive rounding error for W very near 1 (a potential problem in very large samples)
  $ssassx = ($ssa * $ssx) ** 0.5;
  $w1 = ($ssassx - $sax) * ($ssassx + $sax)/($ssa * $ssx);
}

$w = $one - $w1;

# Calculate significance level for W (exact for N=3)

if($n == 3){
  $pw = $pi6 * (asin($w**0.5) - $stqr); 
  return ($pw);
}
$y = log($w1);
$xx = log($an);
$m = $zero;
$s = $one;
if($n <= 11){
  $gamma = poly_(\@g, 2, $an); 
  if($y >= $gamma){
    $pw = $small;
    return ($pw);
  }
  $y = -log($gamma - $y);
  $m = poly_(\@c3, 4, $an);
  $s = exp(poly_(\@c4, 4, $an))
}
else{
  $m = poly_(\@c5, 4, $xx);
  $s = exp(poly_(\@c6, 3, $xx));
}
if($ncens > 0){
  
# Censoring by proportion NCENS/N.  Calculate mean and sd of normal equivalent deviate of W.
  
  $ld = -log($delta);
  $bf = $one + $xx * $bf1;
  $z90f = $z90 + $bf * poly_(\@c7, 2, $xx90 ** $xx) ** $ld;
  $z95f = $z95 + $bf * poly_(\@c8, 2, $xx95 ** $xx) ** $ld;
  $z99f = $z99 + $bf * poly_(\@c9, 2, $xx) ** $ld;
  
# Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get pseudo-mean and pseudo-sd of z as the slope and intercept
  
  $zfm = ($z90f + $z95f + $z99f)/$three;
  $zsd = ($z90*($z90f-$zfm)+$z95*($z95f-$zfm)+$z99*($z99f-$zfm))/$zss;
  $zbar = $zfm - $zsd * $zm;
  $m = $m + $zbar * $s;
  $s = $s * $zsd;
}
$pw = alnorm(($y - $m)/$s, $upper);
return ($pw);

}

############################################### ppnd function############################################################
sub ppnd {

# read input
my ($p) = @_;

# define output
my $ifault;
my $normal_dev;

# Local variables

my $zero = 0.0;
my $one = 1.0;
my $half = 0.5;
my $split1 = 0.425;
my $split2 = 5.0;
my $const1 = 0.180625;
my $const2 = 1.6;
my $q;
my $r;

# Coefficients for P close to 0.5

my $a0 = 3.3871327179;
my $a1 = 5.0434271938*10;
my $a2 = 1.5929113202*10**2;
my $a3 = 5.9109374720*10;
my $b1 = 1.7895169469*10;
my $b2 = 7.8757757664*10;
my $b3 = 6.7187563600*10;
# HASH SUM AB          32.3184577772

# Coefficients for P not close to 0, 0.5 or 1.

my $c0 = 1.4234372777;
my $c1 = 2.7568153900;
my $c2 = 1.3067284816;
my $c3 = 1.7023821103*(10**(-1));
my $d1 = 7.3700164250*(10**(-1));
my $d2 = 1.2021132975*(10**(-1));
# HASH SUM CD          15.7614929821

# Coefficients for P near 0 or 1.

my $e0 = 6.6579051150;
my $e1 = 3.0812263860;
my $e2 = 4.2868294337*(10**(-1));
my $e3 = 1.7337203997*(10**(-2));
my $f1 = 2.4197894225*(10**(-1));
my $f2 = 1.2258202635*(10**(-2));
# HASH SUM EF          19.4052910204

$ifault = 0;
$q = $p - $half;
if(abs($q) <= $split1) { 
  $r = $const1 - $q * $q;
  $normal_dev = $q * ((($a3 * $r + $a2) * $r + $a1) * $r + $a0) / ((($b3 * $r + $b2) * $r + $b1) * $r + $one);
}
else{
  if($q < $zero) {
    $r = $p;}
  else {
    $r = $one - $p;
  }
  if($r <= $zero) {
    $ifault = 1;
    $normal_dev = $zero;
    return ($normal_dev, $ifault);
  }
  $r = (-log($r))**0.5;
  if($r <= $split2) {
    $r = $r - $const2;
    $normal_dev = ((($c3 * $r + $c2) * $r + $c1) * $r + $c0) / (($d2 * $r + $d1) * $r + $one);
  }
  else{
    $r = $r - $split2;
    $normal_dev = ((($e3 * $r + $e2) * $r + $e1) * $r + $e0) / (($f2 * $r + $f1) * $r + $one);
  }
  if($q < $zero) {$normal_dev = - $normal_dev;}
  return ($normal_dev, $ifault);
}
}


################################ alnorm function #####################################################################
#FUNCTION alnorm(x, upper) RESULT(fn_val)

#  Evaluates the tail area of the standardised normal curve from x to infinity if upper is .true. or
#  from minus infinity to x if upper is .false.

sub alnorm {

# get the inputs
my ($x, $upper) = @_;

# define the output
my $fn_val;

my $zero = 0.0;
my $one = 1.0;
my $half = 0.5;
my $con = 1.28;
my $z;
my $y;
my $up;

#!*** machine dependent constants
my $ltone = 7.0;
my $utzero = 18.66;

my $p = 0.398942280444;
my $q = 0.39990348504;
my $r = 0.398942280385;
my $a1 = 5.75885480458;
my $a2 = 2.62433121679; 
my $a3 = 5.92885724438;
my $b1 = -29.8213557807; 
my $b2 = 48.6959930692;
my $c1 = -3.8052*(10**(-8)); 
my $c2 = 3.98064794*(10**(-4));
my $c3 = -0.151679116635; 
my $c4 = 4.8385912808;
my $c5 = 0.742380924027; 
my $c6 = 3.99019417011;
my $d1 = 1.00000615302; 
my $d2 = 1.98615381364;
my $d3 = 5.29330324926;
my $d4 = -15.1508972451;
my $d5 = 30.789933034;

$up = $upper;
$z = $x;
if($z <  $zero){
  if($up == 0){$up = 1;}
  else {$up = 0;}
  $z = -$z;
}
if($z <= $ltone || $up && $z <= $utzero){
  $y = $half*$z*$z;
  if($z > $con) {$fn_val = $r*exp(-$y)/($z+$c1+$d1/($z+$c2+$d2/($z+$c3+$d3/($z+$c4+$d4/($z+$c5+$d5/($z+$c6))))));}
  else {$fn_val = $half - $z*($p - $q*$y/($y+$a1+$b1/($y+$a2+$b2/($y+$a3))));}
}
else{
  $fn_val = $zero;
}
if(!$up) {$fn_val = $one - $fn_val};
return $fn_val;
}

##################################### poly function #########################################################

# FUNCTION poly(c, nord, x) RESULT(fn_val)

# Calculates the algebraic polynomial of order nored-1 with array of coefficients c.  Zero order coefficient is c(1)

sub poly_ {

my ($c, $nord, $x) = @_;

my $fn_val;

my $i;
my $j;
my $n2;
my $p;

$fn_val = @$c[1-1];
if($nord == 1) {return $fn_val;}
$p = $x*@$c[$nord-1];
if($nord == 2) {
  $fn_val = $fn_val + $p; 
}
else {
  $n2 = $nord - 2;
  $j = $n2 + 1;
  for($i = 1; $i <= $n2; $i++){
    $p = ($p + @$c[$j-1])*$x;
    $j = $j - 1;
  }
  $fn_val = $fn_val + $p;
}

return $fn_val;
}

###################################### min function ############################################

# FUNCTION min(a,b) RESULT min of a, b

sub min {

my ($a, $b) = @_;

if($a <= $b) { return $a; }
else { return $b; }
}

###################################### sign function ############################################

# FUNCTION sign(a,b) RESULT sign of a, b

sub sign {

my ($a) = @_;

if($a >= 0) { return 1; }
else { return -1; }
}

###################################### asin function ############################################

# FUNCTION asin(a) RESULT arcsin of a

sub asin {

my ($a) = @_;
if($a > 1 || $a < -1) { 
  print "error: abs sin value > 1";
  return 0; }
else {
  my $b = (1-$a**2)**0.5;
  my $fn_val = atan2($a, (1-$a**2)**0.5);
  my $fn_test = (1-$a**2)**0.5;
  return $fn_val;
}
}


sub close_samtools {
    my ($fh, $pid) = @_;

    #kill samtools view nicely
    stderr_log("Send TERM signal for $pid");
    kill('TERM', $pid);
    sleep 2;

    if (kill(0 , $pid)) {
        stderr_log("samtools pid process $pid is still there...");
        stderr_log("invoking kill -9 on $pid ...");
        kill(9, $pid) or die "[$0] [err] Trouble killing samtools pid: $pid !\n";
    }
    else {
        stderr_log("samtools pid process $pid is now gone");
    }

    stderr_log('Closing samtools process : ', $pid);
    close $fh;

    return 1;
}

sub stderr_log {
    my @msg = @_;
    print STDERR '[', scalar localtime, " $0] ", @msg, "\n";
}
