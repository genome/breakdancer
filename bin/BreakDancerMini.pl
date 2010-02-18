#!/gsc/bin/perl
# This is the BreakDancer script that focuses on small indel detection
# Indels as short as 1 std and as long as (maximum insert size - mean insert size) can be detected
# Copyright (C) 2009 Washington University in St. Louis

use strict;
use warnings;
use Getopt::Std;
use Statistics::Descriptive;
use Math::CDF;

# load Breakdancer Specific libraries from the distribution dir
use FindBin qw($Bin);
use lib "$FindBin::Bin/../lib";

use AlnParser;

my $version="BreakDancerMini-0.0.1r57";
my %opts = (i=>200, m=>100000, q=>9, s=>20, r=>2, b=>100, p=>0.001, v=>0, e=>0, a=>0.0001);
my %opts1;
getopts('o:s:p:m:q:r:b:v:e:tfl:g:', \%opts1);
die("
Usage:   BreakDancerMini.pl <analysis_config.lst>
Options:
         -o STRING Operate on a single chromosome [all chromosome]
         -p INT    Cutoff in the Kolmogorov-Smirnov(KS) P value [$opts{p}]
         -a STRING A priori probability of in an SV [$opts{a}]
         -l STRING Use only the specified library
         -g STRING Output Lib Q values in position-sorted regions in the file [chr start stop]
         -f        Use Fisher's method to combine P values from multiple library
         -m INT    Maximum SV size [$opts{m}]
         -q INT    Minimum alternative mapping quality [$opts{q}]
         -w INT    SV flanking window size [mean + std insert size - mean read length]
         -s INT    Step size [$opts{s}]
         -e INT    SV score method (0: average, 1: max) of the flanking region scores [$opts{e}]
         -r INT    Minimum number of read pairs required to establish a connection [$opts{r}]
         -b INT    Number of regions registered before calling SVs [$opts{b}]
         -v INT    Verbose mode 0-2 [$opts{v}]
Version: $version\n
") unless (@ARGV);

my $options='';
foreach my $opt(keys %opts1){
  $options.=$opt.$opts1{$opt};
  $opts{$opt}=$opts1{$opt};
}

my %SVtype=(
	    '2'=>'DEL',  #deletion
	    '3'=>'INS',  #insertion
	   );

my $RG=&ReadRegions($opts{g}) if(defined $opts{g});
my $AP=new AlnParser();
my %exes;
my %format;
my %mean_insertsize;
my %std_insertsize;
my %readlens;
my $max_readlen=0;
my %winsize;
my $minwinsize=1e10;
my $maxwinsize=0;
my %fmaps;
my %libmaps;
my %readgroup_library;
my %uppercutoff;
my $PI=3.1415926;
my $LZERO=-99;
my $ZERO=exp($LZERO);
my $bufferend=0;
open(CONFIG,"<$ARGV[0]") || die "unable to open $ARGV[0]\n";
while(<CONFIG>){
  next unless (/\S+/);
  chomp;
  my $fh;
  my ($readgroup)=($_=~/group\w*\:(\S+)\t*/);
  my ($fmap)=($_=~/map\w*\:(\S+)\t*/);

  my ($mean)=($_=~/mean\w*\:(\S+)\t*/);
  my ($std)=($_=~/std\w*\:(\S+)\t*/);
  my ($readlen)=($_=~/readlen\w*\:(\S+)\t*/);
  my ($upper)=($_=~/upp\w*\:(\S+)\t*/);
  my ($lib)=($_=~/lib\w*\:(\S+)\t*/);
  ($lib)=($_=~/samp\w*\:(\S+)\t*/) if(!defined $lib);
  die "Please include a library name for each row in the configure file.\n" if(!defined $lib);
  my ($exe)=($_=~/exe\w*\:(.+)\t*/);

  $readgroup=$lib if(!defined $readgroup);
  $readgroup_library{$readgroup}=$lib;
  $libmaps{$lib}=$fmap;
  $winsize{$lib}=$opts{w} || $mean+$std-$readlen;
  $fmaps{$fmap}=$lib;

  if(defined $mean && defined $std && !defined $upper){
    $upper=$mean+$std*5;
  }
  $mean_insertsize{$lib}=$mean;
  $std_insertsize{$lib}=$std;
  $uppercutoff{$lib}=$upper;
  $readlens{$lib}=$readlen;

  if(!defined $exes{$fmap}){
    $exes{$fmap}=$exe;
  }
  else{
    die "Please use identical exe commands to open the same map file.\n" if($exes{$fmap} ne $exe);
  }

  $minwinsize=($minwinsize>$winsize{$lib})?$winsize{$lib}:$minwinsize;
  $maxwinsize=($maxwinsize<$winsize{$lib})?$winsize{$lib}:$maxwinsize;
  $bufferend=($bufferend<$winsize{$lib})?$winsize{$lib}:$bufferend;
  $max_readlen=($max_readlen<$readlen)?$readlen:$max_readlen;
}
close(CONFIG);

my @maps=keys %fmaps;

my %pmfs;
my %n_pmfreads;
my %nreads;
my @format;
my $ref_len=0;
for(my $i=0;$i<=$#maps;$i++){
  my $fh;

  my $reflen=0;
  my $exe=$exes{$maps[$i]};
  if($exe){
    open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
    push @format,($exe=~/samtools/)?'sam':'maq';
  }
  else{
    open($fh,"<$maps[$i]") || die "unable to open $maps[$i]\n";
    push @format,'maq';
  }
  my $p_pos=0;
  my $p_chr='';

  while(<$fh>){
    my $t=$AP->in($_,$format[$i]);
    $reflen+=$t->{pos}-$p_pos if($t->{chr} eq $p_chr);
    $p_pos=$t->{pos};
    $p_chr=$t->{chr};
    my $lib=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[$i]};  #when multiple libraries are in a BAM file
    next unless(defined $lib);
    next if(defined $opts{l} && $opts{l} ne $lib);
    $nreads{$lib}++;

    #using only high quality normally mapped reads to build the null emprical distribution
    next if ($t->{qual}<40 || $t->{flag}!=18 || $t->{dist}<=0);

    my $sepdist=abs($t->{dist});
    if($sepdist<=$uppercutoff{$lib}){  #only use normally mapped reads
      $pmfs{$lib}[$sepdist]++;
      $n_pmfreads{$lib}++;
    }
  }
  close($fh);
  $ref_len=($ref_len<$reflen)?$reflen:$ref_len;
}

my %CMF;  #histogram of spanning distance in each library
my %CMFn;
my $total_phy_cov=0;
my $total_seq_cov=0;
foreach my $lib(keys %pmfs){
  #compute cummulative distribution
  my @pmf=@{$pmfs{$lib}};
  $CMFn{$lib}=$n_pmfreads{$lib}+$#pmf+1;

  my $pmfmean=0;
  my $pmfstd=0;
  my @cumtmp;
  my $cumcount=0;
  for(my $dist=0;$dist<=$#pmf;$dist++){
    my $sudocount=$pmf[$dist] || 0 + 1.0;
    $cumcount+=$sudocount;
    push @cumtmp,$cumcount/$CMFn{$lib};
    $pmf[$dist]= $sudocount / $CMFn{$lib};
    $pmfmean+=$dist*$pmf[$dist];
  }
  $CMF{$lib}=\@cumtmp;

  for(my $dist=0;$dist<=$#pmf;$dist++){
    $pmfstd+=$pmf[$dist]*($dist-$pmfmean)**2;
  }

  #update mean insert size;
  $mean_insertsize{$lib}=$pmfmean;
  $std_insertsize{$lib}=sqrt($pmfstd);

  my $sequence_coverage=$nreads{$lib}*$readlens{$lib}/$ref_len;
  $total_seq_cov+=$sequence_coverage;
  my $physical_coverage=$nreads{$lib}*$mean_insertsize{$lib}/2/$ref_len;
  $total_phy_cov+=$physical_coverage;

  printf "#%s\tmean:%.3f\tstd:%.3f\tupper:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx\n", $libmaps{$lib},$mean_insertsize{$lib},$std_insertsize{$lib},$uppercutoff{$lib},$readlens{$lib},$lib,$ref_len,$sequence_coverage,$physical_coverage;
}

print "#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tAllele_frequency\tVersion\tRun_Param\n";

my ($begins, $beginc, $lasts, $lastc) = ('', -1, '', -1);
my (%regs,   # tabulate regions flanking alignment gaps
    %read,   # assign read pairs to one or more region
    %reg_name    # chromosome, start, end, and length classes of the regions
   );
my %reg_seq;
my $current_chr;  #current chromosome

#warn("-- read the mapview output\n");
my $idx_buff=0;
my $final_buff=0;
my $reg_idx=0;
my $pre_reg_idx=0;
my @FHs;
my %Idxs;
for(my $i=0;$i<=$#maps;$i++){
  my $fh;
  my $exe=$exes{$maps[$i]};
  if($exe){
    open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
  }
  else{
    open($fh,"<$maps[$i]") || die "unable to open $maps[$i]\n";
  }
  push @FHs,$fh;
  $Idxs{$i}=1;
}

my %pmaxQvalue=('+'=>-1,'-'=>-1);
my %turnon=('+'=>0,'-'=>0);
my %cluster;
my %pcluster;
my %pclusterID;
my %bestreg;
my %libpmf;
my @buffer;
my %nreads_lib;
my @cIdxs=keys %Idxs;
while(keys %Idxs){
  my ($t,$lib);
  for(my $i=0;$i<=$#cIdxs;$i++){
    my $idx=$cIdxs[$i];
    my $fh=$FHs[$idx];
    if(eof($fh)){
      delete($Idxs{$idx});
      next;
    }
    do{
      $_=<$fh>;
      chomp;
      $t=$AP->in($_,$format[$i]);
      $buffer[$idx]=$t;
      $lib=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[$i]};
    } until(((!defined $opts{o} || $t->{chr} eq $opts{o}) && (!defined $opts{l} || !defined $lib || $opts{l} eq $lib)) || eof($fh));  #analyze specific regions
  }

  my ($minchr,$minpos)=(chr(255),1e10);
  my $minidx;
  my $library;
  my $min_t;
  foreach my $i(keys %Idxs){
    $t=$buffer[$i];
    next if($t->{chr}gt$minchr || $t->{chr}eq$minchr && $t->{pos}>$minpos);
    $minchr=$t->{chr};
    $minpos=$t->{pos};
    $minidx=$i;
    $min_t=$t;
    $library=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[$i]};
  }

  if(defined $minidx){
    &Segmentation($library, $min_t) if(defined $library);
    @cIdxs=($minidx);
  }
}
$final_buff=1;
&buildConnection();

sub Segmentation{
  my ($lib,$t)=@_;

  #main analysis code
  return if ($t->{qual}<=$opts{q});  #skip low quality mapped reads
  return unless ($t->{flag}==2 || $t->{flag}==18 || $t->{flag}==130);  #analyze only reads in the FR configuration and within max insert size distance

  my $sepdist=abs($t->{dist});
  if($t->{pos} >= $bufferend){
    foreach my $ori('+','-'){

      my $logpvalue=0;
      my $ntotalreads=0;
      my @blibs=sort keys %{$libpmf{$ori}};
      my @lps;

      foreach my $blib(@blibs){
	my $Prob=0;
	if($nreads_lib{$ori}{$blib}>5){  #require at least 5 reads
	  my $maxFnDiff=-1e10;
	  my $cumcount=0;
	  foreach my $x(sort {$a<=>$b} keys %{$libpmf{$ori}{$blib}}){
	    $cumcount+=${$libpmf{$ori}{$blib}}{$x};
	    my $FnpX=$cumcount/$nreads_lib{$ori}{$blib};
	    my $FnDiff=abs($FnpX-${$CMF{$blib}}[$x]);
	    $maxFnDiff=($maxFnDiff<$FnDiff)?$FnDiff:$maxFnDiff;
	  }
	  #compute KS test probabilities based on Beta approximation
	  $Prob=&KS_TestB($maxFnDiff,$nreads_lib{$ori}{$blib});
	}
	$ntotalreads+=$nreads_lib{$ori}{$blib};
	my $lp=(1-$Prob<$ZERO)?$LZERO:log(1-$Prob);
	$logpvalue+=$lp;
	push @lps,$lp;
      }

      if($opts{g} && defined $$RG{$t->{chr}}){
	my @gpos=@{$$RG{$t->{chr}}};
	if(@gpos>0){
	  my ($rs,$re,$indelsize)=split " ",$gpos[0];
	  if($bufferend>=$rs-500){
	    if($bufferend<=$re+500){
	      #print "$bufferend:$ori:$rs\-$re:$indelsize";
	      my $offset;
	      if($bufferend<=$rs){
		$offset=$bufferend-$rs;
	      }
	      elsif($bufferend<$re){
		$offset=0;
	      }
	      else{
		$offset=$bufferend-$re;
	      }
	      printf "%d\-%d\:%d%s\t%d",$rs,$re,$indelsize,$ori,$offset;
	      for(my $tmpid=0;$tmpid<=$#lps;$tmpid++){
		my $lp=$lps[$tmpid];
		my $blib=$blibs[$tmpid];
		my $qvalue=-10*($lp/log(10));
		$qvalue=($qvalue>99)?99:int($qvalue+0.5);
		print "\t$blib\=$qvalue";
	      }
	      print "\n";
	    }
	    else{
	      shift @gpos;
	      $$RG{$t->{chr}}=\@gpos;
	    }
	  }
	}
      }

      if($opts{f}){
	#Fisher's Method
	my $fprob=Math::CDF::pchisq(-2*$logpvalue,2*($#blibs+1));
	my $fisherP=(defined $fprob)?(1-$fprob):1.0;
	$logpvalue=($fisherP>$ZERO)?log($fisherP):$LZERO;
      }
      print "$bufferend:$logpvalue\n" if($opts{v}>2);

      my $Qvalue=-10*($logpvalue/log(10));
      $Qvalue=($Qvalue>99)?99:int($Qvalue+0.5);

      if($logpvalue<log($opts{p})){
	if($pmaxQvalue{$ori}<$Qvalue && $ori eq '+' ||
	   $pmaxQvalue{$ori}<=$Qvalue && $ori eq '-'){
	  $pmaxQvalue{$ori}=$Qvalue;
	  #my $bufferstart=$bufferend-$minwinsize;
	  my $bufferstart=$bufferend-$maxwinsize;

	  $cluster{$ori}=sprintf "%s\t%d\t%d\t%.1f\t%s", $t->{chr},$bufferstart,$bufferend,$Qvalue,$ori;
	  $current_chr=$t->{chr};
	  print "$bufferend\t$Qvalue\n" if($opts{v}==2);
	  my %newreg_seq=%reg_seq;
	  $bestreg{$ori}=\%newreg_seq;
	}
	$turnon{$ori}=1;
      }
      elsif($turnon{$ori}){

	#determine if current cluster overlap with the previous one
	my ($cchr,$cps,$cpe,$cscore,$cori)=split /\t+/,$cluster{$ori};
	my ($pchr,$pps,$ppe,$pscore,$pori);
	if(defined $pcluster{$ori}){
	  ($pchr,$pps,$ppe,$pscore,$pori)=split /\t+/,$pcluster{$ori};
	}
	else{
	  $pchr=''; $ppe=0;
	}
	my $overlap=($cchr eq $pchr)?$ppe-$cps:0;

	if($overlap<=0 || $cscore>$pscore){
	  my $k;
	  if($overlap<=0){  #no overlap
	    #register new
	    $k=$reg_idx++;
	    $idx_buff++;
	  }
	  else{  #overlap and the current one has higher KS stat
	    #replace previous with current
	    $k=$pclusterID{$ori};
	  }
	  my @p;
	  foreach my $blib(keys %{$bestreg{$ori}}) {
	    foreach(@{$bestreg{$ori}{$blib}}){
	      my @s = split;
	      next unless($s[3] eq $ori);  #same orientation
	      push(@p, $_);
	    }
	  }
	  $reg_name{$k}=$cluster{$ori};
	  $regs{$k}=\@p;
	  $pcluster{$ori}=$cluster{$ori};
	  $pclusterID{$ori}=$k;
	
	  print "$k\t$cluster{$ori}\n" if($opts{v}>0);
	}

	$turnon{$ori}=0;
	$pmaxQvalue{$ori}=0;

	if($idx_buff>$opts{b}){
	  &buildConnection();
	  $idx_buff=0;
	  $pre_reg_idx=$reg_idx;
	}
      }
    }

    #recycle the start of the windows
    my $bs1=$bufferend+$opts{s};
    my $bs2=$t->{pos};
    $bufferend=($bs1>$bs2)?$bs1:$bs2;

    my @u;
    foreach my $blib(keys %reg_seq){
      my @tmp_reg_seq=@{$reg_seq{$blib}};
      foreach my $read(@{$reg_seq{$blib}}){
	@u=split " ",$read;
	#last if($u[2]>$bufferend-$winsize{$blib});
	last if($u[2]>$bufferend-$maxwinsize);

	shift @tmp_reg_seq;
	my $bdist=abs($u[4]);
	$bdist=($bdist<=$#{$CMF{$blib}})?$bdist:$#{$CMF{$blib}};

	${$libpmf{$u[3]}{$blib}}{$bdist}--;
	delete ${$libpmf{$u[3]}{$blib}}{$bdist} if(${$libpmf{$u[3]}{$blib}}{$bdist}<=0);
	$nreads_lib{$u[3]}{$blib}--;
      }
      $reg_seq{$blib}=\@tmp_reg_seq;
    }
  }

  #Add reads
  $nreads_lib{$t->{ori}}{$lib}++;
  $t->{readname} =~ s/\/[12]$//;
  push @{$reg_seq{$lib}}, join(" ", $t->{readname},$t->{chr},$t->{pos},$t->{ori},$t->{dist},$t->{flag},$t->{qual},$t->{readlen},$lib);

  $sepdist=($sepdist<=$#{$CMF{$lib}})?$sepdist:$#{$CMF{$lib}};

  ${$libpmf{$t->{ori}}{$lib}}{$sepdist}++;
}

sub buildConnection{

  #assign reads to regions
  for(my $j=$pre_reg_idx;$j<$reg_idx;$j++){
    foreach (@{$regs{$j}}) {
      my @s = split;
      $read{$s[0]}{$j}=1;
    }
  }
  # build connections
  # find paired regions that are supported by paired reads
  #warn("-- link regions\n");
  my %link;
  my %free_nodes;
  foreach my $x (keys %read) {
    my @p = sort {$a<=>$b} keys %{$read{$x}};

    if($#p==0){  #singleton connection
      delete $read{$x}{$p[0]} if(!defined $reg_name{$p[0]});
    }
    elsif($#p>0){   #non-singleton connection
      for(my $i=0;$i<$#p;$i++){
	if(!defined $reg_name{$p[$i]}){   #clean up abolished links
	  delete $read{$x}{$p[$i]};
	  next;
	}
	for(my $j=$i+1;$j<=$#p;$j++){
	  if(!defined $reg_name{$p[$j]}){   #clean up abolished links
	    delete $read{$x}{$p[$j]};
	    next;
	  }
	  my @reg0=split /\t+/,$reg_name{$p[$i]};
	  my @reg1=split /\t+/,$reg_name{$p[$j]};

	  next if($p[$i] ne $p[$j] && $reg0[4] eq $reg1[4]);  #different node can only be paired if in different orientation
	  if (defined $link{$p[$i]}{$p[$j]}) {
	    ++$link{$p[$i]}{$p[$j]};
	    ++$link{$p[$j]}{$p[$i]};
	  } else {
	    $link{$p[$i]}{$p[$j]} = 1;
	    $link{$p[$j]}{$p[$i]} = 1;
	  }
	}
      }
    }
    @p = sort {$a<=>$b} keys %{$read{$x}};
    delete $read{$x} if($#p<0);  #remove died reads
  }
  my %clink=%link;
  # segregate graph, find nodes that have connections
  foreach my $s0(sort {$a<=>$b} keys %clink){
    next unless(defined $clink{$s0});
    #construct a subgraph
    my @tails=($s0);
    while(@tails){
      my @newtails;
      foreach my $tail(@tails){
	next unless(defined $clink{$tail});
	next unless (defined $reg_name{$tail});  # a node must be defined
	my @s1s=keys %{$clink{$tail}};
	foreach my $s1(@s1s){
	  my %nodepair;
	  my $nlinks=$clink{$tail}{$s1};
	  next if($nlinks<$opts{r});  # require sufficient number of pairs
	  next if(defined $nodepair{$s1});  # a node only appear once in a nodepair
	  next unless (defined $reg_name{$s1});  # a node must be defined
	  $nodepair{$tail}{$s1}=$clink{$tail}{$s1};
	  $nodepair{$s1}{$tail}=$clink{$s1}{$tail};
	  delete $clink{$tail}{$s1};  # use a link only once
	  delete $clink{$s1}{$tail};
	  push @newtails,$s1;

	  #Analysis a nodepair
	  my @snodes=sort {$a<=>$b} keys %nodepair;
	  my ($node1,$node2)=@snodes;
	  die "Overlapping node number\n" if(!defined $node2);
	  next unless($nodepair{$node1}{$node2});

	  my $nread_pairs=0;

	  my %type;
	  my %library_readcount;
	  my %library_meanspan;   # diff span distance;
	  my @paired_reads;

	  my @reg0=split /\t+/,$reg_name{$node1};
	  my @reg1=split /\t+/,$reg_name{$node2};
	  printf "try %s\n",join("=>","$node1:".$reg_name{$node1},"$node2:".$reg_name{$node2}) if($opts{v}>0);

	  my @nonpaired_reads1;
	  my @nonpaired_reads2;

	  my %orient_count1; my %orient_count2;
	  my %read_pair;
	  foreach my $y (@{$regs{$node1}}) {
	    my @t = split(" ", $y);
	    $read_pair{$t[0]}=\@t;
	  }

	  my ($start1,$end1,$start2,$end2)=(1e10,0,1e10,0);
	  foreach my $y (@{$regs{$node2}}) {
	    my @t2 = split(" ", $y);

	    if(defined $read_pair{$t2[0]}){  #paired
	      my @t1=@{$read_pair{$t2[0]}};

	      #tighten up the boundary based on paired reads
	      $start1=($start1>$t1[2])?$t1[2]:$start1;
	      $end1=($end1<$t1[2])?$t1[2]:$end1;
	      $start2=($start2>$t2[2])?$t2[2]:$start2;
	      $end2=($end2<$t2[2])?$t2[2]:$end2;

	      my $lib=$t2[8];
	      $orient_count1{pair}++;
	      $orient_count2{pair}++;

	      $library_readcount{$lib}++;
	      $library_meanspan{$lib}+=abs($t2[4]);
	      $nread_pairs++;

	      #free connections
	      delete $read{$t1[0]}{$node1};
	      delete $read{$t2[0]}{$node2};

	      push @paired_reads,$y;
	      delete $read_pair{$t2[0]};
	    }
	    else{
	      push @nonpaired_reads2,$y;
	    }
	  }

	  $orient_count1{all}=$orient_count1{pair};
	  $orient_count2{all}=$orient_count2{pair};
	  foreach my $rdname (keys %read_pair) {
	    my @t=@{$read_pair{$rdname}};
	    next if($t[2]<$start1 || $t[2]>$end1);
	    my $y=join(" ",@t);
	    push @nonpaired_reads1,$y;
	    $orient_count1{all}++;
	  }

	  my @tmparray;
	  foreach my $y(@nonpaired_reads2){
	    my @t2 = split(" ", $y);
	    next if($t2[2]<$start2 || $t2[2]>$end2);
	    push @tmparray,$y;
	    $orient_count2{all}++;
	  }
	  @nonpaired_reads2=@tmparray;

	  my $min_nonpaired=($#nonpaired_reads1>$#nonpaired_reads2)?$#nonpaired_reads2+1:$#nonpaired_reads1+1;
	  $regs{$node1}=\@nonpaired_reads1;
	  $regs{$node2}=\@nonpaired_reads2;

	  my @type_orient_counts=(\%orient_count1,\%orient_count2);

	  if( $nread_pairs>=$opts{r}){
	    my $sptype;
	    my $diffspan=0;
	    foreach my $sp(keys %library_readcount){
	      if($sptype){
		$sptype=join(':',$sptype,$sp . '|'. $library_readcount{$sp});
	      }
	      else{
		$sptype=$sp . '|'. $library_readcount{$sp};
	      }
	      $diffspan+=$library_meanspan{$sp}-$library_readcount{$sp}*$mean_insertsize{$sp};
	    }

	    my $indelsize=int($diffspan/$nread_pairs+0.5);
	    #print out result
	    my ($sv_chr1, $sv_pos1, $sv_ori1, $sv_chr2,$sv_pos2, $sv_ori2);

	    #find outer most coordinates
	    foreach my $node(@snodes){
	      my ($chr,$start,$end,$score,$ori)=split /\t+/, $reg_name{$node};

	      ($start,$end)=($start1,$end1) if($node==$node1);
	      ($start,$end)=($start2,$end2) if($node==$node2);

	      my $ori_readcount=shift @type_orient_counts;
	      if($sv_chr1 && $sv_chr2){
		if($ori eq '-'){
		  $sv_chr1=$sv_chr2;
		  $sv_pos1=$sv_pos1;
		  ($sv_chr2,$sv_pos2)=($chr,$end);
		  $sv_ori2=join('', $$ori_readcount{pair}||'0','/',$$ori_readcount{all}||'0');
		}
		else{  #swap
		  $sv_chr1=$sv_chr2;
		  ($sv_chr2,$sv_pos2)=($chr,$sv_pos2);
		  $sv_pos1=$start;
		  $sv_ori2=$sv_ori1;
		  $sv_ori1=join('', $$ori_readcount{pair}||'0','/',$$ori_readcount{all}||'0');
		}
	      }
	      else{  #singleton region
		$sv_chr1=$chr;
		$sv_chr2=$chr;
		($sv_pos1,$sv_pos2)=($start,$end);
		$sv_ori1=join('', $$ori_readcount{pair}||'0','/',$$ori_readcount{all}||'0');
		$sv_ori2=join('', $$ori_readcount{pair}||'0','/',$$ori_readcount{all}||'0');
	      }
	    }

	    #ROC comparison on venter simulation data suggests that the average score has more robust performance than the max score,
	    #while the max score may boost sensitivity in the case of mixed insert size
	    my $SVscore=0;
	    if($opts{e}==1){
	      $SVscore=($reg0[3]<$reg1[3])?$reg1[3]:$reg0[3];  #take the max
	    }
	    else{
	      $SVscore=($reg0[3]+$reg1[3])*99/200; #take the average
	    }

	    my $flag=($diffspan>0)?2:3;
	    my $AF=($#paired_reads+1)/($#paired_reads+1+$min_nonpaired);  #estimate allele frequency
	    printf "%s\t%d\t%s+\t%s\t%d\t%s-\t%s\t%d\t%d\t%d\t%s\t%.2f\t%s\t%s\n",$sv_chr1,$sv_pos1,$sv_ori1,$sv_chr2,$sv_pos2,$sv_ori2,$SVtype{$flag},$indelsize,$SVscore,$nread_pairs,$sptype,$AF,$version,$options;

	    #record list of nodes that can be potentially freed
	    $free_nodes{$node1}=1;
	    $free_nodes{$node2}=1;
	  }
	}
	delete $clink{$tail};
      }
      @tails=@newtails;
    }
  }

  #free nodes
  foreach my $node(keys %reg_name){
    #remove reads in the regions
    my @reads=@{$regs{$node}};
    my ($chr,$start,$end,$score,$ori)=split /\t+/, $reg_name{$node};
    if($#reads+1<$opts{r} && $free_nodes{$node} ||
       $end<$bufferend-$opts{m} ||  #out of range
       $chr ne $current_chr   #on different chromosome
      ){

      print "remove $node\t$reg_name{$node}\n" if($opts{v}>0);
      #remove regions
      delete $regs{$node};
      delete $reg_name{$node};

    }
  }
}

sub KS_TestB{
  #compute P values of one side KS test based Beta approximation

  my ($D,$n)=@_;
  my $ah=0.003326-6.012/$n+5.52/($n**0.53);
  my $bh=-0.0004245-0.003397/$n+0.3204/($n**0.48);
  my $ph=3.258-3.727/$n+4.607/($n)**1.6;
  my $qh=25-161.2/$n+162.2/($n**1.3);
	
  my $x=($D-$bh)/$ah;
  $x=($x<0)?0:$x;  $x=($x>1)?1:$x;
  my $P=Math::CDF::pbeta($x, $ph, $qh);  #prob of not being an SV
  return $P;
}

sub FishersMethod{
  my @ps=@_;
  my $pvalue;
  if($#ps==0){
    $pvalue=$ps[0];
  }
  else{
    my $s=0;
    my $k=$#ps+1;
    foreach my $p(@ps){
      $s+=($p>$ZERO)?log($p):$LZERO;
    }
    $s*=-2;
    $pvalue=1-Math::CDF::pchisq($s,2*$k);
  }
  return $pvalue;
}

sub ReadRegions{
  my ($f)=@_;
  my %Rgs;
  open(FIN,"<$f") || die "unable to open $f\n";
  while(<FIN>){
    chomp;
    my ($chr,@c)=split;
    push @{$Rgs{$chr}},join(" ",@c);
  }
  close(FIN);
  return \%Rgs;
}
