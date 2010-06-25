#!/gsc/bin/perl
#joint analysis of multiple maps (lib specific)
#  Copyright (C) 2009 Washington University in St. Louis

use strict;
use warnings;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$FindBin::Bin";
#use lib '/gscuser/kchen/1000genomes/analysis/scripts/';
use Poisson;
use AlnParser;
use Statistics::Descriptive;
use Math::CDF;
use IO::File;

my $version="BreakDancerMax-1.0r112";
my %opts = (i=>200, c=>3, m=>1000000, q=>35, s=>7, r=>2, b=>100, p=>0.001, x=>1000);
my %opts1;
getopts('o:s:c:m:q:r:b:ep:tfd:g:lCx:', \%opts1);
#         -C        change system default from Illumina to SOLiD
die("
Usage:   BreakDancerMax.pl <analysis_config.lst>
Options:
         -o STRING operate on a single chromosome [all chromosome]
         -c INT    cutoff in unit of standard deviation [$opts{c}]
         -t 	   only detect transchromosomal rearrangement
         -p FLOAT  prior probability of SV [$opts{p}]
         -e        learn parameters from data before applying to SV detection
         -l        analyze Illumina long insert (mate-pair) library
         -f        use Fisher's method to combine P values from multiple library
         -d STRING prefix of fastq files that SV supporting reads will be saved by library
         -g FILE   dump SVs and supporting reads in BED format for GBrowse
         -m INT    maximum SV size [$opts{m}]
         -q INT    minimum alternative mapping quality [$opts{q}]
         -s INT    minimum length of a region [$opts{s}]
         -b INT    buffer size for building connection [$opts{b}]
         -r INT    minimum number of read pairs required to establish a connection
         -x INT    ignore regions that have greater than [$opts{x}] fold haploid sequence coverage
Version: $version\n
") unless (@ARGV);

my $options='';
foreach my $opt(keys %opts1){
  $options.='|'.$opt.$opts1{$opt};
  $opts{$opt}=$opts1{$opt};
}

my %SVtype;
if($opts{l}){  #Illumina long insert
  %SVtype=(
	   '1'=>'INV',  #inversion
	   '3'=>'INS',  #insertion
	   '4'=>'DEL',  #intra-chromosomal translocation
	   '8'=>'INV',  #inversion
	   '32'=>'CTX'  #inter-chromosomal translocation
	  );
}
else{
  %SVtype=(
	   '1'=>'INV',  #inversion
	   '2'=>'DEL',  #deletion
	   '3'=>'INS',  #insertion
	   '4'=>'ITX',  #intra-chromosomal translocation
	   '8'=>'INV',  #inversion
	   '32'=>'CTX'  #inter-chromosomal translocation
	  );
}

my $AP=new AlnParser(platform=>$opts{C});
my $LZERO=-99;
my $ZERO=exp($LZERO);
my %exes;
my %fmaps;
my %libmaps;
my %mean_insertsize;
my %std_insertsize;
my %uppercutoff;
my %lowercutoff;
my %readlens;
my %mapQual;
my $max_readlen=0;
my %x_readcounts;   #reads spanning longer than expect distance
my %readgroup_library;
my %readgroup_platform;
my %ReadsOut;
my $d=1e10;
open(CONFIG,"<$ARGV[0]") || die "unable to open $ARGV[0]\n";
while(<CONFIG>){
  next unless (/\S+/);
  next if(/^\#/);
  chomp;
  my $fh;
  my ($fmap)=($_=~/map\:(\S+)\b/i);

  my ($mean)=($_=~/mean\w*\:(\S+)\b/i);
  my ($std)=($_=~/std\w*\:(\S+)\b/i);
  my ($readlen)=($_=~/readlen\w*\:(\S+)\b/i);
  my ($upper,$lower);
  ($upper)=($_=~/upp\w*\:(\S+)\b/i);
  ($lower)=($_=~/low\w*\:(\S+)\b/i);
  my ($mqual)=($_=~/map\w*qual\w*\:(\d+)\b/i);
  my ($lib)=($_=~/lib\w*\:(\S+)\b/i);
  ($lib)=($_=~/samp\w*\:(\S+)\b/i) if(!defined $lib);
  $lib='NA' if(!defined $lib);

  my ($readgroup)=($_=~/group\:(\S+)\b/i);
  $readgroup=$lib if(!defined $readgroup);
  $readgroup_library{$readgroup}=$lib;

  my ($platform)=($_=~/platform\:(\S+)\b/i);
  $readgroup_platform{$readgroup}=$platform || (($opts{C})?'solid':'illumina');  #default to illumina

  my ($exe)=($_=~/exe\w*\:(.+)\b/i);
  if(defined $opts{d}){
    open($ReadsOut{$lib.'1'},">$opts{d}.$lib.1.fastq") || die "unable to open $lib.1.fastq, check write permission\n";
    open($ReadsOut{$lib.'2'},">$opts{d}.$lib.2.fastq") || die "unable to open $lib.2.fastq, check write permission\n";
  }
  $libmaps{$lib}=$fmap;
  $mapQual{$lib}=$mqual if(defined $mqual && ($mqual=~/^\d+$/));
  $fmaps{$fmap}=$lib;

  if(defined $mean && defined $std && (!defined $upper || !defined $lower)){
    $upper=$mean+$std*$opts{c};
    $lower=$mean-$std*$opts{c};
    $lower=($lower>0)?$lower:0;
  }
  $max_readlen=($max_readlen<$readlen)?$readlen:$max_readlen;

  $mean_insertsize{$lib}=$mean;
  $std_insertsize{$lib}=$std;
  $uppercutoff{$lib}=$upper;
  $lowercutoff{$lib}=$lower;
  $readlens{$lib}=$readlen;

  if(!defined $exes{$fmap}){
    $exes{$fmap}=$exe || 'cat';
  }
  else{
    die "Please use identical exe commands to open the same map file.\n" if($exes{$fmap} ne $exe);
  }
  my $tmp=$mean - $readlen*2;  #this determines the mean of the max of the SV flanking region
  $d=($d<$tmp)?$d:$tmp;
}
close(CONFIG);
$d=50 if($d<50);

open(BED,">$opts{g}") if (defined $opts{g});
my @maps=keys %fmaps;
my @format;
for(my $i=0;$i<=$#maps;$i++){
  my $exe=$exes{$maps[$i]};
  if($exe=~/maq/){
    push @format,'maq';
  }
  else{
    push @format,'sam';
  }
}
&EstimatePriorParameters() if($opts{e});

my $reference_len=1;
my %nreads;
my %cmds;
my $defined_all_readgroups=1;
for(my $i=0;$i<=$#maps;$i++){
  my $fh;
  my $ref_len;
  my $exe=$exes{$maps[$i]};
  $cmds{$exe}++;

  if($exe=~/samtools/){
    $exe=join(" ", $exe, $maps[$i], $opts{o} || '');
    open($fh,"$exe |") || die "unable to open $maps[$i]\n";
  }
  elsif($exe=~/maq/){
    open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
  }
  else{
    open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
  }

  my $p_pos=0;
  my $p_chr='';
  while(<$fh>){
    chomp;
    my $t=$AP->in($_,$format[$i],\%readgroup_platform);
    next if(!defined $t || defined $opts{o} && $t->{chr} ne $opts{o});  #analyze only one chromosome
    $ref_len+=$t->{pos}-$p_pos if($t->{chr} eq $p_chr);
    $p_pos=$t->{pos};
    $p_chr=$t->{chr};
    my $lib;
    if(defined $t->{readgroup}){
      $lib=$readgroup_library{$t->{readgroup}};
    }
    else{
      $defined_all_readgroups=0;
      $lib=$fmaps{$maps[$i]};  #when multiple libraries are in a BAM file
    }
    next unless(defined $lib);

    $nreads{$lib}++;
    if(defined $mapQual{$lib}){
      next if ($t->{qual}<=$mapQual{$lib});
    }
    else{
      next if ($t->{qual}<=$opts{q});
    }
    next if($t->{flag}==0); #return fragment reads
    next if($opts{t} && ($t->{flag}<32 || $t->{flag}>=64));  #only care flag >32 for CTX

    if($opts{l}){
      $t->{flag}=4 if(abs($t->{dist})>$uppercutoff{$lib} && $t->{flag}==20);
      $t->{flag}=20 if(abs($t->{dist})<$uppercutoff{$lib} && $t->{flag}==4);
      $t->{flag}=3 if(abs($t->{dist})<$lowercutoff{$lib} && $t->{flag}==20);
    }
    else{
      $t->{flag}=2 if(abs($t->{dist})>$uppercutoff{$lib} && $t->{flag}==18);
      $t->{flag}=18 if(abs($t->{dist})<$uppercutoff{$lib} && $t->{flag}==2);
      $t->{flag}=3 if(abs($t->{dist})<$lowercutoff{$lib} && $t->{flag}==18);
    }

    next if ($t->{flag} == 18 || $t->{flag}==20 || $t->{flag}==130);
    $x_readcounts{$t->{flag}}{$lib}++;
  }
  close($fh);
  die "$maps[$i] does not contain legitimate paired end alignment. Please check that you have the correct paths and the map/bam files are properly formated and indexed." if(!defined $ref_len);
  $reference_len=$ref_len if($reference_len<$ref_len);
}
my $merge=((keys %cmds)==1 && $defined_all_readgroups)?1:0;

my $total_phy_cov=0;
my $total_seq_cov=0;
my @recflags=sort{$a<=>$b} keys %x_readcounts;
foreach my $lib(keys %nreads){
  my $sequence_coverage=$nreads{$lib}*$readlens{$lib}/$reference_len;
  $total_seq_cov+=$sequence_coverage;
  my $physical_coverage=$nreads{$lib}*$mean_insertsize{$lib}/2/$reference_len;
  $total_phy_cov+=$physical_coverage;

  my $nread_lengthDiscrepant=$x_readcounts{'2'}{$lib} if($x_readcounts{'2'}{$lib});
  $nread_lengthDiscrepant+=$x_readcounts{'3'}{$lib} if($x_readcounts{'3'}{$lib});
  my $tmp=(defined $nread_lengthDiscrepant && $nread_lengthDiscrepant>0)?$reference_len/$nread_lengthDiscrepant:50;
  $d=($d<$tmp)?$d:$tmp;

  printf STDOUT "#%s\tmean:%.3f\tstd:%.3f\tuppercutoff:%.3f\tlowercutoff:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx", $libmaps{$lib},$mean_insertsize{$lib},$std_insertsize{$lib},$uppercutoff{$lib},$lowercutoff{$lib},$readlens{$lib},$lib,$reference_len, $sequence_coverage,$physical_coverage;
  #foreach my $t(sort keys %SVtype){
  foreach my $t(@recflags){
    printf STDOUT "\t%s\:%d",$t,$x_readcounts{$t}{$lib} || 0;
  }
  print STDOUT "\n";
}
print "#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tAllele_frequency\tVersion\tRun_Param\n";

my ($begins, $beginc, $lasts, $lastc) = ('', -1, '', -1);
my (%regs,   # tabulate regions flanking alignment gaps
    %read,   # assign read pairs to one or more region
    %reg_name    # chromosome, start, end, and length classes of the regions
   );
my @reg_seq;

#warn("-- read the mapview output\n");
my $idx_buff=0;
my $final_buff=0;
my $reg_idx=0;
my $normal_switch=0;
my $nnormal_reads=0;
my $ntotal_nucleotides=0;
$max_readlen=0;

if($merge && (@maps>1) && $format[0] eq 'sam' && !defined $opts{o}){
  # open pipe, improvement made by Ben Oberkfell (boberkfe@genome.wustl.edu)
  # samtools merge - in1.bam in2.bam in3.bam in_N.bam | samtools view - 
  # maq mapmerge

  my $merge_command_line=sprintf("samtools merge - %s | samtools view - ", join(" ", @maps));
  my $merge_fh = IO::File->new($merge_command_line . "|");

  while(<$merge_fh>){
    my $t = $AP->in($_, $format[0],\%readgroup_platform);
    next unless(defined $t);
    my $library=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[0]};  #use statistics from the first library
    next unless(!defined $opts{o} || $t->{chr} eq $opts{o});  #analyze only 1 chromosome
    &Analysis($library, $t) if(defined $library);
  }
  $merge_fh->close();
}
else{  #customized merge & sort
  my @FHs;
  my %Idxs;
  for(my $i=0;$i<=$#maps;$i++){
    my $fh;
    my $exe=$exes{$maps[$i]};

    if($exe=~/samtools/){
      $exe=join(" ", $exe, $maps[$i], $opts{o} || '');
      open($fh,"$exe |") || die "unable to open $maps[$i]\n";
    }
    else{
      open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
    }
    push @FHs,$fh;
    $Idxs{$i}=1;
  }

  my @buffer;
  my @cIdxs=keys %Idxs;
  while(keys %Idxs){
    for(my $i=0;$i<=$#cIdxs;$i++){
      my $idx=$cIdxs[$i];
      my $fh=$FHs[$idx];
      if(eof($fh)){
	delete($Idxs{$idx});
	next;
      }
      my $t;
      do{
	$_=<$fh>;
	chomp;
	$buffer[$idx]=$_;
	$t=$AP->in($_,$format[$i],\%readgroup_platform);
      } until(defined $t && (defined $opts{o} && $t->{chr} eq $opts{o} || !defined $opts{o} || eof($fh)));  #analyze only 1 chromosome
    }

    my ($minchr,$minpos)=(chr(255),1e10);
    my $minidx;
    my $min_t;
    my $library;
    foreach my $i(keys %Idxs){
      my $t=$AP->in($buffer[$i],$format[$i],\%readgroup_platform);
      next if(!defined $t || $t->{chr} gt $minchr || $t->{chr} eq $minchr && $t->{pos} > $minpos);
      $minchr=$t->{chr};
      $minpos=$t->{pos};
      $minidx=$i;
      $min_t=$t;
      $library=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[$i]};
    }
    if(defined $minidx){
      &Analysis($library, $min_t) if(defined $library);
      @cIdxs=($minidx);
    }
  }
}

$final_buff=1;
&buildConnection();

if(defined $opts{d}){
  foreach my $fh(keys %ReadsOut){
    close($fh);
  }
}
close(BED) if (defined $opts{g});

sub Analysis{
  my ($lib,$t)=@_;

  #main analysis code
  #return if($t->{qual}<$opts{q} && $t->{flag}!=64 && $t->{flag}!=192);   #include unmapped reads, high false positive rate
  if(defined $mapQual{$lib}){
    return if ($t->{qual}<=$mapQual{$lib});
  }
  else{
    return if ($t->{qual}<=$opts{q});
  }
  return if($t->{chr} eq '*');  #ignore reads that failed to associate with a reference
  return if($t->{flag}==0); #return fragment reads
  return if($opts{t} && ($t->{flag}<32  || $t->{flag}>=64));  #only care flag 32 for CTX

  if($opts{l}){  #for long insert
    $t->{flag}=4 if(abs($t->{dist})>$uppercutoff{$lib} && $t->{flag}==20);
    $t->{flag}=20 if(abs($t->{dist})<$uppercutoff{$lib} && $t->{flag}==4);
    $t->{flag}=3 if(abs($t->{dist})<$lowercutoff{$lib} && $t->{flag}==20);
  }
  else{
    $t->{flag}=2 if(abs($t->{dist})>$uppercutoff{$lib} && $t->{flag}==18);
    $t->{flag}=18 if(abs($t->{dist})<$uppercutoff{$lib} && $t->{flag}==2);
    $t->{flag}=3 if(abs($t->{dist})<$lowercutoff{$lib} && $t->{flag}==18);
    $t->{flag}=4 if($t->{flag}==20); #if it is RF orientation, then regardless of distance
  }
  $t->{flag}=1 if($t->{flag}==8);  #both flag 8 and 1 indicate inversion

  return if ($t->{flag}<32 && abs($t->{dist})>$opts{m});  #skip read pairs mapped too distantly on the same chromosome

  if ($t->{flag} == 18 || $t->{flag} == 20 || $t->{flag} == 130){
    $nnormal_reads++ if($normal_switch && $t->{dist}>0);
    return;
  }

  if($normal_switch){
    $ntotal_nucleotides+=$t->{readlen};
    $max_readlen=($max_readlen<$t->{readlen})?$t->{readlen}:$max_readlen;
  }
  my $do_break = ($t->{chr} ne $lasts || $t->{pos} - $lastc > $d)? 1 : 0;
  if ($do_break) {  # breakpoint in the assembly
    my $seq_coverage=$ntotal_nucleotides/($lastc-$beginc+1+$max_readlen);
    if (($lastc - $beginc > $opts{s}) && ($seq_coverage<$opts{x})) { # skip short/unreliable flanking supporting regions
      #register reliable region and supporting reads across gaps
      my $k=$reg_idx++;
      $reg_name{$k}="$begins\t$beginc\t$lastc\t$nnormal_reads";

      my @p;
      foreach (@reg_seq) {
	push(@p, $_);
	my @s = split;
	push(@{$read{$s[0]}}, $k);
      }
      $regs{$k}=\@p;
      $idx_buff++;
      if($idx_buff>$opts{b}){
	&buildConnection();
	$idx_buff=0;
      }
    }
    else{
      foreach (@reg_seq) {
	my @s = split;
	delete $read{$s[0]};
      }
    }
    ($begins, $beginc) = ($t->{chr},$t->{pos});
    @reg_seq = ();
    $normal_switch=0;
    $nnormal_reads=0;
    $max_readlen=0;
    $ntotal_nucleotides=0;
  }
  $t->{readname} =~ s/\/[12]$//;
  if(defined $opts{d} && defined $t->{seq} && defined $t->{basequal}){
    push @reg_seq, join(" ", $t->{readname},$t->{chr},$t->{pos},$t->{ori},$t->{dist},$t->{flag},$t->{qual},$t->{readlen},$lib,$t->{seq},$t->{basequal});
  }
  else{
    push @reg_seq, join(" ", $t->{readname},$t->{chr},$t->{pos},$t->{ori},$t->{dist},$t->{flag},$t->{qual},$t->{readlen},$lib);
  }
  if($#reg_seq==0){
    $normal_switch=1;
  }
  ($lasts, $lastc) = ($t->{chr},$t->{pos});
}

sub buildConnection{
  # build connections
  # find paired regions that are supported by paired reads
  #warn("-- link regions\n");
  my %link;
  foreach my $x (sort keys %read) {
    my $p = $read{$x};
    next if (@$p != 2);  #skip singleton read (non read pairs)
    if (defined $link{$p->[0]}{$p->[1]}) {
      ++$link{$p->[0]}{$p->[1]};
      ++$link{$p->[1]}{$p->[0]};
    } else {
      $link{$p->[0]}{$p->[1]} = 1;
      $link{$p->[1]}{$p->[0]} = 1;
    }
  }
  my %clink=%link;
  # segregate graph, find nodes that have connections
  my %free_nodes;
  foreach my $s0(sort {$a<=>$b} keys %clink){
    next unless(defined $clink{$s0});
    #construct a subgraph
    my @tails=($s0);
    while(@tails){
      my @newtails;
      foreach my $tail(@tails){
	next unless(defined $clink{$tail});
	next unless (defined $reg_name{$tail});  # a node must be defined
	my @s1s=sort {$a<=>$b} keys %{$clink{$tail}};
	foreach my $s1(@s1s){
	  my @free_reads;
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
	  $node2=$node1 if(!defined $node2);
	  next unless($nodepair{$node1}{$node2});

	  my $nread_pairs=0;
	  my %read_pair;
	  my %type;
	  my %type_library_readcount;
	  my @type_orient_counts;
	  my %type_library_meanspan;   # diff span distance;
	  my @support_reads;
	  foreach my $node(@snodes){
	    my %orient_count;
	    my @nonsupportives;
	    foreach my $y (@{$regs{$node}}){
	      my @t = split(" ", $y);
	      next unless(defined $read{$t[0]});
	      $orient_count{$t[3]}++;

	      if(!defined $read_pair{$t[0]}){
		$read_pair{$t[0]}=$y;
		push @nonsupportives,$y;
	      }
	      else{
		$type{$t[5]}++;
		$type_library_readcount{$t[5]}{$t[8]}++;
		$type_library_meanspan{$t[5]}{$t[8]}+=abs($t[4]);
		$nread_pairs++;
		push @free_reads,$t[0];
		push @support_reads,$y;
		push @support_reads,$read_pair{$t[0]};
		delete $read_pair{$t[0]};
	      }
	    }
	    $regs{$node}=\@nonsupportives;
	    push @type_orient_counts,\%orient_count;
	  }

	  #clean out supportive reads
	  foreach my $node(@snodes){
	    my @nonsupportives;
	    foreach my $y (@{$regs{$node}}){
	      my @t = split(" ", $y);
	      next unless(defined $read_pair{$t[0]});
	      push @nonsupportives,$y;
	    }
	    $regs{$node}=\@nonsupportives;
	  }

	  my ($score,$bestIndelSize);
	  if( $nread_pairs>=$opts{r}){
	    my $maxscore=0;
	    my $flag=0;
	    my %diffspans;
	    my %sptypes;
	    foreach my $fl(sort {$a<=>$b} keys %type){
	      my $ptype=$type{$fl}/$nread_pairs;
	      if($maxscore<$ptype){
		$maxscore=$ptype;
		$flag=$fl;
	      }
	      my $sptype;
	      my $diffspan=0;
	      foreach my $sp(keys %{$type_library_readcount{$fl}}){
		if($sptype){
		  $sptype=join(':',$sptype,$sp . '|'. $type_library_readcount{$fl}{$sp});
		}
		else{
		  $sptype=$sp . '|'. $type_library_readcount{$fl}{$sp};
		}
		$diffspan+=$type_library_meanspan{$fl}{$sp}-$type_library_readcount{$fl}{$sp}*$mean_insertsize{$sp};
	      }
	      $diffspans{$fl}=int($diffspan/$type{$fl}+0.5);
	      $sptypes{$fl}=$sptype;
	    }

	    if($type{$flag}>=$opts{r}){
	      #print out result
	      my ($sv_chr1, $sv_pos1, $sv_ori1, $sv_chr2,$sv_pos2, $sv_ori2);
	      my $normal_rp;
	      #find inner most positions
	      foreach my $node(@snodes){
		my ($chr,$start,$end,$nrp)=split /\t+/, $reg_name{$node};
		my $ori_readcount=shift @type_orient_counts;
		if($sv_chr1 && $sv_chr2){
		  $sv_chr1=$sv_chr2;
		  $sv_pos1=$sv_pos2;
		  ($sv_chr2,$sv_pos2)=($chr,$start);
		  $sv_ori2=join('', $$ori_readcount{'+'}||'0','+',$$ori_readcount{'-'}||'0','-');
		}
		else{
		  $sv_chr1=$chr;
		  $sv_chr2=$chr;
		  ($sv_pos1,$sv_pos2)=($start,$end);
		  $sv_ori1=join('', $$ori_readcount{'+'}||'0','+',$$ori_readcount{'-'}||'0','-');
		  $sv_ori2=$sv_ori1;
		  $normal_rp=$nrp;
		}
	      }

	      my $LogPvalue=&ComputeProbScore(\@snodes, $type_library_readcount{$flag},$flag);
	      my $PhredQ=-10*($LogPvalue/log(10));
	      $PhredQ=($PhredQ>99)?99:int($PhredQ+0.5);

	      my $AF=$type{$flag}/($type{$flag}+$normal_rp);

	      #$sv_pos2+=$max_readlen;  #apply extra padding to the end coordinates
	      $sv_pos1+=$max_readlen-5;  #apply extra padding to the start coordinates

	      my $SVT=$SVtype{$flag} || 'UN';  # UN stands for unknown
	      printf "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%.2f\t%s\t%s\n",$sv_chr1,$sv_pos1,$sv_ori1,$sv_chr2,$sv_pos2,$sv_ori2,$SVT,$diffspans{$flag},$PhredQ,$type{$flag},$sptypes{$flag},$AF,$version,$options;

	      if($opts{d}){  #print out supporting read pairs
		my %pairing;
		foreach my $y (@support_reads) {
		  my @t = split(" ", $y);
		  next unless($#t==10 && $t[5] eq $flag);
		  my $fh=(defined $pairing{$t[0]})?$ReadsOut{$t[8].'1'}:$ReadsOut{$t[8].'2'};
		  $pairing{$t[0]}=1;
		  printf $fh "\@%s\n",$t[0];
		  #printf $fh "\@%s %s %d %s %d %s %d\n",$t[0],$sv_chr1,$sv_pos1,$SVT,$diffspans{$flag},$t[8],$t[6];
		  printf $fh "%s\n",$t[9];
		  printf $fh "\+\n";
		  printf $fh "%s\n",$t[10];
		}
	      }

	      if($opts{g}){  #print out SV and supporting reads in BED format
		my $trackname=join('_',$sv_chr1,$sv_pos1,$SVT,$diffspans{$flag});
		printf BED "track name=%s  description=\"BreakDancer %s %d %s %d\" useScore=0\n",$trackname,$sv_chr1,$sv_pos1,$SVT,$diffspans{$flag};
		foreach my $y (@support_reads) {
		  my @t = split(" ", $y);
		  next unless($#t>=8 && $t[5] eq $flag);
		  $t[1]=~s/chr//;
		  my $aln_end=$t[2]+$t[7]-1;
		  my $color=($t[3] eq '+')?'255,0,0':'0,0,255';
		  #printf BED "chr%s\t%d\t%d\t%s\t1\t%s\t%d\t%d\t%s\n",$t[1],$t[2],$aln_end,$t[0],$t[3],$t[2],$aln_end,$color;
		  my $aln_score=$t[6]*10;  #mapping quality * 10
		  printf BED "chr%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n",$t[1],$t[2],$aln_end,join('|',$t[0],$t[8]),$aln_score,$t[3],$t[2],$aln_end,$color;
		}
	      }
	    }
	    #free reads
	    foreach my $readname(@free_reads){
	      delete $read{$readname};
	    }
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
  foreach my $node(sort {$a<=>$b} keys %free_nodes){
    #remove reads in the regions
    my @reads=@{$regs{$node}};
    if($#reads+1<$opts{r}){
      foreach my $y (@reads) {
	my @t = split(" ", $y);
	my $readname=$t[0];
	delete $read{$readname};
      }
      #remove regions
      delete $regs{$node};
      delete $reg_name{$node};
    }
  }
}

sub EstimatePriorParameters{
  my %es_means;
  my %es_stds;
  my %es_readlens;
  my %es_uppercutoff;
  my %es_lowercutoff;
  my %insert_stat;
  my %readlen_stat;
  for(my $i=0;$i<=$#maps;$i++){
    my $fh;
    my $exe=$exes{$maps[$i]};
    if($exe){
      open($fh,"$exe $maps[$i] |") || die "unable to open $maps[$i]\n";
    }
    else{
      open($fh,"<$maps[$i]") || die "unable to open $maps[$i]\n";
    }

    while(<$fh>){
      chomp;
      my $t=$AP->in($_,$format[$i],\%readgroup_platform);
      next unless(defined $t);
      my $lib=($t->{readgroup})?$readgroup_library{$t->{readgroup}}:$fmaps{$maps[$i]};  #when multiple libraries are in a BAM file
      next unless(defined $lib);
      next if(defined $opts{o} && $t->{chr} ne $opts{o});  #analysis 1 chromosome
      $readlen_stat{$lib}=Statistics::Descriptive::Sparse->new() if(!defined $readlen_stat{$lib});
      $readlen_stat{$lib}->add_data($t->{readlen});
      next if ($t->{qual}<=$opts{q});  #skip low quality mapped reads
      next if($t->{flag}!=18 && $t->{flag}!=20 || $t->{dist}<=0);
      $insert_stat{$lib}=Statistics::Descriptive::Sparse->new() if(!defined $insert_stat{$lib});
      $insert_stat{$lib}->add_data($t->{dist});
    }
    close($fh);
  }
  foreach my $lib(keys %readlen_stat){
    my $mean=$insert_stat{$lib}->mean();
    my $std=$insert_stat{$lib}->standard_deviation();
    my $uppercutoff=$mean+$std*$opts{c};
    my $lowercutoff=$mean-$std*$opts{c};
    $es_readlens{$lib}=$readlen_stat{$lib}->mean();
    $es_means{$lib}=$mean;
    $es_stds{$lib}=$std;
    $es_uppercutoff{$lib}=$uppercutoff;
    $es_lowercutoff{$lib}=$lowercutoff;

  }
  %mean_insertsize=%es_means;
  %std_insertsize=%es_stds;
  %uppercutoff=%es_uppercutoff;
  %lowercutoff=%es_lowercutoff;
  %readlens=%es_readlens;
}

sub ComputeProbScore{
  my ($rnode,$rlibrary_readcount,$type)=@_;
  my $total_region_size=&PutativeRegion($rnode);

  my $lambda;
  my $logpvalue=0;
  my @libs=keys %{$rlibrary_readcount};
  foreach my $lib(@libs){
    $lambda=$total_region_size*$x_readcounts{$type}{$lib}/$reference_len;
    my $p=new Poisson();
    $logpvalue+=$p->LogPoissonTailProb($$rlibrary_readcount{$lib},$lambda);
  }

  if($opts{f} && $logpvalue<0){
      #Fisher's Method
      my $fisherP=1-Math::CDF::pchisq(-2*$logpvalue,2*($#libs+1));
      $logpvalue=($fisherP>$ZERO)?log($fisherP):$LZERO;
  }
  return $logpvalue;
}

sub PutativeRegion{
  my ($rnode)=@_;
  my $total_region_size=0;
  foreach my $node(@{$rnode}){
    my ($clust_chr,$clust_start,$clust_end,@extra)=split /\t+/,$reg_name{$node};
    $total_region_size += $clust_end - $clust_start + 1;
  }
  return $total_region_size;
}
