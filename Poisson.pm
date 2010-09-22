#!/gsc/bin/perl

use strict;
use warnings;

package Poisson;

my $LZERO=-10000000000;
my $LSMALL=$LZERO/2;
my $SMALL=exp($LSMALL);
my $minLogExp = -log(-$LZERO);
my $PI=3.1415927;

sub new{
  my ($class, %arg) = @_;
  my $self={
	   };
  bless($self, $class || ref($class));
  return $self;
}

sub LogPoissonTailProb{
  my ($self,$n,$lambda)=@_;
  my $logprob=$LZERO;
  my $plogprob;
  do{
    $plogprob=$logprob;
    $logprob=&LAdd($plogprob,$self->LogPoissonPDF($n++,$lambda));
  } until ($logprob-$plogprob<0.01);
  return $logprob;
}

sub LogPoissonPDF{
  my ($self,$k,$lambda)=@_;
  my $logk_factorial=($k==0)?0:$k*log($k)-$k+0.5*log(2*$PI*$k);
  #my $logp=$k*log($lambda)-$lambda-$logk_factorial;
  my $log_lambda=($lambda<=0)?$LSMALL:log($lambda);
  my $logp=$k*$log_lambda-$lambda-$logk_factorial;

  return $logp;
}

sub LAdd{
  my ($x, $y)=@_;
  my ($temp,$diff,$z);
  if ($x<$y) {
    $temp = $x; $x = $y; $y = $temp;
  }
  $diff = $y-$x;
  if ($diff<$minLogExp){
    return  ($x<$LSMALL)?$LZERO:$x;
  }
  else {
    $z = exp($diff);
    return $x+log(1.0+$z);
  }
  return $z;
}

1;
