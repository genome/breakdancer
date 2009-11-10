#!/gsc/bin/perl
#Parse different alignment format and convert it to MAQ format

use strict;
use warnings;

package AlnParser;

sub new{
  my ($class, %arg) = @_;
  my $self={
	   };
  bless($self, $class || ref($class));
  return $self;
}

sub in{
  my ($self,$line,$format,$alt)=@_;
  my $t;
  if($format eq 'maq'){
    my @s=split /\t+/,$line;
    ($t->{readname},$t->{chr},$t->{pos},$t->{ori},$t->{dist},$t->{flag},$t->{qual},$t->{readlen},$t->{seq},$t->{basequal})=@s[0..5,8,13,14,15];
  }
  elsif($format eq 'sam'){
    my @s=split /\t+/,$line;
    my ($flag,$mchr,$mpos);
    ($t->{readname},$flag,$t->{chr},$t->{pos},$t->{qual},$mchr,$mpos,$t->{dist},$t->{seq},$t->{basequal})=@s[0..4,6..10];
    $t->{readlen}=length($t->{seq});
    $t->{ori}=($flag & 0x0010)?'-':'+';
    $t->{flag}=0;
    ($t->{readgroup})=($line=~/(RG\S+)/);
    $t->{readgroup}=~s/.*\://g if(defined $t->{readgroup});
    if(! defined $alt){  #default to alternative mapping quality
      if($line=~/Aq\:i\:(\d+)/){  #if there is an alternative mapping quality flag
	$t->{qual}=$1;
      }
      elsif($line=~/AM\:i\:(\d+)/){
	$t->{qual}=$1;
      }
    }

    #convert to Maq flag
    if($line=~/MF\:i\:(\d+)/){
      $t->{flag}=$1;
    }
    elsif($flag & 0x0001){
      my $ori2=($flag & 0x0020)?'-':'+';
      if($flag & 0x0004){  #read itself unmapped
	$t->{flag}=192;
      }
      elsif($flag & 0x0008){  #mate is unmapped
	$t->{flag}=64;
      }
      elsif($mchr ne '='){  #RP mapped to different chromosome
	$t->{flag}=32;
      }
      elsif($flag & 0x0002){  #read properly mapped
	$t->{flag}=18;
      }
      else{
	if($t->{ori} eq $ori2){   #RP mapped to the same unexpected orientation
	  $t->{flag}=($ori2 eq '+')?1:8;
	}
	elsif($mpos>$t->{pos} && $t->{ori} eq '-' || $t->{pos}>$mpos && $t->{ori} eq '+'){  #larger coordinate read is not on the negative strand
	  $t->{flag}=4;
      	}
	else{
	  $t->{flag}=2;
	}
      }
    }
  }
  else{}

  return $t;
}

1;
