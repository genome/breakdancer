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
    if(defined $alt){
      $t->{qual}=$s[6];
    }
  }
  elsif($format eq 'sam'){
    my @s=split /\t+/,$line;
    my ($flag,$mchr,$mpos);
    ($t->{readname},$flag,$t->{chr},$t->{pos},$t->{qual},$mchr,$mpos,$t->{dist},$t->{seq},$t->{basequal})=@s[0..4,6..10];
    $t->{readlen}=length($t->{seq});
    $t->{ori}=($flag&0x0010 || $flag=~/r/)?'-':'+';  #reverse orientation?
    $t->{flag}=0;  #fragment reads
    ($t->{readgroup})=($line=~/(RG\S+)/);
    $t->{readgroup}=~s/.*\://g if(defined $t->{readgroup});
    #$t->{readgroup}='NA' if(!defined $t->{readgroup});
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
    else{
      if($flag & 0x0001 || $flag=~/p/){  #paired reads
	my $ori2=($flag & 0x0020 || $flag=~/R/)?'-':'+';
	if($flag & 0x0004 || $flag=~/u/){  #read itself unmapped
	  $t->{flag}=192;
	}
	elsif($flag & 0x0008 || $flag=~/U/){  #mate is unmapped
	  $t->{flag}=64;
	}
	elsif($mchr ne '='){  #RP mapped to different chromosome
	  $t->{flag}=32;
	}
	elsif($flag & 0x0002 || $flag=~/P/){  #read properly mapped
	  if($flag & 0x0040){  #first read
	    if($flag & 0x0020){  #mate reversed
	      $t->{flag}=18;
	    }
	    else{  #mate forward
	      $t->{flag}=20;
	    }
	  }
	  else{  #second read
	    if($flag & 0x0010){  #itself reversed
	      $t->{flag}=18;
	    }
	    else{
	      $t->{flag}=20;
	    }
	  }

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
  }
  else{}

  return $t;
}

1;
