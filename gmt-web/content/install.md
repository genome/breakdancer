Download (Ubuntu 8.04+)
=======================

* * *

To install the base system if you haven't already, run this in a terminal:

<p class='terminal' markdown='1'>
wget -q -O - http://apt.genome.wustl.edu/ubuntu/files/install.sh | sudo bash
</p>

To install BreakDancer, run:

<p class='terminal' markdown='1'>
sudo gmt install breakdancer{{ page.version_suffix }}
</p>

Behind the Scenes:
------------------

Genome Modeling System uses the native apt package manager, and links to an apt repository served from apt.genome.wustl.edu for systematic updates.

repository linkage:

<p class='terminal' markdown='1'>
sudo apt-add-respository "deb http://apt.genome.wustl.edu lucid-genome main"<br/>
wget http://apt.genome.wustl.edu/ubuntu/files/genome-center.asc | sudo apt-key add<br/>
sudo apt-get update<br/>
</p>

install BreakDancer:

<p class='terminal' markdown='1'>
sudo apt-get install breakdancer{{ page.version_suffix }}
</p>

Running BreakDancer
===================

To run BreakDancer, first use *bam2cfg.pl* to prepare the required per-invocation config file.  See the usage docs for params.
<p class='terminal' markdown='1'>
/usr/lib/breakdancer-max{{ page.version_suffix }}/bam2cfg.pl bam_files breakdancer_options
</p>

Then run BreakDancer on the config
<p class='terminal' markdown='1'>
breakdancer-max config_file.cfg
</p>



