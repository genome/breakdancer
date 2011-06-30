# EXAMPLE PIPELINE

Please review the [full usage instructions](breakdancer.html) for detailed information about the *bam2cfg.pl* and *breakdancer-max* commands.

## STEP 1
Create a configuration file using bam2cfg.pl:

The precompiled Debian package will install this in /usr/lib/breakdancer-max{{page.version_suffix}} along with a few required Perl modules.

	/usr/lib/breakdancer-max{{page.version_suffix}}/bam2cfg.pl -g -h tumor.bam normal.bam > BRC6.cfg

## STEP 2
Detect inter-chromosomal translocations:

	breakdancer_max -t -q 10 -d BRC6.ctx BRC6.cfg > BRC6.ctx

The -d option dumps CTX supporting read pairs into fastq files (in this case BRC6.ctx) by library.

This step normally takes 12 hours or so for three bam files, 8 hours or so for two bam files for cpp version, around three days for perl version.
