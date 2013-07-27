#include "Options.hpp"

#include "version.h"

#include <cstdio>
#include <getopt.h>

using namespace std;

Options::Options()
    : min_len(7)
    , cut_sd(3)
    , max_sd(1000000000)
    , min_map_qual(35)
    , min_read_pair(2)
    , seq_coverage_lim(1000)
    , buffer_size(100)
    , transchr_rearrange(false)
    , fisher(false)
    , Illumina_long_insert(false)
    , CN_lib(false)
    , print_AF(false)
    , score_threshold(30)
{
}

Options::Options(int argc, char** argv)
        : min_len(7)
        , cut_sd(3)
        , max_sd(1000000000)
        , min_map_qual(35)
        , min_read_pair(2)
        , seq_coverage_lim(1000)
        , buffer_size(100)
        , transchr_rearrange(false)
        , fisher(false)
        , Illumina_long_insert(false)
        , CN_lib(false)
        , print_AF(false)
        , score_threshold(30)

{
    int c;
    while((c = getopt(argc, argv, "o:s:c:m:q:r:x:b:tfd:g:lahy:C:R:")) >= 0) {
        switch(c) {
            case 'C': cache_file = optarg; break;
            case 'R': {
                if (argc != 3) {
                    throw runtime_error("When using -R, no other options are allowed");
                }
                restore_file = optarg;
                return;
                break;
            }
            case 'o': chr = optarg; break;
            case 's': min_len = atoi(optarg); break;
            case 'c': cut_sd = atoi(optarg); break;
            case 'm': max_sd = atoi(optarg); break;
            case 'q': min_map_qual = atoi(optarg); break;
            case 'r': min_read_pair = atoi(optarg); break;
            case 'x': seq_coverage_lim = atoi(optarg); break;
            case 'b': buffer_size = atoi(optarg); break;
            case 't': transchr_rearrange = true; break;
            case 'f': fisher = true; break;
            case 'd': prefix_fastq = optarg; break;
            case 'g': dump_BED = optarg; break;
            case 'l': Illumina_long_insert = true; break;
            case 'a': CN_lib = true; break;
            case 'h': print_AF = true; break;
            case 'y': score_threshold = atoi(optarg); break;
            default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
                exit(1);
        }
    }

    // FIXME: instead of printing out defaults, this will print any partial options
    // specified. Let's try to get clearance to use boost::program_options or something
    // more reasonable.
    if(optind == argc) {
        fprintf(stderr, "\nbreakdancer-max version %s (commit %s)\n\n", __g_prog_version, __g_commit_hash);
        fprintf(stderr, "Usage: breakdancer-max <analysis.config>\n\n");
        fprintf(stderr, "Options: \n");
        fprintf(stderr, "       -o STRING       operate on a single chromosome [all chromosome]\n");
        fprintf(stderr, "       -s INT          minimum length of a region [%d]\n", min_len);
        fprintf(stderr, "       -c INT          cutoff in unit of standard deviation [%d]\n", cut_sd);
        fprintf(stderr, "       -m INT          maximum SV size [%d]\n", max_sd);
        fprintf(stderr, "       -q INT          minimum alternative mapping quality [%d]\n", min_map_qual);
        fprintf(stderr, "       -r INT          minimum number of read pairs required to establish a connection [%d]\n", min_read_pair);
        fprintf(stderr, "       -x INT          maximum threshold of haploid sequence coverage for regions to be ignored [%d]\n", seq_coverage_lim);
        fprintf(stderr, "       -b INT          buffer size for building connection [%d]\n", buffer_size);
        fprintf(stderr, "       -t              only detect transchromosomal rearrangement, by default off\n");
        //fprintf(stderr, "    -f INT    use Fisher's method to combine P values from multiple library [%d]\n", fisher);
        fprintf(stderr, "       -d STRING       prefix of fastq files that SV supporting reads will be saved by library\n");
        fprintf(stderr, "       -g STRING       dump SVs and supporting reads in BED format for GBrowse\n");
        fprintf(stderr, "       -l              analyze Illumina long insert (mate-pair) library\n");
        fprintf(stderr, "       -a              print out copy number and support reads per library rather than per bam, by default off\n");
        fprintf(stderr, "       -h              print out Allele Frequency column, by default off\n");
        fprintf(stderr, "       -y INT          output score filter [%d]\n", score_threshold);
        //fprintf(stderr, "Version: %s\n", version);
        fprintf(stderr, "\n");
        exit(1);
    }

    // define the map SVtype
    if (Illumina_long_insert) {
        SVtype[breakdancer::ARP_FF] = "INV";
        SVtype[breakdancer::ARP_FR_small_insert] = "INS";
        SVtype[breakdancer::ARP_RF] = "DEL";
        SVtype[breakdancer::ARP_RR] = "INV";
        SVtype[breakdancer::ARP_CTX] = "CTX";
    }
    else {
        SVtype[breakdancer::ARP_FF] = "INV";
        SVtype[breakdancer::ARP_FR_big_insert] = "DEL";
        SVtype[breakdancer::ARP_FR_small_insert] = "INS";
        SVtype[breakdancer::ARP_RF] = "ITX";
        SVtype[breakdancer::ARP_RR] = "INV";
        SVtype[breakdancer::ARP_CTX] = "CTX";
    }

    bam_config_path = argv[optind];
}

bool Options::operator==(Options const& rhs) const {
    return chr == rhs.chr
        && min_len == rhs.min_len
        && cut_sd == rhs.cut_sd
        && max_sd == rhs.max_sd
        && min_map_qual == rhs.min_map_qual
        && min_read_pair == rhs.min_read_pair
        && seq_coverage_lim == rhs.seq_coverage_lim
        && buffer_size == rhs.buffer_size
        && transchr_rearrange == rhs.transchr_rearrange
        && fisher == rhs.fisher
        && Illumina_long_insert == rhs.Illumina_long_insert
        && CN_lib == rhs.CN_lib
        && print_AF == rhs.print_AF
        && score_threshold == rhs.score_threshold
        && bam_file == rhs.bam_file
        && prefix_fastq == rhs.prefix_fastq
        && dump_BED == rhs.dump_BED
        && SVtype == rhs.SVtype
    ;
}

bool Options::operator!=(Options const& rhs) const {
    return !(*this == rhs);
}
