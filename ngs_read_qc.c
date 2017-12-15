/*
 * ngs_read_qc.c 
 * 
 * Raw read QC using Heng Li's kseq.h library
 * 
 * I wrote this code due to "Too slow" speed of Python language
 * Author : Sunghoon Heo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>  // reading gz compressed file
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

unsigned long n_total_reads = 0;
unsigned long n_total_bases = 0;
unsigned long n_q20_bases = 0;
unsigned long n_q30_bases = 0;

const char *usage = "./ngs_read_qc fq1 fq2 fq3 ...\nOutput is stdout\tSimple usage example ../ngs_read_qc 1.fq 2.fq > stat.txt\n";
int main(int argc, char **argv)
{
    if(argc < 2)
    {
        printf("%s" , usage);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    char *fname;
    gzFile fp;
    kseq_t *seq;
    char *qual_str;
    int l;
    int j = 0;
    int q_len = 0;
    int qs = 0;
    double q20_perc;
    double q30_perc;
    fprintf(stdout , "#FILE\tTotal_Bases\tTotal_Reads\tQ20_Bases\tQ30_Bases\tQ20_Bases_Perc\tQ30_Based_Perc\n");
    for(i = 1; i < argc; i++)
    {
        fname = argv[i];
        fp = gzopen(fname , "r");
        seq = kseq_init(fp);
        while( ( l = kseq_read(seq) ) >= 0 )
        {
             n_total_bases += strlen(seq->seq.s);
             n_total_reads++;
             qual_str = seq->qual.s;
             q_len = strlen(seq->qual.s);
             for(j = 0; j < q_len; j++)
             {
                 qs = (int)(seq->qual.s[j]);
                 qs = qs - 33;
                 if(qs >= 20)
                 {
                     n_q20_bases++;
                     if(qs >= 30)
                     {
                         n_q30_bases++;
                     }
                 }
             }
        }
        q20_perc = (double)(n_q20_bases) / n_total_bases;
        q30_perc = (double)(n_q20_bases) / n_total_bases;
        fprintf(stdout , "%s\t%ld\t%ld\t%ld\t%ld\t%lf\t%lf\n" , fname , n_total_bases, \
                                                       n_total_reads, n_q20_bases, n_q30_bases,\
                                                       q20_perc, q30_perc);
        kseq_destroy(seq);
        gzclose(fp);
        n_total_bases = 0; n_total_reads = 0;
        n_q20_bases = 0; n_q30_bases = 0;
    }

   return 0;
}
