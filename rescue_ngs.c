/*
 * rescue_ngs.c
 * Rescue undetermined NGS data using index information
 *
 * Code written in C due to fast computing speed
 * Author : Sunghoon Heo
 * 
 * This code requires Heng Li's kseq.h C library for fastq parsing for safety
 * Also requres mg_string library for split string and mg_vector for storage system
 */


#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdlib.h>
#include "kseq.h"

#include "kseq.h"
// Remark
// NGS index combination is always p7 + p5 reverse complement
char *fwd_index = NULL;
char *rev_index = NULL;
int use_rc = 1;
int ed = 1;
char *oprefix = NULL;
char *fwd_fname = NULL;
char *rev_fname = NULL;
char rc_index[9] = {0,};
int strict = 0;
int single_mode = 0;
const char *usage = "./rescue [options] undet1.fq <undet2.fq...>\n"
"[Options]\n"
"    -e     INT    tolerate untill INT mismatches. Default = 1\n"
"    -f     STR    p5_index(fwd in D.Bang's Lab). Required. Always must be set\n"
"    -r     STR    p7_index(rev in D. Bang's lab) \n"
"    -u     BOOL   If set, will not apply reverse complement to reverse illumina index ( dfault = ON )\n"
"    -o     STR    Output prefix. Deault : <fwd_index>-<rev_index>_1.fq and 2.fq\n"
"    -s     BOOL   Strictly compare index information. Default = NO\n"
"    -S     BOOL   Use single index\n";

int hamming(char *s1, char *s2)
{
    assert(strlen(s1) == strlen(s2));
    size_t s1l = strlen(s1);
    size_t idx = 0;
    int diff = 0;
    for(idx = 0; idx < s1l; idx++)
    {
        if(s1[idx] != s2[idx]) diff++;
    }
    return diff;
}


void rc_comp(char *index)
{
    int i = 0;
    char base;
    for( i = 7 ; i != -1 ; i--)
    {
        base = index[i];
        switch(base)
        {
            case 'A': rc_index[7-i] = 'T'; break;
            case 'T': rc_index[7-i] = 'A'; break;
            case 'G': rc_index[7-i] = 'C'; break;
            case 'C': rc_index[7-i] = 'G'; break;
            default : rc_index[7-i] = 'N'; break;
        }
    }
    rc_index[strlen(index)] = '\0';
}
/* Pre-defined functions */

int main_perform_two_sample(int argc, char** argv, int o);
int main_perform_single(int argc , char **argv, int i);

KSEQ_INIT(gzFile, gzread)
int main(int argc , char **argv)
{
	 if(argc < 2)
	 {
		   printf("%s\n" , usage);
		   exit(EXIT_FAILURE);
	 }
    int c;
    int help = 0;
    while( ( c = getopt(argc, argv, "husSf:r:e:o:" )) != -1 )
    {
        switch(c)
        {
            case 'h':
                help = 1;
                break;
            case 'u':
                use_rc = 0;
                break;
            case 'f':
                fwd_index = strdup(optarg);
                break;
            case 'r':
                rev_index = strdup(optarg);
                break;
            case 'e':
                 ed = atoi(optarg);
                 break;
            case 'o':
                 oprefix = strdup(optarg);
                 break;
            case 's':
                 strict = 1;
                 break;
		    case 'S' :
			      single_mode = 1;
				  break;
        }
    }

    if(help)
    {
        printf("%s\n" , usage);
        if(fwd_index != NULL) free(fwd_index);
        if(rev_index != NULL) free(rev_index);
        exit(EXIT_FAILURE);
    }
	int ret = 0;
	printf("fwd index = %s\n" , fwd_index);
    if( single_mode == 0 )
	{
		printf("-S ( single mode ) is not set. \nPerform In PAIRED mode\n");
		ret = main_perform_two_sample(argc, argv, optind);
	}
    else
	{
		printf("-S ( single mode ) is set. \nPerform In SINGLE mode\n");
		ret = main_perform_single(argc , argv, optind);
	}
    return ret;
}

int main_perform_two_sample(int argc, char ** argv, int oid)
{
	
    fwd_fname = argv[oid];
    rev_fname = argv[oid + 1];

    printf("[CommandLine] rescue_ngs ");
    if(use_rc)
    {
        printf("-u ");
    }
    char *fwd_oname;
    char *rev_oname;
    size_t default_oprefix_len = strlen(fwd_index) + strlen(rev_index) + strlen("-") + 1;
    int nreads = 0;
    int nrescued = 0;
    if(oprefix == NULL)
    {
        oprefix = (char*)malloc(default_oprefix_len);
        strncpy(oprefix , fwd_index , strlen(fwd_index));
        strncpy(oprefix + strlen(fwd_index) , "-" , strlen("-"));
        strncpy(oprefix + strlen(fwd_index) + strlen("-"), rev_index , strlen(rev_index));
        oprefix[strlen(oprefix) + strlen(fwd_index) + strlen("-") + strlen(rev_index)] = '\0';
    }

    fwd_oname = (char*)malloc(strlen(oprefix) + strlen("_1.fq") + 1);
    rev_oname = (char*)malloc(strlen(oprefix) + strlen("_2.fq") + 1);

    strncpy(fwd_oname , oprefix, strlen(oprefix));
    strncpy(rev_oname , oprefix, strlen(oprefix));

    strncpy(fwd_oname + strlen(oprefix) , "_1.fq" , strlen("_1.fq"));
    strncpy(rev_oname + strlen(oprefix) , "_2.fq" , strlen("_2.fq"));
	
    fwd_oname[strlen(oprefix) + strlen("_1.fq")] = '\0';
    rev_oname[strlen(oprefix) + strlen("_2.fq")] = '\0';
	
    printf("-o %s ", oprefix);
    printf("-f %s -r %s -e %d ",fwd_index, rev_index, ed);
    printf("%s %s\n\n" , fwd_fname, rev_fname);

	printf("Forward Output name = %s\n", fwd_oname);
	printf("Reverse Output name = %s\n" , rev_oname);
    gzFile fwd_fp; gzFile rev_fp;
    fwd_fp = gzopen(fwd_fname , "r"); rev_fp = gzopen(rev_fname , "r");

    kseq_t *fwd_seq;
    kseq_t *rev_seq;


    fwd_seq = kseq_init(fwd_fp);
    rev_seq =kseq_init(rev_fp);
    FILE *ofp1 = fopen(fwd_oname , "w");
    FILE *ofp2 = fopen(rev_oname , "w");

    // Main algorithm
    int l1 , l2;
    mg_str_vec vec1;
    mg_str_vec vec2;
    char *p5_index = NULL;
    char *p7_index = NULL;
    if(use_rc)
    {
        rc_comp(fwd_index);
        p5_index = rc_index;
    }
    else
    {
        p5_index = fwd_index;
    }

    p7_index = rev_index;

    char *name1 = NULL; char *name2 = NULL;
    char *q1 = NULL; char* q2 = NULL;
    int ham1_1 = 0; int ham2_1 = 0; int ham1_2= 0; int ham2_2 = 0;
    clock_t s_t , e_t;
    s_t = clock();
	char *tmp = NULL;
	char *index_comb1 = NULL;
	char *index_comb2 = NULL;
	mg_str_vec v1;
	mg_str_vec v2;
	
    while( ( ( l1 = kseq_read(fwd_seq) ) >= 0 ) && ( ( l2 = kseq_read(rev_seq) ) >= 0 ) )
    {
        nreads++;
		if(nreads % 1000000 == 0)
		{
			printf("Processed %d reads\n" , nreads);
		}
        q1 = fwd_seq->name.s;
        q2 = rev_seq->name.s;
        if(strcmp(q1,q2) != 0 )
        {
            continue; // Illumina sequencing reads must have same id line by line
        }
		mg_init_vec(&vec1 , 's');
		mg_init_vec(&vec2 , 's');
		mg_init_vec(&v1 , 's');
		mg_init_vec(&v2 , 's');
        name1 = fwd_seq->comment.s;
        name2 = rev_seq->comment.s;

        mg_split2(name1 , ":", &vec1);
        mg_split2(name2 , ":", &vec2);
        index_comb1 = mg_svec_at(&vec1 , 3);
		index_comb2 = mg_svec_at(&vec2 , 3);
        // Compute hamming distance
        // p7-p5
		mg_split2(index_comb1 , "+" , &v1);
		mg_split2(index_comb2 , "+" , &v2);
		tmp = mg_svec_at(&v1,0);
        ham1_1 = hamming(tmp , p7_index);
		tmp = mg_svec_at(&v1,1);
        ham1_2 = hamming(tmp , p5_index);
		tmp = mg_svec_at(&v2,0);
        ham2_1 = hamming(tmp , p7_index);
		tmp = mg_svec_at(&v1,1);
        ham2_2 = hamming(tmp , p5_index);
		tmp = NULL;
        if(strict)
        {
             if( (ham1_1 <= ed ) && (ham1_2 <= ed) && (ham2_1 <= ed ) && (ham2_2 <= ed ) )
             {
                 nrescued++;
                 fprintf(ofp1 , "@%s %s\n%s\n+\n%s\n" , q1 , name1, fwd_seq->seq.s, fwd_seq->qual.s);
                 fprintf(ofp2 , "@%s %s\n%s\n+\n%s\n" , q2 , name2 ,rev_seq->seq.s , rev_seq->qual.s);
				 mg_clear_vec(&vec1 ,'s');
                 mg_clear_vec(&vec2, 's');
		
	           	mg_clear_vec(&v1 , 's'); mg_clear_vec(&v2, 's');
             }
             else
             {
				 mg_clear_vec(&vec1 ,'s');
                 mg_clear_vec(&vec2, 's');
		
            		mg_clear_vec(&v1 , 's'); mg_clear_vec(&v2, 's');
                 continue;
             }
        }
        else
        {
             if ( ( ( ham1_1 <= ed ) || ( ham1_2 <= ed ) ) && ( ( ham2_1 <= ed ) || ( ham2_2 <= ed ) ) )
             {
                 nrescued++;
                 fprintf(ofp1 , "@%s %s\n%s\n+\n%s\n" , q1 , name1, fwd_seq->seq.s, fwd_seq->qual.s);
                 fprintf(ofp2 , "@%s %s\n%s\n+\n%s\n" , q2 , name2 ,rev_seq->seq.s , rev_seq->qual.s);
				 mg_clear_vec(&vec1 ,'s');
                mg_clear_vec(&vec2, 's');
		
	          	mg_clear_vec(&v1 , 's'); mg_clear_vec(&v2, 's');
             }
             else
             {
				 mg_clear_vec(&vec1 ,'s');
                 mg_clear_vec(&vec2, 's');
		
		         mg_clear_vec(&v1 , 's'); mg_clear_vec(&v2, 's');
                 continue;
             }
        }
        // strict
        // All hamming distance must be smaller than setting value
      
    }
    e_t = clock();
    printf("For index combination (p5,p7) = (%s , %s), from %d reads %d reads rescued\n" , fwd_index, rev_index, nreads, nrescued);
    printf("%lf sec elapsed for processing...\n" , (double)(e_t - s_t) / CLOCKS_PER_SEC);
    fclose(ofp1);
    fclose(ofp2);
    kseq_destroy(fwd_seq);
    kseq_destroy(rev_seq);
    gzclose(fwd_fp);
    gzclose(rev_fp);
    free(fwd_index); free(rev_index); free(oprefix);
    free(fwd_oname); free(rev_oname);

//    printf("For index combination (p5,p7) = (%s , %s), from %d reads %d reads rescued\n" , fwd_seq, rev_seq, nreads, nrescued);
    return 0;
}

int main_perform_single(int argc , char **argv, int oid)
{
	fwd_fname = argv[oid];

    printf("[CommandLine] rescue_ngs ");
    if(use_rc)
    {
        printf("-u ");
    }
    char *fwd_oname;
    
    size_t default_oprefix_len = strlen(fwd_index) + 1;
    int nreads = 0;
    int nrescued = 0;
    if(oprefix == NULL)
    {
        oprefix = (char*)malloc(default_oprefix_len);
        strncpy(oprefix , fwd_index , strlen(fwd_index));
        oprefix[strlen(oprefix)] = '\0';
    }

    fwd_oname = (char*)malloc(strlen(oprefix) + strlen(".fq") + 1);
    

    strncpy(fwd_oname , oprefix, strlen(oprefix));
    

    strncpy(fwd_oname + strlen(oprefix) , ".fq" , strlen(".fq"));
    
	
    fwd_oname[strlen(oprefix) + strlen("_1.fq")] = '\0';
    
	
    printf("-o %s ", oprefix);
    printf("-f %s -e %d ",fwd_index, ed);
    printf("%s \n\n" , fwd_fname);

	printf("Forward Output name = %s\n", fwd_oname);
	
    gzFile fwd_fp; 
    fwd_fp = gzopen(fwd_fname , "r");

    kseq_t *fwd_seq;



    fwd_seq = kseq_init(fwd_fp);

    FILE *ofp1 = fopen(fwd_oname , "w");


    // Main algorithm
    int l1 ;
    mg_str_vec vec1;
    
    char *p5_index = NULL;
    
    if(use_rc)
    {
        rc_comp(fwd_index);
        p5_index = rc_index;
    }
    else
    {
        p5_index = fwd_index;
    }

    
	
    char *name1 = NULL; 
    char *q1 = NULL; 
    int ham1_1 = 0;
    clock_t s_t , e_t;
    s_t = clock();
	char *tmp = NULL;
	char *index_comb1 = NULL;
	
	mg_str_vec v1;
	
	
    while( ( ( l1 = kseq_read(fwd_seq) ) >= 0 ) )
    {
        nreads++;
		if(nreads % 1000000 == 0)
		{
			printf("Processed %d reads\n" , nreads);
		}
        q1 = fwd_seq->name.s;
        
        
		mg_init_vec(&vec1 , 's');
		
		mg_init_vec(&v1 , 's');
		
        name1 = fwd_seq->comment.s;
        
        mg_split2(name1 , ":", &vec1);
        
        index_comb1 = mg_svec_at(&vec1 , 3);
		
		
        ham1_1 = hamming(index_comb1 , p5_index);
		
		tmp = NULL;
        if(strict)
        {
             if( (ham1_1 <= ed ) )
             {
                 nrescued++;
                 fprintf(ofp1 , "@%s %s\n%s\n+\n%s\n" , q1 , name1, fwd_seq->seq.s, fwd_seq->qual.s);
                 
				 mg_clear_vec(&vec1 ,'s');
                 
		
	           	mg_clear_vec(&v1 , 's'); 
             }
             else
             {
				 mg_clear_vec(&vec1 ,'s');
                 
		
            		mg_clear_vec(&v1 , 's'); 
                 continue;
             }
        }
        else
        {
             if ( ( ham1_1 <= ed ) )
             {
                 nrescued++;
                 fprintf(ofp1 , "@%s %s\n%s\n+\n%s\n" , q1 , name1, fwd_seq->seq.s, fwd_seq->qual.s);
                 
				 mg_clear_vec(&vec1 ,'s');
                
		
	          	mg_clear_vec(&v1 , 's');
             }
             else
             {
				 mg_clear_vec(&vec1 ,'s');
                 
		
		         mg_clear_vec(&v1 , 's');
                 continue;
             }
        }
        // strict
        // All hamming distance must be smaller than setting value
      
    }
    e_t = clock();
    printf("For index combination (p5,p7) = (%s , %s), from %d reads %d reads rescued\n" , fwd_index, rev_index, nreads, nrescued);
    printf("%lf sec elapsed for processing...\n" , (double)(e_t - s_t) / CLOCKS_PER_SEC);
    fclose(ofp1);
    
    kseq_destroy(fwd_seq);
    
    gzclose(fwd_fp);
    
    free(fwd_index); free(oprefix);
    free(fwd_oname); 

//    printf("For index combination (p5,p7) = (%s , %s), from %d reads %d reads rescued\n" , fwd_seq, rev_seq, nreads, nrescued);
    return 0;
}
