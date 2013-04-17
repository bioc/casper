#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include "header.h"

//********** Functions ***********

int *procCigar(char *cigar, int *cigs){

  char *num;
  cigs[0]=0;
  num=malloc(10*sizeof(int));
  strcpy(num, "\0");
  while(*cigar != '\0'){
    switch(*cigar)
      {
      case 'M': 
	sscanf(num, "%d", &cigs[cigs[0]+1]);
	cigs[0]++;
	strcpy(num, "\0");
	break;

      case 'D': case 'S': case 'P': case 'H': case 'N':
	sscanf(num, "%d", &cigs[cigs[0]+1]);
	cigs[cigs[0]+1]*=-1;
	*num='\0';
	cigs[0]++;
        break;   
	
      case 'I': break;
      
      default:
	strncat(num, cigar, 1);
	break;
      }
    cigar++;
  }
  
  //free(cigar);
  free(num);
  return(cigs);
}
	


/*char *pch, *Mtab="M", *Stab="IDSHP";
  pch=strtok(cigar, tab);
  cigs[0]=0;
  while(pch!=NULL){
    sscanf(pch, "%d", &cigs[cigs[0]+1]);
    pch=strtok(NULL, tab);
    cigs[0]++;
  }
  free(cigar);
  return(cigs);
  }*/


void addRead2Frag(const char *qname, const char *chr, int start, int strand, const char *cigar, int totF, read_t *frags, int read){
	if(read==1){
	  frags[totF].qname = malloc((strlen(qname)+1) * sizeof(char));
	  strcpy(frags[totF].qname, qname);
	  frags[totF].chr_1 = malloc((strlen(chr)+1) * sizeof(char));
	  strcpy(frags[totF].chr_1, chr);
	  frags[totF].st_1=start;
	  frags[totF].cigar_1 = malloc((strlen(cigar)+1) * sizeof(char));
	  strcpy(frags[totF].cigar_1, cigar);
	  frags[totF].strand_1=strand;
	  frags[totF].nreads=1;
	} else {    
	  frags[totF].chr_2 = malloc((strlen(chr)+1) * sizeof(char));
	  strcpy(frags[totF].chr_2, chr);
	  frags[totF].st_2=start; 
	  frags[totF].cigar_2 = malloc((strlen(cigar)+1) * sizeof(char));
	  strcpy(frags[totF].cigar_2, cigar);
	  frags[totF].strand_2 = strand;
	  frags[totF].nreads=2;
	}   
}
