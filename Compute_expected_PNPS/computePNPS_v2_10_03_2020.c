//1a versio lliurada a la Mar el 1 Mar 2017
//2a versio mes documentada
#define MAXLONGNAME 1000
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int num_syn_cod[64] = {1,1,1,1,3,3,3,3,5,1,
                        5,1,2,2,0,2,1,1,1,1,
                        3,3,3,3,5,5,5,5,5,5,
                        5,5,1,1,1,1,3,3,3,3,
                        3,3,3,3,3,3,3,3,2,1,
                        2,1,3,3,3,3,2,1,0,1,
                        5,1,5,1};
//num_syn_cod: number of synonymous codons: codon 0 has 1 synonymous, codon 1 has one synonymous, ...

const int syn_cod[64][64]={{2},{3},{0},{1},{5,6,7},{4,6,7},{4,5,7},{4,5,6},{10,24},{11},
			 {8,26},{9},{13,15},{12,15},{},{12,13},{18},{19},{16},{17},
                         {21,22,23},{20,22,23},{20,21,23},{20,21,22},{8,25,26,27},{24,26,27},{10,24,25,27},{24,25,26},{29,30,31,60},{28,30,31},
                         {28,29,31,62},{28,29,30},{34},{35},{32},{33},{37,38,39},{36,38,39},{36,37,39},{36,37,38},
                         {41,42,43}, {40,42,43},{40,41,43},{40,41,42},{45,46,47},{44,46,47},{44,45,47},{44,45,46},{50,56},{51},
                         {48,56},{49},{53,54,55},{52,54,55},{52,53,55},{52,53,54},{48,50},{59},{},{57},
                         {28,62},{63},{30,60},{61}};
//syn_cod: Matrix of synonymous codons: codon 0 (AAA) is synonymous  of codon 2 (AAG), codon 1 (AAC)  is synonymous  of codon 3 (AAG),...

void compute_prob_from_transitions(long int transitions[4][4],float probmat[4][4]){
  long int sum,i,j;
  float sumf;
printf("\nMatrix of transitions:\n"); 
  for (i=0;i<4;i++){for (j=0;j<4;j++) printf("%8ld  ",transitions[i][j]);printf("\n");}

  sum=0;
  for (i=0;i<4;i++) for (j=0;j<4;j++)    sum += transitions[i][j]; printf("%ld \n",sum);
  for (i=0;i<4;i++) for (j=0;j<4;j++)   if (i!=j) probmat[i][j]=((float)transitions[i][j])/sum; else probmat[i][j]=0;

  for (i=0;i<4;i++){
    sumf=0; for (j=0;j<4;j++) sumf += probmat[i][j];
    probmat[i][i]=1-sumf;
  } 
  printf("\nProbability Matrix:\n"); 
  for (i=0;i<4;i++){for (j=0;j<4;j++) printf("%1.5f  ",probmat[i][j]);printf("\n");}printf("\n");
}

int valor(char c){
  int r;
  switch (c){
  case 'A':r=0;break; 
  case 'C':r=1;break; 
  case 'G':r=2;break; 
  case 'T':r=3;break;
  case 'a':r=0;break; 
  case 'c':r=1;break; 
  case 'g':r=2;break; 
  case 't':r=3;break;
  }
  return r;
}
//check if codgen is synonymous of  codon
int is_synonymous(int codon,int codgen){
  int i=0;
  while(syn_cod[codon][i]!=codgen && i<num_syn_cod[codon]){
    //printf("syn_cod[%d][%d]=%d\n",base,i,syn_cod[base][i]);
    i++;
  }
  //printf("hi_es %d = %d\n",codon,i<num_syn_cod[base]);
  return i<num_syn_cod[codon];
 
}
//Computes the value of PS and PN for each codon and registers them in tpn and  tps vectors
void compute_prob_aminos(float *tpn, float *tps,float mat[4][4]){
  int i,codon,exp,base,vcod,codgen,pos;
  float probps,probpn;
  for (codon=0;codon<64;codon++){//recall that a codon is defined by its value 0 (AAA),1 (AAC) ,2 (AAG),...
    probpn=0;
    probps=0;
    exp=1;
    for (pos=0;pos<3;pos++){//for all three bases of a codon
      base = (codon/exp)%4;//the base in position pos is determined
      vcod=codon-base*exp;
      //printf("codon=%d pos=%d base=%d  vcod=%d\n",codon,pos,base,vcod);
      for (i=0;i<4;i++) //genero els codons mutats 
	if (i!=base){
	  codgen=vcod+i*exp;
	  //printf("is_synonymous(%d,%d)=%d  mat[%d][%d]=%1.3f\n",codon,codgen,is_synonymous(codon,codgen),base,i,mat[base][i]);
	  if (is_synonymous(codon,codgen)) probps += mat[base][i];
	      else probpn += mat[base][i];
	  //printf("                                    probpn=%1.3f probps=%1.3f\n",probpn,probps);
	}
      exp *=4;
    }
    tpn[codon]=probpn;
    tps[codon]=probps;
    //printf("                         pn=%1.3f ps=%1.3f\n",probpn,probps);
  }  
  //make prob. of stop codons to  0
  tpn[48]=0;tpn[50]=0;tpn[56]=0;
  tps[48]=0;tps[50]=0;tps[56]=0;
  //for (codon=0;codon<64;codon++) printf("%4d %1.3f %1.3f  %2.4f \n",codon,tpn[codon],tps[codon],tpn[codon]/tps[codon]);
}

void read_transitions(char fname[MAXLONGNAME],long int transitions[4][4]){
    int i,j;
    FILE *f;
    if(!(f = fopen(fname,"r"))){printf("error: opening data file of transitions\n"); exit(-1);}
    for (i=0;i<4;i++) for (j=0;j<4;j++) fscanf(f,"%ld",&transitions[i][j]);
    fclose(f);
  }

//Llegeix cada sequencia de tres en tres bases, mira a la taula quin PS i PN tenen i ho va sumant
void main(int argc, char * argv []){

  FILE *fseq;
  char c;
  int codon,i,ncodon,error;
  long int transitions[4][4];
  float tps[64],tpn[64],ps,pn;
  char nom[1000];
  float probmat[4][4];
  if (argc==1){
    printf("This software computes, for each codon, the probability of  be sinonymous when only one base is mutated, and then computes the mean value for all codons of a sequence.\n \n");
    printf("\t computePNPS File_with_transitions\n");
    printf("\t computePNPS File_with_transitions DNA_fasta_file\n\n");
    printf("The output is composed by the transitions matrix, the probability matrix and ");
    printf(" a table with the following columns:\n");
    printf("\t identifier of the sequence\n");
    printf(" \t PN\n");
    printf(" \t PS\n");
    printf(" \t PN/PS\n \t PN/(PN+PS)*length of the sequence\n");
    printf(" \t PS/(PN+PS)*length of the sequence\n");
    exit(-1);
  }
  read_transitions(argv[1],transitions);
  compute_prob_from_transitions(transitions,probmat);
  compute_prob_aminos(tpn,tps,probmat);
  if (argc == 3){
    if(!(fseq = fopen(argv[2],"r"))){
      printf("error: opening data file %s  \n",argv[1]);
      exit(-1);
    }
    
    fscanf(fseq,"%c",&c);
    while (!feof(fseq)){
      ncodon=0;
      ps=0;pn=0;
      while (!feof(fseq) && c!='>') fscanf(fseq,"%c",&c);
      fscanf(fseq,"%c",&c);   
      i=0;while (!feof(fseq) && c!='\n') {nom[i]=c;i++;fscanf(fseq,"%c",&c);}
      nom[i]='\0';
      fscanf(fseq,"%c",&c);
      while (!feof(fseq) && c!='\n'){
	//llegir codon
	codon=valor(c);
	error=0;
	for (i=0;i<2 && !error;i++){
	  fscanf(fseq,"%c",&c);
	  if (c != '\n') codon=codon*4+valor(c);
	  else error=1;
	}
	//printf("codon %d ps%f pn%f\n",codon,tpn[codon],tps[codon]);
	if (!error){
	  ps +=tps[codon];
	  pn +=tpn[codon];
	  ncodon ++;
	  fscanf(fseq,"%c",&c);
	}
      }  
      if(!feof(fseq))
	if (!error) printf("%s \t%3.6f  \t%3.6f \t%3.6f \t%ld \t%ld \n",nom,pn/ncodon,ps/ncodon,pn/ps,(long int)(pn/(pn+ps)*ncodon*3),(long int)(ps/(pn+ps)*ncodon*3));
	else printf("%s: cannot be analyzed\n",nom);
      //if(!feof(fseq))printf("\t%s \t%3.6f  \t%3.6f \t%3.6f \t%ld \t%ld \t%ld \n",nom,pn/ncodon,ps/ncodon,pn/ps,(long int)(pn/(pn+ps)*ncodon*3),(long int)(ps/(pn+ps)*ncodon*3),(long int)((pn/(pn+ps)*ncodon*3)+(ps/(pn+ps)*ncodon*3)) );
    }
  }
}
