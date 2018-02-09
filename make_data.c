#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define N 1000

int main(int argc, char **argv) {
  int bitcount = 0;
  short data[N];
  char datacheck[5] = "sato";
  FILE *fp;
  FILE *fq;
  char *filename = argv[1];
  char *filename2 = argv[2];

  if(argc!=3) {
    printf("failed to do\n");
    return 0;
  }
  if((fp=fopen(filename, "rb"))==NULL) {
    printf("failed to open\n");
    return 0;
  }
  if((fq=fopen(filename2, "w"))==NULL) {
    printf("failed to open\n");
    return 0;
  }
  while(1) {
    if(strcmp(datacheck, "data")==0) {
      break;
    }
    fread(datacheck, 4, 1, fp);
    datacheck[4] = '\0';
    printf("%s\n", datacheck);
    bitcount += 1;
    fseek(fp, bitcount, SEEK_SET);
  }
  printf("bitcount = %d\n", bitcount);
  fseek(fp, 7+bitcount, SEEK_SET);
  int n = 100;
  while(n!=0) {
     n = fread(data, 2, N, fp);
     if(n!=0) {
       fwrite(data, 2, n, fq);
     }
  }

  fclose(fp);
  fclose(fq);
  return 0;
}
