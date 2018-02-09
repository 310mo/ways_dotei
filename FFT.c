#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void FFT( double *xr, double *xi, double *Xr, double *Xi, int N )
{
  int i, j, k, n, n2;
  double theta, wr, wi;

  static double *rbuf, *ibuf;
  static int    bufsize = 0;

  /* memory allocation for buffers */
  if( bufsize != N ) {
    bufsize = N;
    rbuf = (double*)calloc( sizeof(double), bufsize );
    ibuf = (double*)calloc( sizeof(double), bufsize );
  }

  /* bit reverse of xr[] & xi[] --> store to rbuf[] and ibuf[] */
  i = j = 0 ;
  //バタフライ演算のために、順番を入れ替えておく
  rbuf[j] = xr[j];  ibuf[j] = xi[j];
  for( j = 1 ; j < N-1 ; j++ ) {
    for( k = N/2 ; k <= i ; k /= 2 )  i -= k;
    i += k;
    rbuf[j] = xr[i];  ibuf[j] = xi[i];
  }
  rbuf[j] = xr[j];  ibuf[j] = xi[j];

  /* butterfly calculation */
  theta = -2.0*M_PI;
  //nはforが回るごとに２倍ずつされていく、バタフライ演算を行う段数だと考えてよい？
  for( n = 1 ; ( n2 = n*2 ) <= N ; n = n2 ) {
    //バタフライ演算を複数回行う時に、角度は1/2になっていく（１段目、２段目、３段目、、、ってやっていく感じ）
    theta *= 0.5;
    for ( i = 0 ; i < n ; i++ ) {
      wr = cos(theta*i);  wi = sin(theta*i);
      //rbufですでに順番をいい感じに入れ替えているので、隣のこうと計算を行えばいいようになっている
      for ( j = i ; j < N ; j += n2 ) {
        k = j + n;
        /*
          X[j] = 1*buf[j] + W*buf[k];
          X[k] = 1*buf[j] - W*buf[k];
          Note : X[], buf[], and W are complex numbers.
          Re{ X[n] } = Xr[n], Im{ X[n] } = Xi[n];
          Re{ buf[n] } = rbuf[n], Im{ buf[n] } = ibuf[n];
          Re{ W } = wr, Im{ W } = wi;
        */
        Xr[j] = rbuf[j] + (rbuf[k]*wr-ibuf[k]*wi);  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xi[j] = ibuf[j] + (rbuf[k]*wi+ibuf[k]*wr);  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xr[k] = rbuf[j] - (rbuf[k]*wr-ibuf[k]*wi);  /* ??????????, using wr, wi, rbuf, and ibuf */
        Xi[k] = ibuf[j] - (rbuf[k]*wi+ibuf[k]*wr);  /* ??????????, using wr, wi, rbuf, and ibuf */
      }
    }

    for( i = 0 ; i < N ; i++ ) {
      rbuf[i] = Xr[i];
      ibuf[i] = Xi[i];
    }
  }
  return;
}

void IFFT( double *Xr, double *Xi, double *xr, double *xi, int N )
{
  int i, j, k, n, n2;
  double theta, wr, wi;

  static double *rbuf, *ibuf;
  static int    bufsize = 0;

  /* memory allocation for buffers */
  if( bufsize != N ) {
    bufsize = N;
    rbuf = (double*)calloc( sizeof(double), bufsize );
    ibuf = (double*)calloc( sizeof(double), bufsize );
  }

  /* bit reverse of Xr[] & Xi[] --> store to rbuf[] and ibuf[] */
  i = j = 0 ;
  rbuf[j] = Xr[j]/N;  ibuf[j] = Xi[j]/N;
  for( j = 1 ; j < N-1 ; j++ ) {
    for( k = N/2 ; k <= i ; k /= 2 )  i -= k;
    i += k;
    rbuf[j] = Xr[i]/N;  ibuf[j] = Xi[i]/N;
  }
  rbuf[j] = Xr[j]/N;  ibuf[j] = Xi[j]/N;

  /* butterfly calculation */
  //ifftはthetaの値の正負がひっくり返っただけ
  theta = 2.0*M_PI;  /* not -2.0*M_PI !!! */
  for( n = 1 ; ( n2 = n*2 ) <= N ; n = n2 ) {
    theta *= 0.5;
    for ( i = 0 ; i < n ; i++ ) {
      wr = cos(theta*i);  wi = sin(theta*i);
      for ( j = i ; j < N ; j += n2 ) {
      k = j + n;
      /*
        x[j] = 1*buf[j] + W*buf[k];
        x[k] = 1*buf[j] - W*buf[k];
        Note : x[], buf[], and W are complex numbers.
        Re{ x[n] } = xr[n], Im{ x[n] } = xi[n];
        Re{ buf[n] } = rbuf[n], Im{ buf[n] } = ibuf[n];
        Re{ W } = wr, Im{ W } = wi;
      */
      xr[j] = rbuf[j] + (rbuf[k]*wr-ibuf[k]*wi);   /* ??????????, using wr, wi, rbuf, and ibuf */
      xi[j] = ibuf[j] + (rbuf[k]*wi+ibuf[k]*wr);   /* ??????????, using wr, wi, rbuf, and ibuf */
      xr[k] = rbuf[j] - (rbuf[k]*wr-ibuf[k]*wi);   /* ??????????, using wr, wi, rbuf, and ibuf */
      xi[k] = ibuf[j] - (rbuf[k]*wi+ibuf[k]*wr);   /* ??????????, using wr, wi, rbuf, and ibuf */
      }
    }

    for( i = 0 ; i < N ; i++ ) {
      rbuf[i] = xr[i];
      ibuf[i] = xi[i];
    }
  }
  return;
}

//コマンドライン引数はDFTの時と同じ、１引数目はファイル名、２引数目は開始する地点、３引数目は窓長
int main( int argc, char *argv[] )
{
  int  i, nskip, framelen;
  double spec;
  short  *sdata;
  double *xr, *xi, *Xr, *Xi;
  FILE *fpDAT;
  
  if( 4 != argc ) {
    fprintf( stderr, "Usage : %s DATfile skip[sample] frame_length[sample]\n", argv[0] );
    exit( 1 );
  }
  int datanum = atoi(argv[3]);
  int num = 1;
  while(num<datanum) {
    num *= 2;
  }
  num /= 2;

  if( NULL == ( fpDAT = fopen( argv[1], "r" ) ) )  exit( 1 );
  if( 0 > ( nskip    = atoi( argv[2] ) ) )  exit( 1 );
  if( 0 > ( framelen = atoi( argv[3] ) ) )  exit( 1 );
  //framelen = num;

  fprintf( stderr, "# DATfile = %s\n", argv[1] );
  fprintf( stderr, "# %d samples are skipped.\n", nskip );
  fprintf( stderr, "# 1 frame contains %d sampels.\n", framelen );

  /* memory allocation */
  sdata = (short*)calloc( sizeof(short), framelen );
  xr = (double*)calloc( sizeof(double), framelen );
  xi = (double*)calloc( sizeof(double), framelen );
  Xr = (double*)calloc( sizeof(double), framelen );
  Xi = (double*)calloc( sizeof(double), framelen );
  if( NULL == sdata || NULL == xi || NULL == xi || NULL == Xr || NULL == Xi )  exit( 1 );

  fseek( fpDAT, nskip*sizeof(short), SEEK_SET );
  fread( sdata, sizeof(short), framelen, fpDAT );
  fclose( fpDAT );

  for( i = 0 ; i < framelen ; i++ ) {
    xr[i] = (0.54-0.46*cos(2*M_PI*i/(framelen-1)))*sdata[i];
    xi[i] = 0.0;
  }

  /*for( i = 0 ; i < framelen ; i++ ) {
    Xr[i] = sdata[i];
    Xi[i] = 0.0;
  }*/

  /* for( i = 0 ; i < 1000 ; i++ ) { */
  FFT( xr, xi, Xr, Xi, framelen );
  /* } */

  for( i = 0 ; i < framelen ; i++ ) {
    Xr[i] = log10( (1.0/framelen)*( Xr[i]*Xr[i]+Xi[i]*Xi[i] ) );
    Xi[i] = 0;
  }

  IFFT( Xr, Xi, xr, xi, framelen );

  /*
  for( i = 0 ; i < framelen ; i++ ) {
    if(i<30 || i>framelen -30) {
      xr[i] = xr[i];
      xi[i] = xi[i];
    }
    else {
      xr[i] = 0.0;
      xi[i] = 0.0;
    }
  }
  */

  for( i = 0 ; i < framelen ; i++ ) {
    printf( "%d %lf\n", i, xr[i] );
  }


  /*FFT( xr, xi, Xr, Xi, framelen );

  for( i = 0 ; i < framelen ; i++ ) {
    printf( "%d %lf\n", i, Xr[i] );
  }*/


  return( 0 );
}
