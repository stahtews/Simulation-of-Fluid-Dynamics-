#include <iostream>
#include "Timer.h"
#include <stdlib.h> // atoi
#include <math.h>
#include <mpi.h>
#include <omp.h>
//constants for Schroedinger’s formula and java program 
int default_size = 100; // the default system size
int defaultCellWidth = 8 ;
double c = 1.0; // wave speed
double dt = 0.1; // time quantum
double dd = 2.0; // change in system
using namespace std;
//------------------Schroedinger’s formula calculation when t =1 z[1]-------------------------
double ZtOne( double center,double down, double up, double right,double left)
 {
 double z ;
 double second ;
 double third;
 double fourth;
 double internal;
 second = pow(c,2)/2;
 internal = 4 * center;
 third = ( down + up + right + left - internal );
 fourth =(double) pow(dt/dd,2);
 z = center + (second * fourth * third) ;
 return z;
 }
//------------------Schroedinger’s formula calculation when t =1 z[1]-------------------------
 double ZtTwo(double previous, double center, double down,double up, double right,double left)
 {
 double first;
 double finalPart;
 double internal;
 double last;
 double constantPart;
 first = 2 * center;
 constantPart = pow(dt/dd,2) ;
 internal = 4 * center;
 last = down + up + right + left - internal;
 finalPart =(first - previous + (constantPart * last));

 return finalPart;
 }
//-----------------Main-----------------------------------------------------------------------
int main( int argc, char *argv[] ) {
 MPI_Request request;
 MPI_Status mStatus;
 int my_rank ;
 my_rank = 0; // used by MPI
 int mpi_size;
 mpi_size= 1; // used by MPI
 // verify arguments
 if ( argc != 5 ) {
 cerr << "usage: Wave2D size max_time interval threads" << endl;
 return -1;
 }
 //get the arguments from command line
 int size = atoi( argv[1] );
 int max_time = atoi( argv[2] );
 int interval = atoi( argv[3] );
 int nThreads;
 nThreads = 1;
 nThreads = atoi(argv[4]);
 omp_set_num_threads(nThreads);

 // display errors when number of arguments is not right
 if ( size < default_size || max_time < 3 || interval < 0 ) {
 cerr << "usage: Wave2D size max_time interval" << endl;
 cerr << " where size >= 100 && time >= 3 && interval >= 0" << endl;
 return -1;
 }
 // create a simulation space
 // open mp parallelization for for loop
 double z[3][size][size];
 #pragma omp parallel for
 for ( int p = 0; p < 3; p++ ) 
 for ( int i = 0; i < size; i++ )
 for ( int j = 0; j < size; j++ )
z[p][i][j] = 0.0; // no wave
 // start a timer
 Timer time;
 time.start( );
 // time = 0;
 // initialize the simulation space: calculate z[0][][]
 // open mp parallelization for for loop
 int weight = size / default_size;
 #pragma omp parallel for
 for( int i = 0; i < size; i++ ) {
 for( int j = 0; j < size; j++ ) {

 if( (i > 40 * weight) && (i < 60 * weight) && (j > 40 * weight) && (j < 60 * weight) ) {
 z[0][i][j] = 20.0;
 } else {
 z[0][i][j] = 0.0;
 }
 }
 }
 //cout << "this was z[0]" << endl;
 // time = 1
 // calculate z[1][][] 
 // cells not on edge
 #pragma omp parallel for
 for( int i = 0; i < size; i++ ) {
 for( int j = 0; j < size; j++ ) {
 if( (i > 0) && (i < (size-1)) && (j > 0) && (j < (size-1)) ) {
 z[1][i][j] = ZtOne(z[0][i][j], z[0][i+1][j], z[0][i-1][j], z[0][i][j+1], z[0][i][j-1]);
 }
 else {
 z[1][i][j] = 0.0;
 }
 }
 }
 MPI_Init( &argc, &argv); // start MPI
 MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
 MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
 int inter;
 int minusOne;
 int minusTwo;

 int stripe = size/mpi_size;
 // display the range of each rank
 cerr << "rank[" << my_rank << "]'s range = " << (my_rank * stripe) << " ~ " << (my_rank * stripe) +
stripe - 1 << endl; // display range of each rank
 // iterations for different times from 2 to max_time
 for(int ti = 2; ti < max_time ; ti++ ){
 if( ti % 3 == 2 ){
 inter = 2;
 minusOne = 1;
 minusTwo = 0;
 }
 else if( ti % 3 == 0){
 inter = 0;
 minusOne = 2;
 minusTwo = 1;
 }
 else {
 inter = 1;
 minusOne = 0;
 minusTwo = 2;
 }

 //send and receive data only if necessary
 if(mpi_size > 1 && interval != 0)
 {
 if(my_rank==0 ){
 for(int rank =1 ;rank < mpi_size ; rank++)
 MPI_Recv(z[inter][rank * stripe], stripe * size ,MPI_DOUBLE, rank , 0 ,
MPI_COMM_WORLD,&mStatus);
 }
 else{
 MPI_Send(z[inter][my_rank * stripe], stripe * size ,MPI_DOUBLE, 0 , 0 , MPI_COMM_WORLD);
 }
 }
 // print received data only if necessary
 if(my_rank ==0 && interval != 0 ){
 if(ti % interval == 0 || max_time - ti == 1)
 {
 cout << ti << endl ;
 for(int i = 0; i < size; i++){
 for(int j = 0; j < size; j++){
 cout << z[inter][i][j]<< " " ;
 }
 cout << endl;
 }
 cout << endl;
 }
 }
 int start=0;
 double previousRow[size];
 double nextRow[size];
 int startS =0;

 //---------------------
// share borders with all the ranks synchronously and takes care of synchronization by
diving into even and odd stripes 
// rank = 0
 if(mpi_size > 1){
 if(my_rank % 2 == 0)
 {
 if(my_rank == 0)
 { //master is always even
 MPI_Send(z[minusOne][stripe-1], size , MPI_DOUBLE, my_rank + 1 , 0 ,
MPI_COMM_WORLD); // send 24 to 0 + 1 rank
 MPI_Recv(z[minusOne][stripe], size , MPI_DOUBLE, my_rank + 1, 0 , MPI_COMM_WORLD,
&mStatus); // recv 25 from 0 + 1 rank
 }

 else if(mpi_size - my_rank == 1)
{ // if its a last process and even
 MPI_Send(z[minusOne][(my_rank * stripe)], size , MPI_DOUBLE, my_rank - 1, 0
,MPI_COMM_WORLD);
 MPI_Recv(z[minusOne][(my_rank * stripe)-1], size , MPI_DOUBLE, my_rank -1, 0 ,
MPI_COMM_WORLD, &mStatus);
}
else
{ // intermediate even process 2 50 - 74
 MPI_Send(z[minusOne][my_rank * stripe], size , MPI_DOUBLE, my_rank - 1 , 0
,MPI_COMM_WORLD ); // send 50 to 2 - 1
 MPI_Send(z[minusOne][(my_rank * stripe) + stripe - 1], size , MPI_DOUBLE, my_rank + 1 , 0
,MPI_COMM_WORLD); // send 74 to 2 + 1
 MPI_Recv(z[minusOne][(my_rank * stripe) - 1], size , MPI_DOUBLE, my_rank - 1 , 0
,MPI_COMM_WORLD, &mStatus); // recv 49 from 2 -1
 MPI_Recv(z[minusOne][(my_rank * stripe) + stripe ], size , MPI_DOUBLE, my_rank + 1 , 0
,MPI_COMM_WORLD, &mStatus); // recv 75 from 2 + 1
}
 }
 else {
 if(mpi_size - my_rank == 1)
 { // if its the last process and odd
 MPI_Recv(z[minusOne][(my_rank * stripe)-1], size , MPI_DOUBLE, my_rank -1, 0 ,
MPI_COMM_WORLD, &mStatus); // recv 74 from 3 - 1
 MPI_Send(z[minusOne][(my_rank * stripe)], size , MPI_DOUBLE, my_rank - 1, 0
,MPI_COMM_WORLD); // send 75 to 3-1
 }
 else { // intermediate odd process 1 25 - 49
 MPI_Recv(z[minusOne][(my_rank * stripe) - 1], size , MPI_DOUBLE, my_rank - 1 , 0
,MPI_COMM_WORLD, &mStatus); // for 1 recv 24 from 1-1
 MPI_Recv(z[minusOne][(my_rank * stripe) + stripe], size , MPI_DOUBLE, my_rank + 1 , 0
,MPI_COMM_WORLD, &mStatus); // rec 50 from 1 + 1
 MPI_Send(z[minusOne][my_rank * stripe], size , MPI_DOUBLE, my_rank - 1 , 0
,MPI_COMM_WORLD ); // send 25 to 1-1
 MPI_Send(z[minusOne][(my_rank * stripe) + stripe - 1], size , MPI_DOUBLE, my_rank + 1 , 0
,MPI_COMM_WORLD); // send 49 to 1 + 1

 }
 }
}
 // pass data for shroedingers formula for all the ranks
 #pragma omp parallel for
 for(int i = my_rank * stripe; i < (my_rank * stripe) + stripe ; i++ ) {

 for( int j = 0; j < size; j++ ) {
 if( (i > 0) && (i < (size-1)) && (j > 0) && (j < (size-1)) )
 {
 z[inter][i][j] = ZtTwo(z[minusTwo][i][j],z[minusOne][i][j], z[minusOne][i+1][j], z[minusOne][i-
1][j], z[minusOne][i][j+1], z[minusOne][i][j-1]);
 }
 else {
 z[inter][i][j] = 0.0;
 }
 }
 }
}


 // finish the timer
 if(my_rank == 0)
 cerr << "Elapsed time = " << time.lap( ) << endl;

 MPI_Finalize();
return 0;
}
