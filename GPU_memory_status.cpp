#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <unistd.h>

int main(){
  long unsigned int free_byte ;

  long unsigned int  total_byte ;

  while(1) {
    int cuda_status = cudaMemGetInfo( &free_byte, &total_byte  ) ;

    if ( cudaSuccess != cuda_status  ){
      printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
      return 1;
    }



    double free_db = (double)free_byte ;

    double total_db = (double)total_byte ;

    double used_db = total_db - free_db ;

    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
        used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
    usleep(10000);
  }

}
