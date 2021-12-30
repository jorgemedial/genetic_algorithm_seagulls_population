#include "lib/random.c"
#include <stdio.h>

int main(){
    for(int i=0; i<100; i++) printf("%lu \n", ULONGran(13U));
    return 0;
}
