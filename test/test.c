#include <stdio.h>
#include <stddef.h>

typedef struct PACK{
    double l;
    double N;
    struct PACK *T;
}PACK;

int main(){
    
    printf("%d\n",offsetof(PACK, N));
    
    return 0;
}
