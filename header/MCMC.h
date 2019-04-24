#define MAXDATASIZE 500
#define MAXMCTIMES 20000
#define INF 1e20
#define LOWBODUND -1e200
#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))

//global par
double lambda=0.001;  //mutation rate
double Nbegin=1e5;
int seqlength=-1;
int flagcount=0;
#include "fasta.h"
#include "tree.h"
#include "randTree.h"
#include "../init/init_p3n.h"
#include "Coleselikelihood.h"
#include "Treelikelihood.h"
#include "paraTree.h"

int checkTree(){
    int i,k;
    int j=randindex;
    for(i=0;i<j;i++){
        if(randlength[i]==randlength[j]){
            int kk=0;
            for(k=0;k<randlength[j];k++)
                if(randstore[i][k]==randstore[j][k])
                    kk++;
            if(kk==randlength[j]) {
                randlength[j]=0;
                return 1;
            }
        }
    }
    return 0;
}
void initcheck(){
    int i;
    for(i=0;i<randindex;i++) randlength[i]=0;
    randindex=0;
}

void MCMC(){
    
    int ncount=0;
    int tcount=0;
    double ALPHA,r;
    double llhood;
    double pllhood=LOWBODUND;
    BiNode *TREE,*TEMPTREE;
    int idnum=getpid();
    char fname[40];
    FILE *fp;
    double tlambda,Ntmp;
    int wrank,wsize;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    sprintf(fname,"result/p3n_%d.txt",idnum);
    
    srand((unsigned)time(NULL)+idnum*wrank+strlen(fname));
    
    clock_t startime,endtime;
    double cpu_time_used;
    startime=clock();
    if (wrank==0) {
        TREE=init();
        fp=fopen(fname,"w+");
    }
    MPI_Bcast(&seqlength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    while (ncount<MAXMCTIMES){
        
        if(wrank==0){
            tlambda=lambda;
            r=rand()/(RAND_MAX+0.0)-0.5;
            lambda*=exp(r);
            
            Ntmp=Nbegin;
            r=rand()/(RAND_MAX+0.0)-0.5;
            Nbegin*=exp(0.2*r);
            
            //boast tlambda N Treeroot
            
        }
        MPI_Barrier(MPI_COMM_WORLD);  //barrier
        MPI_Bcast(&lambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Nbegin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(wrank==0){
            int tmpx;
            BiNode *TempsendTree;
            TempsendTree=NULL;
            initcheck();
            for(tmpx=1;tmpx<wsize;tmpx++) {
                // random change tree
                int checkcount=0;
                do{
                    if(TempsendTree!=NULL) DestroyTree(TempsendTree);
                    TempsendTree=cloneTree(TREE,NULL);
                    TempsendTree=NNItheTree(TempsendTree);
                    checkcount++;
                }while(checkTree()&&checkcount<10);
                randindex++;
                r=rand()/(RAND_MAX+0.0);
                RandTime(TempsendTree,TempsendTree->NodeTime+log(1-r)*10);
                ParaSendTree(TempsendTree,tmpx);
                DestroyTree(TempsendTree);
                TempsendTree=NULL;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);  //barrier

        if(wrank!=0){
            TEMPTREE=malloc(sizeof(BiNode));
            ParaRecvTree(TEMPTREE,0,NULL);
            //get out the time value
            flagcount=0;
            TimeListN=0;
            PrintTime(TEMPTREE);
            sortTimeN=sortTime();
            
            //calculate likelihood
            llhood=LogTreeLikelihood(TEMPTREE);
            //printf("\n%lf ",llhood,wrank);
            llhood+=coallikelihood(sortedTime,sortTimeN);
            //printf("%lf %d\n",llhood,wrank);
            MPI_Send(&llhood, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        int tmpi=-1;
        double tmax=LOWBODUND;
        MPI_Barrier(MPI_COMM_WORLD);
        if (wrank==0){
            int tmpx;
            for(tmpx=1;tmpx<wsize;tmpx++){
                MPI_Recv(&llhood, 1, MPI_DOUBLE, tmpx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (tmax<llhood){
                    tmax=llhood;
                    tmpi=tmpx;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&tmpi, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (wrank==tmpi) ParaSendTree(TEMPTREE,0);
        if (wrank==0){
            if(tmpi!=-1){
                TEMPTREE=malloc(sizeof(BiNode));
                TEMPTREE->parent=NULL;
                ParaRecvTree(TEMPTREE,tmpi,NULL);
            }
            llhood=tmax;
            ALPHA=MIN(0,llhood-pllhood);
            r=log(rand()/(RAND_MAX+0.0));
            if (r<ALPHA && tmpi!=-1){
                pllhood=llhood;
                DestroyTree(TREE);
                TREE=TEMPTREE;
                tcount++;
                endtime=clock();
                cpu_time_used = ((double) (endtime - startime)) / CLOCKS_PER_SEC;
                printf("%lf n=%d time=%lf\n",llhood,tcount,cpu_time_used);
                fprintf(fp,"%e %e\n",lambda,Nbegin);
            }
            else{
                DestroyTree(TEMPTREE);
                lambda=tlambda;
                Nbegin=Ntmp;
            }
        }
        if(wrank!=0) DestroyTree(TEMPTREE);
        ncount++;
    }
    if (wrank==0) {
        fclose(fp);
        printTree(TREE);
    }
    
}

