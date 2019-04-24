
int ParaSendTree(BiNode *T,int wrank){

    MPI_Send(T,sizeof(BiNode),MPI_BYTE,wrank,0,MPI_COMM_WORLD);
    MPI_Send(T->name,10,MPI_BYTE,wrank,0,MPI_COMM_WORLD);
    MPI_Send(&T->NodeTime,1,MPI_DOUBLE,wrank,0,MPI_COMM_WORLD);
    if (T->sequence!=NULL){
        MPI_Send(T->sequence,500,MPI_BYTE,wrank,0,MPI_COMM_WORLD);
        return 1;
    }
    else{
        MPI_Send(T->rchi,sizeof(BiNode),MPI_BYTE,wrank,0,MPI_COMM_WORLD);
        MPI_Send(T->lchi,sizeof(BiNode),MPI_BYTE,wrank,0,MPI_COMM_WORLD);
        ParaSendTree(T->rchi,wrank);
        ParaSendTree(T->lchi,wrank);
    }
    return 1;
}
int ParaRecvTree(BiNode *T,int wrank,BiNode *Tp){
    
    T->parent=Tp;
    MPI_Recv(T,sizeof(BiNode), MPI_BYTE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    T->name=calloc(10,sizeof(char));
    MPI_Recv(T->name,10,MPI_BYTE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&T->NodeTime,1,MPI_DOUBLE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(T->sequence!=NULL){
        T->sequence=calloc(500,sizeof(char));
        MPI_Recv(T->sequence,500,MPI_BYTE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        T->rchi=NULL;
        T->lchi=NULL;
        return 1;
     }
     else{
         T->rchi=malloc(sizeof(BiNode));
         MPI_Recv(T->rchi,sizeof(BiNode), MPI_BYTE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         T->lchi=malloc(sizeof(BiNode));
         MPI_Recv(T->lchi,sizeof(BiNode), MPI_BYTE, wrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         T->sequence=NULL;
         ParaRecvTree(T->rchi,wrank,T);
         ParaRecvTree(T->lchi,wrank,T);
     }
    return 1;
}
