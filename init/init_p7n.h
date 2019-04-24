//
#define DATADAYS 10
double datatime[DATADAYS+1]={0,2306,2369,2726,3028,3495,3992,4503,5092,5450,5869};
int datanum[DATADAYS+1]={0};

BiNode *init()
{

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int   L;
    char *filename="fasta/haplotypes_p7_p17.fasta";  //

	int n=0;
	BiNode* S[MAXDATASIZE];
	int top=-1;
	ffp = OpenFASTA(filename);

    if (ffp==NULL) {

        printf("wrong input!\n");
        return NULL;
    }

	while(ReadFASTA(ffp, &seq, &name, &L)){
        S[n]=calloc(1,sizeof(BiNode));
		S[n]->lchi=NULL;
		S[n]->rchi=NULL;
		S[n]->parent=NULL;
        S[n]->sequence=calloc(500,sizeof(char));
        S[n]->name=calloc(10,sizeof(char));
	    strcpy(S[n]->sequence,seq);
	    strncpy(S[n]->name,name,9);
        //
	    if (strcmp(S[n]->name,"days_2306")==0) {S[n]->NodeTime=datatime[1];datanum[1]++;}
	    else if (strcmp(S[n]->name,"days_2369")==0) {S[n]->NodeTime=datatime[2];datanum[2]++;}
	    else if (strcmp(S[n]->name,"days_2726")==0) {S[n]->NodeTime=datatime[3];datanum[3]++;}
	    else if (strcmp(S[n]->name,"days_3028")==0) {S[n]->NodeTime=datatime[4];datanum[4]++;}
        else if (strcmp(S[n]->name,"days_3495")==0) {S[n]->NodeTime=datatime[5];datanum[5]++;}
	    else if (strcmp(S[n]->name,"days_3992")==0) {S[n]->NodeTime=datatime[6];datanum[6]++;}
	    else if (strcmp(S[n]->name,"days_4503")==0) {S[n]->NodeTime=datatime[7];datanum[7]++;}
	    else if (strcmp(S[n]->name,"days_5092")==0) {S[n]->NodeTime=datatime[8];datanum[8]++;}
	    else if (strcmp(S[n]->name,"days_5450")==0) {S[n]->NodeTime=datatime[9];datanum[9]++;}
	    else if (strcmp(S[n]->name,"days_5869")==0) {S[n]->NodeTime=datatime[10];datanum[10]++;}
	    //else if (strcmp(S[n]->name,"days_2922")==0) {S[n]->NodeTime=datatime[11];datanum[11]++;}
        //else if (strcmp(S[n]->name,"days_2996")==0) {S[n]->NodeTime=datatime[12];datanum[12]++;}
        else {
            printf("no name %s\n",S[n]->name);
            //exit(1);
        }

       //printf("length %d\n",strlen(S[n]->sequence));
        if (strlen(S[n]->sequence)!=402)
            printf("%s\n%d\n",name,strlen(S[n]->sequence));
	    free(name);
	    free(seq);
		n++;
	}
	seqlength=strlen(S[1]->sequence);
	//printf("number of node=%d\n",n);
	top=n-1;
	//printf("length %d\n",seqlength);
	int i,j,len=n;
	int judge[MAXDATASIZE];
	//for (i=0;i<len;i++) printf("S[%d]=%d\n",i,S[i]->data);

	for (i=0;i<MAXDATASIZE;i++) judge[i]=i;
	//srand(time(0));
	int r;
	while(len>1){ //shuffle algorithm to build the Tree
		top++;
		S[top]=calloc(1,sizeof(BiNode));
		S[top]->name=calloc(10,sizeof(char));
		strcpy(S[top]->name,"AYNULL");
		S[top]->sequence=NULL;
		r=rand()%len;
		j=judge[r];
		S[top]->lchi=S[j];
        S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		//printf("j=%d\n",j);

		r=rand()%len;
		j=judge[r];
		S[top]->rchi=S[j];
		S[j]->parent=S[top];
		len--;
		judge[r]=judge[len];
		judge[len]=top;
		len++;

		S[top]->NodeTime=MIN(S[top]->lchi->NodeTime,S[top]->rchi->NodeTime)-20;
		S[top]->parent=NULL;
		//printf("j=%d\n\n",j);

	}
    //printf("top=%d\n",top);
	//for (i=0;i<top+1;i++) printf("S[%d]=%d\n",i,S[i]->data);
	return S[top];
}
