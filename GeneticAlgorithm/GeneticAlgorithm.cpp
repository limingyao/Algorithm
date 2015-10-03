#include<iostream>
#include<stdlib.h> 
#include<math.h> 
#include<time.h> 
#include<fstream> 
#include<stdio.h> 
#include<limits.h> 

using namespace std;

#define nofcity 80 //�������� 80 ��
#define maxgen 300 // maximum generation number 
#define popsize 50 //��ʼ��Ⱥ�� 50 
#define pm 0.05    //������
#define pc 0.85    //������
#define maxdot (nofcity+1)*(nofcity+1) 
#define nofcar 13  //total number of car here  
#define length 94  //Ⱦɫ�峤�� #define ee 30//��������ǰ 30EE ���� 
#define ll 30      //������ӳ��� 30 ���� 
#define v 60       //�ٶ� 60Km/h 
#define M1 5       //������ʱ�䴰�ͷ�����ϵ�� M1=5 
#define M2 1.1     //��ʱ�䴰�ͷ�����ϵ�� M2=1.1  

int  pop[(popsize+1)*length];//=(int*)malloc((popsize+1)*(length)*sizeof(int));
//����Ⱦɫ�峤�� nofcity+nofcar+1�����һ��Ϊ������Ⱦɫ��λ��     
int  tempop[(popsize+1)*length];//=(int*)malloc((popsize+1)*(length)*sizeof(int));
// ����Ⱦɫ�峤�� nofcity+nofcar+1�����һ��Ϊ������Ⱦɫ��λ��     
static int newpopsize;       
int temp[(popsize+1)*nofcity];//;temp=(int*)malloc((popsize+2)*nofcity*sizeof(int));
float bestval=32766;     
float val[popsize+1];       
int bestval_now=32766;        
int popx1[length];//;popx1=(int *)malloc((length)*sizeof(int));        
int popx2[length];//;popx2=(int *)malloc((length)*sizeof(int));         
double popfit[popsize+1];//;popfit=(double *)malloc((popsize+1)*sizeof(double));        
double p[popsize+1];//p=(double *)malloc((popsize+1)*sizeof(double));        
double c[popsize+1];//c=(double *)malloc((popsize+1)*sizeof(double));
//�ۼӸ��ʱ�����������ѡ��        
int aa[(popsize+1)*nofcity];//aa=(int*)malloc((popsize+1)*nofcity*sizeof(int));    
int can[nofcar];//capcity of each car     
int dotneed[nofcity];//need of each dot    
double pp;

static int randnum(int b) // Uniform Distribution ����  0-b+1 ������ 
{    
    double y;   
    if(b<0) {   
        printf("\nThe first parameter should be less than the second!");  
        exit(1);    
    }    
    y = (double)rand()/(RAND_MAX);  
    return (int)(b+1)*y;  
} 
static int randnum1(int b) // Uniform Distribution ����  1~b ������ 
{    
    double y;    
    if(b<0) {  
        printf("\nThe first parameter should be less than the second!");   
        exit(1);    
    }    
    y = (double)rand()/(RAND_MAX);    
    return (int)(b-1)*y+1;  
} 
struct jianju 
{  
    int number;   
    int x;   
    int y;   
    double need; 
}; 
struct jianju dot[maxdot];//��������㼯�� 
double julix(int &a,int &b)//jianju dot[maxdot],//��������֮��ľ��� 
{    
    double len;   
    len=(double)sqrt((double)(dot[a].x-dot[b].x)*(dot[a].x-dot[b].x)+(double)(dot[a].y-dot[b].y)*(dot[a].y-dot[b].y));   
    return len; 
}  

int check(int x[])//����Ⱦɫ���Ƿ�����ȶ���Լ������??????????? 
{   
    int n,sum;      
    int i;   
    sum=0;   n=0;   
    for(i=0;i<length;i++)sum+=x[i];//cout<<endl;      

    if(sum!=3240)//Ⱦɫ����������ܺ�Ϊ 3240���൱����������㶼����һ�� 0-80 ����    
        return -1;   
    sum=0;   
    for(i=0;i<length;i++)    
        sum+=x[i];   
    if(sum<0)    
        return -1;      
    else    
        sum=0;       
    for(i=0;i<nofcar;i++)     
    {      
        n++;      
        while(x[n]!=0)      
        {       
            sum+=dot[x[n]].need;             
            n++;      
        }             

        if(sum>can[i]) return -1;      
        else    
            sum=0;         
    }     
    return 1; 
} 
static  double  objval(int  x[])
//����Ŀ�꺯��ֵ��С��������������;�������� KC,·�̷���֮�͡��Լ�ʱ��ͷ�֮��,����ȡ������ 
{    
    int i;    
    int n,flag;   
    double sum;   
    int j;   
    double sump,sumlate;   
    double value;   
    value=0;   
    flag=1;   
    n=0;     
    sum=0;     
    sump=600;//����ʮ�����������ɷ���Ϊ���ٷ�����Ϊ��ʼ��ʱ��     
    sumlate=0;   //���� 6 ����������ɷ���Ϊ 300 ���ӿ�ʼ   //���������ĳ�ʼ����ʱ�̶�Ϊ 10��00��###########################      
    for(i=1;i<length;i++)//   *******length-1*****   
    {     
        flag=1;    
        sum+=julix(x[i-1],x[i]);      
        if(x[i]==0 )    
        {     
            sump=600;     
            flag=0;     
            sumlate=sumlate+0;    
        }    
        else if (x[i-1]==0)    
        {     
            sump=sump+julix(x[i-1],x[i])/v*60;      
            if(sump>=660 && sump<=720)//�ڹ涨ʱ�����ʹ������������ͷ�             
                sumlate=sumlate+0;             
            if(sump>720 && sump<=750)//���ڵ������Ը��� 1000 �ĳͷ�             
                sumlate=sumlate+M2*(sump-720);            
            if(sump>=630 && sump<660)             
                sumlate=sumlate+M2*(660-sump);            
            if(sump>750)             
                sumlate=sumlate+M1*(sump-750);            
            if(sump<630)             
                sumlate=sumlate+M1*(630-sump);    
        }    
        else    
        {     
            sump=sump+julix(x[i-1],x[i])/v*60;     
            if(sump>=660 && sump<=720)//�ڹ涨ʱ�����ʹ������������ͷ�             
                sumlate=sumlate+0;             
            if(sump>720 && sump<=750)//���ڵ������Ը��� 1000 �ĳͷ�             
                sumlate=sumlate+M2*(sump-720);            
            if(sump>=630 && sump<660)             
                sumlate=sumlate+M2*(660-sump);            
            if(sump>750)             
                sumlate=sumlate+M1*(sump-750); 
            if(sump<630)             
                sumlate=sumlate+M1*(630-sump);    
        }   
    }        
    value=sum+sumlate*0.5; 
    return value; 
} 
static void initialize()// ##################��ʼ��############### 
{      
    int i,j,n,k,l,m,point,sum;      
    int flag,gen;       
    int a[popsize*nofcity];   
    srand((unsigned)time(NULL));      
    rand();   
    cout<<"����Ϊ��Ⱥ��ʼ�����飬�����ĵȺ�~~~~~~��лл��"<<endl;   
    for(j=0;j<nofcity;j++)   
    {    
        do    
        {         
            a[j]=randnum1(nofcity);     
            flag=0;     
            for(n=0;n<j;n++)     
            {         
                if(a[n]==a[j])flag=1;         
            }    
        }while(flag!=0);   
    }//�����ĳ�ʼ����������      
    for(n=0;n<=popsize;n++)   
    {       
        for(i=0;i<nofcity;i++)    
        {       
            aa[n*nofcity+i]=a[i];    
        }   
    }      
    for(n=0;n<popsize;n++)   
    {        
        k=nofcity-1;           
        for(i=0;i<nofcity;i++)     
        {            
            point=randnum(k);             
            if(aa[n*nofcity+point]!=0)     
            {               
                temp[n*nofcity+i]=aa[n*nofcity+point];               
                aa[n*nofcity+point]=0;     
            }            
            m=0;          
            for(j=0;j<nofcity;j++)     
            {           
                if(aa[n*nofcity+j]!=0){aa[n*nofcity+m]=aa[n*nofcity+j];m++;}     
            }            
            k--;     
        }   
    }        
    cout<<endl;        
    for(n=0;n<popsize;n++)     
    {           
        cout<<"��"<<n+1<<"����ʼ����Ϊ��"<<endl;      
        sum=0;         
        for(i=0;i<nofcity;i++)      
        {          
            sum+=temp[n*nofcity+i];      
            // cout<<temp[n*nofcity+i]<<"   ";      
        }            
        // cout<<"TOTAL IS : "<<sum<<endl;     
    }        
    //#####################����Ⱦɫ�壬��ʽ�� 01302450���������� 0 �����ʼ���������??????have promble!!!        
    for(i=0;i<popsize;i++)     
    {        
        m=0;k=0;         
        tempop[i*length+k]=0;         
        for(j=0;j<nofcar-1;j++)    
        {        
            sum=0;        
            flag=0;        
            while(flag!=1)     
            {             
                sum+=dot[temp[i*nofcity+m]].need;              
                if(sum<=can[j])       
                {              
                    k++;              
                    tempop[i*length+k]=temp[i*80+m];             
                    flag=0;              
                    m++;       
                    //  cout<<tempop[i*length+k];       
                }             
                else flag=1;      
            }          
            k++;        
            tempop[i*length+k]=0;    
        }     
        if(j=nofcar-1)     
        {      
            do      
            {       
                k++;       
                tempop[i*length+k]=temp[i*nofcity+m];       
                m++;      
            }while(m<80);             
            tempop[i*length+length-1]=0;     
        }     
    }        
    n=0;        
    cout<<"###########################"<<endl;        
    for(i=0;i<popsize;i++)     
    {           
        for(j=0;j<length;j++)     
        {      
            popx1[j]=tempop[i*length+j];     
            //cout<<popx1[j]<<" ";     
        }//new            
        if(check(popx1)!=1){cout<<"��  "<<i+1<<"����ʼ���򲻿��У�"<<endl;}            
        else        
        {              
            for(l=0;l<length;l++)        
            {               
                pop[n*length+l]=tempop[i*length+l];        
            }           
            n++;      
        }     
    }         
    newpopsize=n;         
    cout<<"###########################"<<endl;         
    cout<<"     ��ʼ���н����:   "<<n<<endl;         
    cout<<"###########################"<<endl;         
    for(i=0;i<newpopsize;i++)      
    {           
        cout<<"��  "<<i+1<<"����ʼ��Ⱦɫ��Ϊ��  ";     
        for(j=0;j<length;j++){popx1[j]=pop[i*length+j];}cout<<"   best_pop: "<<objval(popx1)<<"   ";           
        for(j=0;j<length;j++)     
        {             
            //   cout<<pop[i*length+j]<<" ";     
        }              
        cout<<endl;      
    } 
}   

static double leng(int x[])//����Ⱦɫ���·�̳��� 
{    
    int i,sum=0;       
    for(i=1;i<length;i++)//   *******length-1*****   
    {     
        sum+=julix(x[i-1],x[i]);   
    } 
    return sum; 
} 
static void keepbest() 
{   
    int i,j,n,l;   
    n=0;   
    for(i=0;i<popsize;i++)   
    {       for(j=0;j<length;j++){popx1[j]=pop[i*length+j];}//new          
        if(check(popx1)==1)       
        {          
            if( objval(popx1)<abs(bestval))       
            {                   
                bestval=objval(popx1);                
                for(j=0;j<length;j++){pop[popsize*length+j]=pop[i*length+j];}// ����ѵ�Ⱦɫ����򱣴��� popsize ������       
            }           
            for(l=0;l<length;l++)       
            {            
                tempop[n*length+l]=pop[i*length+l];       
            }           
            n++;       
        }   
    }      
    newpopsize=n; 
}  

//#################################################ѡ�� popsize ����Ⱥ 
static void selection() 
{   
    int i,j,n;   
    //#############################����Ⱦɫ���Ŀ�꺯��ֵ��С      
    for(i=0;i<newpopsize;i++)   
    {        
        for(j=0;j<length;j++){popx1[j]=pop[i*length+j];}           
        val[i]= objval(popx1);//new           
        cout<<"   Ŀ�꺯��ֵΪ��"<<i+1<<"     "<<val[i]<<endl;   
    } 
    //b########################################���浱ǰ��ѵ�Ⱦɫ��      
    for(i=0;i<newpopsize;i++)   
    {        
        if(i==0)       
            bestval_now=val[0];//objval(popx1,last,unloadt);           
        else     
        {           
            if( val[i]<abs(bestval_now))     
            {            
                bestval_now=val[i];            
                if(bestval_now<abs(bestval))      
                {             
                    bestval=bestval_now;                
                    for(j=0;j<length;j++){pop[popsize*length+j]=pop[i*length+j];}// ����ѵ�Ⱦɫ����򱣴��� popsize ������      
                }     
            }       
        }   
    } 
    //cout<<"  ��ǰ���ź���Ŀ�꺯��ֵ��    "<<bestval_now<<"      "<<endl;     
    double sumfit=0;     
    for(i=0;i<newpopsize;i++)     
    {         
        popfit[i]=bestval_now/val[i];//����Ⱦɫ����Ӧ��ֵ              
        sumfit+=popfit[i];//����Ⱥ����Ӧ��֮��     
    }     
    for(i=0;i<newpopsize;i++)          
        p[i]=popfit[i]/sumfit;//�� i ����Ⱥ��ѡ��ĸ���       
    c[0]=p[0];   
    for(i=1;i<newpopsize;i++)   
    {    
        c[i]=c[i-1]+p[i];      
        //  cout<<c[i]<<"   "<<endl;   
    }//c Ϊ�ۼӵ�ѡ����� //##########################################��Ⱥѡ��ѡ�� POPSIZE  ����Ⱥ�����̶ķ�ʽ     
    for(j=0;j<length;j++)//�����ŵ�Ⱦɫ�������һ����   
    {       
        pop[j]=pop[popsize*length+j];   
    }     
    for(i=1;i<popsize;i++)     
    {      
        pp=(double)rand()/(RAND_MAX);//�������һ��������Ϊ���̶ı��   
        if(pp<c[0])    
        {         
            for (j=0;j<length;j++)         
                tempop[i*length+j]=pop[j];    
        }       
        else    
        {        
            for(j=0;j<newpopsize;j++)          
                if(pp>=c[j] && pp<c[j+1])      
                {           
                    for (n=0;n<length;n++)                   
                        tempop[i*length+n]=pop[j*length+n];          
                    break;      
                }    
        }     
    } 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~����ѡ����µ���Ⱥ��� 
    //������������������������������������������������������������������Ⱥ����      
    for(i=1;i<popsize;i++)   
    {       
        for(j=0;j<length;j++)    
        {         
            pop[i*length+j]=tempop[i*length+j];       
            // cout<<pop[i*length+j]<<" ";    
        }    
        //cout<<"No:   "<<i<<"   "<<endl;   
    }  
    //����������������������������������������Ⱥ��� 
}  


static void crossover() 
{ 
    int i,j,men,one,point0,point1,point2,mm,flag1,flag2,flag3,flag4,flag5; 
    int f,t; 
    double x; 
    int first=0;//��¼��ѡ���Ⱦɫ����� 
    for(men=0;men<popsize;men++) 
    {   
        x=(double)rand()/(RAND_MAX);//�������һ�������뽻���� PC ���бȽ�    
        if(x<pc)    
        {      
            first++;      
            if (first%2==0)//������Ⱦɫ�屻ѡ����Ϊ������������һ�ֵĽ�������      
            {       
                flag1=0;flag2=0;       //����Ȳ��� 0-nofcar ֮���һ����������ѡ��һ����·���� 03280����Ϊƥ��ν���       
                point0=randnum(nofcar-1);//���� 0-nofcar ֮�������   �� men ���Ⱦɫ����Ե�  men Ϊ�� 1 one Ϊĸ 1      
                for(i=0;i<length;i++)        
                    if(pop[men*length+i]==0)        
                    {         
                        if(flag1==point0){point1=i;break;}//��¼·����Ⱦɫ�� men �����λ�� point1         
                        else          
                            flag1++; 
                    }       
                for(i=0;i<length;i++)        
                    if(pop[men*length+i]==0)        
                    {         
                        if(flag2==(point0+1)){point2=i;break;}//��¼·����Ⱦɫ�� men ���յ�λ�� point2         
                        else          
                            flag2++;        
                    }          
                for(i=0;i<length;i++)//���游��Ⱦɫ��       
                {        
                    tempop[one*length+i]=pop[one*length+i];        
                    tempop[men*length+i]=pop[men*length+i];       
                }       
                mm=0;        
                for(i=0;i<length;i++)      
                {       
                    if(i<=point1 || i>=point2)//&& i<point2       
                    {        if(pop[men*length+i]!=0)        
                        {         
                            for(j=mm;j<length;j++)         
                            {             
                                flag4=0;          
                                if(pop[one*length+j]!=0)          
                                {             
                                    flag3=1;              
                                    for(f=point1+1;f<point2;f++)               
                                        if(pop[one*length+j]==pop[men*length+f])            
                                        {flag3=0;break;}            
                                    mm=j;            
                                    if(flag3!=0){pop[men*length+i]=pop[one*length+j];flag4=1;}          
                                }          
                                if(flag4==1){mm++;break;}         
                            }        
                        }       
                    }      
                }     
                mm=0;        
                for(i=0;i<length;i++)      
                {       
                    if(i<=point1 || i>=point2)//&& i<point2       
                    {        
                        if(pop[one*length+i]!=0)        
                        {         
                            for(j=mm;j<length;j++)         
                            {             
                                flag4=0;          
                                if(pop[men*length+j]!=0)          
                                {             
                                    flag3=1;              
                                    for(f=point1+1;f<point2;f++)               
                                        if(pop[men*length+j]==pop[one*length+f])            
                                        {flag3=0;break;}            
                                    mm=j; 
                                    if(flag3!=0){pop[men*length+i]=pop[one*length+j];flag4=1;}
                                }          
                                if(flag4==1){mm++;break;}         
                            }        
                        }       
                    }      
                }      
                mm=0;        
                for(i=0;i<length;i++)      
                {       
                    if(i<=point1 || i>=point2)//&& i<point2       
                    {        
                        if(pop[one*length+i]!=0)        
                        {         
                            for(j=mm;j<length;j++)         
                            {             
                                flag4=0;          
                                if(pop[men*length+j]!=0)          
                                {             
                                    flag3=1;              
                                    for(f=point1+1;f<point2;f++)               
                                        if(pop[men*length+j]==pop[one*length+f])            
                                        {flag3=0;break;}            
                                    mm=j;            
                                    if(flag3!=0){pop[one*length+i]=pop[men*length+j];flag4=1;}          
                                }          
                                if(flag4==1){mm++;break;}         
                            }        
                        }       
                    }      
                }      
            }        

            else       
                one=men;    
        } 
    }//##############################�������㽻����� 
} 

static void mutation(void) 
{   
    int i,t;   
    int point1,point2,flag1;   
    double y;   
    for(i=0;i<popsize;i++)//���ѡ��������Ž��н��������������------------2 ��������   
    {        
        y=(double)rand()/(RAND_MAX);//�������һ������������� pm ���бȽ�       
        if(y<pm)    
        {        
            do     
            {         
                flag1=0;         
                point1=randnum1(length-3);//����һ�� 1-15 �����������콻���ĵ�һ��         
                point2=randnum1(length-3);//����һ�� 1-15 �����������콻���ĵڶ���          
                if(pop[i*length+point1]!=0 && pop[i*length+point2]!=0 && point1!=point2)      
                    flag1=1;     
            }while(flag1!=1);        
            if(flag1==1)     
            {          
                t=pop[i*length+point1];         
                pop[i*length+point1]=pop[i*length+point2];         
                pop[i*length+point2]=t;     
            }    
        }   
    } //##############################�����������  
}    

void report(void)   
{    
    int i;    
    cout<<"game over"<<endl;      
    cout<<"the best value is "<<bestval<<endl;    
    cout<<"the num is :";       
    for(i=0;i<length;i++)    
    {popx1[i]=pop[popsize*length+i];}    
    cout<<"  ����·��̵�·��Ϊ��  "<<leng(popx1)<<endl;    
    for(i=0;i<length;i++)cout<<pop[popsize*length+i]<<" ";    
    cout<<endl;    
    for(i=0;i<length-1;i++)      
    {     
        cout<<julix(pop[popsize*length+i],pop[popsize*length+i+1])<<"   ";       
    }    
}     

int main () 
{   
    int ii,gen,ss=0,index=0;   
    double maxgener;     
    int max1=0,max2=0;   
    for(ii=1;ii<nofcar;ii++)//capcity of A car =500   
    {    
        can[1]=500;    
        can[2]=500;    
        can[3]=500;    
        can[4]=500;    
        can[5]=500;    
        can[6]=500;    
        can[7]=500;    
        can[8]=500;    
        can[9]=500;    
        can[10]=500;    
        can[11]=400;    
        can[12]=400;     
    }   
    can[0]=500;    
    FILE *infile;//��ȡ�������Ϣ   
    if((infile=fopen("haodata.txt","r"))==NULL)   
    {    
        fprintf(infile,"\nCan not open input file!\n");    
        exit(1);   
    }   
    for(ii=0;ii<=nofcity;ii++)   
    {    
        fscanf(infile,"%d",&dot[ii].number);    
        fscanf(infile,"%d",&dot[ii].x);    
        fscanf(infile,"%d",&dot[ii].y);    
        fscanf(infile,"%lf",&dot[ii].need);   
    }   
    initialize();   
    keepbest(); 
    //do{ 
    for(maxgener=0;maxgener<100;maxgener++) 
    {    
        for(gen=0;gen<maxgen;gen++)    
        {     
            selection();   
            crossover();    
            mutation();    
            keepbest();    
        }      
        cout<<"��  "<<gen<<"�����Ϊ    "<<bestval<<endl;     
        max1=bestval;     
        initialize();   
        keepbest();    
        for(gen=0;gen<maxgen;gen++)    
        {     
            selection();    
            crossover();    
            mutation();    
            keepbest();    
        }   
        cout<<"���֣���  "<<gen<<"�����Ϊ    "<<bestval<<endl;   
        max2=bestval;   
        if(abs(max1-max2)<=10 && max1<=6900) 
            break; 
    }      
    report();   
    cout<<maxgener<<endl;  
    return 0;
} 
