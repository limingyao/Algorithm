#include<iostream>
#include<stdlib.h> 
#include<math.h> 
#include<time.h> 
#include<fstream> 
#include<stdio.h> 
#include<limits.h> 

using namespace std;

#define nofcity 80 //需求点个数 80 个
#define maxgen 300 // maximum generation number 
#define popsize 50 //初始种群数 50 
#define pm 0.05    //变异率
#define pc 0.85    //交叉率
#define maxdot (nofcity+1)*(nofcity+1) 
#define nofcar 13  //total number of car here  
#define length 94  //染色体长度 #define ee 30//最早能提前 30EE 分钟 
#define ll 30      //最迟能延长到 30 分钟 
#define v 60       //速度 60Km/h 
#define M1 5       //超过软时间窗惩罚因子系数 M1=5 
#define M2 1.1     //软时间窗惩罚因子系数 M2=1.1  

int  pop[(popsize+1)*length];//=(int*)malloc((popsize+1)*(length)*sizeof(int));
//构造染色体长度 nofcity+nofcar+1，最后一个为存放最佳染色体位置     
int  tempop[(popsize+1)*length];//=(int*)malloc((popsize+1)*(length)*sizeof(int));
// 构造染色体长度 nofcity+nofcar+1，最后一个为存放最佳染色体位置     
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
//累加概率便于最后的轮盘选择。        
int aa[(popsize+1)*nofcity];//aa=(int*)malloc((popsize+1)*nofcity*sizeof(int));    
int can[nofcar];//capcity of each car     
int dotneed[nofcity];//need of each dot    
double pp;

static int randnum(int b) // Uniform Distribution 返回  0-b+1 的整数 
{    
    double y;   
    if(b<0) {   
        printf("\nThe first parameter should be less than the second!");  
        exit(1);    
    }    
    y = (double)rand()/(RAND_MAX);  
    return (int)(b+1)*y;  
} 
static int randnum1(int b) // Uniform Distribution 返回  1~b 的整数 
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
struct jianju dot[maxdot];//建立坐标点集合 
double julix(int &a,int &b)//jianju dot[maxdot],//返回两点之间的距离 
{    
    double len;   
    len=(double)sqrt((double)(dot[a].x-dot[b].x)*(dot[a].x-dot[b].x)+(double)(dot[a].y-dot[b].y)*(dot[a].y-dot[b].y));   
    return len; 
}  

int check(int x[])//测试染色体是否满足既定的约束条件??????????? 
{   
    int n,sum;      
    int i;   
    sum=0;   n=0;   
    for(i=0;i<length;i++)sum+=x[i];//cout<<endl;      

    if(sum!=3240)//染色体各个排序总和为 3240，相当于所有需求点都经历一遍 0-80 都有    
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
//计算目标函数值大小，包括三个部分;常数部分 KC,路程费用之和、以及时间惩罚之和,这里取后两者 
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
    sump=600;//早上十点出发，换算成分钟为六百分钟作为起始点时刻     
    sumlate=0;   //早上 6 点出发，换成分钟为 300 分钟开始   //各个车俩的初始出发时刻都为 10：00。###########################      
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
            if(sump>=660 && sump<=720)//在规定时间内送达，理想情况，不惩罚             
                sumlate=sumlate+0;             
            if(sump>720 && sump<=750)//低于底限所以给以 1000 的惩罚             
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
            if(sump>=660 && sump<=720)//在规定时间内送达，理想情况，不惩罚             
                sumlate=sumlate+0;             
            if(sump>720 && sump<=750)//低于底限所以给以 1000 的惩罚             
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
static void initialize()// ##################初始化############### 
{      
    int i,j,n,k,l,m,point,sum;      
    int flag,gen;       
    int a[popsize*nofcity];   
    srand((unsigned)time(NULL));      
    rand();   
    cout<<"正在为种群初始化数组，请耐心等候~~~~~~，谢谢！"<<endl;   
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
    }//需求点的初始随机排列组合      
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
        cout<<"第"<<n+1<<"个初始排列为："<<endl;      
        sum=0;         
        for(i=0;i<nofcity;i++)      
        {          
            sum+=temp[n*nofcity+i];      
            // cout<<temp[n*nofcity+i]<<"   ";      
        }            
        // cout<<"TOTAL IS : "<<sum<<endl;     
    }        
    //#####################构造染色体，形式如 01302450，添加虚拟点 0 插入初始排列组合中??????have promble!!!        
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
        if(check(popx1)!=1){cout<<"第  "<<i+1<<"个初始排序不可行！"<<endl;}            
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
    cout<<"     初始可行解个数:   "<<n<<endl;         
    cout<<"###########################"<<endl;         
    for(i=0;i<newpopsize;i++)      
    {           
        cout<<"第  "<<i+1<<"个初始种染色体为：  ";     
        for(j=0;j<length;j++){popx1[j]=pop[i*length+j];}cout<<"   best_pop: "<<objval(popx1)<<"   ";           
        for(j=0;j<length;j++)     
        {             
            //   cout<<pop[i*length+j]<<" ";     
        }              
        cout<<endl;      
    } 
}   

static double leng(int x[])//计算染色体的路程长度 
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
                for(j=0;j<length;j++){pop[popsize*length+j]=pop[i*length+j];}// 将最佳的染色体基因保存在 popsize 数组中       
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

//#################################################选择 popsize 个种群 
static void selection() 
{   
    int i,j,n;   
    //#############################计算染色体的目标函数值大小      
    for(i=0;i<newpopsize;i++)   
    {        
        for(j=0;j<length;j++){popx1[j]=pop[i*length+j];}           
        val[i]= objval(popx1);//new           
        cout<<"   目标函数值为；"<<i+1<<"     "<<val[i]<<endl;   
    } 
    //b########################################保存当前最佳的染色体      
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
                    for(j=0;j<length;j++){pop[popsize*length+j]=pop[i*length+j];}// 将最佳的染色体基因保存在 popsize 数组中      
                }     
            }       
        }   
    } 
    //cout<<"  当前最优函数目标函数值：    "<<bestval_now<<"      "<<endl;     
    double sumfit=0;     
    for(i=0;i<newpopsize;i++)     
    {         
        popfit[i]=bestval_now/val[i];//计算染色体适应度值              
        sumfit+=popfit[i];//所有群体适应度之和     
    }     
    for(i=0;i<newpopsize;i++)          
        p[i]=popfit[i]/sumfit;//第 i 个种群被选择的概率       
    c[0]=p[0];   
    for(i=1;i<newpopsize;i++)   
    {    
        c[i]=c[i-1]+p[i];      
        //  cout<<c[i]<<"   "<<endl;   
    }//c 为累加的选择概率 //##########################################种群选择，选出 POPSIZE  个种群，轮盘赌方式     
    for(j=0;j<length;j++)//将最优的染色体加入这一带中   
    {       
        pop[j]=pop[popsize*length+j];   
    }     
    for(i=1;i<popsize;i++)     
    {      
        pp=(double)rand()/(RAND_MAX);//随机产生一个数，作为轮盘赌标杆   
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
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~重新选择出新的种群完成 
    //………………………………………………………………………………保存种群操作      
    for(i=1;i<popsize;i++)   
    {       
        for(j=0;j<length;j++)    
        {         
            pop[i*length+j]=tempop[i*length+j];       
            // cout<<pop[i*length+j]<<" ";    
        }    
        //cout<<"No:   "<<i<<"   "<<endl;   
    }  
    //……………………………………………保存种群完毕 
}  


static void crossover() 
{ 
    int i,j,men,one,point0,point1,point2,mm,flag1,flag2,flag3,flag4,flag5; 
    int f,t; 
    double x; 
    int first=0;//记录被选择的染色体个数 
    for(men=0;men<popsize;men++) 
    {   
        x=(double)rand()/(RAND_MAX);//随机产生一个数，与交叉率 PC 进行比较    
        if(x<pc)    
        {      
            first++;      
            if (first%2==0)//有两个染色体被选择作为父代，进行下一轮的交叉运算      
            {       
                flag1=0;flag2=0;       //随机先产生 0-nofcar 之间的一个整数（即选择一个子路径如 03280）作为匹配段交叉       
                point0=randnum(nofcar-1);//产生 0-nofcar 之间的整数   对 men 这个染色体而言的  men 为父 1 one 为母 1      
                for(i=0;i<length;i++)        
                    if(pop[men*length+i]==0)        
                    {         
                        if(flag1==point0){point1=i;break;}//记录路径在染色体 men 的起点位置 point1         
                        else          
                            flag1++; 
                    }       
                for(i=0;i<length;i++)        
                    if(pop[men*length+i]==0)        
                    {         
                        if(flag2==(point0+1)){point2=i;break;}//记录路径在染色体 men 的终点位置 point2         
                        else          
                            flag2++;        
                    }          
                for(i=0;i<length;i++)//保存父代染色体       
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
    }//##############################交叉运算交叉完成 
} 

static void mutation(void) 
{   
    int i,t;   
    int point1,point2,flag1;   
    double y;   
    for(i=0;i<popsize;i++)//随机选择两个编号进行交换（变异操作）------------2 交换变异   
    {        
        y=(double)rand()/(RAND_MAX);//随机产生一个数，与变异率 pm 进行比较       
        if(y<pm)    
        {        
            do     
            {         
                flag1=0;         
                point1=randnum1(length-3);//产生一个 1-15 的整数，变异交换的第一点         
                point2=randnum1(length-3);//产生一个 1-15 的整数，变异交换的第二点          
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
    } //##############################变异运算结束  
}    

void report(void)   
{    
    int i;    
    cout<<"game over"<<endl;      
    cout<<"the best value is "<<bestval<<endl;    
    cout<<"the num is :";       
    for(i=0;i<length;i++)    
    {popx1[i]=pop[popsize*length+i];}    
    cout<<"  该线路最短的路程为：  "<<leng(popx1)<<endl;    
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
    FILE *infile;//读取需求点信息   
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
        cout<<"第  "<<gen<<"代最佳为    "<<bestval<<endl;     
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
        cout<<"二轮：第  "<<gen<<"代最佳为    "<<bestval<<endl;   
        max2=bestval;   
        if(abs(max1-max2)<=10 && max1<=6900) 
            break; 
    }      
    report();   
    cout<<maxgener<<endl;  
    return 0;
} 
