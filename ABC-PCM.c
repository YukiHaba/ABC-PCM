#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

#define PI 3.141592654

int poisson_rnd(double);
int poisson_rnd(double lambda){
    int k;
    
    lambda = lambda+log(drand48());
    k=0;
    while (lambda>0){
        lambda += log(drand48());
        k++;
    }
    return k;
}


// define tree //

// # of nodes including tips = 19 
#define N 19

void molerat_data();
double phenotype[N];
double cal_LL(double, double, double, double, double, double);

// define each node's ancestor and descendant in the tree. b_a = ancestor's node ; b_d = descendant's node //
int b_a[N] = {99, 0, 1, 2, 3, 4, 5, 5, 4, 3, 9, 9, 2, 12, 13, 13, 12, 1, 0};
int b_d[N] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
// branch length (M years)
double bl[N] = {9999, 7.1, 13, 8.7, 8.8, 0.9, 5.5, 5.5, 0.4, 6.0, 6.0, 15.2, 6.4, 16.5, 1, 1, 17.5, 36.9, 44};
// tip position. 1 = the node is a tip; 999 = the node is not a tip
int tip[N] =  {999, 999, 999, 999, 999, 999, 1, 1, 999, 1, 1, 1, 999, 999, 1, 1, 1, 1, 1};





// phenotype data of tips //
double spp_m[N], spp_sd[N];
int spp_soc[N]; // 0 = solitary; 1 = social

// group size 
#define UNK_SD 2.56

void molerat_data(){
    // values updated
    spp_m[6] = 11; spp_sd[6] = 6.26; //F_damarensis
    spp_m[7] = 8.72; spp_sd[7] = 2.15; //F_anselli
    spp_m[9] = 7; spp_sd[9] = 2.42; //F_darlingi
    spp_m[10] = 9.906; spp_sd[10] = 2.493338; //F_mechowi
    spp_m[11] = 5.163; spp_sd[11] = 2.6222; //C_h_hottentotus
    spp_m[14] = 1; spp_sd[14] = UNK_SD; //B_janetta
    spp_m[15] = 1; spp_sd[15] = UNK_SD; //B_suillus
    spp_m[16] = 1; spp_sd[16] = UNK_SD; //G_capensis
    spp_m[17] = 1; spp_sd[17] = UNK_SD; //H_argenteocinereus
    spp_m[18] = 75; spp_sd[18] = 48.65; //H_glaber
}

// sociality 
// 1 = social species; 0 = solitary species
int social[N] = {999, 999, 999, 999, 999, 999, 1, 1, 999, 1, 1, 1, 999, 999, 0, 0, 0, 0, 1};






//// ABC-PCM ////

// prior //
double MRCA;
double MIN_MRCA = 1.001;
double MAX_MRCA = 20;

double ev_rate;
double MIN_ev_rate = 0.001;
double MAX_ev_rate = 30.;

double eu_1;
double MIN_eu_1 = 0.001;
double MAX_eu_1 = 30.;

double eu_2;
double MIN_eu_2 = 0.001;
double MAX_eu_2 = 5.;

double mutation_step = 0.1;


// parameters for link functions. social --> solitary (q),  solitary --> social (p) //
double s_a; // curvature coefficient of p
double MAX_s_a = 1;
double MIN_s_a = 0.1; //arbitary value, but consiering the biological reality

double s_b; // curvature coefficient of q
double MAX_s_b = 1;
double MIN_s_b = 0.1; //arbitary value, but consiering the biological reality

double p(double phenotype, double s_a); //phenotype, a coefficient (m) transition prob. of solitary --> social

double p(double phenotype, double s_a){
    double temp2;
    double exp(double);
    temp2 = -exp(s_a*(phenotype-1))+1; // it's the refrection about y=0.5 of q
    return temp2;
}

double q(double, double); //phenotype, a coefficient (m) transition prob. of social --> solitary

double q(double phenotype, double s_b){
    double temp1;
    double exp(double);
    temp1 = exp(-s_b*(phenotype-1));
    return temp1;
}



// calculate liklihood //
double cal_LL(double MRCA, double ev_rate, double eu_1, double eu_2, double s_a, double s_b){
    
    int i, j, d_plus, d_minus;
    int x;
    int temp_plus, temp_minus, temp_total;
    double lik = 0;
    int temp_sol, temp_soc;
    double temp_phenotype;
    
    for (i=0; i<N; i++) phenotype[i] = 0;
    
    // setting spp i=0
    phenotype[0] = MRCA;
    if (0.5 > drand48()) spp_soc[0] = 0;
    else spp_soc[0] = 1;
    
    
    // passing MRCA to ancestors of two branches
    for (j=0; j<N; j++) {
        if (b_d[0]==b_a[j]) {
            phenotype[j] = phenotype[0];
            spp_soc[j] = spp_soc[0];
            
        }
    }
    
    // simulation
    for (i=1; i<N; i++) {    // note we start from i=1, not 0
        
        // setting # of evolutionary changes
        if (i==18) {
            temp_plus = poisson_rnd( ev_rate * bl[i] * eu_1 );
            temp_minus = poisson_rnd( ev_rate* bl[i] / eu_1 );
        }
        else if (i==6) {
            temp_plus = poisson_rnd( ev_rate * bl[i] * eu_2 );
            temp_minus = poisson_rnd( ev_rate * bl[i] / eu_2 );
        }
        else {
            temp_plus = poisson_rnd( ev_rate * bl[i] );
            temp_minus = poisson_rnd( ev_rate * bl[i] );
        }
        
        temp_total = temp_plus + temp_minus;
        
        // trait evolution at each branch
        for (x=0; x<temp_total; x++) {
            
            temp_phenotype = phenotype[i];
            
            if (spp_soc[i] == 1) {      // BM (social spp)
                
                if (temp_plus>0 && temp_minus>0) {
                    if (drand48()<0.5){
                        d_plus=1; d_minus=0; temp_plus--;
                    } else {
                        d_plus=0; d_minus=1; temp_minus--;
                    }
                } else if (temp_plus==0 && temp_minus>0){
                    d_plus=0; d_minus=1; temp_minus--;
                } else if (temp_plus>0 && temp_minus==0){
                    d_plus=1; d_minus=0; temp_plus--;
                }
                
                phenotype[i] += d_plus*mutation_step - d_minus*mutation_step;
                
                if (phenotype[i] < 1) {
                    
                    phenotype[i] = temp_phenotype;
                    
                }
                
                if (q(phenotype[i], s_b) > drand48()) spp_soc[i] = 0; // whether social spp transits to solitary spp
                
            }
            
            else if (spp_soc[i] == 0) { // BM (solitary spp)
                
                if (temp_plus>0 && temp_minus>0) {
                    if (drand48()<0.5){
                        d_plus=1; d_minus=0; temp_plus--;
                    } else {
                        d_plus=0; d_minus=1; temp_minus--;
                    }
                } else if (temp_plus==0 && temp_minus>0){
                    d_plus=0; d_minus=1; temp_minus--;
                } else if (temp_plus>0 && temp_minus==0){
                    d_plus=1; d_minus=0; temp_plus--;
                }
                
                phenotype[i] += d_plus*mutation_step - d_minus*mutation_step;
                
                if (phenotype[i] < 1) {
                    
                    phenotype[i] = temp_phenotype; // if the group size become less than 1, it remains unchanged 
                    
                }
                
                if (p(phenotype[i], s_a) > drand48()) spp_soc[i] = 1;   // whether solitary spp transits to social spp
                
                
            }
        }
        
        for (j=0; j<N; j++) {
            if (b_d[i]==b_a[j]) {
                phenotype[j] = phenotype[i];
                spp_soc[j] = spp_soc[i];
            }
        }
    }
    // end of the i loop, end of trait evolution
    
    // checking distribution of solitary and social spp.
    temp_sol = 0;
    temp_soc = 0;
    
    //for (i=0; i<N; i++) printf("%d\t", spp_soc[i]);
    
    for (i=0; i<N; i++) {
        if (tip[i]==1 && social[i]==0){
            if (spp_soc[i] == 0) temp_sol++;
        }
        else if (tip[i]==1 && social[i]==1){
            if (spp_soc[i] == 1) temp_soc++;
        }
    }
    
    if (temp_sol == 4 && temp_soc == 6){
        for (i=0; i<N; i++) {
            if (tip[i]==1 && social[i]==1){
                lik += ( -(phenotype[i]-spp_m[i])*(phenotype[i]-spp_m[i]) / (2*spp_sd[i]*spp_sd[i]) )
                ;
                // note "lik" is not a "absolute" but "relative" likelihood
            }
        }
        
    }
    
    
    else lik = -9999999999;
    
    return lik;
    
}
// calculate liklihood - end //

/// ABC-PCM - end ///




// if accepted, write results //
int main(void)
{
    
    molerat_data();
    
    // name of a output file
    struct tm result;
    time_t timep;
    struct timeb timeb;
    FILE *fp;
    char fname[80];
    
    ftime(&timeb);
    timep = timeb.time;
    localtime_r(&timep, &result);
    snprintf(fname, 80, "ABC_MR_%02d%02d%02d_%02d%02d%02d_eu1_2.txt",
             result.tm_year-100,result.tm_mon+1, result.tm_mday, result.tm_hour, result.tm_min, result.tm_sec);
    fp = fopen(fname, "w");
    
    freopen (fname, "a", fp);
    //end
    
    
    // title of this program
    time_t timer;
    struct tm *t_st;
    time(&timer);
    
    printf("\n*** %s*** Phylogenetic model of phenotypic evolution version. 1.0 by Yuki Haba and Nobuyuki Kutsukake ***\n\n", ctime(&timer));
    //fprintf(fp, "\n*** %s*** Phylogenetic model of phenotypic evolution version. 1.0 by Nobuyuki Kutsukake ***\n\n", ctime(&timer));
    //end
    
    
    
    int i;
    srand48(time(NULL));
    double ABC_LL;
    double pro = 0.1; // this is an adujusting value for efficient computation because model likelihood is far less than a random number from U[0,1]
    double rejection;
    double bar;
    int count, count_pro;
    
    printf("MRCA\tMRCAsociality\tev_rate\teu1\teu2\ts_a\ts_b\tABC_LL\tbar\trejection\tTorF\tphenotypes\tsocialities\n");
    fprintf(fp, "MRCA\tMRCAsociality\tev_rate\teu1\teu2\ts_a\ts_b\tABC_LL\tbar\trejection\tTorF\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\tp9\tp10\tp11\tp12\tp13\tp14\tp15\tp16\tp17\tp18\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\ts10\ts11\ts12\ts13\ts14\ts15\ts16\ts17\ts18\n");
    
    for (i=0; ; ) {
        MRCA = (MAX_MRCA-MIN_MRCA)*drand48()+MIN_MRCA;
        ev_rate = (MAX_ev_rate - MIN_ev_rate) * drand48() + MIN_ev_rate;
        eu_1 = (MAX_eu_1 - MIN_eu_1) * drand48() + MIN_eu_1;
        eu_2 = (MAX_eu_2 - MIN_eu_2) * drand48() + MIN_eu_2;
        s_a = (MAX_s_a - MIN_s_a) * drand48() + MIN_s_a;
        s_b = (MAX_s_b - MIN_s_b) * drand48() + MIN_s_b;
        
        
        //printf("%2.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t%1.5lf\t%1.5lf\n", MRCA, ev_rate, eu_1, eu_2, p, s_a);
        
        bar = drand48();
        ABC_LL = cal_LL(MRCA, ev_rate, eu_1, eu_2, s_a,s_b);
        rejection = log(bar)-pro;
        
        //printf("%.20lf\n", record_q);
        
        
        //printf("%2.2lf\n", ABC_LL);
        
        if (ABC_LL != -9999999999){
            if (rejection<ABC_LL) {
                printf("%2.4lf\t%d\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t", MRCA, spp_soc[0], ev_rate, eu_1, eu_2, s_a, s_b, ABC_LL, bar, rejection);
                fprintf(fp, "%2.4lf\t%d\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t", MRCA, spp_soc[0], ev_rate, eu_1, eu_2, s_a, s_b, ABC_LL, bar, rejection);
                
                // print rejection TorF //
                if (-1*pro < ABC_LL) {
                    printf("1\t"); fprintf(fp, "1\t");
                    count_pro++;
                }
                
                else {
                    printf("0\t"); fprintf(fp, "0\t");
                }
                
                //record all phenotype to see intermedeate phenotypes//
                for (i=1; i<N; i++) {
                    fprintf(fp, "%2.4lf\t", phenotype[i]);
                }
                for (i=1; i<N; i++) {
                    printf("%2.4lf\t", phenotype[i]);
                }
                
                //record final phenotype//
                fprintf(fp, "%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t", phenotype[5],phenotype[7],phenotype[9],phenotype[10],phenotype[14],phenotype[15],phenotype[16],phenotype[17],phenotype[20],phenotype[21],phenotype[22],phenotype[23],phenotype[24]);
                 printf("%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t%2.2lf\t", phenotype[5],phenotype[7],phenotype[9],phenotype[10],phenotype[14],phenotype[15],phenotype[16],phenotype[17],phenotype[20],phenotype[21],phenotype[22],phenotype[23],phenotype[24]);
                
                
                // record all branch sociality : to see social evolution process //
                for (i=1; i<N; i++) {
                    fprintf(fp, "%d\t", spp_soc[i]);
                }
                for (i=1; i<N; i++) {
                    printf("%d\t", spp_soc[i]);
                }
                
                printf("\n");
                fprintf(fp, "\n");
                
                count++;
                //if (count_pro==25) break;
            }
        }
        
        // keep running until it accepts 500 successful simulations
        if (count==500) break;
    }
    
    fclose (fp);
    
}