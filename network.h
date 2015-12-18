//テスト

/* !ATTENTION! : This program run only on and over C99. */
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

/* network.h */
#ifndef NETWORK_H
#define NETWORK_H



/* ライブラリ一覧 */
//！注意！
//	☆隣接行列は重みも含めた計算ができるようにdouble型にしています．
//	☆有向グラフにおけるクラスター係数の計算は無向と同じように計算されます（つまりリンクがあればそれを一辺と数えます）．
//	☆「非連結なグラフ」における平均頂点間距離は正しく計算はされません．
//	☆関数vertex_distanceに重み付き隣接行列を引数として渡すと頂点uから頂点vまでリンクの重みを含めた最短経路を計算します．ステップ数は計算されません．
//	ステップを計算したいときは1と0で記述された隣接行列を渡してください．
//	☆グラフ生成関数は全て無向グラフです．
//	☆グラフ描画ファイルは極めて基本的なものしか作成できません．
//	☆次数分散や次数分布は無向グラフのみ計算が保証されています．

//ある頂点に対するクラスター係数計算関数clusterの使い方 -> cluster(調べたい頂点番号，総ノード数，隣接行列の配列へのポインタ)
double cluster(int node_label, int node_number, double adjancency_matrix[node_number][node_number]);

//あるグラフに対するクラスター係数計算関数av_cllusterの使い方 -> av_cluster(総ノード数，隣接行列の配列へのポインタ)
double av_cluster(int node_number, double adjancency_matrix[node_number][node_number]);

//頂点uから頂点vまでの距離を計算する関数vertex_distanceの使い方 -> vertex_distance(頂点u，頂点v，総ノード数，隣接行列の配列へのポインタ)
double vertex_distance(int node_u, int node_v, int node_number, double adjancency_matrix[node_number][node_number]);

//あるグラフに対する平均頂点間距離を計算する関数av_vertex_distanceの使い方 -> av_vertex_distance(総ノード数，隣接行列の配列へのポインタ)
double av_vertex_distance(int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列に正方格子グラフを適応する関数square_latticeの使い方 -> square_lattice(総ノード数，隣接行列の配列へのポインタ)
void square_lattice(int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列にレギュラーグラフを適応する関数regular_latticeの使い方 -> regular_lattice(平均次数，総ノード数，隣接行列の配列へのポインタ)
void regular_lattice(int average_degree, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列にランダムネットワークを適応する関数random_networkの使い方 -> random_network(平均次数，総ノード数，隣接行列の配列へのポインタ)
void random_network(int average_degree, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列にWSモデルを適応する関数WS_modelの使い方 -> WS_model(平均次数，リンクの張替え確率，総ノード数，隣接行列の配列へのポインタ)
void WS_model(int average_degree, double replacing_parameter, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列にBAモデルを適応する関数BA_modelの使い方 -> BA_model(いくつの頂点から始めるか，追加していく頂点数，総ノード数，隣接行列の配列へのポインタ)
void BA_model(int start_node_number, int add_link_number, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列を元にpajekファイルを生成する関数pajek_graphの使い方 -> pajek_graph(文字列(ファイル名)を指すポインタ．文字列(directもしくはundirected)，総ノード数，隣接行列の配列へのポインタ)
void pajek_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列を元にgephiファイルを生成する関数gephi_graphの使い方 -> gephi_graph(文字列(ファイル名)を指すポインタ．文字列(directもしくはundirected)，総ノード数，隣接行列の配列へのポインタ)
void gephi_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列を元にgephiファイルを生成する関数gephi_graphの使い方 -> sociarium_graph(文字列(ファイル名)を指すポインタ．文字列(directもしくはundirected)，総ノード数，隣接行列の配列へのポインタ)
void sociarium_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列を元に次数分散を計算する関数degree_varianceの使い方 -> degree_variance(総ノード数，隣接行列の配列へのポインタ)
double degree_variance(int node_number, double adjancency_matrix[node_number][node_number]);

//渡された隣接行列を元に次数分布(CSV)を出力する関数degree_distributionの使い方 -> degree_distribution(文字列(ファイル名)を指すポインタ，総ノード数，隣接行列の配列へのポインタ)
void degree_distribution(char *file_name, int node_number, double adjancency_matrix[node_number][node_number]);





// MT //
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
// MT end //

// MT //
/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

// MT end //


//calculate cluster confficient function
double cluster(int node_label, int node_number, double adjancency_matrix[node_number][node_number]){
	int j,k;
	int count=0;
	int sum=0;
	int degree_num=degree_count(node_label,node_number,adjancency_matrix);
	double matrix_dummy[node_number][node_number];
	
	for(j=0;j<node_number;j++){
		for(k=0;k<node_number;k++){
			matrix_dummy[j][k] = adjancency_matrix[j][k];
		}
	}
	
	for(j=0;j<node_number;j++){
		for(k=0;k<=j;k++){
			if((adjancency_matrix[j][k] > 0) && (adjancency_matrix[k][j] > 0)){
				matrix_dummy[j][k] = 0;
			}
		}
	}
		
	//Error Output
	if((node_label < 0) || (node_number < node_label)){
		printf("[network.h ERROR] cluster : first argument is over\n");
		exit(1);
	}
		
	if(degree_num<2){	//次数が0か1であればクラスターは形成できないので0を返す
		return 0;
	}
	
	for(j=0;j<node_number;j++){
		if((matrix_dummy[node_label][j] > 0) || (matrix_dummy[j][node_label] > 0)){
			for(k=0;k<node_number;k++){
				if((matrix_dummy[node_label][k]) || (matrix_dummy[k][node_label] > 0)){
					if((matrix_dummy[j][k] > 0) || (matrix_dummy[j][k] > 0)){
						sum++;
					}
				}
			}
		}
	}
		
	return (2.0 / (degree_num * (degree_num - 1))) * sum;
}

int degree_count(int node_i, int node_number, double adjancency_matrix[node_number][node_number]){
	int i;
	int count=0;
	
	for(i=0;i<node_number;i++){
		if((node_i != i) && (adjancency_matrix[node_i][i] > 0) || (adjancency_matrix[i][node_i] > 0)){
			count++;	//count indegree or outdegree
		}
	}
	
	return count;
}

//calculate average of cluster confficient function
double av_cluster(int node_number, double adjancency_matrix[node_number][node_number]){
	int node_i;
	double sum=0, average=0;
	
	for(node_i=0;node_i<node_number;node_i++){
		sum = cluster(node_i, node_number, adjancency_matrix) + sum;
	}
	
	average = sum / node_number;
	
	return average;
}

void int_initialization(int *p,int x){	//pは1次元配列のポインタ，xは配列の大きさ(2次元にも応用可 -> *pに配列の行を渡せば良い)
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

void double_initialization(double *p,int x){	//pは1次元配列のポインタ，xは配列の大きさ(2次元にも応用可 -> *pに配列の行を渡せば良い)
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

double dijkstra(int node_u, int node_v, int *ver, double *dis, int node_number, double adjancency_matrix[node_number][node_number]){	//ダイクストラ法の本計算プログラム
	int i;
	int queue[node_number];	//ある頂点に関する隣接ノードが格納される
	int head = 0, tail = 0;
	double dummy;
	
	while(vertex_list_check(ver,node_number)){
		dummy = DBL_MAX;
		
		for(i=0;i<node_number;i++){
			if(ver[i]==0){
				if(dummy >= dis[i]){
					dummy = dis[i];
				}
			}
		}
		
		if(dummy == DBL_MAX){	//もしも最小値がINT_MAXであればこれ以上探索する意味は無い（->不可能到達点しかない）
			break;
		}
		
		for(i=0;i<node_number;i++){
			if(dummy == dis[i]){
				ver[i] = 1;
				queue[tail] = i;
				tail++;
			}
		}
		
		while(head != tail){
			node_u = queue[head];
			head++;
			
			for(i=0;i<node_number;i++){	//変数 i には隣接ノードが入る
				if(adjancency_matrix[node_u][i]>0){
					if(dis[i] > (dis[node_u] + adjancency_matrix[node_u][i])){
						dis[i] = dis[node_u] + adjancency_matrix[node_u][i];
					}
				}
			}
		}
		//初期化
		int_initialization(queue,node_number);
		head = 0;
		tail = 0;
	}
	
	return dis[node_v];
}

double vertex_distance(int node_u, int node_v, int node_number, double adjancency_matrix[node_number][node_number]){	//頂点uから頂点vまでの最短頂点間距離を計算し，その距離を戻り値とする
	int vertex_list[node_number];
	double distance_list[node_number];
	int prev = node_u;
	int dummy;
	int i;
	
	//初期化
	for(i=0;i<node_number;i++){
		//distance_listの初期化
		if(i==node_u){
			distance_list[i] = 0;
		}
		else{
			distance_list[i] = DBL_MAX;
		}
		
		//vertex_listの初期化
		vertex_list[i] = 0;	//vertex_list[i]==1ならば頂点 i は頂点リストから除外されているものとする(0ならば頂点リストの要素)
	}
	//ここまで初期化
	
	//本計算
	return dijkstra(node_u, node_v, vertex_list, distance_list, node_number, adjancency_matrix);
}

int vertex_list_check(int *p, int array_length){	//配列の中身がすべて1であれば0を返し，それ以外であれば1を返す
	int i;
	
	for(i=0;i<array_length;i++){
		if(*p==0){
			return 1;
		}
		p++;
	}
	return 0;	//*p が全部 1 ならばreturn 0
}

double av_vertex_distance(int node_number, double adjancency_matrix[node_number][node_number]){	//平均頂点間距離を計算するプログラム(d(i,j)がINFの場合は数に入れない)
	int count=0;	//頂点間距離がINFの頂点の数を数える
	double sum=0;
	double d;
	int i,j;
	double av_dis;
	
	for(i=0;i<node_number;i++){
		for(j=0;j<=i;j++){
			d = vertex_distance(i,j,node_number,adjancency_matrix);
			
			if(DBL_MAX == d){
				count++;
			}
			else{
				sum = d + sum;
			}
		}
	}	
	av_dis = 2.0/(node_number*(node_number-1) ) * sum;
	
	return av_dis;
}

void square_lattice(int node_number, double adjancency_matrix[node_number][node_number]){
	int i,j;
	
	if((double)((int)sqrt(node_number)) != sqrt(node_number)){
		printf("[network.h ERROR] square_lattice : second argument is not correct (please square root of node number is integer)\n");
		exit(1);
	}
	
	for(j=0;j<node_number;j++){
		double_initialization(adjancency_matrix[j],node_number);
	}
	
	for(i=0;i<node_number;i++){
		for(j=i;j<node_number;j++){
			if(i == j){
				adjancency_matrix[i][j] = 0;
			}
			else{
				if( ((i+1) == j) && ((i+1)%(int)sqrt(node_number) != 0)){
					adjancency_matrix[i][j] = 1;	//右にリンクを伸ばす
					adjancency_matrix[j][i] = 1;
				}
				if( ((j-i) == (int)sqrt(node_number)) && (i+(int)sqrt(node_number) < node_number)){
					adjancency_matrix[i][i+(int)sqrt(node_number)] = 1;	//下にリンクを伸ばす
					adjancency_matrix[i+(int)sqrt(node_number)][i] = 1;
				}
			}
		}
	}
}

void regular_lattice(int average_degree, int node_number, double adjancency_matrix[node_number][node_number]){
	WS_model(average_degree, 0.0, node_number, adjancency_matrix);
}

void random_network(int average_degree, int node_number, double adjancency_matrix[node_number][node_number]){
	WS_model(average_degree, 1.0, node_number, adjancency_matrix);
}

void WS_model(int average_degree, double replacing_parameter, int node_number, double adjancency_matrix[node_number][node_number]){
	int i,j;
	int count=0;
	int node_i;
	
	init_genrand((unsigned)time(NULL));
	
	if(average_degree%2){
		printf("[network.h ERROR] WS_model : first argument is even \n");
		exit(1);
	}
	
	if(node_number < average_degree){
		printf("[network.h ERROR] WS_model : first argument is over \n");
		exit(1);
	}
	
	if((replacing_parameter < 0) || (1 < replacing_parameter)){
		printf("[network.h ERROR] WS_model : second argument is over (please \" from 0.0 to 1.0\")\n");
		exit(1);
	}
	
	for(j=0;j<node_number;j++){
		double_initialization(adjancency_matrix[j],node_number);
	}
	
	for(i=0;i<node_number;i++){
		for(j=i-(average_degree/2);j<=i+(average_degree/2);j++){
			if(i!=j){
				if(j < 0){
					adjancency_matrix[i][j+node_number] = 1;
					adjancency_matrix[j+node_number][i] = 1;
				}
				else if(node_number-1<j){
					adjancency_matrix[i][j-node_number] = 1;
					adjancency_matrix[j-node_number][i] = 1;
				}
				else{
					adjancency_matrix[i][j] = 1;
					adjancency_matrix[j][i] = 1;
				}
			}
		}
	}
	
	for(i=0;i<node_number;i++){
		for(j=0;j<=i;j++){
			if((adjancency_matrix[i][j]==1) && (genrand_real1()<=replacing_parameter)){
				adjancency_matrix[i][j] = 0;
				adjancency_matrix[j][i] = 0;
				
				if(genrand_real1()<0.5){
					node_i = (int)(genrand_real1()*(double)(node_number-1));
					
					while((adjancency_matrix[i][node_i]==1) || (adjancency_matrix[node_i][i]==1)){
						node_i = (int)(genrand_real1()*(double)(node_number-1));
					}
					
					adjancency_matrix[i][node_i] = 1;
					adjancency_matrix[node_i][i] = 1;
				}
				else{
					node_i = (int)(genrand_real1()*(double)(node_number-1));
					
					while((adjancency_matrix[j][node_i]==1) || (adjancency_matrix[node_i][j]==1)){
						node_i = (int)(genrand_real1()*(double)(node_number-1));
					}
					
					adjancency_matrix[j][node_i] = 1;
					adjancency_matrix[node_i][j] = 1;
				}
			}
		}
	}
}

void BA_model(int start_node_number, int add_link_number, int node_number, double adjancency_matrix[node_number][node_number]){
	int check=0;
	double sum_degree;
	double stretch_of_probability[node_number];
	int vertex_mark[node_number];
	int i,j;
	int kk,mm;
	double stock_rand;
	
	init_genrand((unsigned)time(NULL));
	
	if(start_node_number < 2){
		printf("[network.h ERROR] BA_model : first argument is too few (please \"from 2 to %d \")\n",node_number-1);
		exit(1);
	}
	
	for(j=0;j<node_number;j++){
		double_initialization(adjancency_matrix[j],node_number);
	}
	
	//first graph
	while(check==0){
		for(i=1;i<start_node_number;i++){
			for(j=0;j<i;j++){
				if(genrand_real1()>0.5){
					adjancency_matrix[i][j] = 1;
					adjancency_matrix[j][i] = 1;
				}
			}
		}
		
		for(i=0;i<start_node_number;i++){
			if(degree_count(i,node_number,adjancency_matrix)<1){
				for(j=0;j<start_node_number;j++){
					double_initialization(adjancency_matrix[j],start_node_number);
				}
				break;
			}
			else if(i==start_node_number-1){
				check=1;	//すべての頂点 i の次数が1以上ならばwhile文を抜ける
			}
		}
	}
	
	if((add_link_number < 1) || (start_node_number < add_link_number)){
		printf("[network.h ERROR] BA_model : second argument is over (please \" from 1 to %d\")\n",start_node_number);
		exit(1);
	}
	
	for(i=start_node_number;i<node_number;i++){	//変数 i はこれから加えていく頂点
		//初期化
		sum_degree=0;
		double_initialization(stretch_of_probability,i);
		int_initialization(vertex_mark,i);
		
		for(j=0;j<i;j++){
			sum_degree = degree_count(j,node_number,adjancency_matrix) + sum_degree;
		}
		
		for(j=0;j<i;j++){
			stretch_of_probability[j] = degree_count(j,node_number,adjancency_matrix) / sum_degree;
			
			if(j>0){
				stretch_of_probability[j] = stretch_of_probability[j] + stretch_of_probability[j-1];
			}
		}
		
		while(degree_count(i,node_number,adjancency_matrix) < add_link_number){
			stock_rand = genrand_real1();
			
			for(j=0;j<i;j++){
				if((j==0) && (stock_rand <= stretch_of_probability[j])){
					//0が選ばれました
					adjancency_matrix[i][j] = 1;
					adjancency_matrix[j][i] = 1;
				}
				else if((stretch_of_probability[j-1] < stock_rand) && (stock_rand <= stretch_of_probability[j])){
					//jが選ばれました
					adjancency_matrix[i][j] = 1;
					adjancency_matrix[j][i] = 1;
				}
			}
		}
	}
}

void pajek_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]){
	FILE *fp;
	char string[100];
	int i,j;
	
	if((strcmp(directed_or_undirected,"directed")!=0) && (strcmp(directed_or_undirected,"undirected")!=0)){
		printf("[network.h ERROR] pajek_graph : second argument is \"directed\" or \"undirected\" of string \n");
		exit(1);
	}
	
	sprintf(string,"%s.net",file_name);
	
	/* file open */
	if ((fp = fopen(string, "wt")) == NULL) {
		printf("file open error \n");
		exit(EXIT_FAILURE);	/* エラーの場合は通常、異常終了する */
	}
	
	/* writing */
	fprintf(fp,"*Vertices %d \n",node_number);
	
	for(i=0;i<node_number;i++){
		fprintf(fp,"%d \"agent_%d\" \n",i+1,i+1);
	}
	
	if(strcmp(directed_or_undirected,"directed")==0){
		fprintf(fp,"*Arcs \n");
		
		for(i=0;i<node_number;i++){
			for(j=0;j<node_number;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"%d %d %f \n",i+1,j+1,adjancency_matrix[i][j]);
				}
			}
		}
	}
	else if(strcmp(directed_or_undirected,"undirected")==0){
		fprintf(fp,"*Edges \n");
		
		for(i=0;i<node_number;i++){
			for(j=0;j<=i;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"%d %d %f \n",i+1,j+1,adjancency_matrix[i][j]);
				}
			}
		}
	}
	
	/* file close */
	fclose(fp);
}

void gephi_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]){
	FILE *fp;
	char string[100];
	int edge=0;
	int i,j;
	
	if((strcmp(directed_or_undirected,"directed")!=0) && (strcmp(directed_or_undirected,"undirected")!=0)){
		printf("[network.h ERROR] gephi_graph : second argument is \"directed\" or \"undirected\" of string \n");
		exit(1);
	}
	
	/* gephi file作成 */
	sprintf(string,"%s.gexf",file_name); /*ファイル名整形*/
	
	/* file open */
	if ((fp = fopen(string, "wt")) == NULL) {
		printf("file open error \n");
		exit(EXIT_FAILURE);	/* エラーの場合は通常、異常終了する */
	}
	
	/* writing */
	fprintf(fp,"<gexf \n");
	fprintf(fp,"xmlns:ns0=\"http://www.gexf.net/1.1draft/viz\" \n");
	fprintf(fp,"version=\"1.1\" \n");
	fprintf(fp,"xmlns=\"http://www.gexf.net/1.1draft\" \n");
	fprintf(fp,"xmlns:viz=\"http://www.gexf.net/1.1draft/viz\" \n");
	fprintf(fp,"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n");
	fprintf(fp,"xsi:schemaLocation=\"http://www.w3.org/2001/XMLSchema-instance\"> \n");
	
	if(strcmp(directed_or_undirected,"directed")==0){
		fprintf(fp,"<graph defaultedgetype=\"directed\" idtype=\"string\" type=\"static\">\n");
	}
	else if(strcmp(directed_or_undirected,"undirected")==0){
		fprintf(fp,"<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">\n");
	}
	
	fprintf(fp,"<nodes>\n");
	
	for(i=0;i<node_number;i++){
		fprintf(fp,"<node id=\"%d\" label=\"agent_%d\" >\n",i,i);
		fprintf(fp,"</node>\n");
	}
	
	fprintf(fp,"</nodes>\n");
	fprintf(fp,"<edges>\n");
	
	if(strcmp(directed_or_undirected,"directed")==0){
		for(i=0;i<node_number;i++){
			for(j=0;j<node_number;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%f\"/>\n",edge,i,j,adjancency_matrix[i][j]);
					edge++;
				}
			}
		}
	}
	else if(strcmp(directed_or_undirected,"undirected")==0){
		for(i=0;i<node_number;i++){
			for(j=0;j<=i;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%f\"/>\n",edge,i,j,adjancency_matrix[i][j]);
					edge++;
				}
			}
		}
	}
	
	fprintf(fp,"</edges>\n");
	fprintf(fp,"</graph>\n");
	fprintf(fp,"</gexf>\n");
	
	/* file close */
	fclose(fp);
}

void sociarium_graph(char *file_name, char *directed_or_undirected, int node_number, double adjancency_matrix[node_number][node_number]){
	FILE *fp;
	char string[100];
	int i,j;
	
	sprintf(string,"%s.txt",file_name); /*ファイル名整形*/
	
	/* file open */
	if ((fp = fopen(string, "wt")) == NULL) {
		printf("file open error \n");
		exit(EXIT_FAILURE);	/* エラーの場合は通常、異常終了する */
	}
	
	/* writing */
	fprintf(fp,"# Pajek形式のデータ \n");
	fprintf(fp,"# 参照: http://www.tp.umu.se/~rosvall/code.html \n");
	fprintf(fp,"\n");
	fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
	//fprintf(fp,"@title = The dissemination of culture \n");
	
	if(strcmp(directed_or_undirected,"directed")==0){
		fprintf(fp,"@directed \n");
	}
	else if(strcmp(directed_or_undirected,"undirected")==0){
		fprintf(fp,"@nondirected \n");
	}
	
	fprintf(fp,"\n");
	
	fprintf(fp,"*Vertices %d \n",node_number);
	
	fprintf(fp,"<nodes>\n");
	
	for(i=0;i<node_number;i++){
		fprintf(fp,"%d \"agent_%d\" \n",i,i);
	}
	
	fprintf(fp,"*Arcs \n");
	
	if(strcmp(directed_or_undirected,"directed")==0){
		for(i=0;i<node_number;i++){
			for(j=i;j<node_number;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"%d %d %f \n",i,j,adjancency_matrix[i][j]);
				}
			}
		}
	}
	else if(strcmp(directed_or_undirected,"undirected")==0){
		for(i=0;i<node_number;i++){
			for(j=i;j<=i;j++){
				if(adjancency_matrix[i][j]>0){
					fprintf(fp,"%d %d %f \n",i,j,adjancency_matrix[i][j]);
				}
			}
		}
	}
	
	/* file close */
	fclose(fp);
}

double variance(int node_number, int *array, double average){	//分散を計算する関数．avはaverage
	int i;
	double sum=0;
	
	for(i=0;i<node_number;i++){
		sum = sum + pow(average-*array,2);
		array++;
	}
	
	sum = sum / node_number;
	
	return sum;
}

double degree_variance(int node_number, double adjancency_matrix[node_number][node_number]){
	int degree_array[node_number];
	double average_degree=0;
	int node_i;
	
	for(node_i=0;node_i<node_number;node_i++){
		degree_array[node_i] = degree_count(node_i, node_number, adjancency_matrix);
		average_degree = degree_array[node_i] + average_degree;
	}
	average_degree = average_degree / node_number;
	
	return variance(node_number,degree_array,average_degree);
}

void degree_distribution(char *file_name, int node_number, double adjancency_matrix[node_number][node_number]){
	FILE *fp;
	int degree_array[node_number];
	int degree_distribution_array[node_number];
	char string[100];
	int node_i;
	int deg;
	int i;
	
	int_initialization(degree_array,node_number);
	int_initialization(degree_distribution_array,node_number);
	
	for(node_i=0;node_i<node_number;node_i++){
		degree_array[node_i] = degree_count(node_i, node_number, adjancency_matrix);
	}
	
	for(deg=0;deg<node_number;deg++){
		for(node_i=0;node_i<node_number;node_i++){
			if(deg==degree_array[node_i]){
				degree_distribution_array[deg]++;
			}
		}
	}
	
	sprintf(string,"%s.csv",file_name); /*ファイル名整形*/
	
	/* file open */
	if ((fp = fopen(string, "wt")) == NULL) {
		printf("file open error \n");
		exit(EXIT_FAILURE);	/* エラーの場合は通常、異常終了する */
	}
	
	fprintf(fp,"node degree,number of nodes\n");
	
	for(deg=0;deg<node_number;deg++){
		fprintf(fp,"%d,%d\n",deg,degree_distribution_array[deg]);
	}
	
	/* file close */
	fclose(fp);
}

#endif
