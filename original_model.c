/* ヘッダーファイル */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "MT.h"

/* プリプロセッサ */
#define lg 20	//格子の縦または横の長さ．
#define fe 3	//features．特徴の数
#define tr 20	//traits．特性の範囲
#define ag_num lg*lg	//agent数;
#define tmax 500000	//イベント数
#define wt_num 4	//書き出すファイルの数
#define rgb_conversion 255/tr	//agentにrgbで値を加える
#define make_torus 0	//格子モデルをトーラスにする場合は1にする(しない場合は0)
#define conservativeness_parameter 0.99	//保守性パラメータ(0以上1以下で数字を選ぶ) ->0(保守性を持たない状態)を選ぶとCentolaモデルになる
#define read_data_set 1	//data_setを読み込みdata setのある共同性に対する属性値(colum)の平均に近いグラフを生成する（読み込まない場合は0）
	#define data_set_row 100	//data setの共同性の0から1までの刻みの個数
	#define sigma_option 1	//グラフ生成における条件範囲 標準偏差(sigma)が1->68.27%，2->95.45%，3->99.73%
	#define SimulationStop 100	//シミュレーション回数．この回数を超えた場合強制終了します．
#define sociarium_ptch_width 5000	//sociariumのグラフをアニメーションにする際，時間軸を何単位で進めていきたいか

/* 構造体 */
typedef struct{
	int f[fe];	//特徴
	
} agent;


/* 外部変数 */
int lattice_mdl[ag_num][ag_num];	//格子モデルの隣接行列
int cd_nw[ag_num][ag_num];	//cultural driftするネットワークの隣接行列
int olap[ag_num][ag_num];	;//overlap
int wt[wt_num];	//中身はmain文に書いてある
agent ag[ag_num];
int node_through_list[ag_num];	//dfsで通過したノードは該当する列に1を入れる
int node_through_count;	//dfsで通過したノードをカウントする(dfsは再帰呼び出ししているのでこの変数は外部変数にしなきゃダメ)
int node_list_of_network[2];	//添字0の配列には「何個目の連結グラフか」を，添字1には「すべての連結グラフの中で最も大きいグラフのノード数」を格納
int max_size_network[ag_num][ag_num];	//連結グラフのうち最も多くのagent数をもつグラフを保存する


/* プロトタイプ宣言 */
void lattice_func(void);
void lattice_to_torus(void);
void tr_func(void);
int KD_func(int *p, int*q);
void olap_func(void);
void latcd_cpyfunc(void);
void initialization_1(int *p,int x);
void dfs1(int v);
void dfs2(int v);
int network_island_count(void);
int culture_count_function(void);
void max_size_network_func(void);
void orgnl_mdl(void);
double a_rate_of_bridge_func(void);
int vertex_list_check(int *p, int array_length);
int dijkstra(int node_u,int node_v,int *ver,int *dis,int length);
int shortest_path_func(int node_u, int node_v);
double average_vertex_distance(void);
double cluster_calculate_func(int node_i);
int degree_count(int node_i);
double average_cluster_calculate_func(void);
void read_data_set_func(char *filename, double *cooperation_parameter, double *island_number_av_array, double *agent_number_av_array, double *culture_number_av_array, double *cluster_number_av_array, double *shortest_path_av_array, double *bridge_number_av_array, double *island_number_variance_array, double *agent_number_variance_array, double *culture_number_variance_array, double *cluster_number_variance_array, double *shortest_path_variance_array, double *bridge_number_variance_array, double *island_number_standard_deviation_array,  double *agent_number_standard_deviation_array, double *culture_number_standard_deviation_array, double *cluster_number_standard_deviation_array, double *shortest_path_standard_deviation_array, double *bridge_number_standard_deviation_array);
int network_edges_count(void);


/* 外部関数 */
void lattice_func(void){	//格子モデル作成アルゴリズム
	int i,j;
	
	/* 初期化 */
	for(i=0;i<ag_num;i++){	//行列の縦
		initialization_1(lattice_mdl[i],ag_num);	//行列の横
		initialization_1(cd_nw[i],ag_num);	//行列の横
	}
	
	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(i == j){
				lattice_mdl[i][j] = 0;
			}
			else{
				if( ((i+1) == j) && ((i+1)%lg != 0)){
					lattice_mdl[i][j] = 1;	//右にリンクを伸ばす
					lattice_mdl[j][i] = 1;
				}
				if( ((j-i) == lg) && ((i+lg) < ag_num)){
					lattice_mdl[i][i+lg] = 1;	//下にリンクを伸ばす
					lattice_mdl[i+lg][i] = 1;
				}
			}
		}
	}
	
	if(make_torus){
		lattice_to_torus();
	}
}

void lattice_to_torus(void){
	int i,j;
	
	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(i == j){
				lattice_mdl[i][j] = 0;
			}
			else{
				if( ((j+1-i) == lg) && (i%lg==0)){
					lattice_mdl[i][j] = 1;	//横にリンクを伸ばす
					lattice_mdl[j][i] = 1;
				}
				if( ((j-i) == (lg*lg-lg)) && (i < lg)){
					lattice_mdl[i][j] = 1;	//縦にリンクを伸ばす
					lattice_mdl[j][i] = 1;
				}
			}
		}
	}
}

void tr_func(void){	//個々に特性を割り当てる
	int i,j;
	
	for(i=0;i<ag_num;i++){
		initialization_1(ag[i].f,fe);
	}
	
	for(i=0;i<ag_num;i++){
		
		for(j=0;j<fe;j++){
			ag[i].f[j] = (int)genrand_int31() % tr;	//[0,tr-1]までの間で乱数を取る
		}
	}
}

int KD_func(int *p, int*q){	//2つの配列pとqが互いに共有している文化の数を返す関数
	int i,n=0;
	
	for(i=0;i<fe;i++){
		if(*p == *q){
			n = n + 1;
		}
		p++;
		q++;
	}
	
	return n;
}

void olap_func(void){
	int i,j;
	
	for(i=0;i<ag_num;i++){	//行列の縦
		initialization_1(olap[i],ag_num);	//行列の横
	}
	
	for(i=0;i<ag_num;i++){
		for(j=0;j<ag_num;j++){
			if(lattice_mdl[i][j]==1){
				olap[i][j] = KD_func(ag[i].f,ag[j].f);
			}
		}
	}
}

void latcd_cpyfunc(void){	//レギュラー格子の隣接行列を文化流布して変化するのに使う隣接行列に要素をコピーする
	int i,j;
	
	for(i=0;i<ag_num;i++){
		for(j=0;j<ag_num;j++){
			cd_nw[i][j] = lattice_mdl[i][j];
		}
	}
}

//与えられた配列の中身を全て0にする
void initialization_1(int *p,int x){	//pは1次元配列のポインタ，xは配列の大きさ
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

//深さ優先探索
void dfs1(int v){
	int i;
	
	node_through_list[v] = 1;
	
	for(i=0;i<ag_num;i++){
		if((cd_nw[v][i] == 1 )&& (node_through_list[i] == 0)){
			dfs1(i);
		}
	}
}

//深さ優先探索関数(戻り値なしver(最も大きな連結グラフも同時作成)
void dfs2(int v){
	int i;
	
	node_through_list[v] = 1;
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[v][i] == 1){
			max_size_network[v][i] = 1;
			
			if(node_through_list[i] == 0){
				dfs2(i);
			}
		}
	}
}

int network_island_count(void){
	int count = 0;	//島の数
	int i,j,dummy;
	int af_node_through_list_number = 0;
	int bf_node_through_list_number = 0;
	
	initialization_1(node_through_list,ag_num);	//初期化
	initialization_1(node_list_of_network,2);
	
	for(i=0;i<ag_num;i++){
		if(node_through_list[i] == 0){
			af_node_through_list_number = 0;
			
			dfs1(i);	//深さ優先探索関数
			
			for(j=0;j<ag_num;j++){
				if(node_through_list[j] == 1){
					af_node_through_list_number++;
				}
			}
			
			af_node_through_list_number = af_node_through_list_number - bf_node_through_list_number;
			bf_node_through_list_number = af_node_through_list_number + bf_node_through_list_number;
			
			dummy = af_node_through_list_number;
			
			if(dummy > node_list_of_network[1]){	//最も大きい島だけを記録する
				node_list_of_network[0] = count;	//何回目のdfs1には，
				node_list_of_network[1] = dummy;	//いくつのノードが集まっているか
			}
			
			count++;
		}
	}
	
	return count;
}

int culture_count_function(void){	//文化の種類を数える関数
	agent *kinds_of_culture_list;
	int i,j=0,k,m,count=0;
	
	for(i=0;i<ag_num;i++){				if(i==0){
			kinds_of_culture_list = (agent *)malloc(sizeof(agent)*(j+1));	//ここの時点ではj=0なのでj+1して配列を作る
			
			for(k=0;k<fe;k++){
				kinds_of_culture_list[j].f[k] = ag[i].f[k];
			}
			
			j++;
		}
		else{
			for(m=0;m<j;m++){
				if(KD_func(kinds_of_culture_list[m].f, ag[i].f) != fe){
					count++;	//このcountは2つのfor文の中にあることに注意する
				}
			}
			if(count==j){	//kinds_of_culture_listに含まれない文化ag[i]ならば
				if((kinds_of_culture_list = (agent *)realloc(kinds_of_culture_list,sizeof(agent)*(j+1))) == NULL){
					printf("realloc時にメモリが確保できません\n");
					free(kinds_of_culture_list);  /* 元のkinds_of_culture_listを解放して終了 */
					
					exit(1);
				}
				else{
					for(k=0;k<fe;k++){
						kinds_of_culture_list[j].f[k] = ag[i].f[k];
					}
				}
				j++;
			}
			
			count = 0;
		}
	}
	
	free(kinds_of_culture_list);
	
	return j;
}

void max_size_network_func(void){
	int count = 0;
	int i;
	
	initialization_1(node_through_list,ag_num);	//初期化
	
	for(i=0;i<ag_num;i++){	//行列の縦
		initialization_1(max_size_network[i],ag_num);	//行列の横
	}
	
	for(i=0;i<ag_num;i++){
		if(node_through_list[i] == 0){
			if(node_list_of_network[0] == count){
				dfs2(i);	//深さ優先探索関数(最も大きな連結グラフも同時作成)
			}
			else dfs1(i);	//深さ優先探索関数
			
			count++;
		}
	}
}

void orgnl_mdl(void){
	int i,j,x,y,yy;
	int n=0,m=0;
	int count=0;
	int *nb_lst;
	double a;
	int b;
	
	while(n==0){	//孤独なエージェント(隣接ノードが0)は除いて乱数を取る
		x = (int)(genrand_real1()*(double)(ag_num-1 -0 +1)/(double)ag_num);	//活性化させるagentをランダムに選出
	
		for(i=0;i<ag_num;i++){
			if(cd_nw[x][i] == 1){
				n = n+1;
			}
		}
	}
	
	nb_lst = (int *)malloc(sizeof(int) * n);	//nb_lst: neighbor agent listの略
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[x][i] == 1){
			nb_lst[m] = i;	//nb_lstはagent xの隣接ノードリスト
			m++;
		}
	}
	
	y = (int)genrand_int31() % n;	// 0割に注意
	
	if((0 < KD_func(ag[x].f, ag[ nb_lst[y] ].f)) && (KD_func(ag[x].f, ag[ nb_lst[y] ].f) < fe)){	//相互作用対象が全て自分と同じ文化を持っていた場合何もせずに終わる
		//共有していない文化をランダムに選ぶ
		b = (int)genrand_int31()%fe;
		
		while(ag[x].f[b] == ag[ nb_lst[y] ].f[b]){
			b = (int)genrand_int31()%fe;
		}
		
		//活性化agentのb番目の文化(ag[x].f[b])は他の隣接agentと共有しているか
		for(i=0;i<n;i++){
			if(ag[x].f[b] == ag[ nb_lst[i] ].f[b]){	//共有していた場合
				count++;	//これを使うと個々の保守性パラメータも考えることができる
			}
		}
		
		if(count > 0){	//保守性を含めた相互作用をする
			a = ((double)KD_func(ag[x].f, ag[ nb_lst[y] ].f)*(1-conservativeness_parameter))/(double)fe;	//確率( O(i,j)*保守性) /Fを計算
		}
		else if(count == 0){
			a = (double)KD_func(ag[x].f, ag[ nb_lst[y] ].f)/(double)fe;	//確率O(i,j)/Fを計算
		}
		
		if(a >= genrand_real1()){	//文化類似度による確率が成功した場合
			ag[x].f[b] = ag[ nb_lst[y] ].f[b];
		}
	}
	
	if(KD_func(ag[x].f, ag[ nb_lst[y] ].f) == 0){	//類似性=0のとき
		cd_nw[x][ nb_lst[y] ] = 0;
		cd_nw[ nb_lst[y] ][x] = 0;
		
		yy = (int)(genrand_real1()*(double)(ag_num-1 -0 +1)/(double)ag_num);
		
		while((x == yy) || ( nb_lst[y] == yy) || (cd_nw[x][yy] == 1) || (cd_nw[yy][x] == 1)){
			yy = (int)(genrand_real1()*(double)(ag_num-1 -0 +1)/(double)ag_num);
		}
		
		cd_nw[x][yy] = 1;
		cd_nw[yy][x] = 1;
	}
	
	free(nb_lst);
}

double a_rate_of_bridge_func(void){
	int a_rate_of_bridge=0;
	int all_link_count=0;
	int i,j;
	
	for(i=1;i<ag_num;i++){
		for(j=0;j<i;j++){
			if(cd_nw[i][j]==1){
				all_link_count++;
				
				if((0 < KD_func(ag[i].f,ag[j].f)) && (KD_func(ag[i].f,ag[j].f) < fe)){
					a_rate_of_bridge++;
				}
			}
		}
	}
	
	return a_rate_of_bridge/(double)all_link_count;
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

int dijkstra(int node_u,int node_v,int *ver,int *dis,int length){	//ダイクストラ法の本計算プログラム
	int i;
	int queue[ag_num];	//ある頂点に関する隣接ノードが格納される
	int head = 0, tail = 0;
	int node_adj;
	int dummy;
	
	while(vertex_list_check(ver,ag_num)){
		dummy = INT_MAX;
		
		for(i=0;i<length;i++){
			if(ver[i]==0){
				if(dummy >= dis[i]){
					dummy = dis[i];
				}
			}
		}
		
		if(dummy == INT_MAX){	//もしも最小値がINT_MAXであればこれ以上探索する意味は無い（->不可能到達点しかない）
			break;
		}
		
		for(i=0;i<length;i++){
			if(dummy == dis[i]){
				ver[i] = 1;
				queue[tail] = i;
				tail++;
			}
		}
		
		while(head != tail){
			node_u = queue[head];
			head++;
			
			for(i=0;i<length;i++){	//変数 i には隣接ノードが入る
				if(max_size_network[node_u][i]==1){
					if(dis[i] > (dis[node_u] + max_size_network[node_u][i])){
						dis[i] = dis[node_u] + max_size_network[node_u][i];
					}
				}
			}
		}
		//初期化
		initialization_1(queue,ag_num);
		head = 0;
		tail = 0;
	}
	
	return dis[node_v];
}

int shortest_path_func(int node_u, int node_v){	//頂点uから頂点vまでの最短頂点間距離を計算し，その距離を戻り値とする
	int vertex_list[ag_num];
	int distance_list[ag_num];
	int prev = node_u;
	int dummy;
	int i;
	
	//初期化
	for(i=0;i<ag_num;i++){
		//distance_listの初期化
		if(i==node_u){
			distance_list[i] = 0;
		}
		else{
			distance_list[i] = INT_MAX;
		}
		
		//vertex_listの初期化
		vertex_list[i] = 0;	//vertex_list[i]==1ならば頂点 i は頂点リストから除外されているものとする(0ならば頂点リストの要素)
	}
	//ここまで初期化
	
	//本計算
	return dijkstra(node_u, node_v, vertex_list, distance_list, ag_num);
}

double average_vertex_distance(void){	//平均頂点間距離を計算するプログラム(d(i,j)がINFの場合は数に入れない)
	int count=0;	//頂点間距離がINFの頂点の数を数える
	int sum=0;
	int d;
	int i,j;
	double av_dis;
	
	for(i=1;i<ag_num;i++){
		for(j=0;j<i;j++){
			d = shortest_path_func(i,j);
			
			if(INT_MAX == d){
				count++;
			}
			else{
				sum = d + sum;
			}
		}
	}
	
	av_dis = 2.0/(node_list_of_network[1]*(node_list_of_network[1]-1) ) * sum;
	
	return av_dis;
}

double cluster_calculate_func(int node_i){
	int j,k;
	int sum=0;
	int degree_num = degree_count(node_i);
	
	if(degree_num<2){	//次数が0か1であればクラスターは形成できないので0を返す
		return 0;
	}
	
	for(j=0;j<ag_num;j++){
		if(cd_nw[node_i][j]){
			for(k=0;k<ag_num;k++){
				if(cd_nw[node_i][k]){
					if(cd_nw[j][k]){
						sum++;
					}
				}
			}
		}
	}
	
	sum = sum / 2;	//同じ三角形を重複して数えているので2で割る
	
	return (2.0 / (degree_num * (degree_num - 1))) * sum;
}

int degree_count(int node_i){
	int i;
	int count=0;
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[node_i][i]){
			count++;
		}
	}
	
	return count;
}

double average_cluster_calculate_func(void){
	int node_i;
	double sum=0, average=0;
	int dummy;
	
	for(node_i=0;node_i<ag_num;node_i++){
		sum = cluster_calculate_func(node_i) + sum;
	}
	//"average = sum / ag_num"ではなぜかエラーが起きた
	dummy = ag_num;
	average = sum / dummy;
	
	return average;
}

//共同性に関するdata setを読み込む
void read_data_set_func(char *filename, double *cooperation_parameter, double *island_number_av_array, double *agent_number_av_array, double *culture_number_av_array, double *cluster_number_av_array, double *shortest_path_av_array, double *bridge_number_av_array, double *island_number_variance_array, double *agent_number_variance_array, double *culture_number_variance_array, double *cluster_number_variance_array, double *shortest_path_variance_array, double *bridge_number_variance_array, double *island_number_standard_deviation_array,  double *agent_number_standard_deviation_array, double *culture_number_standard_deviation_array, double *cluster_number_standard_deviation_array, double *shortest_path_standard_deviation_array, double *bridge_number_standard_deviation_array){
	FILE *fp;
	
	fp = fopen( filename, "r" );
	if( fp == NULL ){
		printf( "%sファイルが開けません\n", filename );
		exit(1);
	}
	else{
		while(fscanf( fp, "%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf",  cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array, island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array, island_number_standard_deviation_array,  agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array, bridge_number_standard_deviation_array) != EOF ){
			cooperation_parameter++;
			island_number_av_array++;
			agent_number_av_array++;
			culture_number_av_array++;
			cluster_number_av_array++;
			shortest_path_av_array++;
			bridge_number_av_array++;
			island_number_variance_array++;
			agent_number_variance_array++;
			culture_number_variance_array++;
			cluster_number_variance_array++;
			shortest_path_variance_array++;
			bridge_number_variance_array++;
			island_number_standard_deviation_array++;
			agent_number_standard_deviation_array++;
			culture_number_standard_deviation_array++;
			cluster_number_standard_deviation_array++;
			shortest_path_standard_deviation_array++;
			bridge_number_standard_deviation_array++;
		}
	}
	fclose( fp );
}

int cooperation_parameter_return_index(double *cooperation_parameter){
	int i;
	
	for(i=0;i<data_set_row;i++){
		if(cooperation_parameter[i] == conservativeness_parameter){
			return i;
		}
		else if(cooperation_parameter[i] > conservativeness_parameter){
			printf("data setに対する共同性値(conservativeness_parameter)が正しくないので近い値%fを返しました",cooperation_parameter[i]);
			
			return i;
		}
	}
	
	printf("強制終了します%d line)",__LINE__);
	exit(1);
}

int network_edges_count(void){
	int i,j;
	int count=0;
	
	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(cd_nw[i][j] == 1){
				count++;
			}
		}
	}
	
	return count;
}

/* main文 */
int main(void){
	FILE *fp,*fp2;
	char filename[100];
	char filename1[] = "cltfac_lst.txt";	//cultural factor list
	char filename2[100];
	char filename3[] = "max_size_network.gexf";
	char filename4[] = "max_size_network(sociarium).txt";
	char filename6[] = "cultural_network_animation(sociarium).txt";
	int i,j,k;
	int count;
	int t=0,edg=0;
	int simulation_count=0;
	int island_dummy;
	double distance_dummy;
	double cluster_dummy;
	
	//data setクラスの読み込み
	char filename5[] = "data_set_not_torus20151112.csv";
	char filename7[] = "data_set_torus20151214.csv";
	int index_mark;
	
	/* 共同性を格納する配列 */
	double cooperation_parameter[100];
	
	/* 平均を格納する配列 */
	double island_number_av_array[100];
	double agent_number_av_array[100];
	double culture_number_av_array[100];
	double cluster_number_av_array[100];
	double shortest_path_av_array[100];
	double bridge_number_av_array[100];
	
	/* 分散を格納する配列 */
	double island_number_variance_array[100];
	double agent_number_variance_array[100];
	double culture_number_variance_array[100];
	double cluster_number_variance_array[100];
	double shortest_path_variance_array[100];
	double bridge_number_variance_array[100];
	
	/* 標準偏差を格納する配列 */
	double island_number_standard_deviation_array[100];
	double agent_number_standard_deviation_array[100];
	double culture_number_standard_deviation_array[100];
	double cluster_number_standard_deviation_array[100];
	double shortest_path_standard_deviation_array[100];
	double bridge_number_standard_deviation_array[100];
	
	if(read_data_set){
		if(make_torus == 0){
			read_data_set_func(filename5, cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array,  island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array,  island_number_standard_deviation_array, agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array,  bridge_number_standard_deviation_array);
			index_mark = cooperation_parameter_return_index(cooperation_parameter);
		}
		else if(make_torus == 1){
			read_data_set_func(filename7, cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array,  island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array,  island_number_standard_deviation_array, agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array,  bridge_number_standard_deviation_array);
			index_mark = cooperation_parameter_return_index(cooperation_parameter);
		}
	}
	
	wt[0] = 0;	//書き出すtの時刻(1回目)
	wt[1] = 2500;	//書き出すtの時刻(2回目)
	wt[2] = 25000;	//書き出すtの時刻(3回目)
	wt[3] = 500000;	//書き出すtの時刻(4回目)
	
	init_genrand((unsigned)time(NULL));
	
	while(1){
		t=0;	//t の初期化
		lattice_func();
		tr_func();
		olap_func();
		latcd_cpyfunc();	//レギュラー格子の隣接行列を文化流布して変化するのに使う隣接行列に要素をコピーする
		
		/* ファイルオープン */
		if ((fp2 = fopen(filename6, "wt")) == NULL) {
			fprintf(stderr, "ファイルのオープンに失敗しました．\n");
			return exit(1);
		}
		
		while(t<=tmax){
			for(k=0;k<wt_num;k++){
				if(wt[k] == t){
					for(i=0;i<ag_num;i++){	//書き込む前にolapを正しく計算する
						for(j=0;j<ag_num;j++){
							if (cd_nw[i][j]==1){
								olap[i][j] = KD_func(ag[i].f,ag[j].f);
							}
						}
					}
					
					if(t==0){	//t=0(格子状)のネットワークを書き出す
						edg=0;	//edgの初期化
						
						/* gephi faile生成 */
						sprintf(filename,"nw_%d.gexf",wt[k]); /*ファイル名整形*/
						
						/* ファイルオープン */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "ファイルのオープンに失敗しました．\n");
							return exit(1);
						}
						
						/* 書き込み */
						fprintf(fp,"<gexf \n");
						fprintf(fp,"xmlns:ns0=\"http://www.gexf.net/1.1draft/viz\" \n");
						fprintf(fp,"version=\"1.1\" \n");
						fprintf(fp,"xmlns=\"http://www.gexf.net/1.1draft\" \n");
						fprintf(fp,"xmlns:viz=\"http://www.gexf.net/1.1draft/viz\" \n");
						fprintf(fp,"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n");
						fprintf(fp,"xsi:schemaLocation=\"http://www.w3.org/2001/XMLSchema-instance\"> \n");
						
						fprintf(fp,"<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">\n");
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<ag_num;i++){
							fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
							fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",ag[i].f[0]*rgb_conversion,ag[i].f[1]*rgb_conversion,ag[i].f[2]*rgb_conversion);
							fprintf(fp,"</node>\n");
						}
						
						fprintf(fp,"</nodes>\n");
						fprintf(fp,"<edges>\n");
						
						for(i=0;i<ag_num;i++){
							for(j=i;j<ag_num;j++){
								if(cd_nw[i][j]==1){
									fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)olap[i][j]);
									edg = edg+1;
								}
							}
						}
						
						fprintf(fp,"</edges>\n");
						fprintf(fp,"</graph>\n");
						fprintf(fp,"</gexf>\n");
						
						/* ファイルクローズ */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file生成 */
						sprintf(filename2,"nw_%d(sociarium).txt",wt[k]); /*ファイル名整形*/
						
						/* ファイルオープン */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "ファイルのオープンに失敗しました．\n");
							return exit(1);
						}
						
						/* 書き込み */
						fprintf(fp,"# Pajek形式のデータ \n");
						fprintf(fp,"# 参照: http://www.tp.umu.se/~rosvall/code.html \n");
						fprintf(fp,"\n");
						fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
						fprintf(fp,"@title = The dissemination of culture \n");
						fprintf(fp,"@nondirected \n");
						fprintf(fp,"\n");
						
						fprintf(fp,"*Vertices %d \n",ag_num);
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<ag_num;i++){
							fprintf(fp,"%d \"agent_%d\" \n",i,i);
						}
						
						fprintf(fp,"*Arcs \n");
						
						for(i=0;i<ag_num;i++){
							for(j=i;j<ag_num;j++){
								if(cd_nw[i][j]==1){
									fprintf(fp,"%d %d %0.1f \n",i,j,(double)olap[i][j]);
								}
							}
						}
						
						/* ファイルクローズ */
						fclose(fp);
						/*----------------------------------------------------*/
					}
					else{	//時刻tのネットワークを書き出す
						edg=0;	//edgの初期化
						/* gephi file作成 */
						sprintf(filename,"nw_%d.gexf",wt[k]); /*ファイル名整形*/
						
						/* ファイルオープン */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "ファイルのオープンに失敗しました．\n");
							return exit(1);
						}
						
						/* 書き込み */
						fprintf(fp,"<gexf \n");
						fprintf(fp,"xmlns:ns0=\"http://www.gexf.net/1.1draft/viz\" \n");
						fprintf(fp,"version=\"1.1\" \n");
						fprintf(fp,"xmlns=\"http://www.gexf.net/1.1draft\" \n");
						fprintf(fp,"xmlns:viz=\"http://www.gexf.net/1.1draft/viz\" \n");
						fprintf(fp,"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n");
						fprintf(fp,"xsi:schemaLocation=\"http://www.w3.org/2001/XMLSchema-instance\"> \n");
						
						fprintf(fp,"<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">\n");
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<ag_num;i++){
							fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
							fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",ag[i].f[0]*rgb_conversion,ag[i].f[1]*rgb_conversion,ag[i].f[2]*rgb_conversion);
							fprintf(fp,"</node>\n");
						}
						
						fprintf(fp,"</nodes>\n");
						fprintf(fp,"<edges>\n");
						
						for(i=0;i<ag_num;i++){
							for(j=i;j<ag_num;j++){
								if(cd_nw[i][j]==1){
									fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)olap[i][j]);
									edg = edg+1;
								}
							}
						}
						
						fprintf(fp,"</edges>\n");
						fprintf(fp,"</graph>\n");
						fprintf(fp,"</gexf>\n");
						
						/* ファイルクローズ */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file生成 */
						sprintf(filename2,"nw_%d(sociarium).txt",wt[k]); /*ファイル名整形*/
						
						/* ファイルオープン */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "ファイルのオープンに失敗しました．\n");
							return exit(1);
						}
						
						/* 書き込み */
						fprintf(fp,"# Pajek形式のデータ \n");
						fprintf(fp,"# 参照: http://www.tp.umu.se/~rosvall/code.html \n");
						fprintf(fp,"\n");
						fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
						fprintf(fp,"@title = The dissemination of culture \n");
						fprintf(fp,"@nondirected \n");
						fprintf(fp,"\n");
						
						fprintf(fp,"*Vertices %d \n",ag_num);
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<ag_num;i++){
							fprintf(fp,"%d \"agent_%d\" \n",i,i);
							//fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",ag[i].f[0]*rgb_conversion,ag[i].f[1]*rgb_conversion,ag[i].f[2]*rgb_conversion);
						}
						
						fprintf(fp,"*Arcs \n");
						
						for(i=0;i<ag_num;i++){
							for(j=i;j<ag_num;j++){
								if(cd_nw[i][j]==1){
									fprintf(fp,"%d %d %0.1f \n",i,j,(double)olap[i][j]);
								}
							}
						}
						
						/* ファイルクローズ */
						fclose(fp);
						/*----------------------------------------------------*/
					}
				}
			}
			
			/* sociarium file生成 */	
			if(t%sociarium_ptch_width==0){
				if(t==0){
					/* 書き込み */
					fprintf(fp2,"@module = graph_creation_read_time_series_rect.dll \n");
					fprintf(fp2,"@title = cultural drift %f\n",conservativeness_parameter);
					fprintf(fp2,"@delimiter =\t\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@time_format = Y.M.D # Y[ear], M[onth], D[ay], h[our], m[inute], s[econd] \n");
					fprintf(fp2,"@interval  = %d # ネットワークを作成する時刻の間隔\n",sociarium_ptch_width);
					fprintf(fp2,"@characteristic_time = 1 # 各時刻についてexp(-dt/@characteristic_time)を積算し，@threshold以下のエッジをマスキング \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@threshold = 0\n");
					fprintf(fp2,"@nondirected \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@node_texture = agent.png\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", 0, 0, 0, 0);	//これがないと正しい描画ができない
					
					for(i=0;i<ag_num;i++){
						for(j=i;j<ag_num;j++){
							if(cd_nw[i][j]){
								fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", t+sociarium_ptch_width, i, j, cd_nw[i][j]);
							}
						}
					}
				}
				else{
					for(i=0;i<ag_num;i++){
						for(j=i;j<ag_num;j++){
							if(cd_nw[i][j]){
								fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", t+sociarium_ptch_width, i, j, cd_nw[i][j]);
							}
						}
					}
				}
			}			
			
			/***********プログラムの核(orgnl_mdl)となるもの***********
			---------------------以下具体的な内容-----------------------
			1.活性化agentがランダムに選択される
			2.そのagentの隣接agentをランダムに1つ選択する
			3.類似性を計算する
			3.相互作用を行う
				-隣接agentから確率O(i,j)/Fで文化の影響を受ける
			4.overlapが0ならばrule5に則る
			***********************************************************/
			orgnl_mdl();
			
			t = t+1;
		}
		/* ファイルクローズ */
		fclose(fp2);
		
		max_size_network_func();	//最大の島の隣接行列を生成
		
		island_dummy = network_island_count();
		distance_dummy = average_vertex_distance();
		cluster_dummy = average_cluster_calculate_func();
		
		if(read_data_set){
			if(simulation_count >= SimulationStop){
				printf("シミュレーションを強制終了します(Error : %d line)\n",__LINE__);
				exit(1);
			}			
			//該当のparameterを探す
			
		/* ログファイル生成 */
		/* ファイルオープン */
			if ((fp = fopen("simulation_log.csv", "a")) == NULL) {
				fprintf(stderr, "ファイルのオープンに失敗しました．\n");
				exit(1);
			}
			
			/* 書き込み */
			fprintf(fp,"%d,%f,%f \n",island_dummy,distance_dummy,cluster_dummy);
			
			fclose(fp);
			
			if((island_number_av_array[index_mark] - sigma_option*island_number_standard_deviation_array[index_mark] < island_dummy)&&(island_dummy< island_number_av_array[index_mark] + sigma_option*island_number_standard_deviation_array[index_mark])&&(shortest_path_av_array[index_mark] - sigma_option*shortest_path_standard_deviation_array[index_mark] < distance_dummy)&&(distance_dummy < shortest_path_av_array[index_mark] + sigma_option*shortest_path_standard_deviation_array[index_mark])&&(cluster_number_av_array[index_mark] - sigma_option*cluster_number_standard_deviation_array[index_mark] < cluster_dummy)&&(cluster_dummy < cluster_number_av_array[index_mark] + sigma_option*cluster_number_standard_deviation_array[index_mark])){
				break;
			}
			simulation_count++;
			
			printf("再シミュレーションしています… %d / %d\n",simulation_count,SimulationStop);
		}
		else{
			break;
		}
	}
	
	/* ファイルオープン */
	edg = 0;
	if ((fp = fopen(filename3, "w")) == NULL) {
		fprintf(stderr, "ファイルのオープンに失敗しました．\n");
		return exit(1);
	}
	
	/////////////////////////島のネットワークのうち最も大きなネットワークを描画する/////////////////////////
	//Gephi_file
	/* 書き込み */
	fprintf(fp,"<gexf \n");
	fprintf(fp,"xmlns:ns0=\"http://www.gexf.net/1.1draft/viz\" \n");
	fprintf(fp,"version=\"1.1\" \n");
	fprintf(fp,"xmlns=\"http://www.gexf.net/1.1draft\" \n");
	fprintf(fp,"xmlns:viz=\"http://www.gexf.net/1.1draft/viz\" \n");
	fprintf(fp,"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \n");
	fprintf(fp,"xsi:schemaLocation=\"http://www.w3.org/2001/XMLSchema-instance\"> \n");
	
	fprintf(fp,"<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">\n");
	
	fprintf(fp,"<nodes>\n");
	
	for(i=0;i<ag_num;i++){
		for(j=0;j<ag_num;j++){
			if(max_size_network[i][j] == 1){
				fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
				fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",ag[i].f[0]*rgb_conversion,ag[i].f[1]*rgb_conversion,ag[i].f[2]*rgb_conversion);
				fprintf(fp,"</node>\n");
				
				break;
			}
		}
	}

	fprintf(fp,"</nodes>\n");
	fprintf(fp,"<edges>\n");

	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(max_size_network[i][j]==1){
				fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)olap[i][j]);
				edg = edg+1;
			}
		}
	}

	fprintf(fp,"</edges>\n");
	fprintf(fp,"</graph>\n");
	fprintf(fp,"</gexf>\n");
	
	/* ファイルクローズ */
	fclose(fp);
	
	/* sociarium file生成 */	
	/* ファイルオープン */
	if ((fp = fopen(filename4, "wt")) == NULL) {
		fprintf(stderr, "ファイルのオープンに失敗しました．\n");
		return exit(1);
	}
	
	/* 書き込み */
	fprintf(fp,"# Pajek形式のデータ \n");
	fprintf(fp,"# 参照: http://www.tp.umu.se/~rosvall/code.html \n");
	fprintf(fp,"\n");
	fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
	fprintf(fp,"@title = The dissemination of culture (max size network)\n");
	fprintf(fp,"@nondirected \n");
	fprintf(fp,"\n");
	
	fprintf(fp,"*Vertices %d \n",node_list_of_network[1]);
	
	fprintf(fp,"<nodes>\n");
	
	for(i=0;i<ag_num;i++){
		fprintf(fp,"%d \"agent_%d\" \n",i,i);
	}
	
	fprintf(fp,"*Arcs \n");
	
	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(max_size_network[i][j]==1){
				fprintf(fp,"%d %d %0.1f \n",i,j,(double)olap[i][j]);
			}
		}
	}
	
	/* ファイルクローズ */
	fclose(fp);
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//文化の流布が成功したのか文化リストを作成して確認する
	/* ファイルオープン */
	if ((fp = fopen(filename1, "w")) == NULL) {
		fprintf(stderr, "ファイルのオープンに失敗しました．\n");
		return exit(1);
	}
	
	/* 書き込み */
	fprintf(fp,"島の数 = %3d \n",island_dummy);
	fprintf(fp,"最も大きい島のagent数 = %3d \n",node_list_of_network[1]);
	fprintf(fp,"文化の数 = %3d \n",culture_count_function());
	fprintf(fp,"平均頂点間距離 = %3f \n",distance_dummy);
	fprintf(fp,"平均クラスター係数 = %3f \n",cluster_dummy);
	fprintf(fp,"ブリッジ数(%%) = %f \n",a_rate_of_bridge_func()*100);
	
	fclose(fp);
	
	return 0;
}
