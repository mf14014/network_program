/* プログラミングメモ：
・original_modelAVcal1_3からの変更点
	【拡張】格子モデルをトーラスにする関数の作成
*/


/* ヘッダーファイル */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "MT.h"


/* プリプロセッサ */
#define lg 20	//格子の縦または横の長さ．
#define fe 3	//features．特徴の数
#define tr 20	//traits．特性の範囲
#define ag_num lg*lg	//agent数
#define tmax 500000	//イベント数

//pitch_width 0.1，desired_number 5で計算するとおよそ1分掛かる
#define pitch_width 0.01	//保守性パラメータの増加速度(= 保守性0から1までの刻み幅)
#define pitch_length (int)(1 / pitch_width)
#define desired_number 1000	//1つの保守性パラメータあたりのプログラム実行回数

#define make_torus 0	//格子モデルをトーラスにする場合は1にする(しない場合は0)

#define countIslandStepSize 10000	//一応1000もやる
#define event_length tmax/countIslandStepSize

/* 構造体 */
typedef struct{
	int f[fe];	//特徴
	
} agent;


/* 外部変数 */
int lattice_mdl[ag_num][ag_num];	//格子モデルの隣接行列，agはagentのこと
int cd_nw[ag_num][ag_num];	//cultural driftするネットワークの隣接行列
int olap[ag_num][ag_num];	//overlap
agent ag[ag_num];
int node_through_list[ag_num];	//dfsで通過したノードは該当する列に1を入れる
int node_through_count;	//dfsで通過したノードをカウントする(dfsは再帰呼び出ししているのでこの変数は外部変数にしなきゃダメ)
int node_list_of_network[2];	//添字0の配列には「何個目の連結グラフか」を，添字1には「すべての連結グラフの中で最も大きいグラフのノード数」を格納
int max_size_network[ag_num][ag_num];	//連結グラフのうち最も多くのagent数をもつグラフを保存する
double conservativeness_parameter;	//保守性パラメータ(0以上1以下で数字を選ぶ) ->0(保守性を持たない状態)を選ぶとCentolaモデルになる


/* プロトタイプ宣言 */
void lattice_func(void);
void lattice_to_torus(void);
void tr_func(void);
int KD_func(int *p, int*q);
void olap_func(void);
void latcd_cpyfunc(void);
void initialize_int_array(int *p,int x);
void dfs1(int v);
void dfs2(int v);
int network_island_count(void);
int culture_count_function(void);
void max_size_network_func(void);
void orgnl_mdl(void);
double cluster_calculate_func(int node_i);
int degree_count(int node_i);
double average_cluster_calculate_func(void);
int vertex_list_check(int *p, int array_length);	//配列の中身がすべて1であれば0を返し，それ以外であれば1を返す
int dijkstra(int node_u,int node_v,int *ver,int *dis,int length);	//ダイクストラ法の本計算プログラム
int shortest_path_func(int node_u, int node_v);	//最短頂点間距離を計算する関数
double average_vertex_distance(void);	//平均頂点間距離を求める関数
double a_rate_of_bridge_func(void);	//ブリッジなリンクの数を数える関数


/* 外部関数 */
void lattice_func(void){	//格子モデル作成アルゴリズム
	int i,j;
	
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
		
		for(j=0;j<fe;j++){
			ag[i].f[j] = (int)genrand_int31() % tr;	//[0,tr-1]までの間で乱数を取る
			/*
			while(ag[i].f[j] < 0){
				ag[i].f[j] = (int)genrand_real2()*10000 % tr;
			}*/
			
		}
	}
}

int KD_func(int *p, int *q){	//2つの配列pとqが互いに共有している文化の数を返す関数
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
void initialize_int_array(int *p,int x){	//pは1次元配列のポインタ，xは配列の大きさ(2次元にも応用可 -> *pに配列の行を渡せば良い)
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

//与えられた配列の中身を全て0にする
void initialize_double_array(double *p,int x){	//pは1次元配列のポインタ，xは配列の大きさ(2次元にも応用可 -> *pに配列の行を渡せば良い)
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
		if(cd_nw[v][i] == 1 && node_through_list[i] == 0){
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
			//max_size_network[i][v] = 1;	//多分なくても正常動作する
			
			if(node_through_list[i] == 0){
				dfs2(i);
			}
		}
	}
}

int network_island_count(void){	//島の数を数える関数
	int count = 0;	//島の数
	int i,j,dummy;
	int af_node_through_list_number = 0;
	int bf_node_through_list_number = 0;
	
	initialize_int_array(node_through_list,ag_num);	//初期化
	
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
			
			//printf("%d\n",dfs(i));
			
			if(dummy > node_list_of_network[1]){	//最も大きい島だけを記録する
				node_list_of_network[0] = count;	//何回目のdfs1には，
				node_list_of_network[1] = dummy;	//いくつのノードが集まっているか
				
				//printf("test1\n");
			}
			//printf("test2\n");
			
			count++;
		}
	}
	
	return count;
}

int culture_count_function(void){	//文化の種類を数える関数
	agent *kinds_of_culture_list;
	int i,j=0,k,m,count=0;
	
	for(i=0;i<ag_num;i++){		
		if(i==0){
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
	//printf("test_end\n");	//ここまで行き着かない
	
	free(kinds_of_culture_list);
	
	return j;
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
	//sleep(1);
	//printf("average = %f\n",average);
	return average;
}

void max_size_network_func(void){
	int count = 0;
	int i;
	
	initialize_int_array(node_through_list,ag_num);	//初期化
	
	for(i=0;i<ag_num;i++){
		if(node_through_list[i] == 0){
			if(node_list_of_network[0] == count){
				//printf("node_list_of_network[0] = %d\n",node_list_of_network[0]);
				//printf("test3\n");
				dfs2(i);	//深さ優先探索関数(最も大きな連結グラフも同時作成)
			}
			else dfs1(i);	//深さ優先探索関数
			
			count++;
		}
	}
}

double variance_func(int *p,double av){	//分散を計算する関数．avはaverage
	int i;
	double s=0;
	
	for(i=0;i<desired_number;i++){
		s = s + pow(av-*p,2);
		p++;
	}
	
	s = s / desired_number;
	
	return s;
}

double variance_func2(double *p,double av){	//分散を計算する関数．avはaverage
	int i;
	double s=0;
	
	for(i=0;i<desired_number;i++){
		s = s + pow(av-*p,2);
		p++;
	}
	
	s = s / desired_number;
	
	return s;
}

void orgnl_mdl(void){	//original model（修正モデルのこと）
	int i,x,y,yy;
	int n=0,m=0;
	int count=0;
	int *nb_lst;
	double a;
	int b;
	
	//printf("debag_lineA\n");
	
	while(n==0){	//孤独なエージェント(隣接ノードが0)は除いて乱数を取る
		x = (int)(genrand_real1()*(double)(ag_num-1 -0 +1)/(double)ag_num);	//活性化させるagentをランダムに選出
	
		for(i=0;i<ag_num;i++){
			if(cd_nw[x][i] == 1){
				//printf("test by Centola_AA\n");
				n = n+1;
			}
		}
	}
	
	//printf("debag_lineB\n");
	
	nb_lst = (int *)malloc(sizeof(int) * n);	//nb_lst: neighbor agent listの略
	
	//printf("debag_lineC\n");
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[x][i] == 1){
			nb_lst[m] = i;	//nb_lstはagent xの隣接ノードリスト
			
			m++;
		}
	}
	
	//printf("debag_lineD\n");
	
	y = (int)genrand_int31() % n;	// 0割に注意
	
	//printf("debag_lineE\n");
	
	if((0 < KD_func(ag[x].f, ag[ nb_lst[y] ].f)) && (KD_func(ag[x].f, ag[ nb_lst[y] ].f) < fe)){	//相互作用対象が全て自分と同じ文化を持っていた場合何もせずに終わる
		//printf("test by olgmdl_mdl_A\n");
		
		
		//printf("%d\n",x);
		
		
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
		
		//printf("%f/%f=%f\n",(double)KD_func(ag[x].f, ag[ nb_lst[y] ].f),(double)fe,a);
		//printf("test by Centola_mdl_B\n");
		
		if(a >= genrand_real1()){	//文化類似度による確率が成功した場合
			//printf("test by Centola_mdl_C\n")
			
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
		initialize_int_array(queue,ag_num);
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
	//printf("count->%d, sum->%d\n",count,sum);
	
	av_dis = 2.0/(node_list_of_network[1]*(node_list_of_network[1]-1) ) * sum;
	
	return av_dis;
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
	//printf("link -> %d\n",all_link_count);
	
	return a_rate_of_bridge;
}

void island_variation_output(double islandNumberList[pitch_length][event_length]){
	FILE *countIslandListCsvFilePointer;
	char countIslandListCsv[] = "count_island_list.csv";
	char communalityListString[1000] = ",";
	char dummyStringVariable[100];
	int rowSize = tmax / countIslandStepSize;
	int islandNumberListIndex;
	double communality;
	
	for (communality = 0 ; communality < 1 - pitch_width ; communality += pitch_width) {
		sprintf(dummyStringVariable, "%f, ", communality);	//double型をstring型に変換
		strcat(communalityListString, dummyStringVariable);
	}
	
	if ((countIslandListCsvFilePointer = fopen(countIslandListCsv, "w")) == NULL) {
		fprintf(stderr, "open file failed");
		exit(1);
	} else {
		strcat(communalityListString, "\n");
		fprintf(countIslandListCsvFilePointer, communalityListString);	//列ラベルを書き込む
		
		for (islandNumberListIndex = 0 ; islandNumberListIndex < rowSize ; islandNumberListIndex++) {
			memset(communalityListString, '\0' ,strlen(communalityListString));	//char型の変数を初期化
			sprintf(dummyStringVariable, "%d, ", (islandNumberListIndex+1)*countIslandStepSize);
			strcat(communalityListString, dummyStringVariable);
			
			for (communality = 0 ; communality < 1 - pitch_width ; communality += pitch_width) {
				sprintf(dummyStringVariable, "%f, ", (islandNumberList[(int)(communality*pitch_length)][islandNumberListIndex])/desired_number);	//書き込む際に平均処理も行う
				strcat(communalityListString, dummyStringVariable);
			}
			strcat(communalityListString, "\n");
			fprintf(countIslandListCsvFilePointer, communalityListString);
		}
	}
	
	fclose(countIslandListCsvFilePointer);
}

void count_island_input_list(double islandNumberList[pitch_length][event_length], double  present_communality, int eventTime){
	int islandNumberListIndexRow = (int)(present_communality*pitch_length+0.000001);	//キャストした値が変化することがあったので丸めを利用して回避
	
	if (! tmax / countIslandStepSize) {
		printf("[ERROR] countIslandStepSizeがtmaxを超えて設定されています\n");
		exit(1);
	}
	
	if (! (eventTime % countIslandStepSize) && eventTime != 0) {
		//printf("present_communality*pitch_length = %f -floor-> %f -> %d\n", islandNumberListIndexRow, floor(islandNumberListIndexRow), (int)islandNumberListIndexRow);
		//printf("eventTime ％ countIslandStepSize = %d\n", eventTime % countIslandStepSize);
		islandNumberList[islandNumberListIndexRow][(eventTime / countIslandStepSize) - 1] += (double)(network_island_count());
		//printf("islandNumberList[%d][%d]  = %f\n", islandNumberListIndexRow, (eventTime / countIslandStepSize) - 1 , islandNumberList[islandNumberListIndexRow][(eventTime / countIslandStepSize) - 1] );
	}
}

/* main文 */
int main(void){
	FILE *fp;
	char filename[] = "cltural_factor_average_lst.csv";	//cultural factor list
	double i;
	int j,k;
	int t;
	
	/* 平均 */
	double island_number_av;
	double agent_number_av;
	double culture_number_av;
	double cluster_number_av;
	double shortest_path_av;
	double bridge_number_av;
	
	/* 分散 */
	double island_number_variance;
	double agent_number_variance;
	double culture_number_variance;
	double cluster_number_variance;
	double shortest_path_variance;
	double bridge_number_variance;
	
	int island_number_array[desired_number];
	int agent_number_array[desired_number];
	int culture_number_array[desired_number];
	double cluster_number_array[desired_number];
	double shortest_path_array[desired_number];
	double bridge_number_array[desired_number];
	
	double island_number_list[pitch_length][event_length];
	for(k=0;k<pitch_length;k++){	//行列の縦
		initialize_double_array(island_number_list[k], event_length);
	}
	
	/* ファイルオープン */
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "ファイルのオープンに失敗しました．\n");
		exit(1);
	}
	
	init_genrand((unsigned)time(NULL));
	
	for(i=0;i<=1.0-pitch_width;i+=pitch_width){
		//忘れずに初期化すること
		island_number_av = 0;
		agent_number_av = 0;
		culture_number_av = 0;
		cluster_number_av = 0;
		shortest_path_av = 0;
		bridge_number_av = 0;
		
		for(j=0;j<desired_number;j++){
			//外部変数の初期化
			for(k=0;k<ag_num;k++){	//行列の縦
				initialize_int_array(lattice_mdl[k],ag_num);	//行列の横
			}
			
			for(k=0;k<ag_num;k++){	//行列の縦
				initialize_int_array(cd_nw[k],ag_num);	//行列の横
			}
			
			initialize_int_array(node_through_list,ag_num);
			
			node_through_count = 0;
			
			initialize_int_array(node_list_of_network,2);
			
			for(k=0;k<ag_num;k++){	//行列の縦
				initialize_int_array(max_size_network[k],ag_num);	//行列の横
			}
			
			lattice_func();
			tr_func();
			olap_func();
			latcd_cpyfunc();	//レギュラー格子の隣接行列を文化流布して変化するのに使う隣接行列に要素をコピーする
			
			conservativeness_parameter = i;
			
			for(t=0;t<=tmax;t++){
				//printf("%d event\n",t);
				
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
				count_island_input_list(island_number_list, i, t);
			}
			
			//printf("%d\n",network_island_count());
			
			island_number_array[j] = network_island_count();
			max_size_network_func();
			
			agent_number_array[j] = node_list_of_network[1];
			culture_number_array[j] = culture_count_function();
			cluster_number_array[j] = average_cluster_calculate_func();
			shortest_path_array[j] = average_vertex_distance();
			bridge_number_array[j] = a_rate_of_bridge_func();
			
			/*printf("island_number_array[%d]->%d\n",j,island_number_array[j]);
			printf("agent_number_array[%d]->%d\n",j,agent_number_array[j]);
			printf("culture_number_array[%d]->%d\n",j,culture_number_array[j]);
			printf("cluster_number_array[%d]->%f\n",j,cluster_number_array[j]);
			printf("shortest_path_array[%d]->%f\n",j,shortest_path_array[j]);
			printf("bridge_number_array[%d]->%f\n",j,bridge_number_array[j]);*/
			
			island_number_av = (double)island_number_array[j] + island_number_av;
			agent_number_av = (double)agent_number_array[j] + agent_number_av;
			culture_number_av = (double)culture_number_array[j] + culture_number_av;
			cluster_number_av = cluster_number_array[j] + cluster_number_av;
			shortest_path_av = shortest_path_array[j] + shortest_path_av;
			bridge_number_av = bridge_number_array[j] + bridge_number_av;
		}
		
		island_number_av = island_number_av / desired_number;
		agent_number_av = agent_number_av / desired_number;
		culture_number_av = culture_number_av / desired_number;
		cluster_number_av = cluster_number_av / desired_number;
		shortest_path_av = shortest_path_av / desired_number;
		bridge_number_av = bridge_number_av / desired_number;
		
		island_number_variance = variance_func(island_number_array,island_number_av);
		agent_number_variance = variance_func(agent_number_array,agent_number_av);
		culture_number_variance = variance_func(culture_number_array,culture_number_av);
		cluster_number_variance = variance_func2(cluster_number_array,cluster_number_av);
		shortest_path_variance = variance_func2(shortest_path_array,shortest_path_av);
		bridge_number_variance = variance_func2(bridge_number_array,bridge_number_av);
		
		//文化の流布が成功したのか文化リストを作成して確認する	
		/* 書き込み */
		if(i==0){
			fprintf(fp,",島の数(平均) ,最も大きい島のagent数(平均) , 文化の数(平均), 平均クラスター係数(平均), 平均頂点間距離(平均) ,ブリッジ数(平均) ,島の数(分散) ,最も大きい島のagent数(分散) ,文化の数(分散) ,平均クラスター係数(分散), 平均頂点間距離(分散) ,ブリッジ数(分散) ,島の数(標準偏差) ,最も大きい島のagent数(標準偏差) ,文化の数(標準偏差),平均クラスター係数(標準偏差), 平均頂点間距離(標準偏差), ブリッジ数(標準偏差)\n");
			fprintf(fp,"%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f\n",i ,island_number_av ,agent_number_av ,culture_number_av, cluster_number_av, shortest_path_av, bridge_number_av, island_number_variance ,agent_number_variance , culture_number_variance, cluster_number_variance, shortest_path_variance, bridge_number_variance, sqrt(island_number_variance) ,sqrt(agent_number_variance), sqrt(culture_number_variance), sqrt(cluster_number_variance), sqrt(shortest_path_variance), sqrt(bridge_number_variance));
		}
		else{
			fprintf(fp,"%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f ,%5f\n",i ,island_number_av ,agent_number_av ,culture_number_av, cluster_number_av, shortest_path_av, bridge_number_av, island_number_variance ,agent_number_variance , culture_number_variance, cluster_number_variance, shortest_path_variance, bridge_number_variance, sqrt(island_number_variance) ,sqrt(agent_number_variance), sqrt(culture_number_variance), sqrt(cluster_number_variance), sqrt(shortest_path_variance), sqrt(bridge_number_variance));
		}
		
		printf("%f \n",i);
	}
	
	fclose(fp);
	
	island_variation_output(island_number_list);
	
	return 0;
}
