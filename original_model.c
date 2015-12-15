/* �v���O���~���O�����F
original_model-6a����̕ύX�_
	sociarium�t�@�C���i�A�j���[�V�����j�𐶐�����v���O�����̎���
*/


/* �w�b�_�[�t�@�C�� */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "MT.h"

/* �v���v���Z�b�T */
#define lg 20	//�i�q�̏c�܂��͉��̒����D
#define fe 3	//features�D�����̐�
#define tr 20	//traits�D�����͈̔�
#define ag_num lg*lg	//agent��;
#define tmax 500000	//�C�x���g��
#define wt_num 4	//�����o���t�@�C���̐�
#define rgb_conversion 255/tr	//agent��rgb�Œl��������
#define make_torus 0	//�i�q���f�����g�[���X�ɂ���ꍇ��1�ɂ���(���Ȃ��ꍇ��0)
#define conservativeness_parameter 0.99	//�ێ琫�p�����[�^(0�ȏ�1�ȉ��Ő�����I��) ->0(�ێ琫�������Ȃ����)��I�Ԃ�Centola���f���ɂȂ�
#define read_data_set 1	//data_set��ǂݍ���data set�̂��鋤�����ɑ΂��鑮���l(colum)�̕��ςɋ߂��O���t�𐶐�����i�ǂݍ��܂Ȃ��ꍇ��0�j
	#define data_set_row 100	//data set�̋�������0����1�܂ł̍��݂̌�
	#define sigma_option 1	//�O���t�����ɂ���������͈� �W���΍�(sigma)��1->68.27%�C2->95.45%�C3->99.73%
	#define SimulationStop 100	//�V�~�����[�V�����񐔁D���̉񐔂𒴂����ꍇ�����I�����܂��D
#define sociarium_ptch_width 5000	//sociarium�̃O���t���A�j���[�V�����ɂ���ہC���Ԏ������P�ʂŐi�߂Ă���������

/* �\���� */
typedef struct{
	int f[fe];	//����
	
} agent;


/* �O���ϐ� */
int lattice_mdl[ag_num][ag_num];	//�i�q���f���̗אڍs��
int cd_nw[ag_num][ag_num];	//cultural drift����l�b�g���[�N�̗אڍs��
int olap[ag_num][ag_num];	;//overlap
int wt[wt_num];	//���g��main���ɏ����Ă���
agent ag[ag_num];
int node_through_list[ag_num];	//dfs�Œʉ߂����m�[�h�͊Y��������1������
int node_through_count;	//dfs�Œʉ߂����m�[�h���J�E���g����(dfs�͍ċA�Ăяo�����Ă���̂ł��̕ϐ��͊O���ϐ��ɂ��Ȃ���_��)
int node_list_of_network[2];	//�Y��0�̔z��ɂ́u���ڂ̘A���O���t���v���C�Y��1�ɂ́u���ׂĂ̘A���O���t�̒��ōł��傫���O���t�̃m�[�h���v���i�[
int max_size_network[ag_num][ag_num];	//�A���O���t�̂����ł�������agent�������O���t��ۑ�����


/* �v���g�^�C�v�錾 */
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


/* �O���֐� */
void lattice_func(void){	//�i�q���f���쐬�A���S���Y��
	int i,j;
	
	/* ������ */
	for(i=0;i<ag_num;i++){	//�s��̏c
		initialization_1(lattice_mdl[i],ag_num);	//�s��̉�
		initialization_1(cd_nw[i],ag_num);	//�s��̉�
	}
	
	for(i=0;i<ag_num;i++){
		for(j=i;j<ag_num;j++){
			if(i == j){
				lattice_mdl[i][j] = 0;
			}
			else{
				if( ((i+1) == j) && ((i+1)%lg != 0)){
					lattice_mdl[i][j] = 1;	//�E�Ƀ����N��L�΂�
					lattice_mdl[j][i] = 1;
				}
				if( ((j-i) == lg) && ((i+lg) < ag_num)){
					lattice_mdl[i][i+lg] = 1;	//���Ƀ����N��L�΂�
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
					lattice_mdl[i][j] = 1;	//���Ƀ����N��L�΂�
					lattice_mdl[j][i] = 1;
				}
				if( ((j-i) == (lg*lg-lg)) && (i < lg)){
					lattice_mdl[i][j] = 1;	//�c�Ƀ����N��L�΂�
					lattice_mdl[j][i] = 1;
				}
			}
		}
	}
}

void tr_func(void){	//�X�ɓ��������蓖�Ă�
	int i,j;
	
	for(i=0;i<ag_num;i++){
		initialization_1(ag[i].f,fe);
	}
	
	for(i=0;i<ag_num;i++){
		
		for(j=0;j<fe;j++){
			ag[i].f[j] = (int)genrand_int31() % tr;	//[0,tr-1]�܂ł̊Ԃŗ��������
			/*
			while(ag[i].f[j] < 0){
				ag[i].f[j] = (int)genrand_real2()*10000 % tr;
			}*/
			
		}
	}
}

int KD_func(int *p, int*q){	//2�̔z��p��q���݂��ɋ��L���Ă��镶���̐���Ԃ��֐�
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
	
	for(i=0;i<ag_num;i++){	//�s��̏c
		initialization_1(olap[i],ag_num);	//�s��̉�
	}
	
	for(i=0;i<ag_num;i++){
		for(j=0;j<ag_num;j++){
			if(lattice_mdl[i][j]==1){
				olap[i][j] = KD_func(ag[i].f,ag[j].f);
			}
		}
	}
}

void latcd_cpyfunc(void){	//���M�����[�i�q�̗אڍs��𕶉����z���ĕω�����̂Ɏg���אڍs��ɗv�f���R�s�[����
	int i,j;
	
	for(i=0;i<ag_num;i++){
		for(j=0;j<ag_num;j++){
			cd_nw[i][j] = lattice_mdl[i][j];
		}
	}
}

//�^����ꂽ�z��̒��g��S��0�ɂ���
void initialization_1(int *p,int x){	//p��1�����z��̃|�C���^�Cx�͔z��̑傫��
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

//�[���D��T��
void dfs1(int v){
	int i;
	
	node_through_list[v] = 1;
	
	for(i=0;i<ag_num;i++){
		if((cd_nw[v][i] == 1 )&& (node_through_list[i] == 0)){
			dfs1(i);
		}
	}
}

//�[���D��T���֐�(�߂�l�Ȃ�ver(�ł��傫�ȘA���O���t�������쐬)
void dfs2(int v){
	int i;
	
	node_through_list[v] = 1;
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[v][i] == 1){
			max_size_network[v][i] = 1;
			//max_size_network[i][v] = 1;	//�����Ȃ��Ă����퓮�삷��
			
			if(node_through_list[i] == 0){
				dfs2(i);
			}
		}
	}
}

int network_island_count(void){
	int count = 0;	//���̐�
	int i,j,dummy;
	int af_node_through_list_number = 0;
	int bf_node_through_list_number = 0;
	
	initialization_1(node_through_list,ag_num);	//������
	initialization_1(node_list_of_network,2);
	
	for(i=0;i<ag_num;i++){
		if(node_through_list[i] == 0){
			af_node_through_list_number = 0;
			
			dfs1(i);	//�[���D��T���֐�
			
			for(j=0;j<ag_num;j++){
				if(node_through_list[j] == 1){
					af_node_through_list_number++;
				}
			}
			
			af_node_through_list_number = af_node_through_list_number - bf_node_through_list_number;
			bf_node_through_list_number = af_node_through_list_number + bf_node_through_list_number;
			
			dummy = af_node_through_list_number;
			
			//printf("%d\n",dfs(i));
			
			if(dummy > node_list_of_network[1]){	//�ł��傫�����������L�^����
				node_list_of_network[0] = count;	//����ڂ�dfs1�ɂ́C
				node_list_of_network[1] = dummy;	//�����̃m�[�h���W�܂��Ă��邩
				
				//printf("test1\n");
			}
			//printf("test2\n");
			
			count++;
		}
	}
	
	return count;
}

int culture_count_function(void){	//�����̎�ނ𐔂���֐�
	agent *kinds_of_culture_list;
	int i,j=0,k,m,count=0;
	
	for(i=0;i<ag_num;i++){				if(i==0){
			kinds_of_culture_list = (agent *)malloc(sizeof(agent)*(j+1));	//�����̎��_�ł�j=0�Ȃ̂�j+1���Ĕz������
			
			for(k=0;k<fe;k++){
				kinds_of_culture_list[j].f[k] = ag[i].f[k];
			}
			
			j++;
		}
		else{
			for(m=0;m<j;m++){
				//printf("KD=%d\n",KD_func(kinds_of_culture_list[m].f, ag[i].f));
				if(KD_func(kinds_of_culture_list[m].f, ag[i].f) != fe){
					count++;	//����count��2��for���̒��ɂ��邱�Ƃɒ��ӂ���
				}
			}
			//printf("j=%d,m=%d,test\n",j,m);
			if(count==j){	//kinds_of_culture_list�Ɋ܂܂�Ȃ�����ag[i]�Ȃ��
				if((kinds_of_culture_list = (agent *)realloc(kinds_of_culture_list,sizeof(agent)*(j+1))) == NULL){
					printf("realloc���Ƀ��������m�ۂł��܂���\n");
					free(kinds_of_culture_list);  /* ����kinds_of_culture_list��������ďI�� */
					
					exit(1);
				}
				else{
					for(k=0;k<fe;k++){
						kinds_of_culture_list[j].f[k] = ag[i].f[k];
					}
				}
				j++;
			}
			//printf("test\n");
			
			count = 0;
		}
	}
	//printf("test_end\n");	//�����܂ōs�������Ȃ�
	
	free(kinds_of_culture_list);
	
	return j;
}

void max_size_network_func(void){
	int count = 0;
	int i;
	
	initialization_1(node_through_list,ag_num);	//������
	
	for(i=0;i<ag_num;i++){	//�s��̏c
		initialization_1(max_size_network[i],ag_num);	//�s��̉�
	}
	
	for(i=0;i<ag_num;i++){
		if(node_through_list[i] == 0){
			if(node_list_of_network[0] == count){
				//printf("node_list_of_network[0] = %d\n",node_list_of_network[0]);
				//printf("test3\n");
				dfs2(i);	//�[���D��T���֐�(�ł��傫�ȘA���O���t�������쐬)
			}
			else dfs1(i);	//�[���D��T���֐�
			
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
	
	//printf("debag_lineA\n");
	
	while(n==0){	//�ǓƂȃG�[�W�F���g(�אڃm�[�h��0)�͏����ė��������
		x = (int)(genrand_real1()*(double)(ag_num-1 -0 +1)/(double)ag_num);	//������������agent�������_���ɑI�o
	
		for(i=0;i<ag_num;i++){
			if(cd_nw[x][i] == 1){
				//printf("test by Centola_AA\n");
				n = n+1;
			}
		}
	}
	
	//printf("debag_lineB\n");
	
	nb_lst = (int *)malloc(sizeof(int) * n);	//nb_lst: neighbor agent list�̗�
	
	//printf("debag_lineC\n");
	
	for(i=0;i<ag_num;i++){
		if(cd_nw[x][i] == 1){
			nb_lst[m] = i;	//nb_lst��agent x�̗אڃm�[�h���X�g
			
			m++;
		}
	}
	
	//printf("debag_lineD\n");
	
	y = (int)genrand_int31() % n;	// 0���ɒ���
	
	//printf("debag_lineE\n");
	
	if((0 < KD_func(ag[x].f, ag[ nb_lst[y] ].f)) && (KD_func(ag[x].f, ag[ nb_lst[y] ].f) < fe)){	//���ݍ�p�Ώۂ��S�Ď����Ɠ��������������Ă����ꍇ���������ɏI���
		//printf("test by olgmdl_mdl_A\n");
		
		
		//printf("%d\n",x);
		
		
		//���L���Ă��Ȃ������������_���ɑI��
		b = (int)genrand_int31()%fe;
		
		while(ag[x].f[b] == ag[ nb_lst[y] ].f[b]){
			b = (int)genrand_int31()%fe;
		}
		
		//������agent��b�Ԗڂ̕���(ag[x].f[b])�͑��̗א�agent�Ƌ��L���Ă��邩
		for(i=0;i<n;i++){
			if(ag[x].f[b] == ag[ nb_lst[i] ].f[b]){	//���L���Ă����ꍇ
				count++;	//������g���ƌX�̕ێ琫�p�����[�^���l���邱�Ƃ��ł���
			}
		}
		
		if(count > 0){	//�ێ琫���܂߂����ݍ�p������
			a = ((double)KD_func(ag[x].f, ag[ nb_lst[y] ].f)*(1-conservativeness_parameter))/(double)fe;	//�m��( O(i,j)*�ێ琫) /F���v�Z
		}
		else if(count == 0){
			a = (double)KD_func(ag[x].f, ag[ nb_lst[y] ].f)/(double)fe;	//�m��O(i,j)/F���v�Z
		}
		
		//printf("%f/%f=%f\n",(double)KD_func(ag[x].f, ag[ nb_lst[y] ].f),(double)fe,a);
		//printf("test by Centola_mdl_B\n");
		
		if(a >= genrand_real1()){	//�����ގ��x�ɂ��m�������������ꍇ
			//printf("test by Centola_mdl_C\n")
			
			ag[x].f[b] = ag[ nb_lst[y] ].f[b];
		}
	}
	
	if(KD_func(ag[x].f, ag[ nb_lst[y] ].f) == 0){	//�ގ���=0�̂Ƃ�
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
	//printf("link -> %d\n",all_link_count);
	
	return a_rate_of_bridge/(double)all_link_count;
}

int vertex_list_check(int *p, int array_length){	//�z��̒��g�����ׂ�1�ł����0��Ԃ��C����ȊO�ł����1��Ԃ�
	int i;
	
	for(i=0;i<array_length;i++){
		if(*p==0){
			return 1;
		}
		p++;
	}
	return 0;	//*p ���S�� 1 �Ȃ��return 0
}

int dijkstra(int node_u,int node_v,int *ver,int *dis,int length){	//�_�C�N�X�g���@�̖{�v�Z�v���O����
	int i;
	int queue[ag_num];	//���钸�_�Ɋւ���אڃm�[�h���i�[�����
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
		
		if(dummy == INT_MAX){	//�������ŏ��l��INT_MAX�ł���΂���ȏ�T������Ӗ��͖����i->�s�\���B�_�����Ȃ��j
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
			
			for(i=0;i<length;i++){	//�ϐ� i �ɂ͗אڃm�[�h������
				if(max_size_network[node_u][i]==1){
					if(dis[i] > (dis[node_u] + max_size_network[node_u][i])){
						dis[i] = dis[node_u] + max_size_network[node_u][i];
					}
				}
			}
		}
		//������
		initialization_1(queue,ag_num);
		head = 0;
		tail = 0;
	}
	
	return dis[node_v];
}

int shortest_path_func(int node_u, int node_v){	//���_u���璸�_v�܂ł̍ŒZ���_�ԋ������v�Z���C���̋�����߂�l�Ƃ���
	int vertex_list[ag_num];
	int distance_list[ag_num];
	int prev = node_u;
	int dummy;
	int i;
	
	//������
	for(i=0;i<ag_num;i++){
		//distance_list�̏�����
		if(i==node_u){
			distance_list[i] = 0;
		}
		else{
			distance_list[i] = INT_MAX;
		}
		
		//vertex_list�̏�����
		vertex_list[i] = 0;	//vertex_list[i]==1�Ȃ�Β��_ i �͒��_���X�g���珜�O����Ă�����̂Ƃ���(0�Ȃ�Β��_���X�g�̗v�f)
	}
	//�����܂ŏ�����
	
	//�{�v�Z
	return dijkstra(node_u, node_v, vertex_list, distance_list, ag_num);
}

double average_vertex_distance(void){	//���ϒ��_�ԋ������v�Z����v���O����(d(i,j)��INF�̏ꍇ�͐��ɓ���Ȃ�)
	int count=0;	//���_�ԋ�����INF�̒��_�̐��𐔂���
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

double cluster_calculate_func(int node_i){
	int j,k;
	int sum=0;
	int degree_num = degree_count(node_i);
	
	if(degree_num<2){	//������0��1�ł���΃N���X�^�[�͌`���ł��Ȃ��̂�0��Ԃ�
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
	
	sum = sum / 2;	//�����O�p�`���d�����Đ����Ă���̂�2�Ŋ���
	
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
	//"average = sum / ag_num"�ł͂Ȃ����G���[���N����
	dummy = ag_num;
	average = sum / dummy;
	//sleep(1);
	//printf("average = %f\n",average);
	return average;
}

//�������Ɋւ���data set��ǂݍ���
void read_data_set_func(char *filename, double *cooperation_parameter, double *island_number_av_array, double *agent_number_av_array, double *culture_number_av_array, double *cluster_number_av_array, double *shortest_path_av_array, double *bridge_number_av_array, double *island_number_variance_array, double *agent_number_variance_array, double *culture_number_variance_array, double *cluster_number_variance_array, double *shortest_path_variance_array, double *bridge_number_variance_array, double *island_number_standard_deviation_array,  double *agent_number_standard_deviation_array, double *culture_number_standard_deviation_array, double *cluster_number_standard_deviation_array, double *shortest_path_standard_deviation_array, double *bridge_number_standard_deviation_array){
	FILE *fp;
	
	fp = fopen( filename, "r" );
	if( fp == NULL ){
		printf( "%s�t�@�C�����J���܂���\n", filename );
		exit(1);
	}
	else{
		while(fscanf( fp, "%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf ,%lf",  cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array, island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array, island_number_standard_deviation_array,  agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array, bridge_number_standard_deviation_array) != EOF ){
			//printf("�������F%f\n",*cooperation_parameter);
			
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
/* �z��o�[�W���� */
/*
//�������Ɋւ���data set��ǂݍ���
void read_data_set_func(char *filename, double cooperation_parameter[], double island_number_av_array[], double agent_number_av_array[], double culture_number_av_array[], double cluster_number_av_array[], double shortest_path_av_array[], double bridge_number_av_array[], double island_number_variance_array[], double agent_number_variance_array[], double culture_number_variance_array[], double cluster_number_variance_array[], double shortest_path_variance_array[], double bridge_number_variance_array[], double island_number_standard_deviation_array[],  double agent_number_standard_deviation_array[], double culture_number_standard_deviation_array[], double cluster_number_standard_deviation_array[], double shortest_path_standard_deviation_array[], double bridge_number_standard_deviation_array[]){
	FILE *fp;
	int index = 0;
	
	fp = fopen( filename, "r" );
	if( fp == NULL ){
		printf( "%s�t�@�C�����J���܂���\n", filename );
		exit(1);
	}
	else{		
		while(fscanf( fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",   &cooperation_parameter[index],  &island_number_av_array[index],  &agent_number_av_array[index],  &culture_number_av_array[index],  &cluster_number_av_array[index],  &shortest_path_av_array[index],  &bridge_number_av_array[index],  &island_number_variance_array[index],  &agent_number_variance_array[index],  &culture_number_variance_array[index],  &cluster_number_variance_array[index],  &shortest_path_variance_array[index],  &bridge_number_variance_array[index],  &island_number_standard_deviation_array[index],   &agent_number_standard_deviation_array[index],  &culture_number_standard_deviation_array[index],  &cluster_number_standard_deviation_array[index],  &shortest_path_standard_deviation_array[index],  &bridge_number_standard_deviation_array[index]) != EOF ){
			printf("�������F%f\n",cooperation_parameter[index]);
			index++;
		}
	}
	fclose( fp );
}*/
 
int cooperation_parameter_return_index(double *cooperation_parameter){
	int i;
	
	for(i=0;i<data_set_row;i++){
		//printf("%f\n",cooperation_parameter[data_set_row]);
		
		if(cooperation_parameter[i] == conservativeness_parameter){
			return i;
		}
		else if(cooperation_parameter[i] > conservativeness_parameter){
			printf("data set�ɑ΂��鋤�����l(conservativeness_parameter)���������Ȃ��̂ŋ߂��l%f��Ԃ��܂���",cooperation_parameter[i]);
			
			return i;
		}
	}
	
	printf("�����I�����܂�%d line)",__LINE__);
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

/* main�� */
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
	
	//data set�N���X�̓ǂݍ���
	char filename5[] = "data_set20151112.csv";
	char filename7[] = "data_set_torus20151214.csv";
	int index_mark;
	
	/* ���������i�[����z�� */
	double cooperation_parameter[100];
	
	/* ���ς��i�[����z�� */
	double island_number_av_array[100];
	double agent_number_av_array[100];
	double culture_number_av_array[100];
	double cluster_number_av_array[100];
	double shortest_path_av_array[100];
	double bridge_number_av_array[100];
	
	/* ���U���i�[����z�� */
	double island_number_variance_array[100];
	double agent_number_variance_array[100];
	double culture_number_variance_array[100];
	double cluster_number_variance_array[100];
	double shortest_path_variance_array[100];
	double bridge_number_variance_array[100];
	
	/* �W���΍����i�[����z�� */
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
	
	wt[0] = 0;	//�����o��t�̎���(1���)
	wt[1] = 2500;	//�����o��t�̎���(2���)
	wt[2] = 25000;	//�����o��t�̎���(3���)
	wt[3] = 500000;	//�����o��t�̎���(4���)
	
	init_genrand((unsigned)time(NULL));
	
	while(1){
		t=0;	//t �̏�����
		lattice_func();
		tr_func();
		olap_func();
		latcd_cpyfunc();	//���M�����[�i�q�̗אڍs��𕶉����z���ĕω�����̂Ɏg���אڍs��ɗv�f���R�s�[����
		
		/* �t�@�C���I�[�v�� */
		if ((fp2 = fopen(filename6, "wt")) == NULL) {
			fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
			return exit(1);
		}
		
		while(t<=tmax){
			//printf("%d event\n",t);
			
			for(k=0;k<wt_num;k++){
				if(wt[k] == t){
					for(i=0;i<ag_num;i++){	//�������ޑO��olap�𐳂����v�Z����
						for(j=0;j<ag_num;j++){
							if (cd_nw[i][j]==1){
								olap[i][j] = KD_func(ag[i].f,ag[j].f);
							}
						}
					}
					
					if(t==0){	//t=0(�i�q��)�̃l�b�g���[�N�������o��
						edg=0;	//edg�̏�����
						
						/* gephi faile���� */
						sprintf(filename,"nw_%d.gexf",wt[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							return exit(1);
						}
						
						/* �������� */
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
						
						/*
						if ( count < 0 ) {
							fprintf(stderr, "�t�@�C���̏����݂Ɏ��s���܂���.\n");
							fclose(fp);
							return exit(1);
						}
						*/
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file���� */
						sprintf(filename2,"nw_%d(sociarium).txt",wt[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							return exit(1);
						}
						
						/* �������� */
						fprintf(fp,"# Pajek�`���̃f�[�^ \n");
						fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
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
						
						/*
						if ( count < 0 ) {
							fprintf(stderr, "�t�@�C���̏����݂Ɏ��s���܂���.\n");
							fclose(fp);
							return exit(1);
						}
						*/
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
					}
					else{	//����t�̃l�b�g���[�N�������o��
						edg=0;	//edg�̏�����
						/* gephi file�쐬 */
						sprintf(filename,"nw_%d.gexf",wt[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							return exit(1);
						}
						
						/* �������� */
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
						
						/*
						if ( count < 0 ) {
							fprintf(stderr, "�t�@�C���̏����݂Ɏ��s���܂���.\n");
							fclose(fp);
							return exit(1);
						}
						*/
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file���� */
						sprintf(filename2,"nw_%d(sociarium).txt",wt[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							return exit(1);
						}
						
						/* �������� */
						fprintf(fp,"# Pajek�`���̃f�[�^ \n");
						fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
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
						
						/*
						if ( count < 0 ) {
							fprintf(stderr, "�t�@�C���̏����݂Ɏ��s���܂���.\n");
							fclose(fp);
							return exit(1);
						}
						*/
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
					}
				}
			}
			
			/* sociarium file���� */	
			if(t%sociarium_ptch_width==0){
				if(t==0){
					/* �������� */
					fprintf(fp2,"@module = graph_creation_read_time_series_rect.dll \n");
					fprintf(fp2,"@title = cultural drift %f\n",conservativeness_parameter);
					fprintf(fp2,"@delimiter =\t\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@time_format = Y.M.D # Y[ear], M[onth], D[ay], h[our], m[inute], s[econd] \n");
					fprintf(fp2,"@interval  = %d # �l�b�g���[�N���쐬���鎞���̊Ԋu\n",sociarium_ptch_width);
					fprintf(fp2,"@characteristic_time = 1 # �e�����ɂ���exp(-dt/@characteristic_time)��ώZ���C@threshold�ȉ��̃G�b�W���}�X�L���O \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@threshold = 0\n");
					fprintf(fp2,"@nondirected \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@node_texture = agent.png\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", 0, 0, 0, 0);	//���ꂪ�Ȃ��Ɛ������`�悪�ł��Ȃ�
					
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
			
			/***********�v���O�����̊j(orgnl_mdl)�ƂȂ����***********
			---------------------�ȉ���̓I�ȓ��e-----------------------
			1.������agent�������_���ɑI�������
			2.����agent�̗א�agent�������_����1�I������
			3.�ގ������v�Z����
			3.���ݍ�p���s��
				-�א�agent����m��O(i,j)/F�ŕ����̉e�����󂯂�
			4.overlap��0�Ȃ��rule5�ɑ���
			***********************************************************/
			orgnl_mdl();
			
			t = t+1;
		}
		/* �t�@�C���N���[�Y */
		fclose(fp2);
		
		max_size_network_func();	//�ő�̓��̗אڍs��𐶐�
		
		island_dummy = network_island_count();
		distance_dummy = average_vertex_distance();
		cluster_dummy = average_cluster_calculate_func();
		
		if(read_data_set){
			if(simulation_count >= SimulationStop){
				printf("�V�~�����[�V�����������I�����܂�(Error : %d line)\n",__LINE__);
				exit(1);
			}			
			//�Y����parameter��T��
			
		/* ���O�t�@�C������ */
		/* �t�@�C���I�[�v�� */
			if ((fp = fopen("simulation_log.csv", "a")) == NULL) {
				fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
				exit(1);
			}
			
			/* �������� */
			fprintf(fp,"%d,%f,%f \n",island_dummy,distance_dummy,cluster_dummy);
			
			fclose(fp);
			
			if((island_number_av_array[index_mark] - sigma_option*island_number_standard_deviation_array[index_mark] < island_dummy)&&(island_dummy< island_number_av_array[index_mark] + sigma_option*island_number_standard_deviation_array[index_mark])&&(shortest_path_av_array[index_mark] - sigma_option*shortest_path_standard_deviation_array[index_mark] < distance_dummy)&&(distance_dummy < shortest_path_av_array[index_mark] + sigma_option*shortest_path_standard_deviation_array[index_mark])&&(cluster_number_av_array[index_mark] - sigma_option*cluster_number_standard_deviation_array[index_mark] < cluster_dummy)&&(cluster_dummy < cluster_number_av_array[index_mark] + sigma_option*cluster_number_standard_deviation_array[index_mark])){
				break;
			}
			simulation_count++;
			
			printf("�ăV�~�����[�V�������Ă��܂��c %d / %d\n",simulation_count,SimulationStop);
		}
		else{
			break;
		}
	}
	
	/* �t�@�C���I�[�v�� */
	edg = 0;
	if ((fp = fopen(filename3, "w")) == NULL) {
		fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
		return exit(1);
	}
	
	/////////////////////////���̃l�b�g���[�N�̂����ł��傫�ȃl�b�g���[�N��`�悷��/////////////////////////
	//Gephi_file
	/* �������� */
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
	
	/* �t�@�C���N���[�Y */
	fclose(fp);
	
	/* sociarium file���� */	
	/* �t�@�C���I�[�v�� */
	if ((fp = fopen(filename4, "wt")) == NULL) {
		fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
		return exit(1);
	}
	
	/* �������� */
	fprintf(fp,"# Pajek�`���̃f�[�^ \n");
	fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
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
	
	/* �t�@�C���N���[�Y */
	fclose(fp);
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		//�����̗��z�����������̂��������X�g���쐬���Ċm�F����
	/* �t�@�C���I�[�v�� */
	if ((fp = fopen(filename1, "w")) == NULL) {
		fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
		return exit(1);
	}
	
	/* �������� */
	fprintf(fp,"���̐� = %3d \n",island_dummy);
	fprintf(fp,"�ł��傫������agent�� = %3d \n",node_list_of_network[1]);
	fprintf(fp,"�����̐� = %3d \n",culture_count_function());
	fprintf(fp,"���ϒ��_�ԋ��� = %3f \n",distance_dummy);
	fprintf(fp,"���σN���X�^�[�W�� = %3f \n",cluster_dummy);
	fprintf(fp,"�u���b�W��(%%) = %f \n",a_rate_of_bridge_func()*100);
	
	fclose(fp);
	
	//printf("%d\n",node_list_of_network[0]);
	//printf("%d\n",node_list_of_network[1]);
	
	return 0;
}
