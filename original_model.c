﻿/* �w�b�_�[�t�@�C�� */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "MT.h"

/* �v���v���Z�b�T */
#define TRUE 1
#define FALSE 0
#define latticeSideSize 20	//�i�q�̏c�܂��͉��̒����D
#define features 3	//features�D�����̐�
#define traits 20	//traits�D�����͈̔�
#define agents latticeSideSize*latticeSideSize	//agent��;
#define eventMax 500000	//�C�x���g��
#define writeOutFileNumber 4	//�����o���t�@�C���̐�
#define convertCultureToRgb 255/traits	//agent��rgb�Œl��������
#define latticeToTorus FALSE	//�i�q���f�����g�[���X�ɂ���ꍇ��TRUE�ɂ���(���Ȃ��ꍇ��FALSE)
#define cooperativityParameter 0.95	//�ێ琫�p�����[�^(0�ȏ�1�ȉ��Ő�����I��) ->0(�ێ琫�������Ȃ����)��I�Ԃ�Centola���f���ɂȂ�
#define readDataSet TRUE	//data_set��ǂݍ���data set�̂��鋤�����ɑ΂��鑮���l(colum)�̕��ςɋ߂��O���t�𐶐�����i�ǂݍ��܂Ȃ��ꍇ��FALSE�j
	#define dataSetCooperativityStepSize 100	//data set�̋�������0����1�܂ł̍��݂̌�
	#define sigmaConfficient 1	//�O���t�����ɂ���������͈� �W���΍�(sigma)��1->68.27%�C2->95.45%�C3->99.73%
	#define SimulationStop 100	//�V�~�����[�V�����񐔁D���̉񐔂𒴂����ꍇ�����I�����܂��D
#define sociariumAnimationStepSize 5000	//sociarium�̃O���t���A�j���[�V�����ɂ���ہC���Ԏ������P�ʂŐi�߂Ă���������

/* �\���� */
typedef struct{
	int feature[features];	//����
} createAgentType;


/* �O���ϐ� */
int latticeGraph[agents][agents];	//�i�q���f���̗אڍs��
int evolvingNetwork[agents][agents];	//cultural drift����l�b�g���[�N�̗אڍs��
int networkAddedSimilarity[agents][agents];	//overlap
int writeOutEventTimes[writeOutFileNumber];	//���g��main���ɏ����Ă���
createAgentType agent[agents];
int throughNodeListInDfs[agents];	//dfs�Œʉ߂����m�[�h�͊Y��������1������
int connectedGraphInformation[2];	//�Y��0�̔z��ɂ́u���ڂ̘A���O���t���v���C�Y��1�ɂ́u���ׂĂ̘A���O���t�̒��ōł��傫���O���t�̃m�[�h���v���i�[
int maxSizeNetwork[agents][agents];	//�A���O���t�̂����ł�������agent�������O���t��ۑ�����


/* �v���g�^�C�v�錾 */
void createLatticeGraph(void);
void convertLatticeToTorus(void);
void assignAgentToTrait(void);
int countCommonFeature(int *p, int*q);
void createNetworkAddedSimilarity(void);
void copyLatticeGraphToEvolvingNetwork(void);
void initializeIntArray(int *p,int x);
void dfs1(int v);
void dfs2(int v);
int countIslandNetwork(void);
int countCulture(void);
void createMaxSizeNetwork(void);
void orgnl_mdl(void);
double countNetworkBridge(void);
int reachVertexCheck(int *p, int array_length);
int dijkstra(int node_u,int node_v,int *ver,int *dis,int length);
int calculateVertexDistance(int node_u, int node_v);
double calculateAverageVertexDistance(void);
double calculateClusterCoefficient(int node_i);
int calculateDegree(int node_i);
double calculateAverageClusterCoefficient(void);
void writeDataSetOnMemory(char *filename, double *cooperation_parameter, double *island_number_av_array, double *agent_number_av_array, double *culture_number_av_array, double *cluster_number_av_array, double *shortest_path_av_array, double *bridge_number_av_array, double *island_number_variance_array, double *agent_number_variance_array, double *culture_number_variance_array, double *cluster_number_variance_array, double *shortest_path_variance_array, double *bridge_number_variance_array, double *island_number_standard_deviation_array,  double *agent_number_standard_deviation_array, double *culture_number_standard_deviation_array, double *cluster_number_standard_deviation_array, double *shortest_path_standard_deviation_array, double *bridge_number_standard_deviation_array);
int network_edges_count(void);


/* �O���֐� */
void createLatticeGraph(void){	//�i�q���f���쐬�A���S���Y��
	int i,j;
	
	/* ������ */
	for(i=0;i<agents;i++){	//�s��̏c
		initializeIntArray(latticeGraph[i],agents);	//�s��̉�
		initializeIntArray(evolvingNetwork[i],agents);	//�s��̉�
	}
	
	for(i=0;i<agents;i++){
		for(j=i;j<agents;j++){
			if(i == j){
				latticeGraph[i][j] = 0;
			}
			else{
				if( ((i+1) == j) && ((i+1)%latticeSideSize != 0)){
					latticeGraph[i][j] = 1;	//�E�Ƀ����N��L�΂�
					latticeGraph[j][i] = 1;
				}
				if( ((j-i) == latticeSideSize) && ((i+latticeSideSize) < agents)){
					latticeGraph[i][i+latticeSideSize] = 1;	//���Ƀ����N��L�΂�
					latticeGraph[i+latticeSideSize][i] = 1;
				}
			}
		}
	}
	
	if(latticeToTorus){
		convertLatticeToTorus();
	}
}

void convertLatticeToTorus(void){
	int i,j;
	
	for(i=0;i<agents;i++){
		for(j=i;j<agents;j++){
			if(i == j){
				latticeGraph[i][j] = 0;
			}
			else{
				if( ((j+1-i) == latticeSideSize) && (i%latticeSideSize==0)){
					latticeGraph[i][j] = 1;	//���Ƀ����N��L�΂�
					latticeGraph[j][i] = 1;
				}
				if( ((j-i) == (latticeSideSize*latticeSideSize-latticeSideSize)) && (i < latticeSideSize)){
					latticeGraph[i][j] = 1;	//�c�Ƀ����N��L�΂�
					latticeGraph[j][i] = 1;
				}
			}
		}
	}
}

void assignAgentToTrait(void){	//�X�ɓ��������蓖�Ă�
	int i,j;
	
	for(i=0;i<agents;i++){
		initializeIntArray(agent[i].feature,features);
	}
	
	for(i=0;i<agents;i++){
		
		for(j=0;j<features;j++){
			agent[i].feature[j] = (int)genrand_int31() % traits;	//[0,traits-1]�܂ł̊Ԃŗ��������
		}
	}
}

int countCommonFeature(int *p, int*q){	//2�̔z��p��q���݂��ɋ��L���Ă��镶���̐���Ԃ��֐�
	int i,n=0;
	
	for(i=0;i<features;i++){
		if(*p == *q){
			n = n + 1;
		}
		p++;
		q++;
	}
	
	return n;
}

void createNetworkAddedSimilarity(void){
	int i,j;
	
	for(i=0;i<agents;i++){	//�s��̏c
		initializeIntArray(networkAddedSimilarity[i],agents);	//�s��̉�
	}
	
	for(i=0;i<agents;i++){
		for(j=0;j<agents;j++){
			if(latticeGraph[i][j]==1){
				networkAddedSimilarity[i][j] = countCommonFeature(agent[i].feature,agent[j].feature);
			}
		}
	}
}

void copyLatticeGraphToEvolvingNetwork(void){	//���M�����[�i�q�̗אڍs��𕶉����z���ĕω�����̂Ɏg���אڍs��ɗv�f���R�s�[����
	int i,j;
	
	for(i=0;i<agents;i++){
		for(j=0;j<agents;j++){
			evolvingNetwork[i][j] = latticeGraph[i][j];
		}
	}
}

//�^����ꂽ�z��̒��g��S��0�ɂ���
void initializeIntArray(int *p,int x){	//p��1�����z��̃|�C���^�Cx�͔z��̑傫��
	int i;
	
	for(i=0;i<x;i++){
		*p = 0;
		
		p++;
	}
}

//�[���D��T��
void dfs1(int v){
	int i;
	
	throughNodeListInDfs[v] = 1;
	
	for(i=0;i<agents;i++){
		if((evolvingNetwork[v][i] == 1 )&& (throughNodeListInDfs[i] == 0)){
			dfs1(i);
		}
	}
}

//�[���D��T���֐�(�߂�l�Ȃ�ver(�ł��傫�ȘA���O���t�������쐬)
void dfs2(int v){
	int i;
	
	throughNodeListInDfs[v] = 1;
	
	for(i=0;i<agents;i++){
		if(evolvingNetwork[v][i] == 1){
			maxSizeNetwork[v][i] = 1;
			//maxSizeNetwork[i][v] = 1;	//�����Ȃ��Ă����퓮�삷��
			
			if(throughNodeListInDfs[i] == 0){
				dfs2(i);
			}
		}
	}
}

int countIslandNetwork(void){
	int count = 0;	//���̐�
	int i,j,dummy;
	int af_node_through_list_number = 0;
	int bf_node_through_list_number = 0;
	
	initializeIntArray(throughNodeListInDfs,agents);	//������
	initializeIntArray(connectedGraphInformation,2);
	
	for(i=0;i<agents;i++){
		if(throughNodeListInDfs[i] == 0){
			af_node_through_list_number = 0;
			
			dfs1(i);	//�[���D��T���֐�
			
			for(j=0;j<agents;j++){
				if(throughNodeListInDfs[j] == 1){
					af_node_through_list_number++;
				}
			}
			
			af_node_through_list_number = af_node_through_list_number - bf_node_through_list_number;
			bf_node_through_list_number = af_node_through_list_number + bf_node_through_list_number;
			
			dummy = af_node_through_list_number;
			
			if(dummy > connectedGraphInformation[1]){	//�ł��傫�����������L�^����
				connectedGraphInformation[0] = count;	//����ڂ�dfs1�ɂ́C
				connectedGraphInformation[1] = dummy;	//�����̃m�[�h���W�܂��Ă��邩
			}
			
			count++;
		}
	}
	
	return count;
}

int countCulture(void){	//�����̎�ނ𐔂���֐�
	createAgentType *kinds_of_culture_list;
	int i,j=0,k,m,count=0;
	
	for(i=0;i<agents;i++){				if(i==0){
			kinds_of_culture_list = (createAgentType *)malloc(sizeof(createAgentType)*(j+1));	//�����̎��_�ł�j=0�Ȃ̂�j+1���Ĕz������
			
			for(k=0;k<features;k++){
				kinds_of_culture_list[j].feature[k] = agent[i].feature[k];
			}
			
			j++;
		}
		else{
			for(m=0;m<j;m++){
				if(countCommonFeature(kinds_of_culture_list[m].feature, agent[i].feature) != features){
					count++;	//����count��2��for���̒��ɂ��邱�Ƃɒ��ӂ���
				}
			}
			if(count==j){	//kinds_of_culture_list�Ɋ܂܂�Ȃ�����ag[i]�Ȃ��
				if((kinds_of_culture_list = (createAgentType *)realloc(kinds_of_culture_list,sizeof(createAgentType)*(j+1))) == NULL){
					printf("realloc���Ƀ��������m�ۂł��܂���\n");
					free(kinds_of_culture_list);  /* ����kinds_of_culture_list��������ďI�� */
					
					exit(1);
				}
				else{
					for(k=0;k<features;k++){
						kinds_of_culture_list[j].feature[k] = agent[i].feature[k];
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

void createMaxSizeNetwork(void){
	int count = 0;
	int i;
	
	initializeIntArray(throughNodeListInDfs,agents);	//������
	
	for(i=0;i<agents;i++){	//�s��̏c
		initializeIntArray(maxSizeNetwork[i],agents);	//�s��̉�
	}
	
	for(i=0;i<agents;i++){
		if(throughNodeListInDfs[i] == 0){
			if(connectedGraphInformation[0] == count){
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
	
	while(n==0){	//�ǓƂȃG�[�W�F���g(�אڃm�[�h��0)�͏����ė��������
		x = (int)(genrand_real1()*(double)(agents-1 -0 +1)/(double)agents);	//������������agent�������_���ɑI�o
	
		for(i=0;i<agents;i++){
			if(evolvingNetwork[x][i] == 1){
				n = n+1;
			}
		}
	}
	
	nb_lst = (int *)malloc(sizeof(int) * n);	//nb_lst: neighbor agent list�̗�
	
	for(i=0;i<agents;i++){
		if(evolvingNetwork[x][i] == 1){
			nb_lst[m] = i;	//nb_lst��agent x�̗אڃm�[�h���X�g
			
			m++;
		}
	}
	
	y = (int)genrand_int31() % n;	// 0���ɒ���
	
	if((0 < countCommonFeature(agent[x].feature, agent[ nb_lst[y] ].feature)) && (countCommonFeature(agent[x].feature, agent[ nb_lst[y] ].feature) < features)){	//���ݍ�p�Ώۂ��S�Ď����Ɠ��������������Ă����ꍇ���������ɏI���
	
		//���L���Ă��Ȃ������������_���ɑI��
		b = (int)genrand_int31()%features;
		
		while(agent[x].feature[b] == agent[ nb_lst[y] ].feature[b]){
			b = (int)genrand_int31()%features;
		}
		
		//������agent��b�Ԗڂ̕���(agent[x].feature[b])�͑��̗א�agent�Ƌ��L���Ă��邩
		for(i=0;i<n;i++){
			if(agent[x].feature[b] == agent[ nb_lst[i] ].feature[b]){	//���L���Ă����ꍇ
				count++;	//������g���ƌX�̕ێ琫�p�����[�^���l���邱�Ƃ��ł���
			}
		}
		
		if(count > 0){	//�ێ琫���܂߂����ݍ�p������
			a = ((double)countCommonFeature(agent[x].feature, agent[ nb_lst[y] ].feature)*(1-cooperativityParameter))/(double)features;	//�m��( O(i,j)*�ێ琫) /F���v�Z
		}
		else if(count == 0){
			a = (double)countCommonFeature(agent[x].feature, agent[ nb_lst[y] ].feature)/(double)features;	//�m��O(i,j)/F���v�Z
		}
		
		if(a >= genrand_real1()){	//�����ގ��x�ɂ��m�������������ꍇ
			agent[x].feature[b] = agent[ nb_lst[y] ].feature[b];
		}
	}
	
	if(countCommonFeature(agent[x].feature, agent[ nb_lst[y] ].feature) == 0){	//�ގ���=0�̂Ƃ�
		evolvingNetwork[x][ nb_lst[y] ] = 0;
		evolvingNetwork[ nb_lst[y] ][x] = 0;
		
		yy = (int)(genrand_real1()*(double)(agents-1 -0 +1)/(double)agents);
		
		while((x == yy) || ( nb_lst[y] == yy) || (evolvingNetwork[x][yy] == 1) || (evolvingNetwork[yy][x] == 1)){
			yy = (int)(genrand_real1()*(double)(agents-1 -0 +1)/(double)agents);
		}
		
		evolvingNetwork[x][yy] = 1;
		evolvingNetwork[yy][x] = 1;
	}
	
	free(nb_lst);
}

double countNetworkBridge(void){
	int a_rate_of_bridge=0;
	int all_link_count=0;
	int i,j;
	
	for(i=1;i<agents;i++){
		for(j=0;j<i;j++){
			if(evolvingNetwork[i][j]==1){
				all_link_count++;
				
				if((0 < countCommonFeature(agent[i].feature,agent[j].feature)) && (countCommonFeature(agent[i].feature,agent[j].feature) < features)){
					a_rate_of_bridge++;
				}
			}
		}
	}
	
	return a_rate_of_bridge/(double)all_link_count;
}

int reachVertexCheck(int *p, int array_length){	//�z��̒��g�����ׂ�1�ł����0��Ԃ��C����ȊO�ł����1��Ԃ�
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
	int queue[agents];	//���钸�_�Ɋւ���אڃm�[�h���i�[�����
	int head = 0, tail = 0;
	int node_adj;
	int dummy;
	
	while(reachVertexCheck(ver,agents)){
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
				if(maxSizeNetwork[node_u][i]==1){
					if(dis[i] > (dis[node_u] + maxSizeNetwork[node_u][i])){
						dis[i] = dis[node_u] + maxSizeNetwork[node_u][i];
					}
				}
			}
		}
		//������
		initializeIntArray(queue,agents);
		head = 0;
		tail = 0;
	}
	
	return dis[node_v];
}

int calculateVertexDistance(int node_u, int node_v){	//���_u���璸�_v�܂ł̍ŒZ���_�ԋ������v�Z���C���̋�����߂�l�Ƃ���
	int vertex_list[agents];
	int distance_list[agents];
	int prev = node_u;
	int dummy;
	int i;
	
	//������
	for(i=0;i<agents;i++){
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
	return dijkstra(node_u, node_v, vertex_list, distance_list, agents);
}

double calculateAverageVertexDistance(void){	//���ϒ��_�ԋ������v�Z����v���O����(d(i,j)��INF�̏ꍇ�͐��ɓ���Ȃ�)
	int count=0;	//���_�ԋ�����INF�̒��_�̐��𐔂���
	int sum=0;
	int d;
	int i,j;
	double av_dis;
	
	for(i=1;i<agents;i++){
		for(j=0;j<i;j++){
			d = calculateVertexDistance(i,j);
			
			if(INT_MAX == d){
				count++;
			}
			else{
				sum = d + sum;
			}
		}
	}
	
	av_dis = 2.0/(connectedGraphInformation[1]*(connectedGraphInformation[1]-1) ) * sum;
	
	return av_dis;
}

double calculateClusterCoefficient(int node_i){
	int j,k;
	int sum=0;
	int degree_num = calculateDegree(node_i);
	
	if(degree_num<2){	//������0��1�ł���΃N���X�^�[�͌`���ł��Ȃ��̂�0��Ԃ�
		return 0;
	}
	
	for(j=0;j<agents;j++){
		if(evolvingNetwork[node_i][j]){
			for(k=0;k<agents;k++){
				if(evolvingNetwork[node_i][k]){
					if(evolvingNetwork[j][k]){
						sum++;
					}
				}
			}
		}
	}
	
	sum = sum / 2;	//�����O�p�`���d�����Đ����Ă���̂�2�Ŋ���
	
	return (2.0 / (degree_num * (degree_num - 1))) * sum;
}

int calculateDegree(int node_i){
	int i;
	int count=0;
	
	for(i=0;i<agents;i++){
		if(evolvingNetwork[node_i][i]){
			count++;
		}
	}
	
	return count;
}

double calculateAverageClusterCoefficient(void){
	int node_i;
	double sum=0, average=0;
	int dummy;
	
	for(node_i=0;node_i<agents;node_i++){
		sum = calculateClusterCoefficient(node_i) + sum;
	}
	//"average = sum / agents"�ł͂Ȃ����G���[���N����
	dummy = agents;
	average = sum / dummy;
	
	return average;
}

//�������Ɋւ���data set��ǂݍ���
void writeDataSetOnMemory(char *filename, double *cooperation_parameter, double *island_number_av_array, double *agent_number_av_array, double *culture_number_av_array, double *cluster_number_av_array, double *shortest_path_av_array, double *bridge_number_av_array, double *island_number_variance_array, double *agent_number_variance_array, double *culture_number_variance_array, double *cluster_number_variance_array, double *shortest_path_variance_array, double *bridge_number_variance_array, double *island_number_standard_deviation_array,  double *agent_number_standard_deviation_array, double *culture_number_standard_deviation_array, double *cluster_number_standard_deviation_array, double *shortest_path_standard_deviation_array, double *bridge_number_standard_deviation_array){
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
 
int cooperation_parameter_return_index(double *cooperation_parameter){
	int i;
	
	for(i=0;i<dataSetCooperativityStepSize;i++){
		//printf("%f\n",cooperation_parameter[dataSetCooperativityStepSize]);
		
		if(cooperation_parameter[i] == cooperativityParameter){
			return i;
		}
		else if(cooperation_parameter[i] > cooperativityParameter){
			printf("data set�ɑ΂��鋤�����l(cooperativityParameter)���������Ȃ��̂ŋ߂��l%f��Ԃ��܂���",cooperation_parameter[i]);
			
			return i;
		}
	}
	
	printf("�����I�����܂�%d line)",__LINE__);
	exit(1);
}

int network_edges_count(void){
	int i,j;
	int count=0;
	
	for(i=0;i<agents;i++){
		for(j=i;j<agents;j++){
			if(evolvingNetwork[i][j] == 1){
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
	char filename3[] = "maxSizeNetwork.gexf";
	char filename4[] = "maxSizeNetwork(sociarium).txt";
	char filename6[] = "cultural_network_animation(sociarium).txt";
	int i,j,k;
	int count;
	int t=0,edg=0;
	int simulation_count=0;
	int island_dummy;
	double distance_dummy;
	double cluster_dummy;
	
	//data set�N���X�̓ǂݍ���
	char filename5[] = "data_set_not_torus20151112.csv";
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
	
	if(readDataSet){
		if(latticeToTorus == 0){
			writeDataSetOnMemory(filename5, cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array,  island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array,  island_number_standard_deviation_array, agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array,  bridge_number_standard_deviation_array);
			index_mark = cooperation_parameter_return_index(cooperation_parameter);
		}
		else if(latticeToTorus == 1){
			writeDataSetOnMemory(filename7, cooperation_parameter, island_number_av_array, agent_number_av_array, culture_number_av_array, cluster_number_av_array, shortest_path_av_array, bridge_number_av_array,  island_number_variance_array, agent_number_variance_array, culture_number_variance_array, cluster_number_variance_array, shortest_path_variance_array, bridge_number_variance_array,  island_number_standard_deviation_array, agent_number_standard_deviation_array, culture_number_standard_deviation_array, cluster_number_standard_deviation_array, shortest_path_standard_deviation_array,  bridge_number_standard_deviation_array);
			index_mark = cooperation_parameter_return_index(cooperation_parameter);
		}
	}
	
	writeOutEventTimes[0] = 0;	//�����o��t�̎���(1���)
	writeOutEventTimes[1] = 2500;	//�����o��t�̎���(2���)
	writeOutEventTimes[2] = 25000;	//�����o��t�̎���(3���)
	writeOutEventTimes[3] = 500000;	//�����o��t�̎���(4���)
	
	init_genrand((unsigned)time(NULL));
	
	while(1){
		t=0;	//t �̏�����
		createLatticeGraph();
		assignAgentToTrait();
		createNetworkAddedSimilarity();
		copyLatticeGraphToEvolvingNetwork();	//���M�����[�i�q�̗אڍs��𕶉����z���ĕω�����̂Ɏg���אڍs��ɗv�f���R�s�[����
		
		/* �t�@�C���I�[�v�� */
		if ((fp2 = fopen(filename6, "wt")) == NULL) {
			fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
			exit(1);
		}
		
		while(t<=eventMax){
			//printf("%d event\n",t);
			
			for(k=0;k<writeOutFileNumber;k++){
				if(writeOutEventTimes[k] == t){
					for(i=0;i<agents;i++){	//�������ޑO��olap�𐳂����v�Z����
						for(j=0;j<agents;j++){
							if (evolvingNetwork[i][j]==1){
								networkAddedSimilarity[i][j] = countCommonFeature(agent[i].feature,agent[j].feature);
							}
						}
					}
					
					if(t==0){	//t=0(�i�q��)�̃l�b�g���[�N�������o��
						edg=0;	//edg�̏�����
						
						/* gephi faile���� */
						sprintf(filename,"nw_%d.gexf",writeOutEventTimes[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							exit(1);
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
						
						for(i=0;i<agents;i++){
							fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
							fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",agent[i].feature[0]*convertCultureToRgb,agent[i].feature[1]*convertCultureToRgb,agent[i].feature[2]*convertCultureToRgb);
							fprintf(fp,"</node>\n");
						}
						
						fprintf(fp,"</nodes>\n");
						fprintf(fp,"<edges>\n");
						
						for(i=0;i<agents;i++){
							for(j=i;j<agents;j++){
								if(evolvingNetwork[i][j]==1){
									fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)networkAddedSimilarity[i][j]);
									edg = edg+1;
								}
							}
						}
						
						fprintf(fp,"</edges>\n");
						fprintf(fp,"</graph>\n");
						fprintf(fp,"</gexf>\n");
						
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file���� */
						sprintf(filename2,"nw_%d(sociarium).txt",writeOutEventTimes[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							exit(1);
						}
						
						/* �������� */
						fprintf(fp,"# Pajek�`���̃f�[�^ \n");
						fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
						fprintf(fp,"\n");
						fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
						fprintf(fp,"@title = The dissemination of culture \n");
						fprintf(fp,"@nondirected \n");
						fprintf(fp,"\n");
						
						fprintf(fp,"*Vertices %d \n",agents);
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<agents;i++){
							fprintf(fp,"%d \"agent_%d\" \n",i,i);
						}
						
						fprintf(fp,"*Arcs \n");
						
						for(i=0;i<agents;i++){
							for(j=i;j<agents;j++){
								if(evolvingNetwork[i][j]==1){
									fprintf(fp,"%d %d %0.1f \n",i,j,(double)networkAddedSimilarity[i][j]);
								}
							}
						}
						
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
					}
					else{	//����t�̃l�b�g���[�N�������o��
						edg=0;	//edg�̏�����
						/* gephi file�쐬 */
						sprintf(filename,"nw_%d.gexf",writeOutEventTimes[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							exit(1);
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
						
						for(i=0;i<agents;i++){
							fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
							fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",agent[i].feature[0]*convertCultureToRgb,agent[i].feature[1]*convertCultureToRgb,agent[i].feature[2]*convertCultureToRgb);
							fprintf(fp,"</node>\n");
						}
						
						fprintf(fp,"</nodes>\n");
						fprintf(fp,"<edges>\n");
						
						for(i=0;i<agents;i++){
							for(j=i;j<agents;j++){
								if(evolvingNetwork[i][j]==1){
									fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)networkAddedSimilarity[i][j]);
									edg = edg+1;
								}
							}
						}
						
						fprintf(fp,"</edges>\n");
						fprintf(fp,"</graph>\n");
						fprintf(fp,"</gexf>\n");
						
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
						
						
						/* sociarium file���� */
						sprintf(filename2,"nw_%d(sociarium).txt",writeOutEventTimes[k]); /*�t�@�C�������`*/
						
						/* �t�@�C���I�[�v�� */
						if ((fp = fopen(filename2, "wt")) == NULL) {
							fprintf(stderr, "�t�@�C���̃I�[�v���Ɏ��s���܂����D\n");
							exit(1);
						}
						
						/* �������� */
						fprintf(fp,"# Pajek�`���̃f�[�^ \n");
						fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
						fprintf(fp,"\n");
						fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
						fprintf(fp,"@title = The dissemination of culture \n");
						fprintf(fp,"@nondirected \n");
						fprintf(fp,"\n");
						
						fprintf(fp,"*Vertices %d \n",agents);
						
						fprintf(fp,"<nodes>\n");
						
						for(i=0;i<agents;i++){
							fprintf(fp,"%d \"agent_%d\" \n",i,i);
						}
						
						fprintf(fp,"*Arcs \n");
						
						for(i=0;i<agents;i++){
							for(j=i;j<agents;j++){
								if(evolvingNetwork[i][j]==1){
									fprintf(fp,"%d %d %0.1f \n",i,j,(double)networkAddedSimilarity[i][j]);
								}
							}
						}
						
						/* �t�@�C���N���[�Y */
						fclose(fp);
						/*----------------------------------------------------*/
					}
				}
			}
			
			/* sociarium file���� */	
			if(t%sociariumAnimationStepSize==0){
				if(t==0){
					/* �������� */
					fprintf(fp2,"@module = graph_creation_read_time_series_rect.dll \n");
					fprintf(fp2,"@title = cultural drift %f\n",cooperativityParameter);
					fprintf(fp2,"@delimiter =\t\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@time_format = Y.M.D # Y[ear], M[onth], D[ay], h[our], m[inute], s[econd] \n");
					fprintf(fp2,"@interval  = %d # �l�b�g���[�N���쐬���鎞���̊Ԋu\n",sociariumAnimationStepSize);
					fprintf(fp2,"@characteristic_time = 1 # �e�����ɂ���exp(-dt/@characteristic_time)��ώZ���C@threshold�ȉ��̃G�b�W���}�X�L���O \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@threshold = 0\n");
					fprintf(fp2,"@nondirected \n");
					fprintf(fp2,"\n");
					fprintf(fp2,"@node_texture = agent.png\n");
					fprintf(fp2,"\n");
					fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", 0, 0, 0, 0);	//���ꂪ�Ȃ��Ɛ������`�悪�ł��Ȃ�
					
					for(i=0;i<agents;i++){
						for(j=i;j<agents;j++){
							if(evolvingNetwork[i][j]){
								fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", t+sociariumAnimationStepSize, i, j, evolvingNetwork[i][j]);
							}
						}
					}
				}
				else{
					for(i=0;i<agents;i++){
						for(j=i;j<agents;j++){
							if(evolvingNetwork[i][j]){
								fprintf(fp2,"%d\tagent_%03d\tagent_%03d\t%d\n", t+sociariumAnimationStepSize, i, j, evolvingNetwork[i][j]);
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
		
		createMaxSizeNetwork();	//�ő�̓��̗אڍs��𐶐�
		
		island_dummy = countIslandNetwork();
		distance_dummy = calculateAverageVertexDistance();
		cluster_dummy = calculateAverageClusterCoefficient();
		
		if(readDataSet){
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
			
			if((island_number_av_array[index_mark] - sigmaConfficient*island_number_standard_deviation_array[index_mark] < island_dummy)&&(island_dummy< island_number_av_array[index_mark] + sigmaConfficient*island_number_standard_deviation_array[index_mark])&&(shortest_path_av_array[index_mark] - sigmaConfficient*shortest_path_standard_deviation_array[index_mark] < distance_dummy)&&(distance_dummy < shortest_path_av_array[index_mark] + sigmaConfficient*shortest_path_standard_deviation_array[index_mark])&&(cluster_number_av_array[index_mark] - sigmaConfficient*cluster_number_standard_deviation_array[index_mark] < cluster_dummy)&&(cluster_dummy < cluster_number_av_array[index_mark] + sigmaConfficient*cluster_number_standard_deviation_array[index_mark])){
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
		exit(1);
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
	
	for(i=0;i<agents;i++){
		for(j=0;j<agents;j++){
			if(maxSizeNetwork[i][j] == 1){
				fprintf(fp,"<node id=\"%d\" label=\"agent%d\" >\n",i,i);
				fprintf(fp,"<ns0:color b=\"%d\" g=\"%d\" r=\"%d\" />\n",agent[i].feature[0]*convertCultureToRgb,agent[i].feature[1]*convertCultureToRgb,agent[i].feature[2]*convertCultureToRgb);
				fprintf(fp,"</node>\n");
				
				break;
			}
		}
	}

	fprintf(fp,"</nodes>\n");
	fprintf(fp,"<edges>\n");

	for(i=0;i<agents;i++){
		for(j=i;j<agents;j++){
			if(maxSizeNetwork[i][j]==1){
				fprintf(fp,"<edge id=\"%d\" source=\"%d\" target=\"%d\" weight=\"%0.1f\"/>\n",edg,i,j,(double)networkAddedSimilarity[i][j]);
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
		exit(1);
	}
	
	/* �������� */
	fprintf(fp,"# Pajek�`���̃f�[�^ \n");
	fprintf(fp,"# �Q��: http://www.tp.umu.se/~rosvall/code.html \n");
	fprintf(fp,"\n");
	fprintf(fp,"@module = graph_creation_read_pajek.dll \n");
	fprintf(fp,"@title = The dissemination of culture (max size network)\n");
	fprintf(fp,"@nondirected \n");
	fprintf(fp,"\n");
	
	fprintf(fp,"*Vertices %d \n",connectedGraphInformation[1]);
	
	fprintf(fp,"<nodes>\n");
	
	for(i=0;i<agents;i++){
		fprintf(fp,"%d \"agent_%d\" \n",i,i);
	}
	
	fprintf(fp,"*Arcs \n");
	
	for(i=0;i<agents;i++){
		for(j=i;j<agents;j++){
			if(maxSizeNetwork[i][j]==1){
				fprintf(fp,"%d %d %0.1f \n",i,j,(double)networkAddedSimilarity[i][j]);
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
		exit(1);
	}
	
	/* �������� */
	fprintf(fp,"���̐� = %3d \n",island_dummy);
	fprintf(fp,"�ł��傫������agent�� = %3d \n",connectedGraphInformation[1]);
	fprintf(fp,"�����̐� = %3d \n",countCulture());
	fprintf(fp,"���ϒ��_�ԋ��� = %3f \n",distance_dummy);
	fprintf(fp,"���σN���X�^�[�W�� = %3f \n",cluster_dummy);
	fprintf(fp,"�u���b�W��(%%) = %f \n",countNetworkBridge()*100);
	
	fclose(fp);
	
	return 0;
}
