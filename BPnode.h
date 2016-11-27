#pragma once
#include <math.h>
#include <vector>
#include <complex>
#include <queue>
#include <list>
#include <hash_set>
#include <map>
#include <set>
#include "eigen/Eigen/Dense"
#include "BPrag.h"
using namespace Eigen;
using namespace std;
#define NFEAT 17
#define PI 3.14159
struct _Pixels
{
	int x;
	int y;
	Vector3f value;    //(RGB)
	//Matrix3cf value; 
	void init( int x, int y, float v1, float v2, float v3)
	{
		Vector3f temp(v1,v2,v3);
		this->x = x;
		this->y = y;
		value = temp;   //*temp.adjoint();
	}
};
typedef struct Info{
	int index;
	float features[NFEAT];
}Info;
class BPnode
{
public:
	BPnode(void);
	~BPnode(void);
	BPnode(const _Pixels *pixels,int i, int row, int col, int index);
public:
	set<int> adjacents; //邻域表
    set<int> pixjacents; //像素邻域表
	int pixelnum;
	Vector3f covmat; //平均协方差矩阵，用于计算相似度
	float homogeneity;
	BPnode * parent, *left, *right;
	Info pnodeinfo;
	int regionID;
	int flag; //表示区域是否有效，即已被合并的区域 flag = 0;
public:
	void cal_covmatrix();
	float cal_node_dist(BPnode * temp, int method);
	
    void cal_node_homo(vector<int> pixelind, const _Pixels* pixels);
	void get_node_pixels( vector<int> & pixels);
};

BPnode * bp_mergenode(vector<BPnode*> &vnodes, int r1, int r2, const _Pixels *pixels );

//初始化情况下，按照像素位置增加其邻域，并将邻域对加入 RAGS表中
void bp_init_edges(vector<BPnode*> & vnodes, int row, int col, priority_queue<BPrag*, vector<BPrag*>, cmprag> &RAGS);
//对每个节点node, 按照其邻域表增加对应的边缘到RAGS表中
void bp_add_edges(BPnode * node,vector<BPnode*> & vnodes,  priority_queue<BPrag*, vector<BPrag*>, cmprag> &RAGS);
BPrag * bp_create_rag( BPnode*, BPnode*);
void bp_buildtree(vector<BPnode*> & vnodes, priority_queue<BPrag*, vector<BPrag*>, cmprag> & RAGS, const _Pixels *pixels );
void boundry(BPnode *pnode, vector<BPnode*> & vnodes, vector<int> &bound);
void bp_prun_tree(BPnode* tree, _Pixels* pixels,  float thr, list<BPnode*> & nreg);
bool bp_greater(const BPnode *node1, const BPnode *node2);