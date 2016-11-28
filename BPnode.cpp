#include "BPnode.h"
//node handle cpp
extern int indreg;
extern int method;


BPnode::BPnode(void)
{
	left = NULL;
	right = NULL;
	parent = NULL;
	flag = 1;
	homogeneity = -1;
	for (int i=0;i<NFEAT;i++)
	{
		pnodeinfo.features[i]=0.0;
		pnodeinfo.index=-1;
	}
 
}


BPnode::~BPnode(void)
{
}

BPnode::BPnode(const _Pixels *pixels, int i, int row, int col, int index)
{
	regionID = index;
	pixelnum = 1;
	flag = 1; //表示未使用
	homogeneity = -1;
	left = NULL;
	right = NULL;
	parent = NULL;
	covmat = pixels[i].value;
	for (int i = 0; i<NFEAT; i++)
	{
		pnodeinfo.features[i] = 0.0;
		pnodeinfo.index = -1;
	}
}
void BPnode::get_node_pixels( vector<int> & pixelind)
{
	if( left == NULL) //由于是二叉合并树，因此左右子树要么均存在，要么均不存在
		pixelind.push_back( regionID );
	else
	{
		left->get_node_pixels( pixelind);
		right->get_node_pixels( pixelind);
	}
	
}

void BPnode::cal_node_homo(vector<int> pixelind, const _Pixels* pixels)
{
	int i;
	Vector3f sum, temp;
	float sum1 =0.0, sum2 = 0.0;
	//计算 covmat 的范数
	for( int j = 0; j<3; j++)
				sum2 = sum2 + pow(covmat(j,0),2) ;

	vector<int>::iterator pindex = pixelind.begin();
	while( pindex != pixelind.end() )
	{
		temp = (pixels[(*pindex)].value - covmat);
		for( int j = 0; j<3; j++)
				sum1 = sum1 + pow(temp(j,0),2) ;
		pindex++;
	}
	homogeneity = sum1/sum2/pixelnum;      //异质度，每个节点与其所在树的根节点的差异性均值
}


void BPnode::cal_covmatrix()
{
	Vector3f sum;
	sum = (left->covmat)*left->pixelnum + (right->covmat)*right->pixelnum;
	covmat = sum/(float)pixelnum;
}

float BPnode::cal_node_dist(BPnode * region2, int method)
{
	float D,tmpD;
	vector<int> pix1,pix2;
	Vector3f avg;
	if( method < 0 || method > 5)
		return 0.0;
	switch (method)
	{
	case 1: //三通道差值平方（欧式距离ED）
		   D=pow(covmat(0,0)-region2->covmat(0,0),2)+pow(covmat(1,0)-region2->covmat(1,0),2)+pow(covmat(2,0)-region2->covmat(2,0),2);
		break;
	case 2: //DEDM加权欧式距离
		tmpD=pow(covmat(0,0)-region2->covmat(0,0),2)+pow(covmat(1,0)-region2->covmat(1,0),2)+pow(covmat(2,0)-region2->covmat(2,0),2);
		if (pixelnum>region2->pixelnum)
			D=tmpD*region2->pixelnum;
		else if(pixelnum==region2->pixelnum)
			D=tmpD*(pixelnum+region2->pixelnum)/2.0;
		else
			D=tmpD*pixelnum;
		break;
	case 3: //Ddn
		avg=(pixelnum*covmat+region2->pixelnum*region2->covmat)/(pixelnum+region2->pixelnum);
		get_node_pixels( pix1);
		region2->get_node_pixels(pix2);
		for( int i = 0 ; i< pix1.size(); i++)
		{
			
		}
		D=0.0;
		break;
	default:
		break;
	}

	if( D != D)
		D = 10000.0;
	return D; //之前算的是距离，求反后就将距离转化为相似度了。。。
}

BPnode * bp_mergenode(vector<BPnode*> &vnodes, int r1, int r2, const _Pixels *pixels )
{
	int i = 0;
	BPnode* node1 = vnodes[r1];
	BPnode *node2 = vnodes[r2];
	BPnode * newnode = new BPnode;
	newnode->regionID = indreg;
	newnode->pixelnum = node1->pixelnum + node2->pixelnum;
	newnode->left = node1;
	newnode->right = node2;
	newnode->parent = NULL;
	set<int>::iterator pnode = node1->adjacents.begin();
	while( pnode != node1->adjacents.end() )
	{
		if( *pnode != node2->regionID )
		{
			newnode->adjacents.insert( *pnode);
			vnodes[*pnode]->adjacents.erase( node1->regionID);
			vnodes[*pnode]->adjacents.insert(indreg);
		}
		pnode++;
	}
	pnode = node2->adjacents.begin();
	while( pnode != node2->adjacents.end() )
	{
		if( *pnode != node1->regionID )
		{
			newnode->adjacents.insert( *pnode);
			vnodes[*pnode]->adjacents.erase( node2->regionID);
			vnodes[*pnode]->adjacents.insert(indreg);
		}
		pnode++;
	}
	newnode->cal_covmatrix();
	vnodes[indreg] = newnode;


	node1->flag = 0;
	node1->parent = newnode;
	node2->flag = 0;
	node2->parent = newnode;
	//将 node1 和 node 2 中所有相邻节点的邻域表更新

	indreg++;
	return newnode;
}

BPrag * bp_create_rag( BPnode* pnode1, BPnode* pnode2)
{
	float dist = pnode1->cal_node_dist(pnode2, method);
	//printf("the distance=%f\n",dist);
	BPrag * newrag = new BPrag(pnode1->regionID, pnode2->regionID, dist);
	return newrag;
}

//初始化情况下，按照像素位置增加其邻域，并将邻域对加入 RAGS表中
void bp_init_edges(vector<BPnode*> & vnodes, int row, int col, priority_queue<BPrag*, vector<BPrag*>, cmprag> &RAGS)
{
	int i = 0;
	int p3reg =0, p4reg = col;//指向右边相邻和下边相邻的点
	int NPIXEL = row*col;
	for( i = 0; i< vnodes.size(); i++)
	{
		p3reg = i+1;
		p4reg = i + col;
		if( p3reg < NPIXEL && p3reg % col != 0 )
		{
			vnodes[i]->adjacents.insert(p3reg);
			vnodes[p3reg]->adjacents.insert(i);
            vnodes[i]->pixjacents.insert(p3reg);
			vnodes[p3reg]->pixjacents.insert(i);
			BPrag * newrag = bp_create_rag( vnodes[i], vnodes[p3reg]);
			RAGS.push( newrag);
		}
		if( p4reg < NPIXEL )
		{
			vnodes[i]->adjacents.insert(p4reg);
			vnodes[p4reg]->adjacents.insert(i);
            vnodes[i]->pixjacents.insert(p4reg);
			vnodes[p4reg]->pixjacents.insert(i);
			BPrag * newrag = bp_create_rag( vnodes[i], vnodes[p4reg]);
			RAGS.push( newrag);
		}
	}
}

void bp_add_edges(BPnode * node, vector<BPnode*> & vnodes, priority_queue<BPrag*, vector<BPrag*>, cmprag> &RAGS)
{
	set<int>::iterator padj = node->adjacents.begin();
	while( padj != node->adjacents.end() )
	{
		BPrag * newrag = bp_create_rag( vnodes[*padj], node);
		RAGS.push(newrag);
		padj ++;
	}
}

void bp_buildtree(vector<BPnode*> & vnodes, priority_queue<BPrag*, vector<BPrag*>, cmprag> &RAGS, const _Pixels *pixels )
{
	BPrag * cur_rag;
	BPnode *pnode;
	while( RAGS.size() != 0 )
	{
		cur_rag = RAGS.top(); 
		RAGS.pop();
		if( vnodes[cur_rag->r1]->flag == 0 || vnodes[cur_rag->r2]->flag == 0 )      //合并的前提
		{
			delete cur_rag;
			continue;
		}
		pnode = bp_mergenode(vnodes, cur_rag->r1, cur_rag->r2, pixels);
		delete cur_rag;
		bp_add_edges( pnode, vnodes, RAGS);
	}
	//delete pnode;
}

void bp_prun_tree(BPnode* root,_Pixels *pixels, float thr,  list<BPnode*> & nreg)
{
	if( root == NULL )
		return;
	vector<int> pixelind;
	root->get_node_pixels( pixelind);
	root->cal_node_homo(pixelind, pixels);
	if( root->homogeneity < thr )
	{
		nreg.insert( nreg.end(), root );
		return;
	}
	bp_prun_tree( root->left, pixels, thr, nreg);
	bp_prun_tree( root->right, pixels, thr, nreg);

}

bool bp_greater(const BPnode *node1, const BPnode *node2)
{
	if( node1->homogeneity < node2->homogeneity )
		return true;
	else
		return false;
}
void boundry(BPnode *pnode, vector<BPnode*> & vnodes, vector<int> &bound)
{
	bool flag = false;
	vector<int> pixelindex;
	pnode->get_node_pixels(pixelindex);
	vector<int>::iterator pindex = pixelindex.begin();
	while (pindex != pixelindex.end())
	{
		if (vnodes[(*pindex)]->pixjacents.size()<4)    //the boundry of image is the boundry of a shape
		{
			bound.push_back(*pindex);
			pindex++;
			continue;
		}
		flag = false;
		set<int>::iterator padj = vnodes[(*pindex)]->pixjacents.begin();
		while (padj != vnodes[(*pindex)]->pixjacents.end() && flag == false)
		{
			if (find(pixelindex.begin(), pixelindex.end(), *padj) == pixelindex.end())
			{
				bound.push_back(*pindex);
				flag = true;
			}
			padj++;
		}
		pindex++;
	}
}