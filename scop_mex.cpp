#include <math.h>
#include <mex.h>
#include <limits.h>
#include <matrix.h>
#include "BPnode.h"
#include "BPrag.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <queue>
#include "array2d.h"
using namespace std;

int indreg = 0;
int method = 2;

/*================ Compute the perimeter of one shape ================== */
float cal_primeter(BPnode *pnode, const _Pixels* pixels, int nrow, int ncol)
{
	int num_boun = 0;
	int i, j;
	unsigned char **fimage = malloc_Array2D<unsigned char>(nrow, ncol);
	//unsigned char fimage[256][256];
	//init the flag image zero
	for (i = 0; i<nrow; i++)
	{
		for (j = 0; j<ncol; j++)
		{
			fimage[i][j] = 0;
		}
	}
	//set the node area 1
	vector<int> pixelindex;
	pnode->get_node_pixels(pixelindex);
	vector<int>::iterator pindex = pixelindex.begin();
	while (pindex != pixelindex.end())
	{
		int x = pixels[(*pindex)].x; //col
		int y = pixels[(*pindex)].y; //row
		fimage[y][x] = 1;
		pindex++;
	}

	//calculate the boundary
	for (i = 0; i<nrow; i++)
	{
		for (j = 0; j<ncol - 1; j++)
		{
			if ((i == 0 && fimage[i][j] == 1) || ((i == nrow - 1) && fimage[i][j] == 1))
				num_boun++;
			else if ((j == 0 && fimage[i][j] == 1) || (j == ncol - 1 && fimage[i][j] == 1))
			{
				if (i>0 && i < nrow - 1)num_boun++;
			}
			else if (abs(fimage[i][j] - fimage[i][j + 1]) == 1)
				num_boun++;
		}
	}
	free_Aarray2D((void**)fimage);
	return (float)num_boun;
}
/*========= Compute the boundry rectangle of shapes ========*/
float cal_rectangle(vector<int> ppindex, const _Pixels* pixels, int ncol)
{
	int min_c = INT_MAX;
	int max_c = 0;
	int min_r = INT_MAX;
	int max_r = 0;
	int val = 0;
	vector<int>::iterator pindex = ppindex.begin();
	while (pindex != ppindex.end())
	{
		int x = pixels[(*pindex)].x; //col
		int y = pixels[(*pindex)].y; //row
		val = y*ncol + x;
		if (x < min_c) min_c = x;
		if (x > max_c) max_c = x;
		if (y < min_r) min_r = y;
		if (y > max_r) max_r = y;
		pindex++;
	}
	return (float)(max_c - min_c)*(max_r - min_r);
}
/*========= Compute the orientation, Elongation and Compactness
  (with new definition) of shapes ========*/
void shape_orilam(BPnode* root, vector<int> pixelind, _Pixels *imgin, float *out_ori, float *out_e, float *out_k)
{
	float size;
	float a11, a20, a02, x0, y0, sumx, sumy, lambda1, lambda2;
	int i;
	sumx = 0;
	sumy = 0;

	size = (float)pixelind.size();
	vector<int>::iterator pindex = pixelind.begin();
	while (pindex != pixelind.end())
	{
		sumx += (double)(imgin[(*pindex)].x);
		sumy += (double)(imgin[(*pindex)].y);
		pindex++;
	}

	x0 = sumx / (float)size;
	y0 = sumy / (float)size;

	a11 = 0;
	a20 = 0;
	a02 = 0;
	for (pindex = pixelind.begin(); pindex < pixelind.end(); pindex++)
	{
		a11 += (float)(imgin[(*pindex)].x - x0)*(imgin[(*pindex)].y - y0);
		a20 += (float)(imgin[(*pindex)].x - x0)*(imgin[(*pindex)].x - x0);
		a02 += (float)(imgin[(*pindex)].y - y0)*(imgin[(*pindex)].y - y0);
	}
	a11 = a11 / (float)pow(size, 1);
	a20 = a20 / (float)pow(size, 1) + 1.0 / 12.0;
	a02 = a02 / (float)pow(size, 1) + 1.0 / 12.0;

	/*Orientation*/
	*out_ori = 0.5*atan2(2 * a11, (a20 - a02)) + PI / 2;

	lambda1 = 0.5*(a02 + a20 + sqrt((a20 - a02)*(a20 - a02) + 4 * a11*a11));
	lambda2 = 0.5*(a02 + a20 - sqrt((a20 - a02)*(a20 - a02) + 4 * a11*a11));

	/*Elongation*/
	*out_e = lambda2 / lambda1;

	/*Compactness*/
	*out_k = (size) / (sqrt(lambda2*lambda1) * 4 * PI);
}
/*--------compute the parent scale area-------*/
float parent_scale_mean(BPnode* root, float size, int *nn)
{
	float smean = 0.0;
	int i, nump;
	BPnode* pShapeTemp;
	vector<int> pindex;
	smean = (float)size;
	pShapeTemp = root->parent;
	nump = 1;

	for (i = 0; i < (*nn); i++)
	{
		if (pShapeTemp == NULL)
			break;

		pShapeTemp->get_node_pixels(pindex);
		smean += (float)pindex.size();
		nump++;
		pindex.clear();
		pShapeTemp = pShapeTemp->parent;
	}
	smean /= nump;

	return smean;
}

/*====== Indexing the mn-order parent of the pShape ===*/
BPnode* m_order_parent(BPnode *pnode, int *mn)
{
	int t;
	BPnode *pnodeTemp;
	pnodeTemp = pnode;
	for (t = 0; t < (*mn); t++)
	{
		if (pnodeTemp->parent == NULL)
			break;

		pnodeTemp = pnodeTemp->parent;
	}

	return pnodeTemp;
}
/*--------compute all the features-------*/
void calculate_feature_V1(vector<BPnode*> vnodes, _Pixels *imgin, int min_shape, int max_shape,int nrow,int ncol)
{
	float  fo, fe, fc;
	float val, size, sumx[3],sumxx[3], area;
    vector<int>  sboundry;
	int i, j, npscale = 5;                               ////////////////////////////////////
	float perimeter = 0.0;
	for (i = 0; i < vnodes.size(); i++)
	{
		vector<int> pixelind;
		vnodes[i]->get_node_pixels(pixelind);
		vnodes[i]->cal_node_homo(pixelind, imgin);
		size = (float)pixelind.size();
        if(size<min_shape||size>max_shape)
            continue;

		for (int tt = 0; tt < 3; tt++)
		{
			sumx[tt] = 0.0;
            sumxx[tt]= 0.0;
			vector<int>::iterator pindex = pixelind.begin();
			while (pindex != pixelind.end())
			{
				val = imgin[(*pindex)].value(tt, 0)-vnodes[i]->covmat(tt,0);
				sumx[tt] += val*val;  
                sumxx[tt]+=pow(val,(float)3.0);
				pindex++;
			}
			/* contrast*/
			vnodes[i]->pnodeinfo.features[2 + tt] = sqrt(sumx[tt] / (float)size);
            /*skewness*/
            vnodes[i]->pnodeinfo.features[13 +tt] = pow(sumxx[tt]/(float)size,(float)1.0/3);
		}
		/*====== Compute the attributes of the shape ===*/
		vnodes[i]->pnodeinfo.features[0] = size/ cal_rectangle(pixelind,imgin,256);       //Ô¤Áô¼ÆËã4pi*area/primeter^2
        if(vnodes[i]->pnodeinfo.features[0]>1.0)
            vnodes[i]->pnodeinfo.features[0]=1.0;
		/* scale ratio*/
		vnodes[i]->pnodeinfo.features[1] = size / parent_scale_mean(vnodes[i], size, &npscale);

		/* Orientation, Elongation, Compactness2 */
		shape_orilam(vnodes[i], pixelind, imgin, &fo, &fe, &fc);
		vnodes[i]->pnodeinfo.features[5] = fo / PI;
		vnodes[i]->pnodeinfo.features[6] = fe;
		vnodes[i]->pnodeinfo.features[7] = (fc - 0.3) / 0.7;
		/*homogeneity*/
		vnodes[i]->pnodeinfo.features[8] = vnodes[i]->homogeneity;  // ((float) 4* PI* size)/pow(peri_shape(vnodes[i]),2);
		//if( vnodes[i]->pnodeinfo.features[6] > 1 )
		//	vnodes[i]->pnodeinfo.features[6] = 1;
		/* scale*/
		vnodes[i]->pnodeinfo.features[9] = (log10((float)size) - log10((float)min_shape)) / (log10((float)max_shape) - log10((float)min_shape));
		/*gray level*/
		vnodes[i]->pnodeinfo.features[10] = vnodes[i]->covmat(0, 0) ;
		vnodes[i]->pnodeinfo.features[11] = vnodes[i]->covmat(1, 0) ;
		vnodes[i]->pnodeinfo.features[12] = vnodes[i]->covmat(2, 0) ;
		//vnodes[i]->pnodeinfo.index=vnodes[i]->regionID;
		sboundry.clear();
		boundry(vnodes[i], vnodes, sboundry);
		perimeter = (float)sboundry.size();
		vnodes[i]->pnodeinfo.features[16] = ((float)4 * PI* size) / pow(perimeter,2);       //Ô¤Áô¼ÆËã4pi*area/primeter^2
		if (vnodes[i]->pnodeinfo.features[16]>1.0)
			vnodes[i]->pnodeinfo.features[16] = 1.0;
	}

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int x, y;           //the coordinate (x,y)
	int k, i, j, j2, j3, jc, n_elements, n_structs;
	int numfeature, numfeature2, numfeature3, numfeaturec;
	int max_shape, min_shape, mn_order,mn_order1;
	double Y_double, Cb_double, Cr_double, *Y_data, *Cb_data, *Cr_data;
	float Y_float, Cb_float, Cr_float;
	double	*im_h, *im_w, *temp_data, *nu, *minshape, *maxshape;
	int *featout1_ind, *featout2_ind, *featout3_ind, *featoutc_ind;
	BPnode  *pShape2 = new BPnode;
	BPnode  *pShape3 = new BPnode;
    BPnode  *pShape4 = new BPnode;
	BPnode 	*pShape_c1 = new BPnode;
	_Pixels *pixels;
	int MAX_NUM = 30;
	mxArray  *temp_matrix;
	/* new defination */
	const char *field_names[5];
	char field0[10] = "m_attri";
	char field1[10] = "shapes";
	char field2[10] = "pattern2";
	char field3[10] = "pattern3";
	char field4[10] = "pattern4";

	field_names[0] = field0;
	field_names[1] = field1;
	field_names[2] = field2;
	field_names[3] = field3;
	field_names[4] = field4;

	//pTree = mw_new_shapes();
	numfeature = 0;

	/* check out input parameters */
	if (nrhs != 8)
		mexErrMsgTxt("Five input required.");
	else if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");

	im_h = mxGetPr(prhs[3]);   /* get in image width and height */
	im_w = mxGetPr(prhs[4]);
	nu = mxGetPr(prhs[5]);/*get the m-order*/
	minshape = mxGetPr(prhs[6]);
	maxshape = mxGetPr(prhs[7]);
	/* guarantee: n_elements == im_w * im_h */
	n_elements = mxGetNumberOfElements(prhs[0]);
	if (n_elements != (*im_w)*(*im_h))
		mexErrMsgTxt("Length of image data does not match image width and height.");
	pixels = new _Pixels[n_elements];
	/* initialize image and tree */
	//pImageInput.create((int)(*im_h), (int)(*im_w),uchar);
	//pImagetemp.create((int)(*im_h), (int)(*im_w),uchar);
	Y_data = mxGetPr(prhs[0]);   /*get in image data(double 2 float)  */
	Cb_data = mxGetPr(prhs[1]);
	Cr_data = mxGetPr(prhs[2]);
	for (i = 0; i < n_elements; i++)     /*data type transfermation */
	{
		y = i / (int)(*im_h);    //row
		x = i % (int)(*im_h);    //col
		Y_double = *(Y_data + i);
		Cb_double = *(Cb_data + i);
		Cr_double = *(Cr_data + i);
		Y_float = (float)Y_double;
		Cb_float = (float)Cb_double;
		Cr_float = (float)Cr_double;
		pixels[i].x = x;
		pixels[i].y = y;
		pixels[i].init(x, y, Y_float, Cb_float, Cr_float);
	}
	//mw_alloc_shapes(pTree, (int)(*im_h), (int)(*im_w), pImageInput->gray[0]);
	vector<BPnode*> vnodes;
	priority_queue<BPrag*, vector<BPrag*>, cmprag> RAGS;
	vnodes.resize(2 * n_elements - 1);

	indreg = 0;
	for (int i = 0; i < n_elements; i++)
	{
		BPnode * pnode = new BPnode(pixels, i, (int)(*im_h), (int)(*im_w), indreg++);
		vnodes[i] = pnode;
	}
	bp_init_edges(vnodes, (int)(*im_h), (int)(*im_w), RAGS);     //ÁÚÓòÁ¬½ÓÍ¼
    bp_buildtree(vnodes, RAGS, pixels);
	/* flst operation */
	//flst(NULL, pImagetemp, pTree);
	//flst_pixels(pTree);

	/* allocate space for sorted shapes */
	n_structs = 2 * n_elements;

	/* calculate the features of each shape*/

	featout1_ind = (int*)malloc(sizeof(int)*(n_structs));
	featout2_ind = (int*)malloc(sizeof(int)*(2 * n_structs));
	featout3_ind = (int*)malloc(sizeof(int)*(3 * n_structs));
	featoutc_ind = (int*)malloc(sizeof(int)*(3*n_structs));

	max_shape = *maxshape;
	min_shape = *minshape;
	mn_order = *nu;
	calculate_feature_V1(vnodes, pixels, min_shape, max_shape, (int)(*im_h), (int)(*im_w));

	/*Process Pattern 1*/
	for (i = 0, j = 0; i < vnodes.size(); i++)
	{

		if (vnodes[i] != NULL && vnodes[i]->parent != NULL && vnodes[i]->pixelnum >= min_shape && vnodes[i]->pixelnum <= max_shape
			&& vnodes[i]->pnodeinfo.features[7] >= 0 && vnodes[i]->pnodeinfo.features[2] >= 0 && vnodes[i]->pnodeinfo.features[3] >= 0&&vnodes[i]->pnodeinfo.features[4] >= 0
			&& vnodes[i]->pnodeinfo.features[2] <= 1 && vnodes[i]->pnodeinfo.features[3] <= 1 && vnodes[i]->pnodeinfo.features[4] <= 1)
		{
			vnodes[i]->pnodeinfo.index = j++;
		}

	}
	/*Process Pattern 2, 3, 4*/

	for (i = 0, j2 = 0, j3 = 0, jc = 0; i < vnodes.size(); i++)
	{
		if (vnodes[i]->pnodeinfo.index >= 0)
		{
			/*Pattern 4*/
			//mn_order1 = mn_order - 1;
			pShape2 = vnodes[i]->parent;//m_order_parent(vnodes[i], &mn_order1);
			if (pShape2 != NULL && pShape2->pnodeinfo.index >= 0)
			{
				if (pShape2->left->regionID != vnodes[i]->regionID)
					pShape_c1 = pShape2->left;
				else
					pShape_c1 = pShape2->right;
				if (pShape_c1 != NULL && pShape_c1->pnodeinfo.index >= 0)
				{
					featoutc_ind[jc++] = vnodes[i]->pnodeinfo.index;
					featoutc_ind[jc++] = pShape2->pnodeinfo.index;
					featoutc_ind[jc++] = pShape_c1->pnodeinfo.index;
				}
				
				//pShape_c1 = pShape_c1->next_sibling;
			}
			/*Pattern 2*/
			pShape2 = m_order_parent(vnodes[i], &mn_order);
			if (pShape2->pnodeinfo.index >= 0)
			{
				featout2_ind[j2++] = vnodes[i]->pnodeinfo.index;
				featout2_ind[j2++] = pShape2->pnodeinfo.index;
			}
            /*Pattern 3*/
			mn_order1 = mn_order;
			pShape3 = m_order_parent(pShape2, &mn_order1);
			//mn_order1 = mn_order + 2;
            //pShape4 = m_order_parent(vnodes[i], &mn_order1);
		    if (pShape3->pnodeinfo.index >= 0&&pShape3!=NULL)
				{
					featout3_ind[j3++] = vnodes[i]->pnodeinfo.index;
					featout3_ind[j3++] = pShape2->pnodeinfo.index;
                    featout3_ind[j3++] = pShape3->pnodeinfo.index;
				}
		}
	}
	numfeature = j;
	numfeature2 = j2;
	numfeature3 = j3;
	numfeaturec = jc;

	/*Output Configure*/
	plhs[0] = mxCreateStructMatrix(1, 1, 5, field_names);

	temp_matrix = mxCreateDoubleMatrix(1, 1, mxREAL);
	temp_data = mxGetPr(temp_matrix);
	*temp_data = NFEAT;
	mxSetField(plhs[0], 0, field_names[0], temp_matrix); /*number of Shape attributes*/
	
	temp_matrix = mxCreateDoubleMatrix(1, numfeature*NFEAT, mxREAL);
	temp_data = mxGetPr(temp_matrix);
	k = 0;
	for (i = 0; i < n_structs - 1; i++)
	{
		if (vnodes[i]->pnodeinfo.index < 0)
			continue;
		for (j = 0; j < NFEAT; j++)
			temp_data[k*NFEAT + j] = vnodes[i]->pnodeinfo.features[j];
		k++;
	}
	mxSetField(plhs[0], 0, field_names[1], temp_matrix); /*Shape attributes, Pattern 1*/
	
	temp_matrix = mxCreateDoubleMatrix(1, numfeature2, mxREAL);
	temp_data = mxGetPr(temp_matrix);
	for (i = 0; i < numfeature2; i++)
	{
		temp_data[i] = featout2_ind[i];
	}
	mxSetField(plhs[0], 0, field_names[2], temp_matrix); /*Pattern 2*/

	temp_matrix = mxCreateDoubleMatrix(1, numfeature3, mxREAL);
	temp_data = mxGetPr(temp_matrix);
	for (i = 0; i < numfeature3; i++)
	{
		temp_data[i] = featout3_ind[i];
	}
	mxSetField(plhs[0], 0, field_names[3], temp_matrix); /*Pattern 3*/

	temp_matrix = mxCreateDoubleMatrix(1, numfeaturec, mxREAL);
	temp_data = mxGetPr(temp_matrix);
	for (i = 0; i < numfeaturec; i++)
	{
		temp_data[i] = featoutc_ind[i];
	}
	mxSetField(plhs[0], 0, field_names[4], temp_matrix); /*Pattern 4*/


	/* free pointers */
	delete[] pixels;
    for(i=0;i<vnodes.size();i++)
    {
        delete vnodes[i];
    }
	vnodes.clear();                             //7.11ÐÞ¸Ä
	
	
	free(featout1_ind);
	free(featout2_ind);
	free(featout3_ind);
	free(featoutc_ind);
	//delete pShape2;
	//delete pShape3;
	//delete pShape_c1;
	//free(Y_data);
	//free(Cb_data);
	//free(Cr_data);
	//free(pShape2);
	//free(pShape3);
	//free(pShape_c1);
	//free(pTree);

}