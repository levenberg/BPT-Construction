#pragma once
//two nodes comparison
class BPrag
{
public:
	BPrag(void);
	BPrag(int r1,int r2,float dist):r1(r1),r2(r2),dist(dist){}
	~BPrag(void);
public:
	int r1;
	int r2;
	float dist;
};


struct cmprag
{
	bool operator()( BPrag *rag1, BPrag *rag2)
	{
		return rag1->dist > rag2->dist;  //the small value first
	}
};
