//C�����ж�̬�������ά���� malloc free
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//��̬�����ά����
template <typename T>
T** malloc_Array2D(int row, int col)
{
	int size = sizeof(T);
	int point_size = sizeof(T*);
	//�������ڴ棬����point_size * row��ʾ���row����ָ��
	T **arr = (T **)malloc(point_size * row + size * row * col);
	if (arr != NULL)
	{
		memset(arr, 0, point_size * row + size * row * col);
		T *head = (T*)((int)arr + point_size * row);
		while (row--)
			arr[row] = (T*)((int)head + row * col * size);
	}
	return (T**)arr;
}
//�ͷŶ�ά����
void free_Aarray2D(void **arr)
{
	if (arr != NULL)
		free(arr);
}
/*
int main()
{
	printf("  C�����ж�̬�������ά���� malloc free\n");
	printf(" -- by MoreWindows( http://blog.csdn.net/MoreWindows ) --\n\n");

	printf("����������(�Կո�ֿ�): ");
	int nRow, nCol;
	scanf("%d %d", &nRow, &nCol);

	//��̬���������Ķ�ά����
	int **p = malloc_Array2D<int>(nRow, nCol);

	//Ϊ��ά���鸳ֵ	
	int i, j;
	for (i = 0; i < nRow; i++)
	for (j = 0; j < nCol; j++)
		p[i][j] = i + j;

	//�����ά����	
	for (i = 0; i < nRow; i++)
	{
		for (j = 0; j < nCol; j++)
			printf("%4d ", p[i][j]);
		putchar('\n');
	}

	free_Aarray2D((void**)p);
	return 0;
}*/