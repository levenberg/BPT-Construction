//C语言中动态的申请二维数组 malloc free
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//动态申请二维数组
template <typename T>
T** malloc_Array2D(int row, int col)
{
	int size = sizeof(T);
	int point_size = sizeof(T*);
	//先申请内存，其中point_size * row表示存放row个行指针
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
//释放二维数组
void free_Aarray2D(void **arr)
{
	if (arr != NULL)
		free(arr);
}
/*
int main()
{
	printf("  C语言中动态的申请二维数组 malloc free\n");
	printf(" -- by MoreWindows( http://blog.csdn.net/MoreWindows ) --\n\n");

	printf("请输入行列(以空格分开): ");
	int nRow, nCol;
	scanf("%d %d", &nRow, &nCol);

	//动态申请连续的二维数组
	int **p = malloc_Array2D<int>(nRow, nCol);

	//为二维数组赋值	
	int i, j;
	for (i = 0; i < nRow; i++)
	for (j = 0; j < nCol; j++)
		p[i][j] = i + j;

	//输出二维数组	
	for (i = 0; i < nRow; i++)
	{
		for (j = 0; j < nCol; j++)
			printf("%4d ", p[i][j]);
		putchar('\n');
	}

	free_Aarray2D((void**)p);
	return 0;
}*/