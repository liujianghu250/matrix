#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define MAXN 100

//定义分数结构体
typedef struct
{
    int x,y;//分子。分母

} fraction;
//fraction:分数


//定义以分数为元素的矩阵结构
typedef struct Matrix
{
    //declare Matrix;
    int row,col;   //row:行。col：列。
    fraction value[MAXN][MAXN];
}Matrix;
//matrix:矩阵

// 定义以浮点数为元素的矩阵
typedef struct Matrix_x
{
    int col,row;
    double value[MAXN][MAXN];
}Matrix_x;

void print_fraction(fraction a);//打印分数
int print_matrix(Matrix * mat);//打印分数为元素的矩阵
void print_matrix_x(Matrix_x * mat);
int gcd(int x,int y);//GCD 阿基米德最大公约数法
void reduction(int* x,int* y);//约分
fraction add(fraction* a,fraction*b);//equals "+"
fraction add(fraction* a,fraction*b);//equals "-"
fraction mul(fraction*a, fraction*b);//equals "*"
fraction divide(fraction*a, fraction*b);//equals "/"
Matrix_x  *add_decimal(Matrix_x * mat1,Matrix_x * mat2);//矩阵加法
Matrix  *add_(Matrix * mat1,Matrix * mat2);//分数形式
Matrix_x * sub_decimal(Matrix_x * mat1,Matrix_x * mat2);//矩阵减法，浮点数元素
Matrix * sub_(Matrix * mat1,Matrix * mat2);//矩阵减法，分数元素
Matrix_x * mul_decimal(Matrix_x * mat1,Matrix_x * mat2);//矩阵的乘法，浮点数元素
Matrix * mul_(Matrix * mat1,Matrix * mat2);//矩阵乘法，分数元素
Matrix_x * transpose_decimal(Matrix_x * mat);//矩阵转置，浮点数元素
Matrix * transpose_(Matrix * mat);//矩阵转置，分数元素
Matrix * read(Matrix *  mat);//input Matrix;
Matrix_x * read_decimal(Matrix_x * mat);//读入浮点数为元素的矩阵
void swap_row(Matrix * mat,int row1,int row2);//分数行互换
void swap_row_decimal(Matrix_x * mat,int row1,int row2);//小数行互换
Matrix * cha_into_dia(Matrix * mat);//化成矩阵最简型分数
Matrix_x * cha_into_dia_dicimal(Matrix_x * mat);//小数
int judge_rs(Matrix * mat);//判断分数型矩阵的行数和列数是否相等
int judge_rs_x(Matrix_x * mat);//判断浮点数型矩阵行数与列数是否相等
int judge_symmetric(Matrix * mat);//判断是否是对称矩阵
int judge_symmetric_x(Matrix_x * mat); //判断浮点型矩阵是否为实对称矩阵
void jacbi(Matrix_x *mat ,Matrix_x * dblVect);//计算实对称矩阵的特征值和特征向量 雅克比迭代法
Matrix_x * scalar_decimal(Matrix_x *a,double k);//矩阵数乘
Matrix * scalar(Matrix *a,fraction k);//分数
double  * scalar_(double *a,double b,int n,double *c);//向量数乘
fraction dot_product(fraction *a,fraction * b,int n);//向量内积
double dot_product_decimal(double *a,double * b,int n);
double lenth(double * a,int n);//向量长度
double * add_vector(double * a,double * b,int n,double *c);//向量加法
void schmidt(double *a,int n);//施密特正交化法
fraction  cal_det(Matrix * mat);//计算行列式
double  cal_det_decimal(Matrix_x * mat);//计算行列式
int is_singular(Matrix * mat);//判断是否是奇异矩阵
int is_singular_decimal(Matrix_x * mat);//小数形式
void cal_quad(Matrix_x *mat);//二次型正惯性指数与负惯性指数
Matrix * get_inverse(Matrix * mat);//求逆
Matrix_x * get_inverse_decimal(Matrix_x * mat);
int  rank(Matrix * mat);//求秩
int  rank_decimal(Matrix_x * mat);//求秩
int  syst(Matrix * mat);//求齐次方程的基础解系
int  syst_decimal(Matrix_x * mat);
void linear_ab(Matrix *mat);//齐次线性方程组
void linear_ab_decimal(Matrix_x *mat);
void linear_nor(Matrix *mat);//求解非齐次线性方程组
void linear_nor_decimal(Matrix_x *mat);//线性方程
int  input(double * a);
int  print_tips();//输入提示
int orders(int *n,int *m,Matrix_x *A[],Matrix *A_[],char *name,char*name_);//命令模式
int operations(int n);

int main(int argc, char const *argv[])
{

    while(1)
    {
        int n;
        n = print_tips();
        if(n == 0)
        {
            printf("welcome to orders mode !\n");
            int x=0,y=0;
            int *namen,*namem;
            namen = &x;
            namem = &y;
            Matrix * A_[MAXN] ;
            Matrix_x * A[MAXN];
            char *name = (char * )malloc(sizeof(100));
            char *name_= (char * )malloc(sizeof(100));
            while (1) {

                orders(namen,namem,A,A_,name,name_);
            }
        }
        else
        {
            operations(n);
        }
    }
    system("PAUSE");

    return 0;
}

//打印分数
void print_fraction(fraction a)
{
    if (a.y==1) //如果分母为一
    {
        printf("%d",a.x );//只打印分子
    }

    else 			//否则分母不为一
		printf("%d/%d",a.x,a.y);//分子分母一起打印
}


//打印分数为元素的矩阵
int print_matrix(Matrix * mat)
{
    int i,j;//用作循环用的变量

    printf("MAtrix:\n");//for beauty;
    for (i = 0; i < mat->row; i++)//对行进行循环
	{
        for ( j = 0; j < mat->col; j++)//列循环
		{
        //this part can change into a function;
        if(mat->value[i][j].x>=0)
			printf(" ");//打印空格使得打印出来的矩阵整齐

        print_fraction(mat->value[i][j]);

        if(mat->value[i][j].y==1)//如果分母为一
			printf("    ");

        else
			printf("  ");
        }

        printf("\n");
    }
    return 0;
}


void print_matrix_x(Matrix_x * mat)
{
    int i,j;
    printf("\t");//使输出呈表格形

    for(j=0;j<mat->col;j++)
    {
        printf("第%d列\t\t   ",j+1);//控制输出格式
    }

    printf("\n");

    for(i=0;i<mat->row;i++)
	{
        for(j=0;j<mat->col;j++)
		{
            if(j==0)
			{
                printf("第%d行\t",i+1 );//行首打印当前行数
            }

            printf("%15lf ",mat->value[i][j] );

            if(j==mat->col-1)
			{
                printf("\n");//行尾换行
            }
        }
    }
}


//GCD 阿基米德最大公约数法
int gcd(int x,int y)
{
    return y==0?x:gcd(y,x%y);
}


//约分
void reduction(int* x,int* y)
{
    int cnt=( *x < 0)+(*y < 0);//the number of (<0)
    if (*x < 0)
	{
        *x *= -1;
    }

    if (*y<0)
	{
        *y *= -1;
    }

    int GCD=gcd(*x,*y);//分子分母最大公约数
    *x /= GCD;
    *y /= GCD;

    if (cnt &1)
	{//even &1=0,odd&1=1
        *x *=-1;
    }
}


//equals "+"
fraction add(fraction* a,fraction*b)
{
    fraction result=//分数相加
    {
        a->x*b->y+a->y*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//约分

    return result;
}


//equals "-"
fraction dec(fraction* a,fraction*b)
{
    fraction result={//分数相减
        a->x*b->y-a->y*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//约分

    return result;
}


//equals "*"
fraction mul(fraction*a, fraction*b)
{
    fraction result={//分数相乘
        a->x*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//约分

    return result;
}


//equals "/"
fraction divide(fraction*a, fraction*b)
{
    fraction result ={//分数相除
        a->x*b->y,
        a->y*b->x
    };

    reduction(&result.x,&result.y);//约分

    return result;
}


//矩阵加法
Matrix_x  *add_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j;
    Matrix_x * mat3;
    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//分配内存
    //必须赋值，否则无法输出
    mat3->col=mat1->col;
    mat3->row=mat1->row;

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat3->value[i][j]=mat1->value[i][j] + mat2->value[i][j];//对应元素相加
        }
    }

    return mat3;
}

//分数形式
Matrix  *add_(Matrix * mat1,Matrix * mat2)
{
    int i,j;
    Matrix* mat;
    mat=(Matrix *)malloc(sizeof(Matrix));//分配内存
    //行列赋值
    mat->col=mat1->col;
    mat->row=mat1->row;
    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat->value[i][j] =add(&mat1->value[i][j],&mat2->value[i][j]);//相加
        }
    }

    return mat;
}

//矩阵减法，浮点数元素
Matrix_x * sub_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j;
    Matrix_x *mat3;//输出矩阵

    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//分配空间
    mat3->row=mat1->row;
    mat3->col=mat2->col;//定义输入矩阵行数和列数

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat3->value[i][j] = mat1->value[i][j] - mat2->value[i][j];
        }
    }

    return mat3;//返回输出矩阵
}

//矩阵减法，分数元素
Matrix * sub_(Matrix * mat1,Matrix * mat2)
{
    int i,j;
    Matrix* mat;//输出矩阵

    mat=(Matrix *)malloc(sizeof(Matrix));//分配空间
    mat->col=mat1->col;
    mat->row=mat1->row;//定义行数和列数

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat->value[i][j] = dec(&mat1->value[i][j],&mat2->value[i][j]);
        }
    }

    return mat;
}




//矩阵的乘法，浮点数元素
Matrix_x * mul_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j,k;
    Matrix_x *mat3;//输出矩阵

    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//分配空间

    	mat3->row=mat1->row;
    	mat3->col=mat2->col;//定义行数和列数

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat2->col;j++)
		{
        	mat3->value[i][j]=0;//将输出矩阵元素初始化为0
        	for(k=0;k<mat1->col;k++)
			{
        		mat3->value[i][j]+=(mat1->value[i][k]* mat2->value[k][j]);//矩阵乘法定义
			}

        }
    }

    return mat3;//返回输出矩阵
}


//矩阵乘法，分数元素
Matrix * mul_(Matrix * mat1,Matrix * mat2)
{
    int i,j,k;
    Matrix *mat3;//输出矩阵

    mat3=(Matrix *)malloc(sizeof(Matrix));//分配空间
    	mat3->row=mat1->row;
    	mat3->col=mat2->col;//定义行数和列数

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat2->col;j++)
		{
        	mat3->value[i][j].x=0;//初始化输出矩阵元素分子值
            mat3->value[i][j].y=1;//分母值
        	for(k=0;k<mat1->col;k++)
			{
                fraction b;//临时变量
                b=mul(&mat1->value[i][k],&mat2->value[k][j]);
                mat3->value[i][j]=add(&mat3->value[i][j],&b);
			}

        }
    }

    return mat3;//返回输出矩阵
}



//矩阵转置，浮点数元素
Matrix_x * transpose_decimal(Matrix_x * mat)
{
    int i,j;

    Matrix_x *mat1;//定义输出矩阵
    mat1=(Matrix_x *)malloc(sizeof(Matrix_x));
    mat1->row=mat->col;
    mat1->col=mat->row;

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat1->value[i][j] = mat->value[j][i];
        }
    }

    return mat1;//返回输出矩阵
}


//矩阵转置，分数元素
Matrix * transpose_(Matrix * mat)
{
    int i,j;

    Matrix *mat1;//输出矩阵
    mat1=(Matrix *)malloc(sizeof(Matrix));
    mat1->row=mat->col;
    mat1->col=mat->row;

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat1->value[i][j] = mat->value[j][i];
        }
    }

    return mat1;
}

//input Matrix;
Matrix * read(Matrix *  mat)
{
    int i,j;



    printf("input the numbers of rows and cols\n" );

    scanf("%d%d",&mat->row,&mat->col );//input row and col;

    printf("Input the matrix or determinant\n");
    for ( i = 0; i < mat->row; i++)
	{
        for(  j=0;j<mat->col;j++)
		{
            scanf("%d",&mat->value[i][j].x);
            if (getchar()=='/')
			 {
			 //fraction or not?
                scanf("%d",&mat->value[i][j].y);
            }
            else
			{
                mat->value[i][j].y=1;//denominator=1;
            }
        }
    }

    return mat;
}


//读入浮点数为元素的矩阵
Matrix_x * read_decimal(Matrix_x * mat)
{
    int i,j;

    printf("input the numbers of rows and cols\n" );

    scanf("%d%d",&mat->row,&mat->col );//输入行数和列数

    printf("Input the matrix or determinant\n");

    for ( i = 0; i < mat->row; i++)//输入矩阵的元素
	{
        for(  j=0;j<mat->col;j++)
		{
            double mid;//临时变量

            scanf("%lf",&mid);
            mat->value[i][j]=mid;
        }
    }

    return mat;
}

//分数行互换
void swap_row(Matrix * mat,int row1,int row2)
{
    int i;

    for ( i = 0; i <mat->col; i++)
	 {
        fraction trm;//临时变量
        trm=mat->value[row1][i];

        mat->value[row1][i]=mat->value[row2][i];
        mat->value[row2][i]=trm;//交换某两行的元素
    }

}


//小数行互换
void swap_row_decimal(Matrix_x * mat,int row1,int row2)
 {
    int i;

    for ( i = 0; i <mat->col; i++)
	{
        double trm;
        trm=mat->value[row1][i];

        mat->value[row1][i]=mat->value[row2][i];
        mat->value[row2][i]=trm;// 交换元素
    }

}




//化成矩阵最简型分数
Matrix * cha_into_dia(Matrix * mat)
{
    int m=0;
    int i,k,j;

    // print_matrix(mat);
    for (i = 0; i < mat->row; i++)//行
	 {
        int ok=1;

        if (mat->value[i][m].x==0)//判断元素是否为0
		 {
            ok=0;
            for (k = i+1; k < mat->row; k++)
			 {
                if (mat->value[k][m].x!=0)
				 {
                    ok=1;

                    swap_row(mat,i,k);
                    break;
                }
            }
        }

        if (ok==0)//此列元素i行及以下，均为0
		{
            i--;//下次循环仍在此行
            m++;//已判断列数

            if(m>=mat->col)
			{
                break;
            }
            continue;
        }


        for ( j = 0; j <mat->row; j++)//对i行以外的元素进行初等变换，使该列元素为零
		{
            if (mat->value[j][m].x==0||i==j)//已经为0，或i行不处理
			{
                continue;
            }

            fraction b;

            // printf("%d\n",i);
            // print_fraction(mat->value[j][m]);

            b=divide(&mat->value[j][m],&mat->value[i][m]);//该列 j行 初以 i行

            for ( k = m; k <mat->col; k++)
			 {
                fraction MUL=mul(&b,&mat->value[i][k]);

                mat->value[j][k]=dec(&mat->value[j][k],&MUL);//初等变换
            }
        }

        m++;//已处理列数

        if(m>=mat->col)
		{
            break;
        }
        // print_matrix(mat);
    }

    return mat;
}


//小数
Matrix_x * cha_into_dia_dicimal(Matrix_x * mat)
{
    int m=0;
    int i,k,j;

    // print_matrix(mat);
    for (i = 0; i < mat->row; i++)
	 {
        int ok=1;

        if (mat->value[i][m]==0)
		{
            ok=0;
            for (k = i+1; k < mat->row; k++)
			 {
                if (mat->value[k][m]!=0)
				{
                    ok=1;

                    swap_row_decimal(mat,i,k);
                    break;
                }
            }
        }

        if (ok==0)//已经为0，不用处理
		 {
            i--;
            m++;

            if(m>=mat->col)
			{
                break;
            }

            continue;
        }


        for ( j = 0; j <mat->row; j++)//
		 {
            if (mat->value[j][m]==0||i==j)
			{
                continue;
            }

            double b;

            // printf("%d\n",i);
            // print_fraction(mat->value[j][m]);

            b= mat->value[j][m]/mat->value[i][m];//该列 j行 初以 i行

            for ( k = m; k <mat->col; k++)
			 {
                double MUL=b*mat->value[i][k];
                mat->value[j][k]=mat->value[j][k]-MUL;
            }
        }

        m++;
        if(m>=mat->col)
		{
            break;
        }
        // print_matrix(mat);
    }

    return mat;
}


//判断分数型矩阵的行数和列数是否相等
int judge_rs(Matrix * mat)
{
	if(mat->row==mat->col)//判断条件
		return 1;
	else
		return 0;
 }


//判断浮点数型矩阵行数与列数是否相等
int judge_rs_x(Matrix_x * mat)
{
	if(mat->row==mat->col)//判断条件
		return 1;
	else
		return 0;
 }

//判断分数型矩阵是否为实对称矩阵
int judge_symmetric(Matrix * mat)
{
	int i,j;
	int ok=0;//用于判断的变量

	for(i=0;i<mat->row;i++)
		{
			for(j=i;j<mat->col;j++)
				{
					if((mat->value[i][j].x)!=(mat->value[j][i].x))//先判断分子是否相等
						ok++;
						else if((mat->value[i][j].y)!=(mat->value[j][i].y))//判断分母是否相等
							ok++;
				}
		}

		if(ok>0)
			return 0;//返回为假
		else
			return 1;//返回为真
}

//判断浮点型矩阵是否为实对称矩阵
int judge_symmetric_x(Matrix_x * mat)
{
	int i,j;
	int ok=0;//用于判断的变量

	for(i=0;i<mat->row;i++)
	{
		for(j=i;j<mat->col;j++)
			{
				if((mat->value[i][j])!=(mat->value[j][i]))//判断两个浮点数是否相等
					ok++;
			}
	}

		if(ok>0)
			return 0;//返回为假
		else
			return 1;//返回为真
}


//计算实对称矩阵的特征值和特征向量 雅克比迭代法
void jacbi(Matrix_x *mat ,Matrix_x * dblVect)
{

    int i,j;
    double dbEps = 1e-7;    //  精度
    int nJt = 20;   //最大迭代次数

    //初始化特征向量
    dblVect->col=mat->col;
    dblVect->row=mat->row;

    for(i=0;i<mat->col;i++)
	{
        for(j=0;j<mat->row;j++)
		{
            if(i==j)
			{
                dblVect->value[i][j] = 1.0;
            }
            else dblVect->value[i][j] = 0.0;
        }
    }

    int nCount =0;     //迭代次数
    while(1)
	{
        //在mat的非对角线上找到最大元素
        double dbMax = mat->value[0][1];
        int nRow = 0;//最大元素的行数
        int nCol = 1;//最大元素的列数

        for(i =0 ;i<mat->row;i++)
		{      //行
            for(j=0;j<mat->col;j++)
			{    //列
                double d =fabs(mat->value[i][j]);

                if((i!=j)&&(d>dbMax))
				{
                    dbMax = d;
                    nRow = i;
                    nCol = j;
                }
            }
        }

        if (dbMax<dbEps)//精度
		{
            break;
        }

        if(nCount > nJt)//迭代次数
		{
            break;
        }

        nCount++;
        //p行p列，p行q列，q行q列
        double dbApp = mat->value[nRow][nRow];
        double dbApq = mat->value[nRow][nCol];
        double dbAqq = mat->value[nCol][nCol];


        //计算旋转角度
        double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
        double dbSinTheta = sin(dbAngle);//sin
        double dbCosTheta = cos(dbAngle);//cos
        double dbSin2Theta = sin(2*dbAngle);//sin2
        double dbCos2Theta = cos(2*dbAngle);//cos2


        //矩阵旋转
        mat->value[nRow][nRow] = dbApp*dbCosTheta*dbCosTheta + dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
		mat->value[nCol][nCol] = dbApp*dbSinTheta*dbSinTheta + dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
		mat->value[nRow][nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
		mat->value[nCol][nRow] = mat->value[nRow][nCol];

        //p行或p列或q行或q列
        for (i=0;i<mat->row;i++)
		 {
                if((i!=nCol)&&(i!=nRow))
				{
                    double dbMiddle = mat->value[i][nRow];
                    mat->value[i][nRow] = mat->value[i][nCol]*dbSinTheta + dbMiddle*dbCosTheta;
                    mat->value[i][nCol] = mat->value[i][nCol]*dbCosTheta - dbMiddle*dbSinTheta;

                }
        }
        //q,q以外的元素
        for(j=0;j<mat->col;j++)
		{
            if((j!=nCol)&&(j!=nRow))
			{
                double dbMiddle = mat->value[nRow][j];
                mat->value[nRow][j] = mat->value[nCol][j]*dbSinTheta + dbMiddle*dbCosTheta;
                mat->value[nCol][j] = mat->value[nCol][j]*dbCosTheta - dbMiddle*dbSinTheta;
            }
        }

        //计算特征向量
        for(i=0;i<mat->row;i++)
		{
            double dbMiddle = dblVect->value[i][nRow];
            dblVect->value[i][nRow] = dblVect->value[i][nCol]*dbSinTheta + dbMiddle*dbCosTheta;
            dblVect->value[i][nCol] = dblVect->value[i][nCol]*dbCosTheta - dbMiddle*dbSinTheta;
        }

    }
    //输出特征值和对应的特征向量

    for(i=0;i<mat->row;i++)
	{
        printf("特征值：%f\n",mat->value[i][i] );
        printf("特征向量\n");

        for ( j = 0; j <dblVect->row; j++)
		 {
            printf("%f,",dblVect->value[j][i] );

            if(j==dblVect->row-1)//最后一行
			{
                printf("\n\n");
            }
        }
    }
}


//矩阵数乘
Matrix_x * scalar_decimal(Matrix_x *a,double k)
{
    Matrix_x *b = (Matrix_x *)malloc(sizeof(Matrix_x));
    int i,j;
    b->row=a->row;//行数
    b->col= a->col;//列
    for(i=0;i< a->row;i++)
    {
        for(j=0;j<a->col;j++)
        {
            b->value[i][j] = a->value[i][j] * k;
        }
    }
    return b;
}

//分数
Matrix * scalar(Matrix *a,fraction k)
{
    Matrix *b = (Matrix *)malloc(sizeof(Matrix));
    int i,j;
    b->row=a->row;
    b->col= a->col;
    for(i=0;i< a->row;i++)
    {
        for(j=0;j<a->col;j++)
        {
            b->value[i][j] = mul(&a->value[i][j],&k);
        }
    }
    return b;
}

//向量数乘
double  * scalar_(double *a,double b,int n,double *c)
{
    int i;
    for(i=0;i<n;i++)
    {
        c[i]=a[i] * b;
    }
    return c;
}


//向量内积
fraction dot_product(fraction *a,fraction * b,int n)
{
    int i;
    fraction k, ans;
    for(i=0;i<n;i++)
    {
        k = mul(&a[i],&b[i]);
        ans = add(&k,&ans);
    }
    return k;
}


double dot_product_decimal(double *a,double * b,int n)
{
    int i;
    double k=0;
    for(i=0;i<n;i++)
    {
        k += a[i] * b[i];
    }
    return k;
}

//向量长度
double lenth(double * a,int n)
{
    int i;
    double len=0.0;
    for(i=0;i<n;i++)
    {
        len += a[i]*a[i];
    }
    return sqrt(len);
}


//向量加法
double * add_vector(double * a,double * b,int n,double *c)
{
    int i;
    for(i=0;i<n;i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}



//施密特正交化法
void schmidt(double *a,int n)
{
    int i,j;
    double k[MAXN],b[MAXN][MAXN];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            b[i][j] = a[i*n+j];
        }
        for(j=0;j<i;j++)
        {   double an[MAXN];
            k[j] = - dot_product_decimal(a+i*n,b[i-j-1],n) / dot_product_decimal(b[i-j-1],b[i-j-1],n);
            scalar_(b[i-j-1],k[j],n,an);
            add_vector(a+i*n,an,n,b[i]);

        }
    }

    for(i=0;i<n;i++)
    {
        scalar_(b[i],1/lenth(b[i],n),n,b[i]);
        for(j=0;j<n;j++)//输出
        {
            printf("%15f  ", b[i][j]);
            if(j==n-1)
            {
                printf("\n");//行尾换行
            }
        }
    }
}


//计算行列式
fraction  cal_det(Matrix * mat)
{
    // printf("det?" );
    int i,j;

    Matrix mat_bak;
    for( i=0;i<mat->row;i++)
	{
        for( j=0;j<mat->col;j++)
		{
            mat_bak.value[i][j]=mat->value[i][j];
        }
    }

    mat=cha_into_dia(mat);//here!!!!
    for ( i = 0; i < mat->row; i++)
	 {
        if (mat->value[i][i].x==0)
		 {
            fraction ans;
            ans.x=0;
            ans.y=1;

            return ans;
        }
    }
    // printf("%d",mat_bak.value[0][1].x);
    fraction an;//here is wrong
    an.x=1;
    an.y=1;

    for ( i = 0; i <mat->row; i++)
	{
        an=mul(&an,&mat->value[i][i]);
    }
    for( i=0;i<mat->row;i++)
	{
        for( j=0;j<mat->col;j++)
		{
            mat->value[i][j]=mat_bak.value[i][j];
        }
    }

    return an;
}

//计算行列式
double  cal_det_decimal(Matrix_x * mat)
{
    // printf("det?" );
    int i,j;

    Matrix_x mat_bak;

    for( i=0;i<mat->row;i++)
	{
        for( j=0;j<mat->col;j++)
		{
            mat_bak.value[i][j]=mat->value[i][j];
        }
    }

    mat=cha_into_dia_dicimal(mat);//here!!!!
    for ( i = 0; i < mat->row; i++)
	{
        if (mat->value[i][i]==0)
		 {
            double ans = 0;
            return ans;
        }
    }
    // printf("%d",mat_bak.value[0][1].x);
    double an;//here is wrong
    an=1.0;

    for ( i = 0; i <mat->row; i++)
	 {
        an*=mat->value[i][i];
    }
    for( i=0;i<mat->row;i++)
	{
        for( j=0;j<mat->col;j++)
		{
            mat->value[i][j]=mat_bak.value[i][j];
        }
    }

    return an;
}


//判断是否是奇异矩阵
int is_singular(Matrix * mat)
{

    fraction an=cal_det(mat);

    if (an.x==0)
	{
        // printf("singular wrong?");
        return 1;
    }

    return 0;
}


//小数形式
int is_singular_decimal(Matrix_x * mat)
{

    double an=cal_det_decimal(mat);
    if (an==0)
	{
        // printf("singular wrong?");
        return 1;
    }

    return 0;
}

//二次型正惯性指数与负惯性指数
void cal_quad(Matrix_x *mat) {
	Matrix_x *ans_mat;
	ans_mat=(Matrix_x *)malloc(sizeof(Matrix_x));
	jacbi(mat,ans_mat);

	int i=0,j=0,k;
	for(k=0;k<mat->row;k++)
	{
		if(mat->value[k][k]>0.0)
			i++;

		if(mat->value[k][k]<0.0)
			j++;

	}
	printf("\n该矩阵正惯性指数为%d",i);
	printf("\n负惯性指数为%d\n",j);
	if(j==0)
		{
			if(i==mat->row)
				printf("\n该矩阵是正定矩阵！\n");
			else
				printf("\n该矩阵是半正定矩阵！\n");
		}
	else
		{
			if(i==0)
				{
					if(j==mat->row)
						printf("\n该矩阵为负定矩阵！\n");
					else
						printf("\n该矩阵为半负定矩阵！\n");
				}
			else
				printf("\n该矩阵为不定矩阵！\n\n\n\n");
		}

}


Matrix * get_inverse(Matrix * mat)
{
    //row!=col
    int i,j;

    if (mat->row!=mat->col)
	 {
        printf("row!=col,not a matrix!");

        exit(0);
    }

    //is singular or not?
    if (is_singular(mat))
	{
        printf("This matrix is a singular one which is not inversable.\n" );
        exit(0);
    }

    else
	{//inversable
        //
        for ( i = 0; i <mat->row; i++)
		 {
            for (j = 0; j <mat->col; j++)
			{
                mat->value[i][j+mat->col].y=1;

                if (i==j)
				{
                    mat->value[i][j+mat->col].x=1;
                }

                else mat->value[i][j+mat->col].x=0;
            }
        }

        mat->col*=2;
        // print_matrix(mat);
        mat=cha_into_dia(mat);
        // print_matrix(mat);
        mat->col/=2;

        for ( i = 0; i < mat->row; i++)
		 {
            fraction b;
            b=mat->value[i][i];

            for (j =0; j <mat->col; j++)
			{
                mat->value[i][j]=divide(&mat->value[i][j+mat->col],&b);
            }
        }
    }

    return mat;
}



Matrix_x * get_inverse_decimal(Matrix_x * mat)
{
    //row!=col
    int i,j;

    if (mat->row!=mat->col)
	{
        printf("row!=col,not a matrix!");
        exit(0);
    }

    //is singular or not?
    if (is_singular_decimal(mat))
	{
        printf("This matrix is a singular one which is not inversable.\n" );
        exit(0);
    }
    else
	{//inversable
        //
        for ( i = 0; i <mat->row; i++)
		{
            for (j = 0; j <mat->col; j++)
			{

                if (i==j)
				{
                    mat->value[i][j+mat->col]=1;
                }
                else mat->value[i][j+mat->col]=0;
            }
        }

        mat->col*=2;
        // print_matrix(mat);
        mat=cha_into_dia_dicimal(mat);
        // print_matrix(mat);
        mat->col/=2;

        for ( i = 0; i < mat->row; i++)
		{
            double b;
            b=mat->value[i][i];

            for (j =0; j <mat->col; j++)
			{
                mat->value[i][j] = mat->value[i][j+mat->col] / b;
            }
        }
    }

    return mat;
}


//求秩
int  rank(Matrix * mat)
{
    int i,j;

    mat=cha_into_dia(mat);

    int n=0;
    for( i=0;i<mat->row;i++)
	{
        for ( j = 0; j <mat->col; j++)
		{
            if (mat->value[i][j].x!=0)
			{
                n++;
                break;
            }
        }
    }
    return n;
}


//求秩
int  rank_decimal(Matrix_x * mat)
{
    int i,j;

    mat=cha_into_dia_dicimal(mat);

    int n=0;


    for( i=0;i<mat->row;i++)
	{
        for ( j = 0; j <mat->col; j++)
		{
            if (mat->value[i][j]!=0)
			{
                n++;
                break;
            }
        }
    }
    return n;
}

//求齐次方程的基础解系

#define MAX 100
int  syst(Matrix * mat)
{
    cha_into_dia(mat);

    int i,j;
    fraction x[MAXN][MAX];
    int fuzhi[mat->col],pos[mat->row];

    j=0;
    for(i=0;i<mat->row;i++)
	{
        for(;j<mat->col;j++)
		{
            //非自由未知量
            if(mat->value[i][j].x!=0)
			{
                fuzhi[j]=1;
                pos[j]=i;
                j++;

                break;
            }
            //自由未知量
            else
			{
                fuzhi[j]=0;

            }
        }
    }

    for(i=0;i<mat->col;i++)
	{
        if(fuzhi[i]!=1)
		{
            x[i][i].x=1;
            x[i][i].y=1;
            for(j=0;j<mat->col;j++)
			{
                if (i==j)
				{
                    continue;
                }
                //
                else if(fuzhi[j]==1)
				{
                    x[i][j] = divide(&mat->value[pos[j]][i],&mat->value[pos[j]][j]);
                    x[i][j].x = -x[i][j].x;
                }
                else
				{
                    x[i][j].x = 0;
                    x[i][j].y = 1;
                }
            }
        }
    }
    int m=1;
    printf("基础解系为：\n");
    for(i=0;i<mat->col;i++)
	{
        if(fuzhi[i]==0)
		{
            printf("X%d:(",m);
            m++;

            for(j=0;j<mat->col;j++)
			{
                print_fraction(x[i][j]);
                if(j!=mat->col-1)
				{
                    printf(" , " );
                }
                else
				printf("\n" );
            }
        }
    }
    return m;
}


int  syst_decimal(Matrix_x * mat)
{
    cha_into_dia_dicimal(mat);

    int i,j;
    double x[MAXN][MAX];
    int fuzhi[mat->col],pos[mat->row];

    j=0;
    for(i=0;i<mat->row;i++)
	{
        for(;j<mat->col;j++)
		{
            //非自由未知量
            if(mat->value[i][j]!=0)
			{
                fuzhi[j]=1;
                pos[j]=i;
                j++;

                break;
            }
            //自由未知量
            else
			{
                fuzhi[j]=0;

            }
        }
    }

    for(i=0;i<mat->col;i++)
	{
        if(fuzhi[i]!=1)
		{
            x[i][i]=1;
            for(j=0;j<mat->col;j++)
			{
                if (i==j)
				{
                    continue;
                }
                //
                else if(fuzhi[j]==1)
				{
                    x[i][j] = mat->value[pos[j]][i] / mat->value[pos[j]][j];
                    x[i][j] = -x[i][j];
                }
                else
				{
                    x[i][j] = 0;

                }
            }
        }
    }

    int m=1;
    printf("基础解系为：\n");

    for(i=0;i<mat->col;i++)
	{
        if(fuzhi[i]==0)
		{
            printf("X%d:(",m);
            m++;

            for(j=0;j<mat->col;j++)
			{
                printf("%f",x[i][j]);
                if(j!=mat->col-1)
				{
                    printf(" , " );
                }
                else
					printf(")\n" );
            }
        }
    }

    return m;
}

//齐次线性方程组

void linear_ab(Matrix *mat)
{
    int i;


    //mat为系数矩阵

    if(rank(mat)==mat->col)
	{
        //只有零解
        printf("x = 0\n" );
    }
    else
	{
        cha_into_dia(mat);

        int m;//解向量数量
        m=syst(mat);

        printf("一般解为：\n  X = ");

        for(i=1;i<m;i++)
		{
            printf("k%d X%d",i,i);

            if(i!=m-1)
			{
                printf(" + ");
            }

            else
				printf("\n");
        }
    }
}


void linear_ab_decimal(Matrix_x *mat)
{
    int i;

    //mat为系数矩阵

    if(rank_decimal(mat)==mat->col)
	{
        //只有零解
        printf("x = 0\n" );
    }
    else
	{
        cha_into_dia_dicimal(mat);

        int m;//解向量数量
        m=syst_decimal(mat);

        printf("一般解为：\n  X = ");

        for(i=1;i<m;i++)
		{
            printf("k%d X%d",i,i);
            if(i!=m-1)
			{
                printf(" + ");
            }
            else
				printf("\n");
        }
    }
}


//求解非齐次线性方程组
void linear_nor(Matrix *mat)
{
    int i,j;

    Matrix co_mat;//coefficient matrix;
    co_mat.row=mat->row;
    co_mat.col=mat->col-1;

    for (  i = 0; i <co_mat.row; i++)
	{
        for ( j = 0; j <co_mat.col; j++)
		{
            co_mat.value[i][j]=mat->value[i][j];

        }
    }

    int r1=rank(mat);

    int r2=rank(&co_mat);


    if (r1!=r2)
	{
        printf("no answer!");
    }

    //只有唯一解，b=0时，只有零解
    else if (r1==r2&&r1==co_mat.col)
	{
        cha_into_dia(mat);

        for (i=0;i<mat->col-1;i++)
		{
            fraction ans;
            ans = divide(&mat->value[i][mat->col-1],&mat->value[i][i]);

            printf("x%d = ",i+1 );
            print_fraction(ans);
            printf("\n");
        }

    }

    //无穷多个解
    else
	{
        cha_into_dia(mat);
        //特解
        printf("特解为：\n   X0 = (" );

        fraction x0[mat->col];
        int fuzhi[mat->col];

        memset(fuzhi,0,mat->col);
        x0[0].x=0;
        j=0;

        for(i=0;i<mat->row;i++)
		{
            for(;j<mat->col-1;j++)
			{
                if (mat->value[i][j].x!=0)
				{
                    x0[j]=divide(&mat->value[i][mat->col-1],&mat->value[i][j]);
                    fuzhi[j]=1;
                    j++;

                    print_fraction(x0[j]);

                    if(j!=mat->col-2)
					{
                        printf(" , " );
                    }
                    else
						printf(")\n" );

                    break;
                }
                else
				{
                    x0[j].x=0;
                    x0[j].y=1;
                    fuzhi[j]=0;
                }
                print_fraction(x0[j]);
                //格式控制
                if(j!=mat->col-2)
				{
                    printf(" , " );
                }
                else printf(")\n" );
            }
            if(j==mat->col-2)
				break;
        }

        //AX=0的基础解系
        cha_into_dia(&co_mat);
        int m = syst(&co_mat);

        //一般解
        printf("一般解为：\n  X = x0 ");

        for(i=1;i<m;i++)
		{
            printf("+ k%dX%d",i,i );
            if (i==m-1)
			{
                printf("\n");
            }
        }

    }
}

//线性方程
void linear_nor_decimal(Matrix_x *mat)
{
    int i,j;

    Matrix_x co_mat;//coefficient matrix;
    co_mat.row=mat->row;
    co_mat.col=mat->col-1;

    for (  i = 0; i <co_mat.row; i++)
	{
        for ( j = 0; j <co_mat.col; j++)
		{
            co_mat.value[i][j]=mat->value[i][j];

        }
    }

    int r1=rank_decimal(mat);

    int r2=rank_decimal(&co_mat);


    if (r1!=r2)
	{
        printf("no answer!");
    }

    //只有唯一解，b=0时，只有零解
    else if (r1==r2&&r1==co_mat.col)
	{
        cha_into_dia_dicimal(mat);
        for (i=0;i<mat->col-1;i++)
		{
            double ans;
            ans = mat->value[i][mat->col-1] / mat->value[i][i];

            printf("x%d = ",i+1 );
            printf("%f",ans);
            printf("\n");
        }

    }

    //无穷多个解
    else
	{
        cha_into_dia_dicimal(mat);
        //特解
        printf("特解为：\n   X0 = (" );

        double x0[mat->col];

        int fuzhi[mat->col];

        memset(fuzhi,0,mat->col);
        x0[0]=0;
        j=0;

        for(i=0;i<mat->row;i++)
		{
            for(;j<mat->col-1;j++)
			{
                if (mat->value[i][j]!=0)
				{
                    x0[j]= mat->value[i][mat->col-1] / mat->value[i][j];
                    fuzhi[j]=1;
                    j++;

                    printf("%f",x0[j]);

                    if(j!=mat->col-2)
					{
                        printf(" , " );
                    }
                    else
						printf(")\n" );
                    break;
                }
                else
				{
                    x0[j]=0;
                    fuzhi[j]=0;
                }
                printf("%f",x0[j]);
                if(j!=mat->col-2)
				{
                    printf(" , " );
                }
                else printf(")\n" );
            }
            if(j==mat->col-2)
				break;
        }

        //AX=0的基础解系
        cha_into_dia_dicimal(&co_mat);
        int m = syst_decimal(&co_mat);

        //一般解
        printf("一般解为：\n  X = x0 ");
        for(i=1;i<m;i++)
		{
            printf("+ k%dX%d",i,i );

            if (i==m-1)
			{
                printf("\n");
            }
        }

    }
}

int  input(double * a)
{
    printf("输入向量维数\n" );
    int n;
    scanf("%d",&n );
    int i;
    for(i=0;i<n;i++)
    {
        scanf("%lf", a+i);
    }
    return n;
}

//输入提示
int  print_tips()
{
    printf("This is a Matrix operator.\n");
    printf("0   into the order mode(命令模式)\n");
    printf("1   matrix(矩阵)\n");
    printf("2   caculate detaminate(行列式)\n" );
    printf("3   linear equations(线性方程组)\n" );
    printf("4   向量运算\n");//这部分未完成
    printf("-1  退出程序\n");

    printf("Please input the number before your wanting operation.\n" );
    //waiting to add function;
    int n;
    scanf("%d", &n);
    return n;
}

int operations(int n)
{
    //数字操作
    if (n==1)//矩阵
	{
        printf("1   矩阵求逆\n");
        printf("2   两个矩阵之间运算\n" );
        printf("3   矩阵求秩\n");
        printf("4   矩阵转置\n");
        printf("5   矩阵特征值和特征向量\n");
        printf("Please input the number before your wanting operation.\n" );
        scanf("%d",&n );

        if(n==1)
        {
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            get_inverse(mat);

            print_matrix(mat);
        }

        if (n==2)
        {
            printf("第一个矩阵\n" );
            Matrix_x *mat1 = (Matrix_x *)malloc(sizeof(Matrix_x));
            mat1 = read_decimal(mat1);

            printf("输入你想进行的操作（加“+”，减“-”，乘“*”\n" );
            char sign;
            scanf("%c",&sign );

            printf("第二个矩阵\n" );
            Matrix_x *mat2 = (Matrix_x *)malloc(sizeof(Matrix_x));
            mat2 = read_decimal(mat2);

             if(sign =='+')
             {
                 print_matrix_x(add_decimal(mat1,mat2));
             }

             if (sign == '-')
             {
                 print_matrix_x(sub_decimal(mat1,mat2));
             }

             if(sign == '*')
             {
                     print_matrix_x(mul_decimal(mat1,mat2));
             }


        }
        if(n==3)
        {
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            int ans=rank(mat);
            printf("%d\n",ans );
        }

        if(n==4)
        {
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            print_matrix(transpose_(mat));
        }

        if(n==5)
        {
            Matrix_x *mat;
            mat=(Matrix_x *)malloc(sizeof(Matrix_x));
            mat =read_decimal(mat);
            Matrix_x * ans_mat;
            ans_mat = (Matrix_x * )malloc(sizeof(Matrix_x));

            jacbi(mat,ans_mat);
        }

    }

    if (n==2)//行列式
	{
        Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
        mat=read(mat);
        fraction ans;
        ans=cal_det(mat);

        print_fraction(ans);
    }

    if (n==3)//线性方程组
	{
        printf("非齐次线性方程输1\n齐次线性方程输入2\n");
        int choose;
        scanf("%d",&choose);

        if(choose==1)
		{
            printf("Input the augmented matrix(增广矩阵)\n" );
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            linear_nor(mat);
        }

        else
		{
            printf("输入系数矩阵\n" );
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            linear_ab(mat);

        }
    }

    if (n==4)//未完成
	{
        //向量输入部分未完成
        printf("1   求向量长度\n");//finished lenth()
        printf("2   向量内积\n");
        scanf("%d\n",&n );

        if (n==1)
        {
            double * a = (double * )malloc(MAXN*sizeof(double));
            int l = input(a);
            printf("向量长度为：\t%f",lenth(a,l) );
        }

        if (n==2)
        {
            double * a = (double * )malloc(MAXN*sizeof(double));
            int la = input(a);
            double * b = (double * )malloc(MAXN*sizeof(double));
            int lb = input(b);
            if(la==lb)
            {
                dot_product_decimal(a,b,la);
            }
            else
            {
                    printf("错误操作\n");
            }
        }
    }

    if (n==-1)
	{
        exit(0);
    }


    return 0;
}

//命令模式

int orders(int *n,int *m,Matrix_x *A[],Matrix *A_[],char *name,char*name_)
{
    int i,j;
    char s[MAXN];


    scanf("%s",s);

    //finished
    if(s[0]>='A'&&s[0]<='Z'&&s[1]=='=')//输入矩阵
    {
        name[*m]=s[0];
        printf("默认为小数，分数请使用命令“X_=”\n" );
        A[*m]= (Matrix_x *)malloc(sizeof(Matrix_x));//分配内存
        A[*m] = read_decimal(A[*m]);
        printf("%c=\n",name[*m] );
        print_matrix_x(A[*m]);
        (*m)++;
        return 0;
    }

    else if(s[0]>='A'&&s[0]<='Z'&&s[1]=='_'&&s[2]=='=')//分数
    {
        name_[*n]=s[0];
        A_[*n]= (Matrix *)malloc(sizeof(Matrix));//分配内存
        A_[*n] =read(A_[*n]);
        printf("%c=\n",name_[*n] );
        print_matrix(A_[*n]);
        (*n)++;
        return 0;

    }

    int k=0;

    if(s[0]=='r'&&s[1]=='a'&&s[2]=='n'&&s[3]=='k'&&s[4]=='(')//求秩
    {
        if(s[6]=='_')//分数
        {
            for(i=0;i<(*n);i++)
            {
                if(s[5]==name_[i])

                {
                    k=i;
                    break;
                }
                if(i==*n)
                {
                    printf("No Such Name\n");
                    return 0;
                }
            }
            int ans = rank(A_[k]);
            printf("rank(%c_)\t%d\n",name_[k],ans);
            return 0;
        }
        else//小数
        {
            for(i=0;i<(*m);i++)
            {
                if(s[5]==name[i])
                {
                    k=i;
                    break;
                }
                if(i==*m)//数组内没有此名字
                {
                    printf("No Such Name\n");
                    return 0;
                }

            }

            int ans = rank_decimal(A[k]);
            printf("rank(%c) =\t %d\n",name[k],ans );
            return 0;
        }
    }

    //
    if(s[0]=='d'&&s[1]=='e'&&s[2]=='t'&&s[3]=='(')//det,求行列式
    {
        for(i=0;i<(*m);i++)
        {
            if(s[4]==name[i])
            {
                k=i;
                break;
            }
            if(i==*m)
            {
                printf("No Such Name\n");
                return 0;
            }

        }

        if(s[5]=='_')//小数
        {
            fraction ans = cal_det(A_[k]);
            printf("det(%c)=\n\t\t",name_[k] );
            print_fraction(ans);
        }
        else
        {
            double ans = cal_det_decimal(A[k]);
            printf("det(%c)=\n\t\t%f",name[k],ans);

        }

    }

    if(s[0]=='i'&&s[1]=='n'&&s[2]=='v')//求逆
    {


        if(s[5]=='_')
        {
            for(i=0;i<(*n);i++)
            {
                if(s[4]==name_[i])
                {
                    k=i;
                    break;
                }
                if(i==*n)
                {
                    printf("No Such Name\n");
                    return 0;
                }

            }
            printf("求逆：\n" );
            print_matrix(get_inverse(A_[k]));
        }
        else{
            for(i=0;i<(*m);i++)
            {
                if(s[4]==name[i])
                {
                    k=i;
                    break;
                }
                if(i==*m)
                {
                    printf("No Such Name\n");
                    return 0;
                }

            }
            printf("求逆：\n" );
            print_matrix_x(get_inverse_decimal(A[k]));
            return 0;
        }
    }

    //k*A
    //todo  k值已确定
    if(s[0]>='0'&&s[0]<='9')
    {

        int a = s[0]-'0';//k在符号前的部分，小数的整数部分或分数的分子
        double b;
        fraction c;
        int l=strlen(s);
        for(i=1;i<l;i++)//确定a的值
        {
            j=1;
            if(s[i]>='0'&&s[i]<='9')//数字00
            {
                a*=10;
                a+=(s[i]-'0');
                j=i+1;
            }

            else if(s[i]=='.')//小数点
            {


                double  lo=10;

                b=(s[i+1]-'0')/lo;
                // printf("%f\n",b );
                for(j=i+2;j<l;j++)//小数点后部分
                {
                    if(s[j]>='0'&&s[j]<='9')
                    {
                        lo *= 10;
                        b+=((s[j]-'0')/lo);
                    }
                    else break;
                }

            }

            b+=a;
            // TEST :NOPROBLEM
            if(s[j]=='*')
            {

                int k;
                for(k=0;k<(*m);k++)
                {
                    if(s[j+1]==name[k])
                    {
                        // printf("%c",name[k] );
                        int ak=k;
                        print_matrix_x(scalar_decimal(A[ak],b));
                        return 0;
                    }
                    if(k==*m)
                    {
                        printf("No Such Name\n");
                        return 0;
                    }
                }
            }
            else if(s[i]=='/')//分数
            {
                c.x=a;
                a=s[i+1];
                for(j=i+2;j<l;j++)
                {
                    if(s[j]>='0'&&s[j]<='9')
                    {
                        a*=10;
                        a+=(s[j]-'0');
                    }
                    else break;
                }
                c.y = a;

                int k;
                for(k=0;k<(*m);k++)
                {
                    if(s[j+1]==name[k])
                    {
                        int ak=k;
                        print_matrix(scalar(A_[ak],c));
                        return 0;
                    }
                    if(i==*n)
                    {
                        printf("No Such Name\n");
                        return 0;
                    }
                }
            }

        }
        return 0;
    }



    //schmidt施密特正交化
    if(s[0]=='s'&&s[1]=='m'&&s[2]=='t')//smt(A,B)    A->row=B->row = 1
    {

        int numn,num;
        printf("输入向量个数与向量维数\n");
        scanf("%d%d\n",&numn,&num );
        double a[MAXN];
        for(k=0;k<numn;k++)
        {
            printf("输入第%d个向量\n",k+1 );
            for(j=0;j<num;j++)
            {
                scanf("%lf\n",&a[k*num+j] );
            }
        }

        schmidt(a,num);
    }



    //A+B||A-B||A*B||A^K||A'
    for(i=0;i<(*m);i++)
    {
        if(s[0]==name[i])
        {
            k=i;
            if(s[1]=='+')
            {
                for(j=0;j<(*m);j++)
                {
                    if(s[2]==name[j])
                    {
                        int ak=j;
                        print_matrix_x(add_decimal(A[k],A[ak]));
                        return 0;
                    }
                }
            }
            else if(s[1]=='-')//A-B
            {
                for(j=0;j<(*m);j++)
                {
                    if(s[2]==name[j])
                    {
                        int ak=j;
                        print_matrix_x(sub_decimal(A[k],A[ak]));
                        return 0;
                    }
                }
            }
            else if(s[1]=='*')//A*B
            {
                for(j=0;j<(*m);j++)
                {
                    if(s[2]==name[j])
                    {
                        int ak=j;
                        if(A[k]->col!=A[ak]->row)
                        {
                            printf("无法乘");
                        }
                        else
                            print_matrix_x(mul_decimal(A[k],A[ak]));
                        return 0;
                    }
                }
            }
            //转置

            else if(s[1]=='\'')
            {
                print_matrix_x(transpose_decimal(A[k]));
                return 0;
            }
            //A_小数矩阵运算  //A+B||A-B||A*B||A^K||A'||A.'

            else if(s[1]=='_')
            {
                if(s[2]=='+')
                {
                    for(j=0;j<(*n);j++)
                    {
                        if(s[3]==name_[j])
                        {
                            int ak=j;
                            print_matrix(add_(A_[k],A_[ak]));
                            return 0;
                        }
                    }
                }
                else if(s[2]=='-')
                {
                    for(j=0;j<(*n);j++)
                    {
                        if(s[3]==name_[j])
                        {
                            int ak=j;
                            print_matrix(sub_(A_[k],A_[ak]));
                            return 0;
                        }
                    }
                }
                else if(s[2]=='*')
                {
                    for(j=0;j<(*n);j++)
                    {
                        if(s[3]==name_[j])
                        {
                            int ak=j;
                            if(A_[k]->col!=A_[ak]->row)
                            {
                                printf("错误操作，请检查是否符合矩阵相乘规则\n");
                            }
                            else
                                print_matrix(mul_(A_[k],A_[ak]));
                            return 0;
                        }
                    }
                }
                else if(s[2]=='\'')
                {
                    print_matrix(transpose_(A_[k]));
                }
            }
            return 0;
        }
    }
    printf("无效命令\n");
    return 0;
}
