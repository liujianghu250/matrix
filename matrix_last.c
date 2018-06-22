#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define MAXN 100

//��������ṹ��
typedef struct
{
    int x,y;//���ӡ���ĸ

} fraction;
//fraction:����


//�����Է���ΪԪ�صľ���ṹ
typedef struct Matrix
{
    //declare Matrix;
    int row,col;   //row:�С�col���С�
    fraction value[MAXN][MAXN];
}Matrix;
//matrix:����

// �����Ը�����ΪԪ�صľ���
typedef struct Matrix_x
{
    int col,row;
    double value[MAXN][MAXN];
}Matrix_x;

void print_fraction(fraction a);//��ӡ����
int print_matrix(Matrix * mat);//��ӡ����ΪԪ�صľ���
void print_matrix_x(Matrix_x * mat);
int gcd(int x,int y);//GCD �����׵����Լ����
void reduction(int* x,int* y);//Լ��
fraction add(fraction* a,fraction*b);//equals "+"
fraction add(fraction* a,fraction*b);//equals "-"
fraction mul(fraction*a, fraction*b);//equals "*"
fraction divide(fraction*a, fraction*b);//equals "/"
Matrix_x  *add_decimal(Matrix_x * mat1,Matrix_x * mat2);//����ӷ�
Matrix  *add_(Matrix * mat1,Matrix * mat2);//������ʽ
Matrix_x * sub_decimal(Matrix_x * mat1,Matrix_x * mat2);//���������������Ԫ��
Matrix * sub_(Matrix * mat1,Matrix * mat2);//�������������Ԫ��
Matrix_x * mul_decimal(Matrix_x * mat1,Matrix_x * mat2);//����ĳ˷���������Ԫ��
Matrix * mul_(Matrix * mat1,Matrix * mat2);//����˷�������Ԫ��
Matrix_x * transpose_decimal(Matrix_x * mat);//����ת�ã�������Ԫ��
Matrix * transpose_(Matrix * mat);//����ת�ã�����Ԫ��
Matrix * read(Matrix *  mat);//input Matrix;
Matrix_x * read_decimal(Matrix_x * mat);//���븡����ΪԪ�صľ���
void swap_row(Matrix * mat,int row1,int row2);//�����л���
void swap_row_decimal(Matrix_x * mat,int row1,int row2);//С���л���
Matrix * cha_into_dia(Matrix * mat);//���ɾ�������ͷ���
Matrix_x * cha_into_dia_dicimal(Matrix_x * mat);//С��
int judge_rs(Matrix * mat);//�жϷ����;���������������Ƿ����
int judge_rs_x(Matrix_x * mat);//�жϸ������;��������������Ƿ����
int judge_symmetric(Matrix * mat);//�ж��Ƿ��ǶԳƾ���
int judge_symmetric_x(Matrix_x * mat); //�жϸ����;����Ƿ�Ϊʵ�Գƾ���
void jacbi(Matrix_x *mat ,Matrix_x * dblVect);//����ʵ�Գƾ��������ֵ���������� �ſ˱ȵ�����
Matrix_x * scalar_decimal(Matrix_x *a,double k);//��������
Matrix * scalar(Matrix *a,fraction k);//����
double  * scalar_(double *a,double b,int n,double *c);//��������
fraction dot_product(fraction *a,fraction * b,int n);//�����ڻ�
double dot_product_decimal(double *a,double * b,int n);
double lenth(double * a,int n);//��������
double * add_vector(double * a,double * b,int n,double *c);//�����ӷ�
void schmidt(double *a,int n);//ʩ������������
fraction  cal_det(Matrix * mat);//��������ʽ
double  cal_det_decimal(Matrix_x * mat);//��������ʽ
int is_singular(Matrix * mat);//�ж��Ƿ����������
int is_singular_decimal(Matrix_x * mat);//С����ʽ
void cal_quad(Matrix_x *mat);//������������ָ���븺����ָ��
Matrix * get_inverse(Matrix * mat);//����
Matrix_x * get_inverse_decimal(Matrix_x * mat);
int  rank(Matrix * mat);//����
int  rank_decimal(Matrix_x * mat);//����
int  syst(Matrix * mat);//����η��̵Ļ�����ϵ
int  syst_decimal(Matrix_x * mat);
void linear_ab(Matrix *mat);//������Է�����
void linear_ab_decimal(Matrix_x *mat);
void linear_nor(Matrix *mat);//����������Է�����
void linear_nor_decimal(Matrix_x *mat);//���Է���
int  input(double * a);
int  print_tips();//������ʾ
int orders(int *n,int *m,Matrix_x *A[],Matrix *A_[],char *name,char*name_);//����ģʽ
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

//��ӡ����
void print_fraction(fraction a)
{
    if (a.y==1) //�����ĸΪһ
    {
        printf("%d",a.x );//ֻ��ӡ����
    }

    else 			//�����ĸ��Ϊһ
		printf("%d/%d",a.x,a.y);//���ӷ�ĸһ���ӡ
}


//��ӡ����ΪԪ�صľ���
int print_matrix(Matrix * mat)
{
    int i,j;//����ѭ���õı���

    printf("MAtrix:\n");//for beauty;
    for (i = 0; i < mat->row; i++)//���н���ѭ��
	{
        for ( j = 0; j < mat->col; j++)//��ѭ��
		{
        //this part can change into a function;
        if(mat->value[i][j].x>=0)
			printf(" ");//��ӡ�ո�ʹ�ô�ӡ�����ľ�������

        print_fraction(mat->value[i][j]);

        if(mat->value[i][j].y==1)//�����ĸΪһ
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
    printf("\t");//ʹ����ʱ����

    for(j=0;j<mat->col;j++)
    {
        printf("��%d��\t\t   ",j+1);//���������ʽ
    }

    printf("\n");

    for(i=0;i<mat->row;i++)
	{
        for(j=0;j<mat->col;j++)
		{
            if(j==0)
			{
                printf("��%d��\t",i+1 );//���״�ӡ��ǰ����
            }

            printf("%15lf ",mat->value[i][j] );

            if(j==mat->col-1)
			{
                printf("\n");//��β����
            }
        }
    }
}


//GCD �����׵����Լ����
int gcd(int x,int y)
{
    return y==0?x:gcd(y,x%y);
}


//Լ��
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

    int GCD=gcd(*x,*y);//���ӷ�ĸ���Լ��
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
    fraction result=//�������
    {
        a->x*b->y+a->y*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//Լ��

    return result;
}


//equals "-"
fraction dec(fraction* a,fraction*b)
{
    fraction result={//�������
        a->x*b->y-a->y*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//Լ��

    return result;
}


//equals "*"
fraction mul(fraction*a, fraction*b)
{
    fraction result={//�������
        a->x*b->x,
        a->y*b->y
    };

    reduction(&result.x,&result.y);//Լ��

    return result;
}


//equals "/"
fraction divide(fraction*a, fraction*b)
{
    fraction result ={//�������
        a->x*b->y,
        a->y*b->x
    };

    reduction(&result.x,&result.y);//Լ��

    return result;
}


//����ӷ�
Matrix_x  *add_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j;
    Matrix_x * mat3;
    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//�����ڴ�
    //���븳ֵ�������޷����
    mat3->col=mat1->col;
    mat3->row=mat1->row;

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat3->value[i][j]=mat1->value[i][j] + mat2->value[i][j];//��ӦԪ�����
        }
    }

    return mat3;
}

//������ʽ
Matrix  *add_(Matrix * mat1,Matrix * mat2)
{
    int i,j;
    Matrix* mat;
    mat=(Matrix *)malloc(sizeof(Matrix));//�����ڴ�
    //���и�ֵ
    mat->col=mat1->col;
    mat->row=mat1->row;
    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat->value[i][j] =add(&mat1->value[i][j],&mat2->value[i][j]);//���
        }
    }

    return mat;
}

//���������������Ԫ��
Matrix_x * sub_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j;
    Matrix_x *mat3;//�������

    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//����ռ�
    mat3->row=mat1->row;
    mat3->col=mat2->col;//���������������������

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat3->value[i][j] = mat1->value[i][j] - mat2->value[i][j];
        }
    }

    return mat3;//�����������
}

//�������������Ԫ��
Matrix * sub_(Matrix * mat1,Matrix * mat2)
{
    int i,j;
    Matrix* mat;//�������

    mat=(Matrix *)malloc(sizeof(Matrix));//����ռ�
    mat->col=mat1->col;
    mat->row=mat1->row;//��������������

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat1->col;j++)
		{
            mat->value[i][j] = dec(&mat1->value[i][j],&mat2->value[i][j]);
        }
    }

    return mat;
}




//����ĳ˷���������Ԫ��
Matrix_x * mul_decimal(Matrix_x * mat1,Matrix_x * mat2)
{
    int i,j,k;
    Matrix_x *mat3;//�������

    mat3=(Matrix_x *)malloc(sizeof(Matrix_x));//����ռ�

    	mat3->row=mat1->row;
    	mat3->col=mat2->col;//��������������

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat2->col;j++)
		{
        	mat3->value[i][j]=0;//���������Ԫ�س�ʼ��Ϊ0
        	for(k=0;k<mat1->col;k++)
			{
        		mat3->value[i][j]+=(mat1->value[i][k]* mat2->value[k][j]);//����˷�����
			}

        }
    }

    return mat3;//�����������
}


//����˷�������Ԫ��
Matrix * mul_(Matrix * mat1,Matrix * mat2)
{
    int i,j,k;
    Matrix *mat3;//�������

    mat3=(Matrix *)malloc(sizeof(Matrix));//����ռ�
    	mat3->row=mat1->row;
    	mat3->col=mat2->col;//��������������

    for(i=0;i < mat1->row;i++)
	{
        for(j=0;j<mat2->col;j++)
		{
        	mat3->value[i][j].x=0;//��ʼ���������Ԫ�ط���ֵ
            mat3->value[i][j].y=1;//��ĸֵ
        	for(k=0;k<mat1->col;k++)
			{
                fraction b;//��ʱ����
                b=mul(&mat1->value[i][k],&mat2->value[k][j]);
                mat3->value[i][j]=add(&mat3->value[i][j],&b);
			}

        }
    }

    return mat3;//�����������
}



//����ת�ã�������Ԫ��
Matrix_x * transpose_decimal(Matrix_x * mat)
{
    int i,j;

    Matrix_x *mat1;//�����������
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

    return mat1;//�����������
}


//����ת�ã�����Ԫ��
Matrix * transpose_(Matrix * mat)
{
    int i,j;

    Matrix *mat1;//�������
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


//���븡����ΪԪ�صľ���
Matrix_x * read_decimal(Matrix_x * mat)
{
    int i,j;

    printf("input the numbers of rows and cols\n" );

    scanf("%d%d",&mat->row,&mat->col );//��������������

    printf("Input the matrix or determinant\n");

    for ( i = 0; i < mat->row; i++)//��������Ԫ��
	{
        for(  j=0;j<mat->col;j++)
		{
            double mid;//��ʱ����

            scanf("%lf",&mid);
            mat->value[i][j]=mid;
        }
    }

    return mat;
}

//�����л���
void swap_row(Matrix * mat,int row1,int row2)
{
    int i;

    for ( i = 0; i <mat->col; i++)
	 {
        fraction trm;//��ʱ����
        trm=mat->value[row1][i];

        mat->value[row1][i]=mat->value[row2][i];
        mat->value[row2][i]=trm;//����ĳ���е�Ԫ��
    }

}


//С���л���
void swap_row_decimal(Matrix_x * mat,int row1,int row2)
 {
    int i;

    for ( i = 0; i <mat->col; i++)
	{
        double trm;
        trm=mat->value[row1][i];

        mat->value[row1][i]=mat->value[row2][i];
        mat->value[row2][i]=trm;// ����Ԫ��
    }

}




//���ɾ�������ͷ���
Matrix * cha_into_dia(Matrix * mat)
{
    int m=0;
    int i,k,j;

    // print_matrix(mat);
    for (i = 0; i < mat->row; i++)//��
	 {
        int ok=1;

        if (mat->value[i][m].x==0)//�ж�Ԫ���Ƿ�Ϊ0
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

        if (ok==0)//����Ԫ��i�м����£���Ϊ0
		{
            i--;//�´�ѭ�����ڴ���
            m++;//���ж�����

            if(m>=mat->col)
			{
                break;
            }
            continue;
        }


        for ( j = 0; j <mat->row; j++)//��i�������Ԫ�ؽ��г��ȱ任��ʹ����Ԫ��Ϊ��
		{
            if (mat->value[j][m].x==0||i==j)//�Ѿ�Ϊ0����i�в�����
			{
                continue;
            }

            fraction b;

            // printf("%d\n",i);
            // print_fraction(mat->value[j][m]);

            b=divide(&mat->value[j][m],&mat->value[i][m]);//���� j�� ���� i��

            for ( k = m; k <mat->col; k++)
			 {
                fraction MUL=mul(&b,&mat->value[i][k]);

                mat->value[j][k]=dec(&mat->value[j][k],&MUL);//���ȱ任
            }
        }

        m++;//�Ѵ�������

        if(m>=mat->col)
		{
            break;
        }
        // print_matrix(mat);
    }

    return mat;
}


//С��
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

        if (ok==0)//�Ѿ�Ϊ0�����ô���
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

            b= mat->value[j][m]/mat->value[i][m];//���� j�� ���� i��

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


//�жϷ����;���������������Ƿ����
int judge_rs(Matrix * mat)
{
	if(mat->row==mat->col)//�ж�����
		return 1;
	else
		return 0;
 }


//�жϸ������;��������������Ƿ����
int judge_rs_x(Matrix_x * mat)
{
	if(mat->row==mat->col)//�ж�����
		return 1;
	else
		return 0;
 }

//�жϷ����;����Ƿ�Ϊʵ�Գƾ���
int judge_symmetric(Matrix * mat)
{
	int i,j;
	int ok=0;//�����жϵı���

	for(i=0;i<mat->row;i++)
		{
			for(j=i;j<mat->col;j++)
				{
					if((mat->value[i][j].x)!=(mat->value[j][i].x))//���жϷ����Ƿ����
						ok++;
						else if((mat->value[i][j].y)!=(mat->value[j][i].y))//�жϷ�ĸ�Ƿ����
							ok++;
				}
		}

		if(ok>0)
			return 0;//����Ϊ��
		else
			return 1;//����Ϊ��
}

//�жϸ����;����Ƿ�Ϊʵ�Գƾ���
int judge_symmetric_x(Matrix_x * mat)
{
	int i,j;
	int ok=0;//�����жϵı���

	for(i=0;i<mat->row;i++)
	{
		for(j=i;j<mat->col;j++)
			{
				if((mat->value[i][j])!=(mat->value[j][i]))//�ж������������Ƿ����
					ok++;
			}
	}

		if(ok>0)
			return 0;//����Ϊ��
		else
			return 1;//����Ϊ��
}


//����ʵ�Գƾ��������ֵ���������� �ſ˱ȵ�����
void jacbi(Matrix_x *mat ,Matrix_x * dblVect)
{

    int i,j;
    double dbEps = 1e-7;    //  ����
    int nJt = 20;   //����������

    //��ʼ����������
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

    int nCount =0;     //��������
    while(1)
	{
        //��mat�ķǶԽ������ҵ����Ԫ��
        double dbMax = mat->value[0][1];
        int nRow = 0;//���Ԫ�ص�����
        int nCol = 1;//���Ԫ�ص�����

        for(i =0 ;i<mat->row;i++)
		{      //��
            for(j=0;j<mat->col;j++)
			{    //��
                double d =fabs(mat->value[i][j]);

                if((i!=j)&&(d>dbMax))
				{
                    dbMax = d;
                    nRow = i;
                    nCol = j;
                }
            }
        }

        if (dbMax<dbEps)//����
		{
            break;
        }

        if(nCount > nJt)//��������
		{
            break;
        }

        nCount++;
        //p��p�У�p��q�У�q��q��
        double dbApp = mat->value[nRow][nRow];
        double dbApq = mat->value[nRow][nCol];
        double dbAqq = mat->value[nCol][nCol];


        //������ת�Ƕ�
        double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
        double dbSinTheta = sin(dbAngle);//sin
        double dbCosTheta = cos(dbAngle);//cos
        double dbSin2Theta = sin(2*dbAngle);//sin2
        double dbCos2Theta = cos(2*dbAngle);//cos2


        //������ת
        mat->value[nRow][nRow] = dbApp*dbCosTheta*dbCosTheta + dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
		mat->value[nCol][nCol] = dbApp*dbSinTheta*dbSinTheta + dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
		mat->value[nRow][nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
		mat->value[nCol][nRow] = mat->value[nRow][nCol];

        //p�л�p�л�q�л�q��
        for (i=0;i<mat->row;i++)
		 {
                if((i!=nCol)&&(i!=nRow))
				{
                    double dbMiddle = mat->value[i][nRow];
                    mat->value[i][nRow] = mat->value[i][nCol]*dbSinTheta + dbMiddle*dbCosTheta;
                    mat->value[i][nCol] = mat->value[i][nCol]*dbCosTheta - dbMiddle*dbSinTheta;

                }
        }
        //q,q�����Ԫ��
        for(j=0;j<mat->col;j++)
		{
            if((j!=nCol)&&(j!=nRow))
			{
                double dbMiddle = mat->value[nRow][j];
                mat->value[nRow][j] = mat->value[nCol][j]*dbSinTheta + dbMiddle*dbCosTheta;
                mat->value[nCol][j] = mat->value[nCol][j]*dbCosTheta - dbMiddle*dbSinTheta;
            }
        }

        //������������
        for(i=0;i<mat->row;i++)
		{
            double dbMiddle = dblVect->value[i][nRow];
            dblVect->value[i][nRow] = dblVect->value[i][nCol]*dbSinTheta + dbMiddle*dbCosTheta;
            dblVect->value[i][nCol] = dblVect->value[i][nCol]*dbCosTheta - dbMiddle*dbSinTheta;
        }

    }
    //�������ֵ�Ͷ�Ӧ����������

    for(i=0;i<mat->row;i++)
	{
        printf("����ֵ��%f\n",mat->value[i][i] );
        printf("��������\n");

        for ( j = 0; j <dblVect->row; j++)
		 {
            printf("%f,",dblVect->value[j][i] );

            if(j==dblVect->row-1)//���һ��
			{
                printf("\n\n");
            }
        }
    }
}


//��������
Matrix_x * scalar_decimal(Matrix_x *a,double k)
{
    Matrix_x *b = (Matrix_x *)malloc(sizeof(Matrix_x));
    int i,j;
    b->row=a->row;//����
    b->col= a->col;//��
    for(i=0;i< a->row;i++)
    {
        for(j=0;j<a->col;j++)
        {
            b->value[i][j] = a->value[i][j] * k;
        }
    }
    return b;
}

//����
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

//��������
double  * scalar_(double *a,double b,int n,double *c)
{
    int i;
    for(i=0;i<n;i++)
    {
        c[i]=a[i] * b;
    }
    return c;
}


//�����ڻ�
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

//��������
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


//�����ӷ�
double * add_vector(double * a,double * b,int n,double *c)
{
    int i;
    for(i=0;i<n;i++)
    {
        c[i] = a[i] + b[i];
    }
    return c;
}



//ʩ������������
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
        for(j=0;j<n;j++)//���
        {
            printf("%15f  ", b[i][j]);
            if(j==n-1)
            {
                printf("\n");//��β����
            }
        }
    }
}


//��������ʽ
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

//��������ʽ
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


//�ж��Ƿ����������
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


//С����ʽ
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

//������������ָ���븺����ָ��
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
	printf("\n�þ���������ָ��Ϊ%d",i);
	printf("\n������ָ��Ϊ%d\n",j);
	if(j==0)
		{
			if(i==mat->row)
				printf("\n�þ�������������\n");
			else
				printf("\n�þ����ǰ���������\n");
		}
	else
		{
			if(i==0)
				{
					if(j==mat->row)
						printf("\n�þ���Ϊ��������\n");
					else
						printf("\n�þ���Ϊ�븺������\n");
				}
			else
				printf("\n�þ���Ϊ��������\n\n\n\n");
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


//����
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


//����
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

//����η��̵Ļ�����ϵ

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
            //������δ֪��
            if(mat->value[i][j].x!=0)
			{
                fuzhi[j]=1;
                pos[j]=i;
                j++;

                break;
            }
            //����δ֪��
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
    printf("������ϵΪ��\n");
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
            //������δ֪��
            if(mat->value[i][j]!=0)
			{
                fuzhi[j]=1;
                pos[j]=i;
                j++;

                break;
            }
            //����δ֪��
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
    printf("������ϵΪ��\n");

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

//������Է�����

void linear_ab(Matrix *mat)
{
    int i;


    //matΪϵ������

    if(rank(mat)==mat->col)
	{
        //ֻ�����
        printf("x = 0\n" );
    }
    else
	{
        cha_into_dia(mat);

        int m;//����������
        m=syst(mat);

        printf("һ���Ϊ��\n  X = ");

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

    //matΪϵ������

    if(rank_decimal(mat)==mat->col)
	{
        //ֻ�����
        printf("x = 0\n" );
    }
    else
	{
        cha_into_dia_dicimal(mat);

        int m;//����������
        m=syst_decimal(mat);

        printf("һ���Ϊ��\n  X = ");

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


//����������Է�����
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

    //ֻ��Ψһ�⣬b=0ʱ��ֻ�����
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

    //��������
    else
	{
        cha_into_dia(mat);
        //�ؽ�
        printf("�ؽ�Ϊ��\n   X0 = (" );

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
                //��ʽ����
                if(j!=mat->col-2)
				{
                    printf(" , " );
                }
                else printf(")\n" );
            }
            if(j==mat->col-2)
				break;
        }

        //AX=0�Ļ�����ϵ
        cha_into_dia(&co_mat);
        int m = syst(&co_mat);

        //һ���
        printf("һ���Ϊ��\n  X = x0 ");

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

//���Է���
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

    //ֻ��Ψһ�⣬b=0ʱ��ֻ�����
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

    //��������
    else
	{
        cha_into_dia_dicimal(mat);
        //�ؽ�
        printf("�ؽ�Ϊ��\n   X0 = (" );

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

        //AX=0�Ļ�����ϵ
        cha_into_dia_dicimal(&co_mat);
        int m = syst_decimal(&co_mat);

        //һ���
        printf("һ���Ϊ��\n  X = x0 ");
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
    printf("��������ά��\n" );
    int n;
    scanf("%d",&n );
    int i;
    for(i=0;i<n;i++)
    {
        scanf("%lf", a+i);
    }
    return n;
}

//������ʾ
int  print_tips()
{
    printf("This is a Matrix operator.\n");
    printf("0   into the order mode(����ģʽ)\n");
    printf("1   matrix(����)\n");
    printf("2   caculate detaminate(����ʽ)\n" );
    printf("3   linear equations(���Է�����)\n" );
    printf("4   ��������\n");//�ⲿ��δ���
    printf("-1  �˳�����\n");

    printf("Please input the number before your wanting operation.\n" );
    //waiting to add function;
    int n;
    scanf("%d", &n);
    return n;
}

int operations(int n)
{
    //���ֲ���
    if (n==1)//����
	{
        printf("1   ��������\n");
        printf("2   ��������֮������\n" );
        printf("3   ��������\n");
        printf("4   ����ת��\n");
        printf("5   ��������ֵ����������\n");
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
            printf("��һ������\n" );
            Matrix_x *mat1 = (Matrix_x *)malloc(sizeof(Matrix_x));
            mat1 = read_decimal(mat1);

            printf("����������еĲ������ӡ�+��������-�����ˡ�*��\n" );
            char sign;
            scanf("%c",&sign );

            printf("�ڶ�������\n" );
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

    if (n==2)//����ʽ
	{
        Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
        mat=read(mat);
        fraction ans;
        ans=cal_det(mat);

        print_fraction(ans);
    }

    if (n==3)//���Է�����
	{
        printf("��������Է�����1\n������Է�������2\n");
        int choose;
        scanf("%d",&choose);

        if(choose==1)
		{
            printf("Input the augmented matrix(�������)\n" );
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            linear_nor(mat);
        }

        else
		{
            printf("����ϵ������\n" );
            Matrix *mat = (Matrix *)malloc(sizeof(Matrix));
            mat=read(mat);
            linear_ab(mat);

        }
    }

    if (n==4)//δ���
	{
        //�������벿��δ���
        printf("1   ����������\n");//finished lenth()
        printf("2   �����ڻ�\n");
        scanf("%d\n",&n );

        if (n==1)
        {
            double * a = (double * )malloc(MAXN*sizeof(double));
            int l = input(a);
            printf("��������Ϊ��\t%f",lenth(a,l) );
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
                    printf("�������\n");
            }
        }
    }

    if (n==-1)
	{
        exit(0);
    }


    return 0;
}

//����ģʽ

int orders(int *n,int *m,Matrix_x *A[],Matrix *A_[],char *name,char*name_)
{
    int i,j;
    char s[MAXN];


    scanf("%s",s);

    //finished
    if(s[0]>='A'&&s[0]<='Z'&&s[1]=='=')//�������
    {
        name[*m]=s[0];
        printf("Ĭ��ΪС����������ʹ�����X_=��\n" );
        A[*m]= (Matrix_x *)malloc(sizeof(Matrix_x));//�����ڴ�
        A[*m] = read_decimal(A[*m]);
        printf("%c=\n",name[*m] );
        print_matrix_x(A[*m]);
        (*m)++;
        return 0;
    }

    else if(s[0]>='A'&&s[0]<='Z'&&s[1]=='_'&&s[2]=='=')//����
    {
        name_[*n]=s[0];
        A_[*n]= (Matrix *)malloc(sizeof(Matrix));//�����ڴ�
        A_[*n] =read(A_[*n]);
        printf("%c=\n",name_[*n] );
        print_matrix(A_[*n]);
        (*n)++;
        return 0;

    }

    int k=0;

    if(s[0]=='r'&&s[1]=='a'&&s[2]=='n'&&s[3]=='k'&&s[4]=='(')//����
    {
        if(s[6]=='_')//����
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
        else//С��
        {
            for(i=0;i<(*m);i++)
            {
                if(s[5]==name[i])
                {
                    k=i;
                    break;
                }
                if(i==*m)//������û�д�����
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
    if(s[0]=='d'&&s[1]=='e'&&s[2]=='t'&&s[3]=='(')//det,������ʽ
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

        if(s[5]=='_')//С��
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

    if(s[0]=='i'&&s[1]=='n'&&s[2]=='v')//����
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
            printf("���棺\n" );
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
            printf("���棺\n" );
            print_matrix_x(get_inverse_decimal(A[k]));
            return 0;
        }
    }

    //k*A
    //todo  kֵ��ȷ��
    if(s[0]>='0'&&s[0]<='9')
    {

        int a = s[0]-'0';//k�ڷ���ǰ�Ĳ��֣�С�����������ֻ�����ķ���
        double b;
        fraction c;
        int l=strlen(s);
        for(i=1;i<l;i++)//ȷ��a��ֵ
        {
            j=1;
            if(s[i]>='0'&&s[i]<='9')//����00
            {
                a*=10;
                a+=(s[i]-'0');
                j=i+1;
            }

            else if(s[i]=='.')//С����
            {


                double  lo=10;

                b=(s[i+1]-'0')/lo;
                // printf("%f\n",b );
                for(j=i+2;j<l;j++)//С����󲿷�
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
            else if(s[i]=='/')//����
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



    //schmidtʩ����������
    if(s[0]=='s'&&s[1]=='m'&&s[2]=='t')//smt(A,B)    A->row=B->row = 1
    {

        int numn,num;
        printf("������������������ά��\n");
        scanf("%d%d\n",&numn,&num );
        double a[MAXN];
        for(k=0;k<numn;k++)
        {
            printf("�����%d������\n",k+1 );
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
                            printf("�޷���");
                        }
                        else
                            print_matrix_x(mul_decimal(A[k],A[ak]));
                        return 0;
                    }
                }
            }
            //ת��

            else if(s[1]=='\'')
            {
                print_matrix_x(transpose_decimal(A[k]));
                return 0;
            }
            //A_С����������  //A+B||A-B||A*B||A^K||A'||A.'

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
                                printf("��������������Ƿ���Ͼ�����˹���\n");
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
    printf("��Ч����\n");
    return 0;
}
