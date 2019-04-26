//
//  main.cpp
//  LDA模块
//
//  Created by tsc on 17/4/20.
//  Copyright © 2017年 tsc. All rights reserved.
//

#include <iostream>
#include <unistd.h>
#include <string>
#include <cstring>
using namespace std;

class 数据处理{
public:
	static void 释放二维指针数组内存(int **&二维数组,int 行){
		for(int i=0;i<行;i++){
			delete[] 二维数组[i];
			二维数组[i]=NULL;
		}
		delete[] 二维数组;
		二维数组=NULL;
	}

	static void 释放二维指针数组内存(double **&二维数组,int 行){
		for(int i=0;i<行;i++){
			delete[] 二维数组[i];
			二维数组[i]=NULL;
		}
		delete[] 二维数组;
		二维数组=NULL;
	}

	static void 创建二维指针数组(int **&二维指针数组,int 行,int 列){
		二维指针数组=new int*[行];
		for(int i=0;i<行;i++){
			二维指针数组[i]=new int[列];
			for(int j=0;j<列;j++){
				二维指针数组[i][j]=0;
			}
		}
	}

	static void 创建二维指针数组(double **&二维指针数组,int 行,int 列){
		二维指针数组=new double*[行];
		for(int i=0;i<行;i++){
			二维指针数组[i]=new double[列];
			for(int j=0;j<列;j++){
				二维指针数组[i][j]=0;
			}
		}
	}

	static void 创建一维指针数组(int *&一维指针数组,int 列){
		一维指针数组=new int[列];
		for(int i=0;i<列;i++){
			一维指针数组[i]=0;
		}
	}
};

class 文件读写 { //要保证行列相等，没有的用0代替，所以所有编号都是从1开始
public:
    static void 写入二维矩阵(int **矩阵, int 行数, int 列数, string 地址,string 分隔符="\t"){
        FILE *fp=fopen(地址.c_str(),"w");
        for(int i=0;i<行数;i++) {
            fprintf(fp,"%d",矩阵[i][0]);
            for (int j = 1; j < 列数; j++)
                fprintf(fp,"%s%d", 分隔符.c_str() ,矩阵[i][j]);
            fprintf(fp,"\r\n");
        }
        fclose(fp);
    }

    static void 写入二维矩阵(double **矩阵, int 行数, int 列数, string 地址,string 分隔符="\t"){
        FILE *fp=fopen(地址.c_str(),"w");
        for(int i=0;i<行数;i++) {
            fprintf(fp,"%e",矩阵[i][0]);
            for (int j = 1; j < 列数; j++)
                fprintf(fp,"%s%e", 分隔符.c_str() ,矩阵[i][j]);
            fprintf(fp,"\r\n");
        }
        fclose(fp);
    }

    static void 读取二维矩阵(int **矩阵,int 行数, int 列数, string 地址){
        FILE *fp=fopen(地址.c_str(),"r");
        for(int i=0;i<行数;i++) {
            for (int j = 0; j < 列数; j++) {
                fscanf(fp, "%d", &矩阵[i][j]);
            }
        }
        fclose(fp);
    }

    static void 读取二维矩阵(double **矩阵,int 行数, int 列数, string 地址){
        FILE *fp=fopen(地址.c_str(),"r");
        for(int i=0;i<行数;i++) {
            for (int j = 0; j < 列数; j++)
                fscanf(fp,"%lf" ,&矩阵[i][j]);
        }
        fclose(fp);
    }
};

class LDA_基础{
public:
	int **MN_词号=NULL;
	int **MN_主题=NULL;
	int **VK_词数=NULL;
	int *K_词数=NULL;
	int **MK_词数=NULL;
	int *M_词数=NULL;

    string MN_主题默认输出地址="MN_主题.txt";
    string VK_词数默认输出地址="VK_词数.txt";
    string MK_词数默认输出地址="MK_词数.txt";

	int 主题数K;
	int 词总数V;
	int 文档数M;
	int 最大文档长度N;
	double alpha;
	double beta;

	int 迭代到第几次;
	int 采样到第几篇文档;
	
	static int 累计算法(double 主题概率[],int 主题个数=2147483647){
		srand((unsigned)time(NULL));
		rand();
		double 累计概率=0,骰子概率,总概率=0;
        for(int i=0;i<主题个数;i++){
            总概率+=主题概率[i];
        }
        骰子概率=rand()*1.0/RAND_MAX*总概率;
		int 主题=1;
		for(int i=0;i<主题个数;i++){
			if(累计概率<=骰子概率 && (累计概率+=主题概率[i],累计概率)>=骰子概率)
				break;
			主题++;
		}
		return 主题;
	}
	
	static void 还原矩阵参数(int **MN_词号_原,int **MN_主题_原,int **VK_词数_待,int *K_词数_待,
					   int **MK_词数_待,int *M_词数_待,int 主题数K_原,int 文档数M_原,int 最大文档长度N_原){
		int 主题,第几文档,文档第几词;
		for(第几文档=0;第几文档<文档数M_原;第几文档++){
			for(文档第几词=0;文档第几词<最大文档长度N_原;文档第几词++){
				if(MN_词号_原[第几文档][文档第几词]<=0)
					break;
				主题=MN_主题_原[第几文档][文档第几词];
				VK_词数_待[MN_词号_原[第几文档][文档第几词]-1][主题-1]++;
				K_词数_待[主题-1]++;
				MK_词数_待[第几文档][主题-1]++;
				M_词数_待[第几文档]++;
			}
		}
	}

	static void 计算分布(int **矩阵_词数,int 行,int 列,double 狄利克雷参数,bool 以行为累计,double **&分布结果){
		if(矩阵_词数==NULL){
			printf("变量不能为空！");
			return;
		}
		double 一组累计;
		// 分布结果 永远是按行累计的，比如主题-词分布，所以 矩阵_词数 和 分布结果 可能是转置的关系。
		if(以行为累计){
			if(分布结果==NULL)
				数据处理::创建二维指针数组(分布结果,行,列);
			for(int i=0;i<行;i++){
				一组累计=0;
				for(int j=0;j<列;j++)
					一组累计+=矩阵_词数[i][j];
				for(int j=0;j<列;j++)
					分布结果[i][j]=(矩阵_词数[i][j]+狄利克雷参数)/(一组累计+行*狄利克雷参数);
			}
		}else{
			if(分布结果==NULL)
				数据处理::创建二维指针数组(分布结果,列,行);
			for(int i=0;i<列;i++){
				一组累计=0;
				for(int j=0;j<行;j++)
					一组累计+=矩阵_词数[j][i];
				for(int j=0;j<列;j++)
					分布结果[i][j]=(矩阵_词数[j][i]+狄利克雷参数)/(一组累计+列*狄利克雷参数);
			}
		}
	}

	static void 随机化矩阵参数(int **MN_词号_原,int **MN_主题_待,int **VK_词数_待,int *K_词数_待,
						int **MK_词数_待,int *M_词数_待,int 主题数K_原,int 文档数M_原,int 最大文档长度N_原){
		srand((unsigned)time(NULL));
		rand();
		int 主题,第几文档,文档第几词;
		for(第几文档=0;第几文档<文档数M_原;第几文档++){
			for(文档第几词=0;文档第几词<最大文档长度N_原;文档第几词++){
				if(MN_词号_原[第几文档][文档第几词]<=0)
					break;
				主题=rand()/主题数K_原+1;
				MN_主题_待[第几文档][文档第几词]=主题;
				VK_词数_待[MN_词号_原[第几文档][文档第几词]-1][主题-1]++;
				K_词数_待[主题-1]++;
				MK_词数_待[第几文档][主题-1]++;
				M_词数_待[第几文档]++;
			}
		}
	}
	
	void 初始内部变量设定(int **MN_词号_原,int 主题数K_原,int 词总数V_原,int 文档数M_原,int 最大文档长度N_原,
				  int **MN_主题_原=NULL){ //MN_主题_原 不等于空则可以接着后面训练
		this->MN_词号=MN_词号_原;
		数据处理::创建二维指针数组(this->VK_词数,词总数V_原,主题数K_原);
		数据处理::创建二维指针数组(this->MK_词数,文档数M_原,主题数K_原);
		数据处理::创建一维指针数组(this->K_词数,主题数K_原);
		数据处理::创建一维指针数组(this->M_词数,文档数M_原);
		this->主题数K=主题数K_原;
		this->词总数V=词总数V_原;
		this->文档数M=文档数M_原;
		this->最大文档长度N=最大文档长度N_原;
		this->alpha=50.0/主题数K_原;
		this->beta=0.01;
		if(MN_主题_原==NULL){
			数据处理::创建二维指针数组(this->MN_主题,文档数M_原,最大文档长度N_原);
			LDA_基础::随机化矩阵参数(this->MN_词号,this->MN_主题,this->VK_词数,this->K_词数,this->MK_词数,
							  this->M_词数,主题数K_原,文档数M_原,最大文档长度N_原);
		}else{
			this->MN_主题=MN_主题_原;
			LDA_基础::还原矩阵参数(this->MN_词号,this->MN_主题,this->VK_词数,this->K_词数,this->MK_词数,
							  this->M_词数,主题数K_原,文档数M_原,最大文档长度N_原);
		}
	}

    void 主要文件输出(string MN_主题输出地址="",string VK_词数输出地址="",string MK_词数输出地址=""){
        if(MN_主题输出地址==""){
            MN_主题输出地址=this->MN_主题默认输出地址;
        }
        if(VK_词数输出地址==""){
            VK_词数输出地址=this->VK_词数默认输出地址;
        }
        if(MK_词数输出地址==""){
            MK_词数输出地址=this->MK_词数默认输出地址;
        }
        文件读写::写入二维矩阵(this->MN_主题,this->文档数M,this->最大文档长度N,MN_主题输出地址);
        文件读写::写入二维矩阵(this->VK_词数,this->词总数V,this->主题数K,VK_词数输出地址);
        文件读写::写入二维矩阵(this->MK_词数,this->文档数M,this->主题数K,MK_词数输出地址);
    }

	~LDA_基础(){
		数据处理::释放二维指针数组内存(this->VK_词数,this->词总数V);
		数据处理::释放二维指针数组内存(this->MK_词数,this->文档数M);
		数据处理::释放二维指针数组内存(this->MN_主题,this->文档数M);
		数据处理::释放二维指针数组内存(this->MN_词号,this->文档数M);
		delete[] this->K_词数;
		delete[] this->M_词数;
	}
};

class 吉布斯采样_训练:public LDA_基础{
public:
	double **主题词分布=NULL;
	double **文档主题分布=NULL;
	double **词主题分布=NULL;
	double **主题文档分布=NULL;

    string 主题词分布_默认输出地址="主题词分布.txt";
    string 文档主题分布_默认输出地址="文档主题分布.txt";

	void 开始采样(int 迭代次数,double 阿尔法=-1,double 贝塔=-1){
		if(this->MN_主题==NULL){
			printf("缺少初始化！");
			return;
		}
		int 主题,新主题,词号;
		double 主题概率[this->主题数K];
		if(阿尔法>0) this->alpha=阿尔法;
		if(贝塔>0) this->beta=贝塔;
		for(this->迭代到第几次=0;this->迭代到第几次<迭代次数;this->迭代到第几次++){
			//第几篇、第几词、第几主题都是减一后的编号
			for(this->采样到第几篇文档=0;this->采样到第几篇文档<this->文档数M;this->采样到第几篇文档++){
				for(int 第几词=0;第几词<this->最大文档长度N;第几词++){
					词号=this->MN_词号[this->采样到第几篇文档][第几词];
					if(词号<=0)
						break;
					主题=this->MN_主题[this->采样到第几篇文档][第几词];
					this->K_词数[主题-1]--;
					this->MK_词数[this->采样到第几篇文档][主题-1]--;
					this->VK_词数[词号-1][主题-1]--;
					for(int 第几主题=0;第几主题<this->主题数K;第几主题++){
						主题概率[第几主题]=((this->VK_词数[词号-1][第几主题]+this->beta)/
								(this->K_词数[第几主题]+this->词总数V*this->beta))*
								((this->MK_词数[this->采样到第几篇文档][第几主题]+this->alpha)/
								(this->M_词数[this->采样到第几篇文档]+this->主题数K*this->alpha));
					}
					新主题=累计算法(主题概率,this->主题数K);
					this->MN_主题[this->采样到第几篇文档][第几词]=新主题;
					this->K_词数[新主题-1]++;
					this->MK_词数[this->采样到第几篇文档][新主题-1]++;
					this->VK_词数[词号-1][新主题-1]++;
				}
			}
		}
	}

	double** 计算主题词分布(){
		if(this->VK_词数==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->VK_词数,this->词总数V,this->主题数K,this->beta,false,this->主题词分布);
		return this->主题词分布;
	}

	double** 计算词主题分布(){
		if(this->VK_词数==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->VK_词数,this->词总数V,this->主题数K,this->beta,true,this->词主题分布);
		return this->词主题分布;
	}

	double** 计算文档主题分布(){
		if(this->MK_词数==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->MK_词数,this->文档数M,this->主题数K,this->alpha,true,this->文档主题分布);
		return this->文档主题分布;
	}

	double** 计算主题文档分布(){
		if(this->MK_词数==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->MK_词数,this->文档数M,this->主题数K,this->alpha,false,this->主题文档分布);
		return this->主题文档分布;
	}

    void 主要文件输出(string MN_主题输出地址="",string VK_词数输出地址="",string MK_词数输出地址=""
            ,string 主题词分布_输出地址="",string 文档主题分布_输出地址=""){
        LDA_基础::主要文件输出(MN_主题输出地址,VK_词数输出地址,MK_词数输出地址);
        if(主题词分布_输出地址==""){
            主题词分布_输出地址=this->主题词分布_默认输出地址;
        }
        if(文档主题分布_输出地址==""){
            文档主题分布_输出地址=this->文档主题分布_默认输出地址;
        }
        this->计算主题词分布();
        this->计算文档主题分布();
        文件读写::写入二维矩阵(this->主题词分布,this->主题数K,this->词总数V,主题词分布_输出地址);
        文件读写::写入二维矩阵(this->文档主题分布,this->文档数M,this->主题数K,文档主题分布_输出地址);
    }

	~吉布斯采样_训练(){
		数据处理::释放二维指针数组内存(this->主题词分布,this->主题数K);
		数据处理::释放二维指针数组内存(this->文档主题分布,this->文档数M);
		数据处理::释放二维指针数组内存(this->词主题分布,this->主题数K);
		数据处理::释放二维指针数组内存(this->主题文档分布,this->文档数M);
	}
};

class 吉布斯采样_预测:public LDA_基础{
public:
	int **MN_词号_新=NULL;
	int **MN_主题_新=NULL;
	int **VK_词数_新=NULL;
	int *K_词数_新=NULL;
	int **MK_词数_新=NULL;
	int *M_词数_新=NULL;

	double **主题词分布_新=NULL;
	double **文档主题分布_新=NULL;
	double **词主题分布_新=NULL;
	double **主题文档分布_新=NULL;

    string MN_主题_新默认输出地址="MN_主题_新.txt";
    string VK_词数_新默认输出地址="VK_词数_新.txt";
    string MK_词数_新默认输出地址="MK_词数_新.txt";
    string 主题词分布_新_默认输出地址="主题词分布_新.txt";
    string 文档主题分布_新_默认输出地址="文档主题分布_新.txt";

	int 词总数V_新; // 最小是 this->词总数V
	int 文档数M_新; // MN_词号_新 的文档数
	int 最大文档长度N_新; // MN_词号_新 的最大文档长度

	void 初始内部变量设定(int **MN_词号_原,int 主题数K_原,int 词总数V_原,int 文档数M_原,int 最大文档长度N_原,int **MN_主题_原,
				  int **MN_词号_新_原,int 词总数V_新_原,int 文档数M_新_原,int 最大文档长度N_新_原){
		if(MN_主题_原==NULL){
			printf("没有训练好的数据！");
			return;
		}
		LDA_基础::初始内部变量设定(MN_词号_原,主题数K_原,词总数V_原,文档数M_原,最大文档长度N_原,MN_主题_原);
		this->MN_词号_新=MN_词号_新_原;
		数据处理::创建二维指针数组(this->VK_词数_新,词总数V_新_原,主题数K_原);
		数据处理::创建二维指针数组(this->MK_词数_新,文档数M_新_原,主题数K_原);
		数据处理::创建一维指针数组(this->K_词数_新,主题数K_原);
		数据处理::创建一维指针数组(this->M_词数_新,文档数M_新_原);
		this->词总数V_新=词总数V_新_原;
		this->文档数M_新=文档数M_新_原;
		this->最大文档长度N_新=最大文档长度N_新_原;
		数据处理::创建二维指针数组(this->MN_主题_新,文档数M_新_原,最大文档长度N_新_原);
		LDA_基础::随机化矩阵参数(this->MN_词号_新,this->MN_主题_新,this->VK_词数_新,this->K_词数_新,
						  this->MK_词数_新,this->M_词数_新,主题数K_原,文档数M_新_原,最大文档长度N_新_原);
	}

    void 再指定新文档(int **MN_词号_新_原,int 词总数V_新_原,int 文档数M_新_原,int 最大文档长度N_新_原){
        //释放原空间
        数据处理::释放二维指针数组内存(this->VK_词数_新,this->词总数V_新);
        数据处理::释放二维指针数组内存(this->MK_词数_新,this->文档数M_新);
        数据处理::释放二维指针数组内存(this->MN_主题_新,this->文档数M_新);
        数据处理::释放二维指针数组内存(this->MN_词号_新,this->文档数M_新);
        delete[] this->K_词数_新;
        delete[] this->M_词数_新;
        数据处理::释放二维指针数组内存(this->主题词分布_新,this->主题数K);
        数据处理::释放二维指针数组内存(this->文档主题分布_新,this->文档数M_新);
        数据处理::释放二维指针数组内存(this->词主题分布_新,this->主题数K);
        数据处理::释放二维指针数组内存(this->主题文档分布_新,this->文档数M_新);
        //新文档相关赋值
        this->MN_词号_新=MN_词号_新_原;
        数据处理::创建二维指针数组(this->VK_词数_新,词总数V_新_原,this->主题数K);
        数据处理::创建二维指针数组(this->MK_词数_新,文档数M_新_原,this->主题数K);
        数据处理::创建一维指针数组(this->K_词数_新,this->主题数K);
        数据处理::创建一维指针数组(this->M_词数_新,文档数M_新_原);
        this->词总数V_新=词总数V_新_原;
        this->文档数M_新=文档数M_新_原;
        this->最大文档长度N_新=最大文档长度N_新_原;
        数据处理::创建二维指针数组(this->MN_主题_新,文档数M_新_原,最大文档长度N_新_原);
        LDA_基础::随机化矩阵参数(this->MN_词号_新,this->MN_主题_新,this->VK_词数_新,this->K_词数_新,
                        this->MK_词数_新,this->M_词数_新,this->主题数K,文档数M_新_原,最大文档长度N_新_原);
    }

	void 开始采样(int 迭代次数,double 阿尔法=-1,double 贝塔=-1){
		if(this->MN_主题==NULL || this->MN_主题_新==NULL){
			printf("缺少初始化！");
			return;
		}
		int 主题,新主题,词号,原VK对应值,原K对应值;
		bool 伪计数能力;
		double 主题概率[this->主题数K];
		if(阿尔法>0) this->alpha=阿尔法;
		if(贝塔>0) this->beta=贝塔;
		for(this->迭代到第几次=0;this->迭代到第几次<迭代次数;this->迭代到第几次++){
			//第几篇、第几词、第几主题都是减一后的编号
			for(this->采样到第几篇文档=0;this->采样到第几篇文档<this->文档数M_新;this->采样到第几篇文档++){
				for(int 第几词=0;第几词<this->最大文档长度N_新;第几词++){
					词号=this->MN_词号_新[this->采样到第几篇文档][第几词];
					if(词号<=0)
						break;
					主题=this->MN_主题_新[this->采样到第几篇文档][第几词];
					this->K_词数_新[主题-1]--;
					this->MK_词数_新[this->采样到第几篇文档][主题-1]--;
					this->VK_词数_新[词号-1][主题-1]--;
					伪计数能力=词号<=this->词总数V?true:false;
					for(int 第几主题=0;第几主题<this->主题数K;第几主题++){
						if(伪计数能力){
							原VK对应值=this->VK_词数[词号-1][第几主题];
							原K对应值=this->K_词数[第几主题];
						}else{
							原VK对应值=0;
							原K对应值=0;
						}
						主题概率[第几主题]=((this->VK_词数_新[词号-1][第几主题]+原VK对应值+this->beta)/
									(this->K_词数_新[第几主题]+原K对应值+this->词总数V_新*this->beta))*
								   ((this->MK_词数_新[this->采样到第几篇文档][第几主题]+this->alpha)/
									(this->M_词数_新[this->采样到第几篇文档]+this->主题数K*this->alpha));
					}
					新主题=累计算法(主题概率,this->主题数K);
					this->MN_主题_新[this->采样到第几篇文档][第几词]=新主题;
					this->K_词数_新[新主题-1]++;
					this->MK_词数_新[this->采样到第几篇文档][新主题-1]++;
					this->VK_词数_新[词号-1][新主题-1]++;
				}
			}
		}
	}

	double** 计算主题词分布(){
		if(this->VK_词数_新==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->VK_词数_新,this->词总数V_新,this->主题数K,this->beta,false,this->主题词分布_新);
		return this->主题词分布_新;
	}

	double** 计算词主题分布(){
		if(this->VK_词数_新==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->VK_词数_新,this->词总数V_新,this->主题数K,this->beta,true,this->词主题分布_新);
		return this->词主题分布_新;
	}

	double** 计算文档主题分布(){
		if(this->MK_词数_新==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->MK_词数_新,this->文档数M_新,this->主题数K,this->alpha,true,this->文档主题分布_新);
		return this->文档主题分布_新;
	}

	double** 计算主题文档分布(){
		if(this->MK_词数_新==NULL){
			printf("还未采样！");
			return NULL;
		}
		计算分布(this->MK_词数_新,this->文档数M_新,this->主题数K,this->alpha,false,this->主题文档分布_新);
		return this->主题文档分布_新;
	}

    void 主要文件输出(string MN_主题_新输出地址="",string VK_词数_新输出地址="",string MK_词数_新输出地址=""
            ,string 主题词分布_新_输出地址="",string 文档主题分布_新_输出地址=""){
        if(MN_主题_新输出地址==""){
            MN_主题_新输出地址=this->MN_主题_新默认输出地址;
        }
        if(VK_词数_新输出地址==""){
            VK_词数_新输出地址=this->VK_词数_新默认输出地址;
        }
        if(MK_词数_新输出地址==""){
            MK_词数_新输出地址=this->MK_词数_新默认输出地址;
        }
        if(主题词分布_新_输出地址==""){
            主题词分布_新_输出地址=this->主题词分布_新_默认输出地址;
        }
        if(文档主题分布_新_输出地址==""){
            文档主题分布_新_输出地址=this->文档主题分布_新_默认输出地址;
        }
        this->计算主题词分布();
        this->计算文档主题分布();
        文件读写::写入二维矩阵(this->MN_主题_新,this->文档数M_新,this->最大文档长度N_新,MN_主题_新输出地址);
        文件读写::写入二维矩阵(this->VK_词数_新,this->词总数V_新,this->主题数K,VK_词数_新输出地址);
        文件读写::写入二维矩阵(this->MK_词数_新,this->文档数M_新,this->主题数K,MK_词数_新输出地址);
        文件读写::写入二维矩阵(this->主题词分布_新,this->主题数K,this->词总数V_新,主题词分布_新_输出地址);
        文件读写::写入二维矩阵(this->文档主题分布_新,this->文档数M_新,this->主题数K,文档主题分布_新_输出地址);
    }

	~吉布斯采样_预测(){
		数据处理::释放二维指针数组内存(this->VK_词数_新,this->词总数V_新);
		数据处理::释放二维指针数组内存(this->MK_词数_新,this->文档数M_新);
		数据处理::释放二维指针数组内存(this->MN_主题_新,this->文档数M_新);
		数据处理::释放二维指针数组内存(this->MN_词号_新,this->文档数M_新);
		delete[] this->K_词数_新;
		delete[] this->M_词数_新;
		数据处理::释放二维指针数组内存(this->主题词分布_新,this->主题数K);
		数据处理::释放二维指针数组内存(this->文档主题分布_新,this->文档数M_新);
		数据处理::释放二维指针数组内存(this->词主题分布_新,this->主题数K);
		数据处理::释放二维指针数组内存(this->主题文档分布_新,this->文档数M_新);
	}
};

吉布斯采样_训练 *吉布斯采样_训练实例=new 吉布斯采样_训练();
吉布斯采样_预测 *吉布斯采样_预测实例=new 吉布斯采样_预测();

//python调用方法
extern "C" void LDA训练并写入文件(string MN_词号地址,int 主题数K_原,int 词总数V_原,int 文档数M_原,int 最大文档长度N_原
        ,int 迭代次数,double 阿尔法=-1,double 贝塔=-1,string MN_主题地址=""){
    int **MN_词号=NULL,**MN_主题=NULL;
    数据处理::创建二维指针数组(MN_词号,文档数M_原,最大文档长度N_原);
    文件读写::读取二维矩阵(MN_词号,文档数M_原,最大文档长度N_原,MN_词号地址);
    if(MN_主题地址==""){
        吉布斯采样_训练实例->初始内部变量设定(MN_词号,主题数K_原,词总数V_原,文档数M_原,最大文档长度N_原);
    }else{
        数据处理::创建二维指针数组(MN_主题,文档数M_原,最大文档长度N_原);
        文件读写::读取二维矩阵(MN_主题,文档数M_原,最大文档长度N_原,MN_主题地址);
        吉布斯采样_训练实例->初始内部变量设定(MN_词号,主题数K_原,词总数V_原,文档数M_原,最大文档长度N_原,MN_主题);
    }
    吉布斯采样_训练实例->开始采样(迭代次数,阿尔法,贝塔);
    吉布斯采样_训练实例->主要文件输出();
}

extern "C" void LDA预测并写入文件(string MN_词号地址,int 主题数K_原,int 词总数V_原,int 文档数M_原,int 最大文档长度N_原
        ,string MN_主题地址,string MN_词号_新地址,int 词总数V_新_原,int 文档数M_新_原,int 最大文档长度N_新_原
        ,int 迭代次数,double 阿尔法=-1,double 贝塔=-1){
    int **MN_词号=NULL,**MN_主题=NULL,**MN_词号_新=NULL;
    数据处理::创建二维指针数组(MN_词号,文档数M_原,最大文档长度N_原);
    文件读写::读取二维矩阵(MN_词号,文档数M_原,最大文档长度N_原,MN_词号地址);
    数据处理::创建二维指针数组(MN_主题,文档数M_原,最大文档长度N_原);
    文件读写::读取二维矩阵(MN_主题,文档数M_原,最大文档长度N_原,MN_主题地址);
    数据处理::创建二维指针数组(MN_词号_新,文档数M_新_原,最大文档长度N_新_原);
    文件读写::读取二维矩阵(MN_词号_新,文档数M_新_原,最大文档长度N_新_原,MN_词号_新地址);
    吉布斯采样_预测实例->初始内部变量设定(MN_词号,主题数K_原,词总数V_原,文档数M_原,最大文档长度N_原,MN_主题,
            MN_词号_新,词总数V_新_原,文档数M_新_原,最大文档长度N_新_原);
    吉布斯采样_预测实例->开始采样(迭代次数,阿尔法,贝塔);
    吉布斯采样_预测实例->主要文件输出();
}

int main(int argc, char *argv[]) {
    LDA训练并写入文件("MN_词号.txt"
            ,10,13749,9,4890,1000,-1,-1
            ,"MN_主题.txt");
}