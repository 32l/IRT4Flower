// Time-stamp: <2012-01-26 16:33:24 tenmo>

#include <stdio.h>
#include <stdlib.h>
#include <highgui.h>
#include <cv.h>
#include <cxcore.h>
#include <ml.h>

#ifdef _DEBUG
//Debugモードの場合
#pragma comment(lib,"cv200d.lib")
#pragma comment(lib,"cxcore200d.lib")
#pragma comment(lib,"cvaux200d.lib")
#pragma comment(lib,"highgui200d.lib")
#pragma comment(lib,"ml200d.lib")
#else
//Releaseモードの場合
#pragma comment(lib,"cv200.lib")
#pragma comment(lib,"cxcore200.lib")
#pragma comment(lib,"cvaux200.lib")
#pragma comment(lib,"highgui200.lib")
#pragma comment(lib,"ml200.lib")
#endif

static void read_file(char* file_name, int num_samples, int num_features, CvMat* train_data, CvMat* responses)
{
    FILE* fp = fopen(file_name, "r");
    if (fp == NULL) {
        perror(file_name);
        exit(1);
    }
    for (int i = 0; i < num_samples; i++) {
        double tmp;
        int code = fscanf(fp, "%lf", &tmp);
	if (code != 1) {
	    fprintf(stderr, "%s: Bad format\n", file_name);
	    exit(1);
	}
        cvSetReal1D(responses, i, tmp);
        for (int j = 0; j < num_features; j++) {
            fscanf(fp, "%lf", &tmp);
            cvSetReal2D(train_data, i, j, tmp);
        }
    }
    fclose(fp);
}

int main(int argc, char* argv[])
{
    char* file_name = "C:/flower0126/test.txt";
    int num_classes = 2; //クラスの数
    int num_samples = 5; //学習パターンの数
    int num_features = 4; //特徴の数
    
    CvMat* train_data = cvCreateMat(num_samples, num_features, CV_32FC1);
    CvMat* responses = cvCreateMat(num_samples, 1, CV_32FC1);
    read_file(file_name, num_samples, num_features, train_data, responses);
    
    // ---- k nearest neighbor ----
    CvKNearest knn;
	knn.train(train_data, responses);
    
    int miss = 0;
    for (int i = 0; i < num_samples; i++) {
        printf("%.2lf => ", cvGetReal1D(responses, i));
        CvMat sample;
        cvGetRow(train_data, &sample, i);
		printf("1nn %.2lf\n", knn.find_nearest(&sample, 3));
        if (cvGetReal1D(responses, i) != knn.find_nearest(&sample, 1)) {
	    miss++;
		}
    }
    
    printf("correct rate 1nn %5.2lf\n", (double)(num_samples - miss / (double)num_samples * 100.0));
	getchar();
    exit(0);
}