//赤色を検出するプログラム
/* 構造体操作の詳細は　http://homepage3.nifty.com/mmgames/c_guide/16-02.html　*/

#include <stdio.h>
#include <stdlib.h>
#include <highgui.h>
#include <cv.h>
#include <cxcore.h>
#include "Labeling.h"
#include "LabelingW.h"
/* Labeling.hについて… http://oshiro.bpe.es.osaka-u.ac.jp/people/staff/imura/products/labeling */

#ifdef _DEBUG
//Debugモードの場合
#pragma comment(lib,"cv200d.lib")
#pragma comment(lib,"cxcore200d.lib")
#pragma comment(lib,"cvaux200d.lib")
#pragma comment(lib,"highgui200d.lib")
#else
//Releaseモードの場合
#pragma comment(lib,"cv200.lib")
#pragma comment(lib,"cxcore200.lib")
#pragma comment(lib,"cvaux200.lib")
#pragma comment(lib,"highgui200.lib")
#endif

#define ITERATIONS 2

//キャプチャ画像用IplImage
IplImage *img = NULL;
IplImage *imgResult, *imgH, *imgS, *imgV;

//色の平均値を格納する構造体
typedef struct{
	int h;
	int s;
	int v;
}AveColor;


typedef struct{
	double c_ratio;
	double d_ratio;
}Ratio;

//画像の本体や情報を格納する構造体
typedef struct{
	IplImage *bit;
	IplImage *convex;
	IplImage *label;
	char bitname[20];
	char convexname[20];
	CvScalar color;
}Para;

//
//  花領域を描画し、RGB値が近い部分も花領域に加える
//
//	引数:
//		skinImage       : 花の色抽出画像用IplImage
//		label           : ラベリングした結果
//		convexHullImage : ConvexHull画像用IplImage
//
void drawMaxArea(IplImage *bitImage, IplImage *label, IplImage *convexHullImage ) {
	unsigned char char_h, char_s, char_v;
	unsigned char ex_char_h_x, ex_char_s_x, ex_char_v_x, ex_char_h_y, ex_char_s_y, ex_char_v_y;
	int h, s, v;
	int ex_h_x, ex_s_x, ex_v_x, d_h_x, d_s_x, d_v_x;
	int ex_h_y, ex_s_y, ex_v_y, d_h_y, d_s_y, d_v_y;

	for(int x = 0; x < bitImage->width; x++ ) {
		for( int y = 0; y < bitImage->height; y++ ) {

			if(x < 1 || y < 1) {
				cvSet2D(convexHullImage, y, x, CV_RGB(0, 0, 0));
				cvSetReal2D( bitImage, y, x, 0 );
			} else {

				char_h = (unsigned char)imgH->imageDataOrigin[x+y*imgH->width];
				char_s = (unsigned char)imgS->imageDataOrigin[x+y*imgS->width];
				char_v = (unsigned char)imgV->imageDataOrigin[x+y*imgV->width];

				ex_char_h_x = (unsigned char)imgH->imageDataOrigin[(x-1)+y*imgH->width];
				ex_char_s_x = (unsigned char)imgS->imageDataOrigin[(x-1)+y*imgS->width];
				ex_char_v_x = (unsigned char)imgV->imageDataOrigin[(x-1)+y*imgV->width];

				ex_char_h_y = (unsigned char)imgH->imageDataOrigin[x+(y-1)*imgH->width];
				ex_char_s_y = (unsigned char)imgS->imageDataOrigin[x+(y-1)*imgS->width];
				ex_char_v_y = (unsigned char)imgV->imageDataOrigin[x+(y-1)*imgV->width];


				//それぞれ変数を作成
				//hsv本体
				h = char_h; s = char_s; v = char_v;
				//各hsvの前のxの色
				ex_h_x = ex_char_h_x; ex_s_x = ex_char_s_x; ex_v_x = ex_char_v_x;
				//各hsvの前のyの色
				ex_h_y = ex_char_h_y; ex_s_y = ex_char_s_y; ex_v_y = ex_char_v_y;
				//前のxの色との差異　difference
				d_h_x = h - ex_h_x; d_s_x = s - ex_s_x; d_v_x = v - ex_v_x;
				//前のyの色との差異
				d_h_y = h - ex_h_y; d_s_y = s - ex_s_y; d_v_y = v - ex_v_y;


				if( cvGetReal2D( label, y, x ) == 1) {
					//	最大領域だった場合 面積をインクリメントし、convexを白色にセット
					cvSet2D( convexHullImage, y, x, CV_RGB( 255, 255, 255 ) );
					cvSetReal2D(label, y, x, 1);			//ラベルを付加する(1が一番大きい画像)
					cvSetReal2D(bitImage, y, x, 255);
				} else {
					//でなかった場合、黒にセット
					cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 0 ) );
					cvSetReal2D( bitImage, y, x, 0 );
					//もし近接のRGBが近い値であれば、補正
					//デフォルトの値　abs(d_h_x) <= 2  && abs(d_s_x) <= 13 && abs(d_v_x) <= 26
					if(cvGetReal2D( label, y, x-1 ) == 1 && abs(d_h_x) <= 2 && abs(d_s_x) <= 13 && abs(d_v_x) <= 26) {
						cvSet2D( convexHullImage, y, x, CV_RGB( 255, 0, 0));		//convexを赤で描画
						cvSetReal2D(label, y, x, 1);			//ラベルを付加する(1が一番大きい画像)
						cvSetReal2D(bitImage, y, x, 255);		//bitImageに255を付加して白を描く

					//y方向を一旦お休み　復活させるときは少し近似値を下げる？
					//} else if(cvGetReal2D( label, y-1, x ) == 1 && abs(d_h_y) < 2 && abs(d_s_y) < 6.75 && abs(d_v_y) < 12.5) {
					//	cvSet2D( convexHullImage, y, x, CV_RGB( 0, 255, 0 ));		//convexを緑で描画
					//	cvSetReal2D(label, y, x, 1);			//ラベルを付加する(1が一番大きい画像)
					//	cvSetReal2D(bitImage, y, x, 255);		//bitImageに255を付加して白を描く

					//} else if(cvGetReal2D( label, y, x-1 ) == 1 && s <= 51 && abs(d_s_x) <= 13 && abs(d_v_x) <= 26 ) {
					//	cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 255 ));		//convexを青で描画
					//	cvSetReal2D(label, y, x, 1);			//ラベルを付加する(1が一番大きい画像)
					//	cvSetReal2D(bitImage, y, x, 255);		//bitImageに255を付加して白を描く
					}

				}
			}
		}
	}

}

//
//  色の近似による花領域をカウントする際に誤カウントしている部分を修正する
//
//	引数:
//		skinImage       : 花の色抽出画像用IplImage
//		label           : ラベリングした結果
//		convexHullImage : ConvexHull画像用IplImage
//
void reviseMaxArea(IplImage *bitImage, IplImage *label, IplImage *convexHullImage ) {

	CvScalar color;

	int count = 0, sequence = 0,same = 0;
	int mem = -1;
	int ybuf = bitImage->height;

	//for文開始 
	for( int y = 0; y < bitImage->height; y++ ) {
		for(int x = 0; x < bitImage->width; x++ ) {
			//もしx軸方向に赤色の線(xの補正)が激しく連続していたらカウント
			color = cvGet2D(convexHullImage, y, x);
			if(color.val[0] == 0.0 && color.val[1] == 0.0 && color.val[2] == 255.0) { //BGR
				count++;
				sequence++;
			} else {
				//シーケンスは再度リセット
				sequence = 0;
			}

			if(count > bitImage->width/5) {
				//もし赤の合計数が画像の1/5を突破していたら変数sameをインクリメントし次の行へ
				if(ybuf > y) ybuf = y;
				same++;
				y++;
			} else if(sequence > bitImage->width/8) {
				//もし赤の長さが画像の1/8連続していたらその場で修正する
				for(int m = 0; m < bitImage->width; m++){
					color = cvGet2D(convexHullImage, y, m);
					//convexの色が赤ならば
					if(color.val[0] != 255.0 || color.val[1] != 255.0 || color.val[2] != 255.0) {
						//convexの白い部分 (== 補正をかける前のもともとの花の部分)　以外を全て塗りつぶす
						cvSetReal2D( bitImage, y, m, 0 );
						cvSet2D( convexHullImage, y, m, CV_RGB( 0, 0, 0 ) );
						cvSetReal2D(label, y, m, 0);			//ラベルを付加する(1が一番大きい画像)
					}
				} //for文m終了
				y++;
			} //if文終了
		}//for文x終了
		count = 0;
		sequence = 0;
	} //for文y終了


	//もしsameが10以上あれば、大きすぎるため大幅修正を行う(リセットでもいいのでは？)
	if(same > 10) {
		printf("xの修正が連続しす(same=%d)ぎたため、y=%dから大幅修正\n", same, ybuf);
		//エラー部分を消去して終了
		for( int y = ybuf -1 ; y < bitImage->height ; y++ ) {
			for(int x = 0; x < bitImage->width ; x++ ) {
				color = cvGet2D(convexHullImage, y, x);
				//どれか一つでも255じゃないものがあれば… (== 255 255 255 つまり白でなければ)
				if(color.val[0] != 255.0 || color.val[1] != 255.0 || color.val[2] != 255.0) {
					//convexの白い部分 (== 補正をかける前のもともとの花の部分)　以外を全て塗りつぶす
					//printf("%lf %lf %lf\n", color.val[0], color.val[1], color.val[2]);
					cvSetReal2D( bitImage, y, x, 0 );
					cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 0 ) );
					cvSetReal2D(label, y, x, 0);			//ラベルを付加する(1が一番大きい画像)
				}
			}
		}
	}

	cvDestroyWindow("conv");
}


//
//	花の最大領域の抽出(label==1)を行い、面積をカウントする
//
//	引数:
//		skinImage       : 花の色抽出画像用IplImage
//		label           : ラベリングした結果
//
//	戻り値:
//		花領域の面積
//
int pickupMaxArea(IplImage *bitImage, IplImage *label) {

	int flowerarea = 0;	//	花領域の面積

	for(int x = 0; x < bitImage->width; x++ ) {
		for( int y = 0; y < bitImage->height; y++ ) {
			if(cvGetReal2D(label, y, x) == 1) {
				flowerarea++;
			}
		}
	}

	return flowerarea;
}


//
//	ConvexHullを生成する
//
//	引数:
//		skinImage   : 肌色抽出画像用IplImage
//		flowerarea  : 花領域の面積(点の数)
//		flowerpoint : 花領域内の点の座標配列へのポインタ
//		hull        : ConvexHullの頂点のhandpointにおけるindex番号へのポインタ
//		pointMatrix : 花領域用行列へのポインタ
//		hullMatrix  : ConvexHull用行列へのポインタ
//
void createConvexHull(IplImage *bitImage, IplImage *label,int flowerarea, CvPoint **flowerpoint, int **hull,
					  CvMat *pointMatrix, CvMat *hullMatrix ) {
						  int i=0;

						  //	ConvexHullを計算するために必要な行列を生成する
						  *flowerpoint=( CvPoint * )malloc( sizeof( CvPoint ) * flowerarea);
						  *hull = ( int * )malloc( sizeof( int ) * flowerarea );
						  *pointMatrix = cvMat( 1, flowerarea, CV_32SC2, *flowerpoint );
						  *hullMatrix = cvMat( 1, flowerarea, CV_32SC1, *hull );

						  for( int x = 0; x < bitImage->width; x++ ) {
							  for(  int y = 0; y < bitImage->height; y++ ) {
								  if(cvGetReal2D(label, y, x) == 1) {
									  ( *flowerpoint )[i].x = x;
									  ( *flowerpoint )[i].y = y;
									  i++;
								  }
								} 
						  }

						  //	ConvexHullを生成する
						  cvConvexHull2( pointMatrix, hullMatrix, CV_CLOCKWISE, 0 );

}


//
//	ConvexHullを描画する
//
//	引数:
//		convexHullImage : ConvexHull画像用IplImage
//		flowerpoint     : 花領域内の点の座標配列
//		hull            : ConvexHullの頂点のhandpointにおけるindex番号
//		hullcount       : ConvexHullの頂点の数
//
void drawConvexHull(IplImage *convexHullImage, CvPoint *flowerpoint, int *hull, int hullcount ) {
	CvPoint pt0 = flowerpoint[hull[hullcount-1]];

	for( int i = 0; i < hullcount; i++ ) {
		CvPoint pt = flowerpoint[hull[i]];
		cvLine( convexHullImage, pt0, pt, CV_RGB(255, 255, 0));
		pt0 = pt;
	}
}

//
//	ConvexHull内の面積を求める
//
//	引数:
//		convexHullImage : ConvexHull画像用IplImage
//		flowerpoint       : 花領域内の点の座標配列
//		hull            : ConvexHullの頂点のhandpointにおけるindex番号
//		hullcount       : ConvexHullの頂点の数　　
//
//	戻り値:
//		ConvexHull内の面積
//
int calcConvexHullArea( IplImage *convexHullImage, CvPoint *flowerpoint, int *hull, int hullcount ) {

	//	ConvexHullの頂点からなる行列を生成
	CvPoint *hullpoint = ( CvPoint * )malloc( sizeof( CvPoint ) * hullcount );
	CvMat hMatrix = cvMat( 1, hullcount, CV_32SC2, hullpoint );
	for( int i = 0; i < hullcount; i++ ) {
		hullpoint[i]=flowerpoint[hull[i]];
	}

	//	ConvexHull内の点の数を数える
	int hullarea = 0;
	for( int x = 0; x < convexHullImage->width; x++ ) {
		for( int y = 0;y < convexHullImage->height; y++ ) {
			if( cvPointPolygonTest( &hMatrix, cvPoint2D32f( x, y ), 0 ) > 0) {
				hullarea++;
			}
		}
	}

	free( hullpoint );
	return hullarea;
}

//2点間の距離を求める
int Distance(CvPoint pt1, CvPoint pt2){
	int dis;
	
	dis = int(sqrt(pow(pt2.x - pt1.x, 2.0) + pow(pt2.y - pt1.y, 2.0)));

	return dis;
}

//
//画像の重心を求める関数
//
//	引数:
//		para		    : 画像の構造体
//		flowerpoint     : 花領域内の点の座標配列
//		hull            : ConvexHullの頂点のhandpointにおけるindex番号
//		hullcount       : ConvexHullの頂点の数　　
//
void Gravity(IplImage *bitImage, IplImage *convexHullImage, CvPoint *flowerpoint, int *hull, int hullcount, Ratio *ratiobox ) {

	CvMoments moments;
	CvPoint pt0 = flowerpoint[hull[hullcount-1]];
	CvPoint ptxy;
	int dis, dis_def;
	double sum = 0.0;

	cvMoments(bitImage, &moments, 0);
	double m00 = cvGetSpatialMoment(&moments, 0, 0);
	double m10 = cvGetSpatialMoment(&moments, 1, 0);
	double m01 = cvGetSpatialMoment(&moments, 0, 1);

	double gX = fabs(m10)/fabs(m00);
	double gY = fabs(m01)/fabs(m00);
	ptxy = cvPoint((int)gX, (int)gY);
	cvCircle(convexHullImage, ptxy, 4, CV_RGB(0, 0, 255), 1, 4, 0);

	//printf("gX = %.3f gY = %.3f \n", gX, gY);

	for( int i = 0; i < hullcount; i++ ) {
		if(i == 0) {
			//1本目のみ赤…基準線とする
			dis_def = Distance(ptxy, pt0);
			cvLine( convexHullImage, ptxy, pt0, CV_RGB(255, 0, 0));
		} else {
			cvLine( convexHullImage, ptxy, pt0, CV_RGB(0, 0, 255));
		}
		dis = Distance(ptxy, pt0);
		//printf("dis%d = %3.3f \n", i, (double)dis/dis_def);
		sum += (double)dis/dis_def;

		CvPoint pt = flowerpoint[hull[i]];
		pt0 = pt;
	}

	//各距離の比率の平均値を出す(試験運用中)
	ratiobox->d_ratio = sum / (double)(hullcount-1);
}

//
//	ラベリング部分の色の平均値を求める
//
//	引数:
//		skinImage       : 花の色抽出画像用IplImage
//		para			: 画像の構造体
//		aveColor		: 花の色の平均値を格納する構造体
//
void averageColor(Para *para, AveColor *aveColor) {

	unsigned char char_h, char_s, char_v;
	int h, s, v;
	int h_sum = 0, s_sum = 0, v_sum = 0, in = 0;

	for(int x = 0; x < para->bit->width; x++ ) {
		for( int y=0; y < para->bit->height; y++ ) {
			if( cvGetReal2D( para->label, y, x ) == 1) {
				//	最大領域だった場合 … convexで白だった場合
				char_h = (unsigned char)imgH->imageDataOrigin[x+y*imgH->width];
				char_s = (unsigned char)imgS->imageDataOrigin[x+y*imgS->width];
				char_v = (unsigned char)imgV->imageDataOrigin[x+y*imgV->width];

				h = char_h;
				s = char_s;
				v = char_v;

				h_sum += h;
				s_sum += s;
				v_sum += v;
				in++;
			}
		}		
	}
	aveColor->h = h_sum/in;	//後ろの数字は正規化
	aveColor->s = s_sum/in;
	aveColor->v = v_sum/in;
	printf("Average H: %d(/180) S: %d(/100) V: %d(/100) \n", aveColor->h * 2, (int)aveColor->s * 100 / 255, (int)aveColor->v * 100 / 255);
}

//
//	欠損領域を補間する
//
//	引数:
//		skinImage : 肌色抽出画像用IplImage
//		temp      : 一時保存用IplImage
//
void interpolate( IplImage *skinImage, IplImage *temp, int num) {
	//膨張をnum回行う
	cvDilate( skinImage, temp, NULL, num );

	//収縮をnum回行う
	cvErode( temp, skinImage, NULL, num );
}


//--------------------------------------------------------------- 
//【関数名　】：cv_ColorExtraction 
//【処理概要】：色抽出 
//【引数　　】：src_img        = 入力画像(8bit3ch) 
//　　　　　　：dst_img        = 出力2値画像(8bit1ch) 
//　　　　　　：mask_img        = 出力マスク画像(8bit3ch) 
//　　　　　　：code        = 色空間の指定（CV_BGR2HSV,CV_BGR2Labなど）
//　　　　　　：ch1_lower    = ch1のしきい値(小)
//　　　　　　：ch1_upper    = ch1のしきい値(大)
//　　　　　　：ch2_lower    = ch2のしきい値(小)
//　　　　　　：ch2_upper    = ch2のしきい値(大)
//　　　　　　：ch3_lower    = ch3のしきい値(小)
//　　　　　　：ch3_upper    = ch3のしきい値(大)
//【戻り値　】：なし 
//【備考　　】：lower <= upperの場合、lower以上upper以下の範囲を抽出、
//　　　　　　：lower >  upperの場合、upper以下lower以上の範囲を抽出します。
//---------------------------------------------------------------
void cv_ColorExtraction(IplImage* src_img, IplImage* dst_img, IplImage* mask_img,
                            int code,
                            int ch1_lower, int ch1_upper,
                            int ch2_lower, int ch2_upper,
                            int ch3_lower, int ch3_upper
                        ){

    int i, k;   

    IplImage *Color_img;
    IplImage *ch1_img, *ch2_img, *ch3_img;

    int lower[3];
    int upper[3];
    int val[3];

    CvMat *lut;   

    //codeに基づいたカラー変換
    Color_img = cvCreateImage(cvGetSize(src_img), src_img->depth, src_img->nChannels);
    cvCvtColor(src_img, Color_img, code);
       
    //3ChのLUT作成
    lut    = cvCreateMat(256, 1, CV_8UC3);

    lower[0] = ch1_lower;
    lower[1] = ch2_lower;
    lower[2] = ch3_lower;

    upper[0] = ch1_upper;
    upper[1] = ch2_upper;
    upper[2] = ch3_upper;

    for (i = 0; i < 256; i++){
        for (k = 0; k < 3; k++){
            if (lower[k] <= upper[k]){
                if ((lower[k] <= i) && (i <= upper[k])){
                    val[k] = 255;
                }else{
                    val[k] = 0;
                }
            }else{
                if ((i <= upper[k]) || (lower[k] <= i)){
                    val[k] = 255;
                }else{
                    val[k] = 0;
                }
            }
        }
        //LUTの設定
        cvSet1D(lut, i, cvScalar(val[0], val[1], val[2]));
    }

    //3ChごとのLUT変換（各チャンネルごとに２値化処理）
    cvLUT(Color_img, Color_img, lut);
    cvReleaseMat(&lut);

    //各チャンネルごとのIplImageを確保する
    ch1_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);
    ch2_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);
    ch3_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);

    //チャンネルごとに二値化された画像をそれぞれのチャンネルに分解する
    cvSplit(Color_img, ch1_img, ch2_img, ch3_img, NULL);

    //3Ch全てのANDを取り、マスク画像を作成する。

    cvAnd(ch1_img, ch2_img, dst_img);
    cvAnd(dst_img, ch3_img, dst_img);

    //入力画像(src_img)のマスク領域を出力画像(dst_img)へコピーする
    cvZero(mask_img);
    cvCopy(src_img, mask_img, dst_img);

    //解放
    cvReleaseImage(&Color_img);
    cvReleaseImage(&ch1_img);
    cvReleaseImage(&ch2_img);
    cvReleaseImage(&ch3_img);

}

//
//	上記の関数を用いて計算を行い、convexHull内の花の割合を出し、描画も行う関数
//
//	引数:
//		para :	画像の構造体
//
//	戻り値:
//		ratio :	convex内における花の割合
//

void Func(Para *para, AveColor *aveColor, Ratio *ratiobox) {

	double convexratio = -1.0;
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* find_contour = NULL;

	//収縮、膨張をして画像のノイズを除去	
	//interpolate(para->bit, para->bit, 2);

	//	ラベリングを行う 
	Label *labeling = createLabeling();
	exec( labeling, para->bit, para->label, true, 2000);
	/* exec(ラベリング変数, 入力画像, ラベリング結果を格納する変数, 
	領域が大きさの降順にソートされるか, 最小の領域サイズ) */

	if(getNumOfResultRegions(labeling) > 0 ) {
		//最小の領域サイズよりも大きな領域があった場合
		int flowerarea;		//	花領域の面積
		int hullarea;		//	ConvexHull内の面積
		int hullcount;		//	ConvexHullの頂点の数
		CvPoint *flowerpoint;	//	花領域内の点の座標配列
		int *hull;			//	ConvexHullの頂点のflowerpointにおけるindex番号
		CvMat pointMatrix;	//	花領域用行列
		CvMat hullMatrix;	//	ConvexHull用行列


		//輪郭を検出する
		int contour_num_color = cvFindContours(cvCloneImage(para->bit),
			storage, 
			&find_contour,
			sizeof(CvContour), 
			CV_RETR_EXTERNAL, 
			CV_CHAIN_APPROX_NONE, 
			cvPoint(0,0)
			);

		//printf("check-2 \n");

		//花領域を描画・近似値を補正
		drawMaxArea(para->bit, para->label, para->convex);
		//エラー箇所の修正
		reviseMaxArea(para->bit, para->label, para->convex);

		//bit-2値画像を縮小・圧縮
		interpolate(para->bit, para->bit, 2);
		//輪郭の表示
		cvDrawContours(para->convex, find_contour, para->color, para->color, 1, 2, 8, cvPoint(0,0)); //黄色の場合(2番目のyellow使用しない)]

		//printf("check-3 \n");

		//最大領域を抽出し、面積をカウント
		flowerarea = pickupMaxArea(para->bit, para->label);

		//printf("check-4 flowerarea = %d\n", flowerarea);

		//ConvexHullを生成
		createConvexHull(para->bit,para->label, flowerarea, &flowerpoint, &hull, &pointMatrix, &hullMatrix);
		hullcount = hullMatrix.cols;

		//printf("check-5 \n");

		//ConvexHullを描画
		drawConvexHull(para->convex, flowerpoint, hull, hullcount);

		//printf("check-6 \n");

		//ConvexHull内の面積を求める
		hullarea = calcConvexHullArea(para->convex, flowerpoint, hull, hullcount);

		//printf("check-7 \n");

		//最大領域の色の平均値を求める 構造体ポインタ渡し
		averageColor(para, aveColor);

		//printf("check-8 \n");

		//重心を求め、各頂点まで線を引く
		Gravity(para->bit, para->convex, flowerpoint, hull, hullcount, ratiobox);

		//printf("check-9 \n");

		//if(hullarea = 0) return ratio;
		ratiobox->c_ratio = flowerarea / ( double )hullarea;

		//printf("check-10 \n");

	}

	releaseLabeling(labeling);
	cvReleaseMemStorage(&storage);
}

//
//	各ウィンドウを作成し、画像表示の準備を行ってくれる関数
//
//	引数:
//		img	 :	元となる画像 キャプチャウィンドウ用に準備
//		para :	画像の構造体
//

void Show(Para *para){

	//引数の構造体のウィンドウを作成し
	cvNamedWindow(para->bitname, CV_WINDOW_AUTOSIZE);
	cvNamedWindow(para->convexname, CV_WINDOW_AUTOSIZE);

	//表示
	cvShowImage(para->bitname, para->bit);
	cvShowImage(para->convexname, para->convex);

}

//
//	各ウィンドウを閉じ、メモリの解放もを行ってくれる関数
//
//	引数:
//		img	 :	元となる画像 キャプチャウィンドウ用に準備
//		para :	画像の構造体
//

void Close(Para *para){

	//メモリの解放
	cvReleaseImage(&para->bit);
	cvReleaseImage(&para->convex);
    cvReleaseImage(&para->label);

	//ウィンドウを破棄
	cvDestroyWindow(para->bitname);
	cvDestroyWindow(para->convexname);

}



//
//	main関数
//
int main(int argc, char** argv) {
	int i = 1, key;
	int cls = 7;
	CvCapture *capture = NULL;
	char buf[100];
	FILE *file;
	errno_t err;


	//フルオートの場合以下の行必要
	for(i = 1; i < 10; i++) {
	
	printf("-----------------------------------------------\n");
	printf("version 01/10  i = %d \n", i);

	err = fopen_s(&file, "C:/Users/TK/flower1227/data/r-atsumori-0126.txt", "a");
	if(err != 0) printf("結果格納テキストファイルが見つかりません\n");

	//	画像を読み込む
	sprintf_s(buf, 100, "%s%d%s","C:/Users/TK/flower1227/r-atsumori/r-atsumori (", i, ").jpg");
	//sprintf_s(buf, 100, "%s","C:/Users/TK/flower1227/pencil.png");
	img = cvLoadImage( buf, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR );

	if ( img == NULL){
		//	画像が見つからなかった場合
		printf( "画像が見つかりません\n" );
		return -1;
	}

	const int w = img->width;
	const int h = img->height;	

	//平均値格納構造体
	Ratio ratiobox = {-1.0, 0.0};

	//白色の平均カラー用構造体
	AveColor aveWhite = {0, 0, 0};
	//白色の画像構造体 {bit, convex, label, bit-window, convex-window, color}
	Para paraWhite = {
		cvCreateImage( cvSize(w, h),IPL_DEPTH_8U,1), 
		cvCreateImage( cvSize(w, h),IPL_DEPTH_8U,3 ),
		cvCreateImage( cvSize(w, h),IPL_DEPTH_16S,1),
		"Result-bit","Result-convex", CV_RGB(0, 255, 255)
	};

	//RGB
	imgH = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	imgS = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	imgV = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);

	//マスク画像用(期間限定設置)
	IplImage* mask_img;
    mask_img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);

	//なぜかmainに無いとエラーが起きるので置いとく
	imgResult = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);

	//HSVへ変換
	IplImage *imgHSV = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
	cvCvtColor(img, imgHSV, CV_BGR2HSV); //※RGBではなくBGRということに注意

	// HSVを分解
	cvSplit(imgHSV, imgH, imgS, imgV, NULL);

	//ウィンドウを準備して画像を表示　※フルオート化の場合不要
	cvNamedWindow("Capture", CV_WINDOW_AUTOSIZE);
	cvShowImage("Capture", img);

	//フルオートにする場合、以下の2行カット
	//key = cvWaitKey(0);
	//if(key == 'b'){
		cv_ColorExtraction(img, paraWhite.bit, mask_img, CV_BGR2HSV, 165 ,15 ,76, 255, 127, 255); //H S Vの順
		//cv_ColorExtraction(img, paraWhite.bit, mask_img, CV_BGR2HSV, 15, 40, 80, 255, 0, 255); //H S Vの順

		Func(&paraWhite, &aveWhite, &ratiobox);
	//}

	//ウィンドウを準備して画像を表示　※フルオート化の場合不要
	cvShowImage("Result-mask", mask_img);
	//構造体内を表示
	Show(&paraWhite);

	//データファイルへの書き込み
	if(ratiobox.c_ratio == -1.0){
		printf("error!\n");
	} else {
		fprintf(file,"%d %f %d %d %d\n", cls, ratiobox.c_ratio, aveWhite.h * 2, (int)aveWhite.s * 100 / 255, (int)aveWhite.v * 100 / 255);
		printf("c_ratio = %lf\n", ratiobox.c_ratio);
		printf("d_ratio = %1.3f\n", ratiobox.d_ratio);
	}	
	//ratioを再び-1に初期化
	ratiobox.c_ratio = -1.0;

	//フルオートで行う場合、以下の行もカット　画像変換後の確認のため使用していた
	cvWaitKey(0);	

	//解放軍
	cvReleaseImage(&img);
	cvDestroyWindow("Capture");
	cvReleaseImage(&mask_img);
	cvDestroyWindow("Result-mask");

	cvReleaseImage(&imgHSV);
	cvReleaseImage(&imgH);
	cvReleaseImage(&imgS);
	cvReleaseImage(&imgV);

	//構造体内も片付ける
	Close(&paraWhite);

	fclose(file);
	//for文のケツ
	}

	return 0;
}