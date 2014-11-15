//���F�����o����v���O����
/* �\���̑���̏ڍׂ́@http://homepage3.nifty.com/mmgames/c_guide/16-02.html�@*/

#include <stdio.h>
#include <stdlib.h>
#include <highgui.h>
#include <cv.h>
#include <cxcore.h>
#include "Labeling.h"
#include "LabelingW.h"
/* Labeling.h�ɂ��āc http://oshiro.bpe.es.osaka-u.ac.jp/people/staff/imura/products/labeling */

#ifdef _DEBUG
//Debug���[�h�̏ꍇ
#pragma comment(lib,"cv200d.lib")
#pragma comment(lib,"cxcore200d.lib")
#pragma comment(lib,"cvaux200d.lib")
#pragma comment(lib,"highgui200d.lib")
#else
//Release���[�h�̏ꍇ
#pragma comment(lib,"cv200.lib")
#pragma comment(lib,"cxcore200.lib")
#pragma comment(lib,"cvaux200.lib")
#pragma comment(lib,"highgui200.lib")
#endif

#define ITERATIONS 2

//�L���v�`���摜�pIplImage
IplImage *img = NULL;
IplImage *imgResult, *imgH, *imgS, *imgV;

//�F�̕��ϒl���i�[����\����
typedef struct{
	int h;
	int s;
	int v;
}AveColor;


typedef struct{
	double c_ratio;
	double d_ratio;
}Ratio;

//�摜�̖{�̂�����i�[����\����
typedef struct{
	IplImage *bit;
	IplImage *convex;
	IplImage *label;
	char bitname[20];
	char convexname[20];
	CvScalar color;
}Para;

//
//  �ԗ̈��`�悵�ARGB�l���߂��������ԗ̈�ɉ�����
//
//	����:
//		skinImage       : �Ԃ̐F���o�摜�pIplImage
//		label           : ���x�����O��������
//		convexHullImage : ConvexHull�摜�pIplImage
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


				//���ꂼ��ϐ����쐬
				//hsv�{��
				h = char_h; s = char_s; v = char_v;
				//�ehsv�̑O��x�̐F
				ex_h_x = ex_char_h_x; ex_s_x = ex_char_s_x; ex_v_x = ex_char_v_x;
				//�ehsv�̑O��y�̐F
				ex_h_y = ex_char_h_y; ex_s_y = ex_char_s_y; ex_v_y = ex_char_v_y;
				//�O��x�̐F�Ƃ̍��ف@difference
				d_h_x = h - ex_h_x; d_s_x = s - ex_s_x; d_v_x = v - ex_v_x;
				//�O��y�̐F�Ƃ̍���
				d_h_y = h - ex_h_y; d_s_y = s - ex_s_y; d_v_y = v - ex_v_y;


				if( cvGetReal2D( label, y, x ) == 1) {
					//	�ő�̈悾�����ꍇ �ʐς��C���N�������g���Aconvex�𔒐F�ɃZ�b�g
					cvSet2D( convexHullImage, y, x, CV_RGB( 255, 255, 255 ) );
					cvSetReal2D(label, y, x, 1);			//���x����t������(1����ԑ傫���摜)
					cvSetReal2D(bitImage, y, x, 255);
				} else {
					//�łȂ������ꍇ�A���ɃZ�b�g
					cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 0 ) );
					cvSetReal2D( bitImage, y, x, 0 );
					//�����ߐڂ�RGB���߂��l�ł���΁A�␳
					//�f�t�H���g�̒l�@abs(d_h_x) <= 2  && abs(d_s_x) <= 13 && abs(d_v_x) <= 26
					if(cvGetReal2D( label, y, x-1 ) == 1 && abs(d_h_x) <= 2 && abs(d_s_x) <= 13 && abs(d_v_x) <= 26) {
						cvSet2D( convexHullImage, y, x, CV_RGB( 255, 0, 0));		//convex��Ԃŕ`��
						cvSetReal2D(label, y, x, 1);			//���x����t������(1����ԑ傫���摜)
						cvSetReal2D(bitImage, y, x, 255);		//bitImage��255��t�����Ĕ���`��

					//y��������U���x�݁@����������Ƃ��͏����ߎ��l��������H
					//} else if(cvGetReal2D( label, y-1, x ) == 1 && abs(d_h_y) < 2 && abs(d_s_y) < 6.75 && abs(d_v_y) < 12.5) {
					//	cvSet2D( convexHullImage, y, x, CV_RGB( 0, 255, 0 ));		//convex��΂ŕ`��
					//	cvSetReal2D(label, y, x, 1);			//���x����t������(1����ԑ傫���摜)
					//	cvSetReal2D(bitImage, y, x, 255);		//bitImage��255��t�����Ĕ���`��

					//} else if(cvGetReal2D( label, y, x-1 ) == 1 && s <= 51 && abs(d_s_x) <= 13 && abs(d_v_x) <= 26 ) {
					//	cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 255 ));		//convex��ŕ`��
					//	cvSetReal2D(label, y, x, 1);			//���x����t������(1����ԑ傫���摜)
					//	cvSetReal2D(bitImage, y, x, 255);		//bitImage��255��t�����Ĕ���`��
					}

				}
			}
		}
	}

}

//
//  �F�̋ߎ��ɂ��ԗ̈���J�E���g����ۂɌ�J�E���g���Ă��镔�����C������
//
//	����:
//		skinImage       : �Ԃ̐F���o�摜�pIplImage
//		label           : ���x�����O��������
//		convexHullImage : ConvexHull�摜�pIplImage
//
void reviseMaxArea(IplImage *bitImage, IplImage *label, IplImage *convexHullImage ) {

	CvScalar color;

	int count = 0, sequence = 0,same = 0;
	int mem = -1;
	int ybuf = bitImage->height;

	//for���J�n 
	for( int y = 0; y < bitImage->height; y++ ) {
		for(int x = 0; x < bitImage->width; x++ ) {
			//����x�������ɐԐF�̐�(x�̕␳)���������A�����Ă�����J�E���g
			color = cvGet2D(convexHullImage, y, x);
			if(color.val[0] == 0.0 && color.val[1] == 0.0 && color.val[2] == 255.0) { //BGR
				count++;
				sequence++;
			} else {
				//�V�[�P���X�͍ēx���Z�b�g
				sequence = 0;
			}

			if(count > bitImage->width/5) {
				//�����Ԃ̍��v�����摜��1/5��˔j���Ă�����ϐ�same���C���N�������g�����̍s��
				if(ybuf > y) ybuf = y;
				same++;
				y++;
			} else if(sequence > bitImage->width/8) {
				//�����Ԃ̒������摜��1/8�A�����Ă����炻�̏�ŏC������
				for(int m = 0; m < bitImage->width; m++){
					color = cvGet2D(convexHullImage, y, m);
					//convex�̐F���ԂȂ��
					if(color.val[0] != 255.0 || color.val[1] != 255.0 || color.val[2] != 255.0) {
						//convex�̔������� (== �␳��������O�̂��Ƃ��Ƃ̉Ԃ̕���)�@�ȊO��S�ēh��Ԃ�
						cvSetReal2D( bitImage, y, m, 0 );
						cvSet2D( convexHullImage, y, m, CV_RGB( 0, 0, 0 ) );
						cvSetReal2D(label, y, m, 0);			//���x����t������(1����ԑ傫���摜)
					}
				} //for��m�I��
				y++;
			} //if���I��
		}//for��x�I��
		count = 0;
		sequence = 0;
	} //for��y�I��


	//����same��10�ȏ゠��΁A�傫�����邽�ߑ啝�C�����s��(���Z�b�g�ł������̂ł́H)
	if(same > 10) {
		printf("x�̏C�����A������(same=%d)�������߁Ay=%d����啝�C��\n", same, ybuf);
		//�G���[�������������ďI��
		for( int y = ybuf -1 ; y < bitImage->height ; y++ ) {
			for(int x = 0; x < bitImage->width ; x++ ) {
				color = cvGet2D(convexHullImage, y, x);
				//�ǂꂩ��ł�255����Ȃ����̂�����΁c (== 255 255 255 �܂蔒�łȂ����)
				if(color.val[0] != 255.0 || color.val[1] != 255.0 || color.val[2] != 255.0) {
					//convex�̔������� (== �␳��������O�̂��Ƃ��Ƃ̉Ԃ̕���)�@�ȊO��S�ēh��Ԃ�
					//printf("%lf %lf %lf\n", color.val[0], color.val[1], color.val[2]);
					cvSetReal2D( bitImage, y, x, 0 );
					cvSet2D( convexHullImage, y, x, CV_RGB( 0, 0, 0 ) );
					cvSetReal2D(label, y, x, 0);			//���x����t������(1����ԑ傫���摜)
				}
			}
		}
	}

	cvDestroyWindow("conv");
}


//
//	�Ԃ̍ő�̈�̒��o(label==1)���s���A�ʐς��J�E���g����
//
//	����:
//		skinImage       : �Ԃ̐F���o�摜�pIplImage
//		label           : ���x�����O��������
//
//	�߂�l:
//		�ԗ̈�̖ʐ�
//
int pickupMaxArea(IplImage *bitImage, IplImage *label) {

	int flowerarea = 0;	//	�ԗ̈�̖ʐ�

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
//	ConvexHull�𐶐�����
//
//	����:
//		skinImage   : ���F���o�摜�pIplImage
//		flowerarea  : �ԗ̈�̖ʐ�(�_�̐�)
//		flowerpoint : �ԗ̈���̓_�̍��W�z��ւ̃|�C���^
//		hull        : ConvexHull�̒��_��handpoint�ɂ�����index�ԍ��ւ̃|�C���^
//		pointMatrix : �ԗ̈�p�s��ւ̃|�C���^
//		hullMatrix  : ConvexHull�p�s��ւ̃|�C���^
//
void createConvexHull(IplImage *bitImage, IplImage *label,int flowerarea, CvPoint **flowerpoint, int **hull,
					  CvMat *pointMatrix, CvMat *hullMatrix ) {
						  int i=0;

						  //	ConvexHull���v�Z���邽�߂ɕK�v�ȍs��𐶐�����
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

						  //	ConvexHull�𐶐�����
						  cvConvexHull2( pointMatrix, hullMatrix, CV_CLOCKWISE, 0 );

}


//
//	ConvexHull��`�悷��
//
//	����:
//		convexHullImage : ConvexHull�摜�pIplImage
//		flowerpoint     : �ԗ̈���̓_�̍��W�z��
//		hull            : ConvexHull�̒��_��handpoint�ɂ�����index�ԍ�
//		hullcount       : ConvexHull�̒��_�̐�
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
//	ConvexHull���̖ʐς����߂�
//
//	����:
//		convexHullImage : ConvexHull�摜�pIplImage
//		flowerpoint       : �ԗ̈���̓_�̍��W�z��
//		hull            : ConvexHull�̒��_��handpoint�ɂ�����index�ԍ�
//		hullcount       : ConvexHull�̒��_�̐��@�@
//
//	�߂�l:
//		ConvexHull���̖ʐ�
//
int calcConvexHullArea( IplImage *convexHullImage, CvPoint *flowerpoint, int *hull, int hullcount ) {

	//	ConvexHull�̒��_����Ȃ�s��𐶐�
	CvPoint *hullpoint = ( CvPoint * )malloc( sizeof( CvPoint ) * hullcount );
	CvMat hMatrix = cvMat( 1, hullcount, CV_32SC2, hullpoint );
	for( int i = 0; i < hullcount; i++ ) {
		hullpoint[i]=flowerpoint[hull[i]];
	}

	//	ConvexHull���̓_�̐��𐔂���
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

//2�_�Ԃ̋��������߂�
int Distance(CvPoint pt1, CvPoint pt2){
	int dis;
	
	dis = int(sqrt(pow(pt2.x - pt1.x, 2.0) + pow(pt2.y - pt1.y, 2.0)));

	return dis;
}

//
//�摜�̏d�S�����߂�֐�
//
//	����:
//		para		    : �摜�̍\����
//		flowerpoint     : �ԗ̈���̓_�̍��W�z��
//		hull            : ConvexHull�̒��_��handpoint�ɂ�����index�ԍ�
//		hullcount       : ConvexHull�̒��_�̐��@�@
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
			//1�{�ڂ̂ݐԁc����Ƃ���
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

	//�e�����̔䗦�̕��ϒl���o��(�����^�p��)
	ratiobox->d_ratio = sum / (double)(hullcount-1);
}

//
//	���x�����O�����̐F�̕��ϒl�����߂�
//
//	����:
//		skinImage       : �Ԃ̐F���o�摜�pIplImage
//		para			: �摜�̍\����
//		aveColor		: �Ԃ̐F�̕��ϒl���i�[����\����
//
void averageColor(Para *para, AveColor *aveColor) {

	unsigned char char_h, char_s, char_v;
	int h, s, v;
	int h_sum = 0, s_sum = 0, v_sum = 0, in = 0;

	for(int x = 0; x < para->bit->width; x++ ) {
		for( int y=0; y < para->bit->height; y++ ) {
			if( cvGetReal2D( para->label, y, x ) == 1) {
				//	�ő�̈悾�����ꍇ �c convex�Ŕ��������ꍇ
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
	aveColor->h = h_sum/in;	//���̐����͐��K��
	aveColor->s = s_sum/in;
	aveColor->v = v_sum/in;
	printf("Average H: %d(/180) S: %d(/100) V: %d(/100) \n", aveColor->h * 2, (int)aveColor->s * 100 / 255, (int)aveColor->v * 100 / 255);
}

//
//	�����̈���Ԃ���
//
//	����:
//		skinImage : ���F���o�摜�pIplImage
//		temp      : �ꎞ�ۑ��pIplImage
//
void interpolate( IplImage *skinImage, IplImage *temp, int num) {
	//�c����num��s��
	cvDilate( skinImage, temp, NULL, num );

	//���k��num��s��
	cvErode( temp, skinImage, NULL, num );
}


//--------------------------------------------------------------- 
//�y�֐����@�z�Fcv_ColorExtraction 
//�y�����T�v�z�F�F���o 
//�y�����@�@�z�Fsrc_img        = ���͉摜(8bit3ch) 
//�@�@�@�@�@�@�Fdst_img        = �o��2�l�摜(8bit1ch) 
//�@�@�@�@�@�@�Fmask_img        = �o�̓}�X�N�摜(8bit3ch) 
//�@�@�@�@�@�@�Fcode        = �F��Ԃ̎w��iCV_BGR2HSV,CV_BGR2Lab�Ȃǁj
//�@�@�@�@�@�@�Fch1_lower    = ch1�̂������l(��)
//�@�@�@�@�@�@�Fch1_upper    = ch1�̂������l(��)
//�@�@�@�@�@�@�Fch2_lower    = ch2�̂������l(��)
//�@�@�@�@�@�@�Fch2_upper    = ch2�̂������l(��)
//�@�@�@�@�@�@�Fch3_lower    = ch3�̂������l(��)
//�@�@�@�@�@�@�Fch3_upper    = ch3�̂������l(��)
//�y�߂�l�@�z�F�Ȃ� 
//�y���l�@�@�z�Flower <= upper�̏ꍇ�Alower�ȏ�upper�ȉ��͈̔͂𒊏o�A
//�@�@�@�@�@�@�Flower >  upper�̏ꍇ�Aupper�ȉ�lower�ȏ�͈̔͂𒊏o���܂��B
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

    //code�Ɋ�Â����J���[�ϊ�
    Color_img = cvCreateImage(cvGetSize(src_img), src_img->depth, src_img->nChannels);
    cvCvtColor(src_img, Color_img, code);
       
    //3Ch��LUT�쐬
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
        //LUT�̐ݒ�
        cvSet1D(lut, i, cvScalar(val[0], val[1], val[2]));
    }

    //3Ch���Ƃ�LUT�ϊ��i�e�`�����l�����ƂɂQ�l�������j
    cvLUT(Color_img, Color_img, lut);
    cvReleaseMat(&lut);

    //�e�`�����l�����Ƃ�IplImage���m�ۂ���
    ch1_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);
    ch2_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);
    ch3_img = cvCreateImage(cvGetSize(Color_img), Color_img->depth, 1);

    //�`�����l�����Ƃɓ�l�����ꂽ�摜�����ꂼ��̃`�����l���ɕ�������
    cvSplit(Color_img, ch1_img, ch2_img, ch3_img, NULL);

    //3Ch�S�Ă�AND�����A�}�X�N�摜���쐬����B

    cvAnd(ch1_img, ch2_img, dst_img);
    cvAnd(dst_img, ch3_img, dst_img);

    //���͉摜(src_img)�̃}�X�N�̈���o�͉摜(dst_img)�փR�s�[����
    cvZero(mask_img);
    cvCopy(src_img, mask_img, dst_img);

    //���
    cvReleaseImage(&Color_img);
    cvReleaseImage(&ch1_img);
    cvReleaseImage(&ch2_img);
    cvReleaseImage(&ch3_img);

}

//
//	��L�̊֐���p���Čv�Z���s���AconvexHull���̉Ԃ̊������o���A�`����s���֐�
//
//	����:
//		para :	�摜�̍\����
//
//	�߂�l:
//		ratio :	convex���ɂ�����Ԃ̊���
//

void Func(Para *para, AveColor *aveColor, Ratio *ratiobox) {

	double convexratio = -1.0;
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* find_contour = NULL;

	//���k�A�c�������ĉ摜�̃m�C�Y������	
	//interpolate(para->bit, para->bit, 2);

	//	���x�����O���s�� 
	Label *labeling = createLabeling();
	exec( labeling, para->bit, para->label, true, 2000);
	/* exec(���x�����O�ϐ�, ���͉摜, ���x�����O���ʂ��i�[����ϐ�, 
	�̈悪�傫���̍~���Ƀ\�[�g����邩, �ŏ��̗̈�T�C�Y) */

	if(getNumOfResultRegions(labeling) > 0 ) {
		//�ŏ��̗̈�T�C�Y�����傫�ȗ̈悪�������ꍇ
		int flowerarea;		//	�ԗ̈�̖ʐ�
		int hullarea;		//	ConvexHull���̖ʐ�
		int hullcount;		//	ConvexHull�̒��_�̐�
		CvPoint *flowerpoint;	//	�ԗ̈���̓_�̍��W�z��
		int *hull;			//	ConvexHull�̒��_��flowerpoint�ɂ�����index�ԍ�
		CvMat pointMatrix;	//	�ԗ̈�p�s��
		CvMat hullMatrix;	//	ConvexHull�p�s��


		//�֊s�����o����
		int contour_num_color = cvFindContours(cvCloneImage(para->bit),
			storage, 
			&find_contour,
			sizeof(CvContour), 
			CV_RETR_EXTERNAL, 
			CV_CHAIN_APPROX_NONE, 
			cvPoint(0,0)
			);

		//printf("check-2 \n");

		//�ԗ̈��`��E�ߎ��l��␳
		drawMaxArea(para->bit, para->label, para->convex);
		//�G���[�ӏ��̏C��
		reviseMaxArea(para->bit, para->label, para->convex);

		//bit-2�l�摜���k���E���k
		interpolate(para->bit, para->bit, 2);
		//�֊s�̕\��
		cvDrawContours(para->convex, find_contour, para->color, para->color, 1, 2, 8, cvPoint(0,0)); //���F�̏ꍇ(2�Ԗڂ�yellow�g�p���Ȃ�)]

		//printf("check-3 \n");

		//�ő�̈�𒊏o���A�ʐς��J�E���g
		flowerarea = pickupMaxArea(para->bit, para->label);

		//printf("check-4 flowerarea = %d\n", flowerarea);

		//ConvexHull�𐶐�
		createConvexHull(para->bit,para->label, flowerarea, &flowerpoint, &hull, &pointMatrix, &hullMatrix);
		hullcount = hullMatrix.cols;

		//printf("check-5 \n");

		//ConvexHull��`��
		drawConvexHull(para->convex, flowerpoint, hull, hullcount);

		//printf("check-6 \n");

		//ConvexHull���̖ʐς����߂�
		hullarea = calcConvexHullArea(para->convex, flowerpoint, hull, hullcount);

		//printf("check-7 \n");

		//�ő�̈�̐F�̕��ϒl�����߂� �\���̃|�C���^�n��
		averageColor(para, aveColor);

		//printf("check-8 \n");

		//�d�S�����߁A�e���_�܂Ő�������
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
//	�e�E�B���h�E���쐬���A�摜�\���̏������s���Ă����֐�
//
//	����:
//		img	 :	���ƂȂ�摜 �L���v�`���E�B���h�E�p�ɏ���
//		para :	�摜�̍\����
//

void Show(Para *para){

	//�����̍\���̂̃E�B���h�E���쐬��
	cvNamedWindow(para->bitname, CV_WINDOW_AUTOSIZE);
	cvNamedWindow(para->convexname, CV_WINDOW_AUTOSIZE);

	//�\��
	cvShowImage(para->bitname, para->bit);
	cvShowImage(para->convexname, para->convex);

}

//
//	�e�E�B���h�E����A�������̉�������s���Ă����֐�
//
//	����:
//		img	 :	���ƂȂ�摜 �L���v�`���E�B���h�E�p�ɏ���
//		para :	�摜�̍\����
//

void Close(Para *para){

	//�������̉��
	cvReleaseImage(&para->bit);
	cvReleaseImage(&para->convex);
    cvReleaseImage(&para->label);

	//�E�B���h�E��j��
	cvDestroyWindow(para->bitname);
	cvDestroyWindow(para->convexname);

}



//
//	main�֐�
//
int main(int argc, char** argv) {
	int i = 1, key;
	int cls = 8;
	CvCapture *capture = NULL;
	char buf[100];
	FILE *file;
	errno_t err;


	//�t���I�[�g�̏ꍇ�ȉ��̍s�K�v
	for(i = 1; i < 22; i++) {
	
	printf("-----------------------------------------------\n");
	printf("version 01/10  i = %d \n", i);

	err = fopen_s(&file, "C:/Users/TK/flower1227/data/odamaki-0126.txt", "a");
	if(err != 0) printf("���ʊi�[�e�L�X�g�t�@�C����������܂���\n");

	//	�摜��ǂݍ���
	sprintf_s(buf, 100, "%s%d%s","C:/Users/TK/flower1227/odamaki/odamaki (", i, ").jpg");
	//sprintf_s(buf, 100, "%s","C:/Users/TK/flower1227/pencil.png");
	img = cvLoadImage( buf, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR );

	if ( img == NULL){
		//	�摜��������Ȃ������ꍇ
		printf( "�摜��������܂���\n" );
		return -1;
	}

	const int w = img->width;
	const int h = img->height;	

	//���ϒl�i�[�\����
	Ratio ratiobox = {-1.0, 0.0};

	//���F�̕��σJ���[�p�\����
	AveColor aveWhite = {0, 0, 0};
	//���F�̉摜�\���� {bit, convex, label, bit-window, convex-window, color}
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

	//�}�X�N�摜�p(���Ԍ���ݒu)
	IplImage* mask_img;
    mask_img = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);

	//�Ȃ���main�ɖ����ƃG���[���N����̂Œu���Ƃ�
	imgResult = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);

	//HSV�֕ϊ�
	IplImage *imgHSV = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);
	cvCvtColor(img, imgHSV, CV_BGR2HSV); //��RGB�ł͂Ȃ�BGR�Ƃ������Ƃɒ���

	// HSV�𕪉�
	cvSplit(imgHSV, imgH, imgS, imgV, NULL);

	//�E�B���h�E���������ĉ摜��\���@���t���I�[�g���̏ꍇ�s�v
	cvNamedWindow("Capture", CV_WINDOW_AUTOSIZE);
	cvShowImage("Capture", img);

	//�t���I�[�g�ɂ���ꍇ�A�ȉ���2�s�J�b�g
	key = cvWaitKey(0);
	if(key == 'b'){
		cv_ColorExtraction(img, paraWhite.bit, mask_img, CV_BGR2HSV, 105, 135, 76, 255, 127, 255); //H S V�̏�
		//cv_ColorExtraction(img, paraWhite.bit, mask_img, CV_BGR2HSV, 15, 40, 80, 255, 0, 255); //H S V�̏�

		Func(&paraWhite, &aveWhite, &ratiobox);
	}

	//�E�B���h�E���������ĉ摜��\���@���t���I�[�g���̏ꍇ�s�v
	cvShowImage("Result-mask", mask_img);
	//�\���̓���\��
	Show(&paraWhite);

	//�f�[�^�t�@�C���ւ̏�������
	if(ratiobox.c_ratio == -1.0){
		printf("error!\n");
	} else {
		fprintf(file,"%d %f %d %d %d\n", cls, ratiobox.c_ratio, aveWhite.h * 2, (int)aveWhite.s * 100 / 255, (int)aveWhite.v * 100 / 255);
		printf("c_ratio = %lf\n", ratiobox.c_ratio);
		printf("d_ratio = %1.3f\n", ratiobox.d_ratio);
	}	
	//ratio���Ă�-1�ɏ�����
	ratiobox.c_ratio = -1.0;

	//�t���I�[�g�ōs���ꍇ�A�ȉ��̍s���J�b�g�@�摜�ϊ���̊m�F�̂��ߎg�p���Ă���
	cvWaitKey(0);	

	//����R
	cvReleaseImage(&img);
	cvDestroyWindow("Capture");
	cvReleaseImage(&mask_img);
	cvDestroyWindow("Result-mask");

	cvReleaseImage(&imgHSV);
	cvReleaseImage(&imgH);
	cvReleaseImage(&imgS);
	cvReleaseImage(&imgV);

	//�\���̓����Еt����
	Close(&paraWhite);

	fclose(file);
	//for���̃P�c
	}

	return 0;
}