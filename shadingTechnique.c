//
// sample-zbuf.c
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "imgio.h"
#include "imgioX11.h"
#include "VecMat.h"

#define RESX 512/* 画像サイズ (横) */
#define RESY 512/* 画像サイズ (縦) */

#define FL   256/* 透視投影変換の焦点距離 */

#define COU 20





double zbuf[RESX][RESY];// Z バッファ用配列

//原点を中央に持ってくるような変換
int ConvIJtoXY(int i, int j, double *x, double *y)
{
    *x = (double)(i - RESX/2);
    *y = (double)(RESY/2 - j);
    return 0;
}

// 画素位置と、スクリーン座標との変換
int ConvXYtoIJ(double x, double y, int *i, int *j)
{
    *i = (int)(x + RESX/2);
    *j = (int)(RESY/2 - y);
    return 0;
}

// Z バッファの初期化。右手系を考え、奥行き= Z が負
int InitZbuf()
{
    int i, j;

    for (j=0; j<RESY; j++) {
        for (i=0; i<RESX; i++) {
            zbuf[i][j] = -99999.99;// 最奥値で初期化
        }
    }
    return 0;
}


int InitImage(IMG *img_out, int R, int G, int B)
{
    int i, j;
    for (j=0; j<img_out->height; j++) {
        for (i=0; i<img_out->width; i++) {
            R(*img_out, i, j) = R;
            G(*img_out, i, j) = G;
            B(*img_out, i, j) = B;
        }
    }
    return 0;
}


int DrawPoint(IMG *img_out, int X, int Y, int r, int g, int b)
{
    int x, y;

    x = X + RESX/2; 
    y = RESY/2 -Y;  

    if (x>=0 && x<W(*img_out) && y>=0 && y<H(*img_out)) {
        R(*img_out, x, y) = r;
        G(*img_out, x, y) = g;
        B(*img_out, x, y) = b;
    }
    return 0;
}




// 三角形の表示。三頂点 p0, p1, p2 と変換行列 mat を与える。
int DispOneTriangle(IMG *img, Mat mat, Vec p0, Vec p1, Vec p2, int R, int G, int B, double th)
{
    int i, j;
    Vec q0, q1, q2, q10, q20, qq;
    Vec v0, v1, v2;
    Vec z01, z12, z20;
    Vec n;
    double t;
    Vec amb = {0, 30, 0, 1};//環境光成分
    Vec lv;//平行光線の方向ベクトル
    double a;

    p0[3] = p1[3] = p2[3] = 1;// 同次座標表現のための処理(一応)

    lv[0] = 0.5;
    lv[1] = 0.5;
    lv[2] = -0.5;
    lv[3] = 1;

    //4x4 行列 mat x 4D ベクトルp0を計算しq0に格納
    MultipleVec(mat, p0, q0);
    MultipleVec(mat, p1, q1);
    MultipleVec(mat, p2, q2);
    
    // ベクトルの差 (q10 = q1 - q0)
    SubVec(q1, q0, q10);
    SubVec(q2, q0, q20);
    // ベクトルの外積 q10 × vq20 計算(同次座標表現の第 4 成分は無視)
    OuterProduct(q10, q20, n);
    
    // 画面に対応する領域全体を走査して、三角形の内部か判定
    for (j=0; j< RESY; j++) {
        for (i=0; i<RESX; i++) {

            // レイの生成
            ConvIJtoXY(i, j, &qq[0], &qq[1]);
            qq[0] += 0.5; qq[1] += 0.5; qq[2] = -FL; qq[3] = 1;

            // レイと三角形の乗る平面との交点計算
            // ベクトルの内積 x,y の計算(同次座標表現の第 4 成分は無視)
            t = InnerProduct(q0, n)/InnerProduct(n, qq);
            qq[0] *= t; qq[1] *= t; qq[2] *= t;

            // 交点から3頂点へのベクトルの外積の z 成分を利用して内外判定
            SubVec(q0, qq, v0);
            SubVec(q1, qq, v1);
            SubVec(q2, qq, v2);
            
            // ベクトルの外積 x*y 計算(同次座標表現の第 4 成分は無視)
            OuterProduct(v0, v1, z01);
            OuterProduct(v1, v2, z12);
            OuterProduct(v2, v0, z20);

            // 三角形内部なら
            if ( (z01[2] >0 && z12[2] > 0 && z20[2] > 0) || (z01[2] <0 && z12[2] < 0 && z20[2] < 0)    ){
                //zbuf との比較
                if (zbuf[i][j] < qq[2]) {
                    qq[0] = -FL*qq[0]/qq[2];
                    qq[1] = -FL*qq[1]/qq[2];
                    a = InnerProduct(qq,lv);
                    a /= sqrt(InnerProduct(qq,qq));
                    a /= sqrt(InnerProduct(lv,lv));
                    if (a > 0.0) {
                        R = a * R + amb[0];
                        G = a * G + amb[1];
                        B = a * B + amb[2];
                    } else {
                        R = amb[0];
                        G = amb[1];
                        B = amb[2];
                    }
                    if(R >= 255){
                        R = 255;
                    }
                    if(G >= 255){
                        G = 255;
                    }
                    if(B >= 255){
                        B = 255;
                    }
                    // 投影点を描画
                    DrawPoint(img, qq[0], qq[1], R,G,B); 
                    // zbuf を更新
                    zbuf[i][j] = qq[2];
                }
            }
        }
    }
    return 0;
}

//ラジアン変換
double degrees_to_radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

double calculate_sin(double angle_degrees) {
    double angle_radians = degrees_to_radians(angle_degrees);
    return sin(angle_radians);
}

double calculate_cos(double angle_degrees) {
    double angle_radians = degrees_to_radians(angle_degrees);
    return cos(angle_radians);
}

/* === main === */
int main(int argc, char *argv[])
{
    IMG img_out;
    FILE *fp;
    Mat mat;// 変換
    double th;
    Vec p[20][3];
    int c[20][3];
    int i;
    Vec p0, p1, p2, p3;
    double POA[] = {0,1,1.618};

    double POB[] = {1.618,0,1};
    double POC[] = {0,-1,1.618};
    double POD[] = {-1.618,0,1};
    double POE[] = {-1,1.618,0};
    double POF[] = {1,1.618,0};

    double POG[] = {1.618,0,-1};
    double POH[] = {1,-1.618,0};
    double POI[] = {-1,-1.618,0};
    double POJ[] = {-1.618,0,-1};
    double POK[] = {0,1,-1.618};

    double POL[] = {0,-1,-1.618};

    fprintf(stderr, "START\n");
    if (argc != 2) {
        fprintf(stderr, "USAGE: %s output.ppm\n", argv[0]);
        exit(0);
    }

    imgioX11_InitWindow();// 途中経過表示用の Window を生成
    cIMG(RESX, RESY, &img_out, COLOR);// IMG 型のデータ生成

    //三角面表現用の点xyz
    void create(int x, double xyz0[], double xyz1[], double xyz2[], Vec p[20][3]){
        p[x][0][0] = xyz0[0];
        p[x][0][1] = xyz0[1];
        p[x][0][2] = xyz0[2];
        
        p[x][1][0] = xyz1[0];
        p[x][1][1] = xyz1[1];
        p[x][1][2] = xyz1[2];

        p[x][2][0] = xyz2[0];
        p[x][2][1] = xyz2[1];
        p[x][2][2] = xyz2[2];
    }

    create(0, POA, POB, POC, p);
    create(1, POA, POC, POD, p);
    create(2, POA, POD, POE, p);
    create(3, POA, POE, POF, p);
    create(4, POA, POF, POB, p);
    
    create(5, POH, POB, POC, p);
    create(6, POI, POC, POD, p);
    create(7, POJ, POD, POE, p);
    create(8, POK, POE, POF, p);
    create(9, POG, POF, POB, p);

    create(10, POL, POG, POH, p);
    create(11, POL, POH, POI, p);
    create(12, POL, POI, POJ, p);
    create(13, POL, POJ, POK, p);
    create(14, POL, POK, POG, p);

    create(15, POB, POG, POH, p);
    create(16, POC, POH, POI, p);
    create(17, POD, POI, POJ, p);
    create(18, POE, POJ, POK, p);
    create(19, POF, POK, POG, p);
    



    for (i=0; i<COU; i++) {
        c[i][0] = 0;//100*drand48()+155;//155~255の間の小数,RGB値
        c[i][1] = 255;//100*drand48()+155;
        c[i][2] = 0;100*drand48()+155;
    }

    for (th=0; th<2*M_PI; th=th+0.01) {

        InitImage(&img_out,255, 255, 255);// RGB値で画像内のクリア
        InitZbuf();

        InitMat(mat);// 変換行列の初期化 (単位行列)
        ScaleMat(150, 150, 150, mat);// 拡大
        TranslateMat(0, 0, -512, mat);// 立方体を視野の前方へ
        
        for (i=0; i<COU; i++) {
            DispOneTriangle(&img_out, mat, p[i][0], p[i][1], p[i][2], c[i][0], c[i][1], c[i][2], th);
        }
        
        imgioX11_DisplayImage(img_out);/* 途中経過表示の更新 */

    }
    
    wIMG(img_out, argv[1]);/* 画像の書き出し */
    fprintf(stderr, "SAVED, and WAIT PUSH Q,q \n");
    imgioX11_RedrawLoop(img_out);/* イベント再描画設定。Q, q で終了 */

    
    return 0;
}
