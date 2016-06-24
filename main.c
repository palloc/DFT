#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 256
#define BLOCK_SIZE 16

int main(){
	double f_Real[N][N], F_Real[N][N], F_Im[N][N];
	int Power_spectrum[N][N];
	unsigned char read_file[N][N];
	FILE *fp;
	clock_t start, end;
	int i, j, k, l, u=N/2, v=N/2;
	double x, y, max, min;

	//実行時間計測開始
	start = clock();

	//画像をread_fileに格納する
	fp = fopen("image/lenna.256", "rb");
	fread(read_file, sizeof(unsigned char), N*N, fp);
	fclose(fp);

	//double型の配列に格納し直す&Fの初期化
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			f_Real[i][j] = read_file[i][j];
			F_Real[i][j] = 0;
			F_Im[i][j] = 0;
		}
	}
	
	/*
	  -----------------------------
	  　ここからフーリエ変換の処理
	  -----------------------------
	*/
	/* 実部と虚部をそれぞれ計算する */
	//iはuを回している
	for( i=-N/2; i<N/2; i++ ){
		u = N/2 + i;
		//jはvを回している
		for( j=-N/2; j<N/2; j++ ){
			v = N/2 + j;
			//Σの計算結果
			for( k=0; k<N; k++ ){
				for( l=0; l<N; l++ ){
					//u=0,1,...,N-1 v=0,1,...,N-1
					F_Real[u][v] += f_Real[k][l] * cos(2.0 * M_PI * ((double)i*k+(double)j*l) / N);
					F_Im[u][v] -= f_Real[k][l] * sin(2.0 * M_PI * ((double)i*k+(double)j*l) / N);
				}
			}
		}
	}

	/* 
	-------------------------
	  パワースペクトル生成 
	-------------------------
	*/
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] = log10(pow(F_Real[i][j], 2)+pow(F_Im[i][j], 2));
			read_file[i][j] = (unsigned char)Power_spectrum[i][j];
			printf("%d,",read_file[i][j]);
		}
	}

	
	//画像書き込み
	fp=fopen("lenna_F_Re.raw","wb+"); //フーリエ変換後の実部の画像
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&F_Real[i][j], 1, 1, fp);
		}
	}
	fclose(fp);
	fp=fopen("lenna_F_Im.raw","wb+"); //フーリエ変換後の虚部の画像
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&F_Im[i][j], 1, 1, fp);
		}
	}
	fclose(fp);
	fp=fopen("lenna_Fourier.raw","wb+"); //フーリエ変換後のパワースペクトル
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&read_file[i][j], 1, 1, fp);
		}
	}
	fclose(fp);

	end = clock();
	printf("time = %lu\n",(end-start)/CLOCKS_PER_SEC);
	return 0;
}
