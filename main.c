#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 256
#define BLOCK_SIZE 16

int main(){
	double f_Real[N][N], f_after_R[N][N], f_after_I[N][N], F_Real[N][N], F_Im[N][N];
	int Power_spectrum[N][N];
	unsigned char read_file[N][N];
	FILE *fp;
	clock_t start, end;
	int i, j, k, l, u=N/2, v=N/2;
	double x, y, max[3];
	double ratio;
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
			f_after_R[i][j] = 0;
			f_after_I[i][j] = 0;
		}
	}
	
	/*
	  -----------------------------
	  　ここからフーリエ変換の処理
	  -----------------------------
	*/
	printf("\n----------------------\n");
	printf("      Start DFT       \n");
	printf("----------------------\n");	

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
	max[0] = 0;
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] = log10(pow(F_Real[i][j], 2)+pow(F_Im[i][j], 2));
			//最大値を求めている(あとで正規化に使う)
			if( max[0] < Power_spectrum[i][j] ){
				max[0] = Power_spectrum[i][j];
				max[1] = i;
				max[2] = j;
			}
		}
	}
	/* データの正規化 */
	ratio = 254 / Power_spectrum[(int)max[1]][(int)max[2]];
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] *= ratio;
			read_file[i][j] = (unsigned char)Power_spectrum[i][j];
		}
	}
	
		
    /*
	------------------
	   画像書き込み 
	------------------
	*/
	fp=fopen("image/lenna_F_Re.raw","wb+"); //フーリエ変換後の実部の画像
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&F_Real[i][j], 1, 1, fp);
		}
	}
	fclose(fp);
	fp=fopen("image/lenna_F_Im.raw","wb+"); //フーリエ変換後の虚部の画像
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&F_Im[i][j], 1, 1, fp);
		}
	}
	fclose(fp);
	fp=fopen("image/lenna_Fourier.raw","wb+"); //フーリエ変換後のパワースペクトル
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&read_file[i][j], 1, 1, fp);
		}
	}
	fclose(fp);
	printf("\n----------------------\n");
	printf("        Fin DFT       \n");
	printf("----------------------\n");	

	/*
	  ----------------------------
	     ここからフィルタリング
	  ----------------------------
	*/
	printf("\n----------------------\n");
	printf("    Start Filtering   \n");
	printf("----------------------\n");	
	/* 100*100より外を0にする理想的フィルタリング */
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			if( i<N/2-100 || i>N/2+100 || j<N/2-100 || j>N/2+100 ){
				F_Real[i][j] = 0;
				F_Im[i][j] = 0;
			}
		}
	}
	
	/* 
	-------------------------
	  パワースペクトル生成 
	-------------------------
	*/
	max[0] = 0;
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] = log10(pow(F_Real[i][j], 2)+pow(F_Im[i][j], 2));
			//最大値を求めている(あとで正規化に使う)
			if( max[0] < Power_spectrum[i][j] ){
				max[0] = Power_spectrum[i][j];
				max[1] = i;
				max[2] = j;
			}
		}
	}
	/* データの正規化 */
	ratio = 254 / Power_spectrum[(int)max[1]][(int)max[2]];
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] *= ratio;
			read_file[i][j] = (unsigned char)Power_spectrum[i][j];

		}
	}
	
	fp=fopen("image/After_Filter_lenna_Fourier.raw","wb+"); //フーリエ変換後のパワースペクトル
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&read_file[i][j], 1, 1, fp);
		}
	}
	fclose(fp);

	printf("\n----------------------\n");
	printf("     Fin Filtering    \n");
	printf("----------------------\n");	

	
	/*
	  ----------------------------
	      ここからIDFTの処理
	  ----------------------------
	*/
	/* 実部と虚部をそれぞれ計算する */
	//iはuを回している
	printf("\n----------------------\n");
	printf("      Start IDFT      \n");
	printf("----------------------\n");	
	//iはuを回している
	for( i=0; i<N; i++ ){
		//jはvを回している
		for( j=0; j<N; j++ ){
			//Σの計算結果
			for( k=0; k<N; k++ ){
				for( l=0; l<N; l++ ){
					//u=0,1,...,N-1 v=0,1,...,N-1
					f_after_R[i][j] += F_Real[k][l] * cos(2.0 * M_PI * ((double)i*k + (double)j*l) / N) / (N*N) - F_Im[k][l] * sin(2.0 * M_PI * ((double)i*k + (double)j*l) / N) / (N*N);
				}
			}
		}
	}
	
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			Power_spectrum[i][j] = f_after_R[i][j];
		}
	}
	
	/* 画像出力 */
	fp=fopen("image/After_IDFT.raw","wb+"); //フーリエ変換後のパワースペクトル
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fwrite(&Power_spectrum[i][j], 1, 1, fp);
		}
	}
	
	printf("\n----------------------\n");
	printf("       Fin IDFT       \n");
	printf("----------------------\n");	


	end = clock();
	printf("time = %lus\n",(end-start)/CLOCKS_PER_SEC);
	return 0;
}
