/*
	このヘッダファイルは，
		・Generalization
		・Specificity
	を求めるプログラムが格納されています．

	このヘッダファイルを使用する際は
		・random.h（NumericalRecipesの正規乱数を使用するため）
		・JI.h（一致度算出）
	をincludeしてください．
*/
#include	"other/narimhd.h"

//投影・逆投影
template<class T, class T2>
void	Projection_Back_Projection( T &X, T &myu, T2 &U, T &X_out, int pe, int te){

	//pe:PCA前の軸数
	//te:PCA後の軸数

	//投影
	nari::vector<double> Z( te, 0.0);
#ifdef	_OPENMP
//std::cout << "Projection?" << std::endl;
#pragma omp parallel for schedule(static)
#endif
	for( int t = 0; t < te; t++){
		for( int p = 0; p < pe; p++){
			Z[t] += U( p, t) * ( X[p] - myu[p]);
		}
	}
	//逆投影
	nari::vector<double> UUtX( pe, 0.0);
#ifdef	_OPENMP
//std::cout << "Back Projection?" << std::endl;
#pragma omp parallel for schedule(static)
#endif
	for( int p = 0; p < pe; p++){
		for( int t = 0; t < te; t++){
			UUtX[p] += Z[t] * U( p, t);
		}
		X_out[p] = UUtX[p] + myu[p];
	}
}

template<class T, class T2>
void	Back_Projection( T &Z, T &myu, T2 &U, T &X_out, int pe, int te){

	//pe:PCA前の軸数
	//te:PCA後の軸数
	//逆投影
	nari::vector<double> UUtX( pe, 0.0);	
	for( int p = 0; p < pe; p++){
		for( int t = 0; t < te; t++){
			UUtX[p] += Z[t] * U( p, t);
		}
		X_out[p] = UUtX[p] + myu[p];
	}
}

template<class T, class T2>
void	Projection( T &X, T &myu, T2 &U, T &Z, int pe, int te){

	//pe:PCA前の軸数
	//te:PCA後の軸数

	//投影
	for( int t = 0; t < te; t++){
		for( int p = 0; p < pe; p++){
			Z[t] += U( p, t) * ( X[p] - myu[p]);
		}
	}
}

//ランダムな任意形状の生成
template<class T, class T2, class T3, class T4>
void	Arbitrarily_Shaped_Object( T &myu, T2 &U, T3 &evalue, T4 &X_out, int pe, int te, long *idum){

	//pe:PCA前の軸数
	//te:PCA後の軸数
	nari::vector<double>	Z( te, 0.0);

	for( int t = 0; t < te; t++){
		Z[t] = sqrt(evalue[t]) * normal( idum);
	}

	nari::vector<double> UUtX( pe, 0.0);	
#ifdef	_OPENMP
//std::cout << "Back Projection?" << std::endl;
#pragma omp parallel for schedule(static)
#endif
	for( int p = 0; p < pe; p++){
		for( int t = 0; t < te; t++){
			UUtX[p] += Z[t] * U( p, t);
		}
		X_out[p] = UUtX[p] + myu[p];
	}
}







//generalization
template<class T1, class T2, class T3>
void	generalization( T1 &X, T2 &myu, T3 &U, int pe, int te, int ncase, std::string outDir){
	int se = pe;
	int	cross = 0;
	int join = 0;
	double	ji = 0.0;
	double	generalization = 0.0;
	nari::vector<double>	Xi(se,0.0);
	nari::vector<double>	Xi_out(se,0.0);
	nari::vector<unsigned char>	Xi_label(se,0);
	nari::vector<unsigned char>	Xi_(se,0);
			
	if(!nari::system::directry_is_exist( outDir))
	{
		nari::system::make_directry( outDir);
		std::cout << "Directry was generated : " << outDir << std::endl;
	}

	FILE	*fp;
	std::ostringstream output;
	std::string	i_fn = outDir + "Generalization.csv";
	fp = fopen( i_fn.c_str(), "w");
	fprintf( fp, "症例番号,J.I.\n");

	for( int j = 0; j < ncase; j++){
		for( int p = 0; p < pe; p++)	Xi[p] = X( p, j);

		Projection_Back_Projection( Xi, myu, U, Xi_out, pe, te);
		
		//一致度の算出
		//cross = 0;
		//join = 0;
		//for( int i = 0; i < se; i++){
		//	if( Xi_out[i] < 0.0 || Xi[i] < 0.0){
		//		cross++;
		//		if( Xi_out[i] < 0.0 && Xi[i] < 0.0)	join++;
		//	}
		//}
		//ji = (double)join / (double)cross;

		for( int i = 0; i < se; i++){
			if( Xi_out[i] < 0.0)	Xi_[i] = 1;
			else	Xi_[i] = 0;
			if( Xi[i] < 0.0)	Xi_label[i] = 1;
			else	Xi_label[i] = 0;
		}
		ji = JI( Xi_.ptr(), Xi_label.ptr(), se);


		//逆投影後の形状出力
		output << outDir << "test" << j << ".raw";
		Xi_.save_file_bin( output.str());
		output.str("");

		fprintf( fp, "%d,%f\n", j+1, ji);
		generalization += ji;
		std::cout << "generalization" << j+1 << " = " << ji << std::endl;
	}
	generalization /= (double)ncase;
	fprintf( fp, "Generalization,%f\n", generalization);
	fclose( fp);
}


//specificity
template<class T1, class T2, class T3, class T4>
void	specificity( T1 &X, T2 &myu, T3 &U, T4 &evalue, int pe, int te, int ncase, int ngene, long idum, std::string outDir){
	int se = pe;
	double	ji = 0.0;
	double	case_max_ji = -1.0;
	double	max_ji = -1.0;
	double	min_ji = 2.0;
	double	specificity = 0.0;	
	int	cross = 0;
	int join = 0;
	nari::vector<double>	Xi(se,0.0);
	nari::vector<double>	Xi_out(se,0.0);
	nari::vector<unsigned char>	Xi_label(se,0);
	nari::vector<unsigned char>	Xi_(se,0);

	if(!nari::system::directry_is_exist( outDir))
	{
		nari::system::make_directry( outDir);
		std::cout << "Directry was generated : " << outDir << std::endl;
	}

	std::string	i_fn = outDir + "specificity.csv";
	FILE	*fp;
	fp = fopen( i_fn.c_str(), "w");

	fprintf( fp, "number of times,J.I,specificity\n");

	for( int i = 0; i < ngene; i++){
		Arbitrarily_Shaped_Object( myu, U, evalue, Xi_out, pe, te, &idum);
	
		////任意形状の出力
		//i_fn << outDir << "test" << i << ".raw";
		//Xi_.save_file_bin( i_fn.str());
		//i_fn.str("");

		case_max_ji = 0.0;
		for( int n = 0; n < ncase; n++){
			cross = 0;
			join = 0;
			for( int p = 0; p < pe; p++){
				if( X( p, n) < 0.0 || Xi_out[p] < 0.0){
					join++;
					if( X( p, n) < 0.0 && Xi_out[p] < 0.0)	cross++;
				}
			}
			ji = double(cross)/(double)join;
			if( ji > case_max_ji)	case_max_ji = ji;
		}

		specificity += case_max_ji;
		if( case_max_ji > max_ji)	max_ji = case_max_ji;
		if( case_max_ji < min_ji)	min_ji = case_max_ji;
		//std::cout << "case_max_ji = " << case_max_ji << std::endl;
		fprintf( fp, "%d,%f,%f\n", i + 1 , case_max_ji, specificity / (double)(i+1));
		std::cout << "specificity" << i + 1 << " = " << specificity / (double)(i+1) << std::endl;
	}
	specificity /= (double)ngene;
	fprintf( fp, "max_J.I.,%f\n", max_ji);
	fprintf( fp, "min_J.I.,%f\n", min_ji);
	fprintf( fp, "specificity,%f\n", specificity);
	std::cout << "specificity = " << specificity << std::endl;
	fclose(fp);
}