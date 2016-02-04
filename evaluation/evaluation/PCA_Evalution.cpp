#include <iostream>
#include "naristd.h"
#include "PCA.h"
#include "JI.h"
#include "random.h"
#include "PCA_Evalution.h"
#include "other/narimhd.h"
#include "nariinfocontroller.h"

//入力ファイルの型
typedef double _IN_TYPE_;

//1:テスト症例リスト
//2:固有ベクトル
//3:固有値
//4:平均形状
//5:出力フォルダ
int main(int argc, char *argv[]){
	if(argc != 2){
		std::cout << "usage : .exe [path_list]" << std::endl;
		exit(EXIT_FAILURE);
	}

	//input information
	std::string	settingINI = "setting.ini";
	if(argc == 2) settingINI = argv[1];
	nari::infocontroller infoDir;
	infoDir.load(settingINI);

	//元の肝臓システムからの入力
	const std::string	dirOut = nari::file::add_delim( infoDir["out"]);
	const std::string	caseTxt = infoDir["case"];	//テスト症例リスト
	const std::string	evectTxt = infoDir["evect"];	//使用するモデルの固有ベクトルリスト
	const std::string	dirModel = nari::file::add_delim( infoDir["model"]);	//モデルが格納されているパス
	const std::string	pathEval = dirModel + "eval.vect";		//モデルの固有値
	const std::string	pathAve = dirModel + "average.vect";		//モデルの平均形状



	//Load Filepath by InputText
	//**************テスト症例のロード*************
	std::vector<std::string> inputList;	//テストデータリスト
	int dimension;	//次元数
	const int ncase = nari::file::get_file_list<_IN_TYPE_>( caseTxt, inputList, dimension);	//症例数，次元数，テストデータリスト取得

	mist::matrix <_IN_TYPE_> X(dimension, ncase);
	nari::matrix::load_matrix_bin_col_vector(inputList, X);

	//**************固有ベクトルのロード*************
	const int te = nari::file::get_file_list<_IN_TYPE_>( evectTxt, inputList, dimension);	//症例数，次元数，テストデータリスト取得

	mist::matrix <_IN_TYPE_> U(dimension, te);	//固有ベクトル行列
	nari::matrix::load_matrix_bin_col_vector(inputList, U);

	//**************固有値のロード*************
	nari::vector<double>	evalue( te, 0.0);
	evalue.load_file_bin( pathEval);

	//**************平均形状のロード*************
	nari::vector<double>	myu( dimension, 0.0);
	myu.load_file_bin( pathAve);

	std::cout << "ncase : " << ncase << std::endl;
	std::cout << "Dimension of Input Vector : " << dimension << std::endl;
	std::cout << "軸数" << te << std::endl;

	generalization( X, myu, U, dimension, te, ncase, dirOut);
	specificity( X, myu, U, evalue, dimension, te, ncase, 1000, -1, dirOut);

	////正規乱数のヒストグラムを作成
	//nari::vector<int>	hist( 100,0);
	//long	idum = 1;
	//int j = 0;
	//long	n = 100000000;
	//double	a = 0.0;
	//double	ave = 0.0;
	//double	var = 0.0;
	//double	width = 0.1;
	//double	h_width = 0.05;
	//for( long i = 0; i < n; i++){
	//	 a = normal( &idum);
	//	 ave += a;
	//	// j = 0;
	//	// for( double b = -3.0; b <= 3.0; b += width){
	//	//	if( b - h_width < a && a <= b + h_width)	hist[j]++;	
	//	//	j++;
	//	//}
	//}
	//ave /= (double)n;
	//idum = 1;
	//for( long i = 0; i < n; i++){
	//	a = normal( &idum) - ave;
	//	var += a * a;
	//}
	//var /= (double)n;
	//std::cout << "average : " << ave << "\tvar : " << var << std::endl;
	//hist.save_file_txt("hist.txt");

	
	////任意形状の出力
	//long	idum = 1;
	//nari::vector<double>	X_out( dimension, 0.0);
	//nari::vector<unsigned char>	X_( dimension, 0);
	//std::ostringstream	path;
	//for( int i = 0; i < 10; i++){
	//	Arbitrarily_Shaped_Object( myu, U, evalue, X_out, dimension, te, &idum);
	//	path << outDir << "/shape" << i << ".raw";
	//	for( int p = 0; p < dimension; p++){
	//		if( X_out[p] > 0.0)	X_[p] = 0;
	//		else	X_[p] = 1;
	//	}
	//	X_.save_file_bin( path.str());
	//	path.str("");
	//}
	return	0;
}
