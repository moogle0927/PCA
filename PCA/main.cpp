#include<iostream>
#include"model.h"
#include<fstream>
#include<vector>


void get_list_txt(std::list< std::string > &list_txt, std::string path_list) {
	list_txt.clear();
	std::ifstream ifs(path_list, std::ios::in);
	//fail to open file
	if (!ifs) {
		std::cerr << "Cannot open file: " << path_list << std::endl;
		std::abort();
	}
	else {
		std::string buf;
		while (std::getline(ifs, buf)) {
			//skip empty strings
			if (buf.length()) {
				list_txt.push_back(buf);
			}
		}
	}
	ifs.close();
}

template< class T >
void write_vector(std::vector<T> &v, const std::string filename) {
	FILE *fp;
	if (fopen_s(&fp, filename.c_str(), "wb") != 0) {
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	fwrite(v.data(), sizeof(T), v.size(), fp);
	fclose(fp);
}


template<typename T>
void write_matrix_raw_and_txt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename)
{
	//////////////////////////////////////////////////////////////
	// Wの書き出し												//
	// rowが隠れ層の数，colが可視層の数							//			
	// 重みの可視化を行う場合は，各行を切り出してreshapeを行う  //
	//////////////////////////////////////////////////////////////
	size_t rows = data.rows();
	size_t cols = data.cols();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Data;
	Data = data;
	std::ofstream fs1(filename + ".txt");
	fs1 << "rows = " << rows << std::endl;
	fs1 << "cols = " << cols << std::endl;
	fs1 << typeid(Data).name() << std::endl;
	fs1.close();
	std::vector<T> save_data(rows * cols);
	Data.resize(rows * cols, 1);
	for (size_t i = 0; i <save_data.size(); i++)
		save_data[i] = Data(i, 0);
	write_vector(save_data, filename + ".raw");
	Data.resize(rows, cols);
}

long get_file_size(std::string filename)
{
	FILE *fp;
	struct stat st;
	if (fopen_s(&fp, filename.c_str(), "rb") != 0) {
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	fstat(_fileno(fp), &st);
	fclose(fp);
	return st.st_size;
}

template< class T >
void read_bin(T *p, size_t num, std::string filename) {
	FILE *fp;
	if (fopen_s(&fp, filename.c_str(), "rb") != 0) {
		std::cerr << "Cannot open file" << filename << std::endl;
		std::abort();
	}
	fread(p, sizeof(*p), num, fp);
	fclose(fp);
}

template< class T >
void write_bin(const T *p, size_t num, std::string file_name) {
	FILE *fp;
	fopen_s(&fp, file_name.c_str(), "wb");
	fwrite(p, sizeof(*p), num, fp);
	fclose(fp);
}

template< class T >
void write_text(const T *p, size_t num, std::string file_name) {
	FILE *fp;
	fopen_s(&fp, file_name.c_str(), "w");
	for (int i = 0; i<num; i++) {
		fprintf(fp, "%lf\n", p[i]);
	}
	fclose(fp);
}

int main(int argn, const char *argc[]){
	if(argn!=4){
		std::cerr << "入力がまちがっています" << std::endl;
		return EXIT_FAILURE;
	}

	saito::model<double> pca;
	pca.PCA(argc[1]);
	//pca.WPCA(argc[1],0.5);
	/*Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > X(2,3);
	X << 1,0,1,0,1,1;
	std::cout<<X << std::endl;
	pca.PCA(X);
	pca.disp();*/
	pca.output(argc[3]);

	/* get file list */
	std::list< std::string > file_list;
	get_list_txt(file_list, argc[2]);


	/* load bin */
	const std::size_t num_of_data = file_list.size();
	
	long num_of_dimX = get_file_size(*file_list.begin()) / sizeof(double);
	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > X(num_of_dimX, num_of_data);
	{
		FILE *fp;
		double *x = &X(0, 0);
		for (auto s = file_list.begin(); s != file_list.end(); ++s) {
			fopen_s(&fp, (*s).c_str(), "rb");
			x += fread(x, sizeof(double), num_of_dimX, fp);
			fclose(fp);
		}
	}

	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Y;

	pca.out_of_sample(Y,X);
	std::string o_file = static_cast<std::string>(argc[3]) + "/mat";
	write_matrix_raw_and_txt(Y, o_file);
	std::ofstream mat(o_file + ".csv");
	for (int i = 0; i < Y.rows(); i++) {
		for (int j = 0; j < Y.cols(); j++) {
			mat << Y(i, j) << ",";
		}
		mat << std::endl;
	}
	
	return EXIT_SUCCESS;
	system("pause");

}
