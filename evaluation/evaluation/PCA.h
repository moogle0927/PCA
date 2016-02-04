/*********************************/
/*       By uchida 2010/2/12     */
/*********************************/


#ifndef __NARI_PCA__
#define __NARI_PCA__

#include "naristd.h"

namespace nari
{
	namespace PCA
	{
		template <class IN_TYPE>
		void perform_PCA( const mist::matrix <IN_TYPE> &X, mist::matrix <double> &evect, mist::matrix <double> &eval, mist::matrix <double> &average, std::string outDir){

			const int dimension = X.rows();
			const int ncase = X.cols();

			//����ۂɂ��悤
			eval.clear();
			evect.clear();
			average.resize(dimension,1);

			//�m�F�p
			std::cout << "dimension:" << dimension << "\tncase:" << ncase << std::endl;

			//X_in.clear();
			//Calculate Average Vectors
			std::cout << "Calculate Average" << std::endl;
			for(int j = 0; j < dimension; j++)
			{
				for(int i = 0; i < ncase; i++)
				{
					average( j, 0) += X( j, i);
				}
				average( j, 0) /= (double)ncase;
			}


			nari::matrix::save_matrix_bin_col_vector(outDir + "\\average.vect", average);

			if( dimension < ncase )
			{	//EVD
				std::cout << "EVD" << std::endl;
				//Calculate CovMatrix
				mist::matrix <double> cov_mat(dimension, dimension);
				for(int r=0; r<dimension; r++)
				{
					for(int c=r; c<dimension; c++)
					{
						for(int n=0; n<ncase; n++)
						{
							cov_mat(r, c) += (X(r, n) - average(r, 0)) * (X(c, n) - average(c, 0)) / (ncase-1.0);
						}
						cov_mat(c, r) = cov_mat(r, c);
					}
				}
				//cov_mat = X * X.t() / ncase;

				//Eigenvalue Decomposition
				mist::eigen(cov_mat, eval, evect, mist::matrix_style::sy);

				/*************** Output eval & evec ***************/
				//eval txt
				std::cout << "Eigenvalues" << std::endl;
				{
					FILE *fp_log_txt = nari::file::file_open(outDir + "\\log.txt","wb");
					FILE *fp_eval_txt = nari::file::file_open(outDir + "\\eval.txt","wb");
					double cr=0,ccr=0,total_vari=0;
					for(int i=0; i<eval.rows(); i++)
					{
						total_vari += eval(i,0);
					}
					for(int i=0; i<eval.rows(); i++)
					{
						cr = eval(i,0) / total_vari;
						ccr += cr;
						printf("%4d : %lf (%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
						fprintf(fp_log_txt,"%4d : %lf(%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
						fprintf(fp_eval_txt, "%lf\n", eval(i, 0));
					}
					fclose(fp_log_txt);
					fclose(fp_eval_txt);
				}

				//eval binary
				nari::matrix::save_matrix_bin_col_vector(outDir + "\\eval.vect", eval);

				//eigenvector
				for(int i=0; i<dimension; i++)
				{
					nari::matrix::save_matrix_bin_col_vector(outDir + "\\evect" + nari::string::tostring(i, "%04d") + ".vect", evect, i);
				}
				/*************** END OUTPUT ***************/
			}
			else
			{	//SVD
				std::cout << "SVD" << std::endl;
				mist::matrix <double> inner_mat(ncase,ncase);
				{	//���֍s��̉��Z
					int counter=0;
					for(int i=0; i < ncase; i++)
					{
						counter++;
						printf("%d / %d\r", counter, ncase);
						for(int j=i; j<ncase; j++)
						{
							for(int k=0; k<dimension; k++)
							{
								inner_mat(i,j) += (X(k, i) - average(k, 0)) * (X(k, j) - average(k, 0)) / (ncase-1.0);
							}
							inner_mat(j,i) = inner_mat(i,j);
						}
					}
				}

				//Eigenvalue Decomposition
				std::cout << "EigenValue Decomposition" << std::endl;
				mist::matrix <double> evect_innerp;
				mist::eigen(inner_mat, eval, evect_innerp, mist::matrix_style::sy);

				//Resize
				std::cout << "Resize..." << std::endl;
				{
					evect.resize(dimension, ncase-1);
					int count=0;
					for(int i=0; i<ncase-1; i++)
					{
						printf("%d / %d\r", count+1, ncase-1);
						count++;
						for(int k=0; k<dimension; k++)
						{
							for(int j=0; j<ncase; j++)
							{
								evect(k,i) += (X(k,j) - average(k, 0)) * evect_innerp(j,i);
							}
						}
					}
				}

				//Normalize length
				for(int i=0; i<ncase-1; i++){	//�Ō�̌ŗL�x�N�g���͗v��Ȃ��q�D
					double length = 0;
					for(int j=0; j<dimension; j++){
						length += evect(j,i) * evect(j,i);
					}
					length = sqrt(length);
					for(int j=0; j<dimension; j++){
						evect(j,i) /= length;
					}
					std::cout << length << std::endl;
				}

				/*************** Output eval & evec ***************/
				//log txt
				std::cout << "Eigenvalues" << std::endl;
				{
					FILE *fp_log_txt = nari::file::file_open(outDir + "\\log.txt", "wb");
					double cr=0,ccr=0,total_vari=0;
					for(int i=0; i<ncase-1; i++)
					{
						total_vari += eval(i,0);
					}
					for(int i=0; i<ncase-1; i++)
					{
						cr = eval(i,0) / total_vari;
						ccr += cr;
						printf("%4d : %lf (%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
						fprintf(fp_log_txt,"%4d : %lf(%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
					}
					fclose(fp_log_txt);
				}

				//eval txt
				{
					FILE *fp_eval_txt = nari::file::file_open(outDir + "\\eval.txt","wb");
					for(int i=0; i<ncase-1; i++)
					{
						fprintf(fp_eval_txt,"%lf\n", eval(i,0));
					}
					fclose(fp_eval_txt);
				}

				//eval binary
				nari::matrix::save_matrix_bin_col_vector(outDir + "\\eval.vect", eval);

				//eigenvector
				for(int i=0; i<ncase-1; i++){
					nari::matrix::save_matrix_bin_col_vector(outDir + "\\evect" + nari::string::tostring(i, "%04d") + ".vect", evect, i);
				}
				/*************** END OUTPUT ***************/
			}
		}
	}
}

#endif


		////���e�̂�
		//template <class I>
		//void projection_PCA(const mist::matrix<I> &X, const mist::matrix<double> &average, const mist::matrix<double> &evect, mist::matrix<double> &lowCoor, const std::string outDir, const std::vector<std::string> inputList)
		//{
		//	const int ncase = X.cols();	//�Ǘᐔ
		//	const int dimension = X.rows();	//������
		//	const int basis = evect.cols();	//�厲��
		//	lowCoor.resize(basis, ncase);	//�o�͍s��̊m��

		//	//�Ǘᖈ��
		//	int count = 0;
		//	for(int n=0; n<ncase; n++)
		//	{
		//		printf("%4d / %4d\r", count+1, ncase);	//�i�s��
		//		count++;
		//		mist::matrix<double> X_zeroMean(dimension, 1);
		//		//�厲����
		//		for(int d=0; d<dimension; d++)
		//		{
		//			X_zeroMean(d, 0) = X(d, n) - average(d, 0);	//���ς�����
		//			for(int b=0; b<basis; b++)
		//			{
		//				lowCoor(b, n) += evect(d, b) * X_zeroMean(d, 0);
		//			}
		//		}
		//	}

		//	//txt�֏o��
		//	{
		//		MakeSureDirectoryPathExists((outDir + "\\lowCoor\\lowCoor.txt").c_str());
		//		std::ofstream ofs((outDir + "\\lowCoor\\lowCoor.txt").c_str());
		//		ofs << lowCoor.t() << std::endl;
		//	}

		//	//�o�C�i���t�@�C���֏o��
		//	{
		//		for(int n=0; n<ncase; n++)
		//		{
		//			nari::matrix::save_matrix_bin_col_vector(outDir+ "\\lowCoor\\lowCoor_" + inputList[n].substr(inputList[n].length()-4), lowCoor, n);
		//		}
		//	}
		//}

		////�t���e
		//template <class I>
		//void backProjection_PCA(const mist::matrix<double> &Alpha, const mist::matrix<double> &average, const mist::matrix<double> &evect, mist::matrix<I> &XReconst)
		//{
		//	const int ncase = Alpha.cols();	//�Ǘᐔ
		//	const int basis = Alpha.rows();	//�p����厲��
		//	const int dimension = evect.rows();	//������
		//	if(basis > evect.cols())
		//	{
		//		std::cout << "Error in \"backProjection_PCA\"" << std::endl;
		//		exit(EXIT_FAILURE);
		//	}
		//	XReconst.resize(dimension, ncase);

		//	int count = 0;
		//	for(int n=0; n<ncase; n++)
		//	{
		//		printf("%4d / %4d\r", count+1, ncase);	//�i�s��
		//		count++;
		//		for(int d=0; d<dimension; d++)
		//		{
		//			XReconst(d, n) = average(d, 0);	//���ς��Z�b�g
		//			for(int b=0; b<basis; b++)
		//			{
		//				XReconst(d, n) += Alpha(b, n) * evect(d, b);	//�t���e
		//			}
		//		}
		//	}
		//}

		////���e�E�t���e
		//template <class I>
		//mist::matrix<double> projection_backProjection_PCA(const mist::matrix<I> &X, const mist::matrix<double> &average, const mist::matrix<double> &evect)
		//{
		//	const int ncase = X.cols();	//�Ǘᐔ
		//	const int dimension = X.rows();	//������
		//	const int basis = evect.cols();	//�厲��
		//	mist::matrix<double> reconst(dimension, ncase);	//�o�͍s��̊m��
		//	mist::matrix<double> lowCoor(basis, ncase);	//�听�����_�s��

		//	//�Ǘᖈ��
		//	int count = 0;
		//	for(int n=0; n<ncase; n++)
		//	{
		//		printf("%4d / %4d\r", count+1, ncase);	//�i�s��
		//		count++;
		//		mist::matrix<double> X_zeroMean(dimension, 1);
		//		//�厲����
		//		for(int d=0; d<dimension; d++)
		//		{
		//			X_zeroMean(d, 0) = X(d, n) - average(d, 0);	//���ς�����
		//			for(int b=0; b<basis; b++)
		//			{
		//				lowCoor(b, n) += evect(d, b) * X_zeroMean(d, 0);
		//			}
		//		}
		//	}

		//	count = 0;
		//	for(int n=0; n<ncase; n++)
		//	{
		//		printf("%4d / %4d\r", count+1, ncase);	//�i�s��
		//		count++;
		//		for(int d=0; d<dimension; d++)
		//		{
		//			reconst(d, n) = average(d, 0);	//���ς��Z�b�g
		//			for(int b=0; b<basis; b++)
		//			{
		//				reconst(d, n) += lowCoor(b, n) * evect(d, b);	//�t���e
		//			}
		//		}
		//	}

		//	return reconst;
		//}

		////�ŗL�l�E��^���E�ݐϊ�^�����e�L�X�g�ɕۑ�
		//inline void save_eigenvalue_text(const mist::matrix<double> &eval, const std::string path)
		//{
		//	const int basis = eval.rows();
		//	double sumEval = 0;
		//	for(int b=0; b<basis; b++)
		//	{
		//		sumEval += eval(b, 0);
		//	}

		//	FILE * fp = nari::file::file_open(path, "wb");
		//	double cr=0, ccr=0;
		//	for(int b=0; b<basis; b++)
		//	{
		//		cr = eval(b, 0) / sumEval;
		//		ccr += cr;
		//		fprintf(fp, "%4d : %lf (%lf , %lf)\n", b+1, eval(b, 0), cr, ccr);
		//		printf("%4d : %lf (%lf , %lf)\n", b+1, eval(b, 0), cr, ccr);
		//	}
		//	fclose(fp);
		//}


template <class IN_TYPE>
void perform_PCA( const mist::matrix <IN_TYPE> &X, mist::matrix <double> &evect, mist::matrix <double> &eval, mist::matrix <double> &average, std::string outDir){

	const int dimension = X.rows();
	const int ncase = X.cols();

	//����ۂɂ��悤
	eval.clear();
	evect.clear();
	average.resize(dimension,1);

	//�m�F�p
	std::cout << "dimension:" << dimension << "\tncase:" << ncase << std::endl;

	//X_in.clear();
	//Calculate Average Vectors
	std::cout << "Calculate Average" << std::endl;
	for(int j = 0; j < dimension; j++)
	{
		for(int i = 0; i < ncase; i++)
		{
			average( j, 0) += X( j, i);
		}
		average( j, 0) /= (double)ncase;
	}


	nari::matrix::save_matrix_bin_col_vector(outDir + "\\average.vect", average);

	if( dimension < ncase )
	{	//EVD
		std::cout << "EVD" << std::endl;
		//Calculate CovMatrix
		mist::matrix <double> cov_mat(dimension, dimension);
		for(int r=0; r<dimension; r++)
		{
			for(int c=r; c<dimension; c++)
			{
				for(int n=0; n<ncase; n++)
				{
					cov_mat(r, c) += (X(r, n) - average(r, 0)) * (X(c, n) - average(c, 0)) / ncase;
				}
				cov_mat(c, r) = cov_mat(r, c);
			}
		}
		//cov_mat = X * X.t() / ncase;

		//Eigenvalue Decomposition
		mist::eigen(cov_mat, eval, evect, mist::matrix_style::sy);

		/*************** Output eval & evec ***************/
		//eval txt
		std::cout << "Eigenvalues" << std::endl;
		{
			FILE *fp_log_txt = nari::file::file_open(outDir + "\\log.txt","wb");
			FILE *fp_eval_txt = nari::file::file_open(outDir + "\\eval.txt","wb");
			double cr=0,ccr=0,total_vari=0;
			for(int i=0; i<eval.rows(); i++)
			{
				total_vari += eval(i,0);
			}
			for(int i=0; i<eval.rows(); i++)
			{
				cr = eval(i,0) / total_vari;
				ccr += cr;
				printf("%4d : %lf (%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
				fprintf(fp_log_txt,"%4d : %lf(%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
				fprintf(fp_eval_txt, "%lf\n", eval(i, 0));
			}
			fclose(fp_log_txt);
			fclose(fp_eval_txt);
		}

		//eval binary
		nari::matrix::save_matrix_bin_col_vector(outDir + "\\eval.vect", eval);

		//eigenvector
		for(int i=0; i<dimension; i++)
		{
			nari::matrix::save_matrix_bin_col_vector(outDir + "\\evect" + nari::string::tostring(i, "%04d") + ".vect", evect, i);
		}
		/*************** END OUTPUT ***************/
	}
	else
	{	//SVD
		std::cout << "SVD" << std::endl;
		mist::matrix <double> inner_mat(ncase,ncase);
		{	//���֍s��̉��Z
			int counter=0;
			for(int i=0; i < ncase; i++)
			{
				counter++;
				printf("%d / %d\r", counter, ncase);
				for(int j=i; j<ncase; j++)
				{
					for(int k=0; k<dimension; k++)
					{
						inner_mat(i,j) += (X(k, i) - average(k, 0)) * (X(k, j) - average(k, 0)) / ncase;
					}
					inner_mat(j,i) = inner_mat(i,j);
				}
			}
		}

		//Eigenvalue Decomposition
		std::cout << "EigenValue Decomposition" << std::endl;
		mist::matrix <double> evect_innerp;
		mist::eigen(inner_mat, eval, evect_innerp, mist::matrix_style::sy);

		//Resize
		std::cout << "Resize..." << std::endl;
		{
			evect.resize(dimension, ncase-1);
			int count=0;
			for(int i=0; i<ncase-1; i++)
			{
				printf("%d / %d\r", count+1, ncase-1);
				count++;
				for(int k=0; k<dimension; k++)
				{
					for(int j=0; j<ncase; j++)
					{
						evect(k,i) += (X(k,j) - average(k, 0)) * evect_innerp(j,i);
					}
				}
			}
		}

		//Normalize length
		for(int i=0; i<ncase-1; i++){	//�Ō�̌ŗL�x�N�g���͗v��Ȃ��q�D
			double length = 0;
			for(int j=0; j<dimension; j++){
				length += evect(j,i) * evect(j,i);
			}
			length = sqrt(length);
			for(int j=0; j<dimension; j++){
				evect(j,i) /= length;
			}
			std::cout << length << std::endl;
		}

		/*************** Output eval & evec ***************/
		//log txt
		std::cout << "Eigenvalues" << std::endl;
		{
			FILE *fp_log_txt = nari::file::file_open(outDir + "\\log.txt", "wb");
			double cr=0,ccr=0,total_vari=0;
			for(int i=0; i<ncase-1; i++)
			{
				total_vari += eval(i,0);
			}
			for(int i=0; i<ncase-1; i++)
			{
				cr = eval(i,0) / total_vari;
				ccr += cr;
				printf("%4d : %lf (%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
				fprintf(fp_log_txt,"%4d : %lf(%lf , %lf)\n",i+1 , eval(i,0), cr, ccr);
			}
			fclose(fp_log_txt);
		}

		//eval txt
		{
			FILE *fp_eval_txt = nari::file::file_open(outDir + "\\eval.txt","wb");
			for(int i=0; i<ncase-1; i++)
			{
				fprintf(fp_eval_txt,"%lf\n", eval(i,0));
			}
			fclose(fp_eval_txt);
		}

		//eval binary
		nari::matrix::save_matrix_bin_col_vector(outDir + "\\eval.vect", eval);

		//eigenvector
		for(int i=0; i<ncase-1; i++){
			nari::matrix::save_matrix_bin_col_vector(outDir + "\\evect" + nari::string::tostring(i, "%04d") + ".vect", evect, i);
		}
		/*************** END OUTPUT ***************/
	}
}