#include<fstream>
#include<string>
#include<iostream>
#include<list>
#include<vector>
#include<algorithm>
#include<deque>
#include<queue>
#include<numeric>
#include<iomanip>
#include<sstream>
#include<direct.h>
#include<sys/stat.h>

#include "model.h"

namespace saito{
	template< class T >
	void model< T >::get_list_txt( std::list< std::string > &list_txt, std::string path_list ){
		list_txt.clear();
		std::ifstream ifs( path_list, std::ios::in );
		//fail to open file
		if( !ifs ){
			std::cerr << "Cannot open file: " << path_list << std::endl;
			std::abort();
		}else{
			std::string buf;
			while( std::getline( ifs, buf ) ) {
				//skip empty strings
				if( buf.length() ){
					list_txt.push_back( buf );
				}
			}
		}
		ifs.close();
	}

	template< class T >
	long model< T >::get_file_size( std::string filename )
	{
		FILE *fp;
		struct stat st ;
		if( fopen_s(&fp, filename.c_str(), "rb") != 0 ){
			std::cerr << "Cannot open file: " << filename << std::endl;
			std::abort();
		}
		fstat( _fileno(fp), &st );
		fclose( fp );
		return st.st_size;
	}

	template< class T >
	void model< T >::read_bin( T *p, size_t num, std::string filename ){
		FILE *fp;
		if( fopen_s(&fp, filename.c_str(), "rb") != 0 ){
			std::cerr << "Cannot open file" << filename << std::endl;
			std::abort();
		}
		fread( p, sizeof( *p ), num, fp );
		fclose( fp );
	}

	template< class T >
	void model< T >::write_bin( const T *p, size_t num, std::string file_name ){
		FILE *fp;
		fopen_s(&fp, file_name.c_str(), "wb");
		fwrite( p, sizeof( *p ), num, fp );
		fclose( fp );
	}
	
	template< class T >
	void model< T >::write_text( const T *p, size_t num, std::string file_name ){
		FILE *fp;
		fopen_s(&fp, file_name.c_str(), "w");
		for(int i=0;i<num;i++){
			fprintf(fp,"%lf\n",p[i]);
		}
		fclose( fp );
	}

	template< class T >
	void model< T >::PCA( const std::string &data_list_txt ){

		/* get file list */
		std::list< std::string > file_list;
		get_list_txt( file_list, data_list_txt );

		/* load bin */
		const std::size_t num_of_data = file_list.size();
		num_of_dimX = get_file_size(*file_list.begin()) / sizeof( T );
		Mat_t X(num_of_dimX, num_of_data);
		{
			FILE *fp;
			T *x = &X(0, 0);
			for( auto s = file_list.begin(); s != file_list.end(); ++s ){
				fopen_s(&fp, (*s).c_str(), "rb");
				x += fread( x, sizeof( T ), num_of_dimX, fp );
				fclose( fp );
			}
		}

		//X_mean = mean(X, 2)
		X_mean = X.rowwise().mean();

		//X = bsxfun(@minus, X, M.mean)
		X.colwise() -= X_mean;

		//C = X' * X;
		Mat_t C = X.transpose() * X;

		//[V,D] = EIG(C)
		Mat_t V;
		Vec_t lambda;
		{
			Eigen::SelfAdjointEigenSolver< Mat_t > eig( C );
			const T *it = std::find_if( &eig.eigenvalues()( 1 ), &eig.eigenvalues()( 0 ) + num_of_data, std::bind2nd(std::greater< T >(), static_cast< T >( 0 ) ) );
			num_of_dimY = num_of_data - static_cast< std::size_t >( it - &eig.eigenvalues()( 1 ) ) - 1;
			lambda.resize( num_of_dimY, 1 );
			V.resize( num_of_data, num_of_dimY );
			for( std::size_t s = num_of_data - 1, d = 0; d != num_of_dimY; --s, ++d ){
				lambda( d ) = eig.eigenvalues()( s, 0 );
				V.col( d ) = eig.eigenvectors().col( s );
			}
		}

		//M = X * M.V * invsqrt(D)
		M = X * V;
		M *= lambda.array().sqrt().matrix().asDiagonal().inverse();

		//Compute CCR
		CCR.resize( num_of_dimY, 1 );
		CCR( 0 ) = lambda( 0 );
		for( std::size_t s = 1; s != num_of_dimY; ++s ){
			CCR( s ) = CCR( s - 1 ) + lambda( s );
		}
		for( std::size_t s = 0; s != num_of_dimY; ++s ){
			CCR( s ) /= CCR( num_of_dimY - 1 );
		}
		//Eigen ÇÃ Division assignment ââéZéq (/=) ÇÕÇ«Ç§Ç‚ÇÁèÊéZÇ…íuä∑Ç≥ÇÍÇƒÇ¢ÇÈÇÁÇµÇ¢
		//Ç©ÇÁÇ±Ç§Ç∑ÇÈÇ∆åÎç∑Ç™ëÂÇ´Ç¢
		//CCR /= CCR( num_of_dimY - 1 );

		//unbiased
		const T normalizer = static_cast< T >( num_of_data - 0 );
		Eval = lambda / normalizer;

		isModelExists = true;
		isWeighted = false;

	}
	template< class T >
	void model< T >::PCA( Mat_t X ){

		/* load bin */
		const std::size_t num_of_data = X.cols();
		num_of_dimX = X.rows();

		//X_mean = mean(X, 2)
		X_mean = X.rowwise().mean();

		//X = bsxfun(@minus, X, M.mean)
		X.colwise() -= X_mean;

		//C = X' * X;
		Mat_t C = X.transpose() * X;

		//[V,D] = EIG(C)
		Mat_t V;
		Vec_t lambda;
		{
			Eigen::SelfAdjointEigenSolver< Mat_t > eig( C );
			const T *it = std::find_if( &eig.eigenvalues()( 1 ), &eig.eigenvalues()( 0 ) + num_of_data, std::bind2nd(std::greater< T >(), static_cast< T >( 0 ) ) );
			num_of_dimY = num_of_data - static_cast< std::size_t >( it - &eig.eigenvalues()( 1 ) ) - 1;
			lambda.resize( num_of_dimY, 1 );
			V.resize( num_of_data, num_of_dimY );
			for( std::size_t s = num_of_data - 1, d = 0; d != num_of_dimY; --s, ++d ){
				lambda( d ) = eig.eigenvalues()( s, 0 );
				V.col( d ) = eig.eigenvectors().col( s );
			}
		}

		//M = X * M.V * invsqrt(D)
		M = X * V;
		M *= lambda.array().sqrt().matrix().asDiagonal().inverse();

		//Compute CCR
		CCR.resize( num_of_dimY, 1 );
		CCR( 0 ) = lambda( 0 );
		for( std::size_t s = 1; s != num_of_dimY; ++s ){
			CCR( s ) = CCR( s - 1 ) + lambda( s );
		}
		for( std::size_t s = 0; s != num_of_dimY; ++s ){
			CCR( s ) /= CCR( num_of_dimY - 1 );
		}
		//Eigen ÇÃ Division assignment ââéZéq (/=) ÇÕÇ«Ç§Ç‚ÇÁèÊéZÇ…íuä∑Ç≥ÇÍÇƒÇ¢ÇÈÇÁÇµÇ¢
		//Ç©ÇÁÇ±Ç§Ç∑ÇÈÇ∆åÎç∑Ç™ëÂÇ´Ç¢
		//CCR /= CCR( num_of_dimY - 1 );

		//unbiased
		const T normalizer = static_cast< T >( num_of_data - 0 );
		Eval = lambda / normalizer;

		isModelExists = true;
		isWeighted = false;

	}



	template< class T >
	void model< T >::WPCA( const std::string &data_list_txt, T alpha ){

		/* get file list */
		std::list< std::string > file_list;
		get_list_txt( file_list, data_list_txt );

		/* load bin */
		const std::size_t num_of_data = file_list.size();
		num_of_dimX = get_file_size(*file_list.begin()) / sizeof( T );
		Mat_t X(num_of_dimX, num_of_data);
		{
			FILE *fp;
			T *x = &X(0, 0);
			for( auto s = file_list.begin(); s != file_list.end(); ++s ){
				fopen_s(&fp, (*s).c_str(), "rb");
				x += fread( x, sizeof( T ), num_of_dimX, fp );
				fclose( fp );
			}
		}

		//compute weights 
		//w = mean( 1./(1+exp(alpha.*abs(X))), 2 );
		w = X.unaryExpr([alpha](T elem){return 1/(1+exp(alpha*std::abs(elem)));}).rowwise().mean();

		//compute weights 
		//w = bsxfun(@times, X, w);
		X.cwiseProduct( w.replicate(1,X.cols()) );

		//X_mean = mean(X, 2)
		X_mean = X.rowwise().mean();

		//X = bsxfun(@minus, X, M.mean)
		X.colwise() -= X_mean;

		//C = X' * X;
		Mat_t C = X.transpose() * X;

		//[V,D] = EIG(C)
		Mat_t V;
		Vec_t lambda;
		{
			Eigen::SelfAdjointEigenSolver< Mat_t > eig( C );
			const T *it = std::find_if( &eig.eigenvalues()( 1 ), &eig.eigenvalues()( 0 ) + num_of_data, std::bind2nd(std::greater< T >(), static_cast< T >( 0 ) ) );
			num_of_dimY = num_of_data - static_cast< std::size_t >( it - &eig.eigenvalues()( 1 ) ) - 1;
			lambda.resize( num_of_dimY, 1 );
			V.resize( num_of_data, num_of_dimY );
			for( std::size_t s = num_of_data - 1, d = 0; d != num_of_dimY; --s, ++d ){
				lambda( d ) = eig.eigenvalues()( s, 0 );
				V.col( d ) = eig.eigenvectors().col( s );
			}
		}

		//M = X * M.V * invsqrt(D)
		M = X * V;
		M *= lambda.array().sqrt().matrix().asDiagonal().inverse();

		//Compute CCR
		CCR.resize( num_of_dimY, 1 );
		CCR( 0 ) = lambda( 0 );
		for( std::size_t s = 1; s != num_of_dimY; ++s ){
			CCR( s ) = CCR( s - 1 ) + lambda( s );
		}
		for( std::size_t s = 0; s != num_of_dimY; ++s ){
			CCR( s ) /= CCR( num_of_dimY - 1 );
		}
		//Eigen ÇÃ Division assignment ââéZéq (/=) ÇÕÇ«Ç§Ç‚ÇÁèÊéZÇ…íuä∑Ç≥ÇÍÇƒÇ¢ÇÈÇÁÇµÇ¢
		//Ç©ÇÁÇ±Ç§Ç∑ÇÈÇ∆åÎç∑Ç™ëÂÇ´Ç¢
		//CCR /= CCR( num_of_dimY - 1 );

		//unbiased
		const T normalizer = static_cast< T >( num_of_data - 0 );
		Eval = lambda / normalizer;

		isModelExists = true;
		isWeighted = true;
	}


	template< class T >
	void model< T >::output( const std::string outDir ){
		if(isModelExists){
			std::cout << "mkdir: " << outDir << std::endl;
			_mkdir( outDir.c_str() );
			{
				for( int j = 0; j != num_of_dimY; ++j){
					std::ostringstream name;
					name << outDir << "/vect_" << std::setw(4) << std::setfill('0') << j << ".vect";
					write_bin( &M(0, j), num_of_dimX, name.str().c_str() );
				}
				write_bin( &X_mean(0),   num_of_dimX, outDir + "/mean.vect" );
				write_bin( &Eval(0),     num_of_dimY, outDir + "/eval.vect" );
				write_bin( &CCR(0),      num_of_dimY, outDir + "/CCR.vect" );
				write_text( Eval.data(),     num_of_dimY, outDir + "/eval.txt" );
				write_text( CCR.data(),      num_of_dimY, outDir + "/CCR.txt" );
				if(isWeighted){
					write_bin( &w(0), num_of_dimX, outDir + "/weight.vect" );
				}
				sqrtL = Eval.array().sqrt().matrix();
				invsqrtL = Eval.array().sqrt().inverse().matrix();
			}

		}else{
			std::cerr << "A model is not stored in the class." << std::endl;
		}
	}

	template< class T >//äÓñ{ìIÇ»ì«Ç›çûÇ›
	void model< T >::load( const std::string modelDir ){
		const std::size_t num_of_basis = get_file_size( modelDir + "/eval.vect" ) / sizeof( T );
		load_common( modelDir, num_of_basis );
	}

	template< class T >//ó›êœäÒó^ó¶Ç≈åàÇﬂÇÈéû
	void model< T >::load_with_CCR( const std::string modelDir, const T ccr ){
		const std::size_t num_of_max_basis = get_file_size( modelDir + "/eval.vect" ) / sizeof( T );
		T *CCR_tmp = new T[ num_of_max_basis ];
		read_bin( CCR_tmp, num_of_max_basis, modelDir + "/CCR.vect" );
		T *it = std::find_if( CCR_tmp, CCR_tmp + num_of_max_basis - 1, std::bind2nd(std::greater_equal< T >(), ccr ) );
		const std::size_t num_of_basis = static_cast< std::size_t >( it - CCR_tmp ) + 1;
		load_common( modelDir, num_of_basis );
	}

	template< class T >//éüå≥ñ{êîÇ≈åàÇﬂÇÈéû
	void model< T >::load_with_N( const std::string modelDir, const std::size_t num_of_basis ){
		load_common( modelDir, num_of_basis );
	}

	template< class T >
	void model< T >::load_common( const std::string modelDir, const std::size_t num_of_basis ){
		load_base( modelDir, num_of_basis );
	}

	template< class T >
	void model< T >::load_base( const std::string modelDir, const std::size_t num_of_basis ){
		num_of_dimY = num_of_basis;
		num_of_dimX = get_file_size( modelDir + "/mean.vect" ) / sizeof( T );
		/*resize matrices or vectors*/
		M.resize( num_of_dimX, num_of_dimY );
		X_mean.resize( num_of_dimX );
		Eval.resize( num_of_dimY );
		sqrtL.resize( num_of_dimY) ;
		invsqrtL.resize( num_of_dimY );
		CCR.resize( num_of_dimY );
		/*--------------------------*/
		for(std::size_t j = 0; j != num_of_dimY; ++j ){
			std::ostringstream name;
			name << modelDir << "/vect_" << std::setw(4) << std::setfill('0') << j << ".vect";
			read_bin( &M(0, j), num_of_dimX, name.str() );
		}
		read_bin( &X_mean(0), num_of_dimX, modelDir + "/mean.vect" );
		read_bin( &Eval(0),   num_of_dimY, modelDir + "/eval.vect" );
		read_bin( &CCR(0),    num_of_dimY, modelDir + "/CCR.vect" );
		sqrtL = Eval.array().sqrt().matrix();
		invsqrtL = Eval.array().sqrt().inverse().matrix();
		isModelExists = true;
	}


	template< class T >
	void model< T >::disp( ){
		if(isModelExists){
			std::cout << "------------------------------------" << std::endl;
			std::cout << "Number of dim X: " << num_of_dimX << std::endl;
			std::cout << "Number of dim Y: " << num_of_dimY << std::endl;
			std::cout << "------------------------------------" << std::endl;
			std::cout << "M ="        << std::endl << M        << std::endl << std::endl;
			std::cout << "X_mean ="   << std::endl << X_mean   << std::endl << std::endl;
			std::cout << "Eval ="     << std::endl << Eval     << std::endl << std::endl;
			std::cout << "sqrtL ="    << std::endl << sqrtL    << std::endl << std::endl;
			std::cout << "invsqrtL =" << std::endl << invsqrtL << std::endl << std::endl;
			std::cout << "CCR ="      << std::endl << CCR      << std::endl << std::endl;
		}else{
			std::cerr << "A model is not stored in the class." << std::endl;
		}

	}

	template< class T >
	void model< T >::out_of_sample( Mat_t &Y, const Mat_t &X ){
		Y = M.transpose() * (X.colwise() - X_mean);
	}

	template< class T >
	void model< T >::out_of_sample_normal( Mat_t &Alpha, const Mat_t &X ){
		Alpha = invsqrtL.asDiagonal() * ( M.transpose() * (X.colwise() - X_mean) );
	}

	template< class T >
	void model< T >::pre_image( Mat_t &X, const Mat_t &Y ){
		X = ( M * Y ).colwise() + X_mean;
	}

	template< class T >
	void model< T >::pre_image_normal( Mat_t &X, const Mat_t &Alpha ){
		X = ( M * ( sqrtL.asDiagonal() * Alpha ) ).colwise() + X_mean;
	}

	template< class T >
	void model< T >::reconstruction( Mat_t &Xout, const Mat_t &Xin ){
		Xout = ( M * ( M.transpose() * ( Xin.colwise() - X_mean ) ) ).colwise() + X_mean;
	}

	template class model<double>;
	template class model<float>;

};
