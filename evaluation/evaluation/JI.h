template< class T1, class T2>
void	comp_image_data( T1 *input, T2 *right, double size, double sigma, double lambda, char *path){
	long	err = 0;
	long	cross = 0;
	long	join = 0;
	double	noise = 0.0;
	FILE	*fp;

	fp = fopen( path, "a");

	for( int i = 0; i < size; i++){
		if(input[i] != right[i] && input[i] == 1)	err++;
		if( err > LONG_MAX)	exit(1);

		if( input[i] == 1 || right[i] == 1){
			join++;
			if( input[i] == right[i])	cross++;
		}
	}
	noise /= (double)(size);
	fprintf( fp, "sigma,%f\n", lambda);
	fprintf( fp, "ramda,%f\n", sigma);
	fprintf( fp, "間違い数,%d\n", err);
	fprintf( fp, "一致度,%f[%%]\n\n", (double)cross / (double)join * 100.0);
	fclose( fp);
}

template< class T1, class T2>
double	comp_image_data( T1 *input1, T2 *input2, double size){
	long	cross = 0;
	long	join = 0;

	for( int i = 0; i < size; i++){

		if( input1[i] == 1 || input2[i] == 1){
			join++;
			if( input1[i] == input2[i])	cross++;
		}
	}
	 return	(double)cross / (double)join;
}


//template< class T1, class T2>
//void	comp_image_data( T1 *input, T2 *right, double size, char *path){
//	long	err = 0;
//	long	loss = 0;
//	long	cross = 0;
//	long	join = 0;
//	long	right_area = 0;
//	double	fnr = 0.0;
//	double	fpr = 0.0;
//	double	ji = 0.0;
//	FILE	*fp;
//
//	fp = fopen( path, "a");
//
//	for( int i = 0; i < size; i++){
//		//if(input[i] != right[i] && input[i] == 1)	err++;
//		//if(input[i] != right[i] && input[i] == 0)	loss++;
//		//if( err > LONG_MAX)	exit(1);
//
//		if( input[i] == 1 || right[i] == 1){
//			join++;
//			if( input[i] == right[i])	cross++;
//			else if( input[i] == 1)	err++;
//			else	loss++;
//			if( right[i] == 1)	right_area++;
//		}
//	}
//	fpr = (double)err / (double)right_area;
//	fnr = (double)loss / (double)right_area;
//	ji = (double)cross / (double)join;
//	//fprintf( fp, "間違い数,%d\n", err);
//	//fprintf( fp, "一致度,%f[%%]\n\n", (double)cross / (double)join * 100.0);
//	fprintf( fp, "%f,%f,%f\n", ji, fpr, fnr);
//	std::cout << "J.I. = " << ji << "\tfpr = " << fpr << "\tfnr = " << fnr << std::endl;
//	fclose( fp);
//}

//マスクが無い場合のJI測定
template< class T1, class T2>
void	comp_image_data( T1 *input, T2 *right, int size, std::string &path){
	long	err = 0;
	long	loss = 0;
	long	cross = 0;
	long	join = 0;
	long	right_area = 0;
	long	bkg = 0;
	double	fnr = 0.0;
	double	fpr = 0.0;
	double	ji = 0.0;
	FILE	*fp;

	fp = fopen( path.c_str(), "a");

	for( int i = 0; i < size; i++){
		if( input[i] == 1 || right[i] == 1){
			join++;
			if( input[i] == right[i])	cross++;
			else if( input[i] == 1)	err++;
			else	loss++;
			if( right[i] == 1)	right_area++;
		}
	}
	bkg = size - right_area;
	fpr = (double)err / (double)bkg;
	fnr = (double)loss / (double)right_area;
	ji = (double)cross / (double)join;
	fprintf( fp, "%f,%f,%f\n", ji, fpr, fnr);
	std::cout << "J.I. = " << ji << "\tfpr = " << fpr << "\tfnr = " << fnr << std::endl;
	fclose( fp);
}

//JIをreturnする，マスクなしの場合
template< class T1, class T2>
double	comp_image_data( T1 *input, T2 *right, int size, std::string &path){
	long	err = 0;
	long	loss = 0;
	long	cross = 0;
	long	join = 0;
	long	right_area = 0;
	long	bkg = 0;
	double	fnr = 0.0;
	double	fpr = 0.0;
	double	ji = 0.0;
	FILE	*fp;

	fp = fopen( path.c_str(), "a");

	for( int i = 0; i < size; i++){
		if( input[i] == 1 || right[i] == 1){
			join++;
			if( input[i] == right[i])	cross++;
			else if( input[i] == 1)	err++;
			else	loss++;
			if( right[i] == 1)	right_area++;
		}
	}
	bkg = size - right_area;
	fpr = (double)err / (double)bkg;
	fnr = (double)loss / (double)right_area;
	ji = (double)cross / (double)join;
	fprintf( fp, "%f,%f,%f\n", ji, fpr, fnr);
	std::cout << "J.I. = " << ji << "\tfpr = " << fpr << "\tfnr = " << fnr << std::endl;
	fclose( fp);
	return ji;
}

//純粋にJ.I.を測定する
template< class T1, class T2>
double	JI( T1 *input, T2 *right, int size){

	long	cross = 0;
	long	join = 0;
	double	ji = 0.0;

	for( int i = 0; i < size; i++){
		if( input[i] == 1 || right[i] == 1){
			join++;
			if( input[i] == right[i])	cross++;
		}
	}
	ji = (double)cross / (double)join;
	return	ji;
}


//マスクがある場合のJ.I.,FPR,FNR測定
template< class T1, class T2>
void	comp_image_data( T1 *input, T2 *right, T2 *mask, int size, std::string &path){
	long	err = 0;
	long	loss = 0;
	long	cross = 0;
	long	join = 0;
	long	right_area = 0;
	long	bkg = 0;
	double	fnr = 0.0;
	double	fpr = 0.0;
	double	ji = 0.0;
	FILE	*fp;

	fp = fopen( path.c_str(), "a");

	for( int i = 0; i < size; i++){
		//if(input[i] != right[i] && input[i] == 1)	err++;
		//if(input[i] != right[i] && input[i] == 0)	loss++;
		//if( err > LONG_MAX)	exit(1);

		if( input[i] == 1 || right[i] == 1){
			join++;
			if( input[i] == right[i])	cross++;
			else if( input[i] == 1)	err++;
			else	loss++;
			if( right[i] == 1)	right_area++;
		}
		if( mask[i] == 1 && right[i] == 0)	bkg++;
	}
	fpr = (double)err / (double)bkg;
	fnr = (double)loss / (double)right_area;
	ji = (double)cross / (double)join;
	//fprintf( fp, "間違い数,%d\n", err);
	//fprintf( fp, "一致度,%f[%%]\n\n", (double)cross / (double)join * 100.0);
	fprintf( fp, "%f,%f,%f\n", ji, fpr, fnr);
	std::cout << "J.I. = " << ji << "\tfpr = " << fpr << "\tfnr = " << fnr << std::endl;
	fclose( fp);
}