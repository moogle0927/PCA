#include<list>
#include<string>
#include<Eigen/Core>
#include<Eigen/Eigenvalues>

namespace saito{

  /*=== Class model ========================= ここから =========================*/

  /**
  * PCAによるモデルを扱うためのクラス
  *
  * @author A. Saito
  * @version 1.0
  */
  template< class T >
  class model{
  public:
    typedef typename Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > Mat_t;
    typedef typename Eigen::Matrix< T, Eigen::Dynamic, 1 > Vec_t;
  protected:
    std::size_t num_of_dimX; /**< 入力空間の次元数 ( = D ) */
    std::size_t num_of_dimY; /**< 部分空間の次元数 ( = d ) */
    Vec_t X_mean;            /**< 平均ベクトル ( D x 1 ) */
    Mat_t M;                 /**< 写像行列 ( D x d ) */
    Vec_t Eval;              /**< 固有値 ( d x 1 ) */
    Vec_t sqrtL;             /**< 固有値の平方根 ( d x 1 ) */
    Vec_t invsqrtL;          /**< 固有値の平方根の逆数 ( d x 1 ) */
    Vec_t CCR;               /**< 累積寄与率 ( d x 1 ) */
    Vec_t w;                 /**< 重みベクトル ( d x 1 ) */
    bool isModelExists;      /**< モデルがロード済みであるかどうかを表すフラグ */
    bool isWeighted;      /**< モデルがロード済みであるかどうかを表すフラグ */

    void get_list_txt( std::list< std::string > &list_txt, std::string path_list );
    long get_file_size( std::string filename );
    void read_bin( T *p, std::size_t num, std::string file_name );
    void write_bin( const T *p, std::size_t num, std::string file_name );
	void write_text( const T *p, size_t num, std::string file_name );
    virtual void load_common( std::string modelDir, const std::size_t num_of_basis );
    void load_base( std::string modelDir, const std::size_t num_of_basis );

  public:
    model(){
      isModelExists = false;
      isWeighted = false;
    }
    /**
    * @brief PCAによるモデルの作成
    *
    * 入力データのリストからモデルを作成する．
    * @param[in] path_list 学習データリストのパス
    * @return PCAの実行に成功したかどうかを表すフラグ
    * @retval TRUE PCAの実行に成功
    * @retval FALSE PCAの実行中に失敗
    * @attention 入力データは全て同じ次元数でなければならない
    */
    void PCA( const std::string &path_list );
    /**
    * @brief PCAによるモデルの作成
    *
    * 入力データのリストからモデルを作成する．
    * @param[in] データ行列
    * @return PCAの実行に成功したかどうかを表すフラグ
    * @retval TRUE PCAの実行に成功
    * @retval FALSE PCAの実行中に失敗
    * @attention 入力データは全て同じ次元数でなければならない
    */
	void PCA( Mat_t X );
    /**
    * @brief WPCAによるモデルの作成
    *
    * 入力データのリストからモデルを作成する．
    * @param[in] path_list 学習データリストのパス
    * @param[in] alpha パラメータα (w = 1/(1+exp(α*d)))
    * @return WPCAの実行に成功したかどうかを表すフラグ
    * @retval TRUE PCAの実行に成功
    * @retval FALSE PCAの実行中に失敗
    * @attention 入力データは全て同じ次元数でなければならない
    */
    void WPCA( const std::string &path_list, T alpha );
    /**
    * @brief モデルのファイルへの出力
    *
    * モデルを指定したディレクトリに保存する
    * @param[in] outDir 出力ディレクトリ
    * @attention フォルダがない場合は自動的に作成するが，
    * 親ディレクトリの存在する場所のみ．
    */
    void output( std::string outDir );
    /**
    * @brief モデルの読み込み
    *
    * モデルを指定したディレクトリからロードする．
    * 基底数は全て．
    * @param[in] modelDir モデルのディレクトリ
    */
    void load( std::string modelDir );
    /**
    * @brief モデルの読み込み
    *
    * モデルを指定したディレクトリからロードする．
    * 基底数は指定した数だけ．
    * @param[in] modelDir モデルのディレクトリ
    * @param[in] num_of_basis 読み込む基底の数
    * @attention オーバーロード用にvirtual指定
    */
    void load_with_N( std::string modelDir, const std::size_t num_of_basis );
    /**
    * @brief モデルの読み込み
    *
    * モデルを指定したディレクトリからロードする．
    * 基底数は累積寄与率で指定．
    * @param[in] modelDir モデルのディレクトリ
    * @param[in] num_of_basis 読み込む基底の数
    * @attention 読み込まれる基底数は，累積寄与率が指定した値以上
    * となる最小の本数．
    */
    void load_with_CCR( std::string modelDir, const T ccr );
    /**
    * @brief 未知入力の部分空間への写像 (Out-of-sample)
    *
    * 未知入力を正規化なしで部分空間へ写像
    * @param[out] Y 写像後の主成分得点 (正規化なし)
    * @param[in]  X 写像する座標
    */
    void out_of_sample( Mat_t &Y, const Mat_t &X );
    /**
    * @brief 未知入力の部分空間への写像 (Out-of-sample)
    *
    * 未知入力を正規化して部分空間へ写像
    * @param[out] Alpha 写像後の主成分得点 (正規化あり)
    * @param[in]  X 写像する座標
    */
    void out_of_sample_normal( Mat_t &Alpha, const Mat_t &X );
    /**
    * @brief 部分空間から入力空間への逆写像 (pre_image)
    *
    * 部分空間の座標を入力空間の座標へ逆写像
    * @param[out] X 逆写像後の座標
    * @param[in] Y 逆写像する主成分得点 (正規化なし)
    */
    void pre_image( Mat_t &X, const Mat_t &Y );
    /**
    * @brief 正規化された部分空間から入力空間への逆写像 (pre_image)
    *
    * 正規化された部分空間の座標を入力空間の座標へ逆写像
    * @param[out] X 逆写像後の座標
    * @param[in] Alpha 逆写像する主成分得点 (正規化あり)
    */
    void pre_image_normal( Mat_t &Alpha, const Mat_t &Y );
    /**
    * @brief 写像・逆写像によるデータの再構成
    *
    * @param[out] Xout 写像・逆写像によって再構成されたデータ行列
    * @param[in] Xin 入力空間におけるデータ行列
    */
    void reconstruction( Mat_t &Xout, const Mat_t &Xin );
    /**
    * @brief モデル情報の表示
    * 
    * @attention デバッグ用．大規模なデータは表示しないこと．
    */
    virtual void disp( );
    virtual ~model(){}
  };

  /*=== Class model ========================= ここまで =========================*/


};
