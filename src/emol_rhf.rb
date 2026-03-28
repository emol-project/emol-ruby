#!/usr/bin/ruby
# coding: utf-8

require_relative 'emol_rexml'
require_relative 'emol_util'
require_relative 'emol_consts'
require_relative 'emol_int'
require_relative 'emol_molinfo'
require 'gsl'

# = GSL::Matrix に trace を計算する関数 trace を追加
#
class GSL::Matrix
  def trace
    sum = 0.0
    for i in 0..self.size1-1
      sum += self[i,i]
    end   
    sum
  end
end

# = Array を継承した LinkArray クラス，予めサイズを決めた Array.
# データを追加するために link_push を利用した場合には，予め設定したサイズを越えないようになる
#
class LinkArray < Array
  # 
  def initialize( n )
    @size = n
    []
  end

  # データを末尾に追加する際に，予め決めておいたサイズを越えた場合には
  # 既に保存されたデータを先頭から取り出し，サイズを越えないようにする
  def link_push( a )
    self.push( a )
    if self.size > @size
      self.shift
    end
  end
end

# = DIIS 計算のためのクラス
#
class DIIS
  def initialize
    n = EmolConsts::MAX_DIIS
    @err = LinkArray.new( n )
    @fock = LinkArray.new( n )
  end

  def err_push( a )
    @err.link_push( a )
  end

  def fock_push( a )
    @fock.link_push( a )
  end

  def get_averageFock
    n = @fock.size
    if n < EmolConsts::MIN_DIIS
      return @fock.last
    else
      a = GSL::Matrix.alloc( n+1, n+1 )
      b = GSL::Vector.alloc( n+1 )
      for i in 0..n-1
        for j in 0..i
          a[i,j] = a[j,i] = @err[i] * @err[j].col
        end
        a[i,n] = a[n,i] = -1.0
        b[i] = 0.0
      end
      a[n,n] = 0.0; b[n] = -1.0

#     EmolUtil::print_2dim_ary(a.to_a, "Matrix A", 10, "%12.6f")
#     Emolutil::print_ary(b.to_a, "Matrix B", 10, "%12.6f")

      # 連立方程式 a*x = b を解く
      w = GSL::Linalg::SV.solve( a, b )
      m = @fock[0].size1
      aveFock = @fock[0].scale( w[0] )
      for i in 1..@fock.size-1
        aveFock += @fock[i].scale( w[i] )
      end
    end 
    return aveFock
  end
end

#
#= 2電子積分と密度行列から，G行列を計算
#
#== 入力
# nz_eris:  2電子積分
# d:        密度行列
# nbf:      基底関数数
#
#== 出力
# g:        G行列
#
#== アルゴリズム
# Direct-SCF 計算において i,j,k,l の AO の積分 x が得られている場合にその要素の G 行列への寄与を計算するアルゴリズムを，
# 積分を全て nz_eris に保持している場合に適用する．(nbf * nbf) の全ての要素を計算する．
#
def calc_gmat_nz_eris ( nz_eris, d, nbf )
  g = GSL::Matrix.alloc(nbf, nbf)

  nz_eris.each do |xx|
    i = xx.get_i; j = xx.get_j; k = xx.get_k; l = xx.get_l; 
    x = xx.get_value
    #p xx

    if i == j then; x *= 0.5; end 
    if k == l then; x *= 0.5; end 
    if i == k && j == l then; x *= 0.5; end

    g [ i, j ] += 4 * x * d [ k, l ]
    g [ k, l ] += 4 * x * d [ i, j ]
    g [ i, k ] -= x * d [ j, l ]
    g [ i, l ] -= x * d [ j, k ]
    g [ j, k ] -= x * d [ i, l ]
    g [ j, l ] -= x * d [ i, k ]
#    printf( "%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e : %12.4e%12.4e%12.4e%12.4e%12.4e%12.4e : %12.4e\n",
#    g [ i, j ], g [ k, l ], g [ i, k ], g [ i, l ], g [ j, k ], g[ j, l ],
#    d [ k, l ], d [ i, j ], d [ j, l ], d [ j, k ], d [ i, l ], d[ i, k ], x )
  end

  for i in 0..nbf-1 do
    for j in 0..i-1 do
      tmp = g [ j, i ] + g [ i, j ] 
      g [ j, i ] = g [ i, j ] = tmp / 2.0
    end
  end
#  p g
  g 
end

#
#= 密度行列クラス
#
#== クラスメンバデータ
# n:        密度行列の次元数
# d:        密度行列
# g:        G行列
#
#== アルゴリズム
# Direct-SCF 計算において i,j,k,l の AO の積分 x が得られている場合にその要素の G 行列への寄与を計算するアルゴリズムを，
# 積分を全て nz_eris に保持している場合に適用する．(nbf * nbf) の全ての要素を計算する．
#
class Density
  #== 密度行列の初期化
  #=== 入力
  # c:          MO 係数
  # nocc:       占有軌道
  # occupation: 占有軌道の占有数
  def initialize( c, nocc, occupation )
    @n = c.size1
    @d = GSL::Matrix.alloc( @n, @n )
    (0..@n-1).each do |p|    
      (0..p).each do |q|
        sum = 0.0
        (0..nocc-1).each do |i|
          sum += c[p,i] * c[q,i]
        end
        @d[p,q] = @d[q,p] = sum * occupation
      end
    end  
  end

  #=== 密度行列へのアクセス関数
  # c:          MO 係数
  # nocc:       占有軌道
  # occupation: 占有軌道の占有数
  def get_body
    @d
  end

  def check_idenpotency( s )
        EmolUtil::print_2dim_ary( (@d*s*@d - @d*2.0).to_a, "Check idenpotency of density : PSP - 2P = 0", 10, "%12.6f" )
  end

  def show( os )
    for p in 0...@n do
      for q in 0..p do
        os.printf( "%6d%6d%28.16e\n", p, q, @d[ p, q ] )
      end
    end
  end
end

#
#= RHF 計算関数
#== 入力
# molinfo: MO 係数
# s:       重なり積分
# h:       1電子フォック行列
# nz_g:    全ての2電子積分
def rhf( molinfo, s, h, nz_g )

    nbf = molinfo.get_nbf
    nel = molinfo.get_nel
    #
    # == 正準直交化 (Sazbo, Sec.3.4.5)
    #
    #   Canonical orbital transformation matrix : tm
    #
    # ===重なり積分の対角化
    eval, evec = GSL::Eigen::symmv( s )
    #
    # ===値の大きい順にソート
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_DESC )
       #EmolUtil::print_2dim_ary_tri(s.to_a, "Overlap matrix", 10, "%12.6f")
       #EmolUtil::print_2dim_ary_tri(h.to_a, "Hcore matrix", 10, "%12.6f")
    
    # === EmolConsts::THR_THROWN より小さな値を持つ軌道の数を調べる
    num_thrown = 0; thrown = EmolConsts::THR_THROWN
    eval.each do | x |
      if x < thrown
        num_thrown += 1
      end
    end
    # ===対角化した行列を必要分のみのサイズ: eval.size - num_thrown に変更
    eval = eval.subvector( eval.size - num_thrown )
    evec = evec.submatrix( nil, 0, evec.size[0] - num_thrown )
    for i in 0..eval.size-1
      eval[i] = 1.0 / Math.sqrt( eval[i] )
    end
    # === X = U * s^{-1/2} の計算
    #   EmolUtil::print_2dim_ary_with_value(evec.to_a, eval.to_a, "Rectangular martrix", 10, "%12.6f")
    tm = evec * GSL::Matrix.diagonal( eval )
    #                                                                      X = U s^{-1/2}   Eq(3.172)
    #   EmolUtil::print_2dim_ary((tm.transpose * s * tm).to_a, "Check X(T)SX = 1", 10, "%12.6f")

    # === (X+)*h*X の計算 (Sazbo, Eq.3.172)
    #   make initial guess for mo from core Hamiltonian h
    eval, evec = GSL::Eigen::symmv( tm.transpose * h * tm )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )
    cao = tm * evec
       #EmolUtil::print_2dim_ary_with_value_moname(cao.to_a, eval.to_a, molinfo.get_bs_name, "Orbital Energies and MO Coefficients", 8, "%12.6f")

    # === 密度行列の生成
    dao = Density.new( cao, nel/2, 2.0 )

#
#   scf iteration
#
    diis = DIIS.new
    error = GSL::Vector.alloc( [100.0, 100.0] )
    iter = -1; g = nil
    printf( "\n\nSCF iteration:\n" )
    printf( "    Iter:     Energy                     Error \n" )
#   while error * error.col > EmolConsts::THR_CONV
    while error * error.col > molinfo.get_thresh_hf
      # == 2電子積分と密度行列から2電子フォック行列を計算
      g = calc_gmat_nz_eris( nz_g, dao.get_body, nbf )
      # == 1電子，2電子フォック行列の和
      f = h + g
      iter += 1
      # == 全エネルギー計算 (Szabo, Eq.3.184)
      total_energy = (dao.get_body * (h + f)).trace / 2.0 + molinfo.get_e_nucl
      if iter > 0 then
        printf( "   %3d     %20.16e     %20.16e \n", iter, total_energy, error*error.col )
      else
        printf( "   %3d     %20.16e  \n",            iter, total_energy )
      end

      # == DIIS 計算
      # === error 行列の計算 (Pulay, JCC_2_1982, Eq.2)
      fmo = cao.transpose * f * cao
      err = []
      nocc = nel / 2
      for i in 0..nocc-1
        for j in nocc..nbf-1
          err.push( fmo[i, j] )
        end
      end
      error = GSL::Vector.alloc( err )
      diis.err_push( error )
      diis.fock_push( f )
      #
      #     printf("dmat, gmat, fmat: \n")
      #     ij = -1
      #     for i in 0..g.size1-1
      #       for j in 0..i
      #         ij += 1
      #         printf(" %5d  %5d  %6d  %20.15e  %20.15e  %20.15e  %20.15e \n",
      #                j, i, ij, dao.get_body[j,i], h[j,i], g[j,i], f[j,i])
      #       end
      #     end
      #

      eval, evec = GSL::Eigen::symmv( tm.transpose * diis.get_averageFock * tm )
      GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )
      eval_vector = GSL::Vector.alloc( nbf )
      for i in 0...nbf do
        eval_vector[ i ] = eval[ i ]
      end
      cao = tm * evec
      dao = Density.new( cao, nel/2, 2.0 )
  end
  #
  # 収束軌道にてフォック行列を計算する
  g = calc_gmat_nz_eris( nz_g, dao.get_body, nbf )
  f = h + g
  # 最終全エネルギー計算
  elec_energy  = (dao.get_body * (h + f)).trace / 2.0
  total_energy = elec_energy + molinfo.get_e_nucl

  printf( "\n Nuclear repulsion energy : %28.16e\n", molinfo.get_e_nucl )
  printf( " Electronic energy        : %28.16e\n", elec_energy )
  printf( " Total HF energy          : %28.16e\n", total_energy )

#
  EmolUtil::print_2dim_ary_with_value_moname(cao.submatrix(nil, 0...nocc).to_a, eval.to_a, molinfo.get_bs_name, 
	     "Occupied Orbital Energies and MO Coefficients", 8, "%12.6f")
# EmolUtil::print_2dim_ary_with_value_moname(cao.to_a, eval.to_a, molinfo.get_bs_name, 
#            "Orbital Energies and MO Coefficients", 8, "%12.6f")

  # 全エネルギー，フォック行列，電子フォック行列，係数行列，密度行列データ，Xs^{-1/2}
  return elec_energy, total_energy, eval_vector, f, g, cao, dao.get_body, tm
end
