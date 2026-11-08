#!/usr/bin/ruby
# -*- coding: utf-8 -*-

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative './emol_expression' ## check order !!!
require_relative './emol_trans'
require_relative './emol_liu'
require_relative './emol_consts'
require_relative './emol_futil'

#require 'profile'
require 'gsl'
#
#= H_matrix 生成と対角化するクラス
#
#Authors::   tashi
#Version::   1.0 2015-12-01
#Copyright:: Copyright (C) Sapporo QC, 2015.
#License::   LGPL-3.0-or-later (see COPYING.LESSER and COPYING)
#
#= 要求される中間ファイル
# 1. Data_molinfo: 分子状態情報
# 2. Data_trnint: 分子積分
# 3. Data_expression: エネルギー表式 
#= 出力される中間ファイル
# Data_hmatrix
#= 主なメソッド
# 1. mkhmat : H 行列をエネルギー表式、分子積分から計算し、対角要素(h_diag)、非対角要素(index_i, index_j, value) を生成する。
# 2. get_hc : 複数の CI 係数からなる2次元配列から H 行列との積を計算し、結果を 2 次元配列として返す。Liu などの対角化メソッドから call される。
#== 呼び出されるクラス
# 1. Liu (大規模実対称行列対角化)
class H_matrix
#    csf の個数
    attr_accessor :n
#    molinfo.get_nstate 個の最低固有値 (1 次元配列)
    attr_accessor :eigen_value
#    molinfo.get_nstate 個の最低固有ベクトル(2 次元配列)
    attr_accessor :eigen_vector
#    H 行列の対角要素 (1 次元配列)
    attr_accessor :h_diag
#    H 行列の非零非対角要素の第一, 第二添字 (1 次元配列)
    attr_accessor :index_i, :index_j
#    H 行列の(index_i, index_j)非零非対角要素の値 (1 次元配列)
    attr_accessor :value

    def initialize( obj = nil )
        if obj.nil? then
            @n = 0
            @h_diag = []
            @index_i = []
            @index_j = []
            @value = []
            @eigen_value = []
            @eigen_vector = []
#            @eigen_value = GSL::Vector[]
#            @eigen_vector = GSL::Matrix[]
        end

        if ( obj.instance_of?( Array ) ) then   #
#           molint = load_trnint( obj[ 0 ] )
#           molinfo, s_mo, h_mo, eris_mo = load_trnint( obj[ 0 ] )
            load obj[ 0 ]
            molinfo = Data_molinfo::molinfo_data

            load obj[ 1 ]
#           load MO integrals
            molint = Data_trnint::trnint_data
            core_energy = molint.core_energy
#           one electron integrals (tri)
            h_mo = molint.h_mo
#           two electron integrals (tri, tri, tri)
            eris_mo = molint.eris_mo

#           expr   = load_expression( obj[ 1 ] )
            load obj[ 2 ]
            expr = Data_expression::expr_data
#           expr = Data_super_CI_expression::expr_data
            @n = expr.csf_list.size
#                   H 行列の生成
            mkhmat( expr, h_mo, eris_mo )
            nstate = molinfo.get_nstate
#           thresh = EmolConsts::THR_CONV_WF
	    thresh = molinfo.get_thresh_ci
#                   H 行列の Liu によるアルゴリズムによる対角化
            eigen_val, @eigen_vector = Liu.new( self, nstate, thresh ).solve
#                   凍結殼エネルギーの付加
            @eigen_value = GSL::Vector.alloc( nstate )
            (0...nstate).each do | k |
                @eigen_value[ k ] = eigen_val[ k ] + molint.core_energy
            end
            printf( "\n\n" )
            printf( " Eigenvalues\n" )
            (0...nstate).each do | k |
                printf( "  %5d          ", k + 1 )
            end
            printf( "\n" )

            (0...nstate).each do | k |
                printf( "  %14.10f ", @eigen_value[ k ] )
            end
            printf( "\n\n" )
            
        elsif ( obj.instance_of?( Expression ) ) then
            @ncsf = obj.ncsf
            @h_diag = obj.h_diag.clone
            @index_i = obj.index_i.clone
            @index_j = obj.index_j.clone 
            @value = obj.value.clone 
            @eigen_value = obj.eigen_value.clone
            @eigen_vector = obj.eigen_vector.clone
        else
        end
    end

    def mkhmat( expr, h_mo, eris_mo )
# 機能 ：H 行列要素の生成
#     : 対角要素 ... @h_diag ( 1次元配列)
#     : 非対角要素 ... 非零要素のみ(閾値 EmolConsts::THR_ZERO 以上)
#     :               インデックス(index_i, index_j) と値(value)
#     :                 
# パラメータ：expr ... エネルギー表式インスタンス
#         : h_mo ... 1電子積分
#         ：eris_mo ... 2電子積分
# 戻り値   ：無し
# 備考：エネルギー表式インスタンスは、CSF の組み合わせ ij、 積分通し番号 pqrs、
#      結合係数 coef の配列をクラス変数として持つ。一つの H 要素は複数の積分からの
#      寄与があり、インスタンス生成時に ij に関しては整列化されている。
#      このインスタンスに含まれる H 要素には分子積分の値から非常に絶対値の小さなもの
#      もあり、閾値 EmolConsts::THR_ZERO 以下のものはゼロとみなす。
#                    
# 初期設定
        @h_diag = []
        @index_i = []
        @index_j = []
        @value = []
        i_prev = 1; j_prev = 1; ij_prev = 1; ij_end = 0
        tmp = 0.0
# エネルギー表式ループ開始
        expr.ij.each_with_index do |ij, k|
            num = expr.pqrs[ k ]
            coef = expr.coef[ k ]
# 同じ H 要素の処理を続行
            if ij == ij_prev
                i = i_prev; j = j_prev 
# 処理中のH 要素の書き出し
            else
                if tmp.abs > EmolConsts::THR_ZERO 
                    if i_prev == j_prev
                        @h_diag.push( tmp )
                    else
                        @index_i.push( i_prev )
                        @index_j.push( j_prev )
                        @value.push( tmp )
                    end
                end
# 前表式と異なる H 要素への処理開始、index_i, index_j の生成
                tmp = 0.0
                if ij <= ij_end + i_prev
                    i = i_prev
                    j = ij - ij_end
                else
                    ij_end += i_prev
                    i = i_prev + 1
                    j = ij - ij_end
                end
            end
# 係数*積分値の足し込み
            if num >= h_mo.size
                tmp += coef * eris_mo[ num - h_mo.size ]  
            else
                tmp += coef * h_mo[ num ]  
            end
            ij_prev = ij
            i_prev = i
            j_prev = j
        end
# H 対角要素の書き出し
        if tmp.abs > EmolConsts::THR_ZERO 
            if i_prev == j_prev
                @h_diag.push( tmp )
            else
                @index_i.push( i_prev )
                @index_j.push( j_prev )
                @value.push( tmp )
            end
        end
#        @h_diag.each_with_index do | v, i |
#           printf( " %5d   %12.6f  \n", i, v )
#        end
    end

    def init_vec( _nvec )
        diag_sort = @h_diag.sort
        b = GSL::Matrix.alloc( @n, _nvec )
        for i in 0..._nvec
            index = @h_diag.index( diag_sort[ i ] )
            b[ index, i ] = 1.0
            if index == i 
                b[ mod(i + 1, @n), i ] = 0.00001
            else
                b[ i, i ] = 0.00001
            end
        end
        b
    end

    def get_n
        @n
    end

    def get_diag( i )
#       H 行列の i-番目の対角要素の取得
#       arg[0] ... i
#       戻り値 ... i-番目の対角要素の取得 
        @h_diag[ i ]
    end

    def get_hc( c )
        hc = GSL::Matrix.alloc( c.size1, c.size2 )

        @h_diag.each_with_index do | v, i |
            for p in 0...c.size2
               hc[ i, p ] += v * c[ i, p ]
            end
        end

        @index_i.each_with_index do | i, k |
            j = @index_j[ k ] 
            v = @value[ k ]
            for p in 0...c.size2
                hc[ i - 1, p ] += v * c[ j - 1, p ]
                hc[ j - 1, p ] += v * c[ i - 1, p ]
            end
        end
        hc
    end
#   def print_ci_result( istate, eval, evec, csf_list, n_core, n_active, n_external )
#       if istate == 1
#           printf("\n %2dst CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
#       elsif istate == 2
#           printf("\n %2dnd CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
#       elsif istate == 3
#           printf("\n %2drd CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
#       else
#           printf("\n %2dth CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
#       end
#
#       printf("    CSF       COEF     ( CORE / ACTIVE / EXTERNAL )\n")
#       printf("    ---    ----------   ") 
#       n_core.times do | i |
#           printf("--")
#       end
#       printf("  ")
#       n_active.times do | i |
#           printf("--")
#       end
#       printf("  ")
#       n_external.times do | i |
#           printf("--")
#       end
#       printf("  \n")
#
#       vec_print = {}
#       vec = evec.column( istate - 1 )
#       vec.to_a.each_with_index do |x, i|
#           if x.abs > 0.01
#               vec_print.store( i,  x )
#           end
#       end
#
#       vec_print.sort{ |(k1, v1), (k2, v2)| v2.abs <=> v1.abs }.each do | x |
#           printf("  %5d  ", x[ 0 ] + 1)
#           printf(" %10.6f   ", x[ 1 ])
#           a = csf_list[ x[ 0 ] ]
#
#           k = 0
#           n_core.times do | i |
#               printf("%2s", a[ k ] )
#               k += 1
#           end
#           printf("  ")
#           n_active.times do | i |
#               printf("%2s", a[ k ] )
#               k += 1
#           end
#           printf("  ")
#           n_external.times do | i |
#               printf("%2s", a[ k ] )
#               k += 1
#           end
#           printf("  \n")
#       #   csf_list[ x[ 0 ] ].each do |y|
#       #       printf("%2s", y)
#       #   end
#       end
#   end
end

