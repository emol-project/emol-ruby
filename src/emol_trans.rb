#!/usr/bin/ruby
# -*- coding: utf-8 -*-

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative 'emol_molinfo'
require_relative 'emol_rhf'
require_relative 'emol_util'

require 'Inttrans'
require 'gsl'
#
#= Transformation 原子積分ー分子積分変換
#
#Authors::   tashi
#Version::   1.0 2015-12-01
#Copyright:: Copyright (C) Sapporo QC, 2015.
#License::   LGPL-3.0-or-later (see COPYING.LESSER and COPYING)
#
#= 要求される中間ファイル
# 1. Data_molinfo: 分子状態情報
# 2. Data_molint: 原子積分
# 3. Data_rhf: 分子軌道係数
#= 出力される中間ファイル
# Data_trnint
#= 主なメソッド
# 1. two_one_dimCopy : 原子積分(s, h, nz_g) と分子軌道係数を一次元配列に位置決めコピー
# 2. close_shell_energy : 変換後の分子積分から HF エネルギーを計算 (チェック用)
#== 呼び出される C ライブラリ(Inttrans:trnpp, trnsps, trnrr)
# 吉嶺博士等によって開発された Alchemy：のアルゴリズムに沿って C で作成
# 1. trnpp : 下三角に格納された対称行列を 2 インデックス変換
# 2. trnsps：転置
# 3. trnrr：

class Transformation
    attr_accessor :core_energy, :s_mo, :t_mo, :v_mo, :h_mo, :eris_mo

    def initialize( obj = nil )
        if obj.nil? then
            @core_energy = 0.0
            @s_mo = nil
            @t_mo = nil
            @v_mo = nil
            @h_mo = nil
            @eris_mo = nil
        end

        if ( obj.instance_of?( Array ) ) then     
            load obj[ 0 ]
            molinfo = Data_molinfo::molinfo_data

            load obj[ 1 ]
            s = Data_int::s_data
            t = Data_int::t_data
            v = Data_int::v_data
            h = Data_int::h_data
            nz_g = Data_int::nz_g_data

            load obj[ 2 ]
            cao = Data_rhf::cao_data

            nbf = molinfo.get_nbf
            npq = nbf * (nbf + 1) / 2
            n_frozen = molinfo.get_n_frozen
            nob = molinfo.get_nob

            dao = Density.new( cao, n_frozen, 2 ) 
            g = calc_gmat_nz_eris( nz_g, dao.get_body, nbf )
            f = h + g
            @core_energy = (dao.get_body * (h + f)).trace / 2.0 + molinfo.get_e_nucl
#
            print "\n \n ******T r a n s f o r m a t i o n \n \n"
            printf( "nbf, n_frozen, nob : %5d  %5d  %5d  \n", nbf, n_frozen, nob )
            printf( "frozen core energy : %14.8f \n", @core_energy) 

            cao_trans = GSL::Matrix.alloc( nbf, nob )
            for i in 0...nob
                for j in 0...nbf
                    cao_trans[ j, i ] = cao[ j, n_frozen + i ]
                end
            end

#           EmolUtil::print_2dim_ary(cao.to_a, "Original MO matrix", 8, "%12.6f")
            EmolUtil::print_2dim_ary(cao_trans.to_a, "Transformation matrix", 8, "%12.6f")

            s_work, t_work, v_work, h_work, eris_work, cao_1dim = two_one_dimCopy( nbf, nob, s, t, v, f, nz_g, cao_trans )
            cpq = GSL::Vector.alloc( nbf * nob )
            @s_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
            @t_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
            @v_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
            @h_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
#            Inttrans::trnpp( 1, 1, nbf, nob, s_work, s_mo, 1, cao_1dim, cpq )
#            Inttrans::trnpp( 1, 1, nbf, nob, h_work, h_mo, 1, cao_1dim, cpq )
            Inttrans::trnpp( 1, 1, nbf, nob, s_work, @s_mo, 1, cao_1dim, cpq )
            Inttrans::trnpp( 1, 1, nbf, nob, t_work, @t_mo, 1, cao_1dim, cpq )
            Inttrans::trnpp( 1, 1, nbf, nob, v_work, @v_mo, 1, cao_1dim, cpq )
            Inttrans::trnpp( 1, 1, nbf, nob, h_work, @h_mo, 1, cao_1dim, cpq )

#            EmolUtil::print_ary_tri(s_mo.to_a, "Overlap matrix", 10, "%12.6f")
#            EmolUtil::print_ary_tri(h_mo.to_a, "H matrix", 10, "%12.6f")
            EmolUtil::print_ary_tri(@s_mo.to_a, "Overlap            matrix", 10, "%12.6f")
            EmolUtil::print_ary_tri(@t_mo.to_a, "Kinetic            matrix", 10, "%12.6f")
            EmolUtil::print_ary_tri(@v_mo.to_a, "Nuclear Attraction matrix", 10, "%12.6f")
            EmolUtil::print_ary_tri(@h_mo.to_a, "H                  matrix", 10, "%12.6f")

            mi = 1
            mpq = nob * ( nob + 1 ) / 2
            for n in 1..npq
                Inttrans::trnpp(mi, mi, nbf, nob, eris_work, eris_work, 1, cao_1dim, cpq)
                mi += npq
            end
            Inttrans::trnsps(npq, mpq, eris_work)

            buf = GSL::Vector.alloc( npq )
            @eris_mo = GSL::Vector.alloc( mpq*(mpq + 1)/2 + 1 )
            m = 1; kin = 0
            for n in 1..mpq
                Inttrans::trnrr(m, n, nbf, nob, eris_work, buf, cao_1dim, cpq)
                m += npq
                for k in 0...n
                    @eris_mo[ kin + k ] = buf[ k ]
                end
                kin += n
            end
            EmolUtil::print_ary_tri(eris_mo.to_a, "Eris matrix", 10, "%12.6f")

            printf( "Total energy calculated from mo integrals : %14.8f \n", 
            @core_energy + close_shell_energy(h_mo, eris_mo, molinfo.get_nel - n_frozen * 2) )

        elsif ( obj.instance_of?( Transformation ) ) then
            @core_energy = obj.core_energy
            @s_mo = obj.s_mo.clone
            @t_mo = obj.t_mo.clone
            @v_mo = obj.v_mo.clone
            @h_mo = obj.h_mo.clone
            @eris_mo = obj.eris_mo.clone
        else
        end
    end 

    def ij( i_, j_ )
        i = i_; j = j_
        if i_ < j_
            i = j_; j = i_ 
        end
        i * ( i + 1 ) / 2 + j
    end 

    def two_one_dimCopy( nbf, nob, s, t, v, h, nz_g, cao )
        npq = nbf * (nbf + 1) / 2
        npqrs = npq * (npq + 1) / 2
        s_1dim = GSL::Vector.alloc( npq )
        t_1dim = GSL::Vector.alloc( npq )
        v_1dim = GSL::Vector.alloc( npq )
        h_1dim = GSL::Vector.alloc( npq )
        ij = -1
        for i in 0...nbf
            for j in 0..i
                ij += 1
                s_1dim[ ij ] = s[ i, j ]
                t_1dim[ ij ] = t[ i, j ]
                v_1dim[ ij ] = v[ i, j ]
                h_1dim[ ij ] = h[ i, j ]
            end
        end

        eris_1dim = GSL::Vector.alloc( npq * npq )
        nz_g.each do |x|
            i = x.get_i; j = x.get_j; k = x.get_k; l = x.get_l   
            pq = ij( i, j ); rs = ij( k, l )
            eris_1dim[ pq * npq + rs ] = x.get_value
            eris_1dim[ rs * npq + pq ] = x.get_value
        end

        cao_1dim = GSL::Vector.alloc( nob * nbf )
        ip = -1
        for i in 0...nob
            for p in 0...nbf
                ip += 1
                cao_1dim[ ip ] = cao[ p, i ] 
            end
        end

        return s_1dim, t_1dim, v_1dim, h_1dim, eris_1dim, cao_1dim
    end

    def close_shell_energy( h, eri, nel)
        e = 0.0
        n = nel / 2
        for i in 0...n
            e += 2 * h[ ij( i, i ) ]
        end
        for i in 0...n
            ii = ij( i, i )
            for j in 0...n
                jj = ij( j, j )
                iijj = ij( ii, jj )
                ji = ij( i, j )
                ijij = ij( ji, ji )
                e += 2 * eri[ iijj ] - eri[ ijij ]
            end
        end
        return e
    end
end


