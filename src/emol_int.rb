require_relative "Molint"
require_relative "emol_consts"
require_relative "emol_molinfo"
require "gsl"

def emol_int( molinfo )
  sint = Sint0.new
  tint = Tint0.new
  vint = Vint0.new
  eris = Eri0.new

  charge       = molinfo.get_charge
  numb_shell   = molinfo.get_numb_shell
  numb_atom    = molinfo.get_numb_atm
  numb_prim    = molinfo.get_numb_prim
  angular_type = molinfo.get_angular_type
  shel_lqn     = molinfo.get_shel_lqn
  shel_atm     = molinfo.get_shel_atm
  shel_tem     = molinfo.get_shel_tem
  shel_add     = molinfo.get_shel_add
  shel_ini     = molinfo.get_shel_ini
  atom_charg   = molinfo.get_atom_charg
  atom_xyz     = molinfo.get_atom_xyz
  prim_exp     = molinfo.get_prim_exp
  prim_coe     = molinfo.get_prim_coe
  len_lqn      = molinfo.get_len_lqn
  nbf          = molinfo.get_nbf
  thr_ovch     = molinfo.get_thr_ovch
  level_print  = molinfo.get_level_print

  # --- debug ---
  expected_ini = 0
  numb_shell.times do |ish|
    if shel_ini[ish] != expected_ini
      STDERR.printf("BUG: shel_ini[%d]=%d, expected=%d\n",
                    ish, shel_ini[ish], expected_ini)
    end
    expected_ini += len_lqn[shel_lqn[ish]]
  end
  STDERR.printf("DEBUG: angular=%d nbf=%d expected_nbf=%d\n",
                molinfo.get_angular_type, nbf, expected_ini)
  # --- end debug ---



  sint.init( numb_atom, numb_shell, numb_prim,
             shel_lqn, shel_tem, shel_atm, shel_add,
             atom_charg, atom_xyz, prim_exp, prim_coe,
             thr_ovch, level_print, angular_type )

  tint.init( numb_atom, numb_shell, numb_prim,
             shel_lqn, shel_tem, shel_atm, shel_add,
             atom_charg, atom_xyz, prim_exp, prim_coe,
             thr_ovch, level_print, angular_type )

  vint.init( numb_atom, numb_shell, numb_prim,
             shel_lqn, shel_tem, shel_atm, shel_add,
             atom_charg, atom_xyz, prim_exp, prim_coe,
             thr_ovch, level_print, angular_type )

  eris.init( numb_atom, numb_shell, numb_prim,
             shel_lqn, shel_tem, shel_atm, shel_add,
             atom_charg, atom_xyz, prim_exp, prim_coe,
             thr_ovch, level_print, angular_type )

  s, t, v, h = calc_hmat( numb_shell, shel_lqn, shel_atm, shel_ini, sint, tint, vint, len_lqn, nbf )
  nz_g = mk_eris( numb_shell, shel_lqn, shel_atm, shel_ini, eris, len_lqn, nbf )
  return s, t, v, h, nz_g
end

def calc_hmat ( numb_shell, shel_lqn, shel_atm, shel_ini, sint, tint, vint, len_lqn, nbf )

  s = GSL::Matrix.alloc(nbf, nbf)
  t = GSL::Matrix.alloc(nbf, nbf)
  v = GSL::Matrix.alloc(nbf, nbf)
  h = GSL::Matrix.alloc(nbf, nbf)

  s_int     = [ ]
  t_int     = [ ]
  v_int     = [ ]
  for ish in 0...numb_shell
    for jsh in 0..ish

      inttype, nsize_int, s_int = sint.calc( ish, jsh )
      inttype, nsize_int, t_int = tint.calc( ish, jsh )
      inttype, nsize_int, v_int = vint.calc( ish, jsh )
#           printf "Orb: ( %d %d ): type = %d nsize = %d ncomps = %d %d atoms = %d %d\n", ish, jsh, inttype, nsize_int,
#                  len_lqn[ shel_lqn[ ish ] ], len_lqn[ shel_lqn[ jsh ] ],
#                  shel_atm[ ish ], shel_atm[ jsh ]
#           for i in 0...nsize_int
#                printf " %5d  %23.15e %23.15e %23.15e\n",
#               printf " %5d  %23.15e %23.15e\n",
#                      i, s_int[ i ], t_int[ i ], v_int[ i ]
#           end

      ncomp_i = len_lqn[ shel_lqn[ ish ] ]
      ncomp_j = len_lqn[ shel_lqn[ jsh ] ]
      for ip in 0..ncomp_i - 1 do
        if ( ish ==  jsh )
          last_j = ip + 1
        else
          last_j = ncomp_j
        end
        for jp in 0..last_j - 1 do
#                nn = ip * ncomp_j + jp
          nn = jp * ncomp_i + ip
          i = shel_ini [ ish ] + ip
          j = shel_ini [ jsh ] + jp
          s[ i, j ] = s_int[ nn ]
          t[ i, j ] = t_int[ nn ]
          v[ i, j ] = v_int[ nn ]
          h[ i, j ] = t_int[ nn ] + v_int[ nn ]
          if i != j
            s[ j, i ] = s[ i, j ]
            h[ j, i ] = h[ i, j ]
          end
        end
      end 

    end
  end
  return s, t, v, h
end

def mk_eris ( numb_shell, shel_lqn, shel_atm, shel_ini, eris, len_lqn, nbf )
  nz_eris = []

  for ish in 0...numb_shell
    for jsh in 0..ish
      for ksh in 0..ish
        if ksh < ish then
          l_last = ksh
        else
          l_last = jsh
        end
        for lsh in 0..l_last
          inttype, nsize_int, g_int = eris.calc( ish, jsh, ksh, lsh )
          ash, bsh, csh, dsh        = eris.exchange_order

          #STDERR.printf "Orb: ( %d %d %d %d ): type = %d nsize = %d ncomps = %d %d %d %d atoms = %d %d %d %d\n",
          #    ash, bsh, csh, dsh, inttype, nsize_int,
          #    len_lqn[ shel_lqn[ ash ] ], len_lqn[ shel_lqn[ bsh ] ],
          #    len_lqn[ shel_lqn[ csh ] ], len_lqn[ shel_lqn[ dsh ] ],
          #    shel_atm[ ash ], shel_atm[ bsh ], shel_atm[ csh ], shel_atm[ dsh ]
          #for i in 0..nsize_int-1
          #    STDERR.printf " %5d  %23.15e\n", i, g_int[ i ]
          #end

          thr_zero = EmolConsts.const_get(:THR_ZERO)
          ncomp_i = len_lqn[ shel_lqn[ ash ] ]
          ncomp_j = len_lqn[ shel_lqn[ bsh ] ]
          ncomp_k = len_lqn[ shel_lqn[ csh ] ]
          ncomp_l = len_lqn[ shel_lqn[ dsh ] ]
          flg_ieqj = ( ash == bsh )
          flg_keql = ( csh == dsh )
          flg_same = ( ash == csh ) && ( bsh == dsh )

          ijn = -1
          for ip in 0..ncomp_i - 1 do
            if ( flg_ieqj )
              last_j = ip + 1
            else
              last_j = ncomp_j
            end
            for jp in 0..last_j - 1 do
              ijn += 1; kln = -1
              for kp in 0..ncomp_k - 1 do
                if ( flg_keql )
                  last_l = kp + 1
                else
                  last_l = ncomp_l
                end
                for lp in 0..last_l - 1 do
                  kln += 1
                  if flg_same && kln > ijn then
                    next
                  end

                  nn = ( ( ip * ncomp_j + jp ) * ncomp_k + kp ) * ncomp_l + lp
                  x = g_int [ nn ]
                  if x.abs < thr_zero then
                    next
                  end
                  i = shel_ini [ ash ] + ip
                  j = shel_ini [ bsh ] + jp
                  k = shel_ini [ csh ] + kp
                  l = shel_ini [ dsh ] + lp
#
                  nz_eris.push( NZ_eris.new( i, j, k, l, x ) )
                end  
              end
            end
          end 
        end
      end
    end
  end
  nz_eris
end

#############################
#############################

class NZ_eris
  def initialize(i, j, k, l, value)
    @i = i; @j = j; @k = k; @l = l
    @value = value
  end

  def get_i
    @i
  end

  def get_j
    @j
  end

  def get_k
    @k
  end

  def get_l
    @l
  end

  def get_value
    @value
  end

  def get_index
    [@i, @j, @k, @l]
  end
end

