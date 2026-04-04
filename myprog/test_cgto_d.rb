#!/usr/bin/env ruby
# test_cgto_d.rb
# ==============================================================
# Equivalent to test_ref2_cgto_d.py
#
# System: He(0,0,0) -- H(1.0 Ang, 0, 0)
#   He: contracted d-shell, NPRIM=3
#       alpha = [5.0, 1.5, 0.5], raw_coef = [0.2, 0.5, 0.4]
#   H:  s-shell, NPRIM=1
#       alpha = [1.0],           raw_coef = [1.0]
#
# emol_coe = raw * gto_norm(l, alpha) * sqrt((2l+1)/(4*pi))
#   d (l=2): sqrt(5/(4pi)) = 0.6307831305
#     alpha=5.0: 0.2 * 43.6241415262 * 0.6307831305 = 5.5034745115
#     alpha=1.5: 0.5 *  5.3050452105 * 0.6307831305 = 1.6731665127
#     alpha=0.5: 0.4 *  0.7757591266 * 0.6307831305 = 0.1957343081
#   s (l=0): sqrt(1/(4pi)) = 0.2820947918
#     alpha=1.0: 1.0 *  2.5264751110 * 0.2820947918 = 0.7127054704
#
# Reference values from PySCF (normalize_GTO=True, cart=True):
#   REFS   = 0.7206973446   S(dxx, s)
#   REFT   = 0.6570554634   T(dxx, s)
#   REFVHE = 0.5031533913   +<dxx|1/r_He|s>  (PySCF int1e_rinv)
#   REFVH  = 0.9398107814   +<dxx|1/r_H |s>  (PySCF int1e_rinv)
#   REFERI = 1.5241713338   ERI(dxx dxx|ss)
# ==============================================================

$LOAD_PATH.unshift('.')
require 'Molint'
#require_relative 'Molint'

# --------------------------------------------------------------
# Reference values from PySCF
# --------------------------------------------------------------
REFS   =  0.7206973446
REFT   =  0.6570554634
REFVHE =  0.5031533913   # +<dxx|1/r_He|s>
REFVH  =  0.9398107814   # +<dxx|1/r_H |s>
REFERI =  1.5241713338

# GINT V = sum_A (-Z_A) * <i|1/r_A|j>  = -2*REFVHE - 1*REFVH
REFV = -2.0 * REFVHE - 1.0 * REFVH

# --------------------------------------------------------------
# Molecular data
# --------------------------------------------------------------
ANG2BOHR = 1.8897259886

numb_atom  = 2
numb_shell = 2
numb_prim  = 4

shel_lqn = [2, 0]
shel_tem = [3, 1]
shel_atm = [0, 1]
shel_add = [0, 3]

atom_charg = [2.0, 1.0]
atom_xyz   = [0.0, 0.0, 0.0,
              1.0 * ANG2BOHR, 0.0, 0.0]

prim_exp = [5.0, 1.5, 0.5,
            1.0]

# emol規格化済み係数  (raw * gto_norm * sqrt((2l+1)/(4pi)))
prim_coe = [5.5034745115, 1.6731665127, 0.1957343081,   # d on He
            0.7127054704]                                 # s on H

tol          = 1.0e-10
level_dbg    = 0
angular_type = 0   # Cartesian

ISH_D   = 0
ISH_S   = 1
IDX_DXX = 0   # GAMESS d order: xx,yy,zz,xy,xz,yz -> dxx=0

# --------------------------------------------------------------
# Helper
# --------------------------------------------------------------
def check(label, got, ref, tol = 1.0e-5)
  diff = (got - ref).abs
  ok   = diff < tol ? "OK" : "*** FAIL ***"
  printf("  %-24s got=%+.10f  ref=%+.10f  diff=%.2e  %s\n",
         label, got, ref, diff, ok)
end

# --------------------------------------------------------------
# Init
# --------------------------------------------------------------
args = [numb_atom, numb_shell, numb_prim,
        shel_lqn, shel_tem, shel_atm, shel_add,
        atom_charg, atom_xyz,
        prim_exp, prim_coe,
        tol, level_dbg, angular_type]

sint = Sint0.new; sint.init(*args)
tint = Tint0.new; tint.init(*args)
vint = Vint0.new; vint.init(*args)
eri0 = Eri0.new;  eri0.init(*args)

# --------------------------------------------------------------
# 1E: d(He) -- s(H)
# --------------------------------------------------------------
puts "=== 1E integrals: d(He) -- s(H) ==="

_, _, s_arr = sint.calc(ISH_D, ISH_S)
_, _, t_arr = tint.calc(ISH_D, ISH_S)
_, _, v_arr = vint.calc(ISH_D, ISH_S)

puts "S array (GAMESS order: xx,yy,zz,xy,xz,yz):"
puts "  #{s_arr.map{|x| x.round(8)}}"
puts "T array:"
puts "  #{t_arr.map{|x| x.round(8)}}"
puts "V array:"
puts "  #{v_arr.map{|x| x.round(8)}}"

check("S(dxx,s)",  s_arr[IDX_DXX], REFS)
check("T(dxx,s)",  t_arr[IDX_DXX], REFT)
check("V(dxx,s)",  v_arr[IDX_DXX], REFV)

# --------------------------------------------------------------
# V per nucleus (符号・スケール確認用)
# --------------------------------------------------------------
puts "\n--- V per nucleus (calc_raw) ---"
_, vr_n, vr_arr = vint.calc_raw(ISH_D, ISH_S)
puts "raw V array (#{vr_n} elements):"
puts "  #{vr_arr.map{|x| x.round(8)}}"

# layout: buf[k * nci*ncj + jp*nci + ip]  nci=6, ncj=1
vhe_raw = vr_arr[0 * 6 + IDX_DXX]   # k=0: He
vh_raw  = vr_arr[1 * 6 + IDX_DXX]   # k=1: H
printf("  V_He_raw=%+.10f  PySCF_rinv_He=%+.10f  ratio=%.6f\n",
       vhe_raw, REFVHE, vhe_raw / REFVHE)
printf("  V_H_raw =%+.10f  PySCF_rinv_H =%+.10f  ratio=%.6f\n",
       vh_raw,  REFVH,  vh_raw  / REFVH)
# 期待: GINT raw = -Z * PySCF_rinv なら ratio = -Z
#       He: ratio = -2.0,  H: ratio = -1.0

# --------------------------------------------------------------
# 2E: (dxx dxx | ss)
# --------------------------------------------------------------
puts "\n=== 2E integrals: (d d | s s) ==="

e_type, e_n, e_arr = eri0.calc(ISH_D, ISH_D, ISH_S, ISH_S)
order = eri0.exchange_order
puts "inttype=#{e_type}  nsize=#{e_n}  order=#{order.inspect}"
puts "ERI array (first 6): #{e_arr[0,6].map{|x| x.round(8)}}"

# emol layout: ((ip*ncb+jp)*ncc+kp)*ncd+lp
# nca=ncb=6(d), ncc=ncd=1(s)  -> (dxx dxx|ss): ip=0,jp=0,kp=0,lp=0 -> index 0
check("ERI(dxx,dxx|ss)", e_arr[0], REFERI)

puts "\nDone."
