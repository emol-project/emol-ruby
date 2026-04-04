#!/usr/bin/env ruby
# test_sph_fg.rb
# ==============================================================
# Equivalent to gen_ref_sph_fg.py  (GINT SCALE-free convention)
#
# System: He(0,0,0) -- H(1.0 Ang, 0, 0)
#   ISH=0: f on He  NPRIM=3, alpha=[1.5,0.8,0.4], raw_coef=[0.2,0.5,0.4]
#   ISH=1: g on He  NPRIM=2, alpha=[2.0,0.6],     raw_coef=[0.5,0.7]
#   ISH=2: s on H   NPRIM=1, alpha=[1.0],          raw_coef=[1.0]
#
# angular_type = 1 (Spherical)
#
# Reference values: PySCF raw (no SCALE factor).
#   = gen_ref_sph_fg.py output / SCALE_product
#   where SCALE(l) = sqrt(4*pi/(2l+1))
#   f: SCALE_f = sqrt(4*pi/7),  g: SCALE_g = sqrt(4*pi/9),  s: 1
#
# This convention matches libcint output directly.
# fint_capi.c (GINT backend) divides out the SCALE after GINT_SPH,
# so both backends now return the same values.
# ==============================================================

#require_relative 'Molint'
$LOAD_PATH.unshift('.')
require 'Molint'

# --------------------------------------------------------------
# emol_coe = raw * gto_norm(l, alpha) * sqrt((2l+1)/(4*pi))
# --------------------------------------------------------------
def gto_norm(l, alpha)
  n1 = (2 * l + 3) * 0.5
  gi = Math.gamma(n1) / (2.0 * (2.0 * alpha) ** n1)
  1.0 / Math.sqrt(gi)
end

def emol_coe(l, alpha, raw)
  raw * gto_norm(l, alpha) * Math.sqrt((2.0*l+1.0) / (4.0*Math::PI))
end

alpha_f = [1.5, 0.8, 0.4];  raw_f = [0.2, 0.5, 0.4]
alpha_g = [2.0, 0.6];        raw_g = [0.5, 0.7]
alpha_s = [1.0];              raw_s = [1.0]

coe_f = alpha_f.zip(raw_f).map { |a,r| emol_coe(3, a, r) }
coe_g = alpha_g.zip(raw_g).map { |a,r| emol_coe(4, a, r) }
coe_s = alpha_s.zip(raw_s).map { |a,r| emol_coe(0, a, r) }

puts "--- emol_coe ---"
puts "f: #{coe_f.map{|x| x.round(10)}}"
puts "g: #{coe_g.map{|x| x.round(10)}}"
puts "s: #{coe_s.map{|x| x.round(10)}}"

# --------------------------------------------------------------
# Molecular data
# --------------------------------------------------------------
ANG2BOHR = 1.8897259886

numb_atom  = 2
numb_shell = 3
numb_prim  = 6

shel_lqn = [3, 4, 0]
shel_tem = [3, 2, 1]
shel_atm = [0, 0, 1]
shel_add = [0, 3, 5]

atom_charg = [2.0, 1.0]
atom_xyz   = [0.0, 0.0, 0.0, 1.0*ANG2BOHR, 0.0, 0.0]

prim_exp = alpha_f + alpha_g + alpha_s
prim_coe = coe_f   + coe_g   + coe_s

tol          = 1.0e-10
level_dbg    = 0
angular_type = 1   # Spherical

ISH_F = 0;  ISH_G = 1;  ISH_S = 2

# --------------------------------------------------------------
# Reference values: PySCF raw (no SCALE)
# = gen_ref_sph_fg.py values / SCALE_product
# --------------------------------------------------------------

# 1E f-f diagonal (S and T same for all m)
REF_SFF = 1.0
REF_TFF = 3.189823152748915

REF_VFF = {
  -3 => -1.787437547525, -2 => -1.719251195221, -1 => -1.675749422602,
   0 => -1.715638038646, +1 => -1.795373637185, +2 => -1.747548931481,
  +3 => -1.846121931828
}

# 1E g-g diagonal
REF_SGG = 1.0
REF_TGG = 5.679172505718465

REF_VGG = {
  -4 => -1.882868115664, -3 => -1.822536786701, -2 => -1.799814817489,
  -1 => -1.760203262915,  0 => -1.809347707267, +1 => -1.848262574758,
  +2 => -1.845172324448, +3 => -1.834420358395, +4 => -1.916690811072
}

# 1E f-s [S, T, V]  (nonzero at m=+1, +3 only)
REF_FS = {
  -3 => [0.0, 0.0, 0.0], -2 => [0.0, 0.0, 0.0], -1 => [0.0, 0.0, 0.0],
   0 => [0.0, 0.0, 0.0],
  +1 => [-0.228528755312, -0.538130730561, +0.667709109048],
  +2 => [0.0, 0.0, 0.0],
  +3 => [+0.295029354484, +0.694723785849, -0.862008753151]
}

# 1E g-s [S, T, V]  (nonzero at m=0, +2, +4)
REF_GS = {
  -4 => [0.0, 0.0, 0.0], -3 => [0.0, 0.0, 0.0], -2 => [0.0, 0.0, 0.0],
  -1 => [0.0, 0.0, 0.0],
   0 => [+0.086739382303, +0.286787577778, -0.293735196859],
  +1 => [0.0, 0.0, 0.0],
  +2 => [-0.129303436770, -0.427517679343, +0.437874578374],
  +3 => [0.0, 0.0, 0.0],
  +4 => [+0.171052368680, +0.565552730312, -0.579253619906]
}

# 2E (fm fm|ss)
REF_FF_SS = {
  -3=>0.510231026171, -2=>0.458701131171, -1=>0.416499717813,
   0=>0.449975185083, +1=>0.495792144836, +2=>0.476755558900, +3=>0.520519487279
}

# 2E (gm gm|ss)
REF_GG_SS = {
  -4=>0.520543256594, -3=>0.479832985335, -2=>0.458252672310,
  -1=>0.422262432282,  0=>0.456356233638, +1=>0.494588846475,
  +2=>0.475637886884, +3=>0.484207082428, +4=>0.522865796967
}

# 2E (fm s|fm s) exchange
REF_FS_FS = {
  -3=>0.039145459013, -2=>0.006235656435, -1=>0.003077900579,
   0=>0.015959171448, +1=>0.078146998214, +2=>0.026264188143, +3=>0.126087892733
}

# 2E (gm s|gm s) exchange
REF_GS_GS = {
  -4=>0.035940518080, -3=>0.009394047117, -2=>0.006776567405,
  -1=>0.001509061016,  0=>0.019103740585, +1=>0.023788871966,
  +2=>0.037541965464, +3=>0.018928213520, +4=>0.065552267062
}

# --------------------------------------------------------------
# Helper
# --------------------------------------------------------------
def check(label, got, ref, tol=1.0e-5)
  diff = (got - ref).abs
  ok   = diff < tol ? "OK" : "*** FAIL ***"
  printf("  %-30s got=%+.12f  ref=%+.12f  diff=%.2e  %s\n",
         label, got, ref, diff, ok)
end

# --------------------------------------------------------------
# Init
# --------------------------------------------------------------
args = [numb_atom, numb_shell, numb_prim,
        shel_lqn, shel_tem, shel_atm, shel_add,
        atom_charg, atom_xyz, prim_exp, prim_coe,
        tol, level_dbg, angular_type]

sint = Sint0.new; sint.init(*args)
tint = Tint0.new; tint.init(*args)
vint = Vint0.new; vint.init(*args)
eri0 = Eri0.new;  eri0.init(*args)

# --------------------------------------------------------------
# 1E f-f diagonal
# --------------------------------------------------------------
puts "\n=== 1E f-f (ISH=0 vs ISH=0) ==="
_, _, s_ff = sint.calc(ISH_F, ISH_F)
_, _, t_ff = tint.calc(ISH_F, ISH_F)
_, _, v_ff = vint.calc(ISH_F, ISH_F)

(-3..3).each do |m|
  i = m + 3
  diag = i*7 + i
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("S(f#{ms},f#{ms})", s_ff[diag], REF_SFF)
  check("T(f#{ms},f#{ms})", t_ff[diag], REF_TFF)
  check("V(f#{ms},f#{ms})", v_ff[diag], REF_VFF[m])
end

# --------------------------------------------------------------
# 1E g-g diagonal
# --------------------------------------------------------------
puts "\n=== 1E g-g (ISH=1 vs ISH=1) ==="
_, _, s_gg = sint.calc(ISH_G, ISH_G)
_, _, t_gg = tint.calc(ISH_G, ISH_G)
_, _, v_gg = vint.calc(ISH_G, ISH_G)

(-4..4).each do |m|
  i = m + 4
  diag = i*9 + i
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("S(g#{ms},g#{ms})", s_gg[diag], REF_SGG)
  check("T(g#{ms},g#{ms})", t_gg[diag], REF_TGG)
  check("V(g#{ms},g#{ms})", v_gg[diag], REF_VGG[m])
end

# --------------------------------------------------------------
# 1E f-s
# --------------------------------------------------------------
puts "\n=== 1E f-s (ISH=0 vs ISH=2) ==="
_, _, s_fs = sint.calc(ISH_F, ISH_S)
_, _, t_fs = tint.calc(ISH_F, ISH_S)
_, _, v_fs = vint.calc(ISH_F, ISH_S)

(-3..3).each do |m|
  i = m + 3
  rs, rt, rv = REF_FS[m]
  ms = m>=0 ? "+#{m}" : "#{m}"
  next if rs == 0.0 && rt == 0.0
  check("S(f#{ms},s)", s_fs[i], rs)
  check("T(f#{ms},s)", t_fs[i], rt)
  check("V(f#{ms},s)", v_fs[i], rv)
end
puts "  Zero components (S):"
(-3..3).each do |m|
  next unless REF_FS[m][0] == 0.0
  i = m + 3
  ms = m>=0 ? "+#{m}" : "#{m}"
  printf("    m=%s: S=%.3e  (expect 0)\n", ms, s_fs[i])
end

# --------------------------------------------------------------
# 1E g-s
# --------------------------------------------------------------
puts "\n=== 1E g-s (ISH=1 vs ISH=2) ==="
_, _, s_gs = sint.calc(ISH_G, ISH_S)
_, _, t_gs = tint.calc(ISH_G, ISH_S)
_, _, v_gs = vint.calc(ISH_G, ISH_S)

(-4..4).each do |m|
  i = m + 4
  rs, rt, rv = REF_GS[m]
  ms = m>=0 ? "+#{m}" : "#{m}"
  next if rs == 0.0 && rt == 0.0
  check("S(g#{ms},s)", s_gs[i], rs)
  check("T(g#{ms},s)", t_gs[i], rt)
  check("V(g#{ms},s)", v_gs[i], rv)
end

# --------------------------------------------------------------
# 2E (ff|ss)
# --------------------------------------------------------------
puts "\n=== 2E (f f|s s) ==="
_, _, e_ffss = eri0.calc(ISH_F, ISH_F, ISH_S, ISH_S)
puts "  order: #{eri0.exchange_order.inspect}"

(-3..3).each do |m|
  i = m + 3
  idx = ((i*7+i)*1+0)*1+0
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("(f#{ms} f#{ms}|ss)", e_ffss[idx], REF_FF_SS[m])
end

# --------------------------------------------------------------
# 2E (gg|ss)
# --------------------------------------------------------------
puts "\n=== 2E (g g|s s) ==="
_, _, e_ggss = eri0.calc(ISH_G, ISH_G, ISH_S, ISH_S)
puts "  order: #{eri0.exchange_order.inspect}"

(-4..4).each do |m|
  i = m + 4
  idx = ((i*9+i)*1+0)*1+0
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("(g#{ms} g#{ms}|ss)", e_ggss[idx], REF_GG_SS[m])
end

# --------------------------------------------------------------
# 2E (fs|fs) exchange
# --------------------------------------------------------------
puts "\n=== 2E (f s|f s) exchange ==="
_, _, e_fsfs = eri0.calc(ISH_F, ISH_S, ISH_F, ISH_S)
puts "  order: #{eri0.exchange_order.inspect}"

(-3..3).each do |m|
  i = m + 3
  idx = i*7 + i   # nca=7, ncb=1, ncc=7, ncd=1
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("(f#{ms} s|f#{ms} s)", e_fsfs[idx], REF_FS_FS[m])
end

# --------------------------------------------------------------
# 2E (gs|gs) exchange
# --------------------------------------------------------------
puts "\n=== 2E (g s|g s) exchange ==="
_, _, e_gsgs = eri0.calc(ISH_G, ISH_S, ISH_G, ISH_S)
puts "  order: #{eri0.exchange_order.inspect}"

(-4..4).each do |m|
  i = m + 4
  idx = i*9 + i   # nca=9, ncb=1, ncc=9, ncd=1
  ms = m>=0 ? "+#{m}" : "#{m}"
  check("(g#{ms} s|g#{ms} s)", e_gsgs[idx], REF_GS_GS[m])
end

puts "\nDone."
