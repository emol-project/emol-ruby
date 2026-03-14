#!/usr/bin/env ruby
# test_molint.rb
# Minimal test: H2 molecule with STO-3G basis set
# Compares overlap and kinetic integrals

$LOAD_PATH.unshift('.')
require 'Molint'

# H2 molecule: two H atoms at (0,0,0) and (0,0,1.4) bohr
numb_atom  = 2
atom_charg = [1.0, 1.0]
atom_xyz   = [0.0, 0.0, 0.0,   # H1
              0.0, 0.0, 1.4]   # H2 (1.4 bohr apart)

# STO-3G basis for H: 1s shell = 3 primitives
# exponents and contraction coefficients (already including normalization)
# From standard STO-3G: 
#   alpha:  3.42525091, 0.62391373, 0.16885540
#   coef:   0.15432897, 0.53532814, 0.44463454

sto3g_exp = [3.42525091, 0.62391373, 0.16885540]
sto3g_coe_raw = [0.15432897, 0.53532814, 0.44463454]

# Apply emol-style normalization: coef * (2/pi)^0.75 * alpha^0.75
# For s functions (l=0): normalize factor = (2/pi)^0.75 * alpha^0.75
sto3g_coe = sto3g_coe_raw.map.with_index { |c, i|
  alpha = sto3g_exp[i]
  cons = (2.0 / Math::PI) ** 0.75
  dlqn2 = (0 + 1.5) * 0.5  # l=0
  c * cons * (alpha ** dlqn2)
}

# Two H atoms, each with one s shell (3 primitives)
numb_shell = 2
numb_prim  = 6
shel_lqn   = [0, 0]        # both s shells
shel_tem   = [3, 3]        # 3 primitives each
shel_atm   = [0, 1]        # atom 0 and atom 1
shel_add   = [0, 3]        # offset into prim arrays

prim_exp = sto3g_exp + sto3g_exp
prim_coe = sto3g_coe + sto3g_coe

tol = 1.0e-14
level_print = 0

puts "=" * 60
puts "Test: H2 / STO-3G with libcint-based Molint"
puts "=" * 60
puts

# --- Overlap ---
sint = Sint0.new
sint.init(numb_atom, numb_shell, numb_prim,
          shel_lqn, shel_tem, shel_atm, shel_add,
          atom_charg, atom_xyz, prim_exp, prim_coe,
          tol, level_print)

puts "Overlap integrals:"
for ish in 0...numb_shell
  for jsh in 0..ish
    inttype, nsize, s_int = sint.calc(ish, jsh)
    printf("  S(%d,%d) = %20.14f  (inttype=%d, nsize=%d)\n",
           ish, jsh, s_int[0], inttype, nsize)
  end
end
puts

# --- Kinetic ---
tint = Tint0.new
tint.init(numb_atom, numb_shell, numb_prim,
          shel_lqn, shel_tem, shel_atm, shel_add,
          atom_charg, atom_xyz, prim_exp, prim_coe,
          tol, level_print)

puts "Kinetic energy integrals:"
for ish in 0...numb_shell
  for jsh in 0..ish
    inttype, nsize, t_int = tint.calc(ish, jsh)
    printf("  T(%d,%d) = %20.14f\n", ish, jsh, t_int[0])
  end
end
puts

# --- Nuclear attraction ---
vint = Vint0.new
vint.init(numb_atom, numb_shell, numb_prim,
          shel_lqn, shel_tem, shel_atm, shel_add,
          atom_charg, atom_xyz, prim_exp, prim_coe,
          tol, level_print)

puts "Nuclear attraction integrals:"
for ish in 0...numb_shell
  for jsh in 0..ish
    inttype, nsize, v_int = vint.calc(ish, jsh)
    printf("  V(%d,%d) = %20.14f\n", ish, jsh, v_int[0])
  end
end
puts

# --- ERI ---
eris = Eri0.new
eris.init(numb_atom, numb_shell, numb_prim,
          shel_lqn, shel_tem, shel_atm, shel_add,
          atom_charg, atom_xyz, prim_exp, prim_coe,
          tol, level_print)

puts "Two-electron integrals (selected):"
for ish in 0...numb_shell
  for jsh in 0..ish
    for ksh in 0..ish
      lsh_max = (ksh < ish) ? ksh : jsh
      for lsh in 0..lsh_max
        inttype, nsize, g_int = eris.calc(ish, jsh, ksh, lsh)
        ash, bsh, csh, dsh = eris.exchange_order
        printf("  (%d %d|%d %d) [->(%d %d|%d %d)] = %20.14f\n",
               ish, jsh, ksh, lsh, ash, bsh, csh, dsh, g_int[0])
      end
    end
  end
end

puts
puts "=" * 60
puts "Reference values (Szabo & Ostlund, Table 3.3):"
puts "  S(0,0) = 1.0,  S(1,0) = 0.6593"
puts "  T(0,0) = 0.7600,  T(1,0) = 0.2365"
puts "  V(0,0) = -1.8827,  V(1,0) = -1.1962 (sum of both nuclei)"
puts "  (00|00) = 0.7746"
puts "=" * 60

