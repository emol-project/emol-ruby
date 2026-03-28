# -*- coding: utf-8 -*-
#
#  emol_defaults.rb
#
#  Default values for NameList input used by Mkinput.
#
#  This Hash follows the same structure as the items dictionary
#  returned by SimpleNameList#parse:
#    { "section_name" => { "key_name" => [ "value_string", ... ] } }
#
#  When a section or key is missing from the user's input,
#  the values defined here are automatically filled in.
#
#  Values must be written as arrays of strings to match the
#  output format of the NameList parser.
#
#  How to edit:
#    - To add a new item, add a key under the appropriate section.
#    - To change a default, modify the corresponding value.
#    - If an entire section is absent from the input, the whole
#      section is filled in from these defaults.
#

MKINPUT_DEFAULTS = {

  #---------------------------------------------------------
  #  $runflg section : execution flags for each module (0=skip, 1=run)
  #---------------------------------------------------------
  "runflg" => {
    "int"    => [ "0" ],
    "rhf"    => [ "0" ],
    "trn"    => [ "0" ],
    "expr"   => [ "0" ],
    "ci"     => [ "0" ],
    "mcscf"  => [ "0" ],
  },

  #---------------------------------------------------------
  #  $mol section : basic molecular information
  #---------------------------------------------------------
  "mol" => {
    "charge"       => [ "0" ],        # molecular charge (neutral = 0)
    "thresh_hf"    => [ "1.0e-10" ],  # HF convergence threshold
    "angular_type" => [ "0" ],        # 0 = Cartesian, 1 = spherical
    "atom_xyz"     => [],             # atomic coordinates (empty = no atoms)
    # "basis" has no default (handled as nil downstream)
  },

  #---------------------------------------------------------
  #  $expr section : CI expression settings
  #---------------------------------------------------------
  "expr" => {
    "spin"           => [ "0" ],  # spin quantum number 2S (singlet = 0)
    "fzcore"         => [ "0" ],  # number of frozen-core orbitals
    "core"           => [ "0" ],  # number of core orbitals
    "active"         => [ "0" ],  # number of active orbitals
    "external"       => [ "0" ],  # number of external orbitals
    "flg_use_symbol" => [ "0" ],  # flag: use symbolic expression
    "flg_use_csymb"  => [ "0" ],  # flag: use compressed symbolic expression
  },

  #---------------------------------------------------------
  #  $ci section : CI calculation settings
  #---------------------------------------------------------
  "ci" => {
    "nstate"    => [ "1" ],        # number of states to compute
    "thresh_ci" => [ "1.0e-10" ],  # CI convergence threshold
  },

  #---------------------------------------------------------
  #  $mcscf section : MCSCF calculation settings
  #---------------------------------------------------------
  "mcscf" => {
    "istate"        => [ "1" ],        # target state number
    "max_iteration" => [ "100" ],      # maximum number of iterations
    "thresh_mcscf"  => [ "1.0e-08" ],  # MCSCF convergence threshold
    "energy_shift"  => [ "0.0" ],      # level shift value
  },

}
