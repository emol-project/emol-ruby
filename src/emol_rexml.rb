require_relative './emol_nl'

require "rexml/document"
include REXML

FLG_SORT_SHELLS = nil

class Basis_function
  def initialize(l, zeta, coef)
    @l = l; @zeta = zeta 
    @norm = []
    zeta.each do |x|
      @norm.push(1.0 / Math.sqrt(ovlp(l, x, x)))
    end
    @coef = normalize(l, zeta, coef)  

#    for i in 0...@zeta.length do
#      printf( "L,Zeta,Coef: %6d%28.16e%28.16e\n", l, @zeta[ i ], @coef[ i ] )
#    end
  end

  def normalize(_l, _zeta, _coef)
    __coef = _coef.clone
    pow2   = Array.new( _l + 1 )
    rfact2 = Array.new( _l + 1 )
    pow2[ 0 ] = rfact2[ 0 ] = 1.0e0
    for lqn in 1.._l do
      pow2[   lqn ] = pow2[   lqn - 1 ] * 2.0e0
      rfact2[ lqn ] = rfact2[ lqn - 1 ] * Math.sqrt( ( 1.0e0 / ( 2 * lqn - 1 ) ) )
    end

    pi    = 4.0e0 * Math.atan( 1.0e0 )
    cons  = ( 2.0e0 / pi ) ** 0.75e0
    dlqn2 = ( _l + 1.5e0 ) * 0.5e0
    coefl = cons * pow2[ _l ] * rfact2[ _l ]

    for i in 0..._zeta.length do
      #printf( "%4d: P:%28.16e, l:%28.16e, p:%28.16e, c:%28.16e, coe:%28.16e\n", i, _zeta[ i ], dlqn2, (_zeta[ i ] ** dlqn2), coefl, _coef[ i ] )
      __coef[ i ] *= coefl * ( _zeta[ i ] ** dlqn2 )
    end
    return __coef
#    sum = 0.0
#    (0..zeta.size-1).each do |i|
#      (0..i).each do |j| 
#        if i == j
#          sum += @norm[i] * @norm[i] * coef[i] * coef[i] * ovlp(l, zeta[i], zeta[i]) 
#        else
#          sum += 2.0 * @norm[i] * @norm[j] * coef[i] * coef[j] * ovlp(l, zeta[i], zeta[j]) 
#        end
#      end
#    end
#
#    norm = 1.0 / Math.sqrt(sum)
#    renormalized_coef = []
#    coef.each do |x|
#      renormalized_coef.push(x * norm) 
#    end
#    renormalized_coef
  end

  def ovlp(l, zeta_a, zeta_b)
    x = 2*l + 2
    result = 1.0 
    while x > l + 1
      result *= x
      x -= 1 
    end
    result / 2.0**(2*l + 3) * Math::PI**0.5 / (zeta_a + zeta_b)**(l + 1.5)
  end

  def get_l
    @l
  end

  def get_nterm
    @zeta.size
  end

  def get_zeta
    @zeta
  end

  def get_coef
    @coef
  end

  def get_coef_unnormalized_gtf(l)
    res = []; ang_part = Math.sqrt((2*l + 1)/4.0 / Math::PI)
    @coef.each_with_index do |x, i|
      res.push(x * @norm[i] * ang_part)
    end
    res
  end

  def show
    ang = ["S", "P", "D", "F", "G", "H"]
    printf("*** %3s  %5d \n", ang[@l], @zeta.size)
    for i in 0..@zeta.size-1 do
      printf("%12.6f  %12.6f \n", @zeta[i], @coef[i])
    end
  end
end

class Basis_set
  def initialize(bs)
    @body = bs.clone
    @nbf = bs.size
  end

  def get_body
    @body
  end
 
  def show
    printf("Size of basis set : %2d \n", @nbf)
    @body.each do |b|
      b.show
    end
  end
end

def read_basis_set(bs_name, atom)
  dir_edumolhome = ENV[ "EMOL_HOME" ]
  filename_basis = ""
  if ( bs_name ) then
      if ( dir_edumolhome ) then
         filename_basis = dir_edumolhome + "/basis/" + bs_name + "/" + atom + ".xml"
      else
         filename_basis = "./" + bs_name + "/" + atom + ".xml"
      end
  else
     filename_basis = "./" + atom + ".xml"
  end
  doc = Document.new File.new( filename_basis ) 
  #
  bs = []
  b = XPath.first(doc, "basis_set/body")
  XPath.each(b, "basis") do |bf|
    l = XPath.first(bf, "attribute::l").to_s.to_i
    zeta = []; coef = []
    XPath.each(bf, "primitive") do |p|
      zeta.push(XPath.first(p, "attribute::zeta").to_s.to_f)
      coef.push(XPath.first(p, "attribute::coef").to_s.to_f)
    end
    bs.push(Basis_function.new(l, zeta, coef))
  end
  return Basis_set.new(bs)
end

class Array
  def sum
    s = 0
    self.each do |x|
      s += x
    end
    s  
  end

  def max
    res = self[0]
    self.each do |x|
      if x > res
        res = x
      end
    end
    res
  end
end

class Mkinput
  def initialize(inp)
    #atomic_charge = {"H" => 1, "C" => 6, "N" => 7, "O" => 8}
    atomic_charge = { "H"  =>  1,
                      "He" =>  2, "HE" => 2,
                      "Li" =>  3, "LI" => 3,
                      "Be" =>  4, "BE" => 4,
                      "B"  =>  5, "C"  => 6, "N"  => 7, "O"  => 8, "F" =>  9, 
                      "Ne" => 10, "NE" => 10,
                      "Na" => 11, "NA" => 11,
                      "Mg" => 12, "MG" => 12,
                      "Al" => 13, "AL" => 13,
                      "Si" => 14, "SI" => 14,
                      "P"  => 15,
                      "S"  => 16,
                      "Cl" => 17, "CL" => 17,
                      "Ar" => 18, "AR" => 18 }
    l_name_cart = [["S"],
                  ["X", "Y", "Z"],
                  ["XX", "YY", "ZZ", "XY", "XZ", "YZ"],
                  ["XXX", "YYY", "ZZZ", "XXY", "XXZ", "XYY", "YYZ", "XZZ", "YZZ", "XYZ"],
                  ["XXXX", "YYYY", "ZZZZ", "XXXY", "XXXZ", "XYYY", "YYYZ", "XZZZ", "YZZZ", "XXYY", "XXZZ", "YYZZ", "XXYZ", "XYYZ", "XYZZ"]]
    l_name_sph = [["S"],
                  ["X", "Y", "Z"],
                  ["D-2", "D-1", "D0", "D+1", "D+2"],
                  ["F-3", "F-2", "F-1", "F0", "F+1", "F+2", "F+3"],
                  ["G-4", "G-3", "G-2", "G-1", "G0", "G+1", "G+2", "G+3", "G+4"]]
    l_degeneracy_cart = [1, 3, 6, 10, 15, 21, 28]
    l_degeneracy_sph  = [1, 3, 5,  7,  9, 11, 13]
    @numb_shell = 0; @numb_prim = 0; @numb_atm = 0
    @shel_lqn = []; @shel_atm = []; @shel_tem = []; @shel_add = []; @shel_ini = []
    @atom_charg = []; @atom_xyz = []
    @prim_exp = []; @prim_coe = []

    ########################################################################
    nl = SimpleNameList.new( inp )
    items = nl.parse
    ########################################################################
    if ( items[ "runflg" ] ) then
        @flg_int        = items[ "runflg" ][ "int"   ][ 0 ].to_i
        @flg_rhf        = items[ "runflg" ][ "rhf"   ][ 0 ].to_i
        @flg_trn        = items[ "runflg" ][ "trn"   ][ 0 ].to_i
        @flg_expr       = items[ "runflg" ][ "expr"  ][ 0 ].to_i
        @flg_ci         = items[ "runflg" ][ "ci"    ][ 0 ].to_i
        @flg_mcscf      = items[ "runflg" ][ "mcscf" ][ 0 ].to_i
    end

    if ( items[ "mol" ] ) then
        @charge         = items[ "mol" ][ "charge"    ][ 0 ].to_f
        @thresh_hf      = items[ "mol" ][ "thresh_hf" ][ 0 ].to_f
        if ( items[ "mol" ][ "basis" ] ) then
            @name_basis = items[ "mol" ][ "basis"    ][ 0 ]
        else
            @name_basis = nil
        end
        atom_and_xyz    = items[ "mol" ][ "atom_xyz" ]
        @numb_atm       = atom_and_xyz.length / 4
        @angular_type   = items[ "mol" ][ "angular_type" ][ 0 ].to_i
    end
    l_name              = (@angular_type == 1) ? l_name_sph : l_name_cart
    l_degeneracy        = (@angular_type == 1) ? l_degeneracy_sph : l_degeneracy_cart

    if ( items[ "expr" ] ) then
        @spin           = items[ "expr"  ][ "spin"     ][ 0 ].to_i
        @n_frozen       = items[ "expr"  ][ "fzcore"   ][ 0 ].to_i
        @n_core         = items[ "expr"  ][ "core"     ][ 0 ].to_i
        @n_active       = items[ "expr"  ][ "active"   ][ 0 ].to_i
        @n_external     = items[ "expr"  ][ "external" ][ 0 ].to_i
        @flg_use_symbol = items[ "expr"  ][ "flg_use_symbol" ][ 0 ].to_i
        @flg_use_csymb  = items[ "expr"  ][ "flg_use_csymb"  ][ 0 ].to_i
    end

    if ( items[ "ci" ] ) then
        @nstate         = items[ "ci"  ][ "nstate"   ][ 0 ].to_i
        @thresh_ci      = items[ "ci"  ][ "thresh_ci"    ][ 0 ].to_f
    end

    if ( items[ "mcscf" ] ) then
       @istate         = items[ "mcscf"  ][ "istate"   ][ 0 ].to_i
       @max_iteration  = items[ "mcscf"  ][ "max_iteration" ][ 0 ].to_i
       @thresh_mcscf   = items[ "mcscf"  ][ "thresh_mcscf"  ][ 0 ].to_f
       @energy_shift   = items[ "mcscf"  ][ "energy_shift"  ][ 0 ].to_f
    end
    ########################################################################
    @atom_name = []
    @atom_charg = []
    @atom_xyz = []
    for i in 0...@numb_atm
        @atom_name    << atom_and_xyz[ i * 4 ]
        @atom_charg   << atomic_charge[ @atom_name[ i ] ]
        for ixyz in 0...3
            @atom_xyz << atom_and_xyz[ i * 4 + 1 + ixyz ].to_f
        end
    end

    ########################################################################
    add = 0; ini = 0
    for iatom in 0...@numb_atm
      basis_set = read_basis_set(@name_basis, @atom_name[ iatom ]).get_body
      basis_set.each do |bf|
        @shel_lqn.push(bf.get_l)
        @shel_atm.push(iatom)
        @shel_tem.push(bf.get_nterm)
        @shel_add.push(add); add += bf.get_nterm
        @shel_ini.push(ini); ini += l_degeneracy[bf.get_l]
        exp = bf.get_zeta
        #coe = bf.get_coef_unnormalized_gtf(bf.get_l)
        coe = bf.get_coef
        for i in 0...bf.get_nterm
          @prim_exp.push(exp[i])
          @prim_coe.push(coe[i])
        end
      end
    end

    ########################################################################
    @numb_shell = @shel_lqn.size
    @numb_prim = @shel_tem.sum
    #@numb_atm = @atom_xyz.size / 3
    @maxL = @shel_lqn.max

    ########################################################################
    if ( FLG_SORT_SHELLS ) then
        sort_shells
    end
    #show

    ########################################################################
    @bs_name = []
    @bs_atm = []
    @shel_lqn.each_with_index do |lqn, i|
        atom = @shel_atm[i]
        l_name[ lqn ].each do | xyz |
            @bs_name.push( @atom_name[ atom ] + " " + (@shel_atm[ i ] + 1).to_s + " " + xyz ) 
            @bs_atm.push( atom )
        end
    end
  end
    ########################################################################
#    f = open(inp)
#    line = f.gets.chomp; puts "<input data>";  puts "> " + line
#    line = line.split
#    @charge = line[0].to_f
#    @spin = line[1].to_i
#    line = f.gets.chomp; puts "<input data>";  puts "> " + line
#    line = line.split
#    @n_frozen = line[0].to_i
#    @n_core = line[1].to_i
#    @n_active = line[2].to_i
#    @n_external = line[3].to_i
#    line = f.gets;  puts "> " + line;  printf("\n")
#    basis = line.chomp
#
#    add = 0; ini = 0
#    @atom_name = []
#    while line = f.gets
#      atom = line.split[0]
#      @atom_name.push( atom )
#      @numb_atm += 1; @atom_charg.push(atomic_charge[atom].to_f)
#      @atom_xyz.push(line.split[1].to_f) 
#      @atom_xyz.push(line.split[2].to_f)
#      @atom_xyz.push(line.split[3].to_f)
#
#      basis_set = read_basis_set(name_basis, atom).get_body
#      basis_set.each do |bf|
#        @shel_lqn.push(bf.get_l)
#        @shel_atm.push(@numb_atm - 1)
#        @shel_tem.push(bf.get_nterm)
#        @shel_add.push(add); add += bf.get_nterm
#        @shel_ini.push(ini); ini += l_degeneracy[bf.get_l]
#        exp = bf.get_zeta; coe = bf.get_coef_unnormalized_gtf(bf.get_l)
#        for i in 0..bf.get_nterm-1
#          @prim_exp.push(exp[i])
#          @prim_coe.push(coe[i])
#        end
#      end
#    end
    ########################################################################

  def get_charge
    @charge
  end

  def get_spin
    @spin
  end

  def get_bs_name
     @bs_name
  end

  def get_bs_atm
     @bs_atm
  end


  def get_n_frozen
     @n_frozen
  end

  def get_n_core
     @n_core
  end

  def get_n_active
     @n_active
  end

  def get_n_external
     @n_external
  end

  def get_numb_shell
    @numb_shell
  end

  def get_numb_prim
    @numb_prim
  end

  def get_numb_atm
    @numb_atm
  end

  def get_maxL
    @maxL
  end

  def get_shel_lqn
    @shel_lqn
  end

  def get_shel_atm
    @shel_atm
  end

  def get_shel_tem
    @shel_tem
  end

  def get_shel_add
    @shel_add
  end

  def get_shel_ini
    @shel_ini
  end

  def get_charg
    @atom_charg
  end

  def get_atom_xyz
    @atom_xyz
  end

  def get_prim_exp
    @prim_exp
  end

  def get_prim_coe
    @prim_coe
  end

  def get_angular_type
    @angular_type
  end

  def get_nstate
    @nstate
  end

  def get_istate
    @istate
  end

  def get_max_iteration
    @max_iteration
  end

  def get_thresh_hf
    @thresh_hf
  end

  def get_thresh_ci
    @thresh_ci
  end

  def get_thresh_mcscf
    @thresh_mcscf
  end

  def get_energy_shift
    @energy_shift
  end

  def get_flg_use_symbol
    @flg_use_symbol
  end

  def get_flg_use_csymb
    @flg_use_csymb
  end

  def ncomp_lqn( l )
    n = 0
    if    ( l ==  0 ) then
      n = 1
    elsif ( l ==  1 ) then
      n = 3
    elsif ( l ==  2 ) then
      n = 6
    elsif ( l ==  3 ) then
      n = 10
    else
        STDERR.printf( "ERROR: angular momentum L = %d is not supported.\n", l ) ;
        exit
    end
    return n
  end


  def sort_shells
    ### backup
    @shel_lqn_original = []
    @shel_atm_original = []
    @shel_tem_original = []
    @shel_add_original = []
    @shel_ini_original = []
    @prim_exp_original = []
    @prim_coe_original = []
    tbl_shel_sort      = []
    tbl_shel_sort_inv  = []
    tbl_prim_sort      = []
    tbl_prim_sort_inv  = []
    for i in 0...@numb_shell do
      @shel_lqn_original << @shel_lqn[ i ]
      @shel_atm_original << @shel_atm[ i ]
      @shel_tem_original << @shel_tem[ i ]
      @shel_add_original << @shel_add[ i ]
      @shel_ini_original << @shel_ini[ i ]
      tbl_shel_sort      << -1
      tbl_shel_sort_inv  << -1
    end
    for i in 0...@numb_prim do
      @prim_exp_original << @prim_exp[ i ]
      @prim_coe_original << @prim_coe[ i ]
      tbl_prim_sort      << -1
      tbl_prim_sort_inv  << -1
    end
    @numb_basis = 0
    for i in 0...@numb_shell do
      @numb_basis += ncomp_lqn( @shel_lqn[ i ] )
    end
    #printf( "numb_basis: %6d\n", @numb_basis )

    ### start sorting
    #for ( l = 0 ; l <= MAX_LQN ; l++ ) {
    #    for ( is = 0 ; is < NUMB_SHELL ; is++ ) {
    #        if ( SHEL_LQNt [ is ] == l ) {
    #            TBL_SHEL_SORT     [ is  ] = iss ; 
    #            TBL_SHEL_SORT_INV [ iss ] = is  ; 
    #            ++iss ;
    #        }
    #    }
    #}
    iss = 0
    for l in 0..@maxL do
      for is in 0...@numb_shell do
        if ( @shel_lqn_original[ is ] == l ) then
          tbl_shel_sort[     is  ] = iss
          tbl_shel_sort_inv[ iss ] = is
          iss += 1
        end
      end
    end

    #for ( is = 0 ; is < NUMB_SHELL ; is++ ) {
    #    SHEL_LQN [ TBL_SHEL_SORT [ is ] ] = SHEL_LQNt [ is ] ;
    #    SHEL_TEM [ TBL_SHEL_SORT [ is ] ] = SHEL_TEMt [ is ] ;
    #    SHEL_ATM [ TBL_SHEL_SORT [ is ] ] = SHEL_ATMt [ is ] ;
    #}
    for is in 0...@numb_shell do
      @shel_lqn[ tbl_shel_sort[ is ] ] = @shel_lqn_original[ is ]
      @shel_tem[ tbl_shel_sort[ is ] ] = @shel_tem_original[ is ]
      @shel_atm[ tbl_shel_sort[ is ] ] = @shel_atm_original[ is ]
    end
    
    #for ( is = 0, ip_add = 0, nini = 0 ; is < NUMB_SHELL ; is++ ) {
    #    SHEL_ADD [ is ] = ip_add ;
    #    SHEL_INI [ is ] = nini   ;
    #    ip_add         += SHEL_TEM [ is ] ;
    #    nini           += PSI_ncomp_lqn ( SHEL_LQN [ is ] ) ;
    #}
    ip_add = 0
    nini   = 0
    for is in 0...@numb_shell do
      @shel_add[ is ] = ip_add
      @shel_ini[ is ] = nini
      ip_add         += @shel_tem[ is ]
      nini           += ncomp_lqn( @shel_lqn[ is ] )
    end

    #ip = ipt = 0 ;
    #for ( is = 0 ; is < NUMB_SHELL ; is++ ) {
    #    for ( ip_local = 0 ; ip_local < SHEL_TEMt [ is ] ; ip_local++ ) {
    #        ipt = SHEL_ADDt [ is ] + ip_local ;
    #        ip  = SHEL_ADD  [ is ] + ip_local ;
    #        TBL_PRIM_SORT     [ ip  ] = ipt ; 
    #        TBL_PRIM_SORT_INV [ ipt ] = ip  ; 
    #    }
    #}
    ip  = 0
    ipt = 0
    for is in 0...@numb_shell do
      for ip_local in 0...@shel_tem_original[ is ] do
        ip_original = @shel_add_original[ is ] + ip_local
        ip          = @shel_add[          is ] + ip_local
        tbl_prim_sort[     ip ] = ipt
        tbl_prim_sort_inv[ ip ] = ip
      end
    end

    #prim_exp, prim_coe
    #ip = 0 ;
    #for ( is = 0 ; is < NUMB_SHELL ; is++ ) {
    #    for ( ip_local = 0 ; ip_local < SHEL_TEM [ is ] ; ip_local++, ip++ ) {
    #        PRIM_EXP [ ip ]
    #              = PRIM_EXPt            [ SHEL_ADDt [ TBL_SHEL_SORT_INV [ is ] ] + ip_local ] ;
    #        PRIM_COE_NORMALIZED [ ip ]
    #              = PRIM_COE_NORMALIZEDt [ SHEL_ADDt [ TBL_SHEL_SORT_INV [ is ] ] + ip_local ] ;
    #    }
    #}
    ip = 0
    for is in 0...@numb_shell do
      for ip_local in 0...@shel_tem[ is ] do
        #printf( "%4s, %4s,%4s: %4s, %4s, %4s\n",
        #        ip, is, ip_local, tbl_shel_sort_inv[ is ], @shel_add_original[ tbl_shel_sort_inv[ is ] ] + ip_local,
        #                               @prim_exp_original[ @shel_add_original[ tbl_shel_sort_inv[ is ] ] + ip_local ] )
        @prim_exp[ ip ] = @prim_exp_original[ @shel_add_original[ tbl_shel_sort_inv[ is ] ] + ip_local ]
        @prim_coe[ ip ] = @prim_coe_original[ @shel_add_original[ tbl_shel_sort_inv[ is ] ] + ip_local ]
        ip += 1
      end
    end

#   pow2   = Array.new( @maxL + 1 )
#   rfact2 = Array.new( @maxL + 1 )
#   pow2[ 0 ] = rfact2[ 0 ] = 1.0e0
#   for lqn in 1..@maxL do
#     pow2[   lqn ] = pow2[   lqn - 1 ] * 2.0e0
#     rfact2[ lqn ] = rfact2[ lqn - 1 ] * Math.sqrt( ( 1.0e0 / ( 2 * lqn - 1 ) ) )
#   end
#   pi = 4.0e0 * Math.atan( 1.0e0 )
#   cons = ( 2.0e0 / pi ) ** 0.75e0
#   for is in 0...@numb_shell do
#       lqn    = @shel_lqn[ is ]
#       dlqn2  = ( lqn + 1.5e0 ) * 0.5e0
#       coef   = cons * pow2[ lqn ] * rfact2[ lqn ]
#       ip_add = @shel_add[ is ]
#       for ip_local in 0...@shel_tem[ is ] do
#           ip = ip_add + ip_local
#	printf( "%4d%4d: P:%28.16e, l:%28.16e, p:%28.16e, c:%28.16e, coe:%28.16e\n",
#		is, ip, @prim_exp[ ip ], dlqn2, @prim_exp[ ip ] ** dlqn2, coef, @prim_coe[ ip ] )
#           @prim_coe[ ip ] *= coef * ( @prim_exp[ ip ] ** dlqn2 )
#       end
#   end
  end

  def show
    printf("Input data for integral sub-programs\n")
    printf("       NUMB_SHELL =  %5d \n", @numb_shell)
    printf("       NUMB_PRIM  =  %5d \n", @numb_prim)
    printf("       NUMB_ATOM  =  %5d \n", @numb_atm)
    printf("       MAXL       =  %5d \n \n", @maxL)
    print_ary(@shel_lqn, "SHEL_LQN", 10, "%5d")
    print_ary(@shel_atm, "SHEL_ATM", 10, "%5d")
    print_ary(@shel_tem, "SHEL_TEM", 10, "%5d")
    print_ary(@shel_add, "SHEL_ADD", 10, "%5d")
    print_ary(@shel_ini, "SHEL_INI", 10, "%5d")
    #print_ary(@bs_name,  "BS_NAME ", 10, "%5s")
    print "\n"
    print_ary(@atom_charg, "ATOM_CHARG", 10, "%5.2f")
    print "\n"
    print_ary(@atom_xyz, "ATOM_XYZ", 3, "%16.8f")
    print "\n"
    print_ary(@prim_exp, "PRIM_EXP", 3, "%24.16e")
    print "\n"
    print_ary(@prim_coe, "PRIM_COE", 3, "%24.16e")
    print "\n\n"
    #printf("Run flags\n")
    #printf("int, rhf, trn, expr, ci = %3s %3s %3s %3s %3s\n", @flg_int, @flg_rhf, @flg_trn, @flg_expr, @flg_ci )
    #print "\n\n"
    printf("MO classification \n")
    printf("n_frozen, n_core, n_active, n_externl = %3d %3d %3d %3d \n", @n_frozen, @n_core, @n_active, @n_external)
    printf("istate ( mcscf ) = %5d \n", @istate)
    printf("max_iteration ( mcscf ) = %5d \n", @max_iteration)
    printf("convergence thresholds : hf = %12.6e,  ci = %12.6e,  mcscf = %12.6e \n", @thresh_hf, @thresh_ci, @thresh_mcscf)
    printf("energy_shift ( mcscf ) = %12.6e \n", @energy_shift)
  end

  def print_ary(ary, var, row_size, fmt)
    printf("       %8s = [ ", var)
    if ary.size <= row_size
      print_line(ary, 0, ary.size-2, " ] \n", ary[ary.size-1], fmt)
    else
      from = 0
      while ary.size - from > row_size
        printf("\n            ")
        print_line(ary, from, row_size-2, ", \\", ary[from+row_size-1], fmt)
        from += row_size
      end

      printf("\n            ")
      print_line(ary, from, ary.size-2-from, " ] \n", ary[ary.size-1], fmt)
    end
  end

  def print_line(ary, from, to, br, ary_last, fmt)
    for i in 0..to do
      printf(fmt + ", ", ary[from + i])    
    end
    printf(fmt + br, ary_last)     
  end
end
