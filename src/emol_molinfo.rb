require_relative "emol_rexml"

class Molinfo
    attr_accessor :charge, :spin, :numb_shell, :numb_atm, :numb_prim, :angular_type,
                  :shel_lqn, :shel_atm, :shel_tem, :shel_add, :shel_ini,
                  :atom_charg, :atom_xyz, :prim_exp, :prim_coe, :bs_name, :bs_atm,
                  :thr_ovch, :level_print, :len_lqn, :nel, :nbf, :e_nucl,
                  :n_frozen, :n_core, :n_active, :n_external, :nob, :e_core,
                  :nstate, :flg_use_symbol, :flg_use_csymb,
		  :istate, :max_iteration, :energy_shift,
		  :thresh_hf, :thresh_ci, :thresh_mcscf
  def initialize( obj = nil )
      if obj.nil? then
          @charge        = 0.0e0
          @spin          = 0
          @numb_shell    = 0
          @numb_atm      = 0
          @numb_prim     = 0
          @angular_type  = 0
          @shel_lqn      = []
          @shel_atm      = []
          @shel_tem      = []
          @shel_add      = []
          @shel_ini      = []
          @atom_charg    = []
          @atom_xyz      = []
          @prim_exp      = []
          @prim_coe      = []
          @thr_ovch      = 0.0e0
          @level_print   = 0
          @len_lqn       = []
          @nel           = 0
          @nbf           = 0
          @e_nucl        = 0.0e0
 
          @n_frozen      = 0
          @n_core        = 0
          @n_active      = 0
          @n_external    = 0
          @nob           = 0
          @e_core        = 0.0e0
          @nstate        = 0

	        @istate        = 0
	        @max_iteration = 0
	        @energy_shift  = 0.0
      
	        @thresh_hf     = 0.0
	        @thresh_ci     = 0.0
	        @thresh_mcscf = 0.0
	  
          return
      end
      if ( obj.instance_of?( String ) ) then
          fname = obj
          inp = Mkinput.new( fname )
          inp.show
      
          @charge       = inp.get_charge
          @spin         = inp.get_spin
          @numb_shell   = inp.get_numb_shell
          @numb_atm     = inp.get_numb_atm
          @numb_prim    = inp.get_numb_prim
          @angular_type = inp.get_angular_type
          @shel_lqn     = inp.get_shel_lqn
          @shel_atm     = inp.get_shel_atm
          @shel_tem     = inp.get_shel_tem
          @shel_add     = inp.get_shel_add
          @shel_ini     = inp.get_shel_ini
          @bs_name      = inp.get_bs_name
          @bs_atm       = inp.get_bs_atm
          @atom_charg   = inp.get_charg
          @atom_xyz     = inp.get_atom_xyz
          @prim_exp     = inp.get_prim_exp
          @prim_coe     = inp.get_prim_coe
          @thr_ovch     = 1.0e+30
          @level_print  = 0
      
          @len_lqn      = [1, 3, 6, 10, 15, 21, 28]
          @angular_type = inp.respond_to?(:get_angular_type) ? inp.get_angular_type : 0
          if (@angular_type == 1) then
            @len_lqn    = [1, 3, 5,  7,  9, 11, 13]    # spherical
          else
            @len_lqn    = [1, 3, 6, 10, 15, 21, 28]   # cartesian
          end
          @nel          = (@atom_charg.sum - @charge).to_i
          @nbf          = @shel_ini[ @shel_ini.size - 1 ] + @len_lqn[ @shel_lqn[ @shel_lqn.size - 1 ]]
      
          @e_nucl = 0.0
          for i in 1..@numb_atm-1
            for j in 0..i-1
              @e_nucl  += @atom_charg[i] * @atom_charg[j] / Math.sqrt( (@atom_xyz[i*3  ] - @atom_xyz[j*3  ])**2 +
                                                                       (@atom_xyz[i*3+1] - @atom_xyz[j*3+1])**2 +
                                                                       (@atom_xyz[i*3+2] - @atom_xyz[j*3+2])**2 )
            end
          end

          @n_frozen       = inp.get_n_frozen
          @n_core         = inp.get_n_core
          @n_active       = inp.get_n_active
          @n_external     = inp.get_n_external
          @nob            = @n_core + @n_active + @n_external
          @e_core         = 0.0e0
          @flg_use_symbol = inp.get_flg_use_symbol
          @flg_use_csymb  = inp.get_flg_use_csymb
          @nstate         = inp.get_nstate

	        @thresh_hf      = inp.get_thresh_hf
	        @thresh_ci      = inp.get_thresh_ci
	        @thresh_mcscf   = inp.get_thresh_mcscf
	        @istate         = inp.get_istate
	        @max_iteration  = inp.get_max_iteration
	        @energy_shift   = inp.get_energy_shift

      elsif ( obj.instance_of?( Molinfo ) ) then
          @charge         = obj.charge
          @spin           = obj.spin
          @numb_shell     = obj.numb_shell
          @numb_atm       = obj.numb_atm
          @numb_prim      = obj.numb_prim
          @angular_type   = obj.angular_type
          @shel_lqn       = obj.shel_lqn
          @shel_atm       = obj.shel_atm
          @shel_tem       = obj.shel_tem
          @shel_add       = obj.shel_add
          @shel_ini       = obj.shel_ini
          @atom_charg     = obj.atom_charg
          @atom_xyz       = obj.atom_xyz
          @prim_exp       = obj.prim_exp
          @prim_coe       = obj.prim_coe
          @thr_ovch       = obj.thr_ovch
          @level_print    = obj.level_print
          @len_lqn        = obj.len_lqn
          @bs_name        = obj.bs_name
          @bs_atm         = obj.bs_atm
          @nel            = obj.nel
          @nbf            = obj.nbf
          @e_nucl         = obj.e_nucl
   
          @n_frozen       = obj.n_frozen
          @n_core         = obj.n_core
          @n_active       = obj.n_active
          @n_external     = obj.n_external
          @nob            = obj.nob
          @e_core         = obj.e_core
          @flg_use_symbol = obj.flg_use_symbol
          @flg_use_csymb  = obj.flg_use_csymb 
          @nstate         = obj.nstate
	        @thresh_hf      = obj.thresh_hf
	        @thresh_ci      = obj.thresh_ci
	        @thresh_mcscf   = obj.thresh_mcscf
       
	        @istate         = obj.istate
	        @max_iteration  = obj.max_iteration
	        @energy_shift   = obj.energy_shift
      else  
      end
  end

  def get_charge
    @charge
  end

  def get_spin
    @spin
  end

  def get_numb_shell
    @numb_shell
  end

  def get_numb_atm
    @numb_atm
  end

  def get_numb_prim
    @numb_prim
  end

  def get_angular_type
    @angular_type
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

  def get_atom_charg
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


  def get_thr_ovch
    @thr_ovch
  end

  def get_level_print
    @level_print
  end

  def get_len_lqn
    @len_lqn
  end

  def get_bs_name
    @bs_name
  end

  def get_bs_atm
    @bs_atm
  end

  def get_nel
    @nel
  end

  def get_nbf
    @nbf
  end

  def get_e_nucl
    @e_nucl
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

  def get_nob
    @nob
  end

  def get_e_core
    @e_core
  end

  def set_e_core( en )
    @e_core = en
  end

  def get_nstate
    @nstate
  end

  def get_max_iteration
    @max_iteration
  end

  def get_istate
    @istate
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
end

