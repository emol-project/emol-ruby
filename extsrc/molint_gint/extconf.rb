# extconf_gint.rb  --  mkmf configuration for Molint (GINT backend)
#
# Directory layout:
#   workdir/
#     molint_libgint.c    <- C extension source
#     gint_c_api.h        <- Fortran wrapper declarations
#     extconf_gint.rb     <- this file
#     libgint/
#       blkdta000.f       \
#       fint_api.f         | GINT Fortran sources
#       ...               /
#       Makefile          <- builds ../libgint.a
#
# Usage:
#   ruby extconf_gint.rb
#   make

require 'mkmf'

# ----------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------
WORKDIR     = File.expand_path(__dir__).freeze
LIBGINT_DIR = File.join(WORKDIR, 'libgint').freeze
LIBGINT_MK  = File.join(LIBGINT_DIR, 'Makefile').freeze
LIBGINT_A   = File.join(WORKDIR, 'libgint.a').freeze

unless Dir.exist?(LIBGINT_DIR)
  abort "*** ERROR: #{LIBGINT_DIR} not found."
end
unless File.exist?(LIBGINT_MK)
  abort "*** ERROR: #{LIBGINT_MK} not found."
end

# ----------------------------------------------------------------
# 2. gfortran runtime libraries
# ----------------------------------------------------------------
FC = (ENV['FC'] || 'gfortran').freeze

gfortran_lib = `#{FC} -print-file-name=libgfortran.a 2>/dev/null`.strip
unless gfortran_lib.empty? || gfortran_lib == 'libgfortran.a'
  $LDFLAGS += " -L#{File.dirname(gfortran_lib)}"
end
$libs = append_library($libs, 'gfortran')

quadmath = `#{FC} -print-file-name=libquadmath.a 2>/dev/null`.strip
unless quadmath.empty? || quadmath == 'libquadmath.a'
  $libs = append_library($libs, 'quadmath')
end

$libs = append_library($libs, 'm')

# ----------------------------------------------------------------
# 3. Static-link libgint.a
# ----------------------------------------------------------------
$LOCAL_LIBS += " #{LIBGINT_A}"
$INCFLAGS   += " -I#{WORKDIR}"
$CFLAGS     += ' -O2 -std=c99 -fPIC'

# ----------------------------------------------------------------
# 4. Generate Makefile
# ----------------------------------------------------------------
create_makefile('Molint')

# ----------------------------------------------------------------
# 5. Patch the generated Makefile
#
#    mkmf generates:  all: $(DLLIB)
#    We change it to: all: build_libgint $(DLLIB)
#    and add the build_libgint rule.
#
#    This guarantees libgint.a is built before molint_libgint.c
#    is linked into Molint.so.
# ----------------------------------------------------------------
makefile = File.read('Makefile')

# Replace "all: $(DLLIB)" with "all: build_libgint $(DLLIB)"
makefile.sub!(/^all:[ \t]+\$\(DLLIB\)/) do
  "all: build_libgint $(DLLIB) install_to_src"
end

# Insert build_libgint and clean-libgint rules right after the all: line
insert_rules = <<~MAKE

  # --- libgint.a rules (inserted by extconf_gint.rb) ---
  LIBGINT_DIR = #{LIBGINT_DIR}
  LIBGINT_A   = #{LIBGINT_A}

  build_libgint:
  \t$(MAKE) -C $(LIBGINT_DIR)

  clean-libgint:
  \t$(MAKE) -C $(LIBGINT_DIR) clean

  clean: clean-libgint
  # ---------------------------------------------------

MAKE

# Insert after the first blank line following the all: line
makefile.sub!(/^all:.*\n/) { |m| m + insert_rules }

File.write('Makefile', makefile)

# ----------------------------------------------------------------
# 6. Post-build copy rule: install Molint.so into ../../src/
# ----------------------------------------------------------------
target_dir = File.expand_path('../../src', __dir__)
so_ext = (RUBY_PLATFORM =~ /darwin/) ? 'bundle' : 'so'

open('Makefile', 'a') do |f|
  f.puts
  f.puts "# --- post-build copy (appended by extconf_gint.rb) ---"
  f.puts "TARGET_DIR = #{target_dir}"
  f.puts
  f.puts "install_to_src: $(DLLIB)"
  f.puts "\t@mkdir -p $(TARGET_DIR)"
  f.puts "\t@cp $(DLLIB) $(TARGET_DIR)/"
  f.puts "\t@echo '>>> Molint.#{so_ext} copied to $(TARGET_DIR)/'"
end

puts ""
puts "Generated Makefile."
puts "  libgint dir : #{LIBGINT_DIR}"
puts "  libgint.a   : #{LIBGINT_A}"
puts "Run 'make' to build libgint.a and Molint.so"
