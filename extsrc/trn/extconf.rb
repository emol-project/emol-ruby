# extconf.rb for Molint (GINT backend)
#
# Directory layout:
#   molint/
#     extconf.rb          <- this file
#     molint_libgint.c
#     fint_capi.h         <- copied here by: cd libgint && make
#     libgint.a           <- built by:       cd libgint && make
#     libgint/
#       Makefile
#       fint_capi.h
#       fint_capi.c
#       Fmtmake5.c
#       *.f
#
# Usage:
#   cd molint && ruby extconf.rb
#   make                  # builds libgint.a, Molint.so, copies to ../../src/
#
require 'mkmf'

here = __dir__

libgint_dir = File.join(here, 'libgint')
libgint_a   = File.join(here, 'libgint.a')

# ------------------------------------------------------------------
# 1. Include path: current directory (for fint_capi.h)
# ------------------------------------------------------------------
$INCFLAGS += " -I#{here}"

# ------------------------------------------------------------------
# 2. Link libgint.a (full path)
# ------------------------------------------------------------------
$LOCAL_LIBS += " #{libgint_a}"

# ------------------------------------------------------------------
# 3. Link gfortran runtime and libm
# ------------------------------------------------------------------
def find_gfortran_lib(name)
  path = `gfortran -print-file-name=#{name} 2>/dev/null`.strip
  (path && path != name && File.exist?(path)) ? path : nil
end

if RUBY_PLATFORM =~ /darwin/
  %w[libgfortran.a libquadmath.a].each do |lib|
    path = find_gfortran_lib(lib)
    if path
      $LOCAL_LIBS += " #{path}"
    else
      $libs = append_library($libs, lib.sub(/^lib/, '').sub(/\.a$/, ''))
    end
  end
else
  # Linux: try static first, fall back to -lgfortran
  path = find_gfortran_lib('libgfortran.a')
  if path
    $LOCAL_LIBS += " #{path}"
  else
    $libs = append_library($libs, 'gfortran')
  end
  # quadmath may be needed on some Linux toolchains
  qpath = find_gfortran_lib('libquadmath.a')
  if qpath
    $LOCAL_LIBS += " #{qpath}"
  else
    $libs = append_library($libs, 'quadmath')
  end
end

$libs = append_library($libs, 'm')

# ------------------------------------------------------------------
# 4. Source files and headers
# ------------------------------------------------------------------
$srcs = %w[molint_libgint.c]
$hdrs = %w[fint_capi.h]

# ------------------------------------------------------------------
# 5. Compiler flags
# ------------------------------------------------------------------
$CFLAGS += ' -O2 -std=c99 -fPIC'

# ------------------------------------------------------------------
# 6. Create Makefile
# ------------------------------------------------------------------
create_makefile('Molint')

# ------------------------------------------------------------------
# 7. Patch the generated Makefile:
#    a) Rewrite 'all:' to include build_libgint and install_to_src
#    b) Append the extra rules
# ------------------------------------------------------------------
target_dir = File.expand_path('../../src', here)
so_ext     = (RUBY_PLATFORM =~ /darwin/) ? 'bundle' : 'so'

makefile = File.read('Makefile')

# mkmf generates a line like:  all: $(DLLIB)
# Rewrite to:  all: build_libgint $(DLLIB) install_to_src
makefile.sub!(/^all:(.*)$/) do
  deps = $1.strip.split
              .reject { |d| ['build_libgint', 'install_to_src'].include?(d) }
              .join(' ')
  "all: build_libgint #{deps} install_to_src"
end

File.write('Makefile', makefile)

# Append extra rules
File.open('Makefile', 'a') do |f|
  f.puts
  f.puts '# --- libgint build rules (inserted by extconf.rb) ---'
  f.puts "LIBGINT_DIR = #{libgint_dir}"
  f.puts "LIBGINT_A   = #{libgint_a}"
  f.puts
  f.puts '.PHONY: build_libgint clean-libgint install_to_src'
  f.puts
  f.puts 'build_libgint:'
  f.puts "\t$(MAKE) -C $(LIBGINT_DIR)"
  f.puts
  f.puts 'clean-libgint:'
  f.puts "\t$(MAKE) -C $(LIBGINT_DIR) clean"
  f.puts
  f.puts 'clean: clean-libgint'
  f.puts
  f.puts '# --- post-build copy ---'
  f.puts "TARGET_DIR = #{target_dir}"
  f.puts
  f.puts 'install_to_src: $(DLLIB)'
  f.puts "\t@mkdir -p $(TARGET_DIR)"
  f.puts "\t@cp $(DLLIB) $(TARGET_DIR)/"
  f.puts "\t@echo '>>> Molint.#{so_ext} copied to $(TARGET_DIR)/'"
end
