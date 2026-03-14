# extconf.rb for Molint (libcint backend)
#
# Usage:
#   ruby extconf.rb
#   make          # builds Molint.so and copies it to ../../src/
#
require 'mkmf'

# 1. Search paths for macOS Homebrew
if RUBY_PLATFORM =~ /darwin/
  ['/opt/homebrew/opt', '/usr/local/opt'].each do |brew_path|
    if Dir.exist?(brew_path)
      ['openblas', 'lapack', 'libcint'].each do |pkg|
        dir_config(pkg, "#{brew_path}/#{pkg}/include", "#{brew_path}/#{pkg}/lib")
      end
    end
  end
end

# 2. Find BLAS / LAPACK (prefer OpenBLAS over reference Netlib)
found_linear_algebra = false

if have_library('openblas', 'dgemm_')
  found_linear_algebra = true
  have_library('lapack', 'dgesv_')  # may be needed separately
else
  puts "OpenBLAS not found. Trying reference BLAS and LAPACK..."
  if have_library('blas', 'dgemm_') && have_library('lapack', 'dgesv_')
    found_linear_algebra = true
  end
end

unless found_linear_algebra
  abort "*** ERROR: BLAS/LAPACK not found. ***\n" \
        "Please install: \n" \
        "  Debian/Ubuntu: libopenblas-dev or liblapack-dev\n" \
        "  Fedora: openblas-devel or lapack-devel\n" \
        "  macOS: brew install openblas lapack"
end

# 3. Find libcint (check LIBCINT_DIR env variable first)
if ENV['LIBCINT_DIR']
  dir_config('cint', "#{ENV['LIBCINT_DIR']}/include", "#{ENV['LIBCINT_DIR']}/lib")
  $LDFLAGS += " -Wl,-rpath,#{ENV['LIBCINT_DIR']}/lib"
else
  dir_config('cint')
end

$libs = append_library($libs, 'm')

unless have_library('cint', 'cint2e_cart')
  abort "*** ERROR: libcint not found. ***"
end

# 4. Check headers
unless have_header('cint.h')
  abort "*** ERROR: cint.h not found. ***"
end

have_header('cint_funcs.h')  # present in v6, absent in v3 (optional)

# 5. Compiler flags
$CFLAGS += " -O3 -std=c99 -fPIC"

create_makefile('Molint')

# 6. Append a post-build copy rule to install Molint.so into ../../src/
target_dir = File.expand_path('../../src', __dir__)

so_ext = (RUBY_PLATFORM =~ /darwin/) ? 'bundle' : 'so'

open('Makefile', 'a') do |f|
  f.puts
  f.puts "# --- post-build copy rule (appended by extconf.rb) ---"
  f.puts "TARGET_DIR = #{target_dir}"
  f.puts
  f.puts "all: install_to_src"
  f.puts
  f.puts "install_to_src: $(DLLIB)"
  f.puts "\t@mkdir -p $(TARGET_DIR)"
  f.puts "\t@cp $(DLLIB) $(TARGET_DIR)/"
  f.puts "\t@echo '>>> Molint.#{so_ext} copied to $(TARGET_DIR)/'"
end
