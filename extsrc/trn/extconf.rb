require 'fileutils'
require 'mkmf'
create_makefile( 'Inttrans' )

################################################
################################################
################################################
uname                   = `uname`
uname.chop!
################################################
dir_target              = "../../src/"
fort_basenames          = [ "trnpp", "trnsps", "trnrr" ]
include_additional_dir  = [ "./", "./include/", "./include/trn/", "${HOME}/include", "${HOME}/local/include", "/usr/local/include", "/opt/local/include", "/Library/Ruby/Site/1.8/universal-darwin11.0" ]
library_additional_dir  = [ "./", "./lib", "./lib/trn", "${HOME}/lib", "${HOME}/local/lib", "/usr/local/lib", "/opt/local/lib", "/opt/local/lib/ruby/vendor_ruby/1.8/i686-darwin11/", "/opt/local/lib/ruby/vendor_ruby/1.8/", "/Library/Ruby/Site/1.8/universal-darwin11.0" ]
library_additional_item = [ ]
cc                      = "gcc"
fc                      = "gfortran"
fflags                  = "-fPIC -O3"
echo                    = "echo"
qadd                    = "V = 1\nQ1 = $(V:1=)\nQ = $(Q1:0=@)"
vval                    = "1"

###############
srcfort_add = ""
#fort_basenames.each do
#    |basename|
#    srcfort_add += " " + basename + ".f"
#end
##printf "%s\n", srcfort_add
###############
objfort_add = ""
#fort_basenames.each do
#    |basename|
#    objfort_add += " " + basename + ".o"
#end
##printf "%s\n", objfort_add
###############
inc_add_flag = ""
include_additional_dir.each do
    |dir|
    inc_add_flag += "-I" + dir + " "
end
#printf "%s\n", inc_add_flag
##############
lib_add_flag = ""
library_additional_dir.each do
    |dir|
    lib_add_flag += "-L" + dir + " "
end
#printf "%s\n", inc_add_flag
##############
library_add_flag = ""
library_additional_item.each do
    |item|
    library_add_flag += "-l" + item + " "
end
#printf "%s\n", library_add_flag
##############
file         = "Makefile"
file_tmp     = "Makefile.tmp"
sbst_echo    = "ECHO = "     + echo
sbst_qadd    = ""            + qadd
sbst_vval    = "V = "        + vval
sbst_target  = "-o "         + dir_target + "$@"
sbst_srcf    = "SRCF = "     + srcfort_add
sbst_cc      = "CC = "       + cc
sbst_fc      = "FC = "       + fc
sbst_fflags  = "FFLAGS = "   + fflags
sbst_include = "INCFLAGS = " + inc_add_flag
sbst_libpath = "LIBPATH = "  + lib_add_flag
sbst_library = "LIBS = "     + library_add_flag
sbst_suffix0 = ".f.o:"
sbst_suffix1 = "	$(ECHO) compiling $(<)"
sbst_suffix2 = "	$(Q) $(FC) $(FFLAGS) -c $<"

#printf "%s\n", sbst_include

port_r = open( file,     "r" )
port_w = open( file_tmp, "w" )

port_w.puts sbst_echo
port_w.puts sbst_qadd
begin
    port_r.each_line { |line|
        new_line = line
        if ( uname == "Darwin" ) then
            new_line = new_line.gsub( /-arch[ ][ ]*i386/, "" )
            new_line = new_line.gsub( /-arch[ ][ ]*x86_64/, "" )
            new_line = new_line.gsub( /-arch[ ][ ]*x86_64/, "" )
            new_line = new_line.sub( /-Wshorten-64-to-32/, "" )
            new_line = new_line.sub( /-Wdivision-by-zero/, "" )
            new_line = new_line.sub( /-Wextra-tokens/, "" )
        end
        if new_line.index( /^CC =[ ][ ]*.*$/ ) then
            new_line = new_line.sub( /^CC =[ ][ ]*.*$/,    sbst_cc      )
            port_w.puts new_line
            port_w.puts sbst_fc
        elsif new_line.index( /^CFLAGS[ ][ ]*.*$/ ) then
            port_w.puts new_line
            port_w.puts sbst_fflags
        elsif new_line.index( /^SRCS[ ][ ]*.*$/ ) then
            port_w.puts new_line
            port_w.puts sbst_srcf
        elsif new_line.index( /^.*DLLIB.* Makefile$/ ) then
#            puts new_line
            port_w.puts sbst_suffix0
            port_w.puts sbst_suffix1
            port_w.puts sbst_suffix2
            port_w.puts ""
            port_w.puts new_line
        elsif new_line.index( /^OBJS =[ ][ ]*/ ) then
            new_line = new_line.chop + objfort_add
            port_w.puts new_line
        else
            new_line = new_line.sub( /^V =[ ][ ]*.*/,      sbst_vval    )
            new_line = new_line.sub( /-o[ ][ ]*\$@/,       sbst_target  )
            new_line = new_line.sub( /^INCFLAGS =[ ][ ]*/, sbst_include )
            new_line = new_line.sub( /^LIBPATH =[ ][ ]*/,  sbst_libpath )
            new_line = new_line.sub( /^LIBS =[ ][ ]*/,     sbst_library )
            new_line = new_line.sub( /^LIBS =[ ][ ]*/,     sbst_library )
            port_w.puts new_line
        end
    }
ensure
    port_r.close
    port_w.close
end

FileUtils.mv( file_tmp, file )
