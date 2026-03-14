#!/bin/sh

DIR=$PWD
RUBYLIB=$DIR/../test:$DIR/../src:$RUBYLIB
EDUMOL_HOME="$DIR/.."
export RUBYLIB EDUMOL_HOME
#######################################################################

ruby test_molint.rb

#ruby main_input.rb h2o.dz.inp
#ruby main_rhf.rb h2o.dz.inp
#ruby main_pop.rb
