#!/bin/sh -x

DIR=$PWD
RUBYLIB=$DIR/../src:$DIR/../modules:$RUBYLIB
EMOL_HOME="$DIR/.."
export RUBYLIB EMOL_HOME

echo 'ruby path:' `which ruby`

#ruby mod_input.rb h2o-dz.inp
#ruby mod_input.rb h2o-dzp.inp
#ruby mod_input.rb h2o-6_31G_d.inp
#ruby mod_input.rb h2o-sapTZP.inp
ruby mod_input.rb h2o-sapTZP_sph.inp
#ruby mod_input.rb ch4-6_31G_d.inp
#ruby mod_input.rb ch4-sapTZP.inp
#ruby mod_input.rb ch4-sapTZP_sph.inp
ruby mod_int.rb
ruby mod_rhf.rb
ruby mod_pop.rb
