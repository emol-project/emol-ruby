#!/bin/sh

MODULES='molint_libcint trn wig symbolic'
OPT=
#RUBY=$HOME/ruby-2.4.1/bin/ruby
RUBY=ruby
MAKE=make

if [ x$1 = 'xclean' ] ; then
    OPT=clean
fi

for M in $MODULES
do
    (
        echo $M:
        cd $M
        $RUBY extconf.rb
        $MAKE $OPT
    )
done
