#!/bin/sh

MODULES='molint trn wig symbolic'
#MODULES='molint'
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
