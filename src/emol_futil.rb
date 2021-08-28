require_relative './emol_int'
require 'gsl'

#module EmolFUtil    
    #def write_object( obj, oport )
    #   obj.
    #end

    def typeprint( value, sep, oport )
        #oport.printf "each: value: %s: \"%s\"\n", value.class, value
        if    ( value.instance_of?( Array ) ) then
            oport.printf "[ " 
            for i in 0...(value.length-1)
                typeprint( value[ i ], ", ", oport )
            end
            typeprint( value[ value.length-1 ], " ]", oport )
            oport.printf "%s", sep
        elsif ( value.instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum) ) ) then
            oport.printf "%d%s",     value, sep
        elsif ( value.instance_of?( Float  ) ) then
            oport.printf "%22.15e%s", value, sep
        elsif ( value.instance_of?( String ) ) then
            oport.printf "\"%s\"%s",  value, sep
        elsif ( value.instance_of?( NZ_eris ) ) then
            oport.printf "%8svalue = NZ_eris.new( %8d,%8d,%8d,%8d,%22.15e )\n",
                         "", value.get_i, value.get_j, value.get_k, value.get_l, value.get_value
        else
            STDERR.printf "EMOL FUTIL: Error: class = %s is not supported.\n", value
            exit
        end
    end
    
    def alltypeprint( obj, symbol_obj, blk, sep, oport )
        oport.printf "%s#Class: %s\n", blk, obj.class
        oport.printf "%s%s = ", blk, symbol_obj
        if ( obj.instance_of?( Array )       ||
             obj.instance_of?( GSL::Vector ) ||
             obj.instance_of?( GSL::Matrix ) ) then
            lenprt = 10000000000
            if    ( obj[ 0 ].instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum)) ) then
                lenprt = 10
            elsif ( obj[ 0 ].instance_of?( Float  ) ) then
                lenprt = 4 
            end
            ###
            size = 0
            if    ( obj.instance_of?( Array       ) ) then
                size = obj.length
            elsif ( obj.instance_of?( GSL::Vector ) ) then
                size = obj.size
            elsif ( obj.instance_of?( GSL::Matrix ) ) then
                size = obj.size1 * obj.size2
            else
            end
            ###
            ### write Array
            ###
            if ( obj.instance_of?( Array ) ) then
                if ( obj[ 0 ].instance_of?( NZ_eris ) ) then
                    ###
                    ### NZ_eris Array
                    ###
                    ip = 0
                    oport.printf "[]\n"
                    for i in 0...size do
                        typeprint( obj[ i ], ",", oport )
                        oport.printf "%8sobj.push( value )\n", ""
                    end
                else
                    ###
                    ### General Array
                    ###
                    #ip = 0
                    oport.printf "[ "
                    for i in 0...(size-1) do
                        #oport.printf "all: value: %s: \"%s\"\n", obj[ i ].class, obj[ i ]
                        if ( obj[ i ].instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum)) ) ||
                           ( obj[ i ].instance_of?( Float ) ) ||
                           ( obj[ i ].instance_of?( String ) ) ||
                           ( obj[ i ].instance_of?( Array ) ) then
                            typeprint( obj[ i ], ", ", oport )
                        else
                            lenprt = 5
                            oport.printf( "%s.new ( ), ", obj[ i ].class )
                        end
                        #if ( ( (ip+1) % lenprt ) == 0 ) then
                        #    oport.printf "\n"
                        #end
                        #ip += 1
                    end
                    if ( obj[ obj.length-1 ].instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum)) ) ||
                       ( obj[ obj.length-1 ].instance_of?( Float ) ) ||
                       ( obj[ obj.length-1 ].instance_of?( String ) ) ||
                       ( obj[ obj.length-1 ].instance_of?( Array ) ) then
                        typeprint( obj[ obj.length - 1 ], "", oport )
                    else
                        oport.printf( "%s.new ( ) ", obj[ obj.length - 1 ].class )
                    end

                    oport.printf " ].clone"
                    oport.printf "\n"

                    ###
                    #for i in 0...obj.length do
                    #    if ( obj[ i ].instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum)) ) || ( obj[ i ].instance_of?( Float  ) ) ||
                    #                                             ( obj[ i ].instance_of?( String ) ) then
                    #    else
                    #        p obj, obj.class
                    #        name = obj.instance_variables[ i ]
                    #        printf "name: %d: %s\n", i, name
                    #        oport.printf( "%s.%s =", symbol_obj, name[ 1, name.length - 1 ] )
                    #        typeprint( obj.instance_variable_get( obj.instance_variables[ i ] ), "\n", oport )
                    #    end
                    #end
                end
            ###
            ### write GSL::Vector
            ###
            elsif ( obj.instance_of?( GSL::Vector ) ) then
                ip = 0
                oport.printf "GSL::Vector[\n"
                for i in 0...(obj.size-1) do
                    typeprint( obj[ i ], ",", oport )
                    if ( ( (ip+1) % lenprt ) == 0 ) then
                        oport.printf "\n"
                    end
                    ip += 1
                end
                typeprint( obj[ size - 1 ], "", oport )
                oport.printf " ]"
                oport.printf "\n"
            ###
            ### write GSL::Matrix
            ###
            elsif ( obj.instance_of?( GSL::Matrix ) ) then
                #ip = 0
                oport.printf "GSL::Matrix[\n"
                for i in 0...(obj.size1) do
                    oport.printf "[ "
                    for j in 0...(obj.size2-1) do
                        typeprint( obj[ i, j ], ",", oport )
                        #if ( ( (ip+1) % lenprt ) == 0 ) then
                        #    oport.printf "\n"
                        #end
                        #ip += 1
                    end
                    #typeprint( obj[ i * obj.size2 + ( obj.size2 - 1 ) ], "", oport )
                    typeprint( obj[ i, ( obj.size2 - 1 ) ], "", oport )
                    oport.printf " ]"
                    if ( i != ( obj.size1 - 1 ) ) then
                        oport.printf ",\n"
                    else
                        oport.printf " ]\n"
                    end
                end
            else
            end
        else
            typeprint( obj, "", oport )
            oport.printf "\n"
        end
    end

    def rec_print_classobject( obj, symbol_obj, blk, oport )
        if ( obj.instance_variables.length <= 0 ) then
            # obj does not have member data
            alltypeprint( obj, symbol_obj, blk, ",", oport )
        else
            # obj has member data
            oport.printf "%s%s = %s.new( )\n", blk, symbol_obj, obj.class
            for i in 0...obj.instance_variables.length do
                instance_name_next = obj.instance_variables[ i ]
                symbol_obj_next = symbol_obj + "." +
                                  instance_name_next[ 1, instance_name_next.length - 1 ]
                rec_print_classobject( obj.instance_variable_get( obj.instance_variables[ i ] ),
                                       symbol_obj_next, blk, oport )
                ##obj_instance_name = obj.instance_variables[ i ]
                ##oport.printf "%sobj.%s = %s.new( )\n", blk,
                ##             obj_instance_name[ 1, obj_instance_name.length - 1 ],
                ##             obj.instance_variables[ i ].class
                ##if ( obj.instance_variables[ i ].instance_variable.length > 0 ) then
                ##    rec_print_classobject( obj.instance_variables[ i ], blk + "    ", oport )
                ##else
                ##    alltypeprint( obj.instance_variables[ i ], ",", oport )
                ##end
            end
        end
    end

    def write_object_contents( obj, key, modulename, oport )
        oport.printf "    def %s.%s_data( )\n", modulename, key
        rec_print_classobject( obj, "obj", "        ", oport )

#        if ( obj.instance_variables.length > 0 ) then
#            oport.printf "        obj = %s.new( )\n", obj.class
#            for i in 0...obj.instance_variables.length do
#                obj_instance_name = obj.instance_variables[ i ]
#                oport.printf "        obj.%s = ", obj_instance_name[ 1, obj_instance_name.length - 1 ]
#                object = obj.instance_variable_get( obj.instance_variables[ i ] )
#                alltypeprint( object, ",", oport )
##                if object.methods.index( "length" ) then
##                    oport.printf " [ ", obj.instance_variables[ i ]
##                    for j in 0...(object.length-1) do
##                        typeprint( object[ j ], ",", oport )
##                    end
##                    typeprint( object[ object.length - 1 ], "", oport )
##                    oport.printf " ] ", obj.instance_variables[ i ]
##                    oport.printf "\n"
##                else
##                    value = object
##                    if    ( value.instance_of?( (RUBY_VERSION.to_f >= 2.4 ? Integer : Fixnum)) ) then
##                        oport.printf "%4d\n",     value
##                    elsif ( value.instance_of?( Float )  ) then
##                        oport.printf "%25.16e\n", value
##                    elsif ( value.instance_of?( String ) ) then
##                        oport.printf "\"%s\"\n",  value
##                    else
##                        oport.printf "%s\n",      value
##                    end
##                end
#            end

        oport.printf "        return obj\n"
        oport.printf "    end\n"
#        else
#            oport.printf "        value = "
#            alltypeprint( obj, ",", oport )
#            oport.printf "        return value\n"
#        end
#        oport.printf "    end\n"
    end

    def write_object_header( modulename, oport )
        oport.printf "module %s\n", modulename
    end

    def write_object_footer( modulename, oport )
        oport.printf "end\n"
    end

    def write_multiobject( arrobj, arrkey, modulename, oport )
        write_object_header( modulename, oport )
        for i in 0...arrobj.length
            write_object_contents( arrobj[ i ], arrkey[ i ], modulename, oport )
        end
        write_object_footer( modulename, oport )
    end
    
    def read_object( obj, iport )
      return obj
    end
    
    ##################
    def fwrite_multiobject( arrobj, arrkey, modulename, filename )
        oport = File::open( filename, "w" )
        write_multiobject( arrobj, arrkey, modulename, oport )
        oport.close
    end
    
    def fread_object( obj, filename )
        iport = File::open( filename, "r" )
        read_object( obj, iport )
        iport.close
    end

    ##################
#end
