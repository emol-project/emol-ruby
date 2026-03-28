# -*- coding: utf-8 -*-
#
#  emol_nl.rb
#
#  SimpleNameList : parser for NameList-style input files
#

class SimpleNameList
    def initialize( _filename = nil, _key_field = "$", _end_field = "$end", _septoken = ",", _comments_regexp = /#|!/ )
        @filename = _filename
        @keyfield = _key_field
        @endfield = _end_field
        @septoken = _septoken
        @comments = _comments_regexp
    end
    
    ######################################################
    def readline
        if ( !@filename ) then
            fp = STDIN
        else
            fp = open( @filename )
        end
        str = ""
        while line = fp.gets do
            line = line.chop
            s = line.split( @comments )[ 0 ]
            if ( !s.nil? ) then
                str += s
            end
        end
        return str
    end
    
    ######################################################
    def get_token( str, sep )
#        printf "get_token: str: \"%s\"\n", str
        if ( str.nil? ) then
            return nil, nil
        end
        pos_sep  = str.index( sep )
        if ( pos_sep ) then
            pos_psep = pos_sep + sep.length
#            printf "get_token: psep, ppsep: %d, %d\n", pos_sep, pos_psep
            str_car = str[ 0...pos_psep ]
            if ( pos_psep == str.length ) then
                str_cdr = nil
            else
                str_cdr = str[ pos_psep...str.length ]
            end
        else
            str_car = nil
            str_cdr = str
        end
#        printf "get_token: str_car, str_cdr: \"%s\", \"%s\"\n", str_car, str_cdr
        return str_car, str_cdr
    end

    ######################################################
    def query_fieldname( str )
        name = ( str.split )[ 0 ].downcase
        return name
    end

    ######################################################
    ######################################################
    def parse_field( str_field, item )
        name_field = query_fieldname( str_field )
        strcar = ""
        strcdr = str_field[ name_field.length...str_field.length ].strip
        name = ""
        while strcar do
#            printf "parse_field: cdr:      \"%s\"\n", strcdr
            strcar, strcdr = get_token( strcdr, @septoken )
            if ( strcar ) then
                strcar = strcar[ 0...(strcar.length - @septoken.length ) ].strip
                if ( strcar.index( "=" ) ) then
                    data = []
                    strs = strcar.split( "=" )
                    name = strs[ 0 ].downcase
                    if ( !item[ name ] ) then
                        item[ name ] =  []
                        item[ name ] << strs[ 1 ].strip
                    else
                        item[ name ] << strs[ 0 ].strip
                    end
                else
                    if ( name.length > 0 ) && ( item[ name ] ) then
                        item[ name ] << strcar
                    end
                end
            end
            if ( strcdr ) then
                strcdr = strcdr.strip
#                printf "parse_field: car, cdr: \"%s\", \"%s\"\n", strcar, strcdr
            end
        end
    end

    ######################################################
    def parse_allfield( fields )
        items = Hash.new
        fields.each do |field|
            name = query_fieldname( field )
            if ( !items[ name ] ) then
                items[ name ] = Hash.new
            end
            parse_field( field, items[ name ] )
        end
        return items
    end

    ######################################################
    def split_to_field( str )
        fields = []
        pre_substr = str[ 0...str.length ].split( @keyfield )
        for i in 1...pre_substr.length
            fields << pre_substr[ i ]
        end
        return fields
    end

    ######################################################
    #  parse : read input, parse it, and return the items dictionary
    #          (original behavior, unchanged)
    ######################################################
    def parse
        strfile = readline 
        fields = split_to_field( strfile )
        items  = parse_allfield( fields )
        return items
    end

    ######################################################
    #  parse_with_defaults : parse the input and then merge
    #  default values into the result.  Keys that already
    #  exist in the user's input are kept as-is; only missing
    #  sections and keys are filled in from the defaults.
    #
    #  defaults : a Hash with the same structure as the items
    #             dictionary returned by parse, e.g.
    #             { "section" => { "key" => [ "value", ... ] } }
    #             (see emol_defaults.rb for the concrete definition)
    #
    #  Usage:
    #    nl = SimpleNameList.new( filename )
    #    items = nl.parse_with_defaults( MKINPUT_DEFAULTS )
    #
    ######################################################
    def parse_with_defaults( defaults )
        items = parse
        merge_defaults( items, defaults )
        return items
    end

    ######################################################
    #  nl_fetch (class method)
    #
    #  Safely fetch a value from the parsed items dictionary.
    #  If the section or key is missing (e.g. due to a missing
    #  entry in the defaults definition), print a warning to
    #  stderr and return a safe fallback value instead of
    #  crashing with NoMethodError on nil.
    #
    #  Parameters:
    #    items     - the parsed items dictionary
    #    section   - NameList section name (e.g. "mol")
    #    key       - key name within the section (e.g. "charge")
    #    as_array  - if true, fallback is [] (for array-valued keys)
    #                if false (default), fallback is ["0"]
    #
    #  Returns:
    #    The value array (e.g. ["1.0e-10"]) if found,
    #    or the fallback value if not found.
    #
    #  Usage:
    #    @charge = SimpleNameList.nl_fetch(items, "mol", "charge")[0].to_f
    #    atom_xyz = SimpleNameList.nl_fetch(items, "mol", "atom_xyz", true)
    #
    ######################################################
    def self.nl_fetch( items, section, key, as_array = false )
        if items[ section ] && items[ section ][ key ]
            return items[ section ][ key ]
        end
        if as_array
            fallback = []
        else
            fallback = [ "0" ]
        end
        $stderr.printf( "Warning: NameList key '%s/%s' not found " +
                        "in input or defaults. Using fallback: %s\n",
                        section, key, fallback.inspect )
        return fallback
    end

    private

    ######################################################
    #  merge_defaults : fill missing entries in the items
    #  dictionary with values from the defaults Hash.
    #
    #  Rules:
    #    - If an entire section is absent, add the whole
    #      section from defaults.
    #    - If a section exists but a key is missing, add
    #      only that key from defaults.
    #    - If a key already exists, do NOT overwrite it
    #      (user input always takes priority).
    ######################################################
    def merge_defaults( items, defaults )
        defaults.each do | section, keys |
            if !items[ section ]
                # Section itself is absent -- add the whole section
                items[ section ] = {}
                keys.each do | key, val |
                    items[ section ][ key ] = val.clone
                end
            else
                # Section exists -- fill in only the missing keys
                keys.each do | key, val |
                    if !items[ section ][ key ]
                        items[ section ][ key ] = val.clone
                    end
                end
            end
        end
    end

end
