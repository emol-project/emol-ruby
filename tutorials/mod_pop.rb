
require 'emol_molinfo'
require 'emol_util'
require 'emol_rhf'

load 'Data_molinfo'
load 'Data_int'
load 'Data_rhf'

molinfo = Data_molinfo::molinfo_data
s       = Data_int::s_data
cao     = Data_rhf::cao_data
nao     = molinfo.nbf
natm    = molinfo.numb_atm
nocc    = molinfo.nel / 2
atms    = molinfo.bs_atm
aoname  = molinfo.bs_name

p 'molinfo:',     molinfo
p 'overlap:',     s
p 'rhfvec:',      cao
p 'nao:',         nao
p 'natm:',        natm
p 'nocc:',        nocc
p 'atm[iao]:',    atms
p 'aoname[iao]:', aoname

pop_ao = Array.new( nao )
pop    = Array.new( natm, 0.0e0 )
occ    = Array.new( nocc, 2.0e0 )
for iao in 0...nao do
    pop_ao[ iao ] = 0
    for iob in 0...nocc do
        val_ao = 0
        for jao in 0...nao do
            val_ao += cao[ iao, iob ] * cao[ jao, iob ] * s[ iao, jao ]
        end
        pop_ao[ iao ] += val_ao * occ[ iob ]
    end
end

for iao in 0...nao do
    iatm = atms[ iao ]
    pop[ iatm ] += pop_ao[ iao ]
end

pop_total = 0
for iatm in 0...natm do
    pop_total += pop[ iatm ]
end

printf( "Atomic Orbital Populations:\n" ) ;
for i in 0...nao do
    printf( "%6d %-10s %f\n", i, aoname[ i ], pop_ao[ i ] )
end

printf( "Atomic Populations:\n" ) ;
for i in 0...natm do
    printf( "  %4d %10.4f\n", i, pop[ i ] )
end
printf( "Total: %10.4f\n", pop_total ) ;

print "\n Population Analysis     : completed.\n"
