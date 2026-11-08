// ====================================================================
// 
//      EMOL: Hiroaki HONDA
// 
//      SPDX-License-Identifier: LGPL-3.0-or-later
//      See COPYING.LESSER and COPYING in the top directory.
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_TRN_
#define _INCLUDE_PSI_PROTO_TRN_

typedef int PSI_LI ;
typedef double PSI_LD ;

PSI_LI trnpp_  ( PSI_LI *nin, PSI_LI *nout, PSI_LI *np,  PSI_LI *mp,
                 PSI_LD *a,   PSI_LD *b,    PSI_LI *ndd, PSI_LD *cp, PSI_LD *cpq) ;
//
PSI_LI trnrr_  ( PSI_LI *nin, PSI_LI *ntl,  PSI_LI *np,  PSI_LI *mp,
                 PSI_LD *a,   PSI_LD *b,    PSI_LD *cp,  PSI_LD *cpq) ;
//
PSI_LI trnsps_ ( PSI_LI *nr,  PSI_LI *nc,   PSI_LD *b ) ;

#endif
