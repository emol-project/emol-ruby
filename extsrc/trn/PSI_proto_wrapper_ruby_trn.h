/// ====================================================================
/// 
///      EMOL: Hiroaki HONDA
/// 
///      SPDX-License-Identifier: LGPL-3.0-or-later
///      See COPYING.LESSER and COPYING in the top directory.
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_TRN_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_TRN_

#include "ruby.h"

VALUE inttrans_wrapper_init   ( VALUE self ) ;

VALUE inttrans_wrapper_trnpp  ( VALUE self,
          VALUE _nin, VALUE _nout, VALUE _np,  VALUE _mp,
          VALUE _a,   VALUE _b,    VALUE _ndd, VALUE _cp, VALUE _cpp ) ;

VALUE inttrans_wrapper_trnrr  ( VALUE self,
          VALUE _nin, VALUE _ntl, VALUE _np, VALUE _mp,
          VALUE _a,   VALUE _b,   VALUE _cp, VALUE _cpq ) ;

VALUE inttrans_wrapper_trnsps ( VALUE self,
          VALUE _nr, VALUE _nc, VALUE _b ) ;

#endif
