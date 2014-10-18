/*******************************************************************************
 *	FILE: GNUC_Extensions.hpp
 *	DESC: Compiler-independent extension definitions
 *	DATE: Wed 17 Aug 2011 10:15:14 PM EDT
 *	MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __GNUC_EXTENSIONS_HPP__INCLUDED__
#define __GNUC_EXTENSIONS_HPP__INCLUDED__


// GNU C only!
#ifdef __GNUC__

#	ifndef __nowarn_unused
#		define __nowarn_unused	__attribute__((unused))
#	endif //__nowarn_unused

#	ifndef __deprecated
#		define __deprecated		__attribute__((deprecated))
#	endif //__deprecated

#	ifndef __error
#		define __error(x)		__attribute__((error(x)))
#	endif //__error

#	ifndef __warning
#		define __warning(x)		__attribute__((warning(x)))
#	endif //__warning

#	ifndef __always_inline
#		define __always_inline	__attribute__((always_inline))
#	endif //__always_inline

#	ifndef __flatten
#		define __flatten		__attribute__((flatten))
#	endif //__flatten

#	ifndef __hot
#		define __hot			__attribute__((hot))
#	endif //__hot

#	ifndef __cold
#		define __cold			__attribute__((cold))
#	endif //__cold

#	ifndef __fastcall
#		define __fastcall		__attribute__((fastcall))
#	endif //__fastcall

#	ifndef __thiscall
#		define __thiscall		__attribute__((thiscall))
#	endif //__thiscall

#	ifndef __malloc
#		define __malloc			__attribute__((malloc))
#	endif //__malloc

#   ifndef __likely
#       define __likely(x)      __builtin_expect (!!(x), 1)
#   endif //__likely

#   ifndef __unlikely
#       define __unlikely(x)    __builtin_expect (!!(x), 0)
#   endif //__unlikely

#else //__GNUC__

#	ifndef __nowarn_unused
#		define __nowarn_unused
#	endif //__nowarn_unused

#	ifndef __deprecated
#		define __deprecated
#	endif //__deprecated

#	ifndef __error
#		define __error(x)
#	endif //__error

#	ifndef __warning
#		define __warning(x)
#	endif //__warning

#	ifndef __always_inline
#		define __always_inline
#	endif //__always_inline

#	ifndef __flatten
#		define __flatten
#	endif //__flatten

#	ifndef __hot
#		define __hot
#	endif //__hot

#	ifndef __cold
#		define __cold
#	endif //__cold

#	ifndef __fastcall
#		define __fastcall
#	endif //__fastcall

#	ifndef __thiscall
#		define __thiscall
#	endif //__thiscall

#	ifndef __malloc
#		define __malloc
#	endif //__malloc

#   ifndef __likely
#       define __likely(x)
#   endif //__likely

#   ifndef __unlikely
#       define __unlikely(x)
#   endif //__unlikely

#endif //__GNUC__


#endif //!__GNUC_EXTENSIONS_HPP__INCLUDED__

