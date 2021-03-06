#pragma once

#define _FILE_OFFSET_BITS 64

#if defined(WIN32)
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif defined(__APPLE__) || defined (MACOSX)
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef	off_t			off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef __off64_t off_type;
#endif
