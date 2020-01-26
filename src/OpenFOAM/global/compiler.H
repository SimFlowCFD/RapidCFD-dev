#pragma once

#if defined(__clang__)
  #define RC_DISABLE_CLANG_VIRTUAL_OVERLOAD_WARNING_BEGIN       \
    _Pragma("clang diagnostic push")                            \
    _Pragma("clang diagnostic ignored \"-Woverloaded-virtual\"")
#else
  #define RC_DISABLE_CLANG_VIRTUAL_OVERLOAD_WARNING_BEGIN
#endif

#if defined(__clang__)
  #define RC_DISABLE_CLANG_VIRTUAL_OVERLOAD_WARNING_END          \
    _Pragma("clang diagnostic pop")
#else
  #define RC_DISABLE_CLANG_VIRTUAL_OVERLOAD_WARNING_END
#endif
