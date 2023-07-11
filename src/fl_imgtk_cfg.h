#ifndef __FLIMGTK_CONFIG_H__
#define __FLIMGTK_CONFIG_H__

#ifdef USING_OMP
    #include <omp.h>
#endif /// of USING_OMP

// OpenMP compatibility for M$VC.
#ifdef _MSC_VER
    #define OMPSIZE_T       long
#else
    #define OMPSIZE_T       size_t
#endif 

#define FLIMGTK_MAX_WIDTH   (500000)
#define FLIMGTK_MAX_HEIGHT  (500000)

#endif /// of __FLIMGTK_CONFIG_H__
