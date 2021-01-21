#ifndef __FLIMGTK_CONFIG_H__
#define __FLIMGTK_CONFIG_H__

// if this flag defined, all new Fl_RGB_Image may alloc_array set to 1.
#define FLIMGTK_IMGBUFF_OWNALLOC 

#ifdef USING_OMP
    #include <omp.h>
#endif /// of USING_OMP

// OpenMP compatibility for M$VC.
#ifdef _MSC_VER
    #define OMPSIZE_T       long
#else
    #define OMPSIZE_T       size_t
#endif 

#endif /// of __FLIMGTK_CONFIG_H__