#ifndef __FL_SMIMG_H__
#define __FL_SMIMG_H__

#include <cmath>

#include <FL/Fl_RGB_Image.h>

////////////////////////////////////////////////////////////////////////////////
//
// This class was belongs to below project:
//
// ==========================================================
// FreeImage 3
//
// Design and implementation by
// - Floris van den Berg (flvdberg@wxs.nl)
// - Herv?Drolon (drolon@infonie.fr)
//
// ==========================================================
// Modified for FLTK by rageworx@gmail.com
//

////////////////////////////////////////////////////////////////////////////////
// Filters

class GenericFilter
{
    protected:

        #define FILTER_PI  double (3.1415926535897932384626433832795)
        #define FILTER_2PI double (2.0 * 3.1415926535897932384626433832795)
        #define FILTER_4PI double (4.0 * 3.1415926535897932384626433832795)

        double  m_dWidth;

    public:
        GenericFilter (double dWidth) : m_dWidth (dWidth) {}
        virtual ~GenericFilter() {}

        double GetWidth()                   { return m_dWidth; }
        void   SetWidth (double dWidth)     { m_dWidth = dWidth; }
        virtual double Filter (double dVal) = 0;
};

class BoxFilter : public GenericFilter
{
    public:
        // Default fixed width = 0.5
        BoxFilter() : GenericFilter(0.5) {}
        virtual ~BoxFilter() {}

    public:
        double Filter (double dVal) { return (fabs(dVal) <= m_dWidth ? 1.0 : 0.0); }
};

class BilinearFilter : public GenericFilter
{
    public:
        BilinearFilter () : GenericFilter(1) {}
        virtual ~BilinearFilter() {}

    public:
        double Filter (double dVal)
        {
            dVal = fabs(dVal);
            return (dVal < m_dWidth ? m_dWidth - dVal : 0.0);
        }
};

class BicubicFilter : public GenericFilter
{
    protected:
        // data for parameterized Mitchell filter
        double p0, p2, p3;
        double q0, q1, q2, q3;

    public:
        // Default fixed width = 2
        BicubicFilter (double b = (1/(double)3), double c = (1/(double)3)) : GenericFilter(2)
        {
            p0 = (6 - 2*b) / 6;
            p2 = (-18 + 12*b + 6*c) / 6;
            p3 = (12 - 9*b - 6*c) / 6;
            q0 = (8*b + 24*c) / 6;
            q1 = (-12*b - 48*c) / 6;
            q2 = (6*b + 30*c) / 6;
            q3 = (-b - 6*c) / 6;
        }
        virtual ~BicubicFilter() {}

    public:
        double Filter(double dVal)
        {
            dVal = fabs(dVal);
            if(dVal < 1)
                return (p0 + dVal*dVal*(p2 + dVal*p3));
            if(dVal < 2)
                return (q0 + dVal*(q1 + dVal*(q2 + dVal*q3)));
            return 0;
        }
};

class CatmullRomFilter : public GenericFilter
{
    public:

        // Default fixed width = 2
        CatmullRomFilter() : GenericFilter(2) {}
        virtual ~CatmullRomFilter() {}

    public:
        double Filter(double dVal)
        {
            if(dVal < -2) return 0;
            if(dVal < -1) return (0.5*(4 + dVal*(8 + dVal*(5 + dVal))));
            if(dVal < 0)  return (0.5*(2 + dVal*dVal*(-5 - 3*dVal)));
            if(dVal < 1)  return (0.5*(2 + dVal*dVal*(-5 + 3*dVal)));
            if(dVal < 2)  return (0.5*(4 + dVal*(-8 + dVal*(5 - dVal))));
            return 0;
        }
};

class Lanczos3Filter : public GenericFilter
{
    public:
        // Default fixed width = 3
        Lanczos3Filter() : GenericFilter(3) {}
        virtual ~Lanczos3Filter() {}

    public:
        double Filter(double dVal)
        {
            dVal = fabs(dVal);
            if(dVal < m_dWidth)
            {
                return (sinc(dVal) * sinc(dVal / m_dWidth));
            }
            return 0;
        }

    private:
        double sinc(double value)
        {
            if(value != 0)
            {
                value *= FILTER_PI;
                return (sin(value) / value);
            }
            return 1;
        }
};

class BSplineFilter : public GenericFilter
{
    public:
        // Default fixed width = 2
        BSplineFilter() : GenericFilter(2) {}
        virtual ~BSplineFilter() {}

    public:
        double Filter(double dVal)
        {
            dVal = fabs(dVal);
            if(dVal < 1) return (4 + dVal*dVal*(-6 + 3*dVal)) / 6;
            if(dVal < 2)
            {
                double t = 2 - dVal;
                return (t*t*t / 6);
            }
            return 0;
        }
};

////////////////////////////////////////////////////////////////////////////////
// Resize relations.

class WeightsTable
{
    typedef struct
    {
        double *Weights;
        unsigned Left, Right;
    }Contribution;

    private:
        Contribution*   m_WeightTable;
        unsigned        m_WindowSize;
        unsigned        m_LineLength;

    public:
        WeightsTable(GenericFilter *pFilter, unsigned uDstSize, unsigned uSrcSize);
        ~WeightsTable();

    public:
        double getWeight(unsigned dst_pos, unsigned src_pos);
        unsigned getLeftBoundary(unsigned dst_pos);
        unsigned getRightBoundary(unsigned dst_pos);
};

class ResizeEngine
{
    private:
        GenericFilter* m_pFilter;

    public:
        ResizeEngine( GenericFilter* filter );
        virtual ~ResizeEngine() {}

    public:
        Fl_RGB_Image* scale(Fl_RGB_Image *src, unsigned dst_width, unsigned dst_height);

    public:
        // cindex : 0 = RED, 1 = GREEEN, 2 = BLUE, 3 = ALPHA
        void useSingleChannel( bool f, char cindex );
        bool useSingleChannel() { return useSCh; }
        char refSingleChannel() { return refSCh; }

    private:
        void horizontalFilter( const uchar* src, const unsigned height, const unsigned src_width,
                               const unsigned src_bpp,
                               const unsigned src_offset_x, const unsigned src_offset_y,
                               uchar* dst, const unsigned dst_width);
        void verticalFilter( const uchar* src, const unsigned width, const unsigned src_height,
                             const unsigned src_bpp,
                             const unsigned src_offset_x, const unsigned src_offset_y,
                             uchar* dst, const unsigned dst_width, const unsigned dst_height);

    protected:
        bool useSCh;
        int  refSCh;
};


#endif /// of __FL_SMIMG_H__
