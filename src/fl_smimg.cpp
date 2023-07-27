#ifdef _MSC_VER
    // Holy mother F-hawk M$VC ...
    #pragma warning(disable : 4018)
    #pragma warning(disable : 4068)
    #pragma warning(disable : 4244)
    #pragma warning(disable : 4996)
#endif

#include <FL/fl_draw.H>
#include "fl_smimg.h"
#include "fl_imgtk_minmax.h"
#include "fl_imgtk_cfg.h"

// --- formal RGB and RGBA ---
#define FI_RGBA_RED             0
#define FI_RGBA_GREEN           1
#define FI_RGBA_BLUE            2
#define FI_RGBA_ALPHA           3
// --- related in 1 or 2 depth ---
#define FI_GRAYA_GRAY           0
#define FI_GRAYA_ALPHA          1

/// Clamp function
template <class T> T CLAMP(const T &value, const T &min_value, const T &max_value) 
{
    return ((value < min_value) ? min_value : (value > max_value) ? max_value : value);
}

WeightsTable::WeightsTable( GenericFilter *pFilter, unsigned uDstSize, unsigned uSrcSize )
{
    OMPSIZE_T       u;
    double          dWidth;
    double          dFScale         = 1.0;
    const double    dFilterWidth    = pFilter->GetWidth();

    // scale factor
    const double dScale = double(uDstSize) / double(uSrcSize);

    if(dScale < 1.0)
    {
        // minification
        dWidth  = dFilterWidth / dScale;
        dFScale = dScale;
    }
    else
    {
        // magnification
        dWidth= dFilterWidth;
    }

    // allocate a new line contributions structure
    //
    // window size is the number of sampled pixels
    m_WindowSize = 2 * (int)ceil(dWidth) + 1;
    m_LineLength = uDstSize;
     // allocate list of contributions
    m_WeightTable = new Contribution[ m_LineLength + 1 ];
    for(u = 0 ; u < m_LineLength ; u++)
    {
        // allocate contributions for every pixel
        m_WeightTable[u].Weights = new double[ m_WindowSize + 1 ];
    }

    // offset for discrete to continuous coordinate conversion
    const double dOffset = ( 0.5 / dScale) - 0.5;

    #pragma omp parallel for
    for(u = 0; u < m_LineLength; u++)
    {
        // scan through line of contributions
        const double dCenter = (double)u / dScale + dOffset;   // reverse mapping

        // find the significant edge points that affect the pixel
        int iLeft  = MAX( 0, (int)floor (dCenter - dWidth) );
        int iRight = MIN( (int)ceil (dCenter + dWidth), int(uSrcSize) - 1 );

        // cut edge points to fit in filter window in case of spill-off
        if((iRight - iLeft + 1) > int(m_WindowSize))
        {
            if(iLeft < (int(uSrcSize) - 1 / 2))
            {
                iLeft++;
            }
            else
            {
                iRight--;
            }
        }

        m_WeightTable[u].Left  = iLeft;
        m_WeightTable[u].Right = iRight;

        int iSrc = 0;
        double dTotalWeight = 0;  // zero sum of weights
        for(iSrc = iLeft; iSrc <= iRight; iSrc++)
        {
            // calculate weights
            const double weight = dFScale * pFilter->Filter(dFScale * (dCenter - (double)iSrc));
            m_WeightTable[u].Weights[iSrc-iLeft] = weight;
            dTotalWeight += weight;
        }

        if((dTotalWeight > 0) && (dTotalWeight != 1))
        {
            // normalize weight of neighbouring points
            for(iSrc = iLeft; iSrc <= iRight; iSrc++)
            {
                // normalize point
                m_WeightTable[u].Weights[iSrc-iLeft] /= dTotalWeight;
            }

            // simplify the filter, discarding null weights at the right
            iSrc = iRight - iLeft;
            while(m_WeightTable[u].Weights[iSrc] == 0)
            {
                m_WeightTable[u].Right--;
                iSrc--;
                if(m_WeightTable[u].Right == m_WeightTable[u].Left)
                {
                    break;
                }
            }

        }
    }
}

WeightsTable::~WeightsTable()
{
    for(unsigned u = 0; u < m_LineLength; u++)
    {
        // free contributions for every pixel
        delete[] m_WeightTable[u].Weights;
    }
    // free list of pixels contributions
    delete[] m_WeightTable;
}

double WeightsTable::getWeight(unsigned dst_pos, unsigned src_pos)
{
    if ( dst_pos < m_LineLength )
    {
        int sz = m_WeightTable[dst_pos].Right - m_WeightTable[dst_pos].Left;
        if ( src_pos < m_WindowSize )
        {
            return m_WeightTable[dst_pos].Weights[src_pos];
        }
    }

    return 0.0;
}

unsigned WeightsTable::getLeftBoundary(unsigned dst_pos)
{
    return m_WeightTable[dst_pos].Left;
}

unsigned WeightsTable::getRightBoundary(unsigned dst_pos)
{
    return m_WeightTable[dst_pos].Right;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

ResizeEngine::ResizeEngine( GenericFilter* filter )
 : m_pFilter(filter),
   useSCh(false),
   refSCh(0)
{
}

void ResizeEngine::useSingleChannel( bool f, char cindex )
{
    useSCh = f;

    if ( useSCh == true )
    {
        switch ( cindex )
        {
            default:
            case 0 :
                refSCh = FI_RGBA_RED;
                break;

            case 1 :
                refSCh = FI_RGBA_GREEN;
                break;

            case 2 :
                refSCh = FI_RGBA_BLUE;
                break;

            case 3 :
                refSCh = FI_RGBA_ALPHA;
                break;
        }
    }
    else
    {
        refSCh = 0;
    }
}

Fl_RGB_Image* ResizeEngine::scale( Fl_RGB_Image* src, unsigned dst_width, unsigned dst_height )
{
    if ( src == NULL)
        return NULL;

    if ( ( src->data_w() == 0 ) && ( src->data_h() == 0 ) )
        return NULL;

    if ( ( src->data_w() == dst_width) && ( src->data_h() == dst_height))
    {
        // fire-egg announced this may cause segment fault,
        // let me test it -
        return (Fl_RGB_Image*)src->copy();
    }

    // allocate the dst image
    uchar* dst_buff = new uchar[ ( dst_width * dst_height * 4 ) + 1 ];
    if ( dst_buff == NULL )
    {
        return NULL;
    }

    const uchar* src_buff = (uchar*)src->data()[0];

    if ( dst_width <= src->data_w() )
    {
        uchar* tmp_buff = NULL;

        if ( src->data_w() != dst_width )
        {
            if ( src->data_h() != dst_height )
            {
                tmp_buff = new uchar[ ( dst_width * src->data_h() * src->d() ) + 1 ];
                if ( tmp_buff == NULL )
                {
                    delete[] dst_buff;
                    return NULL;
                }
            }
            else
            {
                tmp_buff = dst_buff;
            }

            horizontalFilter( src_buff, src->data_h(), src->data_w(), src->d(), 0, 0, tmp_buff, dst_width );

        }
        else
        {
            tmp_buff = (uchar*)src_buff;
        }

        if ( src->h() != dst_height )
        {
            verticalFilter( tmp_buff, dst_width, src->data_h(), src->d(), 0, 0,
                            dst_buff, dst_width, dst_height );
        }

        if ( ( tmp_buff != src_buff ) && ( tmp_buff != dst_buff ) )
        {
            delete[] tmp_buff;
            tmp_buff = NULL;
        }

    }
    else    /// == ( dst_width > src->w() )
    {
        uchar* tmp_buff = NULL;

        if ( src->h() != dst_height )
        {
            if ( src->w() != dst_width )
            {
                tmp_buff = new uchar[ src->data_w() * dst_height * src->d() + 1 ];
                if ( tmp_buff == NULL )
                {
                    delete[] dst_buff;
                    return NULL;
                }
            } else {
                tmp_buff = dst_buff;
            }

            verticalFilter( src_buff, src->data_w(), src->data_h(), src->d(),
                            0, 0, tmp_buff, dst_width, dst_height );

        }
        else
        {
            tmp_buff = (uchar*)src_buff;
        }

        if ( src->w() != dst_width )
        {
            horizontalFilter( tmp_buff, dst_height, src->data_w(), src->d(),
                              0, 0, dst_buff, dst_width );
        }

        if ( ( tmp_buff != src_buff ) && ( tmp_buff != dst_buff ) )
        {
            delete[] tmp_buff;
            tmp_buff = NULL;
        }
    }

    if ( dst_buff != NULL )
    {
        Fl_RGB_Image *dst = new Fl_RGB_Image( dst_buff, dst_width, dst_height, src->d() );
        if ( dst == NULL )
        {
            delete[] dst_buff;
            return NULL;
        }
        dst->alloc_array = 1;
        return dst;
    }

    return NULL;
}

void ResizeEngine::horizontalFilter( const uchar* src, const unsigned height, const unsigned src_width, const unsigned src_bpp, const unsigned src_offset_x, const unsigned src_offset_y, uchar* dst, const unsigned dst_width)
{
    // allocate and calculate the contributions
    WeightsTable weightsTable(m_pFilter, dst_width, src_width);

    switch ( src_bpp )
    {
        case 1: /// single Y channel, gray image.
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(x)
            for ( y = 0; y < height; y++)
            {
                const \
                uchar* src_bits = &src[ ( ( y + src_offset_y ) * src_width * src_bpp ) +
                                        ( src_offset_x * src_bpp ) ];
                uchar* dst_bits = &dst[ y * dst_width * src_bpp ];

                // scale each row
                for ( x=0; x<dst_width; x++ )
                {
                    // loop through row
                    // retrieve left boundary
                    const unsigned iLeft  = weightsTable.getLeftBoundary(x);
                    // retrieve right boundary
                    const unsigned iLimit = weightsTable.getRightBoundary(x) - iLeft;   
                    const uchar*   pixel  = src_bits + iLeft * src_bpp;
                    double v = 0;

                    // for(i = iLeft to iRight)
                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(x, i);

                        // 1, 2 depth avoids useSCh flag.
                        v += (weight * (double)pixel[FI_GRAYA_GRAY]);
                        pixel += src_bpp;
                    }

                    // clamp and place result in destination pixel
                    // and single Y channel don't need useSCh.
                    dst_bits[FI_GRAYA_GRAY]   = (uchar)CLAMP<int>((int)(v + 0.5), 0, 0xFF);
                    dst_bits += src_bpp;
                }
            }
        }
        break;

        case 2: /// single Y and Alpha chaneel.
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(x)
            for ( y=0; y<height; y++ )
            {
                const \
                uchar* src_bits = &src[ ( ( y + src_offset_y ) * src_width * src_bpp ) +
                                        ( src_offset_x * src_bpp ) ];
                uchar* dst_bits = &dst[ y * dst_width * src_bpp ];

                // scale each row
                for ( x = 0; x < dst_width; x++)
                {
                    // loop through row
                    const unsigned iLeft = weightsTable.getLeftBoundary(x);             // retrieve left boundary
                    const unsigned iLimit = weightsTable.getRightBoundary(x) - iLeft;   // retrieve right boundary
                    const uchar *pixel = src_bits + iLeft * src_bpp;
                    double v = 0, a = 0;

                    // for(i = iLeft to iRight)
                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(x, i);

                        // 1 and 2 depth avoids useSCh flag.
                        v += (weight * (double)pixel[FI_GRAYA_GRAY]);
                        a += (weight * (double)pixel[FI_GRAYA_ALPHA]);
                        pixel += src_bpp;
                    }

                    // clamp and place result in destination pixel
                    dst_bits[FI_GRAYA_GRAY] = (uchar)CLAMP<int>((int)(v + 0.5), 0, 0xFF);
                    dst_bits[FI_GRAYA_ALPHA] = (uchar)CLAMP<int>((int)(a + 0.5), 0, 0xFF);
                    dst_bits += src_bpp;
                }
            }
        }
        break;
        
        case 3:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(x)
            for ( y = 0; y < height; y++)
            {
                const
                uchar* src_bits = &src[ ( ( y + src_offset_y ) * src_width * src_bpp ) +
                                        ( src_offset_x * src_bpp ) ];
                uchar* dst_bits = &dst[ y * dst_width * src_bpp ];

                // scale each row
                for ( x = 0; x < dst_width; x++)
                {
                    // loop through row
                    const unsigned iLeft  = weightsTable.getLeftBoundary(x);
                    const unsigned iLimit = weightsTable.getRightBoundary(x) - iLeft;
                    const uchar*   pixel  = src_bits + iLeft * src_bpp;
                    double r = 0, g = 0, b = 0;

                    // for(i = iLeft to iRight)
                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(x, i);

                        if ( useSCh == true )
                        {
                            double c = (weight * (double)pixel[refSCh]);
                            r += c;
                            g += c;
                            b += c;
                        }
                        else
                        {
                            r += (weight * (double)pixel[FI_RGBA_RED]);
                            g += (weight * (double)pixel[FI_RGBA_GREEN]);
                            b += (weight * (double)pixel[FI_RGBA_BLUE]);
                        }

                        pixel += src_bpp;
                    }

                    // clamp and place result in destination pixel
                    if ( useSCh == true )
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int)(r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_RED]   = cmpv;
                        dst_bits[FI_RGBA_GREEN] = cmpv;
                        dst_bits[FI_RGBA_BLUE]  = cmpv;
                    }
                    else
                    {
                        dst_bits[FI_RGBA_RED]   = (uchar)CLAMP<int>((int)(r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_GREEN] = (uchar)CLAMP<int>((int)(g + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_BLUE]  = (uchar)CLAMP<int>((int)(b + 0.5), 0, 0xFF);
                    }
                    dst_bits += src_bpp;
                }
            }
        }
        break;

        case 4:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(x)
            for ( y = 0; y < height; y++)
            {
                const
                uchar* src_bits = &src[ ( ( y + src_offset_y ) * src_width * src_bpp ) +
                                        ( src_offset_x * src_bpp ) ];
                uchar* dst_bits = &dst[ y * dst_width * src_bpp ];

                // scale each row
                for ( x=0; x<dst_width; x++)
                {
                    // loop through row
                    const unsigned iLeft = weightsTable.getLeftBoundary(x);             // retrieve left boundary
                    const unsigned iLimit = weightsTable.getRightBoundary(x) - iLeft;   // retrieve right boundary
                    const uchar *pixel = src_bits + iLeft * src_bpp;
                    double r = 0, g = 0, b = 0, a = 0;

                    // for(i = iLeft to iRight)
                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(x, i);

                        if ( useSCh == true )
                        {
                            double c = (weight * (double)pixel[refSCh]);
                            r += c;
                            g += c;
                            b += c;
                            a += (weight * (double)pixel[FI_RGBA_ALPHA]);
                        }
                        else
                        {
                            r += (weight * (double)pixel[FI_RGBA_RED]);
                            g += (weight * (double)pixel[FI_RGBA_GREEN]);
                            b += (weight * (double)pixel[FI_RGBA_BLUE]);
                            a += (weight * (double)pixel[FI_RGBA_ALPHA]);
                        }
                        pixel += src_bpp;
                    }

                    // clamp and place result in destination pixel
                    if ( useSCh == true )
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int)(r + 0.5), 0, 0xFF);

                        dst_bits[FI_RGBA_RED]   = cmpv;
                        dst_bits[FI_RGBA_GREEN] = cmpv;
                        dst_bits[FI_RGBA_BLUE]  = cmpv;
                        dst_bits[FI_RGBA_ALPHA] = (uchar)CLAMP<int>((int)(a + 0.5), 0, 0xFF);
                    }
                    else
                    {
                        dst_bits[FI_RGBA_RED]   = (uchar)CLAMP<int>((int)(r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_GREEN] = (uchar)CLAMP<int>((int)(g + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_BLUE]  = (uchar)CLAMP<int>((int)(b + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_ALPHA] = (uchar)CLAMP<int>((int)(a + 0.5), 0, 0xFF);
                    }
                    dst_bits += src_bpp;
                }
            }
        }
        break;

    } /// of switch()
}

/// Performs vertical image filtering
void ResizeEngine::verticalFilter( const uchar* src, unsigned width, unsigned src_height, const unsigned src_bpp, unsigned src_offset_x, unsigned src_offset_y, uchar* dst, const unsigned dst_width, unsigned dst_height)
{
    // allocate and calculate the contributions
    WeightsTable weightsTable( m_pFilter, dst_height, src_height );

    //unsigned dst_pitch = dst_width * src_bpp;
    unsigned dst_pitch = width * src_bpp;
    uchar*   dst_base  = dst;
    unsigned src_pitch = width * src_bpp;

    switch( src_bpp )
    {
        case 1:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

#pragma omp parallel for private(y)
            for ( x = 0; x < width; x++)
            {
                // work on column x in dst
                const unsigned index = x * src_bpp;
                uchar* dst_bits = dst_base + index;

                // scale each column
                for ( y = 0; y < dst_height; y++)
                {
                    const
                    uchar* src_base  = &src[ ( src_offset_y * width * src_bpp ) +
                                             ( src_offset_y * src_pitch + src_offset_x * src_bpp ) ];
                    // loop through column
                    const unsigned iLeft = weightsTable.getLeftBoundary(y);             // retrieve left boundary
                    const unsigned iLimit = weightsTable.getRightBoundary(y) - iLeft;   // retrieve right boundary
                    const uchar *src_bits = src_base + iLeft * src_pitch + index;
                    double v = 0;

                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(y, i);

                        if (useSCh)
                        {
                            double c = (weight * (double)src_bits[refSCh]);
                            v += c;
                        }
                        else
                        {
                            v += (weight * (double)src_bits[FI_GRAYA_GRAY]);
                        }
                        src_bits += src_pitch;
                    }

                    // clamp and place result in destination pixel
                    if (useSCh)
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int) (v + 0.5), 0, 0xFF);
                        dst_bits[FI_GRAYA_GRAY]   = cmpv;
                    }
                    else
                    {
                        dst_bits[FI_GRAYA_GRAY]   = (uchar)CLAMP<int>((int) (v + 0.5), 0, 0xFF);
                    }
                    dst_bits += dst_pitch;
                }
            }
        }
            break;

        case 2:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

#pragma omp parallel for private(y)
            for ( x = 0; x < width; x++)
            {
                // work on column x in dst
                const unsigned index = x * src_bpp;
                uchar *dst_bits = dst_base + index;

                // scale each column
                for ( y = 0; y < dst_height; y++)
                {
                    const
                    uchar* src_base  = &src[ ( src_offset_y * width * src_bpp ) +
                                             ( src_offset_y * src_pitch + src_offset_x * src_bpp ) ];
                    // loop through column
                    const unsigned iLeft = weightsTable.getLeftBoundary(y);
                    const unsigned iLimit = weightsTable.getRightBoundary(y) - iLeft;
                    const uchar *src_bits = src_base + iLeft * src_pitch + index;
                    double v = 0, a = 0;

                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(y, i);

                        if (useSCh)
                        {
                            double c = (weight * (double)src_bits[refSCh]);
                            v += c;
                            a += (weight * (double)src_bits[1]);
                        }
                        else
                        {
                            v += (weight * (double)src_bits[FI_GRAYA_GRAY]);
                            a += (weight * (double)src_bits[FI_GRAYA_ALPHA]);
                        }
                        src_bits += src_pitch;
                    }

                    // clamp and place result in destination pixel
                    if (useSCh)
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int) (v + 0.5), 0, 0xFF);
                        dst_bits[FI_GRAYA_GRAY]   = cmpv;
                        dst_bits[FI_GRAYA_ALPHA] = (uchar)CLAMP<int>((int) (a + 0.5), 0, 0xFF);
                    }
                    else
                    {
                        dst_bits[FI_GRAYA_GRAY]   = (uchar)CLAMP<int>((int) (v + 0.5), 0, 0xFF);
                        dst_bits[FI_GRAYA_ALPHA] = (uchar)CLAMP<int>((int) (a + 0.5), 0, 0xFF);
                    }
                    dst_bits += dst_pitch;
                }
            }
        }
            break;

        case 3:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(y)
            for ( x = 0; x < width; x++)
            {
                // work on column x in dst
                const unsigned index = x * src_bpp;
                uchar* dst_bits = dst_base + index;

                // scale each column
                for ( y = 0; y < dst_height; y++)
                {
                    const
                    uchar* src_base  = &src[ ( src_offset_y * width * src_bpp ) +
                                             ( src_offset_y * src_pitch + src_offset_x * src_bpp ) ];
                    // loop through column
                    const unsigned iLeft = weightsTable.getLeftBoundary(y);
                    const unsigned iLimit = weightsTable.getRightBoundary(y) - iLeft;
                    const uchar *src_bits = src_base + iLeft * src_pitch + index;
                    double r = 0, g = 0, b = 0;

                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(y, i);

                        if ( useSCh == true )
                        {
                            double c = (weight * (double)src_bits[refSCh]);
                            r += c;
                            g += c;
                            b += c;
                        }
                        else
                        {
                            r += (weight * (double)src_bits[FI_RGBA_RED]);
                            g += (weight * (double)src_bits[FI_RGBA_GREEN]);
                            b += (weight * (double)src_bits[FI_RGBA_BLUE]);
                        }
                        src_bits += src_pitch;
                    }

                    // clamp and place result in destination pixel
                    if ( useSCh == true )
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int) (r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_RED]   = cmpv;
                        dst_bits[FI_RGBA_GREEN] = cmpv;
                        dst_bits[FI_RGBA_BLUE]  = cmpv;
                    }
                    else
                    {
                        dst_bits[FI_RGBA_RED]   = (uchar)CLAMP<int>((int) (r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_GREEN] = (uchar)CLAMP<int>((int) (g + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_BLUE]  = (uchar)CLAMP<int>((int) (b + 0.5), 0, 0xFF);
                    }
                    dst_bits += dst_pitch;
                }
            }
        }
        break;

        case 4:
        {
            OMPSIZE_T x = 0;
            OMPSIZE_T y = 0;

            #pragma omp parallel for private(y)
            for ( x = 0; x < width; x++ )
            {
                // work on column x in dst
                const unsigned index = x * src_bpp;
                uchar *dst_bits = dst_base + index;

                // scale each column
                for ( y = 0; y < dst_height; y++)
                {
                    const
                    uchar* src_base  = &src[ ( src_offset_y * width * src_bpp ) +
                                             ( src_offset_y * src_pitch + src_offset_x * src_bpp ) ];
                    // loop through column
                    const unsigned iLeft = weightsTable.getLeftBoundary(y);
                    const unsigned iLimit = weightsTable.getRightBoundary(y) - iLeft;
                    const uchar *src_bits = src_base + iLeft * src_pitch + index;
                    double r = 0, g = 0, b = 0, a = 0;

                    for (unsigned i = 0; i <= iLimit; i++)
                    {
                        // scan between boundaries
                        // accumulate weighted effect of each neighboring pixel
                        const double weight = weightsTable.getWeight(y, i);

                        if ( useSCh == true )
                        {
                            double c = (weight * (double)src_bits[refSCh]);
                            r += c;
                            g += c;
                            b += c;
                            a += (weight * (double)src_bits[FI_RGBA_ALPHA]);
                        }
                        else
                        {
                            r += (weight * (double)src_bits[FI_RGBA_RED]);
                            g += (weight * (double)src_bits[FI_RGBA_GREEN]);
                            b += (weight * (double)src_bits[FI_RGBA_BLUE]);
                            a += (weight * (double)src_bits[FI_RGBA_ALPHA]);
                        }
                        src_bits += src_pitch;
                    }

                    // clamp and place result in destination pixel
                    if ( useSCh == true )
                    {
                        uchar cmpv = (uchar)CLAMP<int>((int) (r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_RED]   = cmpv;
                        dst_bits[FI_RGBA_GREEN] = cmpv;
                        dst_bits[FI_RGBA_BLUE]  = cmpv;
                        dst_bits[FI_RGBA_ALPHA] = (uchar)CLAMP<int>((int) (a + 0.5), 0, 0xFF);
                    }
                    else
                    {
                        dst_bits[FI_RGBA_RED]   = (uchar)CLAMP<int>((int) (r + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_GREEN] = (uchar)CLAMP<int>((int) (g + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_BLUE]  = (uchar)CLAMP<int>((int) (b + 0.5), 0, 0xFF);
                        dst_bits[FI_RGBA_ALPHA] = (uchar)CLAMP<int>((int) (a + 0.5), 0, 0xFF);
                    }
                    dst_bits += dst_pitch;
                }
            }
        }
        break;
    }
}
