/*
 * ctmf.c - Constant-time median filtering
 * Copyright (C) 2006  Simon Perreault
 *
 * Reference: S. Perreault and P. Hébert, "Median Filtering in Constant Time",
 * IEEE Transactions on Image Processing, September 2007.
 *
 * This program has been obtained from http://nomis80.org/ctmf.html. No patent
 * covers this program, although it is subject to the following license:
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact:
 *  Laboratoire de vision et systèmes numériques
 *  Pavillon Adrien-Pouliot
 *  Université Laval
 *  Sainte-Foy, Québec, Canada
 *  G1K 7P4
 *
 *  perreaul@gel.ulaval.ca
 */

/**
 * This structure represents a two-tier histogram. The first tier (known as the
 * "coarse" level) is 4 bit wide and the second tier (known as the "fine" level)
 * is 8 bit wide. Pixels inserted in the fine level also get inserted into the
 * coarse bucket designated by the 4 MSBs of the fine bucket value.
 *
 * The structure is aligned on 16 bytes, which is a prerequisite for SIMD
 * instructions. Each bucket is 16 bit wide, which means that extra care must be
 * taken to prevent overflow.
 */
class Histogram {

  int[] coarse;
  int[][] fine;
    
  Histogram() {
    coarse = new int[16];
    fine = new int[16][16];
  }

}

#define HOP(h,x,op) \
    h.coarse[x>>4] op; \
    *((uint16_t*) h.fine + x) op;
  int[] coarse;
  int[][] fine;

#define COP(c,j,x,op) \
    h_coarse[ 16*(n*c+j) + (x>>4) ] op; \
    h_fine[ 16 * (n*(16*c+(x>>4)) + j) + (x & 0xF) ] op;
    
/**
 * Adds histograms \a x and \a y and stores the result in \a y. 
 */
void histogram_add(int[] x[16], int[] y[16]) {
    for (int i = 0; i < 16; ++i) {
        y[i] += x[i];
    }
}

/**
 * Subtracts histogram \a x from \a y and stores the result in \a y.
 */
void histogram_sub(int[] x[16], int[] y[16]) {
    for (int i = 0; i < 16; ++i) {
        y[i] -= x[i];
    }
}

void histogram_muladd(int[] a, int[] x[16], int[] y[16]) {
    for (int i = 0; i < 16; ++i) {
        y[i] += a * x[i];
    }
}

void ctmf_helper(unsigned char* src, unsigned char* dst, int width, int height, int src_step, int dst_step, int r, int cn, int pad_left, int pad_right) {
    int m = height, n = width;
    int i, j, k, c;
    unsigned char *p, *q;

    Histogram H[4];
    int[] *h_coarse, *h_fine, luc[4][16];

    assert(src);
    assert(dst);
    assert(r >= 0);
    assert(width >= 2*r+1);
    assert(height >= 2*r+1);
    assert(src_step != 0);
    assert(dst_step != 0);

    h_coarse = (int[]*) calloc( 1 * 16 * n * cn, sizeof(int[]));
    h_fine   = (int[]*) calloc(16 * 16 * n * cn, sizeof(int[]));

    /* First row initialization */
    for (j = 0; j < n; ++j) {
        for (c = 0; c < cn; ++c) {
            COP(c, j, src[cn*j+c], += r+1);
        }
    }
    
    for (i = 0; i < r; ++i) {
        for (j = 0; j < n; ++j) {
            for (c = 0; c < cn; ++c) {
                COP(c, j, src[src_step*i+cn*j+c], ++);
            }
        }
    }

    for (i = 0; i < m; ++i) {
        /* Update column histograms for entire row. */
        p = src + src_step * MAX(0, i-r-1);
        q = p + cn * n;
        
        for (j = 0; p != q; ++j) {
            for (c = 0; c < cn; ++c, ++p) {
                COP(c, j, *p, --);
            }
        }

        p = src + src_step * MIN(m-1, i+r);
        q = p + cn * n;
        
        for (j = 0; p != q; ++j) {
            for (c = 0; c < cn; ++c, ++p) {
                COP(c, j, *p, ++);
            }
        }

        /* First column initialization */
        memset(H, 0, cn*sizeof(H[0]));
        memset(luc, 0, cn*sizeof(luc[0]));
        if (pad_left) {
            for (c = 0; c < cn; ++c) {
                histogram_muladd(r, &h_coarse[16*n*c], H[c].coarse);
            }
        }
        
        for (j = 0; j < (pad_left ? r : 2*r); ++j) {
            for (c = 0; c < cn; ++c) {
                histogram_add(&h_coarse[16*(n*c+j)], H[c].coarse);
            }
        }
        
        for (c = 0; c < cn; ++c) {
            for (k = 0; k < 16; ++k) {
                histogram_muladd(2*r+1, &h_fine[16*n*(16*c+k)], &H[c].fine[k][0]);
            }
        }

        for (j = pad_left ? 0 : r; j < (pad_right ? n : n-r); ++j) {
            for (c = 0; c < cn; ++c) {
                int[] t = 2*r*r + 2*r;
                int[] sum = 0, *segment;
                int b;

                histogram_add(&h_coarse[16*(n*c + MIN(j+r,n-1))], H[c].coarse);

                /* Find median at coarse level */
                for (k = 0; k < 16 ; ++k) {
                    sum += H[c].coarse[k];
                    if (sum > t) {
                        sum -= H[c].coarse[k];
                        break;
                    }
                }
                assert(k < 16);

                /* Update corresponding histogram segment */
                if (luc[c][k] <= j-r) {
                    memset(&H[c].fine[k], 0, 16 * sizeof(int[]));
                    for (luc[c][k] = j-r; luc[c][k] < MIN(j+r+1,n); ++luc[c][k]) {
                        histogram_add(&h_fine[16*(n*(16*c+k)+luc[c][k])], H[c].fine[k]);
                    }
                    if (luc[c][k] < j+r+1) {
                        histogram_muladd(j+r+1 - n, &h_fine[16*(n*(16*c+k)+(n-1))], &H[c].fine[k][0]);
                        luc[c][k] = j+r+1;
                    }
                } else {
                    for (; luc[c][k] < j+r+1; ++luc[c][k]) {
                        histogram_sub(&h_fine[16*(n*(16*c+k)+MAX(luc[c][k]-2*r-1,0))], H[c].fine[k]);
                        histogram_add(&h_fine[16*(n*(16*c+k)+MIN(luc[c][k],n-1))], H[c].fine[k]);
                    }
                }

                histogram_sub(&h_coarse[16*(n*c+MAX(j-r,0))], H[c].coarse);

                /* Find median in segment */
                segment = H[c].fine[k];
                
                for (b = 0; b < 16 ; ++b) {
                    sum += segment[b];
                    if (sum > t) {
                        dst[dst_step*i+cn*j+c] = 16*k + b;
                        break;
                    }
                }
                
                assert(b < 16);
            }
        }
    }

  free(h_coarse);
  free(h_fine);
}

/**
 * \brief Constant-time median filtering
 *
 * This function does a median filtering of an 8-bit image. The source image is
 * processed as if it was padded with zeros. The median kernel is square with
 * odd dimensions. Images of arbitrary size may be processed.
 *
 * To process multi-channel images, you must call this function multiple times,
 * changing the source and destination adresses and steps such that each channel
 * is processed as an independent single-channel image.
 *
 * Processing images of arbitrary bit depth is not supported.
 *
 * The computing time is O(1) per pixel, independent of the radius of the
 * filter. The algorithm's initialization is O(r*width), but it is negligible.
 * Memory usage is simple: it will be as big as the cache size, or smaller if
 * the image is small. For efficiency, the histograms' bins are 16-bit wide.
 * This may become too small and lead to overflow as \a r increases.
 *
 * \param src           Source image data.
 * \param dst           Destination image data. Must be preallocated.
 * \param width         Image width, in pixels.
 * \param height        Image height, in pixels.
 * \param src_step      Distance between adjacent pixels on the same column in
 *                      the source image, in bytes.
 * \param dst_step      Distance between adjacent pixels on the same column in
 *                      the destination image, in bytes.
 * \param r             Median filter radius. The kernel will be a 2*r+1 by
 *                      2*r+1 square.
 * \param cn            Number of channels. For example, a grayscale image would
 *                      have cn=1 while an RGB image would have cn=3.
 * \param memsize       Maximum amount of memory to use, in bytes. Set this to
 *                      the size of the L2 cache, then vary it slightly and
 *                      measure the processing time to find the optimal value.
 *                      For example, a 512 kB L2 cache would have
 *                      memsize=512*1024 initially.
 */
void ctmf(unsigned char* src, unsigned char* dst, int width, int height, int src_step, int dst_step, int r, int cn, long unsigned int memsize) {
    /*
     * Processing the image in vertical stripes is an optimization made
     * necessary by the limited size of the CPU cache. Each histogram is 544
     * bytes big and therefore I can fit a limited number of them in the cache.
     * That number may sometimes be smaller than the image width, which would be
     * the number of histograms I would need without stripes.
     *
     * I need to keep histograms in the cache so that they are available
     * quickly when processing a new row. Each row needs access to the previous
     * row's histograms. If there are too many histograms to fit in the cache,
     * thrashing to RAM happens.
     *
     * To solve this problem, I figure out the maximum number of histograms
     * that can fit in cache. From this is determined the number of stripes in
     * an image. The formulas below make the stripes all the same size and use
     * as few stripes as possible.
     *
     * Note that each stripe causes an overlap on the neighboring stripes, as
     * when mowing the lawn. That overlap is proportional to r. When the overlap
     * is a significant size in comparison with the stripe size, then we are not
     * O(1) anymore, but O(r). In fact, we have been O(r) all along, but the
     * initialization term was neglected, as it has been (and rightly so) in B.
     * Weiss, "Fast Median and Bilateral Filtering", SIGGRAPH, 2006. Processing
     * by stripes only makes that initialization term bigger.
     *
     * Also, note that the leftmost and rightmost stripes don't need overlap.
     * A flag is passed to ctmf_helper() so that it treats these cases as if the
     * image was zero-padded.
     */
    int stripes = (int) ceil((double) (width - 2*r) / (memsize / sizeof(Histogram) - 2*r));
    int stripe_size = (int) ceil((double) (width + stripes*2*r - 2*r) / stripes);

    for (int i = 0; i < width; i += stripe_size - 2*r) {
        int stripe = stripe_size;
        /* Make sure that the filter kernel fits into one stripe. */
        if (i + stripe_size - 2*r >= width || width - (i + stripe_size - 2*r) < 2*r+1) {
            stripe = width - i;
        }

        ctmf_helper(src + cn*i, dst + cn*i, stripe, height, src_step, dst_step, r, cn, i == 0, stripe == width - i);

        if (stripe == width - i) {
            break;
        }
    }
}

}
