/**
 *  @file  GaussKernelSmoother.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/13/2013
 *
 *  @internal
 *     Created :  12/13/2013
 * Last update :  12/13/2013 09:58:00 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef GAUSSKERNELSMOOTHER_H
#define GAUSSKERNELSMOOTHER_H

#include <TH1.h>
#include <algorithm>

class GaussKernelSmoother
{
    public:
        GaussKernelSmoother();
        GaussKernelSmoother(unsigned int ndim);
        ~GaussKernelSmoother();

        TH1* smooth(const TH1* histo);
        void setWidths(const std::vector<TH1*>& widths);
        void setWidthScalingFactor(double widthScalingFactor){m_widthScalingFactor = widthScalingFactor;}

    private:
        inline double weight(const std::vector<double>& x0, const std::vector<double>& xi);
        inline double weight(const std::vector<double>& x0, double dx, int axis);
        inline std::pair<double,double> smoothedValueError(const TH1* histo, const std::vector<double>& x0);
        inline std::pair<double,double> smoothed2DValueError(const TH1* histo, const std::vector<double>& x0);
        inline std::pair<double,double> smoothed3DValueError(const TH1* histo, const std::vector<double>& x0);

        unsigned int m_ndim;
        std::vector<TH1*> m_widths;
        double m_widthScalingFactor;

};


#endif
