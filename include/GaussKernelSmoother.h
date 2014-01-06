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

#include <TH2F.h>
#include <algorithm>

class GaussKernelSmoother
{
    public:
        GaussKernelSmoother();
        ~GaussKernelSmoother();

        TH2F* smooth(const TH2F* histo);
        void setWidths(TH2F* widthx, TH2F* widthy);

    private:
        TH2F* m_widthx;
        TH2F* m_widthy;
        double weight(double x0, double y0, double xi, double yi);
        std::pair<double,double> smoothedValueError(const TH2F* histo, double x0, double y0);


};


#endif
