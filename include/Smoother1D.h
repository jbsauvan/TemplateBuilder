/**
 *  @file  Smoother1D.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    04/12/2014
 *
 *  @internal
 *     Created :  04/12/2014
 * Last update :  04/12/2014 05:47:10 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H

#include <TH1D.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <algorithm>

class Smoother1D
{
    public:
        Smoother1D();
        ~Smoother1D();

        TH1D* smooth(TH1D* rawHisto);


    private:
        double chi2(TH1D* rawHisto, TH1D* interpHisto);
        TGraphAsymmErrors* graphFromHisto(TH1D* histo, TH1D* histoX);
        void rebinHisto();
        void mergeZero();
        void mergeLonelyBins();
        void computeSmoothingWidth();
        void computeSmoothHisto();

        double m_nSigmasRebin;
        double m_nSigmasLonely;
        TH1D* m_rawHisto;
        TH1D* m_rebinHisto;
        TH1D* m_rebinHistoX;
        TGraph* m_smoothWidth;
        TH1D* m_smoothHisto;


};

#endif
