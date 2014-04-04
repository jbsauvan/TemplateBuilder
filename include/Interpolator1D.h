/**
 *  @file  Interpolator1D.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    04/02/2014
 *
 *  @internal
 *     Created :  04/02/2014
 * Last update :  04/02/2014 05:15:35 PM
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

class Interpolator1D
{
    public:
        Interpolator1D();
        ~Interpolator1D();

        TH1D* interpolate(TH1D* rawHisto);


    private:
        double chi2(TH1D* rawHisto, TH1D* interpHisto);
        TH1D* histoFromGraph(TGraphAsymmErrors* graph, TH1D* templateHisto);
        TGraphAsymmErrors* graphFromHisto(TH1D* histo, TH1D* histoX);
        TGraphAsymmErrors* addPointToGraph(TGraphAsymmErrors* graph, double x, double y, double xdo, double xup, double ydo, double yup);
        void rebinHisto();
        void mergeLonelyBins();
        void constrainInterpolation();
        void correctLocalBias();

        bool m_linearInterp;
        double m_nSigmasRebin;
        double m_nSigmasLonely;
        TH1D* m_rawHisto;
        TH1D* m_rebinHisto;
        TH1D* m_rebinHistoX;
        TH1D* m_interpHisto;
        TGraphAsymmErrors* m_interpGraph;

        TFile* m_outputFile;

};

#endif
