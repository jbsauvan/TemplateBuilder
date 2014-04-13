/**
 *  @file  Smoother1D.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    04/12/2014
 *
 *  @internal
 *     Created :  04/12/2014
 * Last update :  04/12/2014 05:38:26 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "Smoother1D.h"

#include <TCanvas.h>
#include <TMath.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

struct pair_comp
{
   bool operator()( const pair<double, int>& lhs, const double& rhs ) const 
   { 
        return lhs.first < rhs;
   }
   bool operator()( const double& lhs, const pair<double, int>& rhs ) const 
   {
        return lhs < rhs.first;
   }
   bool operator()( const pair<double, int>& lhs, const pair<double, int>& rhs ) const 
   {
        return lhs.first < rhs.first;
   }
};


/*****************************************************************/
Smoother1D::Smoother1D():
    m_nSigmasRebin(2.),
    m_nSigmasLonely(5.),
    m_rawHisto(NULL),
    m_rebinHisto(NULL),
    m_rebinHistoX(NULL),
    m_smoothWidth(NULL),
    m_smoothHisto(NULL)
/*****************************************************************/
{
}

/*****************************************************************/
Smoother1D::~Smoother1D()
/*****************************************************************/
{
    // raw histo and smooth histos are used outside this class
    // So do not delete them here
    if(m_rebinHisto) m_rebinHisto->Delete();
    if(m_rebinHistoX) m_rebinHistoX->Delete();
    if(m_smoothWidth) m_smoothWidth->Delete();
}


/*****************************************************************/
double Smoother1D::chi2(TH1D* rawHisto, TH1D* interpHisto)
/*****************************************************************/
{
    double c = 0.;
    int nbins = rawHisto->GetNbinsX();
    for(int b=1;b<=nbins;b++)
    {
        double y1 = rawHisto->GetBinContent(b);
        double y2 = interpHisto->GetBinContent(b);
        double e  = rawHisto->GetBinError(b);
        if(e>0)
        {
            c += (y2-y1)*(y2-y1)/(e*e);
        }
    }
    return c;
}


/*****************************************************************/
TGraphAsymmErrors* Smoother1D::graphFromHisto(TH1D* histo, TH1D* histoX)
/*****************************************************************/
{
    int nbins = histoX->GetNbinsX();
    double* xs   = new double[nbins];
    double* ys   = new double[nbins];
    double* xups = new double[nbins];
    double* xdos = new double[nbins];
    double* yups = new double[nbins];
    double* ydos = new double[nbins];
    for(int b=1;b<=nbins;b++)
    {
        xs[b-1]   = histoX->GetBinContent(b);
        xups[b-1] = histoX->GetXaxis()->GetBinUpEdge(b) - xs[b-1];
        xdos[b-1] = xs[b-1] - histoX->GetXaxis()->GetBinLowEdge(b);
        ys[b-1]   = histo->GetBinContent(b);
        ydos[b-1] = histo->GetBinError(b);
        yups[b-1] = histo->GetBinError(b);
    }
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(nbins, xs, ys, xdos,xups, ydos,yups);
    stringstream gName;
    gName << m_rebinHisto->GetName() << "_g";
    graph->SetName(gName.str().c_str());
    return graph;
}



/*****************************************************************/
void Smoother1D::rebinHisto()
/*****************************************************************/
{
    //cerr<<"Rebinning histo "<<m_rawHisto->GetName()<<"\n";
    int nBins = m_rawHisto->GetNbinsX();
    bool done = false;
    stringstream name, namex;
    name << m_rawHisto->GetName() << "_r";
    namex << m_rawHisto->GetName() << "_xr";
    TH1D* histoRebin = (TH1D*)m_rawHisto->Clone(name.str().c_str());
    TH1D* histoXpos  = (TH1D*)m_rawHisto->Clone(namex.str().c_str());
    // FIXME: this reference binWidth works only with fixed size bins
    double binWidth = m_rawHisto->GetXaxis()->GetBinWidth(1);



    for(int b=1; b<=nBins; b++)
    {
        double x = histoRebin->GetXaxis()->GetBinCenter(b);
        histoXpos->SetBinContent(b, x);
    }

    for(int b=1; b<nBins; b++)
    {
        double value1 = histoRebin->GetBinContent(b);
        double error1 = histoRebin->GetBinError(b);
        double relError1 = (value1>0 ? error1/value1 : 0);
        double value2 = histoRebin->GetBinContent(b+1);
        double error2 = histoRebin->GetBinError(b+1);
        double relError2 = (value2>0 ? error1/value2 : 0);
        double errorRatio = 0;
        if(value1>value2) errorRatio = (error2>0 ? error1/error2: 0.);
        else errorRatio = (error1>0 ? error2/error1: 0.);
        if(value1>0 && value2>0 && error1>0 && error2>0 && (relError1>0.1 || relError2>0.1))
        {
            double ratio = error1/error2*sqrt(value2/value1);
            if(ratio<0.5)
            {
                histoRebin->SetBinError(b, error1/ratio);
            }
            else if(ratio>2)
            {
                histoRebin->SetBinError(b+1, error2*ratio);
            }
        }
    }



    int iteration = 1;
    // Iterative process. Two consecutive bins are merged together at each step.
    // Bins with the less significative differences are merged first. This is redone until all the differences are significant enough.
    while(!done && nBins>1)
    {
        int binDown = 0;
        int binUp = 0;
        double significanceMin = 9.e10;
        double maxi = 0.;
        double binMaxi = 0;
        for(int b=1; b<=nBins; b++)
        {
            double value = histoRebin->GetBinContent(b);
            if(value>maxi)
            {
                maxi = value;
                binMaxi = b;
            }
        }
        // First choose two bins to be merged
        for(int b=1; b<=nBins-1; b++)
        {
            double valueDown = histoRebin->GetBinContent(b);
            double valueUp   = histoRebin->GetBinContent(b+1);
            double errorDown = histoRebin->GetBinError(b);
            double errorUp   = histoRebin->GetBinError(b+1);
            double relWidth  = (histoRebin->GetXaxis()->GetBinWidth(b)+histoRebin->GetXaxis()->GetBinWidth(b+1))/binWidth;

            if(valueUp>0 && valueDown>0 && (errorDown>0 || errorUp>0))
            {
                double diff      = fabs(valueUp-valueDown);
                double errorDiff = sqrt(errorUp*errorUp+errorDown*errorDown);

                double sig = diff/errorDiff;
                //if(sig<m_nSigmasRebin && sig*sqrt(relWidth-1.)<significanceMin) // favor merging of smaller bins
                double sigCutoff = m_nSigmasRebin;
                if(b==binMaxi || b+1==binMaxi) sigCutoff = m_nSigmasRebin/2.; //lower cutoff in the case of the maximum
                if(sig<sigCutoff && sig*(relWidth-1.)<significanceMin)
                {
                    //significanceMin = sig*sqrt(relWidth-1.);
                    significanceMin = sig*(relWidth-1.);
                    binUp = b+1;
                    binDown = b;
                }
            }
        }

        // then merge the bins (the weighted content average is taken as new bin value)
        // the weighted position is also computed
        if(binDown!=0)
        {
            double xDown     = histoXpos->GetXaxis()->GetBinCenter(binDown);
            double xUp       = histoXpos->GetXaxis()->GetBinCenter(binUp);
            double valueUp   = histoRebin->GetBinContent(binUp); 
            double errorUp   = histoRebin->GetBinError(binUp);
            double valueDown = histoRebin->GetBinContent(binDown); 
            double errorDown = histoRebin->GetBinError(binDown);

            if(errorDown==0 && errorUp>0) errorDown = errorUp;
            if(errorDown>0 && errorUp==0) errorUp = errorDown;

            double meanXpos  = (xDown/(errorDown*errorDown) + xUp/(errorUp*errorUp))/(1./(errorDown*errorDown) + 1./(errorUp*errorUp));
            double meanValue = (valueDown/(errorDown*errorDown) + valueUp/(errorUp*errorUp))/(1./(errorDown*errorDown) + 1./(errorUp*errorUp));
            double meanError = sqrt(1./(1./(errorUp*errorUp) + 1./(errorDown*errorDown)));
            histoRebin->SetBinContent(binDown, meanValue);
            histoRebin->SetBinError(binDown, meanError);
            histoRebin->SetBinContent(binUp, meanValue);
            histoRebin->SetBinError(binUp, meanError);
            histoXpos->SetBinContent(binDown, meanXpos);
            histoXpos->SetBinError(binDown, 0.);
            histoXpos->SetBinContent(binUp, meanXpos);
            histoXpos->SetBinError(binUp, 0.);

            //rebin
            //cerr<<" [";
            double* newBins = new double[nBins];
            int b = 0;
            for(int i=1;i<=nBins;i++)
            {
                newBins[b] = histoRebin->GetXaxis()->GetBinLowEdge(i);
                //cerr<<newBins[b]<<",";
                if(i==binDown) i++;
                b++;
            }
            newBins[nBins-1] = histoRebin->GetXaxis()->GetBinUpEdge(nBins);
            //cerr<<newBins[nBins-1]<<"]\n";

            stringstream nametmp, namextmp;
            nametmp << histoRebin->GetName() << "_tmp";
            namextmp << histoXpos->GetName() << "_tmp";
            TH1D* histoRebinTemp = new TH1D(nametmp.str().c_str(), "tmpHisto", nBins-1, newBins);
            TH1D* histoXposTemp  = new TH1D(namextmp.str().c_str(), "tmpHisto", nBins-1, newBins);
            b = 1;
            for(int i=1;i<=nBins;i++)
            {
                histoRebinTemp->SetBinContent(b, histoRebin->GetBinContent(i));
                histoRebinTemp->SetBinError(b, histoRebin->GetBinError(i));
                histoXposTemp->SetBinContent(b, histoXpos->GetBinContent(i));
                histoXposTemp->SetBinError(b, histoXpos->GetBinError(i));
                if(i==binDown) i++;
                b++;
            }
            histoRebin->Delete();
            histoRebin = histoRebinTemp;
            histoRebin->SetName(name.str().c_str());
            histoRebin->SetTitle(name.str().c_str());
            histoXpos->Delete();
            histoXpos = histoXposTemp;
            histoXpos->SetName(namex.str().c_str());
            histoXpos->SetTitle(namex.str().c_str());

            nBins = histoRebin->GetNbinsX();
            iteration++;
        }
        else
            done = true;
    }
    if(m_rebinHisto) 
    {
        m_rebinHisto->Delete();
        m_rebinHisto = NULL;
    }
    if(m_rebinHistoX) 
    {
        m_rebinHistoX->Delete();
        m_rebinHistoX = NULL;
    }
    m_rebinHisto = histoRebin;
    m_rebinHistoX = histoXpos;
}

/*****************************************************************/
void Smoother1D::mergeZero()
/*****************************************************************/
{
    int nBins = m_rebinHisto->GetNbinsX();
    bool done = false;
    stringstream name, namex;
    name << m_rebinHisto->GetName() << "_r";
    namex << m_rebinHistoX->GetName() << "_xr";
    TH1D* histoRebin = (TH1D*)m_rebinHisto->Clone(name.str().c_str());
    TH1D* histoXpos  = (TH1D*)m_rebinHistoX->Clone(namex.str().c_str());



    int iteration = 1;
    while(!done && nBins>1)
    {
        int binDown = 0;
        int binUp = 0;
        double significanceMin = 9.e10;
        // First choose two bins to be merged
        for(int b=1; b<=nBins-1; b++)
        {
            double valueDown = histoRebin->GetBinContent(b);
            double valueUp   = histoRebin->GetBinContent(b+1);

            if(valueUp==0 && valueDown==0)
            {
                binUp = b+1;
                binDown = b;
            }
        }

        if(binDown!=0)
        {
            double xDown     = histoXpos->GetXaxis()->GetBinCenter(binDown);
            double xUp       = histoXpos->GetXaxis()->GetBinCenter(binUp);
            double valueUp   = histoRebin->GetBinContent(binUp); 
            double errorUp   = histoRebin->GetBinError(binUp);
            double valueDown = histoRebin->GetBinContent(binDown); 
            double errorDown = histoRebin->GetBinError(binDown);


            double meanXpos  = (xUp+xDown)/2.;
            double meanValue = 0.;
            double meanError = 0.;
            histoRebin->SetBinContent(binDown, meanValue);
            histoRebin->SetBinError(binDown, meanError);
            histoRebin->SetBinContent(binUp, meanValue);
            histoRebin->SetBinError(binUp, meanError);
            histoXpos->SetBinContent(binDown, meanXpos);
            histoXpos->SetBinError(binDown, 0.);
            histoXpos->SetBinContent(binUp, meanXpos);
            histoXpos->SetBinError(binUp, 0.);

            //rebin
            //cerr<<" [";
            double* newBins = new double[nBins];
            int b = 0;
            for(int i=1;i<=nBins;i++)
            {
                newBins[b] = histoRebin->GetXaxis()->GetBinLowEdge(i);
                //cerr<<newBins[b]<<",";
                if(i==binDown) i++;
                b++;
            }
            newBins[nBins-1] = histoRebin->GetXaxis()->GetBinUpEdge(nBins);
            //cerr<<newBins[nBins-1]<<"]\n";

            stringstream nametmp, namextmp;
            nametmp << histoRebin->GetName() << "_tmp";
            namextmp << histoXpos->GetName() << "_tmp";
            TH1D* histoRebinTemp = new TH1D(nametmp.str().c_str(), "tmpHisto", nBins-1, newBins);
            TH1D* histoXposTemp  = new TH1D(namextmp.str().c_str(), "tmpHisto", nBins-1, newBins);
            b = 1;
            for(int i=1;i<=nBins;i++)
            {
                histoRebinTemp->SetBinContent(b, histoRebin->GetBinContent(i));
                histoRebinTemp->SetBinError(b, histoRebin->GetBinError(i));
                histoXposTemp->SetBinContent(b, histoXpos->GetBinContent(i));
                histoXposTemp->SetBinError(b, histoXpos->GetBinError(i));
                if(i==binDown) i++;
                b++;
            }
            histoRebin->Delete();
            histoRebin = histoRebinTemp;
            histoRebin->SetName(name.str().c_str());
            histoRebin->SetTitle(name.str().c_str());
            histoXpos->Delete();
            histoXpos = histoXposTemp;
            histoXpos->SetName(namex.str().c_str());
            histoXpos->SetTitle(namex.str().c_str());

            nBins = histoRebin->GetNbinsX();
            iteration++;
        }
        else
            done = true;
    }
    if(m_rebinHisto) 
    {
        m_rebinHisto->Delete();
        m_rebinHisto = NULL;
    }
    if(m_rebinHistoX) 
    {
        m_rebinHistoX->Delete();
        m_rebinHistoX = NULL;
    }
    m_rebinHisto = histoRebin;
    m_rebinHistoX = histoXpos;
}

/*****************************************************************/
void Smoother1D::mergeLonelyBins()
/*****************************************************************/
{
    //cerr<<"Merge lonely bins for "<<m_rawHisto->GetName()<<"\n";
    int nBins = m_rebinHisto->GetNbinsX();
    bool done = false;
    stringstream name, namex;
    name << m_rebinHisto->GetName() << "_r";
    namex << m_rebinHistoX->GetName() << "_xr";
    TH1D* histoRebin = (TH1D*)m_rebinHisto->Clone(name.str().c_str());
    TH1D* histoXpos  = (TH1D*)m_rebinHistoX->Clone(namex.str().c_str());

    int iteration = 1;
    while(!done && nBins>1)
    {
        int binDown = 0;
        int binUp   = 0;
        int bin     = 0;
        double significanceMin = 9.e10;
        for(int b=2; b<=nBins-1; b++)
        {
            double valueDown    = histoRebin->GetBinContent(b-1);
            double value        = histoRebin->GetBinContent(b);
            double valueUp      = histoRebin->GetBinContent(b+1);
            double errorDown    = histoRebin->GetBinError(b-1);
            double error        = histoRebin->GetBinError(b);
            double errorUp      = histoRebin->GetBinError(b+1);
            double binWidthDown = histoRebin->GetXaxis()->GetBinWidth(b-1);
            double binWidth     = histoRebin->GetXaxis()->GetBinWidth(b);
            double binWidthUp   = histoRebin->GetXaxis()->GetBinWidth(b+1);
            if(errorDown>0 || errorUp>0 || error>0)
            {
                double diffDown        = value-valueDown;
                double diffUp          = valueUp-value;
                double doubleDiff      = fabs(diffDown)+fabs(diffUp);
                double errorDoubleDiff = sqrt(errorUp*errorUp + errorDown*errorDown + 4.*error*error);

                double sig = doubleDiff/errorDoubleDiff;
                //if(sig<m_nSigmasLonely && 3.*binWidth<=binWidthUp && 3.*binWidth<=binWidthDown && sig<significanceMin)
                bool candidate = false;
                //if(value>0 && 6.*binWidth<=(binWidthUp+binWidthDown)) candidate = true;
                if(3.*binWidth<=binWidthDown && 3.*binWidth<=binWidthUp) candidate = true;
                //if(4.*binWidth<=binWidthDown && 2.*binWidth<=binWidthUp) candidate = true;
                //if(2.*binWidth<=binWidthDown && 4.*binWidth<=binWidthUp) candidate = true;
                if(sig<m_nSigmasLonely && candidate && sig<significanceMin)
                {
                    //cerr<<" sig="<<sig<<",relWidthUp="<<binWidth/binWidthUp<<",relWidthDo="<<binWidth/binWidthDown<<"\n"; 
                    binUp = b+1;
                    binDown = b-1;
                    bin = b;
                    significanceMin = sig;
                }
            }
        }

        if(bin!=0)
        {
            double xUp       = histoXpos->GetXaxis()->GetBinCenter(binUp);
            double valueUp   = histoRebin->GetBinContent(binUp); 
            double errorUp   = histoRebin->GetBinError(binUp);
            double xDown     = histoXpos->GetXaxis()->GetBinCenter(binDown);
            double valueDown = histoRebin->GetBinContent(binDown); 
            double errorDown = histoRebin->GetBinError(binDown);
            double x         = histoXpos->GetXaxis()->GetBinCenter(bin);
            double value     = histoRebin->GetBinContent(bin); 
            double error     = histoRebin->GetBinError(bin);


            if(errorDown==0 ) errorDown = max(errorUp,error);
            if(errorUp==0)    errorUp   = max(errorDown,error);
            if(error==0)      error     = max(errorDown,errorUp);

            double meanXpos  = (xDown/(errorDown*errorDown) + x/(error*error) + xUp/(errorUp*errorUp))/(1./(errorDown*errorDown) + 1./(error*error) + 1./(errorUp*errorUp));
            double meanValue = (valueDown/(errorDown*errorDown) + value/(error*error) + valueUp/(errorUp*errorUp))/(1./(errorDown*errorDown) + 1./(error*error) + 1./(errorUp*errorUp));
            double meanError = sqrt(1./(1./(errorUp*errorUp) + 1./(error*error) +1./(errorDown*errorDown)));

            histoRebin->SetBinContent(binDown, meanValue);
            histoRebin->SetBinError(binDown, meanError);
            histoRebin->SetBinContent(binUp, meanValue);
            histoRebin->SetBinError(binUp, meanError);
            histoRebin->SetBinContent(bin, meanValue);
            histoRebin->SetBinError(bin, meanError);
            histoXpos->SetBinContent(binDown, meanXpos);
            histoXpos->SetBinError(binDown, 0.);
            histoXpos->SetBinContent(binUp, meanXpos);
            histoXpos->SetBinError(binUp, 0.);
            histoXpos->SetBinContent(bin, meanXpos);
            histoXpos->SetBinError(bin, 0.);

            //rebin
            //cerr<<"[";
            double* newBins = new double[nBins-1];
            int b = 0;
            for(int i=1;i<=nBins;i++)
            {
                newBins[b] = histoRebin->GetXaxis()->GetBinLowEdge(i);
                //cerr<<newBins[b]<<",";
                if(i==binDown) i+=2;
                b++;
            }
            newBins[nBins-2] = histoRebin->GetXaxis()->GetBinUpEdge(nBins);
            //cerr<<newBins[nBins-2]<<"]\n";

            stringstream nametmp, namextmp;
            nametmp << histoRebin->GetName() << "_tmp";
            namextmp << histoXpos->GetName() << "_tmp";
            TH1D* histoRebinTemp = new TH1D(nametmp.str().c_str(), "tmpHisto", nBins-2, newBins);
            TH1D* histoXposTemp  = new TH1D(namextmp.str().c_str(), "tmpHisto", nBins-2, newBins);
            b = 1;
            for(int i=1;i<=nBins;i++)
            {
                histoRebinTemp->SetBinContent(b, histoRebin->GetBinContent(i));
                histoRebinTemp->SetBinError(b, histoRebin->GetBinError(i));
                histoXposTemp->SetBinContent(b, histoXpos->GetBinContent(i));
                histoXposTemp->SetBinError(b, histoXpos->GetBinError(i));
                if(i==binDown) i+=2;
                b++;
            }
            histoRebin->Delete();
            histoRebin = histoRebinTemp;
            histoRebin->SetName(name.str().c_str());
            histoRebin->SetTitle(name.str().c_str());
            histoXpos->Delete();
            histoXpos = histoXposTemp;
            histoXpos->SetName(namex.str().c_str());
            histoXpos->SetTitle(namex.str().c_str());

            nBins = histoRebin->GetNbinsX();
            iteration++;
        }
        else
            done = true;
    }
    if(m_rebinHisto) 
    {
        m_rebinHisto->Delete();
        m_rebinHisto = NULL;
    }
    if(m_rebinHistoX) 
    {
        m_rebinHistoX->Delete();
        m_rebinHistoX = NULL;
    }
    m_rebinHisto = histoRebin;
    m_rebinHistoX = histoXpos;

}


/*****************************************************************/
void Smoother1D::computeSmoothingWidth()
/*****************************************************************/
{
    TH1D* binWidths = (TH1D*)m_rebinHisto->Clone("binWidths");
    int nbins = binWidths->GetNbinsX();
    for(int b=1;b<=nbins;b++)
    {
        double width = m_rebinHisto->GetXaxis()->GetBinWidth(b);
        binWidths->SetBinContent(b, width/4.);
        binWidths->SetBinError(b, 0.);
    }
    m_smoothWidth = graphFromHisto(binWidths, m_rebinHistoX);
}

/*****************************************************************/
void Smoother1D::computeSmoothHisto()
/*****************************************************************/
{
    stringstream name;
    name << m_rawHisto->GetName() << "_s";
    m_smoothHisto = (TH1D*)m_rawHisto->Clone(name.str().c_str());
    int nbins = m_rawHisto->GetNbinsX();
    for(int b=1;b<=nbins;b++)
    {
        double x = m_rawHisto->GetXaxis()->GetBinCenter(b);
        double width = m_smoothWidth->Eval(x);
        double sumw=0, sumwy=0;
        for (int i=1; i<=nbins;++i) 
        {
            double xi = m_rawHisto->GetXaxis()->GetBinCenter(i);
            double yi = m_rawHisto->GetBinContent(i);
            double ei = m_rawHisto->GetBinError(i);
            double dx = (x-xi)/width;
            double wi = TMath::Gaus(dx);
            //if (ei>0.) wi *=1./(ei*ei);// weight with error squared
            sumw += wi; sumwy += wi*yi;
        }
        m_smoothHisto->SetBinContent(b, sumwy/sumw);
    }
}


/*****************************************************************/
TH1D* Smoother1D::smooth(TH1D* rawHisto)
/*****************************************************************/
{

    if(m_rawHisto) m_rawHisto->Delete();
    if(m_rebinHisto) m_rebinHisto->Delete();
    if(m_rebinHistoX) m_rebinHistoX->Delete();
    if(m_smoothWidth) m_smoothWidth->Delete();
    if(m_smoothHisto) m_smoothHisto->Delete();
    m_rawHisto = NULL;
    m_rebinHisto = NULL;
    m_rebinHistoX = NULL;
    m_smoothWidth = NULL;
    m_smoothHisto = NULL;


    m_rawHisto = rawHisto;
    rebinHisto();
    mergeZero();
    mergeLonelyBins();
    computeSmoothingWidth();
    computeSmoothHisto();


    double chi2Histo = chi2(m_rawHisto, m_smoothHisto);
    int ndf = m_rawHisto->GetNbinsX()-m_rebinHisto->GetNbinsX();
    double chi2oNdf = (ndf>0 ? chi2Histo/(double)ndf : chi2Histo);
    cout<<"[INFO]     chi2/ndf(ref-raw) = "<<chi2Histo<<"/"<<ndf<<" = "<<chi2oNdf<<"\n";
    if(chi2oNdf>5)
    {
        cout<<"[WARN]     The reference histogram used for reweighting seems not good. Please check the produced control plot in the output file.\n";
    }



    return m_smoothHisto;
}
