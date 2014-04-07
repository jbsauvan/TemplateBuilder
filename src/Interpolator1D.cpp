/**
 *  @file  Interpolator1D.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    04/02/2014
 *
 *  @internal
 *     Created :  04/02/2014
 * Last update :  04/02/2014 05:20:00 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "Interpolator1D.h"

#include <TCanvas.h>

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
Interpolator1D::Interpolator1D():
    m_linearInterp(false),
    m_nSigmasRebin(2.),
    m_nSigmasLonely(5.),
    m_rawHisto(NULL),
    m_rebinHisto(NULL),
    m_rebinHistoX(NULL),
    m_interpHisto(NULL),
    m_interpGraph(NULL)
/*****************************************************************/
{
}

/*****************************************************************/
Interpolator1D::~Interpolator1D()
/*****************************************************************/
{
    // raw histo and interpolated histos are used outside this class
    // So do not delete them here
    if(m_rebinHisto) m_rebinHisto->Delete();
    if(m_rebinHistoX) m_rebinHistoX->Delete();
    if(m_interpGraph) m_interpGraph->Delete();
}


/*****************************************************************/
double Interpolator1D::chi2(TH1D* rawHisto, TH1D* interpHisto)
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
TH1D* Interpolator1D::histoFromGraph(TGraphAsymmErrors* graph, TH1D* templateHisto)
/*****************************************************************/
{
    stringstream name;
    name << graph->GetName() << "_h";
    TH1D* histo = (TH1D*)templateHisto->Clone(name.str().c_str());
    int nbins = histo->GetNbinsX();
    for(int b=1;b<=nbins;b++)
    {
        double x = histo->GetXaxis()->GetBinCenter(b);
        double value = 0.;
        if(m_linearInterp)
        {
            value = graph->Eval(x); // linear interpolation
        }
        else
        {
            value = graph->Eval(x,0, "S"); // Cubic splines interpolation
            if(!std::isfinite(value)) // Linear interpolation backup in case of crazy output
            {
                value = graph->Eval(x);
            }
        }
        //cerr<<" interp histo value ("<<b<<")="<<value<<"\n";
        histo->SetBinContent(b, value);
        histo->SetBinError(b, 0.);
    }
    return histo;
}

/*****************************************************************/
TGraphAsymmErrors* Interpolator1D::graphFromHisto(TH1D* histo, TH1D* histoX)
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
TGraphAsymmErrors* Interpolator1D::addPointToGraph(TGraphAsymmErrors* graph, double x, double y, double xdo, double xup, double ydo, double yup)
/*****************************************************************/
{
    int n = graph->GetN();
    double* ax   = new double[n+1];
    double* ay   = new double[n+1];
    double* axdo = new double[n+1];
    double* aydo = new double[n+1];
    double* axup = new double[n+1];
    double* ayup = new double[n+1];
    int pp = 0;
    bool added = false;
    for(int p=0;p<n-1;p++)
    {
        double xp   = graph->GetX()[p];
        double yp   = graph->GetY()[p];
        double xdop = graph->GetEXlow()[p];
        double ydop = graph->GetEYlow()[p];
        double xupp = graph->GetEXhigh()[p];
        double yupp = graph->GetEYhigh()[p];

        double xpp = graph->GetX()[p+1];
        ax[pp]   = xp;
        ay[pp]   = yp;
        axdo[pp] = xdop;
        aydo[pp] = ydop;
        axup[pp] = xupp;
        ayup[pp] = yupp;
        pp++;
        if(xp<=x && xpp>x)
        {
            ax[pp]   = x;
            ay[pp]   = y;
            axdo[pp] = xdo;
            aydo[pp] = ydo;
            axup[pp] = xup;
            ayup[pp] = yup;
            pp++;
            added = true;
        }
    }
    double xp   = graph->GetX()[n-1];
    double yp   = graph->GetY()[n-1];
    double xdop = graph->GetEXlow()[n-1];
    double ydop = graph->GetEYlow()[n-1];
    double xupp = graph->GetEXhigh()[n-1];
    double yupp = graph->GetEYhigh()[n-1];

    ax[pp]   = xp;
    ay[pp]   = yp;
    axdo[pp] = xdop;
    aydo[pp] = ydop;
    axup[pp] = xupp;
    ayup[pp] = yupp;
    pp++;
    if(!added)
    {
        ax[pp]   = x;
        ay[pp]   = y;
        axdo[pp] = xdo;
        aydo[pp] = ydo;
        axup[pp] = xup;
        ayup[pp] = yup;
        pp++;
        added = true;
    }
    // delete the old graph and create the new one with the additional point
    const char* graphName = graph->GetName();
    graph->Delete();
    graph = new TGraphAsymmErrors(n+1, ax, ay, axdo,axup, aydo,ayup);
    graph->SetName(graphName);
    return graph;
}



/*****************************************************************/
void Interpolator1D::rebinHisto()
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


    double maxRelError = 1.;
    //while(maxRelError>0.5 && nBins%2==0)
    //{
        double meanValue = 0.;
        maxRelError = 0.;
        for(int b=1; b<=nBins; b++)
        {
            double value = m_rawHisto->GetBinContent(b);
            meanValue += value;
        }
        meanValue /= (double)nBins;
        for(int b=1; b<=nBins; b++)
        {
            double value = m_rawHisto->GetBinContent(b);
            double error = m_rawHisto->GetBinError(b);
            double relError = (value>0 ? error/value : 0);
            if(value>0.1*meanValue && relError>maxRelError) maxRelError = relError;
        }
        //cerr<<"Maximum relative error = "<<maxRelError<<"\n";
        if(maxRelError>0.5 && nBins%2==0)
        {
            histoRebin->Rebin(2);
            histoRebin->Scale(1./2.);
            nBins = histoRebin->GetNbinsX();
        }
    //}

    for(int b=1; b<=nBins; b++)
    {
        double x = histoRebin->GetXaxis()->GetBinCenter(b);
        histoXpos->SetBinContent(b, x);
    }


    int iteration = 1;
    // Iterative process. Two consecutive bins are merged together at each step.
    // Bins with the less significative differences are merged first. This is redone until all the differences are significant enough.
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
            double errorDown = histoRebin->GetBinError(b);
            double errorUp   = histoRebin->GetBinError(b+1);
            double relWidth  = (histoRebin->GetXaxis()->GetBinWidth(b)+histoRebin->GetXaxis()->GetBinWidth(b+1))/binWidth;

            if(valueUp>0 && valueDown>0 && (errorDown>0 || errorUp>0))
            {
                double diff      = fabs(valueUp-valueDown);
                double errorDiff = sqrt(errorUp*errorUp+errorDown*errorDown);

                double sig = diff/errorDiff;
                //if(sig<m_nSigmasRebin && sig*sqrt(relWidth-1.)<significanceMin) // favor merging of smaller bins
                if(sig<m_nSigmasRebin && sig*(relWidth-1.)<significanceMin)
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
void Interpolator1D::mergeLonelyBins()
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

                // merge bins only if the 2 neighbors have more than 3 times the width
                double sig = doubleDiff/errorDoubleDiff;
                if(sig<m_nSigmasLonely && 3.*binWidth<=binWidthUp && 3.*binWidth<=binWidthDown && sig<significanceMin)
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
void Interpolator1D::constrainInterpolation()
/*****************************************************************/
{
    //cerr<<"Constraining interpolation for "<<m_rawHisto->GetName()<<"\n";
    TH1D* interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
    double chi2Histo = chi2(m_rawHisto, interpHisto);
    int ndf = m_rawHisto->GetNbinsX() - m_interpGraph->GetN();
    double chi2oNdf = (ndf>0 ? chi2Histo/(double)ndf : chi2Histo);
    interpHisto->Delete();
    vector<int> constrainedBins;
    // if the chi2 is bad, try to constrain the interpolation with additional points
    while(chi2oNdf>2. && ndf>0)
    {
        //cout<<" Starting chi2oNdf="<<chi2oNdf<<"\n";
        int nbins = m_rawHisto->GetNbinsX();
        interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
        vector<pair<double,int> > diffAndPos;
        for(int b=1;b<=nbins;b++)
        {
            double value  = m_rawHisto->GetBinContent(b);
            double error  = m_rawHisto->GetBinError(b);
            double interp = interpHisto->GetBinContent(b);
            if(value>0 && (value-interp)>0 && find(constrainedBins.begin(), constrainedBins.end(), b)==constrainedBins.end())
            {
                diffAndPos.push_back( make_pair(fabs(value-interp)/error, b) );
                //cerr<<"  Can constrain bin "<<b<<" diff = "<<fabs(value-interp)/error<<"\n";
            }
            
            //if(value>0 && fabs(value-interp)/error > maxDiff && find(constrainedBins.begin(), constrainedBins.end(), b)==constrainedBins.end())
            //{
            //    maxDiff = fabs(value-interp)/error;
            //    posMax = b;
            //}
        }
        interpHisto->Delete();
        if(diffAndPos.size()==0) break;
        std::sort(diffAndPos.begin(), diffAndPos.end(), pair_comp()); 

        double minChi2oNdf = 1.e9;
        double bestShift   = 0.;
        int    bestPos     = 0;

        int ncands = diffAndPos.size();
        int ntries = min(5,ncands);
        for(int b=0;b<ntries;b++)
        {
            int pos      = diffAndPos[ncands-b-1].second;
            double value = m_rawHisto->GetBinContent(pos);
            double error = m_rawHisto->GetBinError(pos);
            double x = m_rawHisto->GetXaxis()->GetBinCenter(pos);
            for(int s=-10;s<=10;s++)
            {
                double shiftedValue = value+(double)s/10.*error; // range between +/- 1 sigma
                //cerr<<"  Trying to contrain at x="<<x<<",y="<<shiftedValue<<" (diff="<<diffAndPos[ncands-b-1].first<<"sigma)";
                stringstream nametmp;
                nametmp << m_interpGraph->GetName() << "_tmp";
                TGraphAsymmErrors* graphTmp = (TGraphAsymmErrors*)m_interpGraph->Clone(nametmp.str().c_str());
                graphTmp = addPointToGraph(graphTmp, x, shiftedValue, 0.,0., error, error);
                interpHisto = histoFromGraph(graphTmp, m_rawHisto);
                double c2 = chi2(m_rawHisto, interpHisto);
                ndf = m_rawHisto->GetNbinsX() - graphTmp->GetN();
                if(!std::isfinite(c2)) // Protection against crazy interpolations
                {
                    c2 = 1.e8; 
                }
                double c2oNdf = (ndf>0 ? c2/(double)ndf : c2);
                //cerr<<", c2oNdf="<<c2<<"/"<<ndf<<"="<<c2oNdf<<"\n";
                if(c2oNdf<minChi2oNdf)
                {
                    bestShift = shiftedValue;
                    bestPos = pos;
                    minChi2oNdf = c2oNdf;
                }
                interpHisto->Delete();
                graphTmp->Delete();
            }
        }
        //cerr<<"  Pushing bin "<<bestPos<<" to constrained bins\n";
        if(bestPos==0 && ncands>0) bestPos = diffAndPos[ncands-1].second; // protection if for some reason no best constraint is found
        constrainedBins.push_back(bestPos);
        if(minChi2oNdf<chi2oNdf)
        {
            m_interpGraph = addPointToGraph(m_interpGraph, m_rawHisto->GetXaxis()->GetBinCenter(bestPos), bestShift, 0.,0., m_rawHisto->GetBinError(bestPos), m_rawHisto->GetBinError(bestPos));
            interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
            chi2Histo = chi2(m_rawHisto, interpHisto);
            ndf = m_rawHisto->GetNbinsX() - m_interpGraph->GetN();
            chi2oNdf = (ndf>0 ? chi2Histo/(double)ndf : chi2Histo);
            //cerr<<"  Best constraint="<<m_rawHisto->GetXaxis()->GetBinCenter(bestPos)<<",y="<<bestShift<<"\n";
            //cerr<<" "<<m_interpGraph->GetN()<<" points "<<" ---> nomchi2="<<chi2oNdf<<"\n";
            interpHisto->Delete();
        }

    }
    constrainedBins.clear();
    ndf = m_rawHisto->GetNbinsX() - m_interpGraph->GetN();
    // Try to constrain two consecutive points if the chi2 is still too bad
    while(chi2oNdf>5. && ndf>1)
    {
        int nbins = m_rawHisto->GetNbinsX();
        interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
        vector<pair<double,int> > diffAndPos;
        for(int b=2;b<=nbins-1;b++)
        {
            double value = m_rawHisto->GetBinContent(b);
            double error = m_rawHisto->GetBinError(b);
            double interp = interpHisto->GetBinContent(b);

            if(value>0 && find(constrainedBins.begin(), constrainedBins.end(), b)==constrainedBins.end())
            {
                diffAndPos.push_back( make_pair(fabs(value-interp)/error, b) );
            }
            
        }
        interpHisto->Delete();
        if(diffAndPos.size()<=1) break;
        std::sort(diffAndPos.begin(), diffAndPos.end(), pair_comp()); 
        int ncand     = diffAndPos.size();
        int posMax    = diffAndPos[ncand-1].second;
        int scdPosMax = diffAndPos[ncand-2].second;
        if(abs(posMax-scdPosMax)>1) break;

        constrainedBins.push_back(posMax);
        constrainedBins.push_back(scdPosMax);
        double minChi2oNdf = 1.e9;
        double bestShift1 = 0.;
        double bestShift2 = 0.;
        double value1 = m_rawHisto->GetBinContent(posMax);
        double value2 = m_rawHisto->GetBinContent(scdPosMax);
        double error1 = m_rawHisto->GetBinError(posMax);
        double error2 = m_rawHisto->GetBinError(scdPosMax);

        double x1 = m_rawHisto->GetXaxis()->GetBinCenter(posMax);
        double x2 = m_rawHisto->GetXaxis()->GetBinCenter(scdPosMax);
        for(int s1=-10;s1<=10;s1++)
        {
            for(int s2=-10;s2<=10;s2++)
            {
                double shiftedValue1 = value1+(double)s1/10.*error1;
                double shiftedValue2 = value2+(double)s2/10.*error2;
                //cerr<<"  Trying to contrain at x1="<<x1<<",y1="<<shiftedValue1<<" (diff="<<diffAndPos[ncand-1].first<<"sigma)"<<"\n";
                //cerr<<"                    and x2="<<x2<<",y2="<<shiftedValue2<<" (diff="<<diffAndPos[ncand-2].first<<"sigma)";
                stringstream nametmp;
                nametmp << m_interpGraph->GetName() << "_tmp";
                TGraphAsymmErrors* graphTmp = (TGraphAsymmErrors*)m_interpGraph->Clone(nametmp.str().c_str());
                graphTmp = addPointToGraph(graphTmp, x1, shiftedValue1, 0.,0., error1, error1);
                graphTmp = addPointToGraph(graphTmp, x2, shiftedValue2, 0.,0., error2, error2);
                interpHisto = histoFromGraph(graphTmp, m_rawHisto);
                double c2 = chi2(m_rawHisto, interpHisto);
                ndf = m_rawHisto->GetNbinsX() - graphTmp->GetN();
                if(!std::isfinite(c2)) // Protection against crazy interpolations
                {
                    c2 = 1.e8; 
                }
                double c2oNdf = (ndf>0 ? c2/(double)ndf : c2);
                //cerr<<", c2oNdf="<<c2oNdf<<"\n";
                if(c2oNdf<minChi2oNdf)
                {
                    bestShift1 = shiftedValue1;
                    bestShift2 = shiftedValue2;
                    minChi2oNdf = c2oNdf;
                }
                interpHisto->Delete();
                graphTmp->Delete();
            }
        }
        if(minChi2oNdf<chi2oNdf)
        {

            m_interpGraph = addPointToGraph(m_interpGraph, x1, bestShift1, 0.,0., error1, error1);
            m_interpGraph = addPointToGraph(m_interpGraph, x2, bestShift2, 0.,0., error2, error2);
            interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
            chi2Histo = chi2(m_rawHisto, interpHisto);
            ndf = m_rawHisto->GetNbinsX() - m_interpGraph->GetN();
            chi2oNdf = (ndf>0 ? chi2Histo/(double)ndf : chi2Histo);
            //cerr<<"  Best constraint is x1="<<x1<<",y1="<<bestShift1<<"\n";
            //cerr<<"                 and x2="<<x2<<",y2="<<bestShift2<<"\n";
            //cout<<" "<<m_interpGraph->GetN()<<" points "<<" ---> nomchi2="<<chi2oNdf<<"\n";
            interpHisto->Delete();
        }

    }
}




/*****************************************************************/
void Interpolator1D::correctLocalBias()
/*****************************************************************/
{
    //cerr<<"Correcting bias for "<<m_rawHisto->GetName()<<"\n";
    TH1D* interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
    int nbins = m_rawHisto->GetNbinsX();
    bool done = false;
    vector<int> nPushs;
    for(int b=1;b<=nbins;b++)
    {
        nPushs.push_back(0);
    }
    while(!done)
    {
        bool foundBias = false;
        for(int b1=1;b1<=nbins;b1++)
        {
            if(nPushs[b1]>20) continue;
            int nBiased = 0;
            int signBias = 0;
            double mean = 0.;
            for(int b2=b1;b2<=nbins;b2++)
            {
                double raw = m_rawHisto->GetBinContent(b2);
                double interp = interpHisto->GetBinContent(b2);
                if(signBias*(interp-raw)<0) break; 
                if(signBias==0)
                {
                    signBias = (fabs(interp-raw)>0 ? (interp-raw)/fabs(interp-raw) : 0);
                }
                mean += raw;
                nBiased ++;
            }
            //cerr<<"bin "<<b1<<": "<<nBiased<<" consecutive bins with "<<signBias<<" bias\n";
            mean /= (double)nBiased;
            double xmin = m_rawHisto->GetXaxis()->GetBinLowEdge(b1);
            double xmax = m_rawHisto->GetXaxis()->GetBinUpEdge(b1+nBiased);
            if(nBiased>=4)// more than 3 consecutive bins with bias in the same direction
            {
                nPushs[b1]++;
                foundBias = true;
                int npoints = m_interpGraph->GetN();
                bool pushedPoint = false;
                for(int p=0;p<npoints;p++)
                {
                    double value = m_interpGraph->GetY()[p];
                    double x = m_interpGraph->GetX()[p];
                    if(x>=xmin && x<=xmax)
                    {
                        pushedPoint = true;
                        //double newValue = value - signBias*0.1*fabs(value-mean);
                        double newValue = value - signBias*0.01*value;
                        m_interpGraph->SetPoint(p, x, newValue);
                        //cerr<<"x="<<x<<", old="<<value<<", new="<<newValue<<"\n";
                    }
                }
                if(!pushedPoint)
                {
                    double xMiddle = (xmax+xmin)/2.;
                    m_interpGraph = addPointToGraph(m_interpGraph, xMiddle, mean, 0., 0., 0., 0.);
                }
                interpHisto->Delete();
                interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
                break;
            }
        }
        if(!foundBias) done = true;
    }
}



/*****************************************************************/
TH1D* Interpolator1D::interpolate(TH1D* rawHisto)
/*****************************************************************/
{

    if(m_rawHisto) m_rawHisto->Delete();
    if(m_rebinHisto) m_rebinHisto->Delete();
    if(m_rebinHistoX) m_rebinHistoX->Delete();
    if(m_interpHisto) m_interpHisto->Delete();
    if(m_interpGraph) m_interpGraph->Delete();
    m_rawHisto = NULL;
    m_rebinHisto = NULL;
    m_rebinHistoX = NULL;
    m_interpHisto = NULL;
    m_interpGraph = NULL;


    m_rawHisto = rawHisto;
    rebinHisto();
    mergeLonelyBins();
    if(m_rawHisto->GetNbinsX()>=5*m_rebinHisto->GetNbinsX()) m_linearInterp = true;
    //m_linearInterp = true;
    //cerr<<"Linear interpolation = "<<m_linearInterp<<"\n";
    m_interpGraph = graphFromHisto(m_rebinHisto, m_rebinHistoX);
    constrainInterpolation();
    m_interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);

    double chi2Histo = chi2(m_rawHisto, m_interpHisto);
    int ndf = m_rawHisto->GetNbinsX() - m_interpGraph->GetN();
    double chi2oNdf = (ndf>0 ? chi2Histo/(double)ndf : chi2Histo);
    cout<<"[INFO]     chi2/ndf(ref-raw) = "<<chi2Histo<<"/"<<ndf<<" = "<<chi2oNdf<<"\n";
    if(chi2oNdf>5)
    {
        cout<<"[WARN]     The reference histogram used for reweighting seems not good. Please check the produced control plot in the output file.\n";
    }

    //cerr<<"Diff integral = "<<(m_interpHisto->GetSumOfWeights() - m_rawHisto->GetSumOfWeights())/m_rawHisto->GetSumOfWeights()<<"\n";
    // if difference of normalisation with raw distribution > 5%, correct local bias
    //if(fabs( (m_interpHisto->GetSumOfWeights() - m_rawHisto->GetSumOfWeights())/m_rawHisto->GetSumOfWeights() ) > 0.05)
    //{
    //    correctLocalBias();
    //}
    //m_interpHisto->Delete();
    //m_interpHisto = NULL;
    //m_interpHisto = histoFromGraph(m_interpGraph, m_rawHisto);
    //cerr<<"Diff integral = "<<(m_interpHisto->GetSumOfWeights() - m_rawHisto->GetSumOfWeights())/m_rawHisto->GetSumOfWeights()<<"\n";

    //stringstream fileName;
    //fileName << m_rawHisto->GetName() << ".root";
    //m_outputFile = TFile::Open(fileName.str().c_str(), "RECREATE");
    //stringstream canvasName;
    //canvasName << m_rawHisto->GetName()<< "_c";
    //TCanvas* c = new TCanvas(canvasName.str().c_str(),canvasName.str().c_str(), 700,700);
    //m_rawHisto->SetMarkerStyle(20);
    //m_rebinHisto->SetLineColor(kBlue);
    //m_rebinHisto->SetLineWidth(2);
    //m_rebinHisto->SetMarkerColor(kBlue);
    //m_interpGraph->SetLineColor(kRed);
    //m_interpHisto->SetLineColor(kRed);
    //m_interpHisto->SetMarkerColor(kRed);
    //m_interpHisto->SetMarkerStyle(25);
    //m_rawHisto->Draw();
    //m_rebinHisto->Draw("same");
    //m_interpGraph->Draw("same");
    //m_interpHisto->Draw("same");
    //cerr<<m_interpHisto->GetNbinsX()<<"\n";
    //c->Write();
    //m_outputFile->Close();
    //c->Delete();

    return m_interpHisto;
}
