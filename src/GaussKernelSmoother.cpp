/**
 *  @file  GaussKernelSmoother.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/13/2013
 *
 *  @internal
 *     Created :  12/13/2013
 * Last update :  12/13/2013 02:30:44 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "GaussKernelSmoother.h"

#include <TMath.h>

#include <sstream>
#include <iostream>

using namespace std;



/*****************************************************************/
GaussKernelSmoother::GaussKernelSmoother()
/*****************************************************************/
{
    m_widthx = NULL;
    m_widthy = NULL;
}

/*****************************************************************/
GaussKernelSmoother::~GaussKernelSmoother()
/*****************************************************************/
{
    if(m_widthx) m_widthx->Delete();
    if(m_widthy) m_widthy->Delete();
}

/*****************************************************************/
TH2F* GaussKernelSmoother::smooth(const TH2F* histo)
/*****************************************************************/
{
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    stringstream hName;
    hName << histo->GetName() << "_smooth";
    TH2F* smoothedHisto = dynamic_cast<TH2F*>(histo->Clone(hName.str().c_str()));
    smoothedHisto->SetDirectory(0);
    for(int bx=1;bx<=nbinsx;bx++)
    {
        for(int by=1;by<=nbinsy;by++)
        {
            double x0 = histo->GetXaxis()->GetBinCenter(bx);
            double y0 = histo->GetYaxis()->GetBinCenter(by);
            pair<double,double> valueError = smoothedValueError(histo, x0, y0);
            //cout<<"Bin "<<bx<<","<<by<<": smoothed = "<<valueError.first<<"\n";
            smoothedHisto->SetBinContent(bx, by, valueError.first);
            smoothedHisto->SetBinError(bx, by, valueError.second);
        }
    }
    return smoothedHisto;
}



/*****************************************************************/
void GaussKernelSmoother::setWidths(TH2F* widthx, TH2F* widthy)
/*****************************************************************/
{
    if(m_widthx)
    {
        m_widthx->Delete();
        m_widthx = NULL;
    }
    if(m_widthy)
    {
        m_widthy->Delete();
        m_widthy = NULL;
    }
    m_widthx = (TH2F*)widthx->Clone("smoother_widthx");
    m_widthy = (TH2F*)widthy->Clone("smoother_widthy");
}


/*****************************************************************/
double GaussKernelSmoother::weight(double x0, double y0, double xi, double yi)
/*****************************************************************/
{
    double widthx = 1.;
    double widthy = 1.;
    if(m_widthx)
    {
        int binx = m_widthx->GetXaxis()->FindBin(x0);
        int biny = m_widthx->GetYaxis()->FindBin(y0);
        widthx = m_widthx->GetBinContent(binx, biny)/2.;
    }
    if(m_widthy)
    {
        int binx = m_widthy->GetXaxis()->FindBin(x0);
        int biny = m_widthy->GetYaxis()->FindBin(y0);
        widthy = m_widthy->GetBinContent(binx, biny)/2.;
    }
    double dx = (xi-x0)/widthx;
    double dy = (yi-y0)/widthy;
    double dr = sqrt(dx*dx+dy*dy);
    double wi = TMath::Gaus(dr);

    return wi;
}

/*****************************************************************/
pair<double,double> GaussKernelSmoother::smoothedValueError(const TH2F* histo, double x0, double y0)
/*****************************************************************/
{
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    double widthx = 1.;
    double widthy = 1.;
    if(m_widthx)
    {
        int binx = m_widthx->GetXaxis()->FindBin(x0);
        int biny = m_widthx->GetYaxis()->FindBin(y0);
        widthx = m_widthx->GetBinContent(binx, biny);
    }
    if(m_widthy)
    {
        int binx = m_widthy->GetXaxis()->FindBin(x0);
        int biny = m_widthy->GetYaxis()->FindBin(y0);
        widthy = m_widthy->GetBinContent(binx, biny);
    }
    int binx = histo->GetXaxis()->FindBin(x0);
    int biny = histo->GetYaxis()->FindBin(y0);
    double binWidthx = histo->GetXaxis()->GetBinWidth(binx);
    double binWidthy = histo->GetYaxis()->GetBinWidth(biny);
    int nbinsWidthX = 3.*widthx/binWidthx;
    int nbinsWidthY = 3.*widthy/binWidthy;
    //cout<<"widthx="<<widthx<<" widthy="<<widthy<<" binwidthx="<<binWidthx<<" binwidthy="<<binWidthy<<"\n";
    //cout<<"nbinsWidthX="<<nbinsWidthX<<"\n";
    //cout<<"nbinsWidthY="<<nbinsWidthY<<"\n";

    double sumw = 0.;
    double sumwz = 0.;
    double sumwze = 0.;
    //cout<<"Smoothing "<<x0<<","<<y0<<"\n";
    //for(int bx=max(binx-nbinsWidthX,1);bx<=min(binx+nbinsWidthX,nbinsx);bx++)
    for(int bx=binx-nbinsWidthX;bx<=binx+nbinsWidthX;bx++)
    {
        int bxcp = bx;
        double xshift = 0.;
        if(bx<1) 
        {
            bxcp = 1;
            xshift = (bx-1)*binWidthx;
        }
        else if(bx>nbinsx) 
        {
            bxcp = nbinsx;
            xshift = (bx-nbinsx)*binWidthx;
        }
        //for(int by=max(biny-nbinsWidthY,1);by<=min(biny+nbinsWidthY,nbinsy);by++)
        for(int by=biny-nbinsWidthY;by<=biny+nbinsWidthY;by++)
        {
            int bycp = by;
            double yshift = 0.;
            if(by<1) 
            {
                bycp = 1;
                yshift = (by-1)*binWidthy;
            }
            else if(by>nbinsy) 
            {
                bycp = nbinsy;
                yshift = (by-nbinsy)*binWidthy;
            }
            double xi = histo->GetXaxis()->GetBinCenter(bxcp) + xshift;
            double yi = histo->GetYaxis()->GetBinCenter(bycp) + yshift;
            double wi = weight(x0, y0, xi, yi);
            double zi = histo->GetBinContent(bxcp,bycp);
            double zei = histo->GetBinError(bxcp,bycp);
            sumw += wi;
            sumwz += wi*zi;
            sumwze += wi*zei;
            //cout<<"  Bin "<<bx<<", "<<by<<": weight="<<wi<<" z="<<zi<<"\n";
        }
    }
    double value = 0.;
    double error = 0.;
    if(sumw>0.)
    {
        value = sumwz/sumw;
        error = sumwze/sumw;
    }
    else
    {
        cout<<"sumofweights=0\n";
    }

    return make_pair(value,error);
}
