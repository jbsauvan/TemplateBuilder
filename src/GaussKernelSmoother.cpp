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
#include <TH2F.h>
#include <TH3F.h>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace std;


/*****************************************************************/
GaussKernelSmoother::GaussKernelSmoother():
    m_ndim(2),
    m_widthScalingFactor(1.)
/*****************************************************************/
{
}

/*****************************************************************/
GaussKernelSmoother::GaussKernelSmoother(unsigned int ndim):
    m_ndim(ndim),
    m_widthScalingFactor(1.)
/*****************************************************************/
{
}

/*****************************************************************/
GaussKernelSmoother::~GaussKernelSmoother()
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_widths.size();axis++)
    {
        if(m_widths[axis]) m_widths[axis]->Delete();
    }
    m_widths.clear();
}

/*****************************************************************/
TH1* GaussKernelSmoother::smooth(const TH1* histo)
/*****************************************************************/
{
    if(m_ndim!=2 && m_ndim!=3)
    {
        stringstream error;
        error << "GaussKernelSmoother::smooth(): Smoothing in 2D or 3D only is implemented";
        throw runtime_error(error.str());
    }
    TH1* smoothedHisto = NULL;
    if(m_ndim==2)
    {
        int nbinsx = histo->GetNbinsX();
        int nbinsy = histo->GetNbinsY();
        stringstream hName;
        hName << histo->GetName() << "_smooth";
        TH2F* smoothedHisto2F = dynamic_cast<TH2F*>(histo->Clone(hName.str().c_str()));
        smoothedHisto2F->SetDirectory(0);
        int counter = 0;
        int total = nbinsx*nbinsy;
        for(int bx=1;bx<=nbinsx;bx++)
        {
            for(int by=1;by<=nbinsy;by++)
            {
                if(counter % (total/100) == 0)
                {
                    double ratio  =  (double)counter/(double)total;
                    int   c      =  ratio * 50;
                    cout << "[INFO]   "<< setw(3) << (int)(ratio*100) << "% [";
                    for (int x=0; x<c; x++) cout << "=";
                    for (int x=c; x<50; x++) cout << " ";
                    cout << "]\r" << flush;
                }
                counter++;
                double x0 = histo->GetXaxis()->GetBinCenter(bx);
                double y0 = histo->GetYaxis()->GetBinCenter(by);
                vector<double> point;
                point.push_back(x0);
                point.push_back(y0);
                pair<double,double> valueError = smoothedValueError(histo, point);
                //cout<<"Bin "<<bx<<","<<by<<": smoothed = "<<valueError.first<<"\n";
                smoothedHisto2F->SetBinContent(bx, by, valueError.first);
                smoothedHisto2F->SetBinError(bx, by, valueError.second);
            }
        }
        smoothedHisto = smoothedHisto2F;
    }
    else if(m_ndim==3)
    {
        int nbinsx = histo->GetNbinsX();
        int nbinsy = histo->GetNbinsY();
        int nbinsz = histo->GetNbinsZ();
        stringstream hName;
        hName << histo->GetName() << "_smooth";
        TH3F* smoothedHisto3F = dynamic_cast<TH3F*>(histo->Clone(hName.str().c_str()));
        smoothedHisto3F->SetDirectory(0);
        int counter = 0;
        int total = nbinsx*nbinsy*nbinsz;
        for(int bx=1;bx<=nbinsx;bx++)
        {
            for(int by=1;by<=nbinsy;by++)
            {
                for(int bz=1;bz<=nbinsz;bz++)
                {
                    if(counter % (total/100) == 0)
                    {
                        double ratio  =  (double)counter/(double)total;
                        int   c      =  ratio * 50;
                        cout << "[INFO]   "<< setw(3) << (int)(ratio*100) << "% [";
                        for (int x=0; x<c; x++) cout << "=";
                        for (int x=c; x<50; x++) cout << " ";
                        cout << "]\r" << flush;
                    }
                    counter++;
                    double x0 = histo->GetXaxis()->GetBinCenter(bx);
                    double y0 = histo->GetYaxis()->GetBinCenter(by);
                    double z0 = histo->GetZaxis()->GetBinCenter(bz);
                    vector<double> point;
                    point.push_back(x0);
                    point.push_back(y0);
                    point.push_back(z0);
                    pair<double,double> valueError = smoothedValueError(histo, point);
                    //cout<<"Bin "<<bx<<","<<by<<","<<bz<<": smoothed = "<<valueError.first<<"\n";
                    smoothedHisto3F->SetBinContent(bx, by, bz, valueError.first);
                    smoothedHisto3F->SetBinError(bx, by, bz, valueError.second);
                }
            }
        }
        smoothedHisto = smoothedHisto3F;
    }
    cout << "[INFO]   "<< setw(3) << 100 << "% [";
    for (int x=0; x<50; x++) cout << "=";
    cout << "]" << endl;
    return smoothedHisto;
}



/*****************************************************************/
void GaussKernelSmoother::setWidths(const std::vector<TH1*>& widths)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_widths.size();axis++)
    {
        if(m_widths[axis]) m_widths[axis]->Delete();
    }
    m_widths.clear();
    for(unsigned int axis=0;axis<widths.size();axis++)
    {
        stringstream name;
        name << "smoother_width" << axis;
        m_widths.push_back( (TH1*)widths[axis]->Clone(name.str().c_str()) );

    }
}


/*****************************************************************/
double GaussKernelSmoother::weight(const std::vector<double>& x0, const std::vector<double>& xi)
/*****************************************************************/
{
    // N-dimension weight
    double dr2 = 0.;
    for(unsigned int axis=0;axis<x0.size();axis++)
    {
        double width = 1.;
        if(m_ndim==2)
        {
            int binx = m_widths[axis]->GetXaxis()->FindBin(x0[0]);
            int biny = m_widths[axis]->GetYaxis()->FindBin(x0[1]);
            width = m_widths[axis]->GetBinContent(binx, biny)/2.;
        }
        else if(m_ndim==3)
        {
            int binx = m_widths[axis]->GetXaxis()->FindBin(x0[0]);
            int biny = m_widths[axis]->GetYaxis()->FindBin(x0[1]);
            int binz = m_widths[axis]->GetZaxis()->FindBin(x0[2]);
            width = m_widths[axis]->GetBinContent(binx, biny, binz)/2.;
        }
        width *= m_widthScalingFactor;
        //cout<<"  Width="<<width<<"\n";
        double dx = (xi[axis]-x0[axis])/width;
        dr2 += dx*dx;
        //cout<<"  dr2 = "<<dr2<<"\n";
    }
    double wi = TMath::Gaus( sqrt(dr2) );

    return wi;
}

/*****************************************************************/
double GaussKernelSmoother::weight(const std::vector<double>& x0, double dx, int axis)
/*****************************************************************/
{
    // 1-dimension weight (can be factorized to obtain N-dimension weights)
    double width = 1.;
    if(m_ndim==2)
    {
        int binx = m_widths[axis]->GetXaxis()->FindBin(x0[0]);
        int biny = m_widths[axis]->GetYaxis()->FindBin(x0[1]);
        width = m_widths[axis]->GetBinContent(binx, biny)/2.;
    }
    else if(m_ndim==3)
    {
        int binx = m_widths[axis]->GetXaxis()->FindBin(x0[0]);
        int biny = m_widths[axis]->GetYaxis()->FindBin(x0[1]);
        int binz = m_widths[axis]->GetZaxis()->FindBin(x0[2]);
        width = m_widths[axis]->GetBinContent(binx, biny, binz)/2.;
    }
    width *= m_widthScalingFactor;
    double dxw = dx/width;
    double wi = TMath::Gaus( dxw );
    return wi;
}

/*****************************************************************/
pair<double,double> GaussKernelSmoother::smoothedValueError(const TH1* histo, const std::vector<double>& x0)
/*****************************************************************/
{
    if(m_ndim==2)
    {
        return smoothed2DValueError(histo, x0);
    }
    else if(m_ndim==3)
    {
        return smoothed3DValueError(histo, x0);
    }
    return make_pair<double,double>(0.,0.);
}

/*****************************************************************/
pair<double,double> GaussKernelSmoother::smoothed2DValueError(const TH1* histo, const std::vector<double>& x0)
/*****************************************************************/
{
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    double widthx = 1.;
    double widthy = 1.;
    int binx = m_widths[0]->GetXaxis()->FindBin(x0[0]);
    int biny = m_widths[0]->GetYaxis()->FindBin(x0[1]);
    widthx = m_widths[0]->GetBinContent(binx, biny);
    widthx *= m_widthScalingFactor;
    binx = m_widths[1]->GetXaxis()->FindBin(x0[0]);
    biny = m_widths[1]->GetYaxis()->FindBin(x0[1]);
    widthy = m_widths[1]->GetBinContent(binx, biny);
    widthy *= m_widthScalingFactor;
    double maxWidth = max(widthx, widthy);
    double widthRatiox = widthx/maxWidth;
    double widthRatioy = widthy/maxWidth;

    binx = histo->GetXaxis()->FindBin(x0[0]);
    biny = histo->GetYaxis()->FindBin(x0[1]);
    double binWidthx = histo->GetXaxis()->GetBinWidth(binx);
    double binWidthy = histo->GetYaxis()->GetBinWidth(biny);
    // FIXME: Gaussian is truncated at 2sigma. We may want to be able to configure this cut
    int nbinsWidthX = 2.*widthx/binWidthx;
    int nbinsWidthY = 2.*widthy/binWidthy;
    //cout<<"widthx="<<widthx<<" widthy="<<widthy<<" binwidthx="<<binWidthx<<" binwidthy="<<binWidthy<<"\n";
    //cout<<"nbinsWidthX="<<nbinsWidthX<<"\n";
    //cout<<"nbinsWidthY="<<nbinsWidthY<<"\n";

    // First compute the factorized weights in each direction
    vector<double> weightsX;
    vector<double> weightsY;
    // FIXME: this assumes that all bins have the same size
    for(int dbx=0;dbx<=nbinsWidthX;dbx++)
    {
        double dx = (double)dbx*binWidthx;
        double wxi = weight(x0, dx, 0);
        weightsX.push_back(wxi);
    }
    for(int dby=0;dby<=nbinsWidthY;dby++)
    {
        double dy = (double)dby*binWidthy;
        double wyi = weight(x0, dy, 1);
        weightsY.push_back(wyi);
    }


    double sumw = 0.;
    double sumwv = 0.;
    double sumwe = 0.;
    //cout<<"Smoothing "<<x0<<","<<y0<<"\n";
    //for(int bx=max(binx-nbinsWidthX,1);bx<=min(binx+nbinsWidthX,nbinsx);bx++)
    for(int bx=binx-nbinsWidthX;bx<=binx+nbinsWidthX;bx++)
    {
        int bxcp = bx;
        //double xshift = 0.;
        if(bx<1) 
        {
            bxcp = 1;
            //xshift = (bx-1)*binWidthx;
        }
        else if(bx>nbinsx) 
        {
            bxcp = nbinsx;
            //xshift = (bx-nbinsx)*binWidthx;
        }
        // FIXME: this assumes that all bins have the same size
        int dbx = abs(bx-binx);
        double wi1 = weightsX[dbx];
        //for(int by=max(biny-nbinsWidthY,1);by<=min(biny+nbinsWidthY,nbinsy);by++)
        for(int by=biny-nbinsWidthY;by<=biny+nbinsWidthY;by++)
        {
            int bycp = by;
            //double yshift = 0.;
            if(by<1) 
            {
                bycp = 1;
                //yshift = (by-1)*binWidthy;
            }
            else if(by>nbinsy) 
            {
                bycp = nbinsy;
                //yshift = (by-nbinsy)*binWidthy;
            }
            // FIXME: this assumes that all bins have the same size
            int dby = abs(by-biny);
            double dbr = sqrt( (double)(dbx*dbx)*widthRatiox+(double)(dby*dby)*widthRatioy);
            double wi = wi1*weightsY[dby];
            wi *= 1./(dbr+1.);
            //////////////////////////
            //double xi = histo->GetXaxis()->GetBinCenter(bxcp) + xshift;
            //double yi = histo->GetYaxis()->GetBinCenter(bycp) + yshift;
            //vector<double> xxi;
            //xxi.push_back(xi);
            //xxi.push_back(yi);
            //double wi = weight(x0, xxi);
            double vi = histo->GetBinContent(bxcp,bycp);
            double ei = histo->GetBinError(bxcp,bycp);
            sumw += wi;
            sumwv += wi*vi;
            sumwe += wi*ei;
            //cout<<"  Bin "<<bx<<", "<<by<<": weight="<<wi<<" z="<<zi<<"\n";
        }
    }
    double value = 0.;
    double error = 0.;
    if(sumw>0.)
    {
        value = sumwv/sumw;
        error = sumwe/sumw;
    }
    else
    {
        //cout<<"sumofweights=0\n";
    }

    return make_pair(value,error);
}

/*****************************************************************/
pair<double,double> GaussKernelSmoother::smoothed3DValueError(const TH1* histo, const std::vector<double>& x0)
/*****************************************************************/
{
    int nbinsx = histo->GetNbinsX();
    int nbinsy = histo->GetNbinsY();
    int nbinsz = histo->GetNbinsZ();
    double widthx = 1.;
    double widthy = 1.;
    double widthz = 1.;
    int binx = m_widths[0]->GetXaxis()->FindBin(x0[0]);
    int biny = m_widths[0]->GetYaxis()->FindBin(x0[1]);
    int binz = m_widths[0]->GetZaxis()->FindBin(x0[2]);
    widthx = m_widths[0]->GetBinContent(binx, biny, binz);
    widthx *= m_widthScalingFactor;
    binx = m_widths[1]->GetXaxis()->FindBin(x0[0]);
    biny = m_widths[1]->GetYaxis()->FindBin(x0[1]);
    binz = m_widths[1]->GetZaxis()->FindBin(x0[2]);
    widthy = m_widths[1]->GetBinContent(binx, biny, binz);
    widthy *= m_widthScalingFactor;
    binx = m_widths[2]->GetXaxis()->FindBin(x0[0]);
    biny = m_widths[2]->GetYaxis()->FindBin(x0[1]);
    binz = m_widths[2]->GetZaxis()->FindBin(x0[2]);
    widthz = m_widths[2]->GetBinContent(binx, biny, binz);
    widthz *= m_widthScalingFactor;
    double maxWidth = max(max(widthx, widthy),widthz);
    double widthRatiox = widthx/maxWidth;
    double widthRatioy = widthy/maxWidth;
    double widthRatioz = widthz/maxWidth;

    histo->GetXaxis()->FindBin(x0[0]);
    histo->GetYaxis()->FindBin(x0[1]);
    histo->GetZaxis()->FindBin(x0[2]);
    double binWidthx = histo->GetXaxis()->GetBinWidth(binx);
    double binWidthy = histo->GetYaxis()->GetBinWidth(biny);
    double binWidthz = histo->GetZaxis()->GetBinWidth(binz);
    // FIXME: Gaussian is truncated at 2sigma. We may want to be able to configure this cut
    int nbinsWidthX = 2.*widthx/binWidthx;
    int nbinsWidthY = 2.*widthy/binWidthy;
    int nbinsWidthZ = 2.*widthz/binWidthz;
    //cout<<"widthx="<<widthx<<" widthy="<<widthy<<" widthz="<<widthz<<" binwidthx="<<binWidthx<<" binwidthy="<<binWidthy<<" binwidthz="<<binWidthz<<"\n";
    //cout<<"nbinsWidthX="<<nbinsWidthX<<"\n";
    //cout<<"nbinsWidthY="<<nbinsWidthY<<"\n";
    //cout<<"nbinsWidthZ="<<nbinsWidthZ<<"\n";

    double sumw = 0.;
    double sumwv = 0.;
    double sumwe = 0.;
    //cout<<"Smoothing "<<x0<<","<<y0<<"\n";
    //for(int bx=max(binx-nbinsWidthX,1);bx<=min(binx+nbinsWidthX,nbinsx);bx++)

    // First compute the factorized weights in each direction
    vector<double> weightsX;
    vector<double> weightsY;
    vector<double> weightsZ;
    // FIXME: this assumes that all bins have the same size
    for(int dbx=0;dbx<=nbinsWidthX;dbx++)
    {
        double dx = (double)dbx*binWidthx;
        double wxi = weight(x0, dx, 0);
        weightsX.push_back(wxi);
    }
    for(int dby=0;dby<=nbinsWidthY;dby++)
    {
        double dy = (double)dby*binWidthy;
        double wyi = weight(x0, dy, 1);
        weightsY.push_back(wyi);
    }
    for(int dbz=0;dbz<=nbinsWidthZ;dbz++)
    {
        double dz = (double)dbz*binWidthz;
        double wzi = weight(x0, dz, 2);
        weightsZ.push_back(wzi);
    }


    for(int bx=binx-nbinsWidthX;bx<=binx+nbinsWidthX;bx++)
    {
        int bxcp = bx;
        //double xshift = 0.;
        if(bx<1) 
        {
            bxcp = 1;
            //xshift = (bx-1)*binWidthx;
        }
        else if(bx>nbinsx) 
        {
            bxcp = nbinsx;
            //xshift = (bx-nbinsx)*binWidthx;
        }
        // FIXME: this assumes that all bins have the same size
        int dbx = abs(bx-binx);
        double wi1 = weightsX[dbx];
        //for(int by=max(biny-nbinsWidthY,1);by<=min(biny+nbinsWidthY,nbinsy);by++)
        for(int by=biny-nbinsWidthY;by<=biny+nbinsWidthY;by++)
        {
            int bycp = by;
            //double yshift = 0.;
            if(by<1) 
            {
                bycp = 1;
                //yshift = (by-1)*binWidthy;
            }
            else if(by>nbinsy) 
            {
                bycp = nbinsy;
                //yshift = (by-nbinsy)*binWidthy;
            }
            // FIXME: this assumes that all bins have the same size
            int dby = abs(by-biny);
            double wi2 = wi1*weightsY[dby];
            for(int bz=binz-nbinsWidthZ;bz<=binz+nbinsWidthZ;bz++)
            {
                int bzcp = bz;
                //double zshift = 0.;
                if(bz<1) 
                {
                    bzcp = 1;
                    //zshift = (bz-1)*binWidthz;
                }
                else if(bz>nbinsz) 
                {
                    bzcp = nbinsz;
                    //zshift = (bz-nbinsz)*binWidthz;
                }
                // FIXME: this assumes that all bins have the same size
                int dbz = abs(bz-binz);
                //double dbr = sqrt( (double)(dbx*dbx+dby*dby+dbz*dbz) );
                double dbr = sqrt( (double)(dbx*dbx)*widthRatiox+(double)(dby*dby)*widthRatioy+(double)(dbz*dbz)*widthRatioz);
                double wi = wi2*weightsZ[dbz];
                wi *= 1./((dbr+1.)*(dbr+1.));
                //////////////////////////////////////////////////////////////
                //double xi = histo->GetXaxis()->GetBinCenter(bxcp) + xshift;
                //double yi = histo->GetYaxis()->GetBinCenter(bycp) + yshift;
                //double zi = histo->GetZaxis()->GetBinCenter(bzcp) + zshift;
                //vector<double> xxi;
                //xxi.push_back(xi);
                //xxi.push_back(yi);
                //xxi.push_back(zi);
                //double wi = weight(x0, xxi);
                //double wi = weightsX[dbx]*weightsY[dby]*weightsZ[dbz];
                double vi = histo->GetBinContent(bxcp,bycp,bzcp);
                double ei = histo->GetBinError(bxcp,bycp,bzcp);
                sumw += wi;
                sumwv += wi*vi;
                sumwe += wi*ei;
                //cout<<"  Bin "<<bx<<", "<<by<<","<<bz<<": weight="<<wi<<" vi="<<vi<<"\n";
            }
        }
    }
    double value = 0.;
    double error = 0.;
    if(sumw>0.)
    {
        value = sumwv/sumw;
        error = sumwe/sumw;
    }
    else
    {
        //cout<<"sumofweights=0\n";
    }

    return make_pair(value,error);
}
