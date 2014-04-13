/**
 *  @file  Template.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/09/2013
 *
 *  @internal
 *     Created :  12/09/2013
 * Last update :  12/09/2013 09:12:22 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "Template.h"

#include "TH2F.h"
#include "TH3F.h"

#include <iostream>


using namespace std;

/*****************************************************************/
PostProcessing::PostProcessing(Type type)
/*****************************************************************/
{
    m_type = type;
}


/*****************************************************************/
PostProcessing::~PostProcessing()
/*****************************************************************/
{
}



///////////////////////////////////////////////////////////////////////////////////////////:
/*****************************************************************/
Template::Template():m_template(NULL),
    m_originalSumOfWeights(0.),
    m_conserveSumOfWeights(false)
/*****************************************************************/
{
}

/*****************************************************************/
Template::Template(const Template& tmp)
/*****************************************************************/
{
    stringstream hname;
    hname << tmp.getName() << "_";
    m_template = NULL;
    if(tmp.getTemplate())
    {
        m_template = dynamic_cast<TH2F*>(tmp.getTemplate()->Clone(hname.str().c_str()));
    }
    setWidths(tmp.getWidths());
    setRaw1DTemplates(tmp.getRaw1DTemplates());
    setOriginalSumOfWeights(tmp.originalSumOfWeights());
    m_conserveSumOfWeights = false;

}


/*****************************************************************/
Template::~Template()
/*****************************************************************/
{
    if(m_template)
    {
        delete m_template;
        m_template = NULL;
    }
    for(unsigned int axis=0;axis<m_widths.size();axis++)
    {
        if(m_widths[axis])
        {
            delete m_widths[axis];
        }
    }
    m_widths.clear();
    for(unsigned int axis=0;axis<m_raw1DTemplates.size();axis++)
    {
        if(m_raw1DTemplates[axis])
        {
            delete m_raw1DTemplates[axis];
        }
    }
    m_raw1DTemplates.clear();
    for(unsigned int n=0;n<m_controlPlots.size();n++)
    {
        if(m_controlPlots[n])
        {
            delete m_controlPlots[n];
        }
    }
    m_controlPlots.clear();
}

/*****************************************************************/
const string& Template::getVariable(unsigned int index) const
/*****************************************************************/
{
    assert(index<3);
    return m_variables.at(index);
}


/*****************************************************************/
TH1D* Template::getProjected1DTemplate(unsigned int axis)
/*****************************************************************/
{
    if(axis>=numberOfDimensions())
    {
        stringstream error;
        error << "Template::getProjected1DTemplate(): Projection requested on axis "<<axis<<" for "<<numberOfDimensions()<<"D template '"<<m_name<<"'\n";
        throw runtime_error(error.str());
    }
    unsigned int nbins1 = 0;
    unsigned int nbins2 = 0;
    if(axis==0)      
    {
        nbins1 = m_template->GetNbinsY();
        nbins2 = m_template->GetNbinsZ();
    }
    else if(axis==1) 
    {
        nbins1 = m_template->GetNbinsX();
        nbins2 = m_template->GetNbinsZ();
    }
    else if(axis==2) 
    {
        nbins1 = m_template->GetNbinsX();
        nbins2 = m_template->GetNbinsY();
    }

    stringstream projName;
    projName << m_template->GetName() << "_projFromTmp"<< axis;
    TH1D* projectedTemplate = NULL;

    if(numberOfDimensions()==2)
    {
        TH2F* tmp = dynamic_cast<TH2F*>(m_template);
        if(axis==0)      projectedTemplate = tmp->ProjectionX(projName.str().c_str(), 1, nbins1, "e");
        else if(axis==1) projectedTemplate = tmp->ProjectionY(projName.str().c_str(), 1, nbins1, "e");
    }
    if(numberOfDimensions()==3)
    {
        TH3F* tmp = dynamic_cast<TH3F*>(m_template);
        if(axis==0)      projectedTemplate = tmp->ProjectionX(projName.str().c_str(), 1, nbins1, 1, nbins2, "e");
        else if(axis==1) projectedTemplate = tmp->ProjectionY(projName.str().c_str(), 1, nbins1, 1, nbins2, "e");
        else if(axis==2) projectedTemplate = tmp->ProjectionZ(projName.str().c_str(), 1, nbins1, 1, nbins2, "e");
    }

    return projectedTemplate;
}

/*****************************************************************/
TH2D* Template::getProjected2DTemplate(unsigned int axis1, unsigned int axis2)
/*****************************************************************/
{
    if(axis1>=numberOfDimensions() || axis2>=numberOfDimensions())
    {
        stringstream error;
        error << "Template::getProjected2DTemplate(): Projection requested on axes "<<axis1<<" and "<<axis2<<" for "<<numberOfDimensions()<<"D template '"<<m_name<<"'\n";
        throw runtime_error(error.str());
    }
    if(numberOfDimensions()!=3)
    {
        stringstream error;
        error << "Template::getProjected2DTemplate(): Cannot do a 2D projection for "<<numberOfDimensions()<<"D template '"<<m_name<<"'\n";
        throw runtime_error(error.str());
    }
    unsigned int axis = 0;
    if(axis1==0 && axis2==1)      
    {
        axis = 2;
    }
    else if(axis1==0 && axis2==2) 
    {
        axis = 1;
    }
    else if(axis1==1 && axis2==2) 
    {
        axis = 0;
    }

    stringstream projName;
    projName << m_template->GetName() << "_proj2dFromTmp"<< axis;
    TH2D* projectedTemplate = NULL;

    if(numberOfDimensions()==3)
    {
        TH3F* tmp = dynamic_cast<TH3F*>(m_template);
        if(axis1==0 && axis2==1)      projectedTemplate = dynamic_cast<TH2D*>(tmp->Project3D("xye"));
        else if(axis1==0 && axis2==2) projectedTemplate = dynamic_cast<TH2D*>(tmp->Project3D("xze"));
        else if(axis1==1 && axis2==2) projectedTemplate = dynamic_cast<TH2D*>(tmp->Project3D("yze"));
        projectedTemplate->SetName(projName.str().c_str());
    }

    return projectedTemplate;
}


/*****************************************************************/
void Template::setVariables(const vector<string>& vars)
/*****************************************************************/
{
    m_variables.clear();
    for(unsigned int v=0;v<vars.size();v++)
    {
        m_variables.push_back(vars[v]);
    }
}

/*****************************************************************/
void Template::setTreeName(const std::string& name) 
/*****************************************************************/
{
    m_treeName = name;
    vector< pair<string, string> >::iterator it = m_inputFileAndTreeNames.begin();
    vector< pair<string, string> >::iterator itE = m_inputFileAndTreeNames.end();
    for(;it!=itE;++it)
    {
        it->second = name;
    }
}



/*****************************************************************/
void Template::createTemplate(const vector<unsigned int>& nbins, const vector< pair<double,double> >& minmax)
/*****************************************************************/
{
    // Clear histograms
    if(m_template)
        m_template->Delete();
    if(m_rawTemplate)
        delete m_rawTemplate;
    for(unsigned int axis=0;axis<m_widths.size();axis++)
    {
        if(m_widths[axis])
        {
            delete m_widths[axis];
        }
        m_widths.clear();
    }
    for(unsigned int axis=0;axis<m_raw1DTemplates.size();axis++)
    {
        if(m_raw1DTemplates[axis])
        {
            delete m_raw1DTemplates[axis];
        }
        m_raw1DTemplates.clear();
    }
    for(unsigned int i=0;i<m_raw2DTemplates.size();i++)
    {
        if(m_raw2DTemplates[i])
        {
            delete m_raw2DTemplates[i];
        }
        m_raw2DTemplates.clear();
    }

    m_minmax = minmax;
    // Initialize full ND templates
    if(nbins.size()==2)
    {
        m_template = new TH2F(m_name.c_str(),m_name.c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second);
        stringstream nameRaw;
        nameRaw << m_name << "_raw";
        m_rawTemplate = new TH2F(nameRaw.str().c_str(),nameRaw.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second);
    }
    else if(nbins.size()==3)
    {
        m_template = new TH3F(m_name.c_str(),m_name.c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second,nbins[2],minmax[2].first,minmax[2].second);
        stringstream nameRaw;
        nameRaw << m_name << "_raw";
        m_rawTemplate = new TH3F(nameRaw.str().c_str(),nameRaw.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second,nbins[2],minmax[2].first,minmax[2].second);
    }
    m_template->Sumw2();
    m_rawTemplate->Sumw2();

    // Initialize 1D projections and width maps
    for(unsigned int axis=0;axis<nbins.size();axis++)
    {
        stringstream nameWidth, nameProj;
        nameWidth << m_name << "_width" << axis;
        nameProj << m_name << "_proj" << axis;
        if(nbins.size()==2)
        {
            m_widths.push_back(new TH2F(nameWidth.str().c_str(),nameWidth.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second));
        }
        else if(nbins.size()==3)
        {
            m_widths.push_back(new TH3F(nameWidth.str().c_str(),nameWidth.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second,nbins[2],minmax[2].first,minmax[2].second));
        }
        m_widths.back()->Sumw2();
        m_raw1DTemplates.push_back(new TH1D(nameProj.str().c_str(),nameProj.str().c_str(),nbins[axis],minmax[axis].first,minmax[axis].second));
        m_raw1DTemplates.back()->Sumw2();
    }
    // Initialize 2D projections if ND=3
    if(nbins.size()==3)
    {
        stringstream nameProj1, nameProj2, nameProj3;
        nameProj1 << m_name << "_proj2d01";
        nameProj2 << m_name << "_proj2d02";
        nameProj2 << m_name << "_proj2d12";
        m_raw2DTemplates.push_back(new TH2D(nameProj1.str().c_str(),nameProj1.str().c_str(),nbins[0],minmax[0].first,minmax[0].second, nbins[1],minmax[1].first,minmax[1].second));
        m_raw2DTemplates.back()->Sumw2();
        m_raw2DTemplates.push_back(new TH2D(nameProj2.str().c_str(),nameProj2.str().c_str(),nbins[0],minmax[0].first,minmax[0].second, nbins[2],minmax[2].first,minmax[2].second));
        m_raw2DTemplates.back()->Sumw2();
        m_raw2DTemplates.push_back(new TH2D(nameProj3.str().c_str(),nameProj3.str().c_str(),nbins[1],minmax[1].first,minmax[1].second, nbins[2],minmax[2].first,minmax[2].second));
        m_raw2DTemplates.back()->Sumw2();
    }
}

/*****************************************************************/
void Template::setTemplate(const TH1* histo)
/*****************************************************************/
{
    if(m_template)
        m_template->Delete();
    m_template = dynamic_cast<TH1*>(histo->Clone(m_name.c_str()));
    //m_template->Sumw2();
}

/*****************************************************************/
void Template::setRawTemplate(const TH1* histo)
/*****************************************************************/
{
    if(m_rawTemplate)
        m_rawTemplate->Delete();
    stringstream nameRaw;
    nameRaw << m_name << "_raw";
    m_rawTemplate = dynamic_cast<TH1*>(histo->Clone(nameRaw.str().c_str()));
}

/*****************************************************************/
void Template::setRaw1DTemplates(const vector<TH1D*>& histos)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_raw1DTemplates.size();axis++)
    {
        if(m_raw1DTemplates[axis])
        {
            delete m_raw1DTemplates[axis];
        }
    }
    m_raw1DTemplates.clear();
    for(unsigned int axis=0;axis<histos.size();axis++)
    {
        stringstream name;
        name << m_name << "_proj" << axis;
        m_raw1DTemplates.push_back(dynamic_cast<TH1D*>(histos[axis]->Clone(name.str().c_str())));
    }

}

/*****************************************************************/
void Template::setRaw2DTemplates(const vector<TH2D*>& histos)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_raw2DTemplates.size();axis++)
    {
        if(m_raw2DTemplates[axis])
        {
            delete m_raw2DTemplates[axis];
        }
    }
    m_raw2DTemplates.clear();
    for(unsigned int axis=0;axis<histos.size();axis++)
    {
        stringstream name;
        name << m_name << "_proj2d" << axis;
        m_raw2DTemplates.push_back(dynamic_cast<TH2D*>(histos[axis]->Clone(name.str().c_str())));
    }

}

/*****************************************************************/
void Template::setWidths(const vector<TH1*>& widths)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_widths.size();axis++)
    {
        if(m_widths[axis])
        {
            delete m_widths[axis];
        }
        m_widths.clear();
    }
    for(unsigned int axis=0;axis<widths.size();axis++)
    {
        stringstream nameWidth;
        nameWidth << m_name << "_width" << axis;
        m_widths.push_back(dynamic_cast<TH1*>(widths[axis]->Clone(nameWidth.str().c_str())));
    }
}

/*****************************************************************/
bool Template::inTemplate(const vector<double>& vs)
/*****************************************************************/
{
    for(unsigned int d=0;d<vs.size();d++)
    {
        if(vs[d]>m_minmax[d].second || vs[d]<m_minmax[d].first)
        {
            return false;
        }
    }
    return true;
}

/*****************************************************************/
void Template::store(const vector<double>& vs, double w)
/*****************************************************************/
{
    m_entries.push_back(vs);
    m_weights.push_back(w);
}



/*****************************************************************/
void Template::reweight1D(unsigned int axis, unsigned int bin, double weight)
/*****************************************************************/
{
    if(axis>=numberOfDimensions())
    {
        stringstream error;
        error << "Template::reweight1D(): Reweighting requested on axis "<<axis<<" for "<<numberOfDimensions()<<"D template '"<<m_name<<"'\n";
        throw runtime_error(error.str());
    }
    unsigned int nbins1 = 0;
    unsigned int nbins2 = 0;
    if(axis==0)
    {
        nbins1 = m_template->GetNbinsY();
        nbins2 = m_template->GetNbinsZ();
    }
    else if(axis==1)
    {
        nbins1 = m_template->GetNbinsZ();
        nbins2 = m_template->GetNbinsX();
    }
    else if(axis==2)
    {
        nbins1 = m_template->GetNbinsX();
        nbins2 = m_template->GetNbinsY();
    }
    for(unsigned int b1=1;b1<=nbins1;b1++)
    {
        for(unsigned int b2=1;b2<=nbins2;b2++)
        {
            if(axis==0) 
            {
                double content = m_template->GetBinContent(bin,b1,b2);
                m_template->SetBinContent(bin,b1,b2, content*weight);
            }
            else if(axis==1) 
            {
                double content = m_template->GetBinContent(b2,bin,b1);
                m_template->SetBinContent(b2,bin,b1, content*weight);
            }
            else if(axis==2) 
            {
                double content = m_template->GetBinContent(b1,b2,bin);
                m_template->SetBinContent(b1,b2,bin, content*weight);
            }
        }
    }

}


/*****************************************************************/
void Template::makeProjectionControlPlot(const string& tag)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<numberOfDimensions();axis++)
    {
        stringstream plotName, rawName, projName;
        plotName << "control_" << getName() << "_projAxis" << axis << "_" << tag;
        rawName << "control_" << getName() << "_projAxis" << axis << "_" << tag << "_raw";
        projName << "control_" << getName() << "_projAxis" << axis << "_" << tag << "_proj";
        TCanvas* c = new TCanvas(plotName.str().c_str(),plotName.str().c_str(), 700,700);
        TH1D* raw1D = dynamic_cast<TH1D*>(getRaw1DTemplate(axis)->Clone(rawName.str().c_str()));
        TH1D* proj1D = dynamic_cast<TH1D*>(getProjected1DTemplate(axis));
        proj1D->SetName(projName.str().c_str());
        proj1D->SetLineColor(kRed);
        proj1D->SetLineWidth(2);
        raw1D->SetLineColor(kBlack);
        raw1D->SetMarkerColor(kBlack);
        raw1D->SetMarkerStyle(20);
        raw1D->SetXTitle(getVariable(axis).c_str());
        //raw1D->Draw("hist");
        raw1D->Draw();
        proj1D->Draw("hist same");
        addControlPlot(c);
    }
}

/*****************************************************************/
void Template::makeResidualsControlPlot(const string& tag, unsigned int rebin)
/*****************************************************************/
{
    if(numberOfDimensions()>0 && m_template->GetNbinsX()%rebin!=0) return;
    if(numberOfDimensions()>1 && m_template->GetNbinsY()%rebin!=0) return;
    if(numberOfDimensions()>2 && m_template->GetNbinsZ()%rebin!=0) return;

    stringstream cpName, cpRawName, resMapName, resDistName, relErrDistName;
    cpName << m_name << "_cp";
    cpRawName << m_name << "_rawcp";
    resMapName << m_name << "_resmap_" << tag << "_rebin" << rebin;
    resDistName << m_name << "_resdist_" << tag << "_rebin" << rebin;
    relErrDistName << m_name << "_relerrdist_" << tag << "_rebin" << rebin;
    TH1* cpTmp = NULL;
    TH1* cpRawTmp = NULL;
    if(rebin==1)
    {
        cpTmp = dynamic_cast<TH1*>(m_template->Clone(cpName.str().c_str()));
        cpRawTmp = dynamic_cast<TH1*>(m_rawTemplate->Clone(cpRawName.str().c_str()));
    }
    else
    {
        if(numberOfDimensions()==2)
        {
            cpTmp = dynamic_cast<TH2F*>(m_template)->Rebin2D(rebin, rebin, cpName.str().c_str());
            cpRawTmp = dynamic_cast<TH2F*>(m_rawTemplate)->Rebin2D(rebin, rebin, cpName.str().c_str());
        }
        else if(numberOfDimensions()==3)
        {
            cpTmp = dynamic_cast<TH3F*>(m_template)->Rebin3D(rebin, rebin, rebin, cpName.str().c_str());
            cpRawTmp = dynamic_cast<TH3F*>(m_rawTemplate)->Rebin3D(rebin, rebin, rebin, cpName.str().c_str());
        }
    }
    TH1* resMap = dynamic_cast<TH1*>(cpTmp->Clone(resMapName.str().c_str()));
    TH1D* resDist = new TH1D(resDistName.str().c_str(), resDistName.str().c_str(), 30, -3, 3);
    resDist->StatOverflows();
    TH1D* relErrDist = new TH1D(relErrDistName.str().c_str(), relErrDistName.str().c_str(), 200, -1, 1);
    relErrDist->StatOverflows();
    unsigned int nbins1 = cpTmp->GetNbinsX();
    unsigned int nbins2 = cpTmp->GetNbinsY();
    unsigned int nbins3 = cpTmp->GetNbinsZ();
    for(unsigned int b1=1;b1<=nbins1;b1++)
    {
        for(unsigned int b2=1;b2<=nbins2;b2++)
        {
            for(unsigned int b3=1;b3<=nbins3;b3++)
            {
                double tmpValue = cpTmp->GetBinContent(b1,b2,b3);
                double tmpRawValue = cpRawTmp->GetBinContent(b1,b2,b3);
                double tmpRawError = cpRawTmp->GetBinError(b1,b2,b3);
                if(tmpValue>0. && tmpRawValue>0. && tmpRawError>0.)
                {
                    double res = (tmpRawValue-tmpValue)/tmpRawError;
                    double relErr = (tmpRawValue-tmpValue)/tmpRawValue;
                    if(numberOfDimensions()==2)
                    {
                        resMap->SetBinContent(b1,b2,b3, res);
                        resMap->SetBinError(b1,b2,b3, 0.);
                    }
                    resDist->Fill(res);
                    relErrDist->Fill(relErr);
                }
            }
        }
    }
    stringstream plotMapName, plotDistName, plotErrDistName;
    plotMapName << "control_" << getName() << "_resMap" << "_" << tag << "_rebin" << rebin;
    plotDistName << "control_" << getName() << "_resDist" << "_" << tag << "_rebin" << rebin;
    plotErrDistName << "control_" << getName() << "_relErrDist" << "_" << tag << "_rebin" << rebin;
    if(numberOfDimensions()==2)
    {
        TCanvas* c = new TCanvas(plotMapName.str().c_str(),plotMapName.str().c_str(), 700,700);
        resMap->SetContour(99);
        resMap->SetAxisRange(-3., 3., "z");
        resMap->Draw("color z");
        resMap->SetXTitle(getVariable(0).c_str());
        resMap->SetYTitle(getVariable(1).c_str());
        addControlPlot(c);
    }
    TCanvas* c2 = new TCanvas(plotDistName.str().c_str(),plotDistName.str().c_str(), 700,700);
    resDist->SetLineColor(kBlack);
    resDist->SetLineWidth(2);
    resDist->SetMarkerColor(kBlack);
    resDist->SetMarkerStyle(20);
    resDist->SetXTitle("(raw-template)/error_{raw}");
    resDist->Draw();
    addControlPlot(c2);

    TCanvas* c3 = new TCanvas(plotErrDistName.str().c_str(),plotErrDistName.str().c_str(), 700,700);
    relErrDist->SetLineColor(kBlack);
    relErrDist->SetLineWidth(2);
    relErrDist->SetMarkerColor(kBlack);
    relErrDist->SetMarkerStyle(20);
    relErrDist->SetXTitle("(raw-template)/raw");
    relErrDist->Draw();
    addControlPlot(c3);

    if(cpTmp) delete cpTmp;
    if(cpRawTmp) delete cpRawTmp;
}

