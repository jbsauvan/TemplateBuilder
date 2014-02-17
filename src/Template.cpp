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
Template::Template():m_template(NULL)
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
    if(m_template)
        m_template->Delete();
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
    m_minmax = minmax;
    if(nbins.size()==2)
    {
        m_template = new TH2F(m_name.c_str(),m_name.c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second);
    }
    else if(nbins.size()==3)
    {
        m_template = new TH3F(m_name.c_str(),m_name.c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second,nbins[2],minmax[2].first,minmax[2].second);
    }
    m_template->Sumw2();

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
void Template::setRaw1DTemplates(const vector<TH1D*>& histos)
/*****************************************************************/
{
    for(unsigned int axis=0;axis<m_raw1DTemplates.size();axis++)
    {
        if(m_raw1DTemplates[axis])
        {
            delete m_raw1DTemplates[axis];
        }
        m_raw1DTemplates.clear();
    }
    for(unsigned int axis=0;axis<histos.size();axis++)
    {
        stringstream name;
        name << m_name << "_proj" << axis;
        m_raw1DTemplates.push_back(dynamic_cast<TH1D*>(histos[axis]->Clone(name.str().c_str())));
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
        raw1D->SetLineColor(kRed);
        raw1D->SetLineWidth(2);
        proj1D->SetLineColor(kBlack);
        proj1D->SetMarkerColor(kBlack);
        proj1D->SetMarkerStyle(20);
        raw1D->SetXTitle(getVariable(axis).c_str());
        raw1D->Draw("hist");
        proj1D->Draw("same");
        addControlPlot(c);
    }
}
