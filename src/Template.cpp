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
#include <sstream>

using namespace std;


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
        m_widths.clear();
    }
}

/*****************************************************************/
const string& Template::getVariable(unsigned int index) const
/*****************************************************************/
{
    assert(index<3);
    return m_variables.at(index);
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
        stringstream nameWidth;
        nameWidth << m_name << "_width" << axis;
        if(nbins.size()==2)
        {
            m_widths.push_back(new TH2F(nameWidth.str().c_str(),nameWidth.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second));
        }
        else if(nbins.size()==3)
        {
            m_widths.push_back(new TH3F(nameWidth.str().c_str(),nameWidth.str().c_str(),nbins[0],minmax[0].first,minmax[0].second,nbins[1],minmax[1].first,minmax[1].second,nbins[2],minmax[2].first,minmax[2].second));
        }
        m_widths.back()->Sumw2();
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
