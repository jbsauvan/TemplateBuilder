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

#include <iostream>
#include <sstream>

using namespace std;


/*****************************************************************/
Template::Template():m_template(NULL),
    m_widthx(NULL),
    m_widthy(NULL)
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
    if(m_widthx)
    {
        delete m_widthx;
        m_widthx = NULL;
    }
    if(m_widthy)
    {
        delete m_widthy;
        m_widthy = NULL;
    }
}

/*****************************************************************/
const string& Template::getVariable(unsigned int index) const
/*****************************************************************/
{
    assert(index<=1);
    return m_variables.at(index);
}

/*****************************************************************/
void Template::setVariables(const string& var1, const string& var2)
/*****************************************************************/
{
    m_variables.clear();
    m_variables.push_back(var1);
    m_variables.push_back(var2);
}


/*****************************************************************/
void Template::createTemplate(unsigned int nbin1, double min1, double max1, unsigned int nbin2, double min2, double max2)
/*****************************************************************/
{
    if(m_template)
        m_template->Delete();
    if(m_widthx)
        m_widthx->Delete();
    if(m_widthy)
        m_widthy->Delete();
    m_template = new TH2F(m_name.c_str(),m_name.c_str(),nbin1,min1,max1,nbin2,min2,max2);
    m_template->Sumw2();
    stringstream nameWidthx;
    nameWidthx << m_name << "_widthx";
    m_widthx = new TH2F(nameWidthx.str().c_str(),nameWidthx.str().c_str(),nbin1,min1,max1,nbin2,min2,max2);
    m_widthx->Sumw2();
    stringstream nameWidthy;
    nameWidthy << m_name << "_widthy";
    m_widthy = new TH2F(nameWidthy.str().c_str(),nameWidthy.str().c_str(),nbin1,min1,max1,nbin2,min2,max2);
    m_widthy->Sumw2();
}

/*****************************************************************/
void Template::setTemplate(const TH2F* histo)
/*****************************************************************/
{
    if(m_template)
        m_template->Delete();
    m_template = dynamic_cast<TH2F*>(histo->Clone(m_name.c_str()));
    //m_template->Sumw2();
}

/*****************************************************************/
void Template::setWidths(const TH2F* widthx, const TH2F* widthy)
/*****************************************************************/
{
    if(m_widthx)
        m_widthx->Delete();
    if(m_widthy)
        m_widthy->Delete();
    stringstream nameWidthx;
    nameWidthx << m_name << "_widthx";
    m_widthx = dynamic_cast<TH2F*>(widthx->Clone(nameWidthx.str().c_str()));
    //m_widthx->Sumw2();
    stringstream nameWidthy;
    nameWidthy << m_name << "_widthy";
    m_widthy = dynamic_cast<TH2F*>(widthy->Clone(nameWidthy.str().c_str()));
    //m_widthy->Sumw2();
}

/*****************************************************************/
void Template::store(double v1, double v2, double w)
/*****************************************************************/
{
    vector<double> entry;
    entry.push_back(v1);
    entry.push_back(v2);
    entry.push_back(w);
    m_entries.push_back(entry);
}
