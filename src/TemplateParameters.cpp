/**
 *  @file  TemplateParameters.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 06:05:40 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include "TemplateParameters.h"

#include "json/json.h"

#include <sstream>
#include <fstream>
#include <stdexcept>

using namespace std;

/*****************************************************************/
void TemplateParameters::read(const string& parFile)
/*****************************************************************/
{

    Json::Value root;   
    Json::Reader reader;
    std::ifstream pars(parFile, std::ifstream::binary);
    if(!pars.is_open())
    {
        stringstream error;
        error << "TemplateParameters::read(): Cannot open parameter file '"<<parFile<<"'\n";
        throw runtime_error(error.str());
    }
    bool parsingSuccessful = reader.parse( pars, root );
    if(!parsingSuccessful)
    {
        stringstream error;
        error << "TemplateParameters::read(): Problem while reading parameter file. Parser error follows. \n"<<reader.getFormatedErrorMessages()<<"\n";
        throw runtime_error(error.str());
    }

    m_inputDirectory = root.get("inputDirectory", "./" ).asString();
    m_outputFileName = root.get("outputFile", "templates.root" ).asString();

    const Json::Value templates = root["templates"];
    if(templates.isNull())
    {
        stringstream error;
        error << "TemplateParameters::read(): No template defined";
        throw runtime_error(error.str());
    }
    for(unsigned int index = 0; index < templates.size(); ++index)
    {
        readTemplate(templates[index]);
    }
    pars.close();
}


/*****************************************************************/
void TemplateParameters::readTemplate(const Json::Value& tmp)
/*****************************************************************/
{
    m_templates.push_back(new Template());
    std::string name = tmp.get("name", "template").asString();   
    m_templates.back()->setName(name);
    const Json::Value files = tmp["files"];
    if(files.isNull())
    {
        const Json::Value tmpSum = tmp["templatesum"];
        if(tmpSum.isNull())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): No file nor template sum defined for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        m_templates.back()->setOrigin(Template::Origin::TEMPLATES);
        for(unsigned int index = 0; index < tmpSum.size(); ++index)
        {
            string name = tmpSum[index].get("name", "template").asString();
            double factor = tmpSum[index].get("factor", "1.").asDouble();
            m_templates.back()->addInputTemplate(name, factor);
        }
        
    }
    else
    {
        m_templates.back()->setOrigin(Template::Origin::FILES);
        for(unsigned int index = 0; index < files.size(); ++index)
        {
            m_templates.back()->addInputFile(files[index].asString());
        }
    }
    if(m_templates.back()->getOrigin()==Template::Origin::FILES)
    {
        std::string tree = tmp.get("tree", "SelectedTree").asString(); 
        m_templates.back()->setTreeName(tree);
        const Json::Value variables = tmp["variables"];
        if(variables.isNull())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): No variables defined for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        if(variables.size()!=2)
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): Number of variables !=2 for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        m_templates.back()->setVariables(variables[(unsigned int)0].asString(),variables[(unsigned int)1].asString());
        std::string weight = tmp.get("weight", "").asString(); 
        m_templates.back()->setWeight(weight);
        std::string selection = tmp.get("selection", "").asString(); 
        m_templates.back()->setSelection(selection);
        std::string assertion = tmp.get("assertion", "1").asString(); 
        m_templates.back()->setAssertion(assertion);


        const Json::Value binning = tmp["binning"];
        if(binning.isNull())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): No binning defined for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        std::string binningType = binning.get("type", "fixed").asString();
        Template::BinningType type = Template::BinningType::FIXED;
        if(binningType=="adaptive") type = Template::BinningType::ADAPTIVE;
        m_templates.back()->setBinningType( type );
        unsigned int entriesPerBin = binning.get("entriesperbin", 200).asUInt();
        m_templates.back()->setEntriesPerBin( entriesPerBin );
        const Json::Value bins = binning["bins"]; 
        if(bins.isNull())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): No binning defined for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        if(bins.size()!=6)
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): Binning should be of the form [nbin1,min1,max1,nbin2,min2,max2] for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        m_templates.back()->createTemplate(
                bins[(unsigned int)0].asUInt(),
                bins[(unsigned int)1].asDouble(),
                bins[(unsigned int)2].asDouble(),
                bins[(unsigned int)3].asUInt(),
                bins[(unsigned int)4].asDouble(),
                bins[(unsigned int)5].asDouble()
                );
    }

    // postprocessing TODO: add possibility to specify parameters
    const Json::Value postprocess = tmp["postprocessing"];
    if(!postprocess.isNull())
    {
        for(unsigned int index = 0; index < postprocess.size(); ++index)
        {
            string pp = postprocess[index].asString();
            if(pp=="smooth_k5b") m_templates.back()->addPostProcessing(Template::PostProcessing::SMOOTH_K5B);
            else if(pp=="smooth_adaptive") m_templates.back()->addPostProcessing(Template::PostProcessing::SMOOTH_ADAPTIVE);
            else if(pp=="mirror") m_templates.back()->addPostProcessing(Template::PostProcessing::MIRROR);
            else if(pp=="mirror_inverse") m_templates.back()->addPostProcessing(Template::PostProcessing::MIRROR_INV);
            else if(pp=="floor") m_templates.back()->addPostProcessing(Template::PostProcessing::FLOOR);
            else
            {
                stringstream error;
                error << "TemplateParameters::readTemplate(): Unknown post-processing '" << pp <<"'";
                throw runtime_error(error.str().c_str());
            }
        }
    }

    // Rescaling TODO: maybe merge this with postprocessing
    double scaleFactor = tmp.get("rescaling", 1.).asDouble();
    m_templates.back()->setRescaling( scaleFactor );

}
