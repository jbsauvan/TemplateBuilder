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
    const Json::Value trees = tmp["trees"];
    if(files.isNull() && trees.isNull())
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
    else if(!files.isNull())
    {
        m_templates.back()->setOrigin(Template::Origin::FILES);
        for(unsigned int index = 0; index < files.size(); ++index)
        {
            m_templates.back()->addInputFileAndTree(files[index].asString(), "SelectedTree");
        }
    }
    else
    {
        m_templates.back()->setOrigin(Template::Origin::FILES);
        for(unsigned int index = 0; index < trees.size(); ++index)
        {
            string fileAndTree = trees[index].asString();
            vector<string> tokens;
            tokenize(fileAndTree, tokens, ":");
            if(tokens.size()!=2)
            {
                stringstream error;
                error << "TemplateParameters::readTemplate(): List of trees should be of the form fileName:treeName";
                throw runtime_error(error.str());
            }
            m_templates.back()->addInputFileAndTree(tokens[0], tokens[1]);
        }
    }
    if(m_templates.back()->getOrigin()==Template::Origin::FILES)
    {
        //std::string tree = tmp.get("tree", "SelectedTree").asString(); 
        const Json::Value tree = tmp["tree"];
        if(!tree.isNull())
        {
            std::string treeName = tree.asString();
            m_templates.back()->setTreeName(treeName);
        }
        const Json::Value variables = tmp["variables"];
        if(variables.isNull())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): No variables defined for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        if(variables.size()!=2 && variables.size()!=3)
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): ('"<<name<<"') Only 2 or 3 dimensions templates are possible\n";
            throw runtime_error(error.str());
        }
        vector<string> vars;
        for(unsigned int v=0;v<variables.size();v++)
        {
            vars.push_back(variables[v].asString());
        }
        m_templates.back()->setVariables(vars);
        std::string weight = tmp.get("weight", "").asString(); 
        m_templates.back()->setWeight(weight);
        bool conserveSumOfWeights = tmp.get("conserveSumOfWeights", false).asBool(); 
        m_templates.back()->setConserveSumOfWeights(conserveSumOfWeights);
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
        if(bins.size()!=3*variables.size())
        {
            stringstream error;
            error << "TemplateParameters::readTemplate(): Binning should be of the form [nbin1,min1,max1, ... ] for template '"<<name<<"'";
            throw runtime_error(error.str());
        }
        vector<unsigned int> nbins;
        vector< pair<double,double> > minmax;
        for(unsigned int v=0;v<variables.size();v++)
        {
            nbins.push_back(bins[v*3].asUInt());
            minmax.push_back( make_pair(bins[v*3+1].asDouble(),bins[v*3+2].asDouble()) );
        }
        m_templates.back()->createTemplate(nbins, minmax);
        //  what to do with overflows
        bool fillOverflows = tmp.get("filloverflows", false).asBool();
        m_templates.back()->setFillOverflows(fillOverflows);
    }

    // postprocessing 
    const Json::Value postprocess = tmp["postprocessing"];
    if(!postprocess.isNull())
    {
        for(unsigned int index = 0; index < postprocess.size(); ++index)
        {
            const Json::Value pp = postprocess[index];
            string type = pp["type"].asString();
            if(type=="smooth")
            {
                PostProcessing postproc(PostProcessing::Type::SMOOTH);
                readSmoothingParameters(pp, postproc);
                //string kernel = postproc.getParameter<string>("kernel");
                //if(kernel=="k5b" && m_templates.back()->numberOfDimensions()!=2)
                //{
                //    stringstream error;
                //    error << "TemplateParameters::readTemplate(): Smoothing kernel 'k5b' cannot be applied for "<<variables.size()<<"D template '"<<name<<"'";
                //    throw runtime_error(error.str());
                //}
                m_templates.back()->addPostProcessing(postproc);
            }
            else if(type=="mirror")
            {
                PostProcessing postproc(PostProcessing::Type::MIRROR);
                readMirrorParameters(pp, postproc);
                //unsigned int axisMirror = postproc.getParameter<unsigned int>("axis");
                //if(axisMirror>variables.size()-1)
                //{
                //    stringstream error;
                //    error << "TemplateParameters::readTemplate(): Axis mirror="<<axisMirror<<" for "<<variables.size()<<"D template '"<<name<<"'";
                //    throw runtime_error(error.str());
                //}
                m_templates.back()->addPostProcessing(postproc);
            }
            else if(type=="floor")
            {
                PostProcessing postproc(PostProcessing::Type::FLOOR);
                m_templates.back()->addPostProcessing(postproc);
            }
            else if(type=="rescale")
            {
                PostProcessing postproc(PostProcessing::Type::RESCALE);
                readRescalingParameters(pp, postproc);
                m_templates.back()->addPostProcessing(postproc);
            }
            else if(type=="reweight")
            {
                PostProcessing postproc(PostProcessing::Type::REWEIGHT);
                readReweightingParameters(pp, postproc);
                m_templates.back()->addPostProcessing(postproc);
            }
            else
            {
                stringstream error;
                error << "TemplateParameters::readTemplate(): Unknown post-processing '" << pp <<"'";
                throw runtime_error(error.str().c_str());
            }
        }
    }

}


/*****************************************************************/
void TemplateParameters::readSmoothingParameters(const Json::Value& smooth, PostProcessing& postproc)
/*****************************************************************/
{
    unsigned int entriesPerBin = smooth.get("entriesperbin", 200).asUInt();
    string kernel = smooth.get("kernel", "adaptive").asString();
    if(kernel!="adaptive" && kernel!="k5b")
    {
        stringstream error;
        error << "TemplateParameters::readSmoothingParameters(): Unknown smoothing kernel '"<<kernel<<"'";
        throw runtime_error(error.str());
    }
    postproc.addParameter("kernel", kernel);
    postproc.addParameter("entriesperbin", entriesPerBin);
}

/*****************************************************************/
void TemplateParameters::readMirrorParameters(const Json::Value& mirror, PostProcessing& postproc)
/*****************************************************************/
{
    bool antiMirror = mirror.get("antisymmetric", false).asBool();
    unsigned int axis = mirror.get("axis", 1).asUInt();
    postproc.addParameter("antisymmetric", antiMirror);
    postproc.addParameter("axis", axis);
}

/*****************************************************************/
void TemplateParameters::readRescalingParameters(const Json::Value& rescaling, PostProcessing& postproc)
/*****************************************************************/
{
    double factor = rescaling.get("factor", 1.).asDouble();
    postproc.addParameter("factor", factor);
}

/*****************************************************************/
void TemplateParameters::readReweightingParameters(const Json::Value& reweighting, PostProcessing& postproc)
/*****************************************************************/
{
    vector<unsigned int> axes;
    const Json::Value ax = reweighting["axes"];
    if(!ax.isNull())
    {
        for(unsigned int index = 0; index < ax.size(); ++index)
        {
            unsigned int axis = ax[index].asUInt();
            axes.push_back(axis);
        }
    }
    vector< vector<double> > binss;
    const Json::Value binnings = reweighting["rebinning"];
    if(!binnings.isNull())
    {
        for(unsigned int index = 0; index < binnings.size(); ++index)
        {
            vector< double > bins;
            const Json::Value binning = binnings[index];
            for(unsigned int b = 0; b < binning.size(); ++b)
            {
                double boundary = binning[b].asDouble();
                bins.push_back(boundary);
            }
            binss.push_back(bins);
        }
    }
    postproc.addParameter("axes", axes);
    postproc.addParameter("rebinning", binss);
}

/*****************************************************************/
void TemplateParameters::tokenize(const string& str,
        vector<string>& tokens,
        const string& delimiter)
/*****************************************************************/
{
    string::size_type length = delimiter.size();
    string::size_type lastPos = 0;
    string::size_type pos     = str.find(delimiter, 0);


    while (string::npos != pos)
    {
        // Found a token, add it to the vector.
        if(str.substr(lastPos, pos - lastPos).size()>0)
            tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + length;
        // Find next "non-delimiter"
        pos = str.find(delimiter, lastPos);
    }
    if(str.substr(lastPos).size()>0)
        tokens.push_back(str.substr(lastPos));
}
