/**
 *  @file  TemplateParameters.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/08/2013
 *
 *  @internal
 *     Created :  12/08/2013
 * Last update :  12/08/2013 05:52:48 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TEMPLATEPARAMETERS_H
#define TEMPLATEPARAMETERS_H

#include "Template.h"

#include <string>

namespace Json
{
    class Value;
};

class TemplateParameters
{
    public:
        TemplateParameters(){};
        ~TemplateParameters(){};

        void read(const std::string& parFile);

        const std::string& inputDirectory() const {return m_inputDirectory;}
        const std::string& outputFileName() const {return m_outputFileName;}
        std::vector<Template*>::iterator templateBegin() {return m_templates.begin();}
        std::vector<Template*>::iterator templateEnd() {return m_templates.end();}

    private:
        void readTemplate(const Json::Value& tmp);
        void readSmoothingParameters(const Json::Value& smooth, PostProcessing& postproc);
        void readMirrorParameters(const Json::Value& mirror, PostProcessing& postproc);
        void readRescalingParameters(const Json::Value& rescaling, PostProcessing& postproc);
        void readReweightingParameters(const Json::Value& reweighting, PostProcessing& postproc);
        void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiter = " ");

        std::string m_inputDirectory;
        std::string m_outputFileName;
        std::vector<Template*> m_templates;


};


#endif
